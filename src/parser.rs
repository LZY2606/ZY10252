//! 受限本地格式（`.rxn` DSL）解析器。
//!
//! 行格式（`#` 起始为注释，空行忽略）：
//! ```text
//! molecule <id>
//!   atom <元素|D|T> [质量数] [形式电荷]
//!   bond <u> <v> [键级1..3，默认1]
//! end
//! reaction <id>
//!   reactant <分子id> [系数，默认1，可写 n/d]
//!   product  <分子id> [系数，可写 n/d]
//!   map
//!     <反应物分子>:<原子下标> -> <产物分子>:<原子下标>
//!   end
//!   direction unknown          # 可选：方向未定
//! end
//! network <id>
//!   step <反应id> [unknown]
//! end
//! ```

use crate::chem::{builtin_molecules, parse_ratio, Atom, Molecule, Rat};
use crate::model::{Dataset, MapEdge, Network, Reaction, Term};
use std::collections::BTreeMap;

/// 缩进块内的一行：(原始行号, 去除注释与空白的内容)。
type BlockLine<'a> = (usize, &'a str);

pub fn parse(input: &str) -> Result<Dataset, String> {
    let mut molecules: BTreeMap<String, Molecule> = builtin_molecules()
        .into_iter()
        .map(|m| (m.id.clone(), m))
        .collect();
    let mut reactions = Vec::new();
    let mut networks = Vec::new();

    let lines: Vec<&str> = input.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let raw = lines[i];
        let line = strip_comment(raw).trim();
        i += 1;
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        match parts.next() {
            Some("molecule") => {
                let id = parts.next().ok_or("molecule 缺少 id")?;
                let (mol, ni) = parse_molecule(id, &lines, i)?;
                i = ni;
                molecules.insert(mol.id.clone(), mol);
            }
            Some("reaction") => {
                let id = parts.next().ok_or("reaction 缺少 id")?;
                let (rxn, ni) = parse_reaction(id, &lines, i)?;
                i = ni;
                reactions.push(rxn);
            }
            Some("network") => {
                let id = parts.next().ok_or("network 缺少 id")?;
                let (net, ni) = parse_network(id, &lines, i)?;
                i = ni;
                networks.push(net);
            }
            other => return Err(format!("第{i}行: 无法识别的声明 {other:?}")),
        }
    }

    validate_refs(&molecules, &reactions, &networks)?;
    Ok(Dataset {
        molecules,
        reactions,
        networks,
    })
}

fn strip_comment(line: &str) -> &str {
    match line.find('#') {
        Some(p) => &line[..p],
        None => line,
    }
}

/// 读取缩进块直到 `end`；`start` 为首行的下一行，返回(块内行, end 后下一行)。
#[allow(clippy::type_complexity)]
fn read_block<'a>(lines: &[&'a str], start: usize) -> Result<(Vec<BlockLine<'a>>, usize), String> {
    let mut body: Vec<BlockLine<'a>> = Vec::new();
    let mut i = start;
    while i < lines.len() {
        let raw = strip_comment(lines[i]);
        let trimmed = raw.trim();
        if trimmed == "end" {
            return Ok((body, i + 1));
        }
        if !trimmed.is_empty() {
            body.push((i + 1, trimmed));
        }
        i += 1;
    }
    Err("块缺少 end".to_string())
}

fn parse_molecule(id: &str, lines: &[&str], start: usize) -> Result<(Molecule, usize), String> {
    let (body, next) = read_block(lines, start)?;
    let mut atoms: Vec<Atom> = Vec::new();
    let mut bonds = Vec::new();
    for (ln, text) in body {
        let mut it = text.split_whitespace();
        match it.next() {
            Some("atom") => {
                let element = it
                    .next()
                    .ok_or_else(|| format!("第{ln}行: atom 缺少元素"))?;
                // 字段顺序：元素 [质量数] [形式电荷]；省略的为 0。
                let mut isotope: i64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                let charge: i64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
                if isotope == 0 && (element == "D" || element == "T") {
                    isotope = if element == "D" { 2 } else { 3 };
                }
                atoms.push(Atom {
                    element: element.to_string(),
                    isotope,
                    charge,
                });
            }
            Some("bond") => {
                let u = it
                    .next()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| format!("第{ln}行: bond 参数错误"))?;
                let v = it
                    .next()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| format!("第{ln}行: bond 参数错误"))?;
                let order: u8 = it.next().and_then(|s| s.parse().ok()).unwrap_or(1);
                bonds.push((u, v, order));
            }
            other => return Err(format!("第{ln}行: molecule 内非法声明 {other:?}")),
        }
    }
    for &(u, v, _) in &bonds {
        if u >= atoms.len() || v >= atoms.len() {
            return Err(format!("molecule {id}: 边引用了不存在的原子 {u}/{v}"));
        }
    }
    Ok((
        Molecule {
            id: id.to_string(),
            atoms,
            bonds,
        },
        next,
    ))
}

fn parse_edge(ln: usize, text: &str) -> Result<MapEdge, String> {
    let (left, right) = text
        .split_once("->")
        .ok_or_else(|| format!("第{ln}行: map 边需要 'from -> to'，得到 {text:?}"))?;
    let parse_atom_ref = |s: &str| -> Result<(String, usize), String> {
        let (mol, idx) = s
            .trim()
            .split_once(':')
            .ok_or_else(|| format!("第{ln}行: 原子引用需要 '分子:下标'，得到 {s:?}"))?;
        let idx: usize = idx
            .trim()
            .parse()
            .map_err(|_| format!("第{ln}行: 原子下标非法 {idx:?}"))?;
        Ok((mol.trim().to_string(), idx))
    };
    let (from_mol, from_atom) = parse_atom_ref(left)?;
    let (to_mol, to_atom) = parse_atom_ref(right)?;
    Ok(MapEdge {
        from_mol,
        from_atom,
        to_mol,
        to_atom,
    })
}

fn parse_reaction(id: &str, lines: &[&str], start: usize) -> Result<(Reaction, usize), String> {
    // 直接按全局行号扫描，支持 map ... end 子块嵌套。
    let mut reactants: Vec<Term> = Vec::new();
    let mut products: Vec<Term> = Vec::new();
    let mut maps: Vec<Vec<MapEdge>> = Vec::new();
    let mut direction_unknown = false;

    let mut i = start;
    let mut depth = 1usize; // 已消费 reaction 首行
    while i < lines.len() {
        let trimmed = strip_comment(lines[i]).trim().to_string();
        i += 1;
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "end" {
            depth -= 1;
            if depth == 0 {
                return Ok((
                    Reaction {
                        id: id.to_string(),
                        reactants,
                        products,
                        maps,
                        direction_unknown,
                    },
                    i,
                ));
            }
            continue;
        }
        if trimmed.starts_with("map") {
            depth += 1;
            let mut edges = Vec::new();
            while i < lines.len() {
                let inner = strip_comment(lines[i]).trim().to_string();
                i += 1;
                if inner.is_empty() {
                    continue;
                }
                if inner == "end" {
                    depth -= 1;
                    break;
                }
                edges.push(parse_edge(i, &inner)?);
            }
            maps.push(edges);
            continue;
        }
        let ln = i;
        let mut it = trimmed.split_whitespace();
        match it.next() {
            Some("reactant") | Some("product") => {
                let is_product = trimmed.starts_with("product");
                let mol = it
                    .next()
                    .ok_or_else(|| format!("第{ln}行: 缺少分子 id"))?
                    .to_string();
                let coef = match it.next() {
                    Some(c) => parse_ratio(c).map_err(|e| format!("第{ln}行: {e}"))?,
                    None => Rat::from_integer(1),
                };
                if coef <= Rat::from_integer(0) {
                    return Err(format!("第{ln}行: 计量系数必须为正"));
                }
                let term = Term { mol, coef };
                if is_product {
                    products.push(term);
                } else {
                    reactants.push(term);
                }
            }
            Some("direction") => {
                let val = it
                    .next()
                    .ok_or_else(|| format!("第{ln}行: direction 缺少值"))?;
                if val != "unknown" {
                    return Err(format!("第{ln}行: direction 仅支持 unknown"));
                }
                direction_unknown = true;
            }
            Some("molecule" | "reaction" | "network") => {
                return Err(format!("第{ln}行: reaction 块内出现未闭合的新声明"));
            }
            other => return Err(format!("第{ln}行: reaction 内非法声明 {other:?}")),
        }
    }
    Err(format!("reaction {id}: 缺少 end"))
}

fn parse_network(id: &str, lines: &[&str], start: usize) -> Result<(Network, usize), String> {
    let (body, next) = read_block(lines, start)?;
    let mut steps = Vec::new();
    for (ln, text) in body {
        let mut it = text.split_whitespace();
        match it.next() {
            Some("step") => {
                let rid = it
                    .next()
                    .ok_or_else(|| format!("第{ln}行: step 缺少反应 id"))?;
                let unknown = it.next() == Some("unknown");
                steps.push((rid.to_string(), unknown));
            }
            other => return Err(format!("第{ln}行: network 内非法声明 {other:?}")),
        }
    }
    Ok((
        Network {
            id: id.to_string(),
            steps,
        },
        next,
    ))
}

fn validate_refs(
    molecules: &BTreeMap<String, Molecule>,
    reactions: &[Reaction],
    networks: &[Network],
) -> Result<(), String> {
    for rxn in reactions {
        for t in rxn.reactants.iter().chain(&rxn.products) {
            if !molecules.contains_key(&t.mol) {
                return Err(format!("reaction {}: 未知分子 {}", rxn.id, t.mol));
            }
        }
        for (mi, edges) in rxn.maps.iter().enumerate() {
            for e in edges {
                let rm = molecules.get(&e.from_mol).ok_or_else(|| {
                    format!(
                        "reaction {} map#{}: 未知反应物分子 {}",
                        rxn.id, mi, e.from_mol
                    )
                })?;
                let pm = molecules.get(&e.to_mol).ok_or_else(|| {
                    format!("reaction {} map#{}: 未知产物分子 {}", rxn.id, mi, e.to_mol)
                })?;
                if e.from_atom >= rm.atoms.len() {
                    return Err(format!(
                        "reaction {} map#{}: {} 原子下标越界",
                        rxn.id, mi, e.from_mol
                    ));
                }
                if e.to_atom >= pm.atoms.len() {
                    return Err(format!(
                        "reaction {} map#{}: {} 原子下标越界",
                        rxn.id, mi, e.to_mol
                    ));
                }
                let ra = &rm.atoms[e.from_atom];
                let pa = &pm.atoms[e.to_atom];
                // 原子身份只要求元素一致；同位素/电荷是否守恒由守恒向量独立核对。
                // （如 H 映射到 H@2 表示该氢原子获得同位素标签。）
                if ra.element() != pa.element() {
                    return Err(format!(
                        "reaction {} map#{}: {}:{} -> {}:{} 元素身份不一致",
                        rxn.id, mi, e.from_mol, e.from_atom, e.to_mol, e.to_atom
                    ));
                }
            }
        }
    }
    for net in networks {
        for (rid, _) in &net.steps {
            if !reactions.iter().any(|r| &r.id == rid) {
                return Err(format!("network {}: 未知反应 {rid}", net.id));
            }
        }
    }
    Ok(())
}
