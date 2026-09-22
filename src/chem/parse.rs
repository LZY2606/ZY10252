//! 受限本地格式解析器。
//!
//! ```text
//! MOL id 名称
//!   atom 1 [13C]
//!   bond 1 2 2
//! RXN id
//!   reactant 1/2 O2-18 + D2
//!   product D2O-18
//!   map a1 b1
//! NET id 名称
//!   reaction R1 ->
//!   boundary CO2 1
//! ```

use super::*;
use std::collections::BTreeSet;

#[derive(Debug)]
pub struct ParseError {
    pub line: usize,
    pub msg: String,
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "第 {} 行: {}", self.line, self.msg)
    }
}
impl std::error::Error for ParseError {}

fn err<T>(line: usize, msg: impl Into<String>) -> Result<T, ParseError> {
    Err(ParseError {
        line,
        msg: msg.into(),
    })
}

/// 解析原子令牌：`H`、`[D]`、`[13C]`、`H+`、`[18O-]`、`Mg2+`、`2H-` 视为 [2H-]。
pub fn parse_atom_token(tok: &str) -> Option<AtomKey> {
    let bytes: Vec<char> = tok.trim().chars().collect();
    if bytes.is_empty() {
        return None;
    }
    let mut i = 0usize;
    let mut isotope: Option<u32> = None;
    if bytes[i].is_ascii_digit() {
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        isotope = Some(tok[start..i].parse().ok()?);
    }
    if i >= bytes.len() || !bytes[i].is_ascii_uppercase() {
        return None;
    }
    let elem_start = i;
    i += 1;
    if i < bytes.len() && bytes[i].is_ascii_lowercase() {
        i += 1;
    }
    let element: String = tok[elem_start..i].chars().collect();
    let mut charge = 0i32;
    let rest: String = tok[i..].chars().collect();
    if !rest.is_empty() {
        let (sign, digits) = if let Some(n) = rest.strip_suffix('+') {
            (1, n)
        } else if let Some(n) = rest.strip_suffix('-') {
            (-1, n)
        } else {
            return None;
        };
        let mag: i32 = if digits.is_empty() {
            1
        } else {
            digits.parse().ok()?
        };
        charge = sign * mag;
    }
    Some(AtomKey::new(&element, isotope, charge))
}

fn parse_rat(s: &str) -> Option<Rat> {
    let s = s.trim();
    if let Some((a, b)) = s.split_once('/') {
        let a: i64 = a.parse().ok()?;
        let b: i64 = b.parse().ok()?;
        if b == 0 {
            return None;
        }
        Some(Ratio::new(a, b))
    } else {
        Some(Ratio::from_integer(s.parse().ok()?))
    }
}

/// `1/2 O2-18 + D2` → 项列表；分子 id 允许字母、数字、`-`。
fn parse_term_list(text: &str) -> Result<Vec<Term>, String> {
    let mut terms = Vec::new();
    for piece in text.split('+') {
        let piece = piece.trim();
        if piece.is_empty() {
            return Err("空的参与项".into());
        }
        let (coef_tok, rest) = match piece.split_once(char::is_whitespace) {
            Some((c, r)) => (c, r.trim()),
            None => ("1", piece),
        };
        let coef = parse_rat(coef_tok)
            .ok_or_else(|| format!("无法解析系数 `{coef_tok}`"))?;
        if coef <= Rat::zero() {
            return Err(format!("系数必须为正: `{coef_tok}`"));
        }
        let mol_id = rest.split_whitespace().collect::<Vec<_>>().join(" ");
        if mol_id.is_empty() || !mol_id.chars().all(allowed_id_char) {
            return Err(format!("非法分子 id `{mol_id}`"));
        }
        terms.push(Term { mol_id, coef });
    }
    Ok(terms)
}

fn allowed_id_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '-'
}

fn parse_bond_order(s: &str) -> Result<BondOrder, String> {
    match s {
        "1" => Ok(BondOrder::Single),
        "2" => Ok(BondOrder::Double),
        "3" => Ok(BondOrder::Triple),
        other => Err(format!("键级必须为 1/2/3，得到 `{other}`")),
    }
}

enum Block {
    None,
    Mol(String, String),
    Rxn(String),
    Net(String, String),
}

/// 解析完整 fixture 文本。
pub fn parse_dataset(text: &str) -> Result<Dataset, ParseError> {
    let mut ds = Dataset::default();
    let mut block = Block::None;
    let mut pending_mol: Option<Molecule> = None;
    let mut pending_rxn: Option<Record> = None;
    let mut pending_net: Option<(Network, BTreeSet<String>, BTreeSet<String>)> = None;

    let finish = |ds: &mut Dataset,
                  block: &mut Block,
                  pm: &mut Option<Molecule>,
                  pr: &mut Option<Record>,
                  pn: &mut Option<(Network, BTreeSet<String>, BTreeSet<String>)>|
     -> Result<(), ParseError> {
        match block {
            Block::Mol(_, _) => {
                if let Some(m) = pm.take() {
                    if ds.molecules.insert(m.id.clone(), m).is_some() {
                        return Err(ParseError {
                            line: 0,
                            msg: "重复的 MOL id".into(),
                        });
                    }
                }
            }
            Block::Rxn(_) => {
                if let Some(r) = pr.take() {
                    if ds.reactions.insert(r.id.clone(), r).is_some() {
                        return Err(ParseError {
                            line: 0,
                            msg: "重复的 RXN id".into(),
                        });
                    }
                }
            }
            Block::Net(_, _) => {
                if let Some((mut n, rx, du)) = pn.take() {
                    n.reaction_ids = rx.into_iter().collect();
                    n.direction_unknown = du;
                    if ds.networks.insert(n.id.clone(), n).is_some() {
                        return Err(ParseError {
                            line: 0,
                            msg: "重复的 NET id".into(),
                        });
                    }
                }
            }
            Block::None => {}
        }
        *block = Block::None;
        Ok(())
    };

    for (idx, raw_line) in text.lines().enumerate() {
        let lineno = idx + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (kw, rest) = match line.split_once(char::is_whitespace) {
            Some((k, r)) => (k, r.trim()),
            None => (line, ""),
        };
        match kw {
            "MOL" => {
                finish(&mut ds, &mut block, &mut pending_mol, &mut pending_rxn, &mut pending_net)
                    .map_err(|e| ParseError {
                        line: lineno,
                        msg: e.msg,
                    })?;
                let mut parts = rest.splitn(2, char::is_whitespace);
                let id = parts.next().unwrap_or("").trim().to_string();
                let name = parts.next().unwrap_or(&id).trim().to_string();
                if id.is_empty() || !id.chars().all(allowed_id_char) {
                    return err(lineno, "非法 MOL id");
                }
                pending_mol = Some(Molecule {
                    id: id.clone(),
                    name,
                    atoms: Vec::new(),
                    bonds: Vec::new(),
                });
                block = Block::Mol(id, String::new());
            }
            "RXN" => {
                finish(&mut ds, &mut block, &mut pending_mol, &mut pending_rxn, &mut pending_net)
                    .map_err(|e| ParseError {
                        line: lineno,
                        msg: e.msg,
                    })?;
                let id = rest.trim().to_string();
                if id.is_empty() || !id.chars().all(allowed_id_char) {
                    return err(lineno, "非法 RXN id");
                }
                pending_rxn = Some(Record {
                    id: id.clone(),
                    ..Default::default()
                });
                block = Block::Rxn(id);
            }
            "NET" => {
                finish(&mut ds, &mut block, &mut pending_mol, &mut pending_rxn, &mut pending_net)
                    .map_err(|e| ParseError {
                        line: lineno,
                        msg: e.msg,
                    })?;
                let mut parts = rest.splitn(2, char::is_whitespace);
                let id = parts.next().unwrap_or("").trim().to_string();
                let name = parts.next().unwrap_or(&id).trim().to_string();
                pending_net = Some((
                    Network {
                        id: id.clone(),
                        name,
                        ..Default::default()
                    },
                    BTreeSet::new(),
                    BTreeSet::new(),
                ));
                block = Block::Net(id, String::new());
            }
            "atom" => {
                let mol = pending_mol
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "atom 必须位于 MOL 块".into() })?;
                let mut it = rest.split_whitespace();
                let id_tok = it.next().unwrap_or("");
                let key_tok = it.next().unwrap_or("");
                let expect = mol.atoms.len() + 1;
                if id_tok != expect.to_string() {
                    return err(lineno, format!("原子编号必须连续从 1 开始，期望 {expect}"));
                }
                let key = parse_atom_token(key_tok)
                    .ok_or_else(|| ParseError { line: lineno, msg: format!("非法原子 `{key_tok}`") })?;
                mol.atoms.push(Atom { key });
            }
            "bond" => {
                let mol = pending_mol
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "bond 必须位于 MOL 块".into() })?;
                let parts: Vec<&str> = rest.split_whitespace().collect();
                if parts.len() != 3 {
                    return err(lineno, "bond 需要三个字段: a b order");
                }
                let a: usize = parts[0].parse().map_err(|_| ParseError { line: lineno, msg: "非法原子编号".into() })?;
                let b: usize = parts[1].parse().map_err(|_| ParseError { line: lineno, msg: "非法原子编号".into() })?;
                let order = parse_bond_order(parts[2]).map_err(|m| ParseError { line: lineno, msg: m })?;
                if a == 0 || b == 0 || a > mol.atoms.len() || b > mol.atoms.len() || a == b {
                    return err(lineno, "键引用了不存在或相同的原子");
                }
                mol.bonds.push(Bond { a: a - 1, b: b - 1, order });
            }
            "reactant" | "reactants" => {
                let r = pending_rxn
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "reactant 必须位于 RXN 块".into() })?;
                r.reactants = parse_term_list(rest).map_err(|m| ParseError { line: lineno, msg: m })?;
            }
            "product" | "products" => {
                let r = pending_rxn
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "product 必须位于 RXN 块".into() })?;
                r.products = parse_term_list(rest).map_err(|m| ParseError { line: lineno, msg: m })?;
            }
            "map" => {
                let r = pending_rxn
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "map 必须位于 RXN 块".into() })?;
                let parts: Vec<&str> = rest.split_whitespace().collect();
                if parts.len() != 2 {
                    return err(lineno, "map 需要两个原子槽标识");
                }
                r.declared_maps.push((parts[0].to_string(), parts[1].to_string()));
            }
            "direction" => {
                let r = pending_rxn
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "direction 必须位于 RXN 块".into() })?;
                r.direction_unknown = rest == "?" || rest == "unknown";
            }
            "reaction" => {
                let pn = pending_net
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "reaction 必须位于 NET 块".into() })?;
                let parts: Vec<&str> = rest.split_whitespace().collect();
                if parts.is_empty() {
                    return err(lineno, "reaction 缺少反应 id");
                }
                let id = parts[0].to_string();
                pn.1.insert(id.clone());
                if parts.get(1).is_some_and(|d| *d == "?") {
                    pn.2.insert(id);
                }
            }
            "boundary" => {
                let pn = pending_net
                    .as_mut()
                    .ok_or_else(|| ParseError { line: lineno, msg: "boundary 必须位于 NET 块".into() })?;
                let mut it = rest.split_whitespace();
                let mol = it.next().unwrap_or("");
                let val = it.next().unwrap_or("");
                let v = parse_rat(val).ok_or_else(|| ParseError { line: lineno, msg: "非法边界通量".into() })?;
                if !mol.chars().all(allowed_id_char) {
                    return err(lineno, "非法分子 id");
                }
                pn.0.boundary.insert(mol.to_string(), v);
            }
            other => {
                return err(lineno, format!("未知关键字 `{other}`"));
            }
        }
    }
    finish(&mut ds, &mut block, &mut pending_mol, &mut pending_rxn, &mut pending_net)
        .map_err(|e| ParseError { line: text.lines().count(), msg: e.msg })?;
    validate_dataset(&ds)?;
    Ok(ds)
}

fn validate_dataset(ds: &Dataset) -> Result<(), ParseError> {
    for mol in ds.molecules.values() {
        let mut seen = BTreeSet::new();
        for b in &mol.bonds {
            let key = if b.a < b.b { (b.a, b.b) } else { (b.b, b.a) };
            if !seen.insert(key) {
                return Err(ParseError { line: 0, msg: format!("{} 中存在重复键", mol.id) });
            }
        }
    }
    for rec in ds.reactions.values() {
        for t in rec.reactants.iter().chain(&rec.products) {
            if !ds.molecules.contains_key(&t.mol_id) {
                return Err(ParseError {
                    line: 0,
                    msg: format!("反应 {} 引用了未定义分子 {}", rec.id, t.mol_id),
                });
            }
        }
        if rec.reactants.is_empty() || rec.products.is_empty() {
            return Err(ParseError { line: 0, msg: format!("反应 {} 两侧均不可为空", rec.id) });
        }
    }
    for net in ds.networks.values() {
        for id in &net.reaction_ids {
            if !ds.reactions.contains_key(id) {
                return Err(ParseError { line: 0, msg: format!("网络 {} 引用未定义反应 {}", net.id, id) });
            }
        }
        for m in net.boundary.keys() {
            if !ds.molecules.contains_key(m) {
                return Err(ParseError { line: 0, msg: format!("网络 {} 边界分子 {} 未定义", net.id, m) });
            }
        }
    }
    Ok(())
}
