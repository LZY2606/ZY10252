//! 化学核心：原子、分子图、有理数、自同构（分子对称性）。
//!
//! 数据口径：
//! - 所有计量/守恒计算使用 [`num_rational::Ratio<i64>`]，不用浮点，没有“容差放过一个原子”。
//! - 原子的语义三元组 = (元素, 质量数, 形式电荷)；同位素通过质量数区分（0 表示未标记/天然丰度）。
//! - 分子总电荷 = 各原子形式电荷之和。`e` 是伪元素，用于电子（电荷 -1）。

use num_rational::Ratio;
use std::collections::BTreeMap;
use std::sync::OnceLock;

pub type Rat = Ratio<i64>;

/// 有理数转字符串：整数不带分母，分数形如 `1/4`。
pub fn rat_string(r: Rat) -> String {
    if *r.denom() == 1 {
        r.numer().to_string()
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Atom {
    /// 元素符号（如 H/D 保留原写法用于展示；归一化见 [`Atom::element`]）。
    pub element: String,
    /// 质量数（同位素标签）。0 = 未标记。
    pub isotope: i64,
    /// 形式电荷。
    pub charge: i64,
}

impl Atom {
    /// 原子身份元素：D/T 在身份上属于 H（同位素差异不改变元素身份）。
    pub fn element(&self) -> &str {
        match self.element.as_str() {
            "D" | "T" => "H",
            other => other,
        }
    }
    /// 元素+同位素分量键（不含电荷）。D/T 归一为 H 的同位素分量。
    pub fn composition_key(&self) -> String {
        let el = match self.element.as_str() {
            "D" | "T" => "H",
            other => other,
        };
        if self.isotope != 0 {
            format!("{el}@{}", self.isotope)
        } else {
            el.to_string()
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Molecule {
    pub id: String,
    pub atoms: Vec<Atom>,
    /// 无向边：(原子下标, 原子下标, 键级 1/2/3)。
    pub bonds: Vec<(usize, usize, u8)>,
}

impl Molecule {
    pub fn charge(&self) -> i64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// 元素/同位素计数（未加权）。
    pub fn composition(&self) -> BTreeMap<String, i64> {
        let mut out = BTreeMap::new();
        for a in &self.atoms {
            *out.entry(a.composition_key()).or_insert(0) += 1;
        }
        out
    }

    /// 原子语义签名：(元素, 同位素, 电荷)。
    pub fn signature(&self) -> Vec<(String, i64, i64)> {
        self.atoms
            .iter()
            .map(|a| (a.element().to_string(), a.isotope, a.charge))
            .collect()
    }

    /// 邻接表（带键级颜色）。
    fn adjacency(&self) -> Vec<Vec<(usize, u8)>> {
        let mut adj = vec![Vec::new(); self.atoms.len()];
        for &(u, v, o) in &self.bonds {
            adj[u].push((v, o));
            adj[v].push((u, o));
        }
        adj
    }
}

/// 枚举分子图的全部自同构（保持元素/同位素/电荷与带色边）。
///
/// 对本项目的小分子使用朴素回溯，并对结果数设上限以防指数膨胀。
/// 对称导致的多个编号映射由此识别，业务层据此把原始映射归并为等价类，
/// 绝不把对称重排当作“按原子编号唯一”的不同映射。
pub fn automorphisms(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    if n == 0 {
        return vec![vec![]];
    }
    let adj = mol.adjacency();
    let sig = mol.signature();
    let mut group: Vec<Vec<usize>> = Vec::new();
    let mut perm = vec![usize::MAX; n];
    let mut used = vec![false; n];

    // 初始颜色 = 语义签名；回溯中保持同签名匹配。
    fn rec(
        i: usize,
        n: usize,
        perm: &mut [usize],
        used: &mut [bool],
        sig: &[(String, i64, i64)],
        adj: &[Vec<(usize, u8)>],
        group: &mut Vec<Vec<usize>>,
    ) {
        if group.len() >= 1024 {
            return;
        }
        if i == n {
            group.push(perm.to_vec());
            return;
        }
        for j in 0..n {
            if used[j] || sig[j] != sig[i] {
                continue;
            }
            // 边一致性：i 与已分配的 k 之间的边必须等于 j 与 perm[k] 之间的边。
            let mut ok = true;
            for (k, &pk) in perm.iter().enumerate().take(i) {
                let left = adj[i].iter().find(|(t, _)| *t == k).map(|(_, o)| *o);
                let right = adj[j].iter().find(|(t, _)| *t == pk).map(|(_, o)| *o);
                if left != right {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }
            used[j] = true;
            perm[i] = j;
            rec(i + 1, n, perm, used, sig, adj, group);
            used[j] = false;
            perm[i] = usize::MAX;
        }
    }

    rec(0, n, &mut perm, &mut used, &sig, &adj, &mut group);
    group
}

/// 解析 `n` 或 `n/d` 形式的非负有理数。
pub fn parse_ratio(text: &str) -> Result<Rat, String> {
    let text = text.trim();
    if let Some((a, b)) = text.split_once('/') {
        let n: i64 = a
            .trim()
            .parse()
            .map_err(|_| format!("非法有理数分子: {text}"))?;
        let d: i64 = b
            .trim()
            .parse()
            .map_err(|_| format!("非法有理数分母: {text}"))?;
        if d <= 0 {
            return Err(format!("分母必须为正: {text}"));
        }
        Ok(Ratio::new(n, d))
    } else {
        let n: i64 = text.parse().map_err(|_| format!("非法有理数: {text}"))?;
        Ok(Ratio::from_integer(n))
    }
}

/// 最小非负整数解（用于分数计量的整数展开）。
pub fn lcm(a: i64, b: i64) -> i64 {
    use num_integer::Integer;
    a.lcm(&b)
}

/// 内置小分子库（修复候选与网络池共用）。
pub fn library() -> &'static BTreeMap<String, Molecule> {
    static LIB: OnceLock<BTreeMap<String, Molecule>> = OnceLock::new();
    LIB.get_or_init(|| {
        let mut m = BTreeMap::new();
        for mol in builtin_molecules() {
            m.insert(mol.id.clone(), mol);
        }
        m
    })
}

fn atom(el: &str, iso: i64, ch: i64) -> Atom {
    let iso = match (el, iso) {
        ("D", 0) => 2,
        ("T", 0) => 3,
        _ => iso,
    };
    Atom {
        element: el.to_string(),
        isotope: iso,
        charge: ch,
    }
}

/// 用简写构建分子：`(id, 原子列表[(元素,同位素,电荷)], 边[(u,v,键级)])`。
pub fn builtin_molecules() -> Vec<Molecule> {
    let mk = |id: &str, atoms: Vec<(&str, i64, i64)>, bonds: Vec<(usize, usize, u8)>| -> Molecule {
        Molecule {
            id: id.to_string(),
            atoms: atoms.into_iter().map(|(e, i, c)| atom(e, i, c)).collect(),
            bonds,
        }
    };
    vec![
        mk(
            "H2O",
            vec![("O", 0, 0), ("H", 0, 0), ("H", 0, 0)],
            vec![(0, 1, 1), (0, 2, 1)],
        ),
        mk(
            "HDO",
            vec![("O", 0, 0), ("H", 0, 0), ("D", 0, 0)],
            vec![(0, 1, 1), (0, 2, 1)],
        ),
        mk("D", vec![("H", 2, 0)], vec![]),
        mk("H", vec![("H", 0, 0)], vec![]),
        mk("OH-", vec![("O", 0, -1), ("H", 0, 0)], vec![(0, 1, 1)]),
        mk("H+", vec![("H", 0, 1)], vec![]),
        mk("H2", vec![("H", 0, 0), ("H", 0, 0)], vec![(0, 1, 1)]),
        mk("O2", vec![("O", 0, 0), ("O", 0, 0)], vec![(0, 1, 2)]),
        mk("N2", vec![("N", 0, 0), ("N", 0, 0)], vec![(0, 1, 3)]),
        mk("NO", vec![("N", 0, 0), ("O", 0, 0)], vec![(0, 1, 2)]),
        mk(
            "N2O",
            vec![("N", 0, 0), ("N", 0, 0), ("O", 0, 0)],
            vec![(0, 1, 2), (1, 2, 2)],
        ),
        mk(
            "NO2",
            vec![("N", 0, 0), ("O", 0, 0), ("O", 0, 0)],
            vec![(0, 1, 2), (0, 2, 2)],
        ),
        mk(
            "N2O4",
            vec![
                ("N", 0, 0),
                ("O", 0, 0),
                ("O", 0, 0),
                ("N", 0, 0),
                ("O", 0, 0),
                ("O", 0, 0),
            ],
            vec![(0, 1, 2), (0, 2, 2), (0, 3, 1), (3, 4, 2), (3, 5, 2)],
        ),
        mk("CO", vec![("C", 0, 0), ("O", 0, 0)], vec![(0, 1, 3)]),
        mk(
            "CO2",
            vec![("C", 0, 0), ("O", 0, 0), ("O", 0, 0)],
            vec![(0, 1, 2), (0, 2, 2)],
        ),
        mk("Fe2+", vec![("Fe", 0, 2)], vec![]),
        mk("Fe3+", vec![("Fe", 0, 3)], vec![]),
        mk("e-", vec![("e", 0, -1)], vec![]),
        // 网络池中的抽象物种（无组成），只用于分子级流量平衡。
        mk("A", vec![], vec![]),
        mk("X", vec![], vec![]),
        mk("B", vec![], vec![]),
    ]
}
