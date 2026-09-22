//! 化学核心：受限本地格式解析、精确有理守恒、原子身份映射、修复候选与网络流量。

pub mod network;
pub mod repair;

use std::collections::{BTreeMap, BTreeSet};

pub use num_rational::Ratio;
use num_traits::Zero;

/// 化学计量数一律使用精确有理数，禁止浮点容差。
pub type Rat = Ratio<i64>;

pub fn rat(n: i64) -> Rat {
    Ratio::from_integer(n)
}

/// 原子身份键：类型（元素）+ 同位素质量数 + 电荷。映射必须保持该键。
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct AtomKey {
    pub element: String,
    pub isotope: Option<u32>,
    pub charge: i32,
}

impl AtomKey {
    pub fn new(element: &str, isotope: Option<u32>, charge: i32) -> Self {
        Self {
            element: element.to_string(),
            isotope,
            charge,
        }
    }

    /// 例如 `[13C]`、`[D]`、`[18O-]`、`[H+]`。
    pub fn tag(&self) -> String {
        let mut s = String::new();
        if let Some(iso) = self.isotope {
            s.push('[');
            s.push_str(&iso.to_string());
            s.push_str(&self.element);
        } else {
            s.push_str(&self.element);
        }
        if self.charge != 0 {
            s.push(if self.charge > 0 { '+' } else { '-' });
            if self.charge.abs() > 1 {
                s.push_str(&self.charge.abs().to_string());
            }
        }
        s
    }
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub key: AtomKey,
}

#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

impl BondOrder {
    pub fn label(&self) -> &'static str {
        match self {
            BondOrder::Single => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Bond {
    pub a: usize,
    pub b: usize,
    pub order: BondOrder,
}

#[derive(Clone, Debug)]
pub struct Molecule {
    pub id: String,
    pub name: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl Molecule {
    /// 原子数不超过 3 的分子视为可用于修复枚举的“小分子”。
    pub fn is_small(&self) -> bool {
        self.atoms.len() <= 3
    }

    /// 未守恒向量中的分子式计数：键、电荷、同位素 → 数量。
    pub fn composition(&self) -> BTreeMap<AtomKey, i64> {
        let mut m = BTreeMap::new();
        for a in &self.atoms {
            *m.entry(a.key.clone()).or_insert(0) += 1;
        }
        m
    }
}

/// 反应一侧的一个参与项：分子 × 有理系数。
#[derive(Clone, Debug)]
pub struct Term {
    pub mol_id: String,
    pub coef: Rat,
}

#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum Side {
    Reactants,
    Products,
}

#[derive(Clone, Debug, Default)]
pub struct Record {
    pub id: String,
    pub reactants: Vec<Term>,
    pub products: Vec<Term>,
    /// 原始记录中声明（但未必合法）的原子映射边，仅作为证据展示。
    pub declared_maps: Vec<(String, String)>,
    /// 方向未定：仅影响展示与网络语义，不改变守恒残差。
    pub direction_unknown: bool,
}

impl Record {
    pub fn terms(&self, side: Side) -> &[Term] {
        match side {
            Side::Reactants => &self.reactants,
            Side::Products => &self.products,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Network {
    pub id: String,
    pub name: String,
    /// 参与反应的记录 id，顺序即网络中的列序。
    pub reaction_ids: Vec<String>,
    /// 方向未定的反应 id 集合；其余固定正向。
    pub direction_unknown: BTreeSet<String>,
    /// 外部净通量：正=净进料，负=净出料，缺失=中间体（净 0）。
    pub boundary: BTreeMap<String, Rat>,
}

/// 元素/同位素/电荷守恒残差（产物 - 反应物）。
#[derive(Clone, Debug, Default)]
pub struct Residual {
    pub atoms: BTreeMap<AtomKey, Rat>,
}

impl Residual {
    pub fn is_zero(&self) -> bool {
        self.atoms.values().all(Rat::is_zero)
    }

    pub fn entries(&self) -> Vec<(AtomKey, Rat)> {
        self.atoms
            .iter()
            .filter(|(_, v)| !v.is_zero())
            .map(|(k, v)| (k.clone(), *v))
            .collect()
    }

    pub fn total_charge(&self) -> Rat {
        self.atoms
            .iter()
            .map(|(k, v)| v * rat(i64::from(k.charge)))
            .fold(Rat::zero(), |a, b| a + b)
    }
}

/// 分子库 + 反应记录 + 网络：一次导入的完整数据口径。
#[derive(Clone, Debug, Default)]
pub struct Dataset {
    pub molecules: BTreeMap<String, Molecule>,
    pub reactions: BTreeMap<String, Record>,
    pub networks: BTreeMap<String, Network>,
}

impl Dataset {
    /// 计算一条记录的守恒残差（产物 - 反应物）。
    pub fn residual(&self, rec: &Record) -> Residual {
        let mut acc: BTreeMap<AtomKey, Rat> = BTreeMap::new();
        for term in &rec.products {
            if let Some(mol) = self.molecules.get(&term.mol_id) {
                for (k, n) in mol.composition() {
                    *acc.entry(k).or_insert_with(Rat::zero) += term.coef * rat(n);
                }
            }
        }
        for term in &rec.reactants {
            if let Some(mol) = self.molecules.get(&term.mol_id) {
                for (k, n) in mol.composition() {
                    *acc.entry(k).or_insert_with(Rat::zero) -= term.coef * rat(n);
                }
            }
        }
        Residual { atoms: acc }
    }

    /// 应用一组修复操作后得到候选记录。
    pub fn apply_ops(&self, rec: &Record, ops: &[RepairOp]) -> Record {
        let mut out = rec.clone();
        for op in ops {
            let side = match op.side {
                Side::Reactants => &mut out.reactants,
                Side::Products => &mut out.products,
            };
            if let Some(t) = side.iter_mut().find(|t| t.mol_id == op.mol_id) {
                t.coef += op.coef;
            } else {
                side.push(Term {
                    mol_id: op.mol_id.clone(),
                    coef: op.coef,
                });
            }
            side.retain(|t| t.coef > Rat::zero());
        }
        out
    }
}

/// 一次“显式添加参与小分子”的修复操作。
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct RepairOp {
    pub side: Side,
    pub mol_id: String,
    pub coef: Rat,
}

impl RepairOp {
    pub fn describe(&self) -> String {
        let side = match self.side {
            Side::Reactants => "反应物侧",
            Side::Products => "产物侧",
        };
        format!("在{}添加 {} × {}", side, self.mol_id, self.coef)
    }
}
