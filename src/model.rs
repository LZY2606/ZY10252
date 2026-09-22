//! 数据模型：反应、网络、用户操作（运行记录）。

use crate::chem::Rat;
use serde::{Deserialize, Serialize};

/// 反应一侧的一个参与项：分子 id + 有理数计量系数。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Term {
    pub mol: String,
    pub coef: Rat,
}

/// 一条原始映射边：反应物原子 -> 产物原子（原子在对应分子模板内的下标）。
/// 反应记录若出现多个独立 map 块，每块都是一份完整映射。
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct MapEdge {
    pub from_mol: String,
    pub from_atom: usize,
    pub to_mol: String,
    pub to_atom: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reaction {
    pub id: String,
    pub reactants: Vec<Term>,
    pub products: Vec<Term>,
    pub maps: Vec<Vec<MapEdge>>,
    /// 反应记录标记的方向未定（fixture 默认值，可由用户操作翻转）。
    pub direction_unknown: bool,
}

#[derive(Clone, Debug)]
pub struct Network {
    pub id: String,
    /// (反应 id, 方向是否未定)。
    pub steps: Vec<(String, bool)>,
}

/// 分子图（fixture 内联定义的分子与内置库合并）。
#[derive(Clone, Debug)]
pub struct Dataset {
    pub molecules: std::collections::BTreeMap<String, crate::chem::Molecule>,
    pub reactions: Vec<Reaction>,
    pub networks: Vec<Network>,
}

/// 用户操作类型（运行记录，可导出、清空后重放复核）。
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
pub enum Op {
    /// 锁定一组映射：每反应选择一个等价类（用代表边集合的稳定 id）。
    LockMap { reaction: String, class_id: String },
    /// 解锁某反应的映射选择。
    UnlockMap { reaction: String },
    /// 添加明确的参与小分子（修复候选应用后写入反应）。
    AddParticipant {
        reaction: String,
        side: Side,
        mol: String,
        /// 有理数字符串，如 "1/4"。
        coef: String,
    },
    /// 标记反应方向未定 / 恢复已定。
    SetDirectionUnknown { reaction: String, unknown: bool },
    /// 重置（清空数据库并重新导入 fixture）的标记操作。
    ResetImported { source: String },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Side {
    Reactant,
    Product,
}

impl Reaction {
    pub fn terms(&self, side: Side) -> &[Term] {
        match side {
            Side::Reactant => &self.reactants,
            Side::Product => &self.products,
        }
    }
}
