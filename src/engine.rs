//! 审查引擎：有理数守恒、原子映射等价类、修复候选枚举、网络中间体流量。
//!
//! 全部计量运算基于 [`Rat`]（精确有理），不存在浮点容差。

use crate::chem::{automorphisms, rat_string, Molecule, Rat};
use crate::model::{Dataset, MapEdge, Op, Reaction, Side};
use num_traits::Zero;
/// 一条平铺映射边（反应物平铺下标, 产物平铺下标）。
type FlatEdge = (usize, usize);
/// 一个分子项的自同构群（每元素为一个排列）。
type PermGroup = Vec<Vec<usize>>;

use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

/// 单个参与项的运行时视图。
#[derive(Clone, Debug, Serialize)]
pub struct TermView {
    pub mol: String,
    pub coef: String,
    pub builtin: bool,
}

/// 一条守恒分量（元素@同位素 或 charge）。
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Component {
    pub key: String,
    /// 反应物 - 产物 的有理净值。0 表示该分量守恒。
    pub delta: String,
    pub balanced: bool,
}

/// 原始原子映射的一条边（视图）。
#[derive(Clone, Debug, Serialize)]
pub struct MapEdgeView {
    pub from: String,
    pub to: String,
}

/// 对称等价类：多个“按原子编号不同但分子对称下等价”的映射归并为一个。
#[derive(Clone, Debug, Serialize)]
pub struct MapClass {
    pub class_id: String,
    /// 代表映射的边。
    pub edges: Vec<MapEdgeView>,
    /// 等价类包含的原始编号映射数量（对称重排数）。
    pub multiplicity: usize,
    /// 该等价类下身份是否守恒（语义单射 + 全覆盖）。
    pub identity_ok: bool,
    /// 不满足单射/覆盖时的原因。
    pub issues: Vec<String>,
}

/// 一个修复动作。
#[derive(Clone, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct RepairAction {
    pub side: String,
    pub mol: String,
    pub coef: String,
}

/// 一个有界候选修复。
#[derive(Clone, Debug, Serialize)]
pub struct Candidate {
    pub actions: Vec<RepairAction>,
    /// 证据代价（越小越优；纯整数、只用内置小分子者更便宜）。
    pub cost: i64,
    /// 修复前赤字（未守恒向量）。
    pub deficit_before: Vec<Component>,
    /// 修复后是否全分量守恒（应为 true）。
    pub balances: bool,
    /// 证据说明。
    pub evidence: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ReactionView {
    pub id: String,
    pub reactants: Vec<TermView>,
    pub products: Vec<TermView>,
    pub conservation: Vec<Component>,
    pub element_balanced: bool,
    pub charge_balanced: bool,
    pub map_classes: Vec<MapClass>,
    pub locked_class: Option<String>,
    pub candidates: Vec<Candidate>,
    pub direction_unknown: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct FlowSolution {
    /// 方向分支，如 `R2=forward,R3=reverse`。
    pub branch: String,
    /// 各反应净流量（forward 约定下的有符号有理数）。
    pub fluxes: Vec<(String, String)>,
    /// 中间体流量（仅多步参与的物种），值为净稳态通量 0 时标注停留。
    pub intermediates: Vec<(String, String)>,
}

#[derive(Clone, Debug, Serialize)]
pub struct NetworkView {
    pub id: String,
    pub steps: Vec<String>,
    pub solutions: Vec<FlowSolution>,
}

#[derive(Clone, Debug, Serialize)]
pub struct AtomView {
    pub index: usize,
    pub label: String,
    pub isotope: i64,
    pub charge: i64,
}

#[derive(Clone, Debug, Serialize)]
pub struct BondView {
    pub u: usize,
    pub v: usize,
    pub order: u8,
}

#[derive(Clone, Debug, Serialize)]
pub struct MoleculeView {
    pub id: String,
    pub charge: i64,
    pub atoms: Vec<AtomView>,
    pub bonds: Vec<BondView>,
    pub builtin: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct Snapshot {
    pub title: String,
    pub molecules: Vec<MoleculeView>,
    pub reactions: Vec<ReactionView>,
    pub networks: Vec<NetworkView>,
    pub log: Vec<LogEntry>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LogEntry {
    pub seq: i64,
    pub at: String,
    pub op_json: String,
    pub summary: String,
}

/// 反应运行时状态（由操作重放得到）。
#[derive(Clone, Debug, Default)]
pub struct ReactionState {
    pub extras: Vec<(Side, crate::model::Term)>,
    pub locked_class: Option<String>,
    pub direction_unknown: Option<bool>,
}

pub struct Engine {
    pub dataset: Dataset,
    pub states: BTreeMap<String, ReactionState>,
    pub builtin_ids: BTreeSet<String>,
}

impl Engine {
    pub fn new(dataset: Dataset) -> Self {
        let builtin_ids = crate::chem::builtin_molecules()
            .into_iter()
            .map(|m| m.id)
            .collect();
        let mut states: BTreeMap<String, ReactionState> = BTreeMap::new();
        for r in &dataset.reactions {
            states.entry(r.id.clone()).or_default();
            if r.direction_unknown {
                states.get_mut(&r.id).unwrap().direction_unknown = Some(true);
            }
        }
        Engine {
            dataset,
            states,
            builtin_ids,
        }
    }

    pub fn apply(&mut self, op: &Op) -> Result<(), String> {
        match op {
            Op::LockMap { reaction, class_id } => {
                let exists = self
                    .map_classes(reaction)?
                    .iter()
                    .any(|c| &c.class_id == class_id);
                if !exists {
                    return Err(format!("反应 {reaction} 不存在映射等价类 {class_id}"));
                }
                self.state_mut(reaction)?.locked_class = Some(class_id.clone());
            }
            Op::UnlockMap { reaction } => {
                self.state_mut(reaction)?.locked_class = None;
            }
            Op::AddParticipant {
                reaction,
                side,
                mol,
                coef,
            } => {
                if !self.dataset.molecules.contains_key(mol) {
                    return Err(format!("未知分子 {mol}"));
                }
                let coef = crate::chem::parse_ratio(coef)?;
                if coef <= Rat::from_integer(0) {
                    return Err("系数必须为正".to_string());
                }
                self.state_mut(reaction)?.extras.push((
                    *side,
                    crate::model::Term {
                        mol: mol.clone(),
                        coef,
                    },
                ));
            }
            Op::SetDirectionUnknown { reaction, unknown } => {
                self.state_mut(reaction)?.direction_unknown = Some(*unknown);
            }
            Op::ResetImported { .. } => {}
        }
        Ok(())
    }

    fn state_mut(&mut self, reaction: &str) -> Result<&mut ReactionState, String> {
        if !self.dataset.reactions.iter().any(|r| r.id == reaction) {
            return Err(format!("未知反应 {reaction}"));
        }
        Ok(self.states.entry(reaction.to_string()).or_default())
    }

    fn terms(&self, rxn_id: &str, side: Side) -> Vec<&crate::model::Term> {
        let rxn = self.rxn(rxn_id);
        let mut out: Vec<&crate::model::Term> = rxn.terms(side).iter().collect();
        if let Some(st) = self.states.get(rxn_id) {
            for (s, t) in &st.extras {
                if *s == side {
                    out.push(t);
                }
            }
        }
        out
    }

    fn rxn(&self, id: &str) -> &Reaction {
        self.dataset
            .reactions
            .iter()
            .find(|r| r.id == id)
            .expect("reaction exists")
    }

    fn is_unknown(&self, id: &str) -> bool {
        match self.states.get(id).and_then(|s| s.direction_unknown) {
            Some(v) => v,
            None => self.rxn(id).direction_unknown,
        }
    }
}

/// 守恒向量：反应物加权组成 - 产物加权组成，外加 `charge`。
fn conservation_vector(
    mols: &BTreeMap<String, Molecule>,
    reactants: &[&crate::model::Term],
    products: &[&crate::model::Term],
) -> Vec<Component> {
    let mut acc: BTreeMap<String, Rat> = BTreeMap::new();
    let mut charge = Rat::from_integer(0);
    let add = |terms: &[&crate::model::Term],
               sign: Rat,
               acc: &mut BTreeMap<String, Rat>,
               charge: &mut Rat|
     -> () {
        for t in terms {
            let mol = &mols[&t.mol];
            for a in &mol.atoms {
                if a.element() != "e" {
                    *acc.entry(a.composition_key()).or_insert_with(Rat::zero) += sign * t.coef;
                }
            }
            *charge += sign * t.coef * Rat::from_integer(mol.charge());
        }
    };
    add(reactants, Rat::from_integer(1), &mut acc, &mut charge);
    add(products, Rat::from_integer(-1), &mut acc, &mut charge);
    acc.insert("charge".to_string(), charge);
    acc.into_iter()
        .map(|(key, delta)| Component {
            key,
            balanced: delta == Rat::from_integer(0),
            delta: rat_string(delta),
        })
        .collect()
}

/// fnv1a-64，用于由代表边集合生成稳定的 class id。
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

impl Engine {
    /// 计算反应的对称映射等价类。
    ///
    /// 平铺原子下标 = 该侧各分子项原子顺序拼接。自同构群按分子项分组生成，
    /// 原始映射在 G_R × G_P 群作用下的轨道就是一个等价类。
    #[allow(clippy::type_complexity)]
    pub fn map_classes(&self, rxn_id: &str) -> Result<Vec<MapClass>, String> {
        let rxn = self.rxn(rxn_id);
        let mols = &self.dataset.molecules;

        let build_side = |terms: &[crate::model::Term]|
            -> Result<(Vec<String>, Vec<usize>, Vec<PermGroup>), String> {
            let mut mol_ids = Vec::new();
            let mut offsets: Vec<usize> = vec![0usize];
            let mut groups: Vec<Vec<Vec<usize>>> = Vec::new();
            let mut total = 0usize;
            for t in terms {
                let mol = &mols[&t.mol];
                if !rxn.maps.is_empty() && t.coef != Rat::from_integer(1) {
                    return Err(format!(
                        "reaction {}: 含原子映射的反应要求单位计量（{} 系数 {}）",
                        rxn.id,
                        t.mol,
                        rat_string(t.coef)
                    ));
                }
                mol_ids.push(t.mol.clone());
                total += mol.atoms.len();
                offsets.push(total);
                groups.push(automorphisms(mol));
            }
            Ok((mol_ids, offsets, groups))
        };

        let (r_ids, r_off, r_groups) = build_side(&rxn.reactants)?;
        let (p_ids, p_off, p_groups) = build_side(&rxn.products)?;

        // 每个分子项生成“局部排列+平移”，再做笛卡尔积得到整侧排列。
        let side_perms =
            |ids_len: usize, offsets: &[usize], groups: &[Vec<Vec<usize>>]| -> Vec<Vec<usize>> {
                let mut acc: Vec<Vec<usize>> = vec![vec![]];
                for k in 0..ids_len {
                    let base = offsets[k];
                    let mut next = Vec::new();
                    for prefix in &acc {
                        for local in &groups[k] {
                            let mut q = prefix.clone();
                            q.extend(local.iter().map(|x| x + base));
                            next.push(q);
                        }
                    }
                    acc = next;
                    if acc.len() > 4096 {
                        break;
                    }
                }
                acc
            };

        let r_perms = side_perms(r_ids.len(), &r_off, &r_groups);
        let p_perms = side_perms(p_ids.len(), &p_off, &p_groups);

        let r_total = *r_off.last().unwrap();
        let p_total = *p_off.last().unwrap();

        let locate = |mol: &str, atom: usize, ids: &[String], off: &[usize]| -> Option<usize> {
            ids.iter().position(|m| m == mol).map(|k| off[k] + atom)
        };

        let atom_at = |flat: usize, ids: &[String], off: &[usize]| -> (String, usize) {
            let k = (0..ids.len()).find(|&k| flat < off[k + 1]).unwrap();
            (ids[k].clone(), flat - off[k])
        };

        // 校验一条原始映射并转成平铺边集合。
        let flatten = |edges: &[MapEdge]| -> Result<BTreeSet<(usize, usize)>, Vec<String>> {
            let mut out = BTreeSet::new();
            let mut issues = Vec::new();
            for e in edges {
                let rf = match locate(&e.from_mol, e.from_atom, &r_ids, &r_off) {
                    Some(x) => x,
                    None => {
                        issues.push(format!("边引用了非反应物 {}:{}", e.from_mol, e.from_atom));
                        continue;
                    }
                };
                let pf = match locate(&e.to_mol, e.to_atom, &p_ids, &p_off) {
                    Some(x) => x,
                    None => {
                        issues.push(format!("边引用了非产物 {}:{}", e.to_mol, e.to_atom));
                        continue;
                    }
                };
                // 类型/同位素/电荷语义单射：解析器已校验，这里对用户状态再防御一次。
                let rm = {
                    let (m, a) = atom_at(rf, &r_ids, &r_off);
                    &mols[&m].atoms[a]
                };
                let pmol = {
                    let (m, a) = atom_at(pf, &p_ids, &p_off);
                    &mols[&m].atoms[a]
                };
                if rm.element() != pmol.element() {
                    issues.push(format!(
                        "{}:{} -> {}:{} 的元素身份不一致",
                        e.from_mol, e.from_atom, e.to_mol, e.to_atom
                    ));
                }
                out.insert((rf, pf));
            }
            // 单射：每个端点至多出现一次。
            let mut froms = BTreeSet::new();
            let mut tos = BTreeSet::new();
            for (f, t) in &out {
                if !froms.insert(*f) {
                    issues.push(format!("反应物原子 {f} 被映射到多个产物（非单射）"));
                }
                if !tos.insert(*t) {
                    issues.push(format!("产物原子 {t} 被多个反应物映射（非单射）"));
                }
            }
            // 全覆盖：单位计量映射必须覆盖每个原子（原子身份守恒）。
            if out.len() != r_total {
                issues.push(format!(
                    "覆盖不全：反应物 {r_total} 个原子，仅映射 {} 个",
                    out.len()
                ));
            }
            if out.len() != p_total {
                issues.push(format!(
                    "覆盖不全：产物 {p_total} 个原子，仅映射 {} 个",
                    out.len()
                ));
            }
            if issues.is_empty() {
                Ok(out)
            } else {
                Err(issues)
            }
        };

        // 用规范代表边集合归并等价类，保留 multiplicity。
        // 等价类表：规范边集合 -> (出现的原始映射数, 问题列表, 群轨道)。
        let mut classes: BTreeMap<
            BTreeSet<FlatEdge>,
            (usize, Vec<String>, BTreeSet<BTreeSet<FlatEdge>>),
        > = BTreeMap::new();

        for raw in &rxn.maps {
            let (canon_set, issues) = match flatten(raw) {
                Ok(set) => {
                    // 在群作用下求规范像。
                    let mut best: Option<BTreeSet<(usize, usize)>> = None;
                    for rp in &r_perms {
                        for pp in &p_perms {
                            let img: BTreeSet<(usize, usize)> =
                                set.iter().map(|(f, t)| (rp[*f], pp[*t])).collect();
                            best = Some(match best {
                                Some(b) => {
                                    if img < b {
                                        img
                                    } else {
                                        b
                                    }
                                }
                                None => img,
                            });
                        }
                    }
                    (best.unwrap(), Vec::new())
                }
                Err(issues) => {
                    // 非法映射单独成类，不参与对称归并。
                    let rawset: BTreeSet<(usize, usize)> = raw
                        .iter()
                        .filter_map(|e| {
                            let rf = locate(&e.from_mol, e.from_atom, &r_ids, &r_off)?;
                            let pf = locate(&e.to_mol, e.to_atom, &p_ids, &p_off)?;
                            Some((rf, pf))
                        })
                        .collect();
                    (rawset, issues)
                }
            };
            // multiplicity = 该轨道在当前自同构群下的大小。
            let mut orbit: BTreeSet<BTreeSet<FlatEdge>> = BTreeSet::new();
            if issues.is_empty() {
                for rp in &r_perms {
                    for pp in &p_perms {
                        let img: BTreeSet<(usize, usize)> =
                            canon_set.iter().map(|(f, t)| (rp[*f], pp[*t])).collect();
                        orbit.insert(img);
                    }
                }
            }
            let entry =
                classes
                    .entry(canon_set.clone())
                    .or_insert((0, issues.clone(), BTreeSet::new()));
            entry.0 += 1;
            if entry.1.is_empty() {
                entry.1 = issues;
            }
            entry.2.extend(orbit);
        }

        let mut out = Vec::new();
        for (canon, (_raw_count, issues, orbit)) in classes {
            let mut edge_views = Vec::new();
            for (f, t) in &canon {
                let (fm, fa) = atom_at(*f, &r_ids, &r_off);
                let (tm, ta) = atom_at(*t, &p_ids, &p_off);
                edge_views.push(MapEdgeView {
                    from: format!("{fm}:{fa}"),
                    to: format!("{tm}:{ta}"),
                });
            }
            let id_src = edge_views
                .iter()
                .map(|e| format!("{}>{}", e.from, e.to))
                .collect::<Vec<_>>()
                .join("|");
            let class_id = format!("{:016x}", fnv1a64(id_src.as_bytes()));
            out.push(MapClass {
                class_id,
                edges: edge_views,
                multiplicity: if orbit.is_empty() { 1 } else { orbit.len() },
                identity_ok: issues.is_empty(),
                issues,
            });
        }
        out.sort_by(|a, b| a.class_id.cmp(&b.class_id));
        Ok(out)
    }
}

fn parse_component_delta(c: &Component) -> Rat {
    crate::chem::parse_ratio(&c.delta).expect("内部 delta 始终可解析")
}

impl Engine {
    /// 有界枚举修复候选：向反应物或产物侧补入正系数内置小分子（1 或 2 个动作）。
    ///
    /// 系数网格取 {1/2, 1, 3/2, 2}，分母有界 => 候选空间有界。
    /// 按动作集合去重，保留每个候选的修复动作与证据代价。
    fn candidates_for(&self, _rxn_id: &str, deficit: &[Component]) -> Vec<Candidate> {
        let mols = &self.dataset.molecules;
        let deficit_map: BTreeMap<String, Rat> = deficit
            .iter()
            .map(|c| (c.key.clone(), parse_component_delta(c)))
            .collect();

        // 可补的小分子：内置库中含原子的分子（含电子伪分子）。
        let pool: Vec<&Molecule> = self
            .builtin_ids
            .iter()
            .filter_map(|id| mols.get(id))
            .filter(|m| !m.atoms.is_empty())
            .collect();

        let grid = [
            Rat::new(1, 2),
            Rat::from_integer(1),
            Rat::new(3, 2),
            Rat::from_integer(2),
        ];
        let sides = ["reactant", "product"];

        // 单个动作的有符号贡献：补反应物为 +，补产物为 -。
        // delta = 反应物 - 产物；补入后要求 delta + signed_add = 0。
        let action_vector = |side: &str, id: &str, coef: Rat| -> BTreeMap<String, Rat> {
            let sign = if side == "reactant" {
                Rat::from_integer(1)
            } else {
                Rat::from_integer(-1)
            };
            let mol = &mols[id];
            let mut v = BTreeMap::new();
            for a in &mol.atoms {
                if a.element() != "e" {
                    *v.entry(a.composition_key()).or_insert_with(Rat::zero) += sign * coef;
                }
            }
            *v.entry("charge".to_string()).or_insert_with(Rat::zero) +=
                sign * coef * Rat::from_integer(mol.charge());
            v
        };

        let mut found: BTreeMap<Vec<RepairAction>, Vec<String>> = BTreeMap::new();
        let mut consider = |acts: Vec<(&str, String, Rat)>| {
            let mut add: BTreeMap<String, Rat> = BTreeMap::new();
            for (side, id, coef) in &acts {
                for (k, val) in action_vector(side, id, *coef) {
                    *add.entry(k).or_insert_with(Rat::zero) += val;
                }
            }
            // 不得引入赤字之外、且非零的新分量。
            for (k, v) in &add {
                if *v != Rat::from_integer(0) && !deficit_map.contains_key(k) {
                    return;
                }
            }
            for (k, delta) in &deficit_map {
                let got = add.get(k).copied().unwrap_or_else(Rat::zero);
                if got != -*delta {
                    return;
                }
            }
            let mut repair: Vec<RepairAction> = acts
                .iter()
                .map(|(side, id, coef)| RepairAction {
                    side: (*side).to_string(),
                    mol: id.clone(),
                    coef: rat_string(*coef),
                })
                .collect();
            repair.sort();
            let evidence = acts
                .iter()
                .map(|(side, id, coef)| {
                    let side_cn = if *side == "reactant" {
                        "反应物"
                    } else {
                        "产物"
                    };
                    format!(
                        "向{side_cn}补 {} ×{}（精确有理核对）",
                        id,
                        rat_string(*coef)
                    )
                })
                .collect();
            found.insert(repair, evidence);
        };

        // 1 个动作
        for m in &pool {
            for side in sides {
                for &c in &grid {
                    consider(vec![(side, m.id.clone(), c)]);
                }
            }
        }
        // 2 个动作（允许同分子不同侧；按“侧+分子”去重）
        let mut one: Vec<(&str, &Molecule)> = Vec::new();
        for m in &pool {
            for side in sides {
                one.push((side, *m));
            }
        }
        for i in 0..one.len() {
            for j in i..one.len() {
                let (s1, m1) = &one[i];
                let (s2, m2) = &one[j];
                for &c1 in &grid {
                    for &c2 in &grid {
                        consider(vec![(*s1, m1.id.clone(), c1), (*s2, m2.id.clone(), c2)]);
                    }
                }
            }
        }

        let mut cands: Vec<Candidate> = found
            .into_iter()
            .map(|(actions, evidence)| {
                let mut frac_penalty = 0i64;
                let mut num_penalty = 0i64;
                let mut sides_used = BTreeSet::new();
                for a in &actions {
                    sides_used.insert(a.side.clone());
                    let r = crate::chem::parse_ratio(&a.coef).unwrap();
                    if *r.denom() != 1 {
                        frac_penalty += 1;
                    }
                    num_penalty += if *r.denom() == 1 {
                        *r.numer()
                    } else {
                        *r.numer() + 2 * *r.denom()
                    };
                }
                // 证据代价：跨两侧的“对消式”补入证据更弱，给予大惩罚。
                let cross_side_penalty = if sides_used.len() > 1 { 1000 } else { 0 };
                Candidate {
                    cost: cross_side_penalty
                        + 100 * frac_penalty
                        + 10 * (actions.len() as i64)
                        + num_penalty,
                    actions,
                    deficit_before: deficit.to_vec(),
                    balances: true,
                    evidence,
                }
            })
            .collect();
        cands.sort_by(|a, b| a.cost.cmp(&b.cost).then_with(|| a.actions.cmp(&b.actions)));
        cands
    }

    /// 网络中间体流量：对每条“方向未定”的 step 枚举方向分支，
    /// 在有界有理数网格上求稳态（每个中间体净通量为 0）。
    /// 不同方向分支即使流量向量相同也分别保留。
    fn network_solutions(&self, net: &crate::model::Network) -> Vec<FlowSolution> {
        let steps = &net.steps;

        // 每个 step 的分子级净化学计量：reactants - products（forward 约定）。
        let mut stoich: Vec<BTreeMap<String, Rat>> = Vec::new();
        for (rid, _) in steps {
            let mut v = BTreeMap::new();
            // 净生成向量：产物为正、反应物为负；稳态 sum_j S_j * f_j = 0。
            for t in &self.rxn(rid).reactants {
                *v.entry(t.mol.clone()).or_insert_with(Rat::zero) -= t.coef;
            }
            for t in &self.rxn(rid).products {
                *v.entry(t.mol.clone()).or_insert_with(Rat::zero) += t.coef;
            }
            stoich.push(v);
        }

        // 中间体：在不同 step 中既作反应物又作产物的物种（含无组成池物种）。
        let mut as_react: BTreeSet<String> = BTreeSet::new();
        let mut as_prod: BTreeSet<String> = BTreeSet::new();
        for (rid, _) in steps {
            for t in &self.rxn(rid).reactants {
                as_react.insert(t.mol.clone());
            }
            for t in &self.rxn(rid).products {
                as_prod.insert(t.mol.clone());
            }
        }
        let intermediates: Vec<String> = as_react.intersection(&as_prod).cloned().collect();

        // 固定符号 / 未知标记。
        let unknown: Vec<bool> = steps
            .iter()
            .map(|(rid, marked)| *marked || self.is_unknown(rid))
            .collect();

        // 有界有理数网格：净流取 {0, 1/2, 1}（reverse 分支取负）。
        let grid_nonneg = [Rat::from_integer(0), Rat::new(1, 2), Rat::from_integer(1)];

        // 枚举方向分支：未知 step 取 forward(净流>=0)/reverse(净流<=0)。
        let uk_positions: Vec<usize> = unknown
            .iter()
            .enumerate()
            .filter(|(_, u)| **u)
            .map(|(i, _)| i)
            .collect();

        let mut solutions = Vec::new();
        let branch_count = 1usize << uk_positions.len();
        for b in 0..branch_count {
            let mut forward_for = vec![true; steps.len()];
            let mut labels = Vec::new();
            for (k, &idx) in uk_positions.iter().enumerate() {
                let fwd = (b >> k) & 1 == 0;
                forward_for[idx] = fwd;
                labels.push(format!(
                    "{}={}",
                    steps[idx].0,
                    if fwd { "forward" } else { "reverse" }
                ));
            }
            if labels.is_empty() {
                labels.push("(全部方向已定)".to_string());
            }
            let branch = labels.join(",");

            // 在网格上求净流 f，使对每个中间体 sum_j S[j][m]*f_j = 0。
            let mut f = vec![Rat::from_integer(0); steps.len()];
            self.enumerate_flows(
                0,
                &stoich,
                &intermediates,
                &unknown,
                &forward_for,
                &grid_nonneg,
                &mut f,
                &mut |f: &[Rat]| {
                    let fluxes = steps
                        .iter()
                        .map(|(rid, _)| {
                            (
                                rid.clone(),
                                rat_string(f[steps.iter().position(|(r, _)| r == rid).unwrap()]),
                            )
                        })
                        .collect();
                    let inter = intermediates
                        .iter()
                        .map(|m| (m.clone(), "0".to_string()))
                        .collect();
                    solutions.push(FlowSolution {
                        branch: branch.clone(),
                        fluxes,
                        intermediates: inter,
                    });
                },
            );
        }

        // 稳定排序去重：同分支同向量只留一份。
        let mut seen = BTreeSet::new();
        solutions.retain(|s| {
            let key = (
                s.branch.clone(),
                s.fluxes
                    .iter()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect::<Vec<_>>(),
            );
            seen.insert(key)
        });
        solutions
    }

    #[allow(clippy::too_many_arguments)]
    fn enumerate_flows(
        &self,
        idx: usize,
        stoich: &[BTreeMap<String, Rat>],
        intermediates: &[String],
        unknown: &[bool],
        forward_for: &[bool],
        grid_nonneg: &[Rat; 3],
        f: &mut [Rat],
        emit: &mut dyn FnMut(&[Rat]),
    ) {
        if idx == stoich.len() {
            // 校验所有中间体净通量为零。
            for m in intermediates {
                let mut net = Rat::from_integer(0);
                for (j, s) in stoich.iter().enumerate() {
                    net += s.get(m).copied().unwrap_or_else(Rat::zero) * f[j];
                }
                if net != Rat::from_integer(0) {
                    return;
                }
            }
            emit(f);
            return;
        }

        // 已定方向步骤按 forward 约定只取非负净流；
        // 未定步骤按当前分支取非负(forward)或非正(reverse)。
        let candidates: Vec<Rat> = if unknown[idx] {
            grid_nonneg
                .iter()
                .copied()
                .map(|v| if forward_for[idx] { v } else { -v })
                .collect()
        } else {
            grid_nonneg.to_vec()
        };

        for val in candidates {
            f[idx] = val;
            self.enumerate_flows(
                idx + 1,
                stoich,
                intermediates,
                unknown,
                forward_for,
                grid_nonneg,
                f,
                emit,
            );
        }
    }
}

impl Engine {
    /// 组装给页面的快照视图。
    pub fn snapshot(&self, log: Vec<LogEntry>) -> Snapshot {
        let mut molecules: Vec<MoleculeView> = self
            .dataset
            .molecules
            .values()
            .map(|m| MoleculeView {
                id: m.id.clone(),
                charge: m.charge(),
                atoms: m
                    .atoms
                    .iter()
                    .enumerate()
                    .map(|(i, a)| AtomView {
                        index: i,
                        label: format!("{}[{}]", a.element(), i),
                        isotope: a.isotope,
                        charge: a.charge,
                    })
                    .collect(),
                bonds: m
                    .bonds
                    .iter()
                    .map(|(u, v, o)| BondView {
                        u: *u,
                        v: *v,
                        order: *o,
                    })
                    .collect(),
                builtin: self.builtin_ids.contains(&m.id),
            })
            .collect();
        molecules.sort_by(|a, b| a.id.cmp(&b.id));

        let mut reactions = Vec::new();
        for rxn in &self.dataset.reactions {
            let reactants = self
                .terms(&rxn.id, Side::Reactant)
                .into_iter()
                .map(|t| TermView {
                    mol: t.mol.clone(),
                    coef: rat_string(t.coef),
                    builtin: self.builtin_ids.contains(&t.mol),
                })
                .collect();
            let products = self
                .terms(&rxn.id, Side::Product)
                .into_iter()
                .map(|t| TermView {
                    mol: t.mol.clone(),
                    coef: rat_string(t.coef),
                    builtin: self.builtin_ids.contains(&t.mol),
                })
                .collect();

            let cons = conservation_vector(
                &self.dataset.molecules,
                &self.terms(&rxn.id, Side::Reactant),
                &self.terms(&rxn.id, Side::Product),
            );
            let element_balanced = cons
                .iter()
                .filter(|c| c.key != "charge" && c.key != "e")
                .all(|c| c.balanced);
            let charge_balanced = cons
                .iter()
                .find(|c| c.key == "charge")
                .map(|c| c.balanced)
                .unwrap_or(true);

            let map_classes = self.map_classes(&rxn.id).unwrap_or_default();
            let locked_class = self
                .states
                .get(&rxn.id)
                .and_then(|s| s.locked_class.clone());

            let deficit: Vec<Component> = cons.iter().filter(|c| !c.balanced).cloned().collect();
            let candidates = if deficit.is_empty() {
                Vec::new()
            } else {
                self.candidates_for(&rxn.id, &deficit)
            };

            reactions.push(ReactionView {
                id: rxn.id.clone(),
                reactants,
                products,
                conservation: cons,
                element_balanced,
                charge_balanced,
                map_classes,
                locked_class,
                candidates,
                direction_unknown: self.is_unknown(&rxn.id),
            });
        }

        let networks = self
            .dataset
            .networks
            .iter()
            .map(|n| NetworkView {
                id: n.id.clone(),
                steps: n
                    .steps
                    .iter()
                    .map(|(r, u)| {
                        if *u || self.is_unknown(r) {
                            format!("{r} (方向未定)")
                        } else {
                            r.clone()
                        }
                    })
                    .collect(),
                solutions: self.network_solutions(n),
            })
            .collect();

        Snapshot {
            title: "反应守恒审查台".to_string(),
            molecules,
            reactions,
            networks,
            log,
        }
    }
}
