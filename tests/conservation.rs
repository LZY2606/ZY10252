//! 验收测试：有理守恒、分数计量、同位素、遗漏小分子修复、
//! 对称映射等价类、电荷、身份守恒、网络多方向流量解、清空重导入、HTTP 冒烟。

use conservation_bench::chem::{parse_ratio, rat_string, Rat};
use conservation_bench::engine::Engine;
use conservation_bench::model::{Op, Side};
use conservation_bench::parser;
use conservation_bench::storage::Store;
use conservation_bench::{FIXTURE, FIXTURE_NAME};

fn engine() -> Engine {
    let ds = parser::parse(FIXTURE).expect("fixture 解析成功");
    Engine::new(ds)
}

fn snap(e: &Engine) -> conservation_bench::engine::Snapshot {
    e.snapshot(vec![])
}

fn find_rxn<'a>(
    s: &'a conservation_bench::engine::Snapshot,
    id: &str,
) -> &'a conservation_bench::engine::ReactionView {
    s.reactions.iter().find(|r| r.id == id).unwrap()
}

#[test]
fn rational_arithmetic_is_exact() {
    let half = parse_ratio("1/2").unwrap();
    let sum = half + half;
    assert_eq!(rat_string(sum), "1");
    assert_eq!(rat_string(Rat::new(1, 4) + Rat::new(1, 4)), "1/2");
    // 不允许浮点容差：极小但非零的差仍为未守恒。
    let tiny = Rat::new(1, 1_000_000_000);
    assert_ne!(tiny, Rat::from_integer(0));
}

#[test]
fn r2_fractional_stoichiometry_conserves_elements() {
    let s = snap(&engine());
    let r = find_rxn(&s, "R2");
    assert!(r.element_balanced, "NO + 1/2 O2 -> NO2 应按有理数守恒");
    assert!(r.charge_balanced);
    assert!(r.conservation.iter().all(|c| c.balanced));
}

#[test]
fn r1_isotopes_and_symmetric_maps() {
    let s = snap(&engine());
    let r = find_rxn(&s, "R1");
    // 元素、同位素（H/H@1/H@2）、电荷全部守恒。
    assert!(r.element_balanced);
    assert!(r.charge_balanced);
    let keys: Vec<&str> = r.conservation.iter().map(|c| c.key.as_str()).collect();
    assert!(keys.contains(&"H"), "应出现未标记氢分量 H: {keys:?}");
    assert!(keys.contains(&"H@2"), "应出现氘分量 H@2: {keys:?}");

    // 两个对称映射必须归并为同一个等价类，且 multiplicity = 2。
    assert_eq!(r.map_classes.len(), 1, "对称映射不得按原子编号伪装成多个");
    let c = &r.map_classes[0];
    assert_eq!(c.multiplicity, 2, "水的 H 互换给出两个等价映射");
    assert!(c.identity_ok, "语义单射且全覆盖 => 身份守恒");
    assert_eq!(c.edges.len(), 4);
}

#[test]
fn r3_missing_small_molecule_candidate() {
    let s = snap(&engine());
    let r = find_rxn(&s, "R3");
    assert!(!r.element_balanced);
    let deficit: Vec<_> = r.conservation.iter().filter(|c| !c.balanced).collect();
    assert_eq!(deficit.len(), 1);
    assert_eq!(deficit[0].key, "O");
    assert_eq!(deficit[0].delta, "-1", "反应物缺一个 O（有理 Δ=-1）");

    // 有界枚举应给出补 1/2 O2 的候选，且它是排序第一（唯一）的修复。
    assert!(!r.candidates.is_empty(), "应枚举到修复候选");
    let best = &r.candidates[0];
    assert_eq!(best.actions.len(), 1);
    assert_eq!(best.actions[0].mol, "O2");
    assert_eq!(best.actions[0].coef, "1/2");
    assert!(best.balances);
    assert!(best.cost > 0, "候选带有证据代价");

    // 应用修复后守恒。
    let mut e = engine();
    e.apply(&Op::AddParticipant {
        reaction: "R3".to_string(),
        side: Side::Reactant,
        mol: "O2".to_string(),
        coef: "1/2".to_string(),
    })
    .unwrap();
    let s2 = snap(&e);
    assert!(find_rxn(&s2, "R3").element_balanced);
}

#[test]
fn r4_charge_conservation_and_r5_missing_electron() {
    let s = snap(&engine());
    let r4 = find_rxn(&s, "R4");
    assert!(r4.charge_balanced, "Fe2+ -> Fe3+ + e- 电荷守恒");
    assert!(r4.element_balanced);

    let r5 = find_rxn(&s, "R5");
    assert!(!r5.charge_balanced, "Fe2+ -> Fe3+ 缺电子，电荷不守恒");
    let charge = r5.conservation.iter().find(|c| c.key == "charge").unwrap();
    assert_eq!(charge.delta, "-1");
    // 修复候选：补电子（反应物侧 + e-）。
    let has_electron = r5
        .candidates
        .iter()
        .flat_map(|c| c.actions.iter())
        .any(|a| a.mol == "e-" && a.coef == "1");
    assert!(has_electron, "应枚举补 e- 的修复");
}

#[test]
fn invalid_mapping_is_rejected() {
    let bad = r#"
reaction BAD
  reactant H2O
  product H2O
  map
    H2O:0 -> H2O:0
    H2O:1 -> H2O:1
  end
end
"#;
    // 解析层允许（仅校验存在性），但等价类标记 identity 失败。
    let ds = parser::parse(bad).unwrap();
    let e = Engine::new(ds);
    let s = snap(&e);
    let r = &s.reactions[0];
    assert!(r.map_classes.iter().all(|c| !c.identity_ok));
}

#[test]
fn type_isotope_charge_mismatch_in_map_is_parse_error() {
    let bad = r#"
molecule X1
  atom O 0 0
end
molecule Y1
  atom H 0 0
end
reaction BAD
  reactant X1
  product Y1
  map
    X1:0 -> Y1:0
  end
end
"#;
    // D = H@2，H = H@1，同位素不一致 => 解析期即拒绝。
    assert!(parser::parse(bad).is_err());
}

#[test]
fn lock_map_requires_existing_class() {
    let mut e = engine();
    assert!(e
        .apply(&Op::LockMap {
            reaction: "R1".to_string(),
            class_id: "deadbeef".to_string(),
        })
        .is_err());
    let s = snap(&e);
    let cid = find_rxn(&s, "R1").map_classes[0].class_id.clone();
    e.apply(&Op::LockMap {
        reaction: "R1".to_string(),
        class_id: cid.clone(),
    })
    .unwrap();
}

#[test]
fn network_keeps_multiple_flow_solutions_across_branches() {
    let s = snap(&engine());
    let n = s.networks.iter().find(|n| n.id == "N1").unwrap();
    let branches: std::collections::BTreeSet<&str> =
        n.solutions.iter().map(|x| x.branch.as_str()).collect();
    assert!(
        branches.len() >= 2,
        "不同方向分支必须分别保留: {branches:?}"
    );

    // forward/forward 分支：稳态 f(R6)=f(R7)+f(R8)，有界网格 {0,1/2,1} 上共 6 解。
    let ff = n
        .solutions
        .iter()
        .filter(|x| x.branch.contains("R7=forward") && x.branch.contains("R8=forward"))
        .collect::<Vec<_>>();
    assert_eq!(ff.len(), 6, "稳态 f6=f7+f8 在 {{0,1/2,1}} 上有 6 解");
    for sol in &ff {
        let f6 = frac_num(
            sol.fluxes
                .iter()
                .find(|(k, _)| k == "R6")
                .unwrap()
                .1
                .as_str(),
        );
        let f7 = frac_num(
            sol.fluxes
                .iter()
                .find(|(k, _)| k == "R7")
                .unwrap()
                .1
                .as_str(),
        );
        let f8 = frac_num(
            sol.fluxes
                .iter()
                .find(|(k, _)| k == "R8")
                .unwrap()
                .1
                .as_str(),
        );
        // 全部以半单位表示：f6 = f7 + f8（X 净通量为零）。
        assert_eq!(f6, f7 + f8, "中间体 X 稳态被破坏: {:?}", sol.fluxes);
    }
    // 至少包含一个分数流量解，证明分数被精确保留。
    assert!(ff
        .iter()
        .any(|s| s.fluxes.iter().any(|(_, v)| v.contains('/'))));

    // 所有解中间体 X 净通量都为零。
    for sol in &n.solutions {
        assert!(sol.intermediates.iter().any(|(m, v)| m == "X" && v == "0"));
    }
}

/// 把有理流量统一换算为“半单位”整数，便于断言。
fn frac_num(s: &str) -> i64 {
    if let Some((a, b)) = s.split_once('/') {
        a.parse::<i64>().unwrap() * (2 / b.parse::<i64>().unwrap())
    } else {
        s.parse::<i64>().unwrap() * 2
    }
}

#[test]
fn sqlite_reset_and_replay_reproduces_state() {
    let mut store = Store::in_memory().unwrap();
    store.reset_import(FIXTURE_NAME, FIXTURE).unwrap();

    store
        .append_op(&Op::AddParticipant {
            reaction: "R3".to_string(),
            side: Side::Reactant,
            mol: "O2".to_string(),
            coef: "1/2".to_string(),
        })
        .unwrap();
    let cid = {
        let (e, _, _) = store.rebuild_engine().unwrap();
        find_rxn(&e.snapshot(store.log_entries().unwrap()), "R1").map_classes[0]
            .class_id
            .clone()
    };
    store
        .append_op(&Op::LockMap {
            reaction: "R1".to_string(),
            class_id: cid.clone(),
        })
        .unwrap();

    let (e, name, sha) = store.rebuild_engine().unwrap();
    assert_eq!(name, FIXTURE_NAME);
    assert_eq!(sha.len(), 64);
    let s = e.snapshot(store.log_entries().unwrap());
    assert!(find_rxn(&s, "R3").element_balanced, "重放后 R3 修复仍生效");
    assert_eq!(
        find_rxn(&s, "R1").locked_class.as_deref(),
        Some(cid.as_str())
    );

    // 清空后重新导入：操作日志被重置，R3 回到未守恒。
    store.reset_import(FIXTURE_NAME, FIXTURE).unwrap();
    let (e2, _, _) = store.rebuild_engine().unwrap();
    let s2 = e2.snapshot(store.log_entries().unwrap());
    assert!(!find_rxn(&s2, "R3").element_balanced);
    assert_eq!(s2.log.len(), 1, "重置后只剩 reset 标记操作");
}

#[test]
fn sha256_matches_known_vector() {
    // 空串的标准 SHA-256。
    assert_eq!(
        conservation_bench::storage::sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
}
