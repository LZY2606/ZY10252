//! HTTP 冒烟：页面标题、快照、操作、导出、清空重导入。

use conservation_bench::storage::Store;
use conservation_bench::web::{router, AppState};
use conservation_bench::{FIXTURE, FIXTURE_NAME};
use http_body_util::BodyExt;
use std::sync::{Arc, Mutex};
use tower::ServiceExt;

fn app() -> axum::Router {
    let mut store = Store::in_memory().unwrap();
    store.reset_import(FIXTURE_NAME, FIXTURE).unwrap();
    let state = AppState {
        store: Arc::new(Mutex::new(store)),
        fixture_name: Arc::new(FIXTURE_NAME.to_string()),
        default_fixture: Arc::new(FIXTURE.to_string()),
    };
    router(state)
}

async fn get_text(app: &axum::Router, uri: &str) -> (axum::http::StatusCode, String) {
    let resp = app
        .clone()
        .oneshot(
            axum::http::Request::builder()
                .uri(uri)
                .body(axum::body::Body::empty())
                .unwrap(),
        )
        .await
        .unwrap();
    let status = resp.status();
    let bytes = resp.into_body().collect().await.unwrap().to_bytes();
    (status, String::from_utf8(bytes.to_vec()).unwrap())
}

async fn post_form(app: &axum::Router, uri: &str, form: &str) -> (axum::http::StatusCode, String) {
    let resp = app
        .clone()
        .oneshot(
            axum::http::Request::builder()
                .method("POST")
                .uri(uri)
                .header("content-type", "application/x-www-form-urlencoded")
                .body(axum::body::Body::from(form.to_string()))
                .unwrap(),
        )
        .await
        .unwrap();
    let status = resp.status();
    let bytes = resp.into_body().collect().await.unwrap().to_bytes();
    (status, String::from_utf8(bytes.to_vec()).unwrap())
}

#[tokio::test]
async fn index_has_title() {
    let app = app();
    let (status, body) = get_text(&app, "/").await;
    assert_eq!(status, 200);
    assert!(body.contains("反应守恒审查台"), "页面必须包含标题");
}

#[tokio::test]
async fn snapshot_reports_fixture_and_r1_class() {
    let app = app();
    let (status, body) = get_text(&app, "/api/snapshot").await;
    assert_eq!(status, 200);
    let v: serde_json::Value = serde_json::from_str(&body).unwrap();
    assert_eq!(v["fixture"], FIXTURE_NAME);
    assert_eq!(v["sha256"].as_str().unwrap().len(), 64);
    // R1 恰好一个映射等价类，multiplicity=2。
    let r1 = v["snapshot"]["reactions"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["id"] == "R1")
        .unwrap();
    assert_eq!(r1["map_classes"].as_array().unwrap().len(), 1);
    assert_eq!(r1["map_classes"][0]["multiplicity"], 2);
}

#[tokio::test]
async fn apply_repair_and_reset() {
    let app = app();
    // 修复 R3：补 1/2 O2。
    let (st, body) = post_form(
        &app,
        "/api/ops/add-participant",
        "reaction=R3&side=reactant&mol=O2&coef=1%2F2",
    )
    .await;
    assert_eq!(st, 200, "{body}");
    let (_, snap) = get_text(&app, "/api/snapshot").await;
    let v: serde_json::Value = serde_json::from_str(&snap).unwrap();
    let r3 = v["snapshot"]["reactions"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["id"] == "R3")
        .unwrap();
    assert_eq!(r3["element_balanced"], true);

    // 非法操作返回 4xx。
    let (st_bad, _) = post_form(&app, "/api/ops/lock", "reaction=R1&class_id=nope").await;
    assert!(st_bad.is_client_error());

    // 导出包含日志。
    let (_, exp) = get_text(&app, "/api/export").await;
    let ev: serde_json::Value = serde_json::from_str(&exp).unwrap();
    assert!(ev["log"]
        .as_array()
        .unwrap()
        .iter()
        .any(|e| e["summary"].as_str().unwrap().contains("O2")));

    // 清空重置后 R3 回到未守恒。
    let resp = app
        .clone()
        .oneshot(
            axum::http::Request::builder()
                .method("POST")
                .uri("/api/reset")
                .body(axum::body::Body::empty())
                .unwrap(),
        )
        .await
        .unwrap();
    assert_eq!(resp.status(), 200);
    let (_, snap2) = get_text(&app, "/api/snapshot").await;
    let v2: serde_json::Value = serde_json::from_str(&snap2).unwrap();
    let r3b = v2["snapshot"]["reactions"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["id"] == "R3")
        .unwrap();
    assert_eq!(r3b["element_balanced"], false);
}

#[tokio::test]
async fn reimport_rejects_invalid_dsl() {
    let app = app();
    let resp = app
        .clone()
        .oneshot(
            axum::http::Request::builder()
                .method("POST")
                .uri("/api/reimport")
                .header("content-type", "application/json")
                .body(axum::body::Body::from("{\"content\":\"reaction X\\n\"}"))
                .unwrap(),
        )
        .await
        .unwrap();
    assert!(resp.status().is_client_error());
}
