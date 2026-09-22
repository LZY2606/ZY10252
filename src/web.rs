//! Axum Web 层：页面、快照 JSON、操作接口、导出/导入/重置。

use crate::model::{Op, Side};
use crate::storage::Store;
use axum::{
    extract::State,
    http::{header, StatusCode},
    response::{IntoResponse, Response},
    routing::{get, post},
    Form, Json, Router,
};
use serde::Deserialize;
use std::sync::{Arc, Mutex};

#[derive(Clone)]
pub struct AppState {
    pub store: Arc<Mutex<Store>>,
    pub fixture_name: Arc<String>,
    pub default_fixture: Arc<String>,
}

pub fn router(state: AppState) -> Router {
    Router::new()
        .route("/", get(index))
        .route("/static/app.js", get(app_js))
        .route("/static/style.css", get(style_css))
        .route("/api/snapshot", get(snapshot))
        .route("/api/ops", get(list_ops))
        .route("/api/ops/lock", post(lock_map))
        .route("/api/ops/unlock", post(unlock_map))
        .route("/api/ops/add-participant", post(add_participant))
        .route("/api/ops/direction", post(set_direction))
        .route("/api/reset", post(reset))
        .route("/api/export", get(export))
        .route("/api/reimport", post(reimport))
        .with_state(state)
}

async fn index() -> Response {
    html(include_str!("../static/index.html"))
}
async fn app_js() -> Response {
    js(include_str!("../static/app.js"))
}
async fn style_css() -> Response {
    css(include_str!("../static/style.css"))
}

fn html(body: &str) -> Response {
    (
        StatusCode::OK,
        [(header::CONTENT_TYPE, "text/html; charset=utf-8")],
        body.to_string(),
    )
        .into_response()
}
fn js(body: &str) -> Response {
    (
        StatusCode::OK,
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        body.to_string(),
    )
        .into_response()
}
fn css(body: &str) -> Response {
    (
        StatusCode::OK,
        [(header::CONTENT_TYPE, "text/css; charset=utf-8")],
        body.to_string(),
    )
        .into_response()
}

#[derive(serde::Serialize)]
struct SnapshotResponse {
    fixture: String,
    sha256: String,
    snapshot: crate::engine::Snapshot,
}

async fn snapshot(State(st): State<AppState>) -> Result<Json<SnapshotResponse>, ApiError> {
    let guard = st.store.lock().unwrap();
    let (engine, name, sha) = guard.rebuild_engine()?;
    let log = guard.log_entries()?;
    Ok(Json(SnapshotResponse {
        fixture: name,
        sha256: sha,
        snapshot: engine.snapshot(log),
    }))
}

async fn list_ops(State(st): State<AppState>) -> Result<Json<serde_json::Value>, ApiError> {
    let guard = st.store.lock().unwrap();
    let entries = guard.log_entries()?;
    Ok(Json(serde_json::json!({ "log": entries })))
}

#[derive(Deserialize)]
struct LockReq {
    reaction: String,
    class_id: String,
}
async fn lock_map(
    State(st): State<AppState>,
    Form(q): Form<LockReq>,
) -> Result<Json<serde_json::Value>, ApiError> {
    apply(
        st,
        Op::LockMap {
            reaction: q.reaction,
            class_id: q.class_id,
        },
    )
}

#[derive(Deserialize)]
struct RxnReq {
    reaction: String,
}
async fn unlock_map(
    State(st): State<AppState>,
    Form(q): Form<RxnReq>,
) -> Result<Json<serde_json::Value>, ApiError> {
    apply(
        st,
        Op::UnlockMap {
            reaction: q.reaction,
        },
    )
}

#[derive(Deserialize)]
struct AddReq {
    reaction: String,
    side: String,
    mol: String,
    coef: String,
}
async fn add_participant(
    State(st): State<AppState>,
    Form(q): Form<AddReq>,
) -> Result<Json<serde_json::Value>, ApiError> {
    let side = match q.side.as_str() {
        "reactant" => Side::Reactant,
        "product" => Side::Product,
        other => {
            return Err(ApiError(format!(
                "side 仅支持 reactant/product，得到 {other}"
            )))
        }
    };
    apply(
        st,
        Op::AddParticipant {
            reaction: q.reaction,
            side,
            mol: q.mol,
            coef: q.coef,
        },
    )
}

#[derive(Deserialize)]
struct DirReq {
    reaction: String,
    unknown: String,
}
async fn set_direction(
    State(st): State<AppState>,
    Form(q): Form<DirReq>,
) -> Result<Json<serde_json::Value>, ApiError> {
    let unknown = matches!(q.unknown.as_str(), "true" | "1" | "yes");
    apply(
        st,
        Op::SetDirectionUnknown {
            reaction: q.reaction,
            unknown,
        },
    )
}

fn apply(st: AppState, op: Op) -> Result<Json<serde_json::Value>, ApiError> {
    let mut guard = st.store.lock().unwrap();
    {
        // 先在引擎上试算，拒绝非法操作，再写日志。
        let (engine, _, _) = guard.rebuild_engine()?;
        let mut trial = engine;
        trial.apply(&op)?;
    }
    let seq = guard.append_op(&op)?;
    Ok(Json(serde_json::json!({ "ok": true, "seq": seq })))
}

async fn reset(State(st): State<AppState>) -> Result<Json<serde_json::Value>, ApiError> {
    let mut guard = st.store.lock().unwrap();
    let content = st.default_fixture.as_str().to_string();
    let name = st.fixture_name.as_str().to_string();
    guard.reset_import(&name, &content)?;
    Ok(Json(serde_json::json!({ "ok": true })))
}

#[derive(Deserialize)]
struct ReimportReq {
    content: String,
    #[serde(default)]
    name: Option<String>,
}
async fn reimport(
    State(st): State<AppState>,
    Json(q): Json<ReimportReq>,
) -> Result<Json<serde_json::Value>, ApiError> {
    // 先解析确认合法，再落库重置。
    crate::parser::parse(&q.content)?;
    let name = q
        .name
        .unwrap_or_else(|| st.fixture_name.as_str().to_string());
    let mut guard = st.store.lock().unwrap();
    guard.reset_import(&name, &q.content)?;
    Ok(Json(serde_json::json!({ "ok": true })))
}

#[derive(serde::Serialize)]
struct ExportBundle {
    fixture_name: String,
    fixture_sha256: String,
    fixture: String,
    log: Vec<crate::engine::LogEntry>,
    exported_at: String,
}

async fn export(State(st): State<AppState>) -> Result<Json<ExportBundle>, ApiError> {
    let guard = st.store.lock().unwrap();
    let (_, name, sha) = guard.rebuild_engine()?;
    let content = guard
        .fixtures()?
        .into_iter()
        .find(|(n, _, _)| n == &name)
        .map(|(_, c, _)| c)
        .ok_or_else(|| ApiError("fixture 内容缺失".to_string()))?;
    let log = guard.log_entries()?;
    Ok(Json(ExportBundle {
        fixture_name: name,
        fixture_sha256: sha,
        fixture: content,
        log,
        exported_at: String::new(),
    }))
}

struct ApiError(String);
impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        (
            StatusCode::BAD_REQUEST,
            [(header::CONTENT_TYPE, "application/json; charset=utf-8")],
            serde_json::json!({ "ok": false, "error": self.0 }).to_string(),
        )
            .into_response()
    }
}
impl From<String> for ApiError {
    fn from(value: String) -> Self {
        ApiError(value)
    }
}
