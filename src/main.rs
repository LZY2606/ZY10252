//! 反应守恒审查台 —— 本地服务入口。
//!
//! 演示：`cargo run --locked -- --listen 127.0.0.1:5592`

use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use conservation_bench::storage;
use conservation_bench::web::{self, AppState};
use conservation_bench::{FIXTURE, FIXTURE_NAME};

struct Args {
    listen: SocketAddr,
    db: PathBuf,
}

fn parse_args() -> Result<Args, String> {
    let mut args = Args {
        listen: "127.0.0.1:5592".parse().unwrap(),
        db: PathBuf::from("data/conservation.sqlite"),
    };
    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--listen" => {
                let v = it.next().ok_or("--listen 需要地址")?;
                args.listen = v.parse().map_err(|e| format!("监听地址非法: {e}"))?;
            }
            "--db" => {
                args.db = it.next().ok_or("--db 需要路径")?.into();
            }
            "-h" | "--help" => {
                println!("用法: conservation-bench [--listen 127.0.0.1:5592] [--db data/conservation.sqlite]");
                std::process::exit(0);
            }
            other => return Err(format!("未知参数 {other}")),
        }
    }
    Ok(args)
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = parse_args()?;

    if let Some(parent) = args.db.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }

    let mut store = storage::Store::open(args.db.to_str().unwrap())?;
    // 数据库为空（清空后首次启动）时自动导入固定 fixture。
    if store.fixtures()?.is_empty() {
        store.reset_import(FIXTURE_NAME, FIXTURE)?;
        println!("已导入固定 fixture {FIXTURE_NAME}");
    }
    let (_, name, sha) = store.rebuild_engine()?;
    println!("当前 fixture: {name} (sha256 {sha})");

    let state = AppState {
        store: Arc::new(Mutex::new(store)),
        fixture_name: Arc::new(FIXTURE_NAME.to_string()),
        default_fixture: Arc::new(FIXTURE.to_string()),
    };

    let app = web::router(state);
    let listener = tokio::net::TcpListener::bind(args.listen).await?;
    println!("反应守恒审查台已启动：http://{}", args.listen);
    axum::serve(listener, app).await?;
    Ok(())
}
