//! SQLite 持久层：操作日志（事件溯源）+ 导入的 fixture 文本与校验和。
//!
//! 清空数据库后可用保存的 fixture 文本重新导入并重放全部操作复核。

use crate::engine::{Engine, LogEntry};
use crate::model::Op;
use crate::parser;
use rusqlite::Connection;

pub struct Store {
    pub conn: Connection,
}

const SCHEMA: &str = r#"
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS ops_log (
    seq     INTEGER PRIMARY KEY AUTOINCREMENT,
    at      TEXT NOT NULL,
    op_json TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS fixtures (
    name    TEXT PRIMARY KEY,
    content TEXT NOT NULL,
    sha256  TEXT NOT NULL,
    imported_at TEXT NOT NULL
);
"#;

pub fn sha256_hex(bytes: &[u8]) -> String {
    // 不引第三方哈希库：用 std 手写 SHA-256。
    let k: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];

    let bit_len = (bytes.len() as u64).wrapping_mul(8);
    let mut msg = bytes.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[i * 4],
                chunk[i * 4 + 1],
                chunk[i * 4 + 2],
                chunk[i * 4 + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
            (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(k[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

impl Store {
    pub fn open(path: &str) -> Result<Self, String> {
        let conn = Connection::open(path).map_err(|e| format!("打开数据库失败: {e}"))?;
        conn.execute_batch(SCHEMA)
            .map_err(|e| format!("初始化 schema 失败: {e}"))?;
        Ok(Store { conn })
    }

    /// 内存数据库（测试用）。
    pub fn in_memory() -> Result<Self, String> {
        let conn = Connection::open_in_memory().map_err(|e| e.to_string())?;
        conn.execute_batch(SCHEMA).map_err(|e| e.to_string())?;
        Ok(Store { conn })
    }

    /// 清空并重新导入 fixture，同时清掉旧操作日志（重置复核）。
    pub fn reset_import(&mut self, name: &str, content: &str) -> Result<(), String> {
        self.conn
            .execute_batch("DELETE FROM ops_log; DELETE FROM fixtures; DELETE FROM meta;")
            .map_err(|e| e.to_string())?;
        self.put_fixture(name, content)?;
        let op = Op::ResetImported {
            source: name.to_string(),
        };
        self.append_op(&op)?;
        Ok(())
    }

    pub fn put_fixture(&mut self, name: &str, content: &str) -> Result<(), String> {
        let sha = sha256_hex(content.as_bytes());
        self.conn
            .execute(
                "INSERT INTO fixtures(name, content, sha256, imported_at)
                 VALUES (?1, ?2, ?3, ?4)
                 ON CONFLICT(name) DO UPDATE SET content=excluded.content,
                     sha256=excluded.sha256, imported_at=excluded.imported_at",
                rusqlite::params![name, content, sha, now_string()],
            )
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn fixtures(&self) -> Result<Vec<(String, String, String)>, String> {
        let mut stmt = self
            .conn
            .prepare("SELECT name, content, sha256 FROM fixtures ORDER BY name")
            .map_err(|e| e.to_string())?;
        let rows = stmt
            .query_map([], |r| {
                Ok((
                    r.get::<_, String>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                ))
            })
            .map_err(|e| e.to_string())?;
        let mut out = Vec::new();
        for row in rows {
            out.push(row.map_err(|e| e.to_string())?);
        }
        Ok(out)
    }

    pub fn append_op(&mut self, op: &Op) -> Result<i64, String> {
        let json = serde_json::to_string(op).map_err(|e| e.to_string())?;
        self.conn
            .execute(
                "INSERT INTO ops_log(at, op_json) VALUES (?1, ?2)",
                rusqlite::params![now_string(), json],
            )
            .map_err(|e| e.to_string())?;
        Ok(self.conn.last_insert_rowid())
    }

    pub fn log_entries(&self) -> Result<Vec<LogEntry>, String> {
        let mut stmt = self
            .conn
            .prepare("SELECT seq, at, op_json FROM ops_log ORDER BY seq")
            .map_err(|e| e.to_string())?;
        let rows = stmt
            .query_map([], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                ))
            })
            .map_err(|e| e.to_string())?;
        let mut out = Vec::new();
        for row in rows {
            let (seq, at, json) = row.map_err(|e| e.to_string())?;
            let summary = summarize(&json);
            out.push(LogEntry {
                seq,
                at,
                op_json: json,
                summary,
            });
        }
        Ok(out)
    }

    /// 读取 fixture 文本、解析、重放操作日志，重建引擎。
    pub fn rebuild_engine(&self) -> Result<(Engine, String, String), String> {
        let (name, content, sha) = self
            .fixtures()?
            .into_iter()
            .next()
            .ok_or_else(|| "数据库中没有 fixture，请先导入".to_string())?;
        let dataset = parser::parse(&content)?;
        let mut engine = Engine::new(dataset);
        for entry in self.log_entries()? {
            let op: Op =
                serde_json::from_str(&entry.op_json).map_err(|e| format!("操作日志损坏: {e}"))?;
            engine.apply(&op)?;
        }
        Ok((engine, name, sha))
    }
}

fn summarize(op_json: &str) -> String {
    match serde_json::from_str::<Op>(op_json) {
        Ok(Op::LockMap { reaction, class_id }) => {
            format!("锁定映射：反应 {reaction} 选择等价类 {class_id}")
        }
        Ok(Op::UnlockMap { reaction }) => format!("解锁映射：反应 {reaction}"),
        Ok(Op::AddParticipant {
            reaction,
            side,
            mol,
            coef,
        }) => {
            let side = match side {
                crate::model::Side::Reactant => "反应物",
                crate::model::Side::Product => "产物",
            };
            format!("添加参与小分子：反应 {reaction} 的{side}侧 + {coef} {mol}")
        }
        Ok(Op::SetDirectionUnknown { reaction, unknown }) => format!(
            "反应 {reaction} 方向{}",
            if unknown {
                "标记为未定"
            } else {
                "恢复已定"
            }
        ),
        Ok(Op::ResetImported { source }) => format!("重置并重新导入 fixture：{source}"),
        Err(_) => "（无法解析的操作）".to_string(),
    }
}

fn now_string() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    format!("unix:{secs}")
}
