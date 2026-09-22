# 反应守恒审查台

对**单反应**和**小型反应网络**做元素、电荷、同位素与原子身份守恒审查的本地服务。
技术栈：Rust + Axum + SQLite（`rusqlite/bundled`）+ 原生 Web 页面，无前端构建步骤。

## 安装与演示

```bash
cargo fetch --locked
cargo test --locked
cargo run --locked -- --listen 127.0.0.1:5592
```

浏览器打开 <http://127.0.0.1:5592>，页面标题为「**反应守恒审查台**」。

其它参数：

```bash
cargo run --locked -- --listen 127.0.0.1:5592 --db data/conservation.sqlite
```

数据库文件默认放在 `data/conservation.sqlite`，目录不存在会自动创建。
首次启动（或库被清空后）自动导入随二进制分发的固定 fixture `fixtures/sample.rxn`。

## 数据口径（重要）

- **全程有理数，不用浮点。** 所有计量系数、未守恒向量、网络流量都用
  `num_rational::Ratio<i64>` 表示；分数写成 `n/d`（如 `1/2`），不存在“容差放过一个原子”。
- **原子语义三元组** = `(元素, 质量数, 电荷)`。
  - 同位素用**质量数**区分：`H@2` 是氘。DSL 中可直接写 `D`（默认质量数 2）或 `T`（默认 3），
    它们在守恒分量上归一为氢的同位素（`H@2` / `H@3`），在原子身份上仍属于氢。
  - `atom` 行字段顺序为 `元素 [质量数] [形式电荷]`，省略的字段为 0。
- **电荷**：分子总电荷 = 各原子形式电荷之和，守恒向量单列 `charge` 分量。
  `e-` 是仅用于半反应的**伪分子**：它只计入电荷（−1），不计入任何元素组成。
- **未守恒向量** = 反应物加权组成 − 产物加权组成，按分量给出有理 `Δ`；`Δ=0` 才守恒。
- **原子身份 / 映射**：
  - 映射边形如 `H2O:1 -> HDO:2`（`分子id:原子下标`），原子下标从 0 开始。
  - 映射必须是**单射**：每个反应物/产物原子至多出现一次；含映射的反应要求单位计量且全覆盖。
  - 身份只要求**元素一致**；同位素/电荷是否守恒由守恒向量独立核对（氢可以映射到带质量数的氢）。
  - 分子对称产生的多个编号映射**不会假装唯一**：系统枚举分子图自同构，
    把在 `G_反应物 × G_产物` 群作用下等价的映射归并为一个**等价类**，并给出
    `multiplicity`（等价映射数）。
- **有界候选修复**：对未守恒反应，在有界网格上枚举“向反应物/产物补入正系数内置小分子”
  的 1–2 个动作。系数网格 `{1/2, 1, 3/2, 2}`，候选空间有界；按动作集合去重，
  每个候选保留**修复动作**与**证据代价**（分数、动作数、系数规模、跨侧对消都会抬升代价）。
- **网络中间体流量**：
  - 对每条标记 `direction unknown` 的反应枚举方向分支；
  - 在有界有理网格 `{0, 1/2, 1}`（reverse 分支取负）上求稳态，使每个中间体净通量为 0；
  - **不同方向分支分别保留**所有可行流量解（即使流量向量相同），页面逐分支列出。

## 受限本地格式（`.rxn`）

`#` 起始为注释。缩进只用于阅读，块以 `end` 结束，`map` 块可嵌在 `reaction` 内。

```text
molecule <id>
  atom <元素|D|T> [质量数] [形式电荷]
  bond <u> <v> [键级1..3，默认1]
end

reaction <id>
  reactant <分子id> [系数，默认1，可写 n/d]
  product  <分子id> [系数，可写 n/d]
  map
    <反应物分子>:<原子下标> -> <产物分子>:<原子下标>
  end
  direction unknown          # 可选：方向未定
end

network <id>
  step <反应id> [unknown]
end
```

内置小分子库（修复候选可直接引用）：`H2O HDO D H OH- H+ H2 O2 N2 NO N2O NO2 N2O4
CO CO2 Fe2+ Fe3+ e-`，以及网络抽象池物种 `A X B`（无组成，仅用于分子级流量平衡）。

## 固定 fixture 内容（验收点）

`fixtures/sample.rxn`：

| 反应 | 覆盖点 |
| --- | --- |
| `R1: H2O + D -> HDO + H` | 同位素标签守恒；水的两个 H 对称 => **两个等价映射归并为一个等价类（multiplicity=2）** |
| `R2: NO + 1/2 O2 -> NO2` | **分数化学计量**按有理数核对守恒 |
| `R3: NO -> NO2` | **遗漏小分子**（缺 1 个 O），首选修复为反应物补 `1/2 O2` |
| `R4: Fe2+ -> Fe3+ + e-` | 半反应电荷守恒（电子只计电荷） |
| `R5: Fe2+ -> Fe3+` | 遗漏电子，电荷 Δ=−1，候选含产物补 `e-` |
| `R6/R7/R8` + 网络 `N1` | `R7/R8` **方向未定**；稳态 `f6=f7+f8`，逐分支保留多个有理流量解 |

## 操作与可重放性

页面支持：

- **锁定一组映射**：在反应卡片中选择某个对称等价类（也可解锁）；
- **添加明确的参与小分子**：一键应用候选修复，或等价地走操作接口；
- **标记方向未定 / 恢复已定**：反应卡片按钮，或网络 `step ... unknown`；
- **导出运行记录**：`GET /api/export` 导出 JSON（fixture 文本、SHA-256、带序号的操作日志）；
- **清空后重新导入复核**：页面「清空并重新导入」或 `POST /api/reset`；
  也可在「自定义重新导入」粘贴 `.rxn` 文本（`POST /api/reimport`）。

持久化采用**事件溯源**：SQLite 只保存 fixture 文本、SHA-256 与操作日志。
重建状态 = 重新解析 fixture + 按序重放全部操作，因此清空库后导入同一份 fixture
即可复现；导出的 JSON 可离线核对每一步。

### HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/api/snapshot` | 完整审查快照（分子图、守恒向量、映射等价类、候选、网络流量、日志） |
| POST | `/api/ops/lock` | `reaction, class_id` 锁定映射等价类 |
| POST | `/api/ops/unlock` | `reaction` 解锁 |
| POST | `/api/ops/add-participant` | `reaction, side(reactant|product), mol, coef` |
| POST | `/api/ops/direction` | `reaction, unknown(true|false)` |
| GET | `/api/export` | 导出 fixture + 校验和 + 操作日志 |
| POST | `/api/reset` | 清空并重新导入内置固定 fixture |
| POST | `/api/reimport` | JSON `{ "content": "...", "name": "..." }` 自定义重导入 |

## 数据库与清空复核

```bash
# 清空数据库后重新启动，会自动重新导入固定 fixture
python3 -c "import shutil; shutil.rmtree('data', ignore_errors=True)"
cargo run --locked -- --listen 127.0.0.1:5592
```

## 测试

```bash
cargo test --locked
```

- `tests/conservation.rs`：有理精确性、分数计量、同位素与对称映射等价类、
  遗漏小分子/电子的候选、身份单射、网络多分支流量解、SQLite 清空重放、SHA-256。
- `tests/http_api.rs`：页面标题、快照、应用修复、非法操作拒绝、导出、重置、重导入。

## 目录结构

```text
fixtures/sample.rxn   固定 fixture
src/chem.rs           原子/分子图、有理数、图自同构（对称性）、内置小分子库
src/parser.rs         受限本地格式解析
src/model.rs          反应/网络/操作模型
src/engine.rs         守恒向量、映射等价类、身份守恒、候选枚举、网络流量
src/storage.rs        SQLite 事件溯源、SHA-256、重建
src/web.rs            Axum 路由与操作接口
static/               操作页面（HTML/CSS/JS，随二进制内嵌）
```
