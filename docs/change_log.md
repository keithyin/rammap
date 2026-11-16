# Changelog

## 2026-11-16 更名为 bandmap 并发布 2.0.0

### 概述

项目在上游 [jwanglab/rammap](https://github.com/jwanglab/rammap) 1.2.0 的基础上扩展，因 `rammap`/`rammap-core` 包名已被上游占用，更名为 **bandmap**并以 2.0.0 首次发布到 crates.io：

- 库 crate：`bandmap-core`（lib 名 `bandmap`）；CLI crate 与二进制：`bandmap`
- workspace 版本统一为 `2.0.0`，`repository` 指向 [keithyin/rammap](https://github.com/keithyin/rammap)

### 破坏性变更

- **环境变量**：`RAMMAP_FORCE_SCALAR`、`RAMMAP_FORCE_SSE`、`RAMMAP_FORCE_AVX2`、`RAMMAP_FORCE_AVX512`、`RAMMAP_FORCE_SCALAR_CHAIN`、`RAMMAP_COMPARE_SCALAR`、`RAMMAP_DP_CACHE_CAP_MB` 全部更名为 `BANDMAP_*` 前缀（同名语义不变）
- **二进制名**：`rammap` → `bandmap`，SAM `@PG` 头的 `ID`/`PN` 字段同步变更
- **依赖写法**：`rammap-core = "1.2.0"` → `bandmap-core = "2.0.0"`，代码导入 `use rammap::` → `use bandmap::`
- 索引文件格式（RMMI magic、`.rmmi`/`.mmi`/`.idx` 扩展名）保持不变，已有索引文件仍可加载

### 配套变更

- 新增 GitHub Actions CI（全平台 build/test/lint/package）与 Release（tag 触发、5 目标二进制、按序发布 bandmap-core → bandmap 到 crates.io）
- WASM 构建产物更名为 `pkg/bandmap.js` / `pkg/bandmap_bg.wasm`，web demo 与构建脚本同步更新

## 2026-11-16 对角线直方图 top-K anchor 过滤

### 概述

将原型 `prototype/anchor_filter/diag_topk.py`（`diagonal_filter`，topk 模式）实现到 bandmap-core，插入 anchor（seed hit）收集之后、DP 链化之前，缓解低复杂/重复区域高频 seed 导致的 anchor 组合爆炸（典型场景：pacbio CCS subreads 到 draft 的比对，比对集中在单一对角线上）。

### 变更文件

- **[NEW]** `bandmap-core/src/align/diag_filter.rs`
  - `DiagFilterParams`：过滤参数（默认：关闭、bucket_size=100、topk=2、neighbor=±1）
  - `filter_anchors_diagonal()`：按 `ref_id_strand()` 连续分组 → `d = ref_pos - query_pos` 分桶 → 计数降序取 top-K（平局按桶号升序）→ 保留选桶 ±neighbor 内 anchor；in-place 复制保持 x 排序，输出可直接喂给链化
  - 7 个单元测试：top-K 过滤、邻桶扩展、平局选桶、分组独立性（ref+strand）、负对角线 floor 语义、空/禁用边界、默认值对照原型
- **[MODIFY]** `bandmap-core/src/align/mod.rs`：声明 `pub mod diag_filter`
- **[MODIFY]** `bandmap-core/src/align/map.rs`
  - `ChainingParams` 新增 `diag_filter` 字段，`MapOptions::default()` 补默认值
  - 仅在 `map_query` / `map_query_multi` 两条主路径链化调用前插入过滤
  - 6 处测试 `ChainingParams` 构造点补字段（chain.rs / chain_rmq.rs / chain_simple.rs / chain_simd.rs ×3）
- **[MODIFY]** `bandmap/src/main.rs`：新增 4 个 CLI 参数

### CLI 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--diag-filter` | 关闭 | 启用对角线直方图 top-K anchor 过滤 |
| `--diag-bucket-size` | 100 | 对角线分桶宽度 B（bp） |
| `--diag-topk` | 2 | 每个 ref/strand 分组保留的 top-K 桶数 |
| `--diag-neighbor` | 1 | 选桶 ±neighbor 邻桶内的 anchor 一并保留 |

### 设计决策

- **分组粒度**：按 `ref_id + strand` 分组独立选桶（多参考基因组语义）；实际位布局为 strand 位于 x 的 bit 63（`seed.rs` 负链时 `x |= 1<<63`），分组键 `ref_id_strand()`（x >> 32）天然涵盖二者。
- **插入路径范围**：仅两条主链化路径过滤。high-occ re-chain 与 RMQ rescue 路径不过滤——前者是救援路径（mid_occ 过滤导致 segment 未链上时用更高 max_occ 重收集 anchor），叠加过滤与救援意图相悖；后者输入为已过滤 anchor 子集。额外收益：主路径误伤导致 segment 链不上时，re-chain 用全量 anchor 重链化，天然形成 segment 级兜底（等价原型 `--fallback` 思路）。
- **负对角线分桶**：`d` 可能为负，用 `div_euclid` 对齐 Python `//` 的 floor 语义（Rust `/` 对负数是截断）。
- **排序保持**：anchors 输入契约要求按 x 排序（两种 seed 收集路径均保证），过滤用 in-place 复制保持相对顺序与分组连续性。
- **平局规则**：计数相同时按桶号升序选桶，与原型行为一致、可复现。
- **默认关闭**：不改变既有行为，需 CLI 显式开启。

### 验证

- `cargo build --workspace` 无错误/警告；`cargo test --workspace` 全部 144 个测试通过（含 7 个新增单元测试）。
- `test_data/low_complexity_*`（61 reads）端到端对比：`--diag-filter` 开启后 62/62 主比对坐标与 baseline 完全一致，仅个别次优链分数（PAF `s2` 字段）下降，符合过滤掉小对角线噪声的预期。
- 已知限制：大缺失导致对角线跨桶跳变时可能截断比对（1/61 reads 左端被截），`--diag-neighbor 2` 可完全恢复；可通过 CLI 调参缓解。
