# diag-filter：对角线直方图 top-K anchor 过滤

## 概述

`diag-filter` 用于缓解低复杂/重复区域高频 seed 导致的 anchor（seed hit）组合爆炸。典型场景：pacbio CCS subreads 到 draft 的比对，大量 anchor 集中在单一对角线上，噪声对角线上的 anchor 会让 DP 链化的输入规模失控。

核心思路：按 `ref_id + strand` 分组统计 anchor 所在对角线桶的直方图，**只保留 top-K 个最"拥挤"对角线桶（及其 ±neighbor 邻桶）内的 anchor**，其余丢弃。

- 原型来源：`prototype/anchor_filter/diag_topk.py`（`diagonal_filter`，topk 模式）
- Rust 实现：`bandmap-core/src/align/diag_filter.rs`
- 插入位置：anchor 收集（`collect_seed_hits`）之后、DP 链化（`chain_anchors` / `chain_anchors_rmq`）之前

## CLI 使用

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--diag-filter` | 关闭 | 启用对角线直方图 top-K anchor 过滤（默认关闭，不改变既有行为） |
| `--diag-bucket-size` | 100 | 对角线分桶宽度 B（bp） |
| `--diag-topk` | 2 | 每个 ref/strand 分组保留的 top-K 桶数 |
| `--diag-neighbor` | 1 | 选桶 ±neighbor 邻桶内的 anchor 一并保留 |

示例：

```bash
# 显式启用过滤（默认参数：bucket_size=100, topk=2, neighbor=±1）
bandmap --diag-filter -x map-hifi ref.fa reads.fq > out.paf

# 调参：桶宽 500，保留 top-3 桶，邻桶范围 ±2
bandmap --diag-filter --diag-bucket-size 500 --diag-topk 3 --diag-neighbor 2 ref.fa reads.fq > out.paf
```

调参建议：

- 大缺失（deletion）会使真实比对的对角线跨桶跳变，可能被截断；增大 `--diag-neighbor`（如 2）可恢复。
- 重复序列较多时适当增大 `--diag-topk`，避免把次要但真实的比对方向过滤掉。

## 算法流程

对每个 `ref_id + strand` 分组（输入 anchors 已按 `x` 排序，因此分组在数组中连续）：

1. **计算对角线**：`d = ref_pos - query_pos`。query 坐标在 seed 收集时已对负链做过镜像，因此正负链公式一致。
2. **分桶计数**：`b = d.div_euclid(B)`。使用 `div_euclid` 对齐 Python `//` 的 floor 语义（Rust `/` 对负数是截断），负对角线才能正确分桶（如 `d = -150` → 桶 `-2`）。
3. **选 top-K 桶**：计数降序排序，平局按桶号升序（与原型一致，保证可复现），取前 K 个桶。
4. **保留 anchor**：保留所在桶 ∈ 选桶 ± `neighbor` 邻桶集合内的 anchor。

过滤为 **in-place 压缩**（双指针复制），保持 anchor 的相对顺序，输出仍按 `x` 排序、分组仍连续，可直接喂给链化。返回值为保留的 anchor 数。

边界行为：

- `enabled = false`、`topk = 0` 或 `bucket_size <= 0` 时为 no-op，直接返回原长度。
- `neighbor <= 0` 时不做邻桶扩展，仅保留选桶本身。

## 库 API

```rust
use bandmap::align::diag_filter::{filter_anchors_diagonal, DiagFilterParams};
// anchors: Vec<Minimizer>，由 collect_seed_hits 等 seed 收集路径产出（按 x 排序）
let mut anchors = collect_seed_hits(...);

let p = DiagFilterParams {
    enabled: true,
    bucket_size: 100, // 对角线分桶宽度 B（bp）
    topk: 2,          // 每个 ref/strand 分组保留的 top-K 桶数
    neighbor: 1,      // 选桶 ±neighbor 邻桶一并保留
};
let kept = filter_anchors_diagonal(&mut anchors, &p); // 返回保留数量，anchors 原地过滤
```

`DiagFilterParams::default()` 与原型一致：`enabled = false`、`bucket_size = 100`、`topk = 2`、`neighbor = 1`。

## 管线集成

### 配置入口

`ChainingParams`（`bandmap-core/src/align/map.rs`）新增 `diag_filter: DiagFilterParams` 字段，默认值为 `DiagFilterParams::default()`（即关闭）。

CLI（`bandmap/src/main.rs`）在组装 `MapOptions` 时把 4 个命令行参数写入：

```rust
opt.chaining.diag_filter.enabled = cli.diag_filter;
opt.chaining.diag_filter.bucket_size = cli.diag_bucket_size;
opt.chaining.diag_filter.topk = cli.diag_topk;
opt.chaining.diag_filter.neighbor = cli.diag_neighbor;
```

### 插入点

仅在两条主比对路径的链化调用前各插入一次过滤（`map.rs`）：

```rust
collect_seed_hits(opt, mi, ..., &mut ctx.anchors, ...);

// Diagonal histogram top-K anchor filter (--diag-filter; no-op unless enabled).
filter_anchors_diagonal(&mut ctx.anchors, &opt.chaining.diag_filter);

// 链化：按 flags 走 chain_anchors（默认 DP）或 chain_anchors_rmq（RMQ）
let (u, chains) = chain_anchors(...);
```

- `map_query_multi`（多 query segment 合并路径，"strong pairing"，用 `qlen_sum`）与 `map_query`（单 query 路径，用 `qlen`）均覆盖。

### 不过滤的路径（设计决策）

- **high-occ re-chain**：mid_occ 过滤导致 segment 未链上时，用更高 `max_occ` 重收集 anchor 重新链化，属于救援路径；若叠加过滤会与救援意图相悖。
- **RMQ rescue**：输入已是过滤后 anchor 的子集（已链化结果），无需再过滤。

额外收益：主路径过滤若误伤导致某个 segment 链不上，re-chain 会用全量 anchor 重新链化，天然形成 segment 级兜底（等价于原型的 `--fallback` 思路）。

## 输入契约

- anchors 必须按 `x` 排序（`collect_seed_hits` / `collect_seed_hits_heap` 两条 seed 收集路径均保证）。
- 分组依赖 `ref_id_strand()`（`x >> 32`，包含 strand bit 63 + ref_id）的 run 连续性；由 x 排序自然保证。
- 过滤不改变相对顺序，输出仍满足链化的排序契约。

## 测试与验证

`bandmap-core/src/align/diag_filter.rs` 内含 7 个单元测试：

| 测试 | 覆盖点 |
|---|---|
| `keeps_topk_buckets_only` | 只保留 top-K 桶，顺序保持（x 排序不变） |
| `neighbor_expansion_keeps_adjacent_buckets` | 邻桶 ±neighbor 扩展保留相邻桶 anchor |
| `ties_broken_by_smaller_bucket` | 计数平局时按桶号升序选桶 |
| `groups_are_independent_per_ref_and_strand` | ref + strand 分组直方图互相独立 |
| `negative_diagonals_use_floor_semantics` | 负对角线 floor 分桶（对齐 Python `//`） |
| `no_op_when_disabled_or_empty` | 禁用 / 空输入 / 非法 bucket_size 均为 no-op |
| `default_params_match_prototype` | 默认参数与原型一致 |

端到端验证（`test_data/low_complexity_*`，61 reads）：

- 开启 `--diag-filter` 后主比对坐标与 baseline 完全一致，仅个别次优链分数（PAF `s2` 字段）下降，符合过滤掉小对角线噪声的预期。
- 已知限制：1/61 reads 因大缺失导致对角线跨桶跳变，比对左端被截断；`--diag-neighbor 2` 可完全恢复。

## 相关文件

| 文件 | 角色 |
|---|---|
| `bandmap-core/src/align/diag_filter.rs` | 核心实现 + 单元测试 |
| `bandmap-core/src/align/map.rs` | `ChainingParams.diag_filter` 字段、两条主路径插入过滤 |
| `bandmap-core/src/align/mod.rs` | `pub mod diag_filter` 声明 |
| `bandmap/src/main.rs` | 4 个 CLI 参数及配置写入 |
| `prototype/anchor_filter/diag_topk.py` | 原型（行为参照） |
