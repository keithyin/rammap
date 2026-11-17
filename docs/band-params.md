# band 参数调研：`bandwidth` 与 `bandwidth_long` 的设计与实现分析

> 调研日期：2026-11-17
> 调研范围：`bandmap-core/src/align/` 中 band 相关参数（`bandwidth` / `bandwidth_long`）在 seeding→chaining→alignment 全链路的实际用途；并与上游 minimap2（master `align.c`，2026-11-17 抓取）逐行对照。
> 起因问题：① anchor 之间的比对是否使用了 banded alignment？② 这组 band 参数设计得是否合理？

## TL;DR

1. **内核层面是 banded 的**：anchor 间 gap-fill 调用 ksw2 风格的带状 SIMD DP（`dp::extend_single_affine` / `extend_dual_affine`），带 banded traceback 数组。
2. **但默认长读长参数下 band 实际不绑定**：gap-fill 传入的带宽是 `bandwidth_long × 1.5 + 1 = 30001`（map-ont），而相邻 anchor 间距受 `max_gap=5000` 约束——band 恒大于段尺寸，等效全矩阵。真正吃到 band 收益的是**左右端部扩展**（band=751）和 **sr preset**。
3. **这是 minimap2 的忠实移植，不是本地缺陷**：上游 `mm_align1()` 同样用 `bw_long` 做 gap-fill band、`bw` 做端部扩展、`×1.5+1` 系数、long-join 全宽。上游给宽 band 有机制性理由（见 §5.1）。
4. **设计欠账集中在参数治理而非算法结构**：`bw_long` 一个旋钮同时管 gap-fill band 和 RMQ rescue 重链化（收紧一个静默关闭另一个）；CLI `-r` 单值语义不对称；文档承诺的 O(n×bandwidth) 复杂度在 gap-fill 路径名不副实。

---

## 1. 背景

bandmap 是 minimap2 系（seed→chain→gap-fill alignment）的长读长/短读长比对器。链化输出的 anchor（minimizer）链之间存在未比对区域（gap），alignment 阶段逐段填充。band width 是 banded DP 的核心参数，直接决定 gap-fill 的计算量与召回上限。

本仓库有两个 band 参数：

```rust
// bandmap-core/src/align/map.rs:110-111
pub bandwidth: i32,
pub bandwidth_long: i32,
// 默认：500 / 20000（map.rs:205-206）
// architecture.md:622-623 注释：bandwidth="Max diagonal distance in chain"，
//                              bandwidth_long="Bandwidth for long gaps"
```

## 2. 两个参数的用途映射

### 2.1 `bandwidth`（小值）——"局部连接尺度"

| 用途 | 位置 | 说明 |
|---|---|---|
| chaining 对角 gate | `chain_rmq.rs:709,723` | 链边须满足 `|ref_diff − query_diff| ≤ bandwidth` |
| chaining 距离窗 | `chain_rmq.rs:636` | `max_dist = max(max_gap, bandwidth)` |
| 左/右端部扩展 DP band | `extend.rs:1562, 1712-1714, 2002-2004` | `bw = bandwidth × 1.5 + 1`（500→751） |
| `fix_bad_ends` 边界修剪 | `extend.rs:1332` | 以 `bandwidth` 度量链端可信长度 |

### 2.2 `bandwidth_long`（大值）——"长事件尺度"

| 用途 | 位置 | 说明 |
|---|---|---|
| **anchor 间 gap-fill DP band** | `extend.rs:1765, 1826-1877` | `bw1 = bandwidth_long × 1.5 + 1`（20000→30001） |
| long-join gap 全宽 | `extend.rs:1766-1768` | `bw1 = max(gap_q_len, gap_t_len)`，即不设限 |
| **RMQ rescue 重链化带宽** | `map.rs:1047-1071, 1373-1398` | 首链覆盖异常时以 `bandwidth: bandwidth_long` 重跑 `chain_anchors_rmq`；触发前提 `bandwidth_long > bandwidth` 且未设 `SPLICE\|SHORT_READ\|NO_LJOIN` |

守恒约束：`bw_long < bw` 时强制拉齐（`extend.rs:1564`）。

### 2.3 各 preset 取值

| Preset | bandwidth | bandwidth_long | 倍数 | 出处 |
|---|---|---|---|---|
| **默认**（map-ont / map-pb / lr:hq / map-hifi / lr:hqae / map-iclr 均继承） | 500 | 20000 | **40×** | `map.rs:205-206` |
| **asm5/asm10/asm20** | 1000 | 100000 | **100×** | `api.rs:891` |
| sr / short | 100 | 100 | 1× | `api.rs:930` |
| splice / cdna（含 splice:hq / splice:sr） | 200000 | 200000 | 1× | `api.rs:946` |
| ava-ont | 2000 | 2000 | 1× | `api.rs:974` |
| ava-pb | 500（默认） | = bandwidth | 1× | `api.rs:984` |

规律：**短读长 / 全对 / splice 场景两者相等；长读长与组装场景拉开 40–100×**。1× 的 preset 里 `bandwidth_long` 实际退化为冗余参数。

### 2.4 CLI

- `-r <bw>,<bw_long>`（`--bw`，`main.rs:82-84`）：逗号分隔双值。
  ⚠️ **只给一个值时只改 `bandwidth`，`bandwidth_long` 保留 preset 值**（`main.rs:577-578`）。map-ont 下 `-r 2000` 并不会收紧 gap-fill band（仍是 30001）。
- `-G`（max intron len）将两者连同 `max_gap_ref` 一起设为同值（`main.rs:619`）。

## 3. banded alignment 的实际生效情况

### 3.1 gap-fill：banded 内核 + 失效的 band

gap-fill 的 DP 段调用链（`extend.rs:1826-1877`）：

- `dp::extend_single_affine` / `dp::extend_dual_affine`：ksw2 风格带状 SIMD 内核，`_n_col_ = min(min(qlen,tlen), bandwidth+1)`，banded traceback 数组（p/off/off_end，`dp/single.rs:58-65`）——**内核是真的 banded**；
- 传入带宽 `bw1 = bandwidth_long × 1.5 + 1`（`extend.rs:1765`）。

但 map-ont 默认下：段尺寸 ≤ `max_gap` 量级（5000），band = 30001 > 段尺寸 → **band 永不裁剪，等效全矩阵**。gap-fill 的实际最坏复杂度是 **O(max_gap²) ≈ 25M cells/段**（SIMD + z-drop 缓解），而不是参数故事暗示的 O(n × bandwidth)。

### 3.2 真正被 band 约束的路径

- **左右端部扩展**：`bw = bandwidth × 1.5 + 1 = 751`（`extend.rs:1562, 1712-1714, 2002-2004`）——band 有效，因为扩展段一侧有锚定 seed；
- **sr preset**：gap-fill band = 151（`api.rs:930` 的 100×1.5+1），且等长 gap 先走无 gap 直通路径（`extend.rs:1803-1832`），gapped DP 仅在直通分过低时触发；
- **bw_long < 段尺寸的用户自定义场景**（如 `-r 500,1000`）。

### 3.3 不使用 band 的路径

- **splice**：`dp::extend_splice` 完全忽略带宽（`dp/splice.rs:400, 1112, 1807` 注释 "Boundaries - NO bandwidth for splice"），全矩阵，靠 z-drop 截断；
- **chaining 阶段**：不是 DP，是对角差 gate（`chain_rmq.rs:709,723`）；
- **参考实现 `NWAligner`**（`align_simple.rs:1-25`）：无 band 的 O(n×m) NW，文件头对照表注明 prod aligner 为 "O(n × bandwidth)"。**未接入默认 pipeline**，仅文件内测试引用。

## 4. 上游对照：minimap2 的忠实移植

minimap2 master `align.c` `mm_align1()`（2026-11-17 核对）：

```c
bw      = (int)(opt->bw * 1.5 + 1.);
bw_long = (int)(opt->bw_long * 1.5 + 1.);
if (bw_long < bw) bw_long = bw;
...
/* gap filling */
int j, bw1 = bw_long, zdrop_code;
if (a[as1+i].y & MM_SEED_LONG_JOIN)
    bw1 = qe - qs > re - rs ? qe - qs : re - rs;   /* long-join 全宽 */
...
/* 左/右端部扩展 */
mm_align_pair(..., bw, opt->end_bonus, ...);
```

| 设计点 | minimap2 | bandmap | 一致性 |
|---|---|---|---|
| gap-fill band = `bw_long × 1.5 + 1` | ✓ | `extend.rs:1765` | ✓ |
| 端部扩展 band = `bw × 1.5 + 1` | ✓ | `extend.rs:1562` | ✓ |
| long-join 全宽 | ✓ | `extend.rs:1766-1768` | ✓ |
| `bw_long < bw` 守恒 | ✓ | `extend.rs:1564` | ✓ |
| RMQ rescue（`bw_long > bw` 时） | ✓（chain_rmq.c） | `map.rs:1047, 1373` | ✓ |
| SR ungapped 直通 | ✓ | `extend.rs:1803-1822` | ✓ |

**结论：`bw_long` 的存在本身来自上游**（`mm_mapopt_t.bw_long`），gap-fill band 失效的现象在 minimap2 同样存在，是上游作者的自觉取舍，不是移植走样。

## 5. 设计评估

### 5.1 站得住的部分：为什么 gap-fill 需要 40× band

链化的对角 gate 是**逐链边**的（每边 `|ref_diff − query_diff| ≤ 500`），但 **gap-fill 的 DP 段不等于一条链边**。fill 循环（`extend.rs:1727-1924`）存在三种"段累积"机制：

1. **锚点跳过**：`SEED_IGNORE | SEED_TANDEM` 的锚点被跳过且 `gap_rs/gap_qs` 不前进（`extend.rs:1738-1739`），前后段合并成一段——对角漂移是**多条链边的累加**，可远超单边的 500；
2. **小 gap 不触发 DP**：`should_align = is_last || is_long_join || (两维均 ≥ min_dp_len)`（`extend.rs:1759-1760`），未触发时 `gap_rs/gap_qs` 同样不前进（仅 `extend.rs:1923-1924` 在成功 DP 后前进）；
3. **坐标调整偏移**：HPC / k-mer 中点调整（`k_half`）使段两角的对角基线偏移。

因此 band 取 500 会裁剪合法比对；`bw_long=20000` 是"覆盖多边累积漂移的宽裕上限"。真正的成本与正确性兜底是 **z-drop 截断**（`extend.rs:1865-1879`）而非 band。这一分层（`bw`=局部连接，`bw_long`=长事件填充+救援）语义自洽，且经过 minimap2 十年实战检验。

### 5.2 欠账的部分

| # | 问题 | 证据 | 影响 |
|---|---|---|---|
| 1 | **语义过载**：`bw_long` 同时是 gap-fill DP band 和 RMQ rescue 链化带宽，两种算法共用一个旋钮 | `extend.rs:1765` vs `map.rs:1396` | 用户收紧 `-r 500,500` 想限 DP 成本，却因 rescue 触发条件 `bandwidth_long > bandwidth`（`map.rs:1373`）**静默关闭 rescue**；反之放宽 rescue 也无意识放大 gap-fill |
| 2 | **CLI 单值语义不对称** | `main.rs:577-578` | `-r 2000` 只改 `bandwidth`，`bandwidth_long` 保留 preset 值，用户以为收紧了 gap-fill 实际没有 |
| 3 | **复杂度承诺名不副实** | `align_simple.rs:11`（"prod aligner: O(n × bandwidth)"） | 该承诺只对端部扩展成立；gap-fill 在默认 lr 参数下是 O(max_gap²)。读代码者（含 `--bw` 文档）会期待它未提供的保证 |
| 4 | **魔法系数无文档** | `×1.5 + 1`（`extend.rs:1562-1563`） | 移植自上游但本地无注释说明其含义（对角漂移 + k-mer 偏移的余量） |

## 6. 改进建议

### P1：参数治理（低成本，建议先做）

- 将 `bandwidth_long` 拆为两个独立字段：`rescue_bandwidth`（chaining）与 `gap_fill_bandwidth`（alignment），preset 单独赋值；短期可保留 `bandwidth_long` 作为兼容别名同时驱动两者，但文档须写明联动。
- 修正 `-r` 单值语义：单值时同时设置两个参数（与 `-G` 行为一致），或在 help 中明确提示。
- 在 `extend.rs:1765` 与 `map.rs` rescue 处加注释，互相引用对方对 `bw_long` 的用途。

### P2：gap-fill band 封顶实验（需 benchmark，收益不确定）

map-ont 下尝试 `-r 500,5000`（gap-fill band = 7501 ≈ max_gap）：

- **理论依据**：单条链边的对角漂移 ≤ 500，封顶 5000 可无损覆盖任何单边 gap；只有"多边累积"的合并段可能被裁剪，而其中最极端者已被 `SEED_LONG_JOIN` 全宽路径（`extend.rs:1766-1768`）和 z-drop 分裂覆盖；
- **预期收益**：最长 gap 段的 DP 从 25M cells 降至 band×段长量级；
- **风险**：ON reads 的大 indel / 结构变异召回可能下降；
- **验证**：用含大 indel 的模拟 ONT 读集 + 真实数据对照 `dp_score` 分布与召回率（对照基准：当前默认 `-r 500,20000`）。

### 不建议改动

- splice / sr / ava 的 1× 设置（`bandwidth_long` 在这些场景已无独立作用）；
- `×1.5 + 1` 系数与 long-join 全宽规则（上游验证过的保守设计）；
- `NWAligner` 保留为参考实现，不建议接入生产路径（无 SIMD、无 z-drop、O(n×m)）。

## 7. 关键代码索引

| 主题 | 位置 |
|---|---|
| 参数定义与默认值 | `bandmap-core/src/align/map.rs:110-111, 205-206` |
| band 换算（×1.5+1） | `bandmap-core/src/align/extend.rs:1562-1568` |
| gap-fill 主循环（跳过/累积/`bw1`） | `bandmap-core/src/align/extend.rs:1727-1924`（跳过 1738-1739，`bw1` 1765-1768，DP 分发 1826-1877，前进 1923-1924） |
| 端部扩展 | `bandmap-core/src/align/extend.rs:1669-1724`（左）、2000 附近（右） |
| banded SIMD 内核 | `bandmap-core/src/align/dp/single.rs:58-65` |
| splice 无 band | `bandmap-core/src/align/dp/splice.rs:400, 1112, 1807` |
| chaining 对角 gate | `bandmap-core/src/align/chain_rmq.rs:589-613, 709, 723` |
| RMQ rescue | `bandmap-core/src/align/map.rs:1047-1071, 1373-1398` |
| preset 覆盖 | `bandmap-core/src/api.rs:891, 930, 946, 974, 984` |
| CLI | `bandmap/src/main.rs:82-84, 577-578, 619` |
| 无 band 参考实现 | `bandmap-core/src/align/align_simple.rs` |
| 既有文档 | `docs/map-ont-pipeline.md:35, 126, 193`；`docs/architecture.md:334, 622-623, 893` |

## 附录：证据链

- 本仓库结论均来自对 `bandmap-core` / `bandmap` 源码的直接阅读，行号以 2026-11-17 的 `main`（HEAD `168ea61`）为准；
- 上游对照基于 minimap2 master 分支 `align.c`（`https://raw.githubusercontent.com/lh3/minimap2/master/align.c`，2026-11-17 抓取），`mm_align1()` 中的 band 计算、gap-filling `bw1 = bw_long`、long-join 全宽、端部扩展 `bw` 均为原文引用；
- "fill 段累积多条链边"机制已在本仓库核实：`gap_rs/gap_qs` 仅在 `extend.rs:1923-1924` 成功 DP 后前进，锚点跳过（1738-1739）与未触发 `should_align`（1759-1760）均不前进。
