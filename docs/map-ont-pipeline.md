# map-ont 执行链路

本文梳理 `bandmap` 在 `map-ont` 预设下的一条完整执行链：从命令行 `bandmap -x map-ont ref.fa reads.fq` 或库调用 `Aligner::from_fasta(path, Preset::MapOnt)` 开始，到输出 SAM / PAF 记录结束。文中所有函数名与行号均以当前源码为准（工作区根目录 `/root/projects/rammap`）；与 [`architecture.md`](architecture.md) 中已过时的描述不一致之处在文末列出。

---

## 目录

1. [map-ont 的参数画像](#1-map-ont-的参数画像)
2. [总览：两条入口，一条主线](#2-总览两条入口一条主线)
3. [阶段一：选项装配（preset + finalize）](#3-阶段一选项装配preset--finalize)
4. [阶段二：参考索引构建 / 加载](#4-阶段二参考索引构建--加载)
5. [阶段三：查询调度（读线程 + rayon）](#5-阶段三查询调度读线程--rayon)
6. [阶段四：单条 read 的 mapping（map.rs）](#6-阶段四单条-read-的-mappingmaprs)
7. [阶段五：延伸与比对（pipeline.rs + extend.rs + dp/）](#7-阶段五延伸与比对pipeliners--extendrs--dp)
8. [阶段六：比对后过滤与 MAPQ](#8-阶段六比对后过滤与-mapq)
9. [阶段七：输出（SAM / PAF）](#9-阶段七输出sam--paf)
10. [map-ont 特性开关一览](#10-map-ont-特性开关一览)
11. [与 architecture.md 的差异](#11-与-architecturemd-的差异)

---

## 1. map-ont 的参数画像

`map-ont` 是 minimap2 的"默认长读"配置。关键点是：**它几乎没有专用参数**，绝大多数值来自 `MapOptions::default()`，预设本身只改 `k`、`w` 并重置 SIMD 上限。

| 项目 | 值 | 来源 |
|---|---|---|
| k / w | 15 / 10 | `api.rs:852-854`（`"lr" \| "map-ont"` 分支） |
| HPC | 否 | `is_hpc` 不修改（`api.rs:852-854`） |
| SIMD cap | `Auto`（允许 AVX-512） | `api.rs:850` |
| 计分 A/B/O/E/O2/E2 | 2/4/4/2/24/1 | `map.rs:218-230`（默认值） |
| transition / 模糊碱基罚分 | 0 / 1 | `map.rs:225-226` |
| 种子 | `mid_occ=0`（自动校准）、`max_occ=0`、`occ_dist=500`、`q_occ_frac=0.01` | `map.rs:189-197` |
| 链化 | `min_cnt=3`、`min_chain_score=40`、`max_gap=5000`、`bw=500`、`bw_long=20000`、`max_chain_skip=25`、`chain_gap_scale=0.8` | `map.rs:198-217` |
| 比对 | `zdrop=400`、`zdrop_inv=200`、`min_dp_max=80`、`min_dp_len=200` | `map.rs:231-241` |
| 过滤 | `best_n=5`、`pri_ratio=0.8`、`mask_level=0.5`、`seed=11` | `map.rs:242-253` |
| flags | `AlignFlags::empty()`（无 RMQ_CHAIN / SPLICE / SHORT_READ / HEAP_SORT / NO_DUAL 等） | `map.rs:259` |

`finalize_options()`（`api.rs:256-259`）在预设应用后重算链化罚分：`chn_pen_gap = 0.8 × 0.01 × k`（k=15 时为 0.12），`chn_pen_skip = chain_skip_scale × 0.01 × k`（为 0）。

> 注意：`k=15, w=10` 恰好也是 CLI 的初始默认值（`main.rs:562-563`），因此显式 `-x map-ont` 与不传预设得到的是同一套参数。`apply_preset_str()` 还会把 SIMD cap 重置为 `Auto`（`api.rs:850`）；由于全局默认值本身就是 `Auto`（`dp/mod.rs:45`），这一步对单次 CLI 运行没有影响，仅保证库调用方在同一进程内切换预设时行为可预测（`api.rs:846-849`）。

---

## 2. 总览：两条入口，一条主线

```
CLI 入口  bandmap/src/main.rs
  main() → Cli::parse() → run()
      └─ apply_preset() ──► api::apply_preset_str()
      └─ Index::build() / Index::load_part()
      └─ cal_mid_occ() 校准 mid_occ
      └─ map_one_part()
           ├─ 预读线程（sync_channel，批量读 FASTQ）
           └─ rayon par_iter → align_and_format_query()

库入口   bandmap-core/src/api.rs
  Aligner::from_fasta(path, Preset::MapOnt) → AlignerBuilder::resolve()/finish()
  Aligner::map_seq() → map_seq_with() → pipeline::process_query()

两条入口在 process_query() 汇合，共享同一条主线：
  process_query()
    → map_query()                       ① sketch → ② seed → ③ chain → ④ chains 过滤
    → process_query_core()
        → align_single_mapping()        ⑤ anchor 延伸 + 双仿射 DP + CIGAR
        → filter_alignment_results()    ⑥ 阈值过滤
        → assign_parents_and_select()   ⑦ primary / secondary 判定
        → compute_mapping_qualities()   ⑧ MAPQ
    → format_output()                   ⑨ SAM / PAF
```

---

## 3. 阶段一：选项装配（preset + finalize）

**CLI**（`bandmap/src/main.rs`）：

- `main()`（`main.rs:451-463`）：`Cli::parse()` → `run(args.align)`。注意 CLI 结构体用 `#[command(flatten)] align: AlignArgs`（`main.rs:32-33`），并非直接 `AlignArgs::parse()`。
- `run()` 中先给出默认值 `k=15, w=10, is_hpc=cli.is_hpc, opt = MapOptions::default()`（`main.rs:562-565`），若指定了 `-x` 则调用 `apply_preset()`（`main.rs:567-569`）。
- `apply_preset()`（`main.rs:513-515`）是薄封装，转发到 `bandmap::api::apply_preset_str()`（`api.rs:845-993`）。该函数先 `set_simd_cap(SimdCap::Auto)`（`api.rs:850`），`map-ont` 分支只设 `k=15, w=10`（`api.rs:852-854`）。
- 之后应用命令行显式覆盖项（`main.rs:571-742`），最后 `finalize_options(&mut opt, k)`（`main.rs:744`）。

**库**（`bandmap-core/src/api.rs`）：

- `Preset` 枚举定义在 `api.rs:39-56`，`MapOnt` 是其中一个变体，`as_str()` 返回 `"map-ont"`。
- `Aligner::from_fasta(path, Preset::MapOnt)` → `AlignerBuilder::from_fasta()`（`api.rs:347-350`）→ `resolve()`（内部调用 `apply_preset_str`，`api.rs:400`）→ 读 FASTA、`Index::build()` → `finish()`（`api.rs:427-433`）。
- `finish()` 同样做 `finalize_options` 与 `mid_occ` 校准，并在库默认配置下置 `OUT_CIGAR`（库默认 `cigar=true`，`api.rs:314`；这与 CLI 不同——CLI 只有指定 `-c/-a/--cs/--md/-ds/--write-junc` 时才输出 CIGAR，`main.rs:525`）。
- `Aligner::map_seq()` → `map_seq_with()`（`api.rs:536-561`）创建 `AlignmentContext` / `MapContext`，调用 `pipeline::process_query()`（`api.rs:550-553`）。

---

## 4. 阶段二：参考索引构建 / 加载

- 预构建索引：`Index::load_part()`（`main.rs:845`）；库侧 `Index::load` 支持 RMMI / MMI v2 等格式。
- FASTA 构建：`Index::build(all_seqs, w, k, is_hpc, usize::MAX)`（`main.rs:925`；流式/多分卷路径见 `main.rs:982, 1063`），定义在 `index.rs:583`。流程为：packed 序列 → `sketch_sequence` 生成 minimizer → `radix_sort_pair` 按哈希排序 → 装入 per-bucket 开放寻址哈希表（`index_bucket.rs`）。
- `mid_occ` 校准：`if opt.seeding.mid_occ == 0 { opt.seeding.mid_occ = mi.cal_mid_occ(...) }`（`main.rs:1277-1279`，默认比例 `2e-4`，`main.rs:809`；库侧 `api.rs:263-268`）。`map-ont` 的 `mid_occ=0`，因此该步骤必然执行。

---

## 5. 阶段三：查询调度（读线程 + rayon）

`map_one_part()`（`main.rs:1258-1600`）是单个索引分卷上的查询驱动：

1. 校准 `mid_occ`（见上）。
2. 开启并行时，用一个预读线程把 FASTQ 记录攒成 `Vec<SeData>`，经 `sync_channel` 送给工作线程；channel 容量 1，内存占用上限约两批（`main.rs:1483-1511`）。另有写线程负责 stdout（`main.rs:1513-1524`）。
3. rayon 对每批 `chunk_data.par_iter().map_init(...)`（`main.rs:1529-1546`），每个线程独立持有一个 `AlignmentContext` 和 `MapContext`（通过 `map_init` 复用缓冲区），调用 `bandmap::align::pipeline::align_and_format_query()`（`main.rs:1534-1537`）；单条 read 用 `catch_unwind` 包裹，避免一条坏 read 中断整个批次。
4. 无 `parallel` feature 时走串行回退（`main.rs:1555-1596`）。

`align_and_format_query()`（`pipeline.rs:2098-2126`）是 CLI 与后续管线的衔接点：`process_query()` → （splice 模式下 jump 延伸，`map-ont` 不触发）→ `format_output()`。

---

## 6. 阶段四：单条 read 的 mapping（map.rs）

`process_query()`（`pipeline.rs:1160-1180`）先做 `max_qlen` 检查，然后 `map_query()`（`pipeline.rs:1177`）得到 `(regs, rep_len, stats, squeezed)`，再交给 `process_query_core()`（`pipeline.rs:1179`）。

`map_query()`（`map.rs:1285-1765`）的内部顺序如下（这是整条链的核心）：

1. **读名哈希**：`compute_read_hash()`（`map.rs:1295`），用于确定性 seed 与排序 tie-break。
2. **sketch**：`sketch_sequence()`（`map.rs:1299`；实现 `sketch.rs:174`），滚动哈希 + 窗口最小化，产出 `Vec<Minimizer>` 到 `ctx.minimizers`。
3. **高频 query minimizer 过滤**：`filter_minimizers_by_occ()`（`map.rs:1303`；实现 `seed.rs:615`），条件是 `q_occ_frac>0 && mid_occ>0 && minimizers.len()>mid_occ`。
4. **种子收集**：非 HEAP_SORT 路径走 `collect_seed_hits()`（`map.rs:1316`；实现 `seed.rs:432-443`）→ `collect_seed_hits_with_occ()`（`seed.rs:316-429`），按哈希查索引，产出 anchor，并在 `seed.rs:427` 调用 `radix_sort_128x()` 按 `(ref_id, ref_pos)` 排序。`map-ont` 无 `HEAP_SORT` 标志，不走 heap 模式。
5. **对角过滤**：`filter_anchors_diagonal()`（`map.rs:1322`），仅 `--diag-filter` 时生效，默认 no-op。
6. **主链化**：`RMQ_CHAIN` 未设置，进入 else 分支（`map.rs:1344-1360`）；`PAR_CHAIN` 也未设置，因此调用 `chain_anchors()`（`map.rs:1356`；实现 `chain.rs:224-278`）。`chain_anchors()` 在 anchor 数 ≥ 32 且非 cdna / 单段时走 SIMD 分派（`chain.rs:234-255`）：AVX-512 → AVX2 → SSE；因为 SIMD cap 为 `Auto`，AVX-512 在支持的 CPU 上可用。低于 32 个 anchor 回退标量实现。
7. **RMQ 救援重链化**：`map-ont` 默认也会进入（`map.rs:1369-1405`）。条件为 `bandwidth_long(20000) > bandwidth(500)`、链数 > 1、且未设置 `SPLICE|SHORT_READ|NO_LJOIN`；若覆盖度不足或过长则重新 `radix_sort_128x` 并调用 `chain_anchors_rmq()`（`map.rs:1396-1399`，使用 `bandwidth_long`）。注意：这是"救援"阶段，`map-ont` 的主链化仍是普通 DP/SIMD 链化。
8. **高频 seed 重链化**：`map.rs:1407-1429`。`map-ont` 默认 `max_occ=0`，该分支实际不触发。
9. **构造 `Vec<Mapping>`**：`map.rs:1434-1481`，按 `min_chain_score` 过滤（`map.rs:1484-1486`）、链过滤按链得分与读名哈希排序（`map.rs:1496-1511`）。
10. **链后过滤（chain_post 内联块）**：`map.rs:1517-1764`。依次做：
    - `ParentState::assign_parents()`（`map.rs:1622`；实现 `filter.rs:99`），按 query 重叠度把链归到 primary 或 secondary；
    - 预比对 subscore / n_sub 计算与两遍 `check_secondary_filter()` 过滤（`map.rs:1624-1696`）、紧凑化（`map.rs:1698-1714`）；
    - `estimate_divergence()` 估算 `dv:f`（`map.rs:1720`）；
    - 构造 squeezed anchor 数组并计算附近种子边界（`map.rs:1734-1760`）。

    返回值：`(filtered_regs, rep_len, stats, squeezed)`（`map.rs:1764`）。

> 与 `architecture.md` 不同：源码中并不存在名为 `chain_post` / `select_sub` 的独立函数，相关逻辑直接内联在 `map_query()` 尾部。

---

## 7. 阶段五：延伸与比对（pipeline.rs + extend.rs + dp/）

`process_query_core()`（`pipeline.rs:1660-1860`）接收上一步的结果：

1. 重排 regs（`pipeline.rs:1683-1697`），一次性把 query 的正向与反向互补序列编码成 nt4（`pipeline.rs:1727-1728`）。
2. 进入工作队列循环（`pipeline.rs:1730-1833`）。`map-ont` 不是 splice/dual-splice 模式，走 `else` 分支（`pipeline.rs:1788-1796`）调用 `align_single_mapping()`（`pipeline.rs:1790`；实现 `pipeline.rs:427`）。
3. `align_single_mapping()` 计算延伸窗口（坏端修剪、k-mer 中心偏移，`pipeline.rs:464-585`），在 `pipeline.rs:742` 与 `751` 调用 `align_anchors()`（实现 `extend.rs:1514`）。
4. `align_anchors()` 完成左延伸 → 间隙填充 → 右延伸，生成原始 CIGAR。`map-ont` 的 DP 分派逻辑在 `extend.rs:1709-1715`（左延伸）、`1825-1878`（间隙填充）、`1999-2005`（右延伸）：
   - `is_splice=false`，不调用 `extend_splice`；
   - `q=4, q2=24`（`q != q2`），因此调用 **`dp::extend_dual_affine()`**（`extend.rs:1714, 1830, 1845, 1861, 1877, 2004`），即双仿射 gap 模型；
   - SIMD 由 `dp/dual.rs:111-134` 分派：AVX-512 → AVX2 → SSE4.1 → 标量。`map-ont` 的 `SimdCap::Auto` 允许 AVX-512。
5. **z-drop**：`test_zdrop()`（`extend.rs:694-797`）在 `extend.rs:1834-1835` 与 `1866-1867` 被调用；触发后用 `zdrop2` 重跑 DP 并标记 `split_right_anchors`（`extend.rs:1907-1908`）与 `split_inv`（`extend.rs:1910-1911`）。
6. **CIGAR 统计**：DP kernel 通过 `push_cigar` 生成 CIGAR；`CigarStats::from_cigar()`（`pipeline.rs:149-234`）在 `pipeline.rs:873` / `1096` 被调用，产出 matches、NM、block_len、divergence 等。
7. 若 z-drop 产生 split，循环回到队列处理剩余片段；`NO_INV` 未设置，因此还会做倒位检测（`pipeline.rs:1821-1832`）。

---

## 8. 阶段六：比对后过滤与 MAPQ

仍在主线 `process_query_core()` 尾部（`pipeline.rs:1835-1859`）：

1. `filter_alignment_results()`（`pipeline.rs:1340-1393`；调用点 `pipeline.rs:1839`，仅 `do_cigar` 时）按 `min_cnt`、`min_chain_score`、`min_dp_max`、`max_clip_ratio` 做阈值过滤。
2. `assign_parents_and_select()`（`pipeline.rs:1398-`；调用点 `pipeline.rs:1843`）做 primary/secondary 判定、secondary 选择（保留 `best_n`、比值 `pri_ratio`）与 `dp_max` 排名。非 CIGAR 路径会保留阶段四中 chain_post 的父子关系（`pipeline.rs:1421-1438`）。
3. 设置 SAM primary 标志（`pipeline.rs:1845-1855`）。
4. `compute_mapping_qualities()`（`pipeline.rs:1228-`；调用点 `pipeline.rs:1857`）按最优/次优 `dp_max` 的对数似然比计算 MAPQ，clamp 到 `[0,60]`。
5. 返回 `ProcessedQuery { results, mapqs, sam_pri, parent_indices, rep_len, stats }`（`pipeline.rs:1859`）。

> 注意：函数名是 `compute_mapping_qualities`，而非 `architecture.md` 中的 `compute_mapqs`。

---

## 9. 阶段七：输出（SAM / PAF）

`format_output()`（`pipeline.rs:2702-2773`）按结果逐条输出：

- 无命中时输出 unmapped 记录（`pipeline.rs:2734-2736`）；
- `-a` / SAM 模式走 `format_sam_record()`（`pipeline.rs:2413-2621`；调用点 `2759`）；否则走 `format_paf_record()`（`pipeline.rs:2624-2697`；调用点 `2766`）。

**SAM 记录**（`pipeline.rs:2440-2618`）：FLAG（`:2440-2458`）、带 soft clip 的 CIGAR（`:2463-2496`）、11 个主字段（`:2535-2570`），随后是可选标签：`NM:i`（`:2572`）、`ms:i`（`:2573`）、`AS:i`（`:2574`）、`nn:i`（`:2575`）、`tp:A`（`:2581`）、`cm:i`（`:2582`）、`s1:i`（`:2583`）、`s2:i`（`:2584`）、`de:f`（`:2587`）、`zd:i`（`:2592`）、`SA:Z`（`:2596`）、`RG:Z`（`:2600`）、`MD/ds/cs`（`:2603`）、`rl:i`（`:2618`）等。splice 专有标签 `ts:A` 在 `:2576-2580`，`map-ont` 不输出。

**PAF 记录**（`pipeline.rs:2648-2694`）：12 个必需列（`:2648-2651`），CIGAR 模式下附 `NM/ms/AS/nn`（`:2653-2663`）、`tp/cm/s1/s2`（`:2664-2669`）、`de:f` 或 `dv:f`（`:2670-2682`）、`zd`（`:2683-2685`）、`rl`（`:2687`）、`cg:Z`（`:2688-2690`）、`MD/ds/cs`（`:2691-2694`）。

`align_and_format_query()` 把输出累积到字符串返回；CLI 侧的写线程负责最终写到 stdout（`main.rs:1517-1519`）。

---

## 10. map-ont 特性开关一览

| 特性 | 是否启用 | 依据 |
|---|---|---|
| 主链化：DP / SIMD 链化 | **是** | `RMQ_CHAIN` 未设置；`map.rs:1344-1360`，`chain.rs:234-255` |
| 主链化：RMQ | 否 | 仅 `lr:hqae` / `asm*` 设置 `RMQ_CHAIN` |
| RMQ 救援重链化 | **是（条件触发）** | `bw_long > bw`、链数 > 1、无 `SPLICE\|SHORT_READ\|NO_LJOIN`，`map.rs:1372-1373` |
| SIMD 链化（AVX-512/AVX2/SSE） | **是**（anchor ≥ 32） | `chain.rs:238-255`，cap `Auto`（`api.rs:850`） |
| 双仿射 DP | **是** | `q=4 ≠ q2=24`，`extend.rs:1714` 等；`NO_DUAL` 未设置 |
| 单仿射 DP | 否 | 仅 `q==q2 && e==e2` 时选，`extend.rs:1711-1712` |
| 剪接感知 DP / jump 延伸 | 否 | `SPLICE` 未设置；`extend.rs:1709`、`pipeline.rs:2113-2120` |
| 短读优化 / HEAP_SORT / ALL_CHAINS | 否 | `flags = empty()`（`map.rs:259`） |
| AVX-512 DP | 允许 | `SimdCap::Auto`（`api.rs:850`） |
| HPC 种子 | 否 | `is_hpc=false`（`api.rs:852-854`） |
| 倒位检测 | **是** | `NO_INV` 未设置，`pipeline.rs:1821-1832` |
| 父子链 / secondary 过滤 | **是** | `map.rs:1622`、`pipeline.rs:1843` |

---

## 11. 与 architecture.md 的差异

`architecture.md` 描述的是整体架构，其中若干细节与当前源码不符，梳理 `map-ont` 链路时需注意：

| `architecture.md` 的说法 | 实际源码 |
|---|---|
| `chain_post()` 是 `map_query()` 调用的独立函数（`architecture.md:298, 307`） | 不存在该函数；逻辑内联在 `map_query()` 尾部（`map.rs:1517-1764`） |
| `select_sub()` 参与链后过滤（`:298`） | 不存在该函数；实际是两遍 `check_secondary_filter()`（`map.rs:1656-1696`） |
| RMQ 救援在 chain_post **之后**（`:308`） | RMQ 救援在 `map.rs:1369-1405`，早于链后过滤块（1517+） |
| `chain_simd.rs` 仅 "AVX2 + NEON + WASM SIMD128"（`:48`） | 另有 AVX-512 与 SSE2 实现 |
| 模块抽象为单个 `dp.rs`（`:162`） | 已拆分为 `dp/{single,dual,splice,lw,common}.rs`，分派在 `dp/dual.rs:111-134` |
| MAPQ 函数名为 `compute_mapqs` | 实际为 `compute_mapping_qualities`（`pipeline.rs:1228`） |
| CLI `AlignArgs::parse()` | 实际为 `Cli::parse()` + `#[command(flatten)]`（`main.rs:32-33, 452`） |
