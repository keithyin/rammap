#!/usr/bin/env python3
"""对角线直方图 top-K anchor 过滤原型（覆盖 anchor 构建 → chaining）。

管线五步（挂在 anchor 构建之后、DP 之前）：
  ① 对每个 anchor 算 d = r - q（query 已按链方向镜像，正/负链公式相同）
  ② 分桶计数 b = d // B
  ③ 每条链位取计数 top-K 的桶（或按覆盖率自适应）
  ④ 只保留入选桶（±neighbor 邻桶）内的 anchor
  ⑤ 照常 radix sort → DP 链化

DP chaining 忠实移植 bandmap-core/src/align/chain.rs 的 scalar 版本：
compute_chain_score、visited/skip-break、backtrack(max_drop)。
省略了 best_anchor_idx 二次检查（原型只关心过滤效果，不影响结论方向）。

运行：
    python3 prototype/diag_topk.py                    # 跑默认实验（topk, K=2）
    python3 prototype/diag_topk.py --K 1              # 演示 K=1 的失败模式
    python3 prototype/diag_topk.py --mode adaptive    # 覆盖率自适应选桶
    python3 prototype/diag_topk.py --fallback         # 开启"过滤后分数骤降→重链化"兜底
"""

from __future__ import annotations

import argparse
import random
import zlib
from collections import Counter
from dataclasses import dataclass
from math import log2

RC = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]


# ---------------------------------------------------------------- anchors ----

@dataclass
class Anchor:
    r: int        # target(ref) 起点
    q: int        # query 起点，负链时已镜像: q = qlen - (pos_rc + k)
    strand: int   # +1 / -1
    span: int     # k


def minimizers(seq: str, k: int, w: int) -> list[tuple[int, int]]:
    """(hash, pos)；长度为 w 的窗口内取 hash 最小的 k-mer，平局取最左。

    同一 minimizer 在连续窗口重复胜出时只发射一次（mm2 行为）。
    """
    if len(seq) < k:
        return []
    kmers = [zlib.crc32(seq[i:i + k].encode()) for i in range(len(seq) - k + 1)]
    out = []
    last_pos = -1
    for i in range(0, len(kmers) - w + 1):
        win = kmers[i:i + w]
        h, p = min((h, -j) for j, h in enumerate(win))
        pos = i - p  # 平局时选窗口内最靠左的 k-mer
        if pos != last_pos:
            out.append((h, pos))
            last_pos = pos
    return out


def build_target_index(target: str, k: int, w: int) -> dict[int, list[int]]:
    idx: dict[int, list[int]] = {}
    for h, p in minimizers(target, k, w):
        idx.setdefault(h, []).append(p)
    return idx


def build_anchors(read: str, tindex: dict[int, list[int]], k: int, w: int) -> list[Anchor]:
    """query 位置按链方向镜像，使 d = r - q 对两条链语义一致。"""
    qlen = len(read)
    anchors: list[Anchor] = []
    for strand, seq in ((+1, read), (-1, revcomp(read))):
        for h, p in minimizers(seq, k, w):
            for tp in tindex.get(h, ()):
                if strand > 0:
                    q = p
                else:
                    q = qlen - (p + k)  # 镜像回正向坐标
                anchors.append(Anchor(r=tp, q=q, strand=strand, span=k))
    return anchors


# ------------------------------------------------- diagonal histogram filter --

def diagonal_filter(anchors: list[Anchor], B: int = 50, K: int = 2,
                    neighbor: int = 1, mode: str = "topk",
                    coverage: float = 0.8, min_bucket: int = 2) -> list[Anchor]:
    """按链位分组 → 直方图选峰 → 只留入选桶 ±neighbor 内的 anchor。"""
    out: list[Anchor] = []
    for strand in (+1, -1):
        group = [a for a in anchors if a.strand == strand]
        if not group:
            continue
        counts = Counter((a.r - a.q) // B for a in group)
        if mode == "adaptive":
            # 按计数降序累加，直到覆盖该链位 ~coverage 或计数低于阈值
            chosen = []
            cum = 0
            for b, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
                chosen.append(b)
                cum += c
                if cum >= coverage * len(group) or c < min_bucket:
                    break
        else:
            chosen = [b for b, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:K]]
        sel = set()
        for b in chosen:
            sel.update(range(b - neighbor, b + neighbor + 1))
        out.extend(a for a in group if (a.r - a.q) // B in sel)
    return out


# --------------------------------------------------------------- DP chaining --

@dataclass
class Chain:
    score: int
    n_anchor: int
    anchors: list[Anchor]  # 链上 anchor（按 ref 升序）


@dataclass
class ChainingParams:
    min_cnt: int = 3
    min_chain_score: int = 40
    max_dist_x: int = 500
    max_dist_y: int = 500
    bandwidth: int = 500
    max_chain_skip: int = 25
    max_chain_iter: int = 50
    chn_pen_gap: float = 0.5
    chn_pen_skip: float = 0.05


def _chain_score(ai: Anchor, aj: Anchor, p: ChainingParams) -> int:
    """移植 chain.rs::compute_chain_score（is_cdna=False, n_seg=1）。"""
    qd = ai.q - aj.q
    if qd <= 0 or qd > p.max_dist_x:
        return -(1 << 30)
    rd = ai.r - aj.r
    if rd == 0 or qd > p.max_dist_y:
        return -(1 << 30)
    gap = rd - qd if rd > qd else qd - rd
    if gap > p.bandwidth:
        return -(1 << 30)
    min_d = rd if rd < qd else qd
    span = aj.span
    sc = span if span < min_d else min_d
    if gap > 0 or min_d > span:
        lin = p.chn_pen_gap * gap + p.chn_pen_skip * min_d
        log_pen = log2(gap + 1) if gap >= 1 else 0.0
        sc -= int(lin + 0.5 * log_pen)
    return sc


@dataclass
class ChainStats:
    dp_scans: int = 0  # 前向 DP 内层 j 循环总步数


def chain_anchors(anchors: list[Anchor], p: ChainingParams,
                  stats: ChainStats | None = None) -> list[Chain]:
    """按链位分组后各自 DP（等价于 chain_anchors_partitioned 的单分区路径）。"""
    chains: list[Chain] = []
    for strand in (+1, -1):
        a = sorted((x for x in anchors if x.strand == strand), key=lambda x: x.r)
        if a:
            chains.extend(_chain_partition(a, p, stats))
    return chains


def _chain_partition(a: list[Anchor], p: ChainingParams,
                     stats: ChainStats | None) -> list[Chain]:
    n = len(a)
    score = [0] * n
    pred = [-1] * n
    visited = [0] * n
    ws = 0
    for i in range(n):
        best = a[i].span
        best_pred = -1
        skip = 0
        while ws < i and a[i].r - a[ws].r > p.max_dist_x:
            ws += 1
        if i - ws > p.max_chain_iter:
            ws = i - p.max_chain_iter
        end_j = ws - 1
        for j in range(i - 1, ws - 1, -1):
            if stats is not None:
                stats.dp_scans += 1
            sc = _chain_score(a[i], a[j], p)
            if sc == -(1 << 30):
                continue
            tot = sc + score[j]
            if tot > best:
                best = tot
                best_pred = j
                if skip > 0:
                    skip -= 1
            elif visited[j] == i:  # 前驱已被更近的 j 覆盖 → skip-break
                skip += 1
                if skip > p.max_chain_skip:
                    end_j = j
                    break
            if pred[j] >= 0:
                visited[pred[j]] = i
        score[i] = best
        pred[i] = best_pred

    # backtrack（移植 chain_backtrack / chain_backtrack_end，max_drop = bw）
    cands = sorted(((s, i) for i, s in enumerate(score) if s >= p.min_chain_score),
                   key=lambda t: (-t[0], t[1]))
    chains: list[Chain] = []
    v2 = [0] * n
    for _, k in cands:
        if v2[k]:
            continue
        # chain_backtrack_end：带 max_drop 的回溯找链终点 max_i
        i = k
        max_i = k
        max_s = 0
        while True:
            v2[i] = 2
            i = pred[i]
            s = score[k] if i < 0 else score[k] - score[i]
            if s > max_s:
                max_s = s
                max_i = i
            elif max_s - s > p.bandwidth:
                break
            if i < 0 or v2[i]:
                break
        # 收集 k → max_i 之间的 pred 链（不含 max_i），链分 = score[k]-score[max_i]
        idxs: list[int] = []
        i = k
        while i >= 0 and i != max_i:
            idxs.append(i)
            v2[i] = 1
            i = pred[i]
        sc = score[k] if max_i < 0 else score[k] - score[max_i]
        if sc >= p.min_chain_score and len(idxs) >= p.min_cnt:
            anchors = [a[j] for j in reversed(idxs)]
            chains.append(Chain(score=sc, n_anchor=len(idxs), anchors=anchors))
    chains.sort(key=lambda c: -c.score)
    return chains


# ------------------------------------------------------------ read 映射管线 --

def map_read(read: str, tindex: dict[int, list[int]], k: int, w: int,
             filt: dict | None, p: ChainingParams, fallback: bool = False) -> dict:
    """anchor 构建 → 对角线过滤（filt=None 则跳过）→ 链化。返回统计与最优链。"""
    anchors = build_anchors(read, tindex, k, w)
    st_all = ChainStats()
    chains_all = chain_anchors(anchors, p, st_all)
    best_all = chains_all[0].score if chains_all else 0

    if filt is None:  # baseline：不过滤
        chains, n_scan, n_anchor = chains_all, st_all.dp_scans, len(anchors)
        used_fallback, kept = False, anchors
    else:
        kept = diagonal_filter(anchors, **filt)
        st_f = ChainStats()
        chains_f = chain_anchors(kept, p, st_f)
        best_f = chains_f[0].score if chains_f else 0
        chains, n_scan, n_anchor = chains_f, st_f.dp_scans, len(kept)
        used_fallback = False
        if fallback and chains_all and best_f < 0.6 * best_all:
            # 兜底：过滤后最优链分数骤降 → 用未过滤 anchor 重链化
            chains, n_scan, n_anchor = chains_all, st_all.dp_scans, len(anchors)
            used_fallback = True

    return {
        "n_anchor": len(anchors), "n_kept": len(kept),
        "dp_scans": n_scan, "best": chains[0] if chains else None,
        "fallback": used_fallback,
    }


# ------------------------------------------------------------ 合成数据实验 --

def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = []
    for b in seq:
        if rng.random() < rate:
            out.append(rng.choice([c for c in "ACGT" if c != b]))
        else:
            out.append(b)
    return "".join(out)


def gen_target(rng: random.Random, size: int = 100_000) -> str:
    """随机背景 + 6 个散置重复家族（8 拷贝、拷贝间 3% 歧化 → 多对角线噪声）
    + 1 个串联重复（occ≈8 精确重复 → 次优链/MAPQ 场景）。"""
    seq = [rng.choice("ACGT") for _ in range(size)]
    for _ in range(10):
        unit = "".join(rng.choice("ACGT") for _ in range(rng.randint(500, 700)))
        loci = rng.sample(range(size - len(unit)), 8)
        for p in loci:
            seq[p:p + len(unit)] = list(_mutate(unit, 0.03, rng))
    unit = "".join(rng.choice("ACGT") for _ in range(300))
    p = rng.randrange(size - 8 * len(unit))
    seq[p:p + 8 * len(unit)] = list(unit * 8)
    return "".join(seq)


def sample_read(target: str, rng: random.Random, length: int = 1700,
                sub: float = 0.08, ins: float = 0.02, dele: float = 0.02,
                big_del_prob: float = 0.4) -> tuple[str, tuple[int, int]]:
    start = rng.randrange(0, len(target) - length)
    end = start + length
    seg = list(target[start:end])
    if rng.random() < big_del_prob:  # 100-200bp 缺失 → 对角线跳变跨桶
        dp = rng.randrange(200, len(seg) - 250)
        seg[dp:dp + rng.randint(100, 200)] = []
    out = []
    for b in seg:
        u = rng.random()
        if u < sub:
            out.append(rng.choice([c for c in "ACGT" if c != b]))
        elif u < sub + ins:
            out.append(rng.choice("ACGT"))
            out.append(b)
        elif u < sub + ins + dele:
            continue
        else:
            out.append(b)
    return "".join(out), (start, end)


def truth_recall(chain: Chain | None, truth: tuple[int, int], k: int) -> float:
    if chain is None or not chain.anchors:
        return 0.0
    r0 = min(a.r for a in chain.anchors)
    r1 = max(a.r + k for a in chain.anchors)
    ov = min(r1, truth[1]) - max(r0, truth[0])
    return max(0.0, ov) / (truth[1] - truth[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", type=int, default=60)
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--w", type=int, default=5)
    ap.add_argument("--B", type=int, default=100)
    ap.add_argument("--K", type=int, default=2)
    ap.add_argument("--neighbor", type=int, default=1)
    ap.add_argument("--mode", choices=["topk", "adaptive"], default="topk")
    ap.add_argument("--fallback", action="store_true",
                    help="开启'过滤后最优链分数 < 60%% 原值 → 重链化'兜底")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    target = gen_target(rng)
    tindex = build_target_index(target, args.k, args.w)
    reads_truth = [sample_read(target, rng) for _ in range(args.reads)]

    filt = {"B": args.B, "K": args.K, "neighbor": args.neighbor, "mode": args.mode}
    p = ChainingParams()

    configs = [("baseline", None)] + [
        (f"K={args.K}({args.mode})", filt),
    ]
    print(f"target={len(target)}bp  reads={args.reads}  k={args.k} w={args.w} "
          f"B={args.B} neighbor=±{args.neighbor} fallback={args.fallback}")
    print(f"{'config':<16} {'anchors':>8} {'kept':>7} {'dp_scans':>10} "
          f"{'scan/anchor':>11} {'best score':>10} {'best n':>7} {'recall':>7}")
    for name, f in configs:
        ta = tk = ts = 0
        nbest = nscore = nreads = 0
        recall_ok = 0
        for read, truth in reads_truth:
            r = map_read(read, tindex, args.k, args.w, f, p, args.fallback)
            ta += r["n_anchor"]
            tk += r["n_kept"] if f else r["n_anchor"]
            ts += r["dp_scans"]
            if r["best"]:
                nbest += r["best"].n_anchor
                nscore += r["best"].score
                nreads += 1
                if truth_recall(r["best"], truth, args.k) >= 0.9:
                    recall_ok += 1
        print(f"{name:<16} {ta:>8} {tk:>7} {ts:>10} "
              f"{ts / ta:>11.1f} {nscore / max(1, nreads):>10.0f} "
              f"{nbest / max(1, nreads):>7.1f} {recall_ok}/{len(reads_truth)}")

    # 单 read 明细：最优链的 score / anchor 数对照
    read, truth = reads_truth[0]
    for name, f in configs:
        r = map_read(read, tindex, args.k, args.w, f, p, args.fallback)
        b = r["best"]
        s = f"{b.score}/{b.n_anchor}" if b else "none"
        print(f"read#0 {name:<12} anchors {r['n_anchor']}→{r['n_kept'] or '-'} "
              f"dp_scans {r['dp_scans']:>6} best {s} recall "
              f"{truth_recall(b, truth, args.k):.2f}")


if __name__ == "__main__":
    main()
