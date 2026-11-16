//! Diagonal histogram top-K anchor filtering.
//!
//! Ported from the prototype `prototype/anchor_filter/diag_topk.py`
//! (`diagonal_filter`, topk mode). Inserted between seed collection and
//! DP chaining to curb anchor blow-up caused by repetitive/low-complexity
//! regions: only anchors on the top-K most populated diagonals (per
//! reference + strand group) are kept.
//!
//! Pipeline steps (per `ref_id + strand` group):
//!   1. compute `d = ref_pos - query_pos` for every anchor (query positions
//!      are already mirrored for the reverse strand at seed-collection time,
//!      so the formula is identical for both strands)
//!   2. bucket count: `b = d.div_euclid(B)` (floor semantics, matching
//!      Python `//` for negative diagonals)
//!   3. take the K buckets with the highest counts (ties broken by smaller
//!      bucket index, matching the prototype)
//!   4. keep only anchors whose bucket is in the chosen set expanded by
//!      ±`neighbor` neighboring buckets
//!
//! Anchors must be sorted by `x` (as produced by `collect_seed_hits` /
//! `collect_seed_hits_heap`); grouping relies on `ref_id_strand()` runs being
//! contiguous and the filter preserves relative order within each group.

use std::collections::HashMap;

use crate::align::sketch::Minimizer;

/// Parameters of the diagonal histogram top-K anchor filter.
///
/// Defaults mirror the prototype: `B = 100`, `K = 2`, `neighbor = +-1`,
/// disabled unless explicitly turned on (CLI `--diag-filter`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagFilterParams {
    /// Master switch; filtering is a no-op when false.
    pub enabled: bool,
    /// Diagonal bucket width `B` (in bases).
    pub bucket_size: i64,
    /// Number of top buckets `K` to keep per group.
    pub topk: usize,
    /// Number of neighboring buckets kept around each chosen bucket.
    pub neighbor: i64,
}

impl Default for DiagFilterParams {
    fn default() -> Self {
        DiagFilterParams {
            enabled: false,
            bucket_size: 100,
            topk: 2,
            neighbor: 1,
        }
    }
}

/// Diagonal value `d = ref_pos - query_pos` of an anchor.
#[inline]
fn diagonal(a: &Minimizer) -> i64 {
    a.ref_pos() as i64 - a.query_pos() as i64
}

/// Expand chosen buckets by +-`neighbor` into a sorted, deduplicated list.
fn expand_selected(chosen: &[i64], neighbor: i64) -> Vec<i64> {
    if neighbor <= 0 {
        let mut v = chosen.to_vec();
        v.sort_unstable();
        v.dedup();
        return v;
    }
    let mut v = Vec::with_capacity(chosen.len() * (2 * neighbor as usize + 1));
    for &b in chosen {
        for k in (b - neighbor)..=(b + neighbor) {
            v.push(k);
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

/// Apply the diagonal histogram top-K filter in place.
///
/// Anchors are grouped by contiguous `ref_id_strand()` runs (the input is
/// sorted by `x`), each group selects its own top-K buckets, and anchors
/// outside the selected buckets are removed. Relative order is preserved so
/// the output stays sorted by `x` and remains directly usable by chaining.
///
/// Returns the number of anchors kept.
pub fn filter_anchors_diagonal(anchors: &mut Vec<Minimizer>, p: &DiagFilterParams) -> usize {
    if !p.enabled || p.topk == 0 || p.bucket_size <= 0 {
        return anchors.len();
    }
    let b = p.bucket_size;
    let n = anchors.len();
    let mut counts: HashMap<i64, u32> = HashMap::new();
    let mut kept = 0usize;
    let mut start = 0usize;
    while start < n {
        let key = anchors[start].ref_id_strand();
        let mut end = start + 1;
        while end < n && anchors[end].ref_id_strand() == key {
            end += 1;
        }

        // 1-2: histogram over diagonal buckets.
        counts.clear();
        for a in &anchors[start..end] {
            *counts.entry(diagonal(a).div_euclid(b)).or_insert(0u32) += 1;
        }

        // 3: top-K buckets, count descending, bucket index ascending on ties.
        let mut ranked: Vec<(i64, u32)> = counts.drain().collect();
        ranked.sort_unstable_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
        ranked.truncate(p.topk);
        let chosen: Vec<i64> = ranked.iter().map(|&(bk, _)| bk).collect();

        // 4: keep anchors inside chosen buckets expanded by +-neighbor.
        let selected = expand_selected(&chosen, p.neighbor);
        for i in start..end {
            let bucket = diagonal(&anchors[i]).div_euclid(b);
            if selected.binary_search(&bucket).is_ok() {
                if kept != i {
                    anchors[kept] = anchors[i];
                }
                kept += 1;
            }
        }
        start = end;
    }
    anchors.truncate(kept);
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an anchor with a given ref_id, strand, ref/query position.
    /// Bit layout matches seed collection: x = strand << 63 | rid << 32 | ref_pos.
    fn mk(rid: u64, strand: u64, ref_pos: i32, query_pos: i32) -> Minimizer {
        let x = (strand << 63) | (rid << 32) | (ref_pos as u32 as u64);
        let y = (5u64 << 32) | (query_pos as u32 as u64); // span = 5
        Minimizer { x, y }
    }

    fn params(topk: usize, neighbor: i64) -> DiagFilterParams {
        DiagFilterParams {
            enabled: true,
            bucket_size: 100,
            topk,
            neighbor,
        }
    }

    /// 10 anchors on diagonal 0 (bucket 0), 3 on diagonal 500 (bucket 5),
    /// 1 on diagonal 1000 (bucket 10). With K=2, bucket 10 must be dropped.
    #[test]
    fn keeps_topk_buckets_only() {
        let mut anchors: Vec<Minimizer> = (0..10).map(|i| mk(0, 0, i * 100, i * 100)).collect();
        anchors.extend((0..3).map(|i| mk(0, 0, 500 + i * 100, i * 100)));
        anchors.push(mk(0, 0, 1000, 0));
        anchors.sort_by_key(|a| a.x); // filter contract: input sorted by x

        let kept = filter_anchors_diagonal(&mut anchors, &params(2, 0));
        assert_eq!(kept, 13);
        assert!(anchors.iter().all(|a| diagonal(a).div_euclid(100) <= 5));
        assert_eq!(anchors.len(), 13);
        // Order preserved (sorted by x).
        assert!(anchors.windows(2).all(|w| w[0].x <= w[1].x));
    }

    #[test]
    fn neighbor_expansion_keeps_adjacent_buckets() {
        let mut anchors: Vec<Minimizer> = (0..10).map(|i| mk(0, 0, i * 100, i * 100)).collect();
        anchors.push(mk(0, 0, 150, 0)); // diagonal 150 -> bucket 1 (neighbor of top bucket 0)

        filter_anchors_diagonal(&mut anchors, &params(1, 1));
        assert_eq!(anchors.len(), 11); // bucket 1 anchor kept via +-1 expansion
        assert!(anchors.iter().any(|a| diagonal(a) == 150));
    }

    #[test]
    fn ties_broken_by_smaller_bucket() {
        // Buckets 0 and 7 both have 3 anchors; K=1 keeps bucket 0 only.
        let mut anchors: Vec<Minimizer> = (0..3).map(|i| mk(0, 0, i * 100, i * 100)).collect();
        anchors.extend((0..3).map(|i| mk(0, 0, 700 + i * 100, i * 100)));

        filter_anchors_diagonal(&mut anchors, &params(1, 0));
        assert_eq!(anchors.len(), 3);
        assert!(anchors.iter().all(|a| diagonal(a).div_euclid(100) == 0));
    }

    #[test]
    fn groups_are_independent_per_ref_and_strand() {
        // ref 0: heavy diagonal at bucket 0, noise at bucket 20.
        // ref 1: its own diagonal at bucket 20 must survive.
        let mut anchors: Vec<Minimizer> = (0..5).map(|i| mk(0, 0, i * 100, i * 100)).collect();
        anchors.push(mk(0, 0, 2000, 0));
        anchors.extend((0..2).map(|i| mk(1, 0, 2000 + i * 100, i * 100)));
        // Negative-strand group of ref 0: independent histogram.
        anchors.extend((0..2).map(|i| mk(0, 1, 3000 + i * 100, i * 100)));
        anchors.push(mk(0, 1, 5000, 0));

        filter_anchors_diagonal(&mut anchors, &params(1, 0));
        assert_eq!(anchors.len(), 5 + 2 + 2);
        // ref 0 forward: noise dropped; ref 1 forward: its diagonal kept.
        assert!(anchors.iter().any(|a| a.ref_id() == 1 && diagonal(a).div_euclid(100) == 20));
        assert!(anchors.iter().all(|a| !(a.ref_id() == 0 && a.ref_id_strand() == (0u64 << 1) && diagonal(a).div_euclid(100) == 20)));
    }

    #[test]
    fn negative_diagonals_use_floor_semantics() {
        // d = -150 -> bucket -2 (Python floor: -150 // 100 == -2).
        // Top bucket must be -2, not -1 (truncating division would give -1).
        let mut anchors = vec![
            mk(0, 0, 0, 150),     // d = -150 -> bucket -2
            mk(0, 0, 10, 160),    // d = -150 -> bucket -2
            mk(0, 0, 20, 170),    // d = -150 -> bucket -2
            mk(0, 0, 300, 50),    // d = 250 -> bucket 2
        ];
        filter_anchors_diagonal(&mut anchors, &params(1, 0));
        assert_eq!(anchors.len(), 3);
        assert!(anchors.iter().all(|a| diagonal(a) == -150));
    }

    #[test]
    fn no_op_when_disabled_or_empty() {
        let mut anchors = vec![mk(0, 0, 0, 0), mk(0, 0, 1000, 0)];
        let mut p = params(1, 1);
        p.enabled = false;
        filter_anchors_diagonal(&mut anchors, &p);
        assert_eq!(anchors.len(), 2);

        let mut empty: Vec<Minimizer> = Vec::new();
        assert_eq!(filter_anchors_diagonal(&mut empty, &params(2, 1)), 0);

        let mut p_bad = params(2, 1);
        p_bad.bucket_size = 0;
        filter_anchors_diagonal(&mut anchors, &p_bad);
        assert_eq!(anchors.len(), 2);
    }

    #[test]
    fn default_params_match_prototype() {
        let p = DiagFilterParams::default();
        assert!(!p.enabled);
        assert_eq!(p.bucket_size, 100);
        assert_eq!(p.topk, 2);
        assert_eq!(p.neighbor, 1);
    }
}
