use indices::LibSaisSa;
use packed_seq::{PackedSeqVec, SeqVec};
use s_index::SIndex;
use sketchers::{IdentityParams, MinimizerParams};
use tracing::trace;
use utils::{read_chromosomes, Timer};

use super::*;

#[test]
fn test_identity_simple() {
    let seq = PackedSeqVec::from_ascii(b"ACGTACGTACGTACGT");
    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: true,
        par: false,
    };
    let uindex = UIndex::build(&seq, sketcher, ms_index);
    let query = PackedSeqVec::from_ascii(b"ACGT");
    let mut occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
    occ.sort();
    assert_eq!(occ, vec![0, 4, 8, 12]);
}
#[test]
fn test_identity_positive() {
    let seq = PackedSeqVec::random(1000000);
    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: true,
        par: false,
    };
    let uindex = UIndex::build(&seq, sketcher, ms_index);
    for _ in 0..100 {
        let len: usize = rand::random_range(..100);
        let pos = rand::random_range(..seq.len() - len);
        let query = seq.slice(pos..pos + len);
        let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
        assert!(occ.contains(&pos));
        for &pos in &occ {
            assert_eq!(seq.slice(pos..pos + len), query);
        }
    }
}
#[test]
fn test_identity_negative() {
    let seq = PackedSeqVec::random(1000000);
    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: true,
        par: false,
    };
    let uindex = UIndex::build(&seq, sketcher, ms_index);
    for _ in 0..100 {
        let len = 32;
        let query = PackedSeqVec::random(len);
        let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
        assert_eq!(occ.len(), 0);
    }
}
#[test]
fn test_minspace_positive() {
    let seq = PackedSeqVec::random(1000000);

    let ms_index = &LibSaisSa {
        store_ms_seq: true,
        par: false,
    };

    for remap in [false, true] {
        for l in [10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = &MinimizerParams {
                    l,
                    k,
                    remap,
                    cacheline_ef: false,
                    skip_zero: false,
                };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..100 {
                    let len: usize = rand::random_range(l..l + 100);
                    let pos = rand::random_range(..seq.len() - len);
                    let query = seq.slice(pos..pos + len);

                    let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                    assert!(
                        uindex_occ.contains(&pos),
                        "l {l} k {k} remap {remap} pos {pos} query {query:?}"
                    );
                }
            }
        }
    }
}
#[test]
fn test_minspace_negative() {
    let seq = PackedSeqVec::random(1000000);

    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: true,
        par: false,
    };
    let index = UIndex::build(&seq, sketcher, ms_index);

    for remap in [false, true] {
        for l in [10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = &MinimizerParams {
                    l,
                    k,
                    remap,
                    cacheline_ef: false,
                    skip_zero: false,
                };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..100 {
                    let len: usize = rand::random_range(l..l + 100);
                    let query = PackedSeqVec::random(len);

                    let mut index_occ = index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                    let mut uindex_occ =
                        uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                    index_occ.sort();
                    uindex_occ.sort();
                    assert_eq!(index_occ, uindex_occ, "l {l} k {k} remap {remap}");
                }
            }
        }
    }
}
#[test]
fn test_identity_positive_noms() {
    let seq = PackedSeqVec::random(1000000);
    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: false,
        par: false,
    };
    let uindex = UIndex::build(&seq, sketcher, ms_index);
    for _ in 0..100 {
        let len = rand::random_range(..100usize);
        let pos = rand::random_range(..seq.len() - len);
        let query = seq.slice(pos..pos + len);
        let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
        assert!(occ.contains(&pos));
        for &pos in &occ {
            assert_eq!(seq.slice(pos..pos + len), query);
        }
    }
}
#[test]
fn test_identity_negative_noms() {
    let seq = PackedSeqVec::random(1000000);
    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: false,
        par: false,
    };
    let uindex = UIndex::build(&seq, sketcher, ms_index);
    for _ in 0..100 {
        let len = 16;
        let query = PackedSeqVec::random(len);
        let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
        assert_eq!(occ.len(), 0);
    }
}
#[test]
fn test_minspace_positive_noms() {
    let seq = PackedSeqVec::random(1000000);

    let ms_index = &LibSaisSa {
        store_ms_seq: false,
        par: false,
    };

    for remap in [false] {
        for l in [10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = &MinimizerParams {
                    l,
                    k,
                    remap,
                    cacheline_ef: false,
                    skip_zero: false,
                };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..100 {
                    let len = l + rand::random_range(..100usize);
                    let pos = rand::random_range(..seq.len() - len);
                    let query = seq.slice(pos..pos + len);

                    let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                    assert!(
                        uindex_occ.contains(&pos),
                        "l {l} k {k} remap {remap} pos {pos} query {query:?}"
                    );
                }
            }
        }
    }
}
#[test]
fn test_minspace_negative_noms() {
    let seq = PackedSeqVec::random(1000000);

    let sketcher = &IdentityParams { skip_zero: false };
    let ms_index = &LibSaisSa {
        store_ms_seq: false,
        par: false,
    };
    let index = UIndex::build(&seq, sketcher, ms_index);

    for remap in [false] {
        for l in [10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = &MinimizerParams {
                    l,
                    k,
                    remap,
                    cacheline_ef: false,
                    skip_zero: false,
                };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..100 {
                    let len = l + rand::random_range(..100usize);
                    let query = PackedSeqVec::random(len);

                    let mut index_occ = index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                    let mut uindex_occ =
                        uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                    index_occ.sort();
                    uindex_occ.sort();
                    assert_eq!(index_occ, uindex_occ, "l {l} k {k} remap {remap}");
                }
            }
        }
    }
}

#[test]
fn test_s_index() {
    let seq = PackedSeqVec::random(1000000);

    for l in [1, 10, 100] {
        for k in [1, 2, 3, 4, 5, 6, 7, 8] {
            if k > l {
                continue;
            }
            let sindex = SIndex::build(&seq, k, l);
            for it in 0..100 {
                let len = l + rand::random_range(..100usize);
                let pos = rand::random_range(..seq.len() - len);
                let query = seq.slice(pos..pos + len);

                let uindex_occ = sindex.query(query).unwrap().collect::<Vec<_>>();
                assert!(
                    uindex_occ.contains(&pos),
                    "l {l} k {k} it {it} pos {pos} occurrences {}\n{uindex_occ:?}\n query {query:?}",
                    uindex_occ.len()
                );
            }
        }
    }
}

#[test]
#[ignore = "needs human-genome.fa"]
fn human_genome() {
    let (seq, ranges) = read_chromosomes::<PackedSeqVec>(1);

    let ql = 256;
    for (k, l) in [(8, 32), (16, 64)] {
        for remap in [false, true] {
            for store_seq in [false, true] {
                let mut timer = Timer::new("Build");
                trace!("remap {remap} l {l} k {k}");
                let sketcher = &MinimizerParams {
                    l,
                    k,
                    remap,
                    cacheline_ef: false,
                    skip_zero: false,
                };
                let ms_index = &LibSaisSa {
                    store_ms_seq: store_seq,
                    par: false,
                };
                let uindex = UIndex::build_with_ranges(&seq, &ranges, sketcher, ms_index);
                timer.next("Query");
                for _ in 0..100000 {
                    let pos = rand::random_range(..seq.len() - ql);
                    let query = seq.slice(pos..pos + ql);

                    let bad_ranges_before = uindex.query_stats.borrow().bad_ranges;
                    let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                    let bad_ranges_after = uindex.query_stats.borrow().bad_ranges;
                    assert!(uindex_occ.contains(&pos) || bad_ranges_after != bad_ranges_before);
                }
            }
        }
    }
}
