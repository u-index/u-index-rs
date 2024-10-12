use std::collections::HashMap;

use indices::IndexBuilderEnum;
use mem_dbg::SizeFlags;
use packed_seq::SeqVec;
use pyo3::pyclass;
use sketchers::SketcherBuilderEnum;
use tracing::trace;
use utils::{Timer, INIT_TRACE};

pub mod indices;
pub mod sketchers;
pub mod traits;
mod u_index;
pub mod utils;

pub use traits::*;
pub use u_index::*;

mod bench;
#[cfg(feature = "python-bindings")]
mod py;

// Terminology and variables:
// - Index: a datastructure that returns all locations where a pattern matches.
// - k-mer and minimizer are the used interchangeably.
// - k: the k-mer length of the minimizers.
// - l: the length of guaranteed preserved matches.
// - w = l - k + 1: the number of k-mers in a preserved window.
// - seq: the input sequence.
//   NOTE: Avoid 'text', always use 'seq'.
// - MS: MinSpace
// - ms_seq: the minimizer-space compressed text.
// - plain seq: the original sequence.

// TODO: Generalize the alphabet.
// - Currently the underlying index is assumed to only support `u8` alphabet.
// - Currently we only take `u8` encoded 2-bit alphabet as input. That should be extended to
//   2-bit encoded data and full 8-bit alphabets.
// TODO: Generalize underlying KmerVal and Pos types.
// TODO: Streaming input
// TODO: Rolling hash for O(1) checking.
// TODO: True constant time checking?
// TODO: Consider whether `Index` should own the input.
// TODO: minmers?
// TODO: Randomize the minimizer order
// TODO: Perfect hash the minimizers
// TODO: Reserve ms ID for sequence ends.
// TODO: Count false positives and other events. Add tracing or return them directly.
//
// TODO: Drop the minimizer space string, and instead compare against plain-space minimizers.

/// A minimizer space sequence.
pub struct MsSequence(Vec<u8>);

pub fn read_chromosomes<SV: SeqVec>(cnt_max: usize) -> SV {
    *INIT_TRACE;
    let _timer = Timer::new("Reading");
    let Ok(mut reader) = needletail::parse_fastx_file("human-genome.fa") else {
        panic!("Did not find human-genome.fa. Add/symlink it to test runtime on it.");
    };
    let mut seq = SV::default();
    let mut cnt = 0;
    while let Some(r) = reader.next() {
        seq.push_ascii(&r.unwrap().seq());
        cnt += 1;
        if cnt == cnt_max {
            break;
        }
    }
    trace!(
        "Read human genome: {cnt} chromosomes of total length {}Mbp and size {}MB",
        seq.len() / 1000000,
        seq.mem_size(SizeFlags::default()) / 1000000
    );
    seq
}

pub fn read_human_genome<SV: SeqVec>() -> SV {
    read_chromosomes(usize::MAX)
}

#[cfg(test)]
mod test {
    use indices::DivSufSortSa;
    use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
    use sketchers::{IdentityParams, MinimizerParams};
    use tracing::trace;

    use super::*;

    #[test]
    fn test_identity_simple() {
        let seq = AsciiSeqVec::from_ascii(b"ACGTACGTACGTACGT");
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq, sketcher, ms_index);
        let query = AsciiSeqVec::from_ascii(b"ACGT");
        let mut occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
        occ.sort();
        assert_eq!(occ, vec![0, 4, 8, 12]);
    }
    #[test]
    fn test_identity_positive() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
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
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = 32;
            let query = AsciiSeqVec::random(len);
            let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive() {
        let seq = PackedSeqVec::random(1000000);

        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });

        for remap in [false, true] {
            for l in [10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
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
        let seq = AsciiSeqVec::random(1000000);

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let index = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false, true] {
            for l in [10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let query = AsciiSeqVec::random(len);

                        let mut index_occ =
                            index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
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
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
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
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = 16;
            let query = AsciiSeqVec::random(len);
            let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive_noms() {
        let seq = PackedSeqVec::random(1000000);

        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });

        for remap in [false] {
            for l in [10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
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
        let seq = AsciiSeqVec::random(1000000);

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let index = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false] {
            for l in [10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let query = AsciiSeqVec::random(len);

                        let mut index_occ =
                            index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
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
    #[ignore = "needs human-genome.fa"]
    fn human_genome() {
        let seq: PackedSeqVec = read_chromosomes(1);

        let ql = 256;
        let compress = true;
        for (k, l) in [(8, 32), (16, 64)] {
            for remap in [false, true] {
                for store_seq in [false, true] {
                    let mut timer = Timer::new("Build");
                    trace!("remap {remap} l {l} k {k}");
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
                        store_ms_seq: store_seq,
                        compress,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    timer.next("Query");
                    for _ in 0..100000 {
                        let pos = rand::random::<usize>() % (seq.len() - ql);
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
}
