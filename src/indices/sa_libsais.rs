use std::collections::HashMap;

use packed_seq::SeqVec;
use serde_json::Value;
use tracing::trace;

use super::suffix_array::SuffixArray;
use crate::{
    utils::{Stats, Timer},
    Index, IndexBuilder,
};

/// Build a 32-bit suffix array using `libsais`.
#[derive(Clone, Copy, Debug)]
pub struct LibSaisSa {
    pub store_ms_seq: bool,
    pub par: bool,
}

impl<SV: SeqVec> IndexBuilder<SV> for LibSaisSa {
    fn build_with_stats(
        &self,
        mut ms_seq: Vec<u8>,
        width: usize,
        stats: &Stats,
    ) -> Box<dyn Index<SV>> {
        stats.set_val("index", Value::String("SA".to_string()));
        stats.set("index_width", width);
        stats.set("index_store_ms_seq", self.store_ms_seq as u64);
        stats.set("index_par", self.par as u64);

        trace!("MS sequence length {}", ms_seq.len());
        stats.set("sequence length", ms_seq.len());

        let _timer = Timer::new_stats("Building suffix array", stats);

        // If we do not store the ms_seq, first invert the byte order of each minimizer to make sorting aligned with sorting packed u64's.

        let mut sa;
        match width {
            0 => panic!("Width 0 is not allowed"),
            1 => {
                sa = vec![0; ms_seq.len()];
                if self.par {
                    #[cfg(feature = "openmp")]
                    libsais_rs::par::suffix_array_u8(&ms_seq, &mut sa).expect("suffix array");
                    #[cfg(not(feature = "openmp"))]
                    panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
                } else {
                    libsais_rs::suffix_array_u8(&ms_seq, &mut sa).expect("suffix array");
                }
            }
            2 => {
                trace!("Building 16-bit suffix array");
                trace!("Transmuting..");
                let (head, ms_seq, tail) = unsafe { ms_seq.as_mut_slice().align_to_mut::<u16>() };
                assert!(head.is_empty());
                assert!(tail.is_empty());
                // TODO: Instead prefer to write in the right way directly during sketching.
                trace!("Reversing byte order..");
                for x in ms_seq.iter_mut() {
                    *x = x.swap_bytes();
                }

                sa = vec![0; ms_seq.len()];
                trace!("Building suffix array");
                if self.par {
                    #[cfg(feature = "openmp")]
                    libsais_rs::par::suffix_array_u16(&ms_seq, &mut sa).expect("suffix array");
                    #[cfg(not(feature = "openmp"))]
                    panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
                } else {
                    libsais_rs::suffix_array_u16(&ms_seq, &mut sa).expect("suffix array");
                }
                trace!("Reversing byte order back..");
                for x in ms_seq.iter_mut() {
                    *x = x.swap_bytes();
                }
                trace!("Suffix array built");
                for x in sa.iter_mut() {
                    *x *= 2;
                }
            }
            3.. => {
                trace!("Building suffix array on {width}-byte input. First remapping to small i32 values.");
                let minimizers_vals = ms_seq.chunks(width).map(|x| {
                    let mut val = [0u8; 8];
                    val[8 - width..].copy_from_slice(x);
                    usize::from_be_bytes(val)
                });

                let mut vals_map = HashMap::new();
                for x in minimizers_vals.clone() {
                    vals_map.insert(x, 0i32);
                }
                let mut vals = vals_map.iter().map(|x| *x.0).collect::<Vec<_>>();
                vals.sort_unstable();
                trace!("MIN VAL {:?}", vals.iter().min().unwrap());
                trace!("MAX VAL {:?}", vals.iter().max().unwrap());
                for (i, x) in vals.iter().enumerate() {
                    *vals_map.get_mut(x).unwrap() = i as i32;
                }
                let alphabet_size = vals.len() + 1;
                drop(vals);
                trace!("alphabet size {alphabet_size}");
                assert!(alphabet_size < i32::MAX as usize);
                let mut remapped_minimizer_vals: Vec<i32> = minimizers_vals
                    .map(|x| *vals_map.get(&x).unwrap())
                    .collect();
                drop(vals_map);
                trace!("Building suffix array");
                sa = vec![0; ms_seq.len() / width];
                if self.par {
                    #[cfg(feature = "openmp")]
                    libsais_rs::par::suffix_array_i32(
                        &mut remapped_minimizer_vals,
                        &mut sa,
                        alphabet_size,
                    )
                    .expect("suffix array");
                    #[cfg(not(feature = "openmp"))]
                    panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
                } else {
                    libsais_rs::suffix_array_i32(
                        &mut remapped_minimizer_vals,
                        &mut sa,
                        alphabet_size as i32,
                    )
                    .expect("suffix array");
                }

                trace!("Suffix array built");
                for x in sa.iter_mut() {
                    *x *= width as i32;
                }
            } // 4 => {
              //     trace!("Building 32-bit suffix array");
              //     trace!("Transmuting..");
              //     let (head, ms_seq, tail) = unsafe { ms_seq.as_mut_slice().align_to_mut::<i32>() };
              //     assert!(head.is_empty(), "Head has size {}", head.len());
              //     assert!(tail.is_empty(), "Tail has size {}", tail.len());
              //     sa = vec![0; ms_seq.len()];
              //     trace!("Reversing byte order..");
              //     // FIXME: What exactly does alphabet size mean?
              //     let alphabet_size = i32::MAX;
              //     for x in ms_seq.iter_mut() {
              //         *x = x.swap_bytes();
              //     }
              //     trace!("Building suffix array");
              //     if self.par {
              //         #[cfg(feature = "openmp")]
              //         libsais_rs::par::suffix_array_i32(ms_seq, &mut sa, alphabet_size)
              //             .expect("suffix array");
              //         #[cfg(not(feature = "openmp"))]
              //         panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
              //     } else {
              //         libsais_rs::suffix_array_i32(ms_seq, &mut sa, alphabet_size)
              //             .expect("suffix array");
              //     }
              //     trace!("Reversing byte order back..");
              //     for x in ms_seq.iter_mut() {
              //         *x = x.swap_bytes();
              //     }
              //     trace!("Suffix array built");
              //     for x in sa.iter_mut() {
              //         *x *= 4;
              //     }
              // }
              // 8 => {
              //     trace!("Building 64-bit suffix array");
              //     trace!("Transmuting..");
              //     let (head, ms_seq, tail) = unsafe { ms_seq.as_mut_slice().align_to_mut::<i64>() };
              //     assert!(head.is_empty(), "Head has size {}", head.len());
              //     assert!(tail.is_empty(), "Tail has size {}", tail.len());
              //     let mut sa_64 = vec![0; ms_seq.len()];
              //     trace!("Reversing byte order..");
              //     for x in ms_seq.iter_mut() {
              //         *x = x.swap_bytes();
              //     }
              //     // FIXME: What exactly does alphabet size mean?
              //     let alphabet_size = i64::MAX;
              //     trace!("Building suffix array");
              //     if self.par {
              //         #[cfg(feature = "openmp")]
              //         libsais_rs::par::long_suffix_array_i64(ms_seq, &mut sa_64, alphabet_size)
              //             .expect("suffix array");
              //         #[cfg(not(feature = "openmp"))]
              //         panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
              //     } else {
              //         libsais_rs::long_suffix_array_i64(ms_seq, &mut sa_64, alphabet_size)
              //             .expect("suffix array");
              //     }
              //     trace!("Reversing byte order back..");
              //     for x in ms_seq.iter_mut() {
              //         *x = x.swap_bytes();
              //     }
              //     trace!("Copy to 32-bit SA");
              //     sa = sa_64.iter().map(|x| *x as i32).collect();
              //     trace!("Suffix array built");
              //     for x in sa_64.iter_mut() {
              //         *x *= 8;
              //     }
              // }
              // _ => {
              //     // TODO: For remaining alphabet sizes, it's probably better to re-code them to the next power of 2 size.
              //     sa = vec![0; ms_seq.len()];
              //     if self.par {
              //         #[cfg(feature = "openmp")]
              //         libsais_rs::par::suffix_array_u8(&ms_seq, &mut sa).expect("suffix array");
              //         #[cfg(not(feature = "openmp"))]
              //         panic!("Parallel suffix array construction with libsais required the 'openmp' feature.");
              //     } else {
              //         libsais_rs::suffix_array_u8(&ms_seq, &mut sa).expect("suffix array");
              //     }
              //     timer.next("Compress SA");
              //     sa.retain(|x| *x % width as i32 == 0);
              // }
        }

        Box::new(SuffixArray {
            sa,
            ms_seq: self.store_ms_seq.then(|| ms_seq),
        })
    }
}
