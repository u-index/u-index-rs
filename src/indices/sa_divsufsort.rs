use tracing::trace;

use super::SuffixArray;
use crate::{
    utils::{Stats, Timer},
    IndexBuilder,
};

/// Build a 32-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSa {
    pub store_ms_seq: bool,
    pub compress: bool,
}

impl IndexBuilder for DivSufSortSa {
    type Index = SuffixArray;

    fn build_with_stats(&self, ms_seq: Vec<u8>, width: usize, stats: &Stats) -> Self::Index {
        let mut timer = Timer::new_stats("Building suffix array", stats);
        trace!("MS sequence length {}", ms_seq.len());
        stats.set("sequence length", ms_seq.len());

        // If we do not store the ms_seq, first invert the byte order of each minimizer to make sorting aligned with sorting packed u64's.

        let mut sa = libdivsufsort_rs::divsufsort(&ms_seq).expect("suffix array");

        if self.compress {
            timer.next("Compress SA");
            sa.retain(|x| *x % width as i32 == 0);
        }

        SuffixArray {
            sa,
            ms_seq: self.store_ms_seq.then(|| ms_seq),
        }
    }
}
