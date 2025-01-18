use packed_seq::SeqVec;
use serde_json::Value;
use tracing::trace;

use super::suffix_array::SuffixArray;
use crate::{
    utils::{Stats, Timer},
    Index, IndexBuilder,
};

/// Build a 32-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy, Debug)]
pub struct DivSufSortSa {
    pub store_ms_seq: bool,
    pub compress: bool,
}

impl<SV: SeqVec> IndexBuilder<SV> for DivSufSortSa {
    fn build_with_stats(&self, ms_seq: Vec<u8>, width: usize, stats: &Stats) -> Box<dyn Index<SV>> {
        stats.set_val("index", Value::String("divsufsort".to_string()));
        stats.set("index_width", width);
        stats.set("index_store_ms_seq", self.store_ms_seq as u64);
        stats.set("index_compress", self.compress as u64);

        trace!("MS sequence length {}", ms_seq.len());
        stats.set("sequence length", ms_seq.len());
        let mut timer = Timer::new_stats("Building suffix array", stats);

        // If we do not store the ms_seq, first invert the byte order of each minimizer to make sorting aligned with sorting packed u64's.

        let mut sa = libdivsufsort_rs::divsufsort(&ms_seq).expect("suffix array");

        if self.compress {
            timer.next("Compress SA");
            sa.retain(|x| *x % width as i32 == 0);
        }

        Box::new(SuffixArray {
            sa,
            ms_seq: self.store_ms_seq.then(|| ms_seq),
        })
    }
}
