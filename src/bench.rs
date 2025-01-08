use std::hint::black_box;

use super::*;
use packed_seq::{Seq, SeqVec};
use s_index::SIndex;
use tracing::trace;

/// Take `count` random substrings with length `len` and time querying them.
pub fn gen_query_positions<'i>(seq: impl Seq<'i>, len: usize, count: usize) -> Vec<(usize, usize)> {
    (0..count)
        .map(|_| {
            let pos = rand::random::<usize>() % (seq.len() - len);
            (pos, pos + len)
        })
        .collect::<Vec<_>>()
}

impl UIndex {
    /// Take `count` random substrings with length `len` and time querying them.
    pub fn bench_positive(&self, queries: &[(usize, usize)]) -> f64 {
        let start = std::time::Instant::now();

        let mut num_matches = 0;
        let mut i = 0usize;
        for &(s, e) in queries {
            num_matches += self.query(self.seq.slice(s..e)).unwrap().count();
            i += 1;
            if i.is_power_of_two() {
                trace!("Processed {i} queries with {num_matches} matches");
            }
        }
        black_box(num_matches);

        start.elapsed().as_secs_f64()
    }
}

// Yes yes, this is just a plain copy of the corresponding UIndex code. I know.
// TODO: Clean up into a single trait.
impl<SV: SeqVec> SIndex<SV> {
    /// Take `count` random substrings with length `len` and time querying them.
    pub fn bench_positive(&self, queries: &[(usize, usize)]) -> f64 {
        let start = std::time::Instant::now();

        let mut num_matches = 0;

        let mut i = 0usize;
        for &(s, e) in queries {
            num_matches += self.query(self.seq.slice(s..e)).unwrap().count();
            i += 1;
            if i.is_power_of_two() {
                trace!("Processed {i} queries with {num_matches} matches");
            }
        }
        black_box(num_matches);

        start.elapsed().as_secs_f64()
    }
}
