use std::hint::black_box;

use super::*;
use packed_seq::SeqVec;
use s_index::SIndex;
use tracing::trace;

pub fn gen_queries<'s, SV: SeqVec>(seq: &'s SV, len: usize, count: usize) -> Vec<SV> {
    (0..count)
        .map(|_| {
            let pos = rand::random_range(..seq.len() - len);
            let mut sv = SV::default();
            sv.push_seq(seq.slice(pos..pos + len));
            sv
        })
        .collect::<Vec<_>>()
}

impl<'s, SV: SeqVec + 'static> UIndex<'s, SV> {
    /// Take `count` random substrings with length `len` and time querying them.
    pub fn bench(&self, queries: &[SV]) -> f64 {
        let start = std::time::Instant::now();

        let mut num_matches = 0;
        let mut i = 0usize;
        for seq in queries {
            if let Some(locate_it) = self.query(seq.as_slice()) {
                num_matches += locate_it.count();
            }
            i += 1;
            if i.is_power_of_two() {
                if i.trailing_zeros() % 4 == 0 {
                    trace!("Processed {i} queries with {num_matches} matches");
                }
            }
        }
        black_box(num_matches);

        start.elapsed().as_secs_f64()
    }
}

// Yes yes, this is just a plain copy of the corresponding UIndex code. I know.
// TODO: Clean up into a single trait.
impl<'s, SV: SeqVec> SIndex<'s, SV> {
    /// Take `count` random substrings with length `len` and time querying them.
    pub fn bench(&self, queries: &[SV]) -> f64 {
        let start = std::time::Instant::now();

        let mut num_matches = 0;

        let mut i = 0usize;
        for seq in queries {
            if let Some(locate_it) = self.query(seq.as_slice()) {
                num_matches += locate_it.count();
            }
            i += 1;
            if i.is_power_of_two() {
                if i.trailing_zeros() % 4 == 0 {
                    trace!("Processed {i} queries with {num_matches} matches");
                }
            }
        }
        black_box(num_matches);

        start.elapsed().as_secs_f64()
    }
}
