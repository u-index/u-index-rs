use std::hint::black_box;

use super::*;
use packed_seq::SeqVec;

impl<SV: SeqVec> UIndex<SV> {
    /// Take `count` random substrings with length `len` and time querying them.
    pub fn gen_query_positions(&self, len: usize, count: usize) -> Vec<(usize, usize)> {
        let n = self.seq.len();

        (0..count)
            .map(|_| {
                let pos = rand::random::<usize>() % (n - len);
                (pos, pos + len)
            })
            .collect::<Vec<_>>()
    }

    /// Take `count` random substrings with length `len` and time querying them.
    pub fn bench_positive(&self, queries: &[(usize, usize)]) -> f64 {
        let start = std::time::Instant::now();

        let mut num_matches = 0;
        for &(s, e) in queries {
            num_matches += self.query(self.seq.slice(s..e)).unwrap().count();
        }
        black_box(num_matches);

        start.elapsed().as_secs_f64()
    }
}
