use std::sync::Arc;

use bio::alphabets::Alphabet;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use bio::data_structures::suffix_array::{suffix_array, SampledSuffixArray, SuffixArray};
use itertools::Itertools;
use mem_dbg::MemSize;
use packed_seq::PackedSeq;
use tracing::{trace, warn};

use crate::utils::{Stats, Timer};
use crate::{Index, IndexBuilder};

#[derive(Debug, Clone, Copy)]
pub struct FmBioParams {
    pub occ_sampling: u32,
    pub sa_sampling: usize,
}

pub struct FmBio {
    fm: FMIndex<Arc<BWT>, Arc<Less>, Arc<Occ>>,
    sampled_sa: SampledSuffixArray<Arc<BWT>, Arc<Less>, Arc<Occ>>,
}

impl MemSize for FmBio {
    fn mem_size(&self, flags: mem_dbg::SizeFlags) -> usize {
        // let fm_size = self.fm.mem_size(flags);
        let bwt_size = self.sampled_sa.bwt().mem_size(flags);
        let less_size = self.sampled_sa.less().mem_size(flags);
        warn!("OCC size not implemented");
        warn!("Sampled SA size not implemented");
        let occ_size = 0;
        let sampled_sa_size = 0;
        // let occ_size = self.sampled_sa.occ().mem_size(flags);
        // let sampled_sa_size = self.sampled_sa.mem_size(flags);
        bwt_size + less_size + occ_size + sampled_sa_size
    }
}

impl FmBio {
    #[allow(unused)]
    pub fn log_sizes(&self, stats: &Stats) {
        let bwt_size = self.sampled_sa.bwt().mem_size(Default::default()) as f32 / 1000000.;
        stats.add("bwt_size_MB", bwt_size);
        trace!("bwt size:   {bwt_size:>8.3} MB",);

        let less_size = self.sampled_sa.less().mem_size(Default::default()) as f32 / 1000000.;
        stats.add("less_size_MB", less_size);
        trace!("less size:  {less_size:>8.3} MB",);

        warn!("OCC size not implemented");
        warn!("Sampled SA size not implemented");
        // let occ_size = self.sampled_sa.occ().mem_size(Default::default()) as f32 / 1000000.;
        // stats.add("occ_size_MB", occ_size);
        // trace!("occ size:   {occ_size:>8.3} MB",);

        // let ssa_size = self.sampled_sa.mem_size(Default::default()) as f32 / 1000000.;
        // stats.add("ssa_size_MB", ssa_size);
        // trace!("ssa size:   {ssa_size:>8.3} MB",);
    }
}

impl IndexBuilder for FmBioParams {
    fn build_with_stats(
        &self,
        mut text: Vec<u8>,
        _width: usize,
        stats: &crate::utils::Stats,
    ) -> Box<dyn Index> {
        // Rust-bio expects a sentinel character at the end of the text.
        text.push(0);

        let alphabet = (0u8..=255).collect_vec();
        let alphabet = Alphabet::new(alphabet);
        let mut timer = Timer::new_stats("FM: suffix array", stats);
        let sa = suffix_array(&text);
        timer.next("FM: BWT");
        let bwt = Arc::new(bwt(&text, &sa));
        timer.next("FM: less");
        let less = Arc::new(less(&bwt, &alphabet));
        timer.next("FM: occ");
        let occ = Arc::new(Occ::new(&bwt, self.occ_sampling, &alphabet));
        timer.next("FM: sampled SA");
        let sampled_sa = sa.sample(
            &text,
            bwt.clone(),
            less.clone(),
            occ.clone(),
            self.sa_sampling,
        );
        timer.next("FM: index");
        let fm = FMIndex::new(bwt, less, occ);
        Box::new(FmBio { fm, sampled_sa })
    }
}

impl Index for FmBio {
    fn query(
        &self,
        pattern: &[u8],
        _seq: PackedSeq,
        _sketcher: &dyn crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize>> {
        let bsr = self.fm.backward_search(pattern.iter());
        let positions = match bsr {
            BackwardSearchResult::Complete(sai) => sai.occ(&self.sampled_sa),
            _ => Vec::<usize>::new(),
            // BackwardSearchResult::Partial(sai, _l) => sai.occ(&self.sampled_sa),
            // BackwardSearchResult::Absent => Vec::<usize>::new(),
        };
        Box::new(positions.into_iter())
    }
}
