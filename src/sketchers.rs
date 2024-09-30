mod identity;
mod minimizers;

use identity::Identity;
use mem_dbg::{MemDbg, MemSize};
use minimizers::MinimizerSketcher;
use packed_seq::Seq;

use crate::{utils::Stats, MsSequence, SketchError, Sketcher, SketcherBuilder};

pub use identity::IdentityParams;
pub use minimizers::MinimizerParams;

#[derive(Clone, Copy)]
pub enum SketcherBuilderEnum {
    IdentityParams(IdentityParams),
    Minimizer(MinimizerParams),
}

#[derive(MemSize, MemDbg)]
pub enum SketcherEnum {
    Identity(Identity),
    Minimizer(MinimizerSketcher),
}

impl SketcherBuilder for SketcherBuilderEnum {
    type Sketcher = SketcherEnum;

    fn sketch_with_stats<'s, S: Seq<'s>>(
        &self,
        seq: S,
        stats: &Stats,
    ) -> (Self::Sketcher, MsSequence) {
        match self {
            SketcherBuilderEnum::IdentityParams(identity) => {
                let (sketcher, ms_seq) = identity.sketch_with_stats(seq, stats);
                (SketcherEnum::Identity(sketcher), ms_seq)
            }
            SketcherBuilderEnum::Minimizer(params) => {
                let (sketcher, ms_seq) = params.sketch_with_stats(seq, stats);
                (SketcherEnum::Minimizer(sketcher), ms_seq)
            }
        }
    }
}

impl Sketcher for SketcherEnum {
    fn width(&self) -> usize {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.width(),
            SketcherEnum::Minimizer(sketcher) => sketcher.width(),
        }
    }

    fn k(&self) -> usize {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.k(),
            SketcherEnum::Minimizer(sketcher) => sketcher.k(),
        }
    }

    fn len(&self) -> usize {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.len(),
            SketcherEnum::Minimizer(sketcher) => sketcher.len(),
        }
    }

    fn sketch<'s>(&self, seq: impl Seq<'s>) -> Result<(MsSequence, usize), SketchError> {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.sketch(seq),
            SketcherEnum::Minimizer(sketcher) => sketcher.sketch(seq),
        }
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.ms_pos_to_plain_pos(ms_pos),
            SketcherEnum::Minimizer(sketcher) => sketcher.ms_pos_to_plain_pos(ms_pos),
        }
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.get_ms_minimizer(ms_seq, ms_pos),
            SketcherEnum::Minimizer(sketcher) => sketcher.get_ms_minimizer(ms_seq, ms_pos),
        }
    }

    fn get_ms_minimizer_via_plaintext<'s>(
        &self,
        seq: impl Seq<'s>,
        ms_pos: usize,
    ) -> Option<usize> {
        match self {
            SketcherEnum::Identity(sketcher) => {
                sketcher.get_ms_minimizer_via_plaintext(seq, ms_pos)
            }
            SketcherEnum::Minimizer(sketcher) => {
                sketcher.get_ms_minimizer_via_plaintext(seq, ms_pos)
            }
        }
    }
}
