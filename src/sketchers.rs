mod minimizers;

use mem_dbg::{MemDbg, MemSize};
pub use minimizers::MinimizerParams;
use minimizers::MinimizerSketcher;

use crate::{utils::Stats, MsSequence, Seq, SketchError, Sketcher, SketcherBuilder};

/// 'Sketch' the sequence to itself.
/// Convenient for testing purposes.
#[derive(Clone, Copy)]
pub struct IdentityParams;

#[derive(MemSize, MemDbg)]
pub struct Identity;

impl SketcherBuilder for IdentityParams {
    type Sketcher = Identity;

    fn sketch_with_stats(&self, seq: Seq, _stats: &Stats) -> (Self::Sketcher, MsSequence) {
        (Identity, MsSequence(seq.to_vec()))
    }
}

impl Sketcher for Identity {
    fn sketch(&self, seq: Seq) -> Result<(MsSequence, usize), SketchError> {
        Ok((MsSequence(seq.to_vec()), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        Some(ms_pos)
    }
}

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

    fn sketch_with_stats(&self, seq: Seq, stats: &Stats) -> (Self::Sketcher, MsSequence) {
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
    fn sketch(&self, seq: Seq) -> Result<(MsSequence, usize), SketchError> {
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
}
