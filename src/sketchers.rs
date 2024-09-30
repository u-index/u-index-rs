mod minimizers;

use mem_dbg::{MemDbg, MemSize};
pub use minimizers::MinimizerParams;
use minimizers::MinimizerSketcher;
use packed_seq::{AsciiSeq, PackedSeq, Seq};

use crate::{utils::Stats, MsSequence, SketchError, Sketcher, SketcherBuilder};

/// 'Sketch' the sequence to itself.
/// Convenient for testing purposes.
#[derive(Clone, Copy)]
pub struct IdentityParams;

#[derive(MemSize, MemDbg)]
pub struct Identity;

impl SketcherBuilder<AsciiSeq<'_>> for IdentityParams {
    type Sketcher = Identity;

    fn sketch_with_stats(&self, seq: AsciiSeq, _stats: &Stats) -> (Self::Sketcher, MsSequence) {
        (Identity, MsSequence(seq.to_vec().seq))
    }
}

impl Sketcher<AsciiSeq<'_>> for Identity {
    fn width(&self) -> usize {
        1
    }

    fn sketch(&self, seq: AsciiSeq) -> Result<(MsSequence, usize), SketchError> {
        Ok((MsSequence(seq.to_vec().seq), 0))
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

impl SketcherBuilder<AsciiSeq<'_>> for SketcherBuilderEnum {
    type Sketcher = SketcherEnum;

    fn sketch_with_stats(&self, seq: AsciiSeq, stats: &Stats) -> (Self::Sketcher, MsSequence) {
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

impl Sketcher<AsciiSeq<'_>> for SketcherEnum {
    fn width(&self) -> usize {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.width(),
            SketcherEnum::Minimizer(sketcher) => Sketcher::<AsciiSeq>::width(sketcher),
        }
    }

    fn sketch(&self, seq: AsciiSeq<'_>) -> Result<(MsSequence, usize), SketchError> {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.sketch(seq),
            SketcherEnum::Minimizer(sketcher) => sketcher.sketch(seq),
        }
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        match self {
            SketcherEnum::Identity(sketcher) => sketcher.ms_pos_to_plain_pos(ms_pos),
            SketcherEnum::Minimizer(sketcher) => {
                Sketcher::<AsciiSeq>::ms_pos_to_plain_pos(sketcher, ms_pos)
            }
        }
    }
}

impl SketcherBuilder<PackedSeq<'_>> for SketcherBuilderEnum {
    type Sketcher = SketcherEnum;

    fn sketch_with_stats(&self, seq: PackedSeq, stats: &Stats) -> (Self::Sketcher, MsSequence) {
        match self {
            SketcherBuilderEnum::Minimizer(params) => {
                let (sketcher, ms_seq) = params.sketch_with_stats(seq, stats);
                (SketcherEnum::Minimizer(sketcher), ms_seq)
            }
            _ => unimplemented!(),
        }
    }
}

impl Sketcher<PackedSeq<'_>> for SketcherEnum {
    fn width(&self) -> usize {
        match self {
            SketcherEnum::Minimizer(sketcher) => Sketcher::<PackedSeq>::width(sketcher),
            _ => unimplemented!(),
        }
    }

    fn sketch(&self, seq: PackedSeq<'_>) -> Result<(MsSequence, usize), SketchError> {
        match self {
            SketcherEnum::Minimizer(sketcher) => sketcher.sketch(seq),
            _ => unimplemented!(),
        }
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        match self {
            SketcherEnum::Minimizer(sketcher) => {
                Sketcher::<PackedSeq>::ms_pos_to_plain_pos(sketcher, ms_pos)
            }
            _ => unimplemented!(),
        }
    }
}
