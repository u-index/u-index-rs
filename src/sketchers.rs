mod minimizers;

use mem_dbg::{MemDbg, MemSize};
pub use minimizers::MinimizerParams;
use minimizers::MinimizerSketcher;
use packed_seq::{AsciiSeq, Seq};

use crate::{utils::Stats, MsSequence, SketchError, Sketcher, SketcherBuilder, S};

/// 'Sketch' the sequence to itself.
/// Convenient for testing purposes.
#[derive(Clone, Copy)]
pub struct IdentityParams;

#[derive(MemSize, MemDbg)]
pub struct Identity {
    len: usize,
}

impl SketcherBuilder for IdentityParams {
    type Sketcher = Identity;

    fn sketch_with_stats(&self, seq: AsciiSeq, _stats: &Stats) -> (Self::Sketcher, MsSequence) {
        let seq = seq.to_vec().seq;
        (Identity { len: seq.len() }, MsSequence(seq))
    }
}

impl Sketcher for Identity {
    fn width(&self) -> usize {
        1
    }

    fn k(&self) -> usize {
        1
    }

    fn len(&self) -> usize {
        self.len
    }

    fn sketch(&self, seq: S) -> Result<(MsSequence, usize), SketchError> {
        Ok((MsSequence(seq.to_vec().seq), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        Some(ms_pos)
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        Some(ms_seq[ms_pos] as usize)
    }

    fn get_ms_minimizer_via_plaintext(&self, seq: S, ms_pos: usize) -> Option<usize> {
        Some(seq.0[ms_pos] as usize)
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

    fn sketch(&self, seq: AsciiSeq<'_>) -> Result<(MsSequence, usize), SketchError> {
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

    fn get_ms_minimizer_via_plaintext(&self, seq: S, ms_pos: usize) -> Option<usize> {
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
