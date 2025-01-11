mod identity;
mod minimizers;

use mem_dbg::MemSize;
use packed_seq::Seq;

use crate::{utils::Stats, MsSequence, SketchError, Sketcher, SketcherBuilder};

pub use identity::IdentityParams;
pub use minimizers::MinimizerParams;
