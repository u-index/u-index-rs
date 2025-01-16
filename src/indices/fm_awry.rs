use awry::alphabet::SymbolAlphabet;
use awry::fm_index::{FmBuildArgs, FmIndex};
use itertools::Itertools;
use mem_dbg::MemSize;
use packed_seq::PackedSeq;
use serde_json::Value;
use tracing::{info, trace, warn};

use crate::{Index, IndexBuilder};

#[derive(Debug, Clone, Copy)]
pub struct FmAwryParams {
    pub sa_sampling: usize,
}

#[derive(MemSize)]
pub struct FmAwry {
    fm: awry::fm_index::FmIndex,
    explode: bool,
}

impl IndexBuilder for FmAwryParams {
    fn try_build_with_stats(
        &self,
        text: Vec<u8>,
        width: usize,
        stats: &crate::utils::Stats,
    ) -> Option<Box<dyn Index>> {
        // AWRY does not support generic ASCII alphabet, so we 'explode' each byte into 4 DNA characters.
        let unpacked = PackedSeq {
            seq: &text,
            offset: 0,
            len: 4 * text.len(),
        }
        .unpack();
        stats.set_val("index", Value::String("FM-awry".to_string()));
        stats.set("index_width", width);
        stats.set("index_sa_sampling", self.sa_sampling as u64);

        let max = text.iter().copied().max().unwrap();
        trace!("Max value in text: {}", max);
        if max == 4 {
            warn!("AWRY takes packed input, or packed values in 0..4. `skip_zero: true` does not make sense here.");
            return None;
        }
        let explode = max >= 4;
        trace!("Explode: {explode}");

        let mut fasta = b">seq\n".to_vec();
        if !explode {
            fasta.extend(text.iter().copied().map(|x| packed_seq::unpack(x)));
        } else {
            // AWRY does not support generic ASCII alphabet, so we 'explode' each byte into 4 DNA characters.
            let unpacked = PackedSeq {
                seq: &text,
                offset: 0,
                len: 4 * text.len(),
            }
            .unpack();

            info!("Build AWRY on length {}", unpacked.len());

            fasta.extend(unpacked);
        }

        fasta.push(b'\n');
        let path = "/tmp/input.fa";
        // Write text to input file.
        std::fs::write(path, fasta).unwrap();

        let build_args = FmBuildArgs {
            input_file_src: path.into(),
            suffix_array_output_src: None,
            suffix_array_compression_ratio: Some(self.sa_sampling.try_into().unwrap()),
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: true,
        };

        let fm = FmIndex::new(&build_args).unwrap();
        std::fs::remove_file(path).unwrap();
        Some(Box::new(FmAwry { fm, explode }))
    }
}

impl Index for FmAwry {
    fn query(
        &self,
        pattern: &[u8],
        _seq: PackedSeq,
        _sketcher: &dyn crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize>> {
        // AWRY does not support generic ASCII alphabet, so we 'explode' each byte into 4 DNA characters.
        let unpacked = unsafe {
            String::from_utf8_unchecked(if !self.explode {
                pattern
                    .iter()
                    .copied()
                    .map(|x| packed_seq::unpack(x))
                    .collect_vec()
            } else {
                PackedSeq {
                    seq: &pattern,
                    offset: 0,
                    len: 4 * pattern.len(),
                }
                .unpack()
            })
        };

        let positions = self.fm.locate_string(&unpacked);
        let explode = self.explode;
        Box::new(
            positions
                .into_iter()
                .map(|x| x.local_position())
                .filter_map(move |x| {
                    if !explode {
                        Some(x)
                    } else {
                        if x % 4 == 0 {
                            Some(x / 4)
                        } else {
                            None
                        }
                    }
                }),
        )
    }
}
