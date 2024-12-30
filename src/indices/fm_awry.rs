use std::path::Path;

use awry::alphabet::SymbolAlphabet;
use awry::fm_index::{FmBuildArgs, FmIndex};
use mem_dbg::MemSize;

use crate::{Index, IndexBuilder};

#[derive(Debug, Clone, Copy)]
pub struct FmAwryParams {
    pub sa_sampling: usize,
}

pub struct FmAwry {
    fm: awry::fm_index::FmIndex,
}

impl MemSize for FmAwry {
    fn mem_size(&self, _flags: mem_dbg::SizeFlags) -> usize {
        let path = Path::new("/tmp/fm_index.awry");
        self.fm.save(path).unwrap();
        let size = path.metadata().unwrap().len() as usize;
        std::fs::remove_file(path).unwrap();
        size
    }
}

impl FmAwry {}

impl IndexBuilder for FmAwryParams {
    type Index = FmAwry;

    fn build_with_stats(
        &self,
        text: Vec<u8>,
        _width: usize,
        _stats: &crate::utils::Stats,
    ) -> Self::Index {
        let path = Path::new("/tmp/input.fa");
        // Write text to input file.
        std::fs::write(path, text).unwrap();

        let build_args = FmBuildArgs {
            input_file_src: "test.fasta".to_owned(),
            suffix_array_output_src: None,
            suffix_array_compression_ratio: Some(self.sa_sampling.try_into().unwrap()),
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: true,
        };

        let fm = FmIndex::new(&build_args).unwrap();
        std::fs::remove_file(path).unwrap();
        FmAwry { fm }
    }
}

impl Index for FmAwry {
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        _seq: impl packed_seq::Seq<'s>,
        _sketcher: &impl crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's> {
        // convert bytes to String
        let pattern = unsafe { String::from_utf8_unchecked(pattern.to_vec()) };
        let positions = self.fm.locate_string(&pattern);
        Box::new(positions.into_iter().map(|x| x.local_position()))
    }
}
