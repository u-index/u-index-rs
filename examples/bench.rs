use packed_seq::{PackedSeqVec, SeqVec};
use sdsl_lite_fm::*;
use uindex::{
    bench::gen_query_positions,
    indices::{FmAwryParams, FmBioParams, FmSdslParams, LibSaisSa},
    sketchers::MinimizerParams,
    utils::{read_chromosomes, Timer},
    IndexBuilder, UIndex,
};

fn main() {
    color_backtrace::install();

    let seq: PackedSeqVec = read_chromosomes(1);

    let sketch_params = Box::new(MinimizerParams {
        k: 28,
        l: 128,
        remap: true,
        cacheline_ef: true,
        skip_zero: false,
    });

    let queries = gen_query_positions(seq.as_slice(), 256, 10000);

    let index_params: Vec<Box<dyn IndexBuilder>> = vec![
        Box::new(FmSdslParams::<FmIndexByte32Ptr, _>::new()),
        Box::new(FmSdslParams::<FmIndexInt32Ptr, _>::new()),
        Box::new(FmSdslParams::<FmIndexByte64Ptr, _>::new()),
        Box::new(FmSdslParams::<FmIndexInt64Ptr, _>::new()),
        Box::new(FmAwryParams { sa_sampling: 64 }),
        Box::new(LibSaisSa {
            store_ms_seq: true,
            par: false,
        }),
        Box::new(FmBioParams {
            occ_sampling: 64,
            sa_sampling: 64,
        }),
        // Box::new(LibSaisSa {
        //     store_ms_seq: true,
        //     par: true,
        // }),
        // Box::new(DivSufSortSa {
        //     store_ms_seq: true,
        //     compress: true,
        // }),
    ];
    for p in index_params {
        let u = UIndex::build(seq.clone(), sketch_params.clone(), p);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }
}
