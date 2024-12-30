use packed_seq::{PackedSeqVec, SeqVec};
use uindex::{
    bench::gen_query_positions,
    indices::{DivSufSortSa, FmAwryParams, FmBioParams, IndexBuilderEnum, LibSaisSa},
    sketchers::{MinimizerParams, SketcherBuilderEnum},
    utils::{read_chromosomes, Timer},
    UIndex,
};

fn main() {
    let seq: PackedSeqVec = read_chromosomes(1);

    let sketch_params = SketcherBuilderEnum::Minimizer(MinimizerParams {
        k: 8,
        l: 16,
        remap: false,
        cacheline_ef: false,
    });

    let mut queries = vec![];

    // SAIS
    if true {
        let index_params = IndexBuilderEnum::LibSaisSa(LibSaisSa {
            store_ms_seq: true,
            par: false,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        queries = gen_query_positions(seq.as_slice(), 256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    // Fm Bio
    if true {
        let index_params = IndexBuilderEnum::FmBio(FmBioParams {
            occ_sampling: 64,
            sa_sampling: 64,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        // queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    // Fm Awry
    if true {
        let index_params = IndexBuilderEnum::FmAwry(FmAwryParams { sa_sampling: 64 });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        // queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    // SAIS parallel
    if false {
        let index_params = IndexBuilderEnum::LibSaisSa(LibSaisSa {
            store_ms_seq: true,
            par: true,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        // queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    // DivSufSort
    if false {
        let index_params = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        // let queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }
}
