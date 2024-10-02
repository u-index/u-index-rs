use packed_seq::PackedSeqVec;
use uindex::{
    indices::{DivSufSortSa, IndexBuilderEnum, LibSaisSa},
    read_chromosomes, read_human_genome,
    sketchers::{MinimizerParams, SketcherBuilderEnum},
    utils::Timer,
    UIndex,
};

fn main() {
    let seq: PackedSeqVec = read_chromosomes(1);

    let sketch_params = SketcherBuilderEnum::Minimizer(MinimizerParams {
        k: 29,
        l: 128,
        remap: false,
        cacheline_ef: false,
    });

    let queries;

    {
        let index_params = IndexBuilderEnum::LibSaisSa(LibSaisSa {
            store_ms_seq: true,
            par: true,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    {
        let index_params = IndexBuilderEnum::LibSaisSa(LibSaisSa {
            store_ms_seq: true,
            par: false,
        });

        let u = UIndex::build(seq.clone(), sketch_params, index_params);

        // queries = u.gen_query_positions(256, 1000000);

        let _t = Timer::new("bench_positive");
        u.bench_positive(&queries);
    }

    {
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
