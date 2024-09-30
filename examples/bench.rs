use uindex::{
    indices::{DivSufSortSa, IndexBuilderEnum},
    read_human_genome,
    sketchers::{MinimizerParams, SketcherBuilderEnum},
    utils::Timer,
    UIndex,
};

fn main() {
    let seq = read_human_genome();

    let sketch_params = SketcherBuilderEnum::Minimizer(MinimizerParams {
        k: 16,
        l: 64,
        remap: false,
        cacheline_ef: false,
    });
    let index_params = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
        store_ms_seq: true,
        compress: true,
    });

    let u = UIndex::build(seq, sketch_params, index_params);

    let queries = u.gen_query_positions(256, 1000000);

    let t = Timer::new("bench_positive");
    u.bench_positive(&queries);
    drop(t);
}
