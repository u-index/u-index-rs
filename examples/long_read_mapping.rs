#![allow(unused)]
use std::{any::type_name, cell::Cell, collections::HashMap, path::PathBuf, sync::LazyLock};

use clap::Parser;
use packed_seq::{unpack_base, AsciiSeqVec, PackedSeqVec, Seq, SeqVec};
use serde_json::{from_slice, Number, Value};
use tracing::{info, warn};
use uindex::{
    bench::gen_queries,
    indices::LibSaisSa,
    s_index::SIndex,
    sketchers::{IdentityParams, MinimizerParams},
    utils::{read_chromosomes, read_fastq, Timer, INIT_TRACE},
    IndexBuilder, SketcherBuilder, UIndex,
};

#[derive(Parser)]
struct Args {
    /// Reads. One per line. Split
    #[clap(long)]
    reads: PathBuf,
    /// Split each read into queries of this length.
    #[clap(long)]
    ql: usize,
    /// Number of reads to use.
    #[clap(long)]
    num_reads: Option<usize>,

    /// Minimizer-k to use.
    #[clap(short)]
    k: usize,
    /// Minimizer-ell to use.
    #[clap(short)]
    l: usize,

    /// Output json file with statistics.
    #[clap(long)]
    output: Option<PathBuf>,

    /// Max number of matches to consider for each pattern.
    #[clap(long)]
    max_matches: Option<usize>,

    /// Limit number of sketch-space matches to consider.
    #[clap(long)]
    lim: Option<usize>,
}

static ARGS: LazyLock<Args> = LazyLock::new(|| Args::parse());

fn main() {
    color_backtrace::install();
    *INIT_TRACE;

    // Read input.
    let (seq, ranges) = read_chromosomes::<PackedSeqVec>(usize::MAX);

    let rss0 = max_rss();

    // Build index.
    let index = {
        let min_remap = MinimizerParams {
            k: ARGS.k,
            l: ARGS.l,
            remap: true,
            cacheline_ef: true,
            skip_zero: false,
        };
        let sais_ms = LibSaisSa {
            store_ms_seq: true,
            par: false,
        };

        UIndex::<PackedSeqVec>::try_build_with_ranges(&seq, &ranges, &min_remap, &sais_ms).unwrap()
    };

    let rss1 = max_rss();

    let reads = if ARGS.reads.extension() == Some("fastq".as_ref()) {
        read_fastq(&ARGS.reads)
    } else {
        let reads = std::fs::read(&ARGS.reads).unwrap();
        reads
            .trim_ascii_end()
            .split(|&c| c == b'\n')
            .map(PackedSeqVec::from_ascii)
            .collect::<Vec<_>>()
    };

    let mut matches = HashMap::<usize, usize>::new();
    let mut mismatches = HashMap::<usize, usize>::new();

    let mut old_matches = 0;
    let mut reads_with_matches = 0;
    let mut patterns_with_matches = 0;
    let start = std::time::Instant::now();
    let mut bigq = 0;
    for (ri, read) in reads
        .iter()
        .take(ARGS.num_reads.unwrap_or(usize::MAX))
        .enumerate()
    {
        // Split read into queries of length ql.
        let queries = (0..read.len() - ARGS.ql + 1)
            .step_by(ARGS.ql)
            .map(|i| read.slice(i..i + ARGS.ql).to_vec())
            .collect::<Vec<_>>();

        let start = std::time::Instant::now();

        let mut big = 0;
        let ql = queries.len();
        for (qi, seq) in queries.iter().enumerate() {
            if let Some(locate_it) = index.query_lim(seq.as_slice(), ARGS.lim) {
                let mismatches_1 = index.query_stats.borrow().mismatches;
                let cnt = locate_it
                    .take(ARGS.max_matches.unwrap_or(usize::MAX))
                    .count();
                let mismatches_2 = index.query_stats.borrow().mismatches;
                let mms = mismatches_2 - mismatches_1;
                *matches.entry(cnt).or_default() += 1;
                *mismatches.entry((mms + 1).ilog2() as usize).or_default() += 1;
                if cnt > 0 {
                    patterns_with_matches += 1;
                }
                if cnt > 200 {
                    big += 1;
                    // if big > 1 {
                    //     continue;
                    // }
                    let seq = seq
                        .as_slice()
                        .iter_bp()
                        .map(unpack_base)
                        .collect::<Vec<_>>();
                    let seq = String::from_utf8(seq).unwrap();
                    eprintln!("{ri:>4} {qi:>4}/{ql:>4} {cnt:>5}: {seq}");
                }
            }
        }
        if big > 0 {
            bigq += 1;
            eprintln!("{ri:>4}: {big:>} hot matches");
        }

        let new_matches = index.query_stats.borrow().matches;
        if new_matches > old_matches {
            reads_with_matches += 1;
        }
        old_matches = new_matches;
    }
    if bigq > 0 {
        eprintln!("{bigq:>}/{} hot reads", reads.len());
    }
    let query_time = start.elapsed().as_secs_f64();
    warn!(
        "Reads with matches: {} / {}",
        reads_with_matches,
        reads.len()
    );
    warn!(
        "Patterns with matches: {} / {}",
        patterns_with_matches,
        index.query_stats.borrow().queries
    );

    let mut matches = matches.into_iter().collect::<Vec<_>>();
    matches.sort_unstable_by_key(|(k, _)| *k);
    for (k, v) in matches {
        eprintln!("{}: {}", k, v);
    }
    let mut mismatches = mismatches.into_iter().collect::<Vec<_>>();
    mismatches.sort_unstable_by_key(|(k, _)| *k);
    for (k, v) in mismatches {
        eprintln!("{}: {}", k, v);
    }

    let mut stats = index.stats();
    stats.insert("rss0", Value::Number(Number::from(rss0)));
    stats.insert("rss1", Value::Number(Number::from(rss1)));
    stats.insert("rss2", Value::Number(Number::from(max_rss())));
    stats.insert(
        "query_time",
        Value::Number(Number::from_f64(query_time).unwrap()),
    );

    let stats_string = serde_json::to_string_pretty(&stats).unwrap();
    // info!("Stats:\n{stats_string}");
    let path = PathBuf::from("stats.json");
    let output = ARGS.output.as_ref().unwrap_or(&path);
    std::fs::write(output, stats_string).unwrap();
}

fn max_rss() -> usize {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    // On linux, the returned value is in kB.
    rusage.ru_maxrss as usize * 1024
}
