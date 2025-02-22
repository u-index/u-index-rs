use std::{
    cell::Cell,
    collections::HashMap,
    ops::Range,
    path::Path,
    sync::{LazyLock, Mutex},
};

use mem_dbg::{MemSize, SizeFlags};
use packed_seq::SeqVec;
use serde_json::{value::Value, Number};
use tracing::{info, trace};

thread_local! {
    static TIMER_DEPTH: Cell<usize> = Cell::new(0);
}

pub struct Timer<'s> {
    name: &'static str,
    stats: Option<&'s Stats>,
    start: std::time::Instant,
    depth: usize,
    info: bool,
}

impl<'s> Timer<'s> {
    pub fn new(name: &'static str) -> Self {
        let depth = TIMER_DEPTH.with(|d| d.get());
        TIMER_DEPTH.with(|d| d.set(depth + 1));
        Self {
            name,
            stats: None,
            start: std::time::Instant::now(),
            depth,
            info: false,
        }
    }
    pub fn new_stats(name: &'static str, stats: &'s Stats) -> Self {
        let depth = TIMER_DEPTH.with(|d| d.get());
        TIMER_DEPTH.with(|d| d.set(depth + 1));
        Self {
            name,
            stats: Some(stats),
            start: std::time::Instant::now(),
            depth,
            info: false,
        }
    }
    pub fn info(mut self) -> Self {
        self.info = true;
        self
    }
    pub fn next(&mut self, name: &'static str) {
        self.log();
        self.name = name;
        self.start = std::time::Instant::now();
    }

    fn log(&self) {
        let elapsed = self.start.elapsed();
        if let Some(stats) = self.stats {
            stats.add(self.name, elapsed.as_secs_f32());
        }
        let mut prefix = String::new();
        for _ in 0..self.depth {
            prefix.push_str(" ");
        }
        if self.info {
            info!("{prefix} {:<30}: {:.3?}", self.name, elapsed,);
        } else {
            trace!("{prefix} {:<30}: {:.3?}", self.name, elapsed,);
        }
    }
}

impl<'s> Drop for Timer<'s> {
    fn drop(&mut self) {
        self.log();
        TIMER_DEPTH.with(|d| d.set(d.get() - 1));
    }
}

fn init_trace() {
    use tracing::level_filters::LevelFilter;
    use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

    tracing_subscriber::registry()
        .with(tracing_subscriber::fmt::layer().with_writer(std::io::stderr))
        .with(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(LevelFilter::TRACE.into())
                .from_env_lossy(),
        )
        .init();
}

pub static INIT_TRACE: LazyLock<()> = LazyLock::new(init_trace);

#[derive(Default)]
pub struct Stats {
    stats: Mutex<HashMap<&'static str, Value>>,
}

impl MemSize for Stats {
    fn mem_size(&self, _flags: SizeFlags) -> usize {
        0
    }
}

impl Stats {
    pub fn set_val(&self, name: &'static str, value: Value) {
        self.stats.lock().unwrap().insert(name, value);
    }
    pub fn set(&self, name: &'static str, value: impl num_traits::ToPrimitive) {
        self.stats.lock().unwrap().insert(
            name,
            Value::Number(Number::from_f64(value.to_f64().unwrap()).unwrap()),
        );
    }
    pub fn add(&self, name: &'static str, value: f32) {
        let mut lock = self.stats.lock().unwrap();
        let entry = lock.entry(name).or_default();
        let old = entry
            .as_number()
            .map_or(0.0, |n| n.as_f64().unwrap_or_default());
        let new = old + value as f64;
        *entry = Value::Number(Number::from_f64(new).unwrap());
    }
    pub fn into(self) -> HashMap<&'static str, Value> {
        self.stats.into_inner().unwrap()
    }
    pub fn clone(&self) -> Self {
        Self {
            stats: Mutex::new(self.stats.lock().unwrap().clone()),
        }
    }
}

pub fn read_fastq<SV: SeqVec>(path: &Path) -> Vec<SV> {
    *INIT_TRACE;
    let _timer = Timer::new("Reading");
    let mut reads = vec![];
    let mut reader = needletail::parse_fastx_file(path).unwrap();
    while let Some(r) = reader.next() {
        let text = r
            .unwrap()
            .seq()
            .iter()
            .filter_map(|&b| if b == b'N' { None } else { Some(b) })
            .collect::<Vec<_>>();
        reads.push(SV::from_ascii(&text));
    }
    reads
}

pub fn read_chromosomes<SV: SeqVec>(cnt_max: usize) -> (SV, Vec<Range<usize>>) {
    *INIT_TRACE;
    let _timer = Timer::new("Reading");
    let Ok(mut reader) = needletail::parse_fastx_file("human-genome.fa") else {
        panic!("Did not find human-genome.fa. Add/symlink it to test runtime on it.");
    };
    let mut seq = SV::default();
    let mut ranges = vec![];
    let mut cnt = 0;
    while let Some(r) = reader.next() {
        let text = r
            .unwrap()
            .seq()
            .iter()
            .filter_map(|&b| if b == b'N' { None } else { Some(b) })
            .collect::<Vec<_>>();
        ranges.push(seq.push_ascii(&text));
        cnt += 1;
        if cnt == cnt_max {
            break;
        }
    }
    trace!(
        "Read human genome: {cnt} chromosomes of total length {}Mbp and size {}MB",
        seq.len() / 1000000,
        seq.mem_size(SizeFlags::default()) / 1000000
    );
    (seq, ranges)
}

pub fn read_human_genome<SV: SeqVec>() -> (SV, Vec<Range<usize>>) {
    read_chromosomes(usize::MAX)
}
