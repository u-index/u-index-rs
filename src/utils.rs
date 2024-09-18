use std::{cell::Cell, sync::LazyLock};

use tracing::trace;

thread_local! {
    static TIMER_DEPTH: Cell<usize> = Cell::new(0);
}

pub struct Timer {
    name: &'static str,
    start: std::time::Instant,
    depth: usize,
}

impl Timer {
    pub fn new(name: &'static str) -> Self {
        let depth = TIMER_DEPTH.with(|d| d.get());
        TIMER_DEPTH.with(|d| d.set(depth + 1));
        Self {
            name,
            start: std::time::Instant::now(),
            depth,
        }
    }
    pub fn next(&mut self, name: &'static str) {
        let elapsed = self.start.elapsed();
        let mut prefix = String::new();
        for _ in 0..self.depth {
            prefix.push_str(" ");
        }
        trace!("{prefix} {:<30}: {:.3?}", self.name, elapsed,);
        self.name = name;
        self.start = std::time::Instant::now();
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed();
        let mut prefix = String::new();
        for _ in 0..self.depth {
            prefix.push_str(" ");
        }
        trace!("{prefix} {:<30}: {:.3?}", self.name, elapsed,);
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
