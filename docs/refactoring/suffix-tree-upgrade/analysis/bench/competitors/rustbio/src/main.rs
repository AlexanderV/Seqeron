// rust-bio (crate `bio`) suffix array + FM-index benchmark harness.
// Same output protocol as every other harness:  R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes
//   usage: rustbio-bench <text-file> <query-prefix>
// Index = suffix_array() (SA-IS) + bwt + less + Occ(k=3) = FMIndex; locate uses the FULL suffix array (8 B/char, usize).
use bio::alphabets::Alphabet;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable, Interval};
use bio::data_structures::suffix_array::suffix_array;
use std::fs;
use std::time::Instant;

fn emit(op: &str, m: &str, v: f64, u: &str, n: &str) { println!("R\t{}\t{}\t{}\t{}\t{}", op, m, v, u, n); }
fn emiti(op: &str, m: &str, v: i64, u: &str, n: &str) { println!("R\t{}\t{}\t{}\t{}\t{}", op, m, v, u, n); }

fn status_kb(key: &str) -> i64 {
    let s = fs::read_to_string("/proc/self/status").unwrap_or_default();
    for l in s.lines() {
        if l.starts_with(key) {
            return l[key.len()..].split_whitespace().next().unwrap().parse::<i64>().unwrap() * 1024;
        }
    }
    0
}

fn iv(r: BackwardSearchResult) -> Option<Interval> {
    match r { BackwardSearchResult::Complete(i) => Some(i), _ => None }
}

fn read_lines(p: &str) -> Vec<Vec<u8>> {
    match fs::read(p) {
        Ok(b) => b.split(|&c| c == b'\n').filter(|l| !l.is_empty()).map(|l| l.to_vec()).collect(),
        Err(_) => vec![],
    }
}

fn passes<F: Fn(&[Vec<u8>]) -> i64>(f: &F, qs: &[Vec<u8>]) -> (f64, usize, i64) {
    let mut chk = f(qs);
    let t0 = Instant::now();
    let mut n = 0;
    loop {
        chk = f(qs);
        n += 1;
        if t0.elapsed().as_secs_f64() >= 1.0 { break; }
    }
    (t0.elapsed().as_secs_f64() / n as f64, n, chk)
}

fn concurrent<F: Fn(&[Vec<u8>]) -> i64 + Sync>(op: &str, qs: &[Vec<u8>], f: &F, threads: usize, unit: &str) {
    let t0 = Instant::now();
    let expect = f(qs);
    let single = t0.elapsed().as_secs_f64();
    let reps = std::cmp::min(400, (1.0 / single.max(1e-6)) as usize + 1);
    let barrier = std::sync::Barrier::new(threads + 1);
    let (res, el): (Vec<i64>, f64) = std::thread::scope(|s| {
        let hs: Vec<_> = (0..threads).map(|_| {
            let b = &barrier;
            s.spawn(move || { b.wait(); let mut r = 0; for _ in 0..reps { r = f(qs); } r })
        }).collect();
        barrier.wait();
        let start = Instant::now();
        let out: Vec<i64> = hs.into_iter().map(|h| h.join().unwrap()).collect();
        (out, start.elapsed().as_secs_f64())
    });
    let ok = res.iter().all(|&r| r == expect);
    emit(op, "throughput", threads as f64 * qs.len() as f64 * reps as f64 / el, unit,
         &format!("threads={} reps={} each thread runs the full query set; one shared index", threads, reps));
    emiti(op, "checksum_ok", ok as i64, "bool", "every thread's result equals the single-thread result");
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut text = fs::read(&a[1]).expect("read text");
    let qp = &a[2];
    let n = text.len();
    let mut syms: Vec<u8> = text.clone();
    syms.sort_unstable(); syms.dedup();
    text.push(b'$');
    syms.insert(0, b'$');
    let alphabet = Alphabet::new(&syms);
    let rss0 = status_kb("RssAnon:");

    let t0 = Instant::now();
    let sa = suffix_array(&text);
    let t_sa = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let bwt_v = bwt(&text, &sa);
    let less_v = less(&bwt_v, &alphabet);
    let occ = Occ::new(&bwt_v, 3, &alphabet);
    let t_fm = t1.elapsed().as_secs_f64();
    let fm = FMIndex::new(&bwt_v, &less_v, &occ);

    emiti("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build");
    emit("build", "time", t_sa + t_fm, "s", "suffix_array() [SA-IS] + bwt + less + Occ(k=3), single thread");
    emit("build", "time_suffix_array_only", t_sa, "s", "bio::data_structures::suffix_array::suffix_array");
    emit("build", "throughput", n as f64 / (t_sa + t_fm), "chars/s", "");
    let sa_bytes = (sa.len() * std::mem::size_of::<usize>()) as i64;
    let anon_delta = status_kb("RssAnon:") - rss0;
    let est = sa_bytes + n as i64 + (((n / 3) + 1) * syms.len() * std::mem::size_of::<usize>()) as i64;
    emiti("build", "index_bytes", est, "bytes", "analytic: SA Vec<usize> 8n + BWT n + Occ(k=3) (n/3+1)*sigma*8 ; excludes the n-byte input text");
    emit("build", "bytes_per_char", est as f64 / n as f64, "bytes/char", "analytic (see index_bytes); with k=32 the Occ part would be ~10x smaller");
    emiti("build", "rss_anon_growth", anon_delta, "bytes", "measured RssAnon growth over construction incl. freed-but-unreturned temporaries (upper estimate)");
    emiti("build", "sa_bytes_exact", sa_bytes, "bytes", "SA as Vec<usize> (8 B/char) - used for locate");
    emit("build", "bytes_per_char_incl_text", (est as f64 + n as f64) / n as f64, "bytes/char", "analytic + 1 B/char text");

    let present = read_lines(&format!("{}.present20.txt", qp));
    let absent = read_lines(&format!("{}.absent20.txt", qp));
    let loc8 = read_lines(&format!("{}.locate8.txt", qp));
    let loc4 = read_lines(&format!("{}.locate4.txt", qp));

    let count_fn = |qs: &[Vec<u8>]| -> i64 {
        let mut s = 0i64;
        for q in qs { if let Some(i) = iv(fm.backward_search(q.iter())) { s += i.upper as i64 - i.lower as i64; } }
        s
    };
    let contains_fn = |qs: &[Vec<u8>]| -> i64 {
        let mut s = 0i64;
        for q in qs { if iv(fm.backward_search(q.iter())).is_some() { s += 1; } }
        s
    };
    let locate_hits = std::sync::atomic::AtomicI64::new(0);
    let locate_fn = |qs: &[Vec<u8>]| -> i64 {
        let mut sum = 0i64; let mut h = 0i64;
        for q in qs {
            if let Some(i) = iv(fm.backward_search(q.iter())) { let pos = i.occ(&sa); h += pos.len() as i64; for p in pos { sum += p as i64; } }
        }
        locate_hits.store(h, std::sync::atomic::Ordering::Relaxed);
        sum
    };

    let rate = |op: &str, qs: &Vec<Vec<u8>>, f: &dyn Fn(&[Vec<u8>]) -> i64| {
        let (spp, np, chk) = { let g = |x: &[Vec<u8>]| f(x); passes(&g, qs) };
        emit(op, "throughput", qs.len() as f64 / spp, "queries/s", &format!("passes={}", np));
        emit(op, "time_per_pass", spp, "s", "");
        emiti(op, "checksum", chk, "count", "sum of results (cross-tool oracle)");
    };
    rate("count_present20", &present, &count_fn);
    rate("count_absent20", &absent, &count_fn);
    rate("contains_present20", &present, &contains_fn);
    rate("contains_absent20", &absent, &contains_fn);
    let loc_rate = |op: &str, qs: &Vec<Vec<u8>>| {
        let (spp, np, chk) = passes(&locate_fn, qs);
        let hits = locate_hits.load(std::sync::atomic::Ordering::Relaxed);
        emit(op, "throughput", qs.len() as f64 / spp, "patterns/s", &format!("passes={}", np));
        emit(op, "positions_per_s", hits as f64 / spp, "positions/s", "");
        emiti(op, "total_hits", hits, "count", "occurrences returned per pass (oracle)");
        emiti(op, "checksum", chk, "count", "sum of positions (oracle)");
    };
    loc_rate("locate8", &loc8);
    if !loc4.is_empty() { loc_rate("locate4", &loc4); }
    {
        let mut d: Vec<f64> = Vec::with_capacity(present.len());
        for q in &present { let a = Instant::now(); let _ = fm.backward_search(q.iter()); d.push(a.elapsed().as_nanos() as f64); }
        d.sort_by(|x, y| x.partial_cmp(y).unwrap());
        emit("count_present20", "latency_p50", d[d.len() / 2], "ns", "per-query Instant::now (incl. ~25ns overhead)");
        emit("count_present20", "latency_p99", d[(d.len() as f64 * 0.99) as usize], "ns", "");
        emit("count_present20", "latency_mean", d.iter().sum::<f64>() / d.len() as f64, "ns", "");
    }
    for th in [1usize, 4] {
        concurrent(&format!("count_present20_{}t", th), &present, &count_fn, th, "queries/s");
        let lf = |qs: &[Vec<u8>]| -> i64 {
            let mut sum = 0i64;
            for q in qs { if let Some(i) = iv(fm.backward_search(q.iter())) { for p in i.occ(&sa) { sum += p as i64; } } }
            sum
        };
        concurrent(&format!("locate8_{}t", th), &loc8, &lf, th, "patterns/s");
    }
    emiti("process", "peak_rss_hwm", status_kb("VmHWM:"), "bytes", "VmHWM at end of process");
}
