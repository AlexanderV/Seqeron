// sdsl-lite benchmark harness: compressed suffix arrays / FM-indexes (csa_*) and a compressed suffix tree (cst_sct3).
// Same output protocol as every other harness:  R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes
//   usage: sdsl_bench <config> <text-file> <query-prefix> [lrs]
//   configs: csa_wt_huff   = csa_wt<wt_huff<>,32,64>                     (FM-index, plain bitvectors, fast)
//            csa_wt_rrr    = csa_wt<wt_huff<rrr_vector<127>>,512,1024>    (FM-index, RRR bitvectors, small)
//            csa_sada      = csa_sada<>                                   (Psi-array compressed SA)
//            csa_bc        = csa_bitcompressed<>                          (plain SA, ceil(log2 n) bits/entry)
//            cst_sct3      = cst_sct3<>                                   (compressed suffix tree)
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <thread>
#include <atomic>
#include <algorithm>
#include <string>
#include <vector>
using namespace sdsl;
using Clock = std::chrono::steady_clock;
using QV = std::vector<std::string>;
using PassFn = std::function<long long(const QV &)>;

static double nowS() { return std::chrono::duration<double>(Clock::now().time_since_epoch()).count(); }
static void emit(const char *op, const char *m, double v, const char *u, const std::string &n = "") { printf("R\t%s\t%s\t%.17g\t%s\t%s\n", op, m, v, u, n.c_str()); }
static void emiti(const char *op, const char *m, long long v, const char *u, const std::string &n = "") { printf("R\t%s\t%s\t%lld\t%s\t%s\n", op, m, v, u, n.c_str()); }
static long long status_kb(const char *key) {
  std::ifstream f("/proc/self/status"); std::string l; size_t kl = strlen(key);
  while (std::getline(f, l)) if (l.compare(0, kl, key) == 0) return atoll(l.c_str() + kl) * 1024;
  return 0;
}
static QV read_lines(const std::string &p) {
  QV v; std::ifstream f(p); std::string l;
  while (std::getline(f, l)) if (!l.empty()) v.push_back(l);
  return v;
}
static std::string read_file(const std::string &p) { std::ifstream f(p, std::ios::binary); std::stringstream ss; ss << f.rdbuf(); return ss.str(); }

static long long g_hits = 0;
static void passes(const PassFn &fn, const QV &qs, double &spp, int &np, long long &chk) {
  chk = fn(qs); double t0 = nowS(); int n = 0;
  for (;;) { chk = fn(qs); n++; if (nowS() - t0 >= 1.0) break; }
  spp = (nowS() - t0) / n; np = n;
}
static void rate(const char *op, const QV &qs, const PassFn &fn) {
  double spp; int n; long long chk; passes(fn, qs, spp, n, chk);
  emit(op, "throughput", qs.size() / spp, "queries/s", "passes=" + std::to_string(n));
  emit(op, "time_per_pass", spp, "s");
  emiti(op, "checksum", chk, "count", "sum of results (cross-tool oracle)");
}
static bool concurrent(const char *op, const QV &qs, const PassFn &fn, int threads, const char *unit) {
  double t0 = nowS(); long long expect = fn(qs); double single = nowS() - t0;
  int reps = std::min(400, (int)(1.0 / std::max(single, 1e-6)) + 1);
  std::vector<long long> res(threads, 0); std::atomic<int> ready(0); std::atomic<bool> go(false);
  std::vector<std::thread> th;
  for (int i = 0; i < threads; i++) th.emplace_back([&, i] { ready++; while (!go) {} long long s = 0; for (int r = 0; r < reps; r++) s = fn(qs); res[i] = s; });
  while (ready < threads) {}
  double t1 = nowS(); go = true; for (auto &t : th) t.join(); double el = nowS() - t1;
  bool ok = true; for (auto r : res) if (r != expect) ok = false;
  emit(op, "throughput", (double)threads * qs.size() * reps / el, unit, "threads=" + std::to_string(threads) + " reps=" + std::to_string(reps) + " each thread runs the full query set; one shared index");
  emiti(op, "checksum_ok", ok, "bool", "every thread's result equals the single-thread result");
  return ok;
}

struct Sets { QV present, absent, loc8, loc4; };

static void query_suite(const char *name, const Sets &S, const PassFn &count_fn, const PassFn &contains_fn, const PassFn &locate_fn, std::atomic<long long> &hits_out) {
  (void)name;
  rate("count_present20", S.present, count_fn);
  rate("count_absent20", S.absent, count_fn);
  rate("contains_present20", S.present, contains_fn);
  rate("contains_absent20", S.absent, contains_fn);
  auto loc_rate = [&](const char *op, const QV &qs) {
    double spp; int n; long long chk; passes(locate_fn, qs, spp, n, chk); long long hits = hits_out.load();
    emit(op, "throughput", qs.size() / spp, "patterns/s", "passes=" + std::to_string(n));
    emit(op, "positions_per_s", hits / spp, "positions/s");
    emiti(op, "total_hits", hits, "count", "occurrences returned per pass (oracle)");
    emiti(op, "checksum", chk, "count", "sum of positions (oracle)");
  };
  loc_rate("locate8", S.loc8);
  if (!S.loc4.empty()) loc_rate("locate4", S.loc4);
  { std::vector<double> d; double sm = 0;
    for (auto &q : S.present) { double a = nowS(); count_fn(QV{q}); double e = (nowS() - a) * 1e9; d.push_back(e); sm += e; }
    std::sort(d.begin(), d.end());
    emit("count_present20", "latency_p50", d[d.size() / 2], "ns", "per-query (incl. ~100ns QV construction+timer overhead)");
    emit("count_present20", "latency_p99", d[(size_t)(d.size() * 0.99)], "ns"); emit("count_present20", "latency_mean", sm / d.size(), "ns"); }
  for (int th : {1, 4}) {
    bool ok = concurrent(("count_present20_" + std::to_string(th) + "t").c_str(), S.present, count_fn, th, "queries/s");
    if (!ok) { emit("locate8_4t", "skipped", 1, "bool", "concurrent reads of this sdsl structure returned wrong results (not thread-safe); further concurrent tests skipped to avoid a crash"); break; }
    concurrent(("locate8_" + std::to_string(th) + "t").c_str(), S.loc8, [&](const QV &q) { return locate_fn(q); }, th, "patterns/s");
  }
}

template <class CSA>
int run_csa(const std::string &textfile, const Sets &S) {
  std::string text = read_file(textfile);
  CSA csa;
  double t0 = nowS();
  construct_im(csa, text.c_str(), 1);
  double el = nowS() - t0;
  emiti("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build (includes construction temporaries and the 1 B/char input)");
  emit("build", "time", el, "s", "sdsl::construct_im (libdivsufsort inside), single-thread");
  emit("build", "throughput", text.size() / el, "chars/s");
  long long bytes = (long long)size_in_bytes(csa);
  emiti("build", "index_bytes", bytes, "bytes", "sdsl::size_in_bytes(csa); includes sampled SA + sampled ISA (text extraction) ; no separate text copy needed");
  emit("build", "bytes_per_char", (double)bytes / text.size(), "bytes/char");
  std::atomic<long long> hits(0);
  auto count_fn = [&](const QV &qs) { long long s = 0; for (auto &q : qs) s += (long long)count(csa, (const unsigned char *)q.data(), (const unsigned char *)q.data() + q.size()); return s; };
  auto contains_fn = [&](const QV &qs) { long long s = 0; for (auto &q : qs) s += count(csa, (const unsigned char *)q.data(), (const unsigned char *)q.data() + q.size()) > 0; return s; };
  auto locate_fn = [&](const QV &qs) { long long sum = 0, h = 0; for (auto &q : qs) { auto r = locate(csa, (const unsigned char *)q.data(), (const unsigned char *)q.data() + q.size()); h += r.size(); for (size_t i = 0; i < r.size(); i++) sum += r[i]; } hits = h; return sum; };
  query_suite("csa", S, count_fn, contains_fn, locate_fn, hits);
  return 0;
}

// CST: descend with child(v,c) and edge(v,k); count = size(v) (leaf count), locate = csa[lb..rb]
using CST = cst_sct3<>;
static bool cst_find(const CST &cst, const std::string &p, CST::node_type &out) {
  auto v = cst.root(); size_t i = 0, m = p.size();
  while (i < m) {
    auto w = cst.child(v, (CST::char_type)(unsigned char)p[i]);
    if (w == cst.root()) return false;
    size_t dw = cst.depth(w), lim = std::min(dw, m);
    for (size_t k = i + 2; k <= lim; ++k) if (cst.edge(w, k) != (CST::char_type)(unsigned char)p[k - 1]) return false;
    i = lim; v = w;
  }
  out = v; return true;
}
int run_cst(const std::string &textfile, const Sets &S, bool lrs, const std::string &qprefix) {
  std::string text = read_file(textfile);
  CST cst; double t0 = nowS();
  construct_im(cst, text.c_str(), 1);
  double el = nowS() - t0;
  emiti("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build");
  emit("build", "time", el, "s", "sdsl::construct_im(cst_sct3<>), single-thread");
  emit("build", "throughput", text.size() / el, "chars/s");
  long long bytes = (long long)size_in_bytes(cst);
  emiti("build", "index_bytes", bytes, "bytes", "sdsl::size_in_bytes(cst_sct3<>) (csa_sada + lcp_dac + bp)");
  emit("build", "bytes_per_char", (double)bytes / text.size(), "bytes/char");
  std::atomic<long long> hits(0);
  auto count_fn = [&](const QV &qs) { long long s = 0; for (auto &q : qs) { CST::node_type v; if (cst_find(cst, q, v)) s += cst.size(v); } return s; };
  auto contains_fn = [&](const QV &qs) { long long s = 0; for (auto &q : qs) { CST::node_type v; s += cst_find(cst, q, v); } return s; };
  auto locate_fn = [&](const QV &qs) { long long sum = 0, h = 0; for (auto &q : qs) { CST::node_type v; if (cst_find(cst, q, v)) { for (auto i = cst.lb(v); i <= cst.rb(v); ++i) { sum += cst.csa[i]; h++; } } } hits = h; return sum; };
  query_suite("cst", S, count_fn, contains_fn, locate_fn, hits);
  if (lrs) {
    double a = nowS(); size_t best = 0; size_t nodes = 0;
    for (auto it = cst.begin(); it != cst.end(); ++it) { if (it.visit() == 1 && !cst.is_leaf(*it)) { size_t d = cst.depth(*it); if (d > best) best = d; nodes++; } }
    emit("lrs", "first_call_time", nowS() - a, "s", "full DFS over internal nodes taking max string depth (user code on cst_sct3)");
    emiti("lrs", "length", best, "chars"); emiti("lrs", "internal_nodes", nodes, "count");
  }
  (void)qprefix;
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: sdsl_bench config text qprefix [lrs]\n"); return 2; }
  std::string cfg = argv[1], textfile = argv[2], qp = argv[3]; bool lrs = argc > 4 && std::string(argv[4]) == "lrs";
  Sets S{read_lines(qp + ".present20.txt"), read_lines(qp + ".absent20.txt"), read_lines(qp + ".locate8.txt"), read_lines(qp + ".locate4.txt")};
  int rc;
  if (cfg == "csa_wt_huff") rc = run_csa<csa_wt<wt_huff<>, 32, 64>>(textfile, S);
  else if (cfg == "csa_wt_rrr") rc = run_csa<csa_wt<wt_huff<rrr_vector<127>>, 512, 1024>>(textfile, S);
  else if (cfg == "csa_sada") rc = run_csa<csa_sada<>>(textfile, S);
  else if (cfg == "csa_bc") rc = run_csa<csa_bitcompressed<>>(textfile, S);
  else if (cfg == "cst_sct3") rc = run_cst(textfile, S, lrs, qp);
  else { fprintf(stderr, "unknown config\n"); return 2; }
  emiti("process", "peak_rss_hwm", status_kb("VmHWM:"), "bytes", "VmHWM at end of process");
  return rc;
}
