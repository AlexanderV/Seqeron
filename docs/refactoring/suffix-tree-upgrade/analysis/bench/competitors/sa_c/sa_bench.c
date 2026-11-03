/* Suffix-ARRAY benchmark harness for libsais and libdivsufsort (C, int32 SA, byte alphabet).
 * Same output protocol as every other harness:  R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes
 *   build with -DUSE_LIBSAIS   (+ -DLIBSAIS_OMP for the OpenMP build row)   or   -DUSE_DIVSUFSORT
 *   usage: sa_bench <text-file> <query-prefix>
 * count / contains / locate are plain binary searches over the suffix array (memcmp on the text):
 * a suffix ARRAY has no count/locate API of its own, so this is the standard minimal user code.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#ifdef USE_LIBSAIS
#include "libsais.h"
#endif
#ifdef USE_DIVSUFSORT
#include "divsufsort.h"
#endif

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static void emit(const char *op, const char *m, double v, const char *u, const char *n) { printf("R\t%s\t%s\t%.17g\t%s\t%s\n", op, m, v, u, n); }
static void emiti(const char *op, const char *m, long long v, const char *u, const char *n) { printf("R\t%s\t%s\t%lld\t%s\t%s\n", op, m, v, u, n); }

static long long status_kb(const char *key) {
  FILE *f = fopen("/proc/self/status", "r"); char line[256]; long long v = 0; size_t kl = strlen(key);
  if (!f) return 0;
  while (fgets(line, sizeof line, f)) if (!strncmp(line, key, kl)) { v = atoll(line + kl); break; }
  fclose(f); return v * 1024;
}

static uint8_t *T; static int32_t N; static int32_t *SA;

typedef struct { char **q; size_t n; } QS;

static QS read_lines(const char *path) {
  QS r = {0, 0};
  FILE *f = fopen(path, "rb"); if (!f) { r.q = NULL; return r; }
  fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  char *buf = malloc(sz + 1); if (fread(buf, 1, sz, f) != (size_t)sz) exit(3); buf[sz] = 0; fclose(f);
  size_t cap = 1024; r.q = malloc(cap * sizeof(char *));
  char *p = buf;
  while (*p) {
    char *e = strchr(p, '\n'); if (!e) e = p + strlen(p);
    if (e > p) { if (r.n == cap) { cap *= 2; r.q = realloc(r.q, cap * sizeof(char *)); } r.q[r.n++] = strndup(p, e - p); }
    p = *e ? e + 1 : e;
  }
  free(buf); return r;
}

/* first SA index whose suffix is >= pattern (compare on min(len) prefix) */
static inline int32_t lower_bound(const uint8_t *p, size_t m) {
  int32_t lo = 0, hi = N;
  while (lo < hi) {
    int32_t mid = lo + ((hi - lo) >> 1); size_t rem = (size_t)(N - SA[mid]); size_t l = rem < m ? rem : m;
    int c = memcmp(T + SA[mid], p, l);
    if (c < 0 || (c == 0 && rem < m)) lo = mid + 1; else hi = mid;
  }
  return lo;
}
/* first SA index whose suffix has prefix > pattern */
static inline int32_t upper_bound(const uint8_t *p, size_t m) {
  int32_t lo = 0, hi = N;
  while (lo < hi) {
    int32_t mid = lo + ((hi - lo) >> 1); size_t rem = (size_t)(N - SA[mid]); size_t l = rem < m ? rem : m;
    int c = memcmp(T + SA[mid], p, l);
    if (c <= 0) lo = mid + 1; else hi = mid;
  }
  return lo;
}
static inline long long count_one(const char *q) { size_t m = strlen(q); return (long long)(upper_bound((const uint8_t *)q, m) - lower_bound((const uint8_t *)q, m)); }
static inline int contains_one(const char *q) {
  size_t m = strlen(q); int32_t i = lower_bound((const uint8_t *)q, m);
  return i < N && (size_t)(N - SA[i]) >= m && !memcmp(T + SA[i], q, m);
}

static long long count_pass(QS *qs) { long long s = 0; for (size_t i = 0; i < qs->n; i++) s += count_one(qs->q[i]); return s; }
static long long contains_pass(QS *qs) { long long s = 0; for (size_t i = 0; i < qs->n; i++) s += contains_one(qs->q[i]); return s; }
static long long g_hits;
static long long locate_pass(QS *qs) {
  long long sum = 0, hits = 0;
  for (size_t i = 0; i < qs->n; i++) {
    const uint8_t *p = (const uint8_t *)qs->q[i]; size_t m = strlen(qs->q[i]);
    int32_t lo = lower_bound(p, m), hi = upper_bound(p, m);
    hits += hi - lo; for (int32_t k = lo; k < hi; k++) sum += SA[k];
  }
  g_hits = hits; return sum;
}

typedef long long (*PassFn)(QS *);
static void passes(PassFn fn, QS *qs, double *spp, int *np, long long *chk) {
  *chk = fn(qs);
  double t0 = now(); int n = 0;
  for (;;) { *chk = fn(qs); n++; if (now() - t0 >= 1.0) break; }
  *spp = (now() - t0) / n; *np = n;
}

static void rate(const char *op, QS *qs, PassFn fn) {
  double spp; int n; long long chk; char note[64];
  passes(fn, qs, &spp, &n, &chk); snprintf(note, sizeof note, "passes=%d", n);
  emit(op, "throughput", qs->n / spp, "queries/s", note);
  emit(op, "time_per_pass", spp, "s", "");
  emiti(op, "checksum", chk, "count", "sum of results (cross-tool oracle)");
}
static void loc_rate(const char *op, QS *qs) {
  double spp; int n; long long chk; char note[96];
  passes(locate_pass, qs, &spp, &n, &chk); long long hits = g_hits; snprintf(note, sizeof note, "passes=%d", n);
  emit(op, "throughput", qs->n / spp, "patterns/s", note);
  emit(op, "positions_per_s", hits / spp, "positions/s", "");
  emiti(op, "total_hits", hits, "count", "occurrences returned per pass (oracle)");
  emiti(op, "checksum", chk, "count", "sum of positions (oracle)");
}

typedef struct { QS *qs; PassFn fn; int reps; long long res; } TA;
static pthread_barrier_t bar;
static void *thr(void *a) { TA *t = a; long long s = 0; pthread_barrier_wait(&bar); for (int r = 0; r < t->reps; r++) s = t->fn(t->qs); t->res = s; return NULL; }
static long long locate_pass_ts(QS *qs) { /* thread-safe variant (no global) */
  long long sum = 0;
  for (size_t i = 0; i < qs->n; i++) {
    const uint8_t *p = (const uint8_t *)qs->q[i]; size_t m = strlen(qs->q[i]);
    int32_t lo = lower_bound(p, m), hi = upper_bound(p, m);
    for (int32_t k = lo; k < hi; k++) sum += SA[k];
  }
  return sum;
}
static void concurrent(const char *op, QS *qs, PassFn fn, int threads, const char *unit) {
  double t0 = now(); long long expect = fn(qs); double single = now() - t0;
  int reps = (int)(1.0 / (single > 1e-6 ? single : 1e-6)) + 1; if (reps > 400) reps = 400;
  pthread_t th[8]; TA ta[8]; pthread_barrier_init(&bar, NULL, threads + 1);
  for (int i = 0; i < threads; i++) { ta[i] = (TA){qs, fn, reps, 0}; pthread_create(&th[i], NULL, thr, &ta[i]); }
  pthread_barrier_wait(&bar); double t1 = now();
  for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
  double el = now() - t1; int ok = 1; for (int i = 0; i < threads; i++) if (ta[i].res != expect) ok = 0;
  char note[128]; snprintf(note, sizeof note, "threads=%d reps=%d each thread runs the full query set; one shared SA", threads, reps);
  emit(op, "throughput", (double)threads * qs->n * reps / el, unit, note);
  emiti(op, "checksum_ok", ok, "bool", "every thread's result equals the single-thread result");
  pthread_barrier_destroy(&bar);
}

static int cmp_d(const void *a, const void *b) { double x = *(double *)a, y = *(double *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: sa_bench text qprefix\n"); return 2; }
  FILE *f = fopen(argv[1], "rb"); if (!f) return 3;
  fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  T = malloc(sz); if (fread(T, 1, sz, f) != (size_t)sz) return 3; fclose(f); N = (int32_t)sz;
  SA = malloc(sizeof(int32_t) * (size_t)N);
  char path[1024];
#if defined(USE_LIBSAIS)
  const char *impl = "libsais";
  double t0 = now(); int32_t rc = libsais(T, SA, N, 0, NULL); double el = now() - t0;
#elif defined(USE_DIVSUFSORT)
  const char *impl = "libdivsufsort";
  double t0 = now(); int32_t rc = divsufsort(T, SA, N); double el = now() - t0;
#endif
  if (rc != 0) { fprintf(stderr, "SA construction failed rc=%d\n", rc); return 4; }
  char note[160]; snprintf(note, sizeof note, "%s single-thread SA construction, int32 SA, byte alphabet", impl);
  emiti("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build");
  emit("build", "time", el, "s", note);
  emit("build", "throughput", N / el, "chars/s", "");
  emiti("build", "index_bytes", 4LL * N, "bytes", "SA only (4 B/char); text (1 B/char) is extra");
  emit("build", "bytes_per_char", 4.0, "bytes/char", "SA only");
  emit("build", "bytes_per_char_incl_text", 5.0, "bytes/char", "SA + 1 B/char text");
#if defined(USE_LIBSAIS) && defined(LIBSAIS_OMP)
  {
    int32_t *SA2 = malloc(sizeof(int32_t) * (size_t)N);
    for (int th = 2; th <= 4; th *= 2) {
      double a = now(); libsais_omp(T, SA2, N, 0, NULL, th); double e = now() - a;
      char nm[32]; snprintf(nm, sizeof nm, "build_omp%d", th);
      emit(nm, "time", e, "s", "libsais_omp"); emiti(nm, "identical_to_serial", !memcmp(SA, SA2, sizeof(int32_t) * (size_t)N), "bool", "");
    }
    free(SA2);
  }
#endif
  snprintf(path, sizeof path, "%s.present20.txt", argv[2]); QS present = read_lines(path);
  snprintf(path, sizeof path, "%s.absent20.txt", argv[2]); QS absent = read_lines(path);
  snprintf(path, sizeof path, "%s.locate8.txt", argv[2]); QS loc = read_lines(path);
  snprintf(path, sizeof path, "%s.locate4.txt", argv[2]); QS loc4 = read_lines(path);
  rate("count_present20", &present, count_pass);
  rate("count_absent20", &absent, count_pass);
  rate("contains_present20", &present, contains_pass);
  rate("contains_absent20", &absent, contains_pass);
  loc_rate("locate8", &loc);
  if (loc4.q) loc_rate("locate4", &loc4);
  { /* latency */
    double *d = malloc(sizeof(double) * present.n); double sm = 0;
    for (size_t i = 0; i < present.n; i++) { double a = now(); (void)count_one(present.q[i]); d[i] = (now() - a) * 1e9; sm += d[i]; }
    qsort(d, present.n, sizeof(double), cmp_d);
    emit("count_present20", "latency_p50", d[present.n / 2], "ns", "per-query clock_gettime (incl. ~25ns timer overhead)");
    emit("count_present20", "latency_p99", d[(size_t)(present.n * 0.99)], "ns", "");
    emit("count_present20", "latency_mean", sm / present.n, "ns", ""); free(d);
  }
  for (int th = 1; th <= 4; th += 3) {
    char op[64];
    snprintf(op, sizeof op, "count_present20_%dt", th); concurrent(op, &present, count_pass, th, "queries/s");
    snprintf(op, sizeof op, "locate8_%dt", th); concurrent(op, &loc, locate_pass_ts, th, "patterns/s");
  }
  emiti("process", "peak_rss_hwm", status_kb("VmHWM:"), "bytes", "VmHWM at end of process");
  return 0;
}
