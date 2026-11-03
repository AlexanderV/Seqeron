#!/usr/bin/env python3
"""Python `suffix-trees` package (pure-Python Ukkonen suffix tree) harness, small n only.
usage: stree_bench.py <text-file> <query-prefix> <n-prefix>
Builds the tree over the first <n-prefix> characters of the text.  Queries are re-sampled DETERMINISTICALLY from that
prefix (random.Random('stree-<n>')) because the global present20 set is sampled from the full text: present = 20k
20-mers inside the prefix, absent = first 20k lines of the global absent20 file, locate = global locate8 set.
Same R-line protocol as all harnesses."""
import random
import sys
import time


def emit(op, m, v, u, n=""):
    print(f"R\t{op}\t{m}\t{v!r}\t{u}\t{n}", flush=True)


def status_kb(key):
    for l in open("/proc/self/status"):
        if l.startswith(key):
            return int(l.split()[1]) * 1024
    return 0


def lines(p):
    try:
        return [l.decode() for l in open(p, "rb").read().split(b"\n") if l]
    except OSError:
        return []


def main():
    text_path, qp, nprefix = sys.argv[1], sys.argv[2], int(sys.argv[3])
    from suffix_trees import STree
    text = open(text_path, "rb").read().decode()[:nprefix]
    n = len(text)
    r = random.Random(f"stree-{n}")
    present = [text[p:p + 20] for p in (r.randrange(n - 19) for _ in range(20000))]
    absent = lines(qp + ".absent20.txt")[:20000]
    loc8 = lines(qp + ".locate8.txt")
    t0 = time.perf_counter()
    st = STree.STree(text)
    el = time.perf_counter() - t0
    emit("build", "time", el, "s", "suffix_trees.STree.STree(text), pure Python Ukkonen")
    emit("build", "throughput", n / el, "chars/s")
    emit("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build")
    emit("build", "bytes_per_char_rss", (status_kb("VmHWM:")) / n, "bytes/char", "peak RSS / n (interpreter baseline included)")

    def timed(fn, qs):
        t = time.perf_counter()
        s = fn(qs)
        return time.perf_counter() - t, s

    def find_fn(qs):
        return sum(1 for q in qs if st.find(q) != -1)

    def findall_fn(qs):
        return sum(len(list(st.find_all(q))) for q in qs)

    e, s = timed(find_fn, present)
    emit("contains_present20", "throughput", len(present) / e, "queries/s", f"STree.find; {len(present)} queries from the prefix")
    emit("contains_present20", "checksum", s, "count")
    e, s = timed(find_fn, absent)
    emit("contains_absent20", "throughput", len(absent) / e, "queries/s", f"STree.find; {len(absent)} queries")
    emit("contains_absent20", "checksum", s, "count")
    e, s = timed(findall_fn, loc8)
    emit("locate8", "throughput", len(loc8) / e, "patterns/s", "STree.find_all")
    emit("locate8", "total_hits", s, "count", "occurrences returned (prefix text only)")
    emit("process", "peak_rss_hwm", status_kb("VmHWM:"), "bytes", "VmHWM at end of process")


if __name__ == "__main__":
    main()
