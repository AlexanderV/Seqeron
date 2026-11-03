#!/usr/bin/env python3
"""pydivsufsort (PyPI; libdivsufsort bindings) suffix-array harness.  Same R-line protocol as all harnesses.
usage: pydivsufsort_bench.py <text-file> <query-prefix>
count/contains/locate use pydivsufsort.sa_search (C) called from a Python loop, so q/s here is dominated by
Python call overhead (that is the realistic cost for a Python user)."""
import resource
import sys
import time

import numpy as np
import pydivsufsort


def emit(op, m, v, u, n=""):
    print(f"R\t{op}\t{m}\t{v!r}\t{u}\t{n}", flush=True)


def status_kb(key):
    for l in open("/proc/self/status"):
        if l.startswith(key):
            return int(l.split()[1]) * 1024
    return 0


def lines(p):
    try:
        return [l for l in open(p, "rb").read().split(b"\n") if l]
    except OSError:
        return []


def passes(fn, qs):
    chk = fn(qs)
    t0 = time.perf_counter()
    n = 0
    while True:
        chk = fn(qs)
        n += 1
        if time.perf_counter() - t0 >= 1.0:
            break
    return (time.perf_counter() - t0) / n, n, chk


def main():
    text_path, qp = sys.argv[1], sys.argv[2]
    T = np.fromfile(text_path, dtype=np.uint8)  # writable array (pydivsufsort rejects read-only buffers)
    n = len(T)
    anon0 = status_kb("RssAnon:")
    t0 = time.perf_counter()
    SA = pydivsufsort.divsufsort(T)
    el = time.perf_counter() - t0
    emit("build", "time", el, "s", "pydivsufsort.divsufsort(np.uint8 array), single thread")
    emit("build", "throughput", n / el, "chars/s")
    emit("build", "peak_rss_hwm_after_build", status_kb("VmHWM:"), "bytes", "VmHWM right after build (includes numpy text copy)")
    emit("build", "index_bytes", int(SA.nbytes), "bytes", f"SA dtype={SA.dtype} ({SA.dtype.itemsize} B/char); text extra")
    emit("build", "bytes_per_char", SA.nbytes / n, "bytes/char", "SA only")
    emit("build", "bytes_per_char_incl_text", (SA.nbytes + n) / n, "bytes/char")
    present, absent, loc8, loc4 = (lines(qp + s) for s in (".present20.txt", ".absent20.txt", ".locate8.txt", ".locate4.txt"))
    search = pydivsufsort.sa_search

    def count_fn(qs):
        s = 0
        for q in qs:
            s += search(T, SA, q)[0]
        return s

    def contains_fn(qs):
        s = 0
        for q in qs:
            s += search(T, SA, q)[0] > 0
        return s

    hits = [0]

    def locate_fn(qs):
        s = 0
        h = 0
        for q in qs:
            c, i = search(T, SA, q)
            if c:
                h += c
                s += int(SA[i:i + c].sum())
        hits[0] = h
        return s

    for op, qs, fn in (("count_present20", present, count_fn), ("count_absent20", absent, count_fn),
                       ("contains_present20", present, contains_fn), ("contains_absent20", absent, contains_fn)):
        spp, np_, chk = passes(fn, qs)
        emit(op, "throughput", len(qs) / spp, "queries/s", f"passes={np_}; Python loop over C sa_search")
        emit(op, "time_per_pass", spp, "s")
        emit(op, "checksum", int(chk), "count", "sum of results (cross-tool oracle)")
    for op, qs in (("locate8", loc8), ("locate4", loc4)):
        if not qs:
            continue
        spp, np_, chk = passes(locate_fn, qs)
        emit(op, "throughput", len(qs) / spp, "patterns/s", f"passes={np_}")
        emit(op, "positions_per_s", hits[0] / spp, "positions/s")
        emit(op, "total_hits", hits[0], "count", "occurrences returned per pass (oracle)")
        emit(op, "checksum", int(chk), "count", "sum of positions (oracle)")
    emit("process", "peak_rss_hwm", status_kb("VmHWM:"), "bytes", "VmHWM at end of process")
    _ = (resource, anon0)


if __name__ == "__main__":
    main()
