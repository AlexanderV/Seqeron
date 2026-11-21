#!/usr/bin/env python3
"""The time that a resume saves (work package 1b.9e): for one boundary of every phase the build is stopped there and killed (`kill -9`), the dirty pages that it left are written out (sync), and then the work is done
twice in two separate trials: once as a resume (`--resume true`) and once as the restart that the build had before the package (`--resume false`: the stale temporary file and scratch are deleted and everything is built
again). Both run with the same number of threads (no variation: this is a measurement of seconds), `--reps` times each, the trials alternating. The bytes of every result are compared with the reference.
usage: run_1b9e_timing.py --bench <dll> --input <text> --kind <label> --work <directory> --rows <tsv> [--reps 3] [--threads 4]
Rows: kind, marker (where the build was killed), mode (resume | restart), repetition, seconds inside Build, wall seconds of the process, hash equal."""
import argparse
import os
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_1b9e_kill as k  # the process wrapper and the helpers of the kill driver

MARKERS = [("ph0", "created"), ("ph1", "rank-written"), ("ph23", "ph23"), ("ph4", "b2"), ("bwt", "lcp-child-built"), ("psi", "bwt-built"), ("ph7", "info-written")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--kind", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    target = os.path.join(a.work, "out.sesa")

    def args(resume, extra=()):
        return ["--input", a.input, "--work", a.work, "--resume", resume, "--threads", str(a.threads), "--psi", "true", "--fsync", "true", "--markers", "true", "--hash", "true", "--keep", "true", "--label", a.kind] + list(extra)

    k.clean(a.work)
    os.sync()
    ref = k.Run(a.bench, args("false"))
    assert ref.finish() == 0
    rr = ref.rows()
    ref_hash, fresh = rr["sha256"], float(rr["total_s"])
    k.clean(a.work)
    out = open(a.rows, "a")
    results = {}
    good = True
    for phase, marker in MARKERS:
        for rep in range(a.reps):
            for mode in ("true", "false"):
                k.clean(a.work)
                os.sync()
                victim = k.Run(a.bench, args("true", ["--stop-at", marker]))
                assert victim.wait_mark(marker, 3600) is not None
                time.sleep(0.05)
                victim.kill()
                os.sync()
                time.sleep(0.3)
                t0 = time.monotonic()
                fin = k.Run(a.bench, args(mode))
                rc = fin.finish()
                wall = time.monotonic() - t0
                fr = fin.rows()
                # a restart without Resume does not know the lock file that the killed build left (an empty file named <target>.lock): it is not a leftover of the restart, and it is allowed there
                left = k.entries(a.work)
                ok = rc == 0 and fr.get("sha256") == ref_hash and (left == ["out.sesa"] or (mode == "false" and left == ["out.sesa", "out.sesa.lock"]))
                good &= ok
                name = "resume" if mode == "true" else "restart"
                out.write("\t".join(str(x) for x in [a.kind, marker, name, rep, fr.get("total_s", ""), f"{wall:.3f}", int(ok), fr.get("skipped", "")]) + "\n")
                out.flush()
                results.setdefault((phase, marker, name), []).append(float(fr["total_s"]))
                print(f"{'ok' if ok else 'FAILED'} {a.kind} killed at {marker} {name} {rep}: {fr.get('total_s')} s, skipped [{fr.get('skipped', '')}]", flush=True)
                k.clean(a.work)
    print(f"# {a.kind} uninterrupted build without Resume (before any kill): {fresh:.2f} s")
    print("| input | killed at (phase) | restart s (median [min-max]) | resume s (median [min-max]) | saved s | saved % of the restart |")
    print("| --- | --- | --- | --- | --- | --- |")
    for phase, marker in MARKERS:
        rs, rm = sorted(results[(phase, marker, "restart")]), sorted(results[(phase, marker, "resume")])
        saved = statistics.median(rs) - statistics.median(rm)
        print(f"| {a.kind} | {marker} ({phase}) | {statistics.median(rs):.2f} [{rs[0]:.2f}-{rs[-1]:.2f}] | {statistics.median(rm):.2f} [{rm[0]:.2f}-{rm[-1]:.2f}] | {saved:.2f} | {100 * saved / statistics.median(rs):.0f} % |")
    out.close()
    print("ALL OK" if good else "SOME FAILED")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
