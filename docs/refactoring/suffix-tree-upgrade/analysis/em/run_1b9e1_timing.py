#!/usr/bin/env python3
"""What the normalisation of the reopened container costs, next to what a resume saves (work package 1b.9e.1).

For every marker the build is killed there once per repetition (`--stop-at`, then SIGKILL; a marker written `name+0.5` is not a stop: the build runs on for 0.5 s after the marker line and is then killed, which leaves the partial output of a
phase in the file), the dirty pages are written out (sync), the work directory is copied, and then the SAME state is finished in separate processes, each from its own copy
(the order rotates with the repetition):
  new       the resume of this package (the normalisation of 1b.9e.1: every byte outside the ranges of the kept stages is read, the ones that are not zero are written as zeros)
  old       the resume of the commit before it (`--bench-old`: no normalisation; right on a state that nothing damaged, which is the case here)
  restart   the build that the package 1b.9e replaced (`--resume false`: the temporary file and the scratch are deleted and everything is built again)
  garbage   the resume of this package, with every byte outside the kept ranges overwritten with 0xA5 first: the worst case of the normalisation (every unit has to be written)
Seconds are the build's own (`total_s`); `dev_mb` is what the block device wrote from the start of the process until the end of a `sync` after it. The bytes of every result are compared with the uninterrupted build.
usage: run_1b9e1_timing.py --bench-new <dll> --bench-old <dll> --input <text> --kind <label> --work <directory> --rows <tsv> [--reps 3] [--threads 4] [--markers a,b,c]
Standard library only."""
import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import time

HEX = re.compile(r"\b[0-9a-f]{64}\b")
MARKERS = ["packed", "rank-written", "ph23", "sa-built", "b2", "lcp-child-built", "bwt-built", "info-written"]


def cmd(bench, inp, work, threads, resume, stop_at=None):
    c = ["dotnet", bench, "sesa", "resume", "build", "--input", inp, "--work", work, "--threads", str(threads), "--resume", resume, "--psi", "true", "--fsync", "true",
         "--markers", "true", "--hash", "true", "--keep", "true"]
    if stop_at:
        c += ["--stop-at", stop_at]
    return c


def rows_of(text):
    d = {}
    for line in text.splitlines():
        p = line.split("\t")
        if len(p) >= 5 and p[0] == "R" and p[1] == "resume":
            d[p[2]] = p[3]
    return d


def clean(path):
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)


def kill_at(bench, inp, work, threads, marker):
    """Blocks the build at `marker` (--stop-at) and kills it; with `marker+seconds` lets it run on that long after the marker line, then kills it. False when the marker was not reached."""
    name, _, delay = marker.partition("+")
    p = subprocess.Popen(cmd(bench, inp, work, threads, "true", None if delay else name), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL, text=True)
    reached = False
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] == "M" and f[1] == name:
            reached = True
            break
    time.sleep(float(delay) if delay else 0.05)
    landed = p.poll() is None
    if landed:
        p.kill()
    p.wait()
    return reached and landed and not os.path.exists(os.path.join(work, "out.sesa")) and os.path.exists(os.path.join(work, "out.sesa.tmp"))


def kept_complement(work):
    st = json.loads(open(os.path.join(work, "out.sesa.tmp.d", "state.json"), "rb").read())
    kept = sorted((f["offset"], f["length"]) for s in st["stages"] for f in s["files"] if f.get("container"))
    length = os.path.getsize(os.path.join(work, "out.sesa.tmp"))
    gaps, at = [], 0
    for off, ln in kept:
        if off > at:
            gaps.append((at, min(off, length) - at))
        at = max(at, off + ln)
    if at < length:
        gaps.append((at, length - at))
    return [(o, l) for o, l in gaps if l > 0]


def fill_garbage(work):
    block = b"\xa5" * (1 << 22)
    with open(os.path.join(work, "out.sesa.tmp"), "r+b") as f:
        for off, ln in kept_complement(work):
            f.seek(off)
            while ln > 0:
                k = min(ln, len(block))
                f.write(block[:k])
                ln -= k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-new", required=True)
    ap.add_argument("--bench-old", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--kind", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--markers", default=",".join(MARKERS))
    a = ap.parse_args()
    markers = a.markers.split(",")
    base = a.work
    victim_dir, snap_dir, run_dir = base + "_victim", base + "_snap", base + "_run"
    clean(victim_dir)
    os.sync()
    p = subprocess.run(cmd(a.bench_new, a.input, victim_dir, a.threads, "false"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, stdin=subprocess.DEVNULL)
    ref = rows_of(p.stdout)
    ref_hash = ref["sha256"]
    print(f"# {a.kind}: uninterrupted build without Resume {float(ref['total_s']):.2f} s, sha256 {ref_hash[:16]}")
    out = open(a.rows, "a")
    results = {}
    good = True
    modes = ["new", "old", "restart", "garbage"]
    for marker in markers:
        for rep in range(a.reps):
            clean(victim_dir)
            os.sync()
            if not kill_at(a.bench_new, a.input, victim_dir, a.threads, marker):
                print("marker not reached", marker)
                continue
            os.sync()
            shutil.rmtree(snap_dir, ignore_errors=True)
            shutil.copytree(victim_dir, snap_dir, symlinks=True)
            shutil.rmtree(victim_dir, ignore_errors=True)
            order = modes[rep % len(modes):] + modes[:rep % len(modes)]
            for mode in order:
                shutil.rmtree(run_dir, ignore_errors=True)
                shutil.copytree(snap_dir, run_dir, symlinks=True)
                if mode == "garbage":
                    fill_garbage(run_dir)
                os.sync()
                time.sleep(0.3)
                bench = a.bench_old if mode == "old" else a.bench_new
                t0 = time.monotonic()
                p = subprocess.run(cmd(bench, a.input, run_dir, a.threads, "false" if mode == "restart" else "true"), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, stdin=subprocess.DEVNULL)
                wall = time.monotonic() - t0
                r = rows_of(p.stdout)
                left = sorted(os.listdir(run_dir))
                ok = p.returncode == 0 and r.get("sha256") == ref_hash and (left == ["out.sesa"] or (mode == "restart" and left == ["out.sesa", "out.sesa.lock"]))
                good &= ok
                dev_mb = float(r.get("dev_sectors_written_synced", "0")) * 512 / 1e6
                row = [a.kind, marker, mode, rep, r.get("total_s", ""), f"{wall:.3f}", f"{dev_mb:.0f}", r.get("verify_s", ""), r.get("clear_s", ""), r.get("clear_scanned_bytes", ""), r.get("clear_bytes", ""),
                       int(ok), r.get("skipped", "")]
                out.write("\t".join(str(x) for x in row) + "\n")
                out.flush()
                results.setdefault((marker, mode), []).append(row)
                print(f"{'ok' if ok else 'FAILED'} {a.kind} {marker} {mode} {rep}: {r.get('total_s')} s, dev {dev_mb:.0f} MB, clear {r.get('clear_s', '')} s ({r.get('clear_bytes', '')} of {r.get('clear_scanned_bytes', '')} bytes), skipped [{r.get('skipped', '')}]", flush=True)
                shutil.rmtree(run_dir, ignore_errors=True)
    shutil.rmtree(snap_dir, ignore_errors=True)

    def med(rows, i):
        v = [float(r[i]) for r in rows if r[i] != ""]
        return statistics.median(v) if v else float("nan")

    print("| input | killed at | restart s | old resume s | new resume s | garbage resume s | normalisation s (new) | scanned MB | written MB (new) | written MB (garbage) | device MB restart / old / new / garbage |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for marker in markers:
        if (marker, "new") not in results:
            continue
        rs, ol, nw, gb = (results[(marker, m)] for m in ("restart", "old", "new", "garbage"))
        print(f"| {a.kind} | {marker} | {med(rs, 4):.2f} | {med(ol, 4):.2f} | {med(nw, 4):.2f} | {med(gb, 4):.2f} | {med(nw, 8):.3f} | {med(nw, 9) / 1e6:.0f} | {med(nw, 10) / 1e6:.1f} | {med(gb, 10) / 1e6:.0f} | "
              f"{med(rs, 6):.0f} / {med(ol, 6):.0f} / {med(nw, 6):.0f} / {med(gb, 6):.0f} |")
    print("ALL OK" if good else "SOME FAILED")
    raise SystemExit(0 if good else 1)


if __name__ == "__main__":
    main()
