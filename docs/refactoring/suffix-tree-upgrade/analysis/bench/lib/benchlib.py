#!/usr/bin/env python3
"""Shared measurement driver.  Pure Python 3 stdlib.

Every benchmark subject (C# harness, Go harness, C/C++/Rust/Python harnesses of the hand-off tools)
prints result lines on stdout:   R<TAB>operation<TAB>metric<TAB>value<TAB>unit<TAB>notes
This driver runs a subject several times, adds externally measured process metrics
(wall time, peak RSS via wait4/ru_maxrss, user+sys CPU), appends one raw row per value/repeat to a
raw file and can aggregate (median over repeats) into the final bench.tsv.

CLI
  benchlib.py measure --raw RAW --tool T --version V --input NAME [--repeats 3] [--drop-caches]
                      [--op-prefix cold_] [--timeout S] [--notes TXT] [--data-dir DIR] -- CMD ARGS...
  benchlib.py aggregate --raw RAW --out bench.tsv

bench.tsv columns: tool version input n operation metric value unit repeats notes
"""
import argparse
import os
import statistics
import subprocess
import sys
import time

HEADER = ["tool", "version", "input", "n", "operation", "metric", "value", "unit", "repeats", "notes"]
RAW_HEADER = ["tool", "version", "input", "n", "operation", "metric", "value", "unit", "repeat", "notes"]
HERE = os.path.dirname(os.path.abspath(__file__))


def clean(s):
    return str(s).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def manifest_n(data_dir, name):
    with open(os.path.join(data_dir, "MANIFEST.tsv")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == name:
                return int(p[1])
    return -1


def drop_caches():
    os.sync()
    try:
        with open("/proc/sys/vm/drop_caches", "w") as f:
            f.write("3\n")
        return True
    except OSError:
        return False


def run_once(cmd, timeout=None, env=None):
    t0 = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    out_chunks = []
    timed_out = []
    timer = None
    if timeout:
        import threading as _t
        def _kill():
            timed_out.append(1)
            try:
                p.kill()
            except OSError:
                pass
        timer = _t.Timer(timeout, _kill)
        timer.start()
    try:
        # read stdout to EOF (stderr drained by a thread would be safer; harnesses write little to stderr)
        import threading
        err = []
        th = threading.Thread(target=lambda: err.append(p.stderr.read()))
        th.start()
        out_chunks.append(p.stdout.read())
        th.join()
        _, status, ru = os.wait4(p.pid, 0)
    except Exception:
        p.kill()
        raise
    if timer:
        timer.cancel()
    wall = time.time() - t0
    rc = os.waitstatus_to_exitcode(status)
    if timed_out:
        rc = 124
        err = [("TIMEOUT after %.0fs (killed)" % timeout).encode()]
    out = b"".join(out_chunks).decode(errors="replace")
    rows = []
    for line in out.splitlines():
        if line.startswith("R\t"):
            f = line.split("\t")
            if len(f) >= 5:
                rows.append((f[1], f[2], f[3], f[4], f[5] if len(f) > 5 else ""))
    ext = [
        ("process", "ext_peak_rss", str(ru.ru_maxrss * 1024), "bytes", "wait4 ru_maxrss of the whole process (independent of harness)"),
        ("process", "ext_wall", repr(wall), "s", "process wall time incl. startup, text read, queries"),
        ("process", "ext_cpu_user_sys", repr(ru.ru_utime + ru.ru_stime), "s", "user+sys CPU seconds"),
    ]
    return rc, rows + ext, (b"".join(err)).decode(errors="replace")


def cmd_measure(a):
    n = a.n if a.n else manifest_n(a.data_dir, a.input)
    new = not os.path.exists(a.raw) or os.path.getsize(a.raw) == 0
    ok_runs = 0
    with open(a.raw, "a") as raw:
        if new:
            raw.write("\t".join(RAW_HEADER) + "\n")
        for rep in range(1, a.repeats + 1):
            if a.drop_caches and not drop_caches():
                print("WARN: could not drop caches", file=sys.stderr)
            rc, rows, err = run_once(a.cmd, a.timeout)
            if rc != 0:
                print(f"[{a.tool} {a.input}] repeat {rep} FAILED rc={rc}: {err[-400:]}", file=sys.stderr)
                raw.write("\t".join(map(clean, [a.tool, a.version, a.input, n, a.op_prefix + "run", "failed_rc", rc, "code", rep,
                                              "BLOCKED: " + err[-200:]])) + "\n")
                continue
            ok_runs += 1
            for op, metric, val, unit, notes in rows:
                opn = (a.phase + "_process") if op == "process" else a.op_prefix + op
                raw.write("\t".join(map(clean, [a.tool, a.version, a.input, n, opn, metric, val, unit, rep,
                                              (notes + " " + a.notes).strip()])) + "\n")
            raw.flush()
            print(f"[{a.tool} {a.input} {a.op_prefix}] repeat {rep}/{a.repeats} ok", file=sys.stderr, flush=True)
    return 0 if ok_runs else 1


def fmt(x):
    if isinstance(x, float):
        return repr(round(x, 6)) if abs(x) < 1e6 else repr(round(x, 1))
    return str(x)


def cmd_aggregate(a):
    groups = {}
    order = []
    with open(a.raw) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            key = tuple(p[:6]) + (p[7],)
            if key not in groups:
                groups[key] = {"vals": [], "notes": p[9], "reps": set()}
                order.append(key)
            groups[key]["vals"].append(p[6])
            groups[key]["reps"].add(p[8])
    with open(a.out, "w") as o:
        o.write("\t".join(HEADER) + "\n")
        for key in order:
            g = groups[key]
            tool, ver, inp, n, op, metric, unit = key
            vals = g["vals"]
            notes = g["notes"]
            try:
                nums = [float(v) for v in vals]
                val = statistics.median(nums)
                if all(x == int(x) for x in nums) and abs(val) < 2 ** 62:
                    val = int(val) if val == int(val) else val
                if len(set(nums)) > 1 and metric in ("checksum", "total_hits", "matches", "length", "node_count", "leaf_count", "serialized_bytes", "index_bytes_file"):
                    notes += f" [WARN values differ across repeats: {sorted(set(nums))[:3]}]"
                if len(nums) > 1:
                    notes += f" [min={fmt(min(nums))} max={fmt(max(nums))}]"
                vs = fmt(val) if not isinstance(val, int) else str(val)
            except ValueError:
                vs = vals[0]
            o.write("\t".join(map(clean, [tool, ver, inp, n, op, metric, vs, unit, len(g["reps"]) if g["reps"] else len(vals), notes.strip()])) + "\n")
    print("wrote", a.out)
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="c", required=True)
    m = sub.add_parser("measure")
    m.add_argument("--raw", required=True)
    m.add_argument("--tool", required=True)
    m.add_argument("--version", required=True)
    m.add_argument("--input", required=True)
    m.add_argument("--repeats", type=int, default=3)
    m.add_argument("--drop-caches", action="store_true")
    m.add_argument("--op-prefix", default="")
    m.add_argument("--phase", default="run", help="names the per-process rows: <phase>_process (peak RSS, wall, cpu)")
    m.add_argument("--timeout", type=float, default=0)
    m.add_argument("--notes", default="")
    m.add_argument("--n", type=int, default=0, help="override n (prefix runs)")
    m.add_argument("--data-dir", default=os.path.join(HERE, "..", "data", ".."))
    m.add_argument("cmd", nargs=argparse.REMAINDER)
    g = sub.add_parser("aggregate")
    g.add_argument("--raw", required=True)
    g.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.c == "measure":
        if a.cmd and a.cmd[0] == "--":
            a.cmd = a.cmd[1:]
        sys.exit(cmd_measure(a))
    sys.exit(cmd_aggregate(a))


if __name__ == "__main__":
    main()
