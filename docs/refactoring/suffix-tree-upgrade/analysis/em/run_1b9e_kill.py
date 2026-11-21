#!/usr/bin/env python3
"""Kill campaign of work package 1b.9e (ADR-004 section 5a "1b.9e"): `kill -9` of a separate process that builds a container with Resume on, at every boundary of the build and at random moments
inside every phase, then a resume in a NEW process, and the check that the container is byte-identical (SHA-256) to the uninterrupted build of the same input, that nothing but the target is left in the
directory, and that a kill during the publication leaves either no target or a complete one.

usage: run_1b9e_kill.py --bench <SuffixTree.Esa.Bench.dll> --input <text> --kind <label> --work <directory> --rows <tsv>
                        [--per-phase 3] [--seed 1] [--boundaries yes|no] [--threads 4] [--chain 0.25] [--psi true] [--lock-check yes]

One row per kill (tab separated): kind, n, id, mode (boundary | random), aimed (marker or seconds), last marker seen, phase, killed (1 = the signal landed, 0 = the build finished first), target after the kill
(none | complete | damaged), the stages that the resume skipped, seconds of the resume (inside Build), seconds of an uninterrupted build, seconds saved, hash equal, leftover ok.
The work directory is emptied before every kill; the input is read in place. Needs only the Python standard library."""
import argparse
import hashlib
import os
import random
import shutil
import signal
import subprocess
import sys
import threading
import time

PHASES = [  # (name, first marker, last marker): the phase of a kill is the interval between two markers of an uninterrupted run
    ("ph0", "start", "packed"),
    ("ph1", "packed", "ph23"),
    ("ph23", "ph23", "sa-built"),
    ("ph4", "sa-built", "lcp-child-built"),
    ("bwt", "lcp-child-built", "bwt-built"),
    ("psi", "bwt-built", "psi-built"),
    ("ph7", "psi-built", "magic-written"),
    ("seal", "magic-written", None),
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


class Run:
    """One process of the Bench: stdout is read by a thread so that the markers are seen as they are printed."""

    def __init__(self, bench, args, env=None):
        self.lines = []
        self.marks = []  # (arrival time, name, seconds since the start of the build)
        self.t_launch = time.monotonic()
        self.proc = subprocess.Popen(["dotnet", bench, "sesa", "resume", "build"] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, env=env)
        self.cv = threading.Condition()
        self.eof = False
        self.err = []
        threading.Thread(target=self._read, daemon=True).start()
        threading.Thread(target=lambda: self.err.extend(self.proc.stderr.readlines()), daemon=True).start()

    def _read(self):
        for line in self.proc.stdout:
            line = line.rstrip("\n")
            with self.cv:
                self.lines.append(line)
                parts = line.split("\t")
                if parts[0] == "M":
                    self.marks.append((time.monotonic(), parts[1], float(parts[2])))
                self.cv.notify_all()
        with self.cv:
            self.eof = True
            self.cv.notify_all()

    def wait_mark(self, name, timeout=3600):
        end = time.monotonic() + timeout
        with self.cv:
            while True:
                for m in self.marks:
                    if m[1] == name:
                        return m
                if self.eof or time.monotonic() >= end:
                    return None
                self.cv.wait(0.1)

    def kill(self):
        killed = False
        if self.proc.poll() is None:
            try:
                self.proc.send_signal(signal.SIGKILL)
                killed = True
            except ProcessLookupError:
                killed = False
        self.proc.wait()
        time.sleep(0.05)
        killed = killed and self.proc.returncode == -signal.SIGKILL
        return killed

    def finish(self):
        self.proc.wait()
        time.sleep(0.05)
        return self.proc.returncode

    def rows(self):
        out = {}
        for line in self.lines:
            p = line.split("\t")
            if p[0] == "R" and len(p) >= 5:
                out[p[2]] = p[3]
        return out

    def last_mark(self):
        return self.marks[-1][1] if self.marks else "(none)"


def clean(work):
    for name in os.listdir(work):
        path = os.path.join(work, name)
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def entries(work):
    return sorted(os.listdir(work))


def phase_of(last):
    """The phase that a kill was in, from the last marker that the process printed before it died."""
    order = ["(none)", "start", "scanned", "created"]
    if last in order:
        return "ph0"
    if last in ("packed", "ph1", "rank-written"):
        return "ph1"
    if last == "ph23":
        return "ph23"
    if last in ("sa-built", "b1", "b2", "b3", "lcp-placed"):
        return "ph4"
    if last == "lcp-child-built":
        return "bwt"
    if last == "bwt-built":
        return "psi"
    if last in ("psi-built", "info-written", "checksummed", "header-written", "header-synced"):
        return "ph7"
    return "seal"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--kind", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--per-phase", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--boundaries", default="yes")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chain", type=float, default=0.25, help="probability that the resume is itself killed once at a random moment before the run that finishes")
    ap.add_argument("--psi", default="true")
    ap.add_argument("--lock-check", default="no")
    ap.add_argument("--fsync", default="true")
    a = ap.parse_args()
    rnd = random.Random(a.seed)
    os.makedirs(a.work, exist_ok=True)
    n = os.path.getsize(a.input)
    target = os.path.join(a.work, "out.sesa")

    def args(resume, threads=None, extra=()):
        return ["--input", a.input, "--work", a.work, "--resume", resume, "--threads", str(threads or a.threads), "--psi", a.psi, "--fsync", a.fsync, "--markers", "true", "--hash", "true", "--keep", "true", "--label", a.kind] + list(extra)

    # ---- the reference: the build without Resume, and the build with Resume that nothing interrupts (its timeline decides where the random kills aim)
    clean(a.work)
    ref = Run(a.bench, args("false"))
    rc = ref.finish()
    assert rc == 0, ("the reference build failed", ref.err[-3:], ref.lines[-3:])
    rrows = ref.rows()
    ref_hash, ref_total = rrows["sha256"], float(rrows["total_s"])
    assert entries(a.work) == ["out.sesa"], entries(a.work)
    clean(a.work)
    tl = Run(a.bench, args("true"))
    assert tl.finish() == 0, ("the timeline build failed", tl.err[-3:])
    trows = tl.rows()
    assert trows["sha256"] == ref_hash, "a build with Resume and no interruption must have the bytes of the build without"
    assert entries(a.work) == ["out.sesa"], entries(a.work)
    marks = {name: secs for _, name, secs in tl.marks}
    order = [name for _, name, _ in tl.marks]
    end = float(trows["total_s"])
    print(f"# {a.kind} n={n} reference sha256 {ref_hash[:16]}.. total {ref_total:.2f} s (Resume off), {end:.2f} s (Resume on); markers {len(order)}", flush=True)
    clean(a.work)

    out = open(a.rows, "a")

    def log(mode, aimed, last, killed, target_state, skipped, resume_total, hash_ok, leftover_ok, kid, extra=""):
        phase = phase_of(last)
        saved = ref_total - resume_total if resume_total is not None else float("nan")
        out.write("\t".join(str(x) for x in [a.kind, n, kid, mode, aimed, last, phase, int(killed), target_state, skipped, "" if resume_total is None else f"{resume_total:.3f}", f"{ref_total:.3f}", f"{saved:.3f}", int(hash_ok), int(leftover_ok), extra]) + "\n")
        out.flush()

    def one(kid, mode, aimed, victim, chain):
        """A kill of `victim` (a Run), then the resume(s), then the checks. Returns True when everything held."""
        killed = victim.kill()
        last = victim.last_mark()
        target_state = "none"
        if os.path.exists(target):
            target_state = "complete" if sha256_file(target) == ref_hash else "damaged"
        ok_hash = target_state != "damaged"
        extra = ""
        if chain and rnd.random() < chain:
            # the resume is itself killed once, at a random moment of its own run
            mid = Run(a.bench, args("true", rnd.choice([1, 2, a.threads])))
            t = rnd.uniform(0.0, 0.8 * end)
            if mid.wait_mark("start", 60):
                time.sleep(t)
            k2 = mid.kill()
            extra = f"chain-killed@{mid.last_mark()}" if k2 else "chain-finished-first"
            if os.path.exists(target) and sha256_file(target) != ref_hash:
                ok_hash = False
                target_state = "damaged"
        fin = Run(a.bench, args("true", rnd.choice([1, 2, a.threads])))
        rc = fin.finish()
        fr = fin.rows()
        h = fr.get("sha256", "")
        hash_equal = rc == 0 and h == ref_hash and ok_hash
        leftover_ok = entries(a.work) == ["out.sesa"]
        skipped = fr.get("skipped", "")
        total = float(fr["total_s"]) if "total_s" in fr else None
        log(mode, aimed, last, killed, target_state, skipped, total, hash_equal, leftover_ok, kid, extra)
        status = "ok" if hash_equal and leftover_ok else "FAILED"
        print(f"{status} {a.kind} {kid} {mode} aimed {aimed} last {last} phase {phase_of(last)} killed {int(killed)} target {target_state} skipped [{skipped}] resume {total if total is None else round(total, 2)} s {extra}", flush=True)
        if not (hash_equal and leftover_ok):
            print("   exit", rc, "hash", h[:16], "entries", entries(a.work), "stderr", fin.err[-2:], flush=True)
        clean(a.work)
        return hash_equal and leftover_ok

    good = True
    kid = 0
    if a.boundaries == "yes":
        for name in order:
            kid += 1
            clean(a.work)
            victim = Run(a.bench, args("true", extra=["--stop-at", name]))
            if victim.wait_mark(name, 3600) is None:
                print("marker never seen:", name, flush=True)
                victim.kill()
                good = False
                continue
            time.sleep(0.05)
            good &= one(kid, "boundary", name, victim, 0.0)

    for phase, first, last in PHASES:
        lo = marks.get(first, 0.0)
        hi = marks[last] if last in marks else end
        for _ in range(a.per_phase):
            kid += 1
            clean(a.work)
            t = rnd.uniform(lo, hi)
            victim = Run(a.bench, args("true"))
            m = victim.wait_mark("start", 120)
            if m is None:
                good = False
                victim.kill()
                continue
            # the markers are timed from the start of the build; the signal is sent t seconds after the "start" marker arrived
            deadline = m[0] + t
            while time.monotonic() < deadline and victim.proc.poll() is None:
                time.sleep(0.002)
            good &= one(kid, "random", f"{phase}@{t:.3f}", victim, a.chain)

    if a.lock_check == "yes":
        clean(a.work)
        holder = Run(a.bench, args("true", extra=["--stop-at", "sa-built"]))
        assert holder.wait_mark("sa-built", 3600) is not None
        second = Run(a.bench, args("true"))
        rc = second.finish()
        refused = rc == 3 and "lock" in " ".join(second.lines)
        before = entries(a.work)
        holder.kill()
        fin = Run(a.bench, args("true"))
        ok = fin.finish() == 0 and fin.rows().get("sha256") == ref_hash and entries(a.work) == ["out.sesa"]
        print(f"{'ok' if refused and ok else 'FAILED'} {a.kind} lock-check: the second build exited {rc} ({'refused with the lock message' if refused else 'NOT refused'}), the directory held {before} meanwhile, the build after the kill of the holder {'finished identical' if ok else 'FAILED'}", flush=True)
        good &= refused and ok
        clean(a.work)

    out.close()
    print("ALL OK" if good else "SOME FAILED", flush=True)
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
