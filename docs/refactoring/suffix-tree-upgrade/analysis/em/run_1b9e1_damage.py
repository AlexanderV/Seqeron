#!/usr/bin/env python3
"""Work package 1b.9e.1: a resumed build must not trust any byte of the temporary container that no kept stage vouches for.

The driver builds the input once without Resume (the reference), then, for every marker at which `sesa resume build` can be stopped (every SesaHooks boundary and the points inside PH1, PH3 and PH4), starts a
build with Resume, lets it block at the marker (`--stop-at`), kills it with SIGKILL, takes a snapshot of the work directory, and for every variant copies the snapshot, damages the temporary file *outside the ranges that the
state vouches for*, resumes in a NEW process and compares the SHA-256 of the finished container (and `cmp -l`, to name the bytes that differ) with the reference.

variants (what is written into the part of the file that no kept stage has a CRC for):
  fill_a5      every byte outside the kept ranges, from the header page to the end of the file, becomes 0xA5
  fill_rnd     the same with pseudo-random bytes (a second pattern)
  tail_a5      only the part after the end of the last kept range (the orchestrator's second case)
  gap_first / gap_mid / gap_last   one byte of every gap of the layout (the padding between the sections, the room of BUILD_INFO behind its text) that is not inside a kept range,
               in one resume (`--gap-mode all`) or one resume per gap and byte (`--gap-mode each`)
  header_first / header_mid / header_last   one byte of the header page
  child_slot   the last 4-byte word of the CHILD range that the reference has as zero (a slot that the child writer never writes) becomes 0xA5A5A5A5, wherever the marker is
  child_last   one bit of the last byte of the CHILD range is flipped (the orchestrator's first case, at the marker lcp-child-built the range is recorded and the stage must be redone)

modes:
  repro   <bench> <input> <work> <rows>            the two cases of the orchestrator on one input (flip the last byte of the recorded CHILD range; 0xA5 from the end of the last kept range to the end of the file at the
                                                   markers packed and sa-built, and at the markers after the LCP stage), printing the number of bytes that differ from the reference
  matrix  <bench> <input> <work> <rows> [--gap-mode all|each] [--markers a,b,..] [--variants a,b,..] [--threads 4]
Standard library only; every row is one resume: label, marker, variant, exit code, stages skipped, SHA-256 equal, bytes differing, files left, seconds."""
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time

HEX = re.compile(r"\b[0-9a-f]{64}\b")
ALL_MARKERS = ["created", "packed", "ph1", "rank-written", "ph23", "sa-built", "lcp-placed", "b1", "b2", "b3", "lcp-child-built", "bwt-built", "psi-built",
               "info-written", "checksummed", "header-written", "header-synced", "magic-written"]
ALL_VARIANTS = ["fill_a5", "fill_rnd", "tail_a5", "gap_first", "gap_mid", "gap_last", "header_first", "header_mid", "header_last", "child_slot", "child_last"]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def fresh(work):
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)


def cmd(bench, inp, work, threads, resume, stop_at=None, hash_=True):
    c = ["dotnet", bench, "sesa", "resume", "build", "--input", inp, "--work", work, "--threads", str(threads), "--resume", "true" if resume else "false",
         "--keep", "true", "--hash", "true" if hash_ else "false", "--psi", "true"]
    if stop_at:
        c += ["--stop-at", stop_at]
    return c


def rows_of(stdout):
    d = {}
    for line in stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 5 and p[0] == "R" and p[1] == "resume":
            d[p[2]] = p[3]
    return d


def finish(bench, inp, work, threads, resume=True):
    t = time.monotonic()
    p = subprocess.run(cmd(bench, inp, work, threads, resume), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, stdin=subprocess.DEVNULL, timeout=3600)
    m = HEX.findall(p.stdout)
    return p.returncode, (m[-1] if m else None), time.monotonic() - t, rows_of(p.stdout)


def stop_and_kill(bench, inp, work, threads, marker, timeout=900):
    """Starts a resumable build that blocks at `marker`, waits for the marker line, kills the process. False when the build ended without reaching it."""
    p = subprocess.Popen(cmd(bench, inp, work, threads, True, stop_at=marker, hash_=False), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL, text=True)
    reached = False
    t0 = time.monotonic()
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] == "M" and f[1] == marker:
            reached = True
            break
        if time.monotonic() - t0 > timeout:
            break
    if p.poll() is None:
        p.kill()
    p.wait()
    return reached


def state_of(work):
    path = os.path.join(work, "out.sesa.tmp.d", "state.json")
    try:
        return json.loads(open(path, "rb").read())
    except (OSError, ValueError):
        return None


def kept_ranges(state):
    """The container ranges that the stages of the state vouch for (offset, length)."""
    out = []
    for s in (state or {}).get("stages", []):
        for f in s.get("files", []):
            if f.get("container"):
                out.append((f["offset"], f["length"]))
    return sorted(out)


def complement(ranges, length):
    gaps, at = [], 0
    for off, ln in sorted(ranges):
        if off > at:
            gaps.append((at, min(off, length) - at))
        at = max(at, off + ln)
        if at >= length:
            break
    if at < length:
        gaps.append((at, length - at))
    return [(o, l) for o, l in gaps if l > 0]


_NONZERO = bytes([1]) + bytes(range(1, 256))  # translate table: 0 -> 1, every other byte stays


def fill(path, spans, pattern, seed=1):
    """Writes `pattern` (a byte value, or "rnd" for pseudo-random non-zero bytes) over the spans of the file."""
    rnd = random.Random(seed)
    with open(path, "r+b") as f:
        for off, ln in spans:
            f.seek(off)
            left = ln
            while left > 0:
                k = min(left, 1 << 22)
                f.write(rnd.randbytes(k).translate(_NONZERO) if pattern == "rnd" else bytes([pattern]) * k)
                left -= k


def poke(path, offset, value):
    with open(path, "r+b") as f:
        f.seek(offset)
        f.write(bytes([value]))


def diff_bytes(a, b):
    r = subprocess.run(["cmp", "-l", a, b], capture_output=True, text=True)
    lines = r.stdout.splitlines()
    if r.returncode == 0:
        return 0, []
    if not lines:  # cmp: EOF on a shorter file
        return -1, [r.stderr.strip()]
    return len(lines), [int(x.split()[0]) - 1 for x in lines[:5]]


def layout_of(bench, inp, work, threads):
    """All ranges of the finished layout and the length of the finished container, from the state of a build that is killed after PSI (every section is recorded then)."""
    fresh(work)
    if not stop_and_kill(bench, inp, work, threads, "info-written"):
        raise SystemExit("the build ended before info-written")
    st = state_of(work)
    ranges = {}
    for s in st["stages"]:
        for f in s["files"]:
            if f.get("container"):
                ranges[f["role"]] = (f["offset"], f["length"])
    return ranges, os.path.getsize(os.path.join(work, "out.sesa.tmp")), st


def main():
    a = [x for x in sys.argv[1:]]
    if len(a) < 5 or a[0] not in ("repro", "matrix"):
        print(__doc__)
        raise SystemExit(2)
    mode, bench, inp, work, rows = a[:5]
    opts = {a[i]: a[i + 1] for i in range(5, len(a) - 1, 2) if a[i].startswith("--")}
    threads = int(opts.get("--threads", "4"))
    label = os.path.basename(inp)
    refdir = work + "_ref"
    snap = work + "_snap"
    var = work + "_var"
    fresh(refdir)
    rc, ref_hash, _, _ = finish(bench, inp, refdir, threads, resume=False)
    ref_file = os.path.join(refdir, "out.sesa")
    assert rc == 0 and ref_hash, "the reference build failed"
    print("reference", ref_hash, os.path.getsize(ref_file))
    out = open(rows, "a")

    def run_variant(marker, variant, damage):
        shutil.rmtree(var, ignore_errors=True)
        shutil.copytree(snap, var, symlinks=True)
        damage(os.path.join(var, "out.sesa.tmp"))
        rc, h, dt, r = finish(bench, inp, var, threads)
        files = sorted(os.listdir(var))
        eq = h == ref_hash
        nd, first = (0, []) if eq else diff_bytes(os.path.join(var, "out.sesa"), ref_file)
        left = files == ["out.sesa"]
        ok = rc == 0 and eq and left
        out.write("\t".join(map(str, [label, marker, variant, rc, r.get("skipped", ""), "eq" if eq else "DIFF", nd if not eq else 0, "clean" if left else "LEFT:" + ",".join(files),
                                     round(dt, 2), r.get("clear_bytes", ""), r.get("clear_scanned_bytes", ""), r.get("clear_s", "")])) + "\n")
        out.flush()
        return ok, nd, r, first

    layout, final_len, _ = layout_of(bench, inp, work, threads)
    print("layout", {k: v for k, v in layout.items()}, "final length", final_len)
    child_off, child_len = layout["child"]
    refbytes = open(ref_file, "rb").read(child_off + child_len)
    zero_word = None
    for at in range(child_off + child_len - 4, child_off, -4):
        if refbytes[at:at + 4] == b"\0\0\0\0":
            zero_word = at
            break
    ranges_all = sorted(layout.values())
    gaps_all = complement(ranges_all, final_len)  # includes the header page, the padding, the room of BUILD_INFO behind its text, and (at the end) nothing
    total = good = 0
    if mode == "repro":
        cases = [("lcp-child-built", "child_last"), ("packed", "tail_a5"), ("sa-built", "tail_a5"), ("lcp-child-built", "tail_a5"), ("psi-built", "tail_a5"),
                 ("lcp-child-built", "fill_a5"), ("psi-built", "fill_a5"), ("info-written", "fill_a5")]
    else:
        markers = opts["--markers"].split(",") if "--markers" in opts else ALL_MARKERS
        variants = opts["--variants"].split(",") if "--variants" in opts else ALL_VARIANTS
        cases = [(m, v) for m in markers for v in variants]
    gap_mode = opts.get("--gap-mode", "all")
    by_marker = {}
    for m, v in cases:
        by_marker.setdefault(m, []).append(v)
    for marker, variants in by_marker.items():
        fresh(work)
        if not stop_and_kill(bench, inp, work, threads, marker):
            print("marker not reached", marker)
            continue
        shutil.rmtree(snap, ignore_errors=True)
        shutil.copytree(work, snap, symlinks=True)
        st = state_of(snap)
        tmp_len = os.path.getsize(os.path.join(snap, "out.sesa.tmp"))
        kept = kept_ranges(st)
        comp = complement(kept, tmp_len)
        end_kept = max((o + l for o, l in kept), default=0)
        for v in variants:
            if v == "fill_a5":
                damages = [("", lambda p: fill(p, comp, 0xA5))]
            elif v == "fill_rnd":
                damages = [("", lambda p: fill(p, comp, "rnd", 7))]
            elif v == "tail_a5":
                damages = [("", lambda p: fill(p, [(end_kept, tmp_len - end_kept)] if tmp_len > end_kept else [], 0xA5))]
            elif v.startswith("gap_"):
                where = v[4:]
                spans = [(o, l) for o, l in gaps_all if o > 0 and o < tmp_len and not any(ko <= o and o + 1 <= ko + kl for ko, kl in kept)]
                if not spans:
                    continue
                pick = (lambda o, l: o) if where == "first" else (lambda o, l: o + l // 2) if where == "mid" else (lambda o, l: o + l - 1)
                at = [min(pick(o, l), tmp_len - 1) for o, l in spans]
                if gap_mode == "each":
                    damages = [(f"@{x}", (lambda x: lambda p: poke(p, x, 0xA5))(x)) for x in at]
                else:
                    damages = [("", lambda p: [poke(p, x, 0xA5) for x in at])]
            elif v.startswith("header_"):
                x = {"first": 8, "mid": 2048, "last": 4095}[v[7:]]
                damages = [("", (lambda x: lambda p: poke(p, x, 0xA5))(x))]
            elif v == "child_slot":
                if zero_word is None or zero_word + 4 > tmp_len:
                    continue
                damages = [("", lambda p: fill(p, [(zero_word, 4)], 0xA5))]
            elif v == "child_last":
                def flip(p):
                    with open(p, "r+b") as f:
                        f.seek(child_off + child_len - 1)
                        b = f.read(1)
                        f.seek(child_off + child_len - 1)
                        f.write(bytes([b[0] ^ 1]))
                damages = [("", flip)] if child_off + child_len <= tmp_len else []
            else:
                raise SystemExit("unknown variant " + v)
            for suffix, d in damages:
                ok, nd, r, first = run_variant(marker, v + suffix, d)
                total += 1
                good += 1 if ok else 0
                print(f"{label} {marker:16s} {v + suffix:22s} {'ok' if ok else 'BAD'} skipped={r.get('skipped', '')} differing_bytes={nd} first={first}", flush=True)
    print(label, "cases", total, "good", good, "bad", total - good)
    shutil.rmtree(var, ignore_errors=True)
    shutil.rmtree(snap, ignore_errors=True)
    shutil.rmtree(refdir, ignore_errors=True)


if __name__ == "__main__":
    main()
