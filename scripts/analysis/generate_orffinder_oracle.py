#!/usr/bin/env python3
"""GENOMIC-ORF-001 (A1-18 / A1-19) oracle: converts the returned output of the NCBI ORFfinder 0.4.3 hand-off package into the
compact embedded test oracle `orffinder_cli_oracle.tsv`.

The oracle is the REAL NCBI standalone program (https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz, version
0.4.3, sha256 f30a0806...09f2), run by the repository owner in Docker (linux/amd64, Debian 12) over the 629 deterministic inputs and
the 1750-run matrix of `docs/Validation/review-2026-09/handoff/A1-18-orffinder/` (package README / INTEGRATION.md). The returned
archive is `.../returned/A1-18-orffinder-out.zip`; this script reads it (or an extracted `out/` directory) as UNTRUSTED DATA: text only,
nothing in it is executed.

  python3 -I scripts/analysis/generate_orffinder_oracle.py [--out FILE] [--zip ZIP | --dir OUTDIR] [--stats]

Formats actually observed (all asserted while converting):
  -outfmt 0 : `>lcl|ORF<n>_<seqid>:<from>:<to> unnamed protein product[, partial]` + protein FASTA; <from>:<to> are 0-BASED
              INCLUSIVE forward-strand coordinates of the ORF including its stop codon; from > to marks the minus strand; ORFs of
              one sequence are numbered 1..n in output order (plus strand first, then minus; each by frame of the strand-coordinate
              start, then start); `partial` = 5'- and/or 3'-partial.
  -outfmt 1 : CDS FASTA `>lcl|<seqid>:<from+1>-<to+1>` (1-based; minus strand `c<hi+1>-<lo+1>`) `ORF<n>_...` + the CDS nucleotides.
  -outfmt 2 : text ASN.1; `fuzz-from lim lt` / `fuzz-to lim gt` on the location say which END (low / high coordinate) is open.
  -outfmt 3 : feature table `>Feature lcl|ORF<n>_<seqid>:<from>:<to>`, `from<TAB>to<TAB>CDS` (0-based, no `<` / `>` markers).
  stderr    : per sequence `Error (<seqid>) : <message>` (`No ORFs found for the specified parameters`, `The sequence length is N nt.
              It's too short for the minimal ORF length.`, `CFastaReader: ...invalid residue(s)...`, region errors, unknown genetic code).

Output TSV (tab separated, three line types, `#` lines are comments):
  S  seqid  sequence                                  every distinct input sequence (the id is `<set>:<fasta id>`)
  C  cfgid  g  s  ml  n  strand  c  b  e              one ORFfinder option set; c is `0`, `1` or `01` (observed byte-identical)
  R  cfgid  seqid  orfs  proteinhash                  one (option set, sequence) result;
       orfs = `;`-joined `strand,lo,hi,flag` (0-based inclusive forward coordinates, stop included, ORFfinder output order) or `-`;
              flag: `0` complete, `P` partial (outfmt 0, 5'/3' not distinguished), `5` / `3` / `B` 5'-only / 3'-only / both (outfmt 2),
              `?` unknown (outfmt 3 only run);
       proteinhash = first 8 hex digits of SHA-256 of the `|`-joined proteins exactly as printed by ORFfinder (`-` for outfmt 3 only);
       or `ERR:<kind>` when ORFfinder rejected the option set for this sequence (region start/end errors).
  E  cfgid  message                                   option sets ORFfinder rejects for every sequence (unknown genetic code ...)
Sequences ORFfinder could not read (`CFastaReader` reader heuristics: very short / C-rich / protein-like / U / gap) have NO R line.
"""
import collections
import hashlib
import io
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PKG = os.path.join(REPO, "docs", "Validation", "review-2026-09", "handoff", "A1-18-orffinder")
DEFAULT_ZIP = os.path.join(PKG, "returned", "A1-18-orffinder-out.zip")
DEFAULT_OUT = os.path.join(REPO, "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData", "Analysis", "orffinder_cli_oracle.tsv")
BOOL = {"true": "1", "t": "1", "yes": "1", "TRUE": "1", "1": "1", "false": "0", "f": "0", "no": "0"}
BAD_TABLES_MESSAGE = "Unable to find genetic code number"
HDR = re.compile(r"^>lcl\|ORF(\d+)_(.*):(-?\d+):(-?\d+) (unnamed protein product(?:, partial)?)$")
ASN = re.compile(r'product whole local str "ORF(\d+)_(.*?):(-?\d+):(-?\d+)",\s*location int \{\s*from (-?\d+),\s*to (-?\d+),\s*strand (plus|minus),'
                 r'\s*id local str "([^"]*)"((?:,\s*fuzz-(?:from|to) lim (?:lt|gt))*)', re.S)
FEAT = re.compile(r"^>Feature lcl\|ORF(\d+)_(.*):(-?\d+):(-?\d+)$")


class Out:
    """Read access to the returned out/ tree (zip or directory)."""

    def __init__(self, zip_path=None, directory=None):
        self.zf = zipfile.ZipFile(zip_path) if zip_path else None
        self.root = directory
        self.names = {}
        if self.zf:
            for n in self.zf.namelist():
                i = n.find("out/")
                self.names[n[i + 4:] if i >= 0 else n] = n

    def read(self, rel):
        if self.zf:
            n = self.names.get(rel)
            return self.zf.read(n).decode("utf-8", "replace") if n else None
        p = os.path.join(self.root, rel)
        return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else None


def read_fasta(text):
    """{id: sequence}; ids that occur more than once map to None (ambiguous: not used)."""
    d, cur = collections.OrderedDict(), None
    for line in text.splitlines():
        if line.startswith(">"):
            h = line[1:].strip()
            cur = h.split()[0] if h else ""
            d[cur] = None if cur in d else []
        elif cur is not None and d[cur] is not None:
            d[cur].append(line.strip())
    return collections.OrderedDict((k, "".join(v) if v is not None else None) for k, v in d.items())


def parse_fmt0(text):
    """{seqid: [(n, a, b, partial, protein)]}; a, b as printed (0-based, a > b = minus strand)."""
    out, cur = collections.OrderedDict(), None
    for line in text.splitlines():
        if line.startswith(">"):
            m = HDR.match(line)
            assert m, "unrecognised outfmt 0 header: " + line[:120]
            n, sid, a, b, rest = m.groups()
            cur = [int(n), int(a), int(b), rest.endswith("partial"), []]
            out.setdefault(sid, []).append(cur)
        elif cur is not None:
            cur[4].append(line.strip())
    return {k: [(n, a, b, p, "".join(pr)) for n, a, b, p, pr in v] for k, v in out.items()}


def parse_fmt3(text):
    out, cur = collections.OrderedDict(), None
    lines = text.splitlines()
    for i, line in enumerate(lines):
        m = FEAT.match(line)
        if m:
            n, sid, a, b = m.groups()
            fr, to, kind = lines[i + 1].split("\t")
            assert kind == "CDS" and int(fr) == int(a) and int(to) == int(b), (line, lines[i + 1])
            out.setdefault(sid, []).append((int(n), int(a), int(b)))
    return out


def parse_fmt2(text):
    out = collections.OrderedDict()
    for m in ASN.finditer(text):
        n, sid, a, b, fr, to, strand, sid2, fz = m.groups()
        assert sid == sid2
        lo_open, hi_open = "fuzz-from lim lt" in fz, "fuzz-to lim gt" in fz
        assert "fuzz-from lim gt" not in fz and "fuzz-to lim lt" not in fz
        out.setdefault(sid, []).append((int(n), int(fr), int(to), "+" if strand == "plus" else "-", lo_open, hi_open))
    return out


def parse_fmt1(text):
    """{seqid: [(n, a, b, cds)]} a, b 0-based from the `ORF<n>_<id>:<a>:<b>` tail of the header."""
    out, cur = collections.OrderedDict(), None
    for line in text.splitlines():
        if line.startswith(">"):
            m = re.match(r"^>lcl\|(.*?):(c?)(\d+)-(\d+) ORF(\d+)_(.*):(-?\d+):(-?\d+)$", line)
            assert m, "unrecognised outfmt 1 header: " + line[:120]
            sid, c, x, y, n, sid2, a, b = m.groups()
            assert sid == sid2
            if c:
                assert int(x) - 1 == int(a) and int(y) - 1 == int(b), line
            else:
                assert int(x) - 1 == int(a) and int(y) - 1 == int(b), line
            cur = [int(n), int(a), int(b), []]
            out.setdefault(sid, []).append(cur)
        elif cur is not None:
            cur[3].append(line.strip())
    return {k: [(n, a, b, "".join(c)) for n, a, b, c in v] for k, v in out.items()}


def stderr_messages(text):
    """{seqid: message} from `Error (<seqid>) : <message>` lines ('' = the sequence-less errors)."""
    res = {}
    for line in (text or "").splitlines():
        m = re.match(r"^Error \((.*?)\) : (.*)$", line)
        if m:
            res.setdefault(m.group(1), m.group(2))
    return res


def parse_args(a):
    t = a.split()
    d = dict(g="1", s="1", ml="75", n="false", strand="both", c="false", b="1", e="0", outfmt="0")
    i = 0
    while i < len(t):
        k = t[i].lstrip("-")
        if k in d and i + 1 < len(t):
            d[k] = t[i + 1]
            i += 2
        else:
            i += 1
    return d


def load_inputs(pkg):
    """{(set, fasta id): sequence}; the set is the input file (designed, random, circular, region, edge/<file>, circular_edge/<file>)."""
    seq = {}
    for name in ("designed", "random", "circular", "region"):
        for k, v in read_fasta(open(os.path.join(pkg, "inputs", name + ".fa"), encoding="utf-8").read()).items():
            seq[(name, k)] = v
    for sub in ("edge", "circular_edge"):
        d = os.path.join(pkg, "inputs", sub)
        for fn in sorted(os.listdir(d)):
            for k, v in read_fasta(open(os.path.join(d, fn), encoding="utf-8", newline="").read().replace("\r", "")).items():
                if k and v is not None and v != "":
                    seq[(sub + "/" + fn[:-3], k)] = v
    return seq


def input_of(infile):
    """(set, single record id or None) of a job input path."""
    if infile.startswith("work/single/"):
        return "designed", os.path.basename(infile)[:-3]
    if infile.startswith("inputs/edge/"):
        return "edge/" + os.path.basename(infile)[:-3], None
    if infile.startswith("inputs/circular_edge/"):
        return "circular_edge/" + os.path.basename(infile)[:-3], None
    return os.path.basename(infile)[:-3], None


def phash(proteins):
    return hashlib.sha256("|".join(proteins).encode()).hexdigest()[:8]


def flag_of(p5, p3):
    return "B" if p5 and p3 else "5" if p5 else "3" if p3 else "0"


def merge_results(a, b, what):
    if isinstance(a, str) or isinstance(b, str):
        assert a == b, ("conflict", what, a, b)
        return a
    la = [x.split(",") for x in a[0].split(";")] if a[0] != "-" else []
    lb = [x.split(",") for x in b[0].split(";")] if b[0] != "-" else []
    assert [x[:3] for x in la] == [x[:3] for x in lb], ("conflict", what, a, b)
    merged = []
    for x, y in zip(la, lb):
        fx, fy = x[3], y[3]
        if fx == fy or fy == "?":
            f = fx
        elif fx == "?":
            f = fy
        elif fx == "P" and fy in "53B":
            f = fy
        elif fy == "P" and fx in "53B":
            f = fx
        else:
            raise AssertionError(("flag conflict", what, a, b))
        merged.append(",".join(x[:3] + [f]))
    ha = a[1] if a[1] != "-" else b[1]
    assert a[1] == b[1] or "-" in (a[1], b[1]), ("protein hash conflict", what, a, b)
    return (";".join(merged) if merged else "-", ha)


def convert(out, pkg):
    """-> (cfgs, rows, errs, seqs, nruns, checks)
    cfgs: set of (g, s, ml, n, strand, b, e) ; rows: {(cfg, c, seqkey): result}; errs: {(g,): message}."""
    jobs = [j + [""] * (6 - len(j)) for j in (l.split("\t") for l in out.read("jobs.tsv").splitlines())]
    byid = {j[0]: j for j in jobs}
    seqs = load_inputs(pkg)
    rows, errs, cfgs = {}, {}, set()
    checks = collections.Counter()
    nruns = 0
    for jid, group, infile, tmo, useout, args in jobs:
        if group == "misc" or useout != "1":
            continue
        is3only = group == "circedge"            # these genomes were run with -outfmt 3 only
        if not (jid.endswith("fmt0") or (is3only and jid.endswith("fmt3"))):
            continue
        a = parse_args(args)
        if a["n"] not in BOOL or a["c"] not in BOOL or a["strand"] not in ("both", "plus", "minus"):
            continue                              # invalid option spellings are behaviour rows, not oracle rows
        nruns += 1
        g = int(a["g"])
        cfg = (g, int(a["s"]), int(a["ml"]), BOOL[a["n"]], a["strand"], int(a["b"]), int(a["e"]))
        c = BOOL[a["c"]]
        sname, only = input_of(infile)
        err = stderr_messages(out.read("runs/%s/stderr.txt" % jid))
        bad = [m for m in err.values() if m.startswith(BAD_TABLES_MESSAGE)]
        if bad:
            errs[(g,)] = bad[0]
            continue
        cfgs.add(cfg)
        text = out.read("runs/%s/orfs.txt" % jid) or ""
        fmt0 = jid.endswith("fmt0")
        got = parse_fmt0(text) if fmt0 else parse_fmt3(text)
        stem = jid[:-1]
        t3 = out.read("runs/%s3/orfs.txt" % stem) if fmt0 and (stem + "3") in byid else None
        t2 = out.read("runs/%s2/orfs.txt" % stem) if fmt0 and (stem + "2") in byid else None
        t1 = out.read("runs/%s1/orfs.txt" % stem) if fmt0 and (stem + "1") in byid else None
        f3 = parse_fmt3(t3) if t3 else None
        f2 = parse_fmt2(t2) if t2 else None
        f1 = parse_fmt1(t1) if t1 else None
        keys = [k for (n, k) in seqs if n == sname and (only is None or k == only)]
        for sid in keys:
            m = err.get(sid)
            skey = "%s:%s" % (sname, sid)
            if m is not None and m.startswith("CFastaReader"):
                checks["reader-rejected"] += 1
                continue
            if m is not None and ("Interval Start" in m or "out of range" in m):
                res = "ERR:start>end" if "Interval Start" in m else "ERR:start-out-of-range"
            elif sid in got:
                lst = got[sid]
                assert [x[0] for x in lst] == list(range(1, len(lst) + 1)), (jid, sid, "ORF numbering")
                if fmt0:
                    items = []
                    spans = []
                    for n, x, y, p, prot in lst:
                        strand, lo, hi = ("+", x, y) if x <= y else ("-", y, x)
                        spans.append((strand, lo, hi, p, prot))
                    if f3 is not None:
                        assert [(n, a_, b_) for n, a_, b_, *_ in lst] == f3.get(sid, []), (jid, sid, "outfmt 3 disagrees")
                        checks["fmt3-agrees"] += 1
                    if f1 is not None:
                        assert [(n, a_, b_) for n, a_, b_, _ in f1.get(sid, [])] == [(n, a_, b_) for n, a_, b_, *_ in lst], (jid, sid, "outfmt 1 disagrees")
                        checks["fmt1-agrees"] += 1
                    flags = ["P" if p else "0" for _, _, _, p, _ in spans]
                    if f2 is not None:
                        l2 = f2.get(sid, [])
                        assert [(n, s_, lo, hi) for n, lo, hi, s_, _, _ in l2] == [(i + 1, s[0], s[1], s[2]) for i, s in enumerate(spans)], (jid, sid, "outfmt 2 disagrees")
                        flags = [flag_of(lo_open if s_ == "+" else hi_open, hi_open if s_ == "+" else lo_open) for _, lo, hi, s_, lo_open, hi_open in l2]
                        for (_, _, _, p, _), fl in zip(spans, flags):
                            assert p == (fl != "0"), (jid, sid, "partial flag disagrees with outfmt 2")
                        checks["fmt2-agrees"] += 1
                    res = (";".join("%s,%d,%d,%s" % (s_[0], s_[1], s_[2], fl) for s_, fl in zip(spans, flags)), phash([s_[4] for s_ in spans]))
                else:
                    res = (";".join("%s,%d,%d,?" % (("+", x, y) if x <= y else ("-", y, x)) for n, x, y in lst), "-")
            else:
                assert m is not None and (m.startswith("No ORFs found") or "too short" in m or m.startswith("No output produced")), (jid, sid, m)
                res = ("-", phash([]) if fmt0 else "-")
            key = (cfg, c, skey)
            old = rows.get(key)
            if old is None:
                rows[key] = res
            else:
                # the same (option set, sequence) measured by several jobs (factorial / minlength / circular / single ...): must agree;
                # an exact 5'/3' flag (outfmt 2) refines the combined `P` flag of outfmt 0
                rows[key] = merge_results(old, res, (key, jid))
                checks["duplicate-agrees"] += 1
    return cfgs, rows, errs, seqs, nruns, checks


def write(cfgs, rows, errs, seqs, fh):
    """Collapse c: (cfg, c, seq) results that are identical for c = 0 and 1 become one `01` config."""
    used = sorted({k[2] for k in rows})
    seqmap = {"%s:%s" % (n, k): v for (n, k), v in seqs.items()}
    cfg_rows = collections.defaultdict(dict)       # (cfg, c) -> {seq: res}
    for (cfg, c, sk), res in rows.items():
        cfg_rows[(cfg, c)][sk] = res
    merged = collections.defaultdict(dict)         # (cfg, c-label) -> {seq: res}
    for (cfg, c), d in cfg_rows.items():
        if c == "1":
            continue
        for sk, res in d.items():
            twin = cfg_rows.get((cfg, "1"), {}).get(sk)
            if twin is not None:
                res = merge_results(res, twin, ("-c true differs from -c false", cfg, sk))
                d[sk] = res
            merged[(cfg, "01" if twin is not None else "0")][sk] = res
    for (cfg, c), d in cfg_rows.items():
        if c == "1":
            for sk, res in d.items():
                if sk not in merged.get((cfg, "01"), {}):
                    merged[(cfg, "1")][sk] = res
    fh.write("# GENOMIC-ORF-001 oracle: NCBI ORFfinder 0.4.3 standalone (see scripts/analysis/generate_orffinder_oracle.py); generated, do not edit\n")
    fh.write("# S seqid sequence | C cfgid g s ml n strand c b e | R cfgid seqid orfs proteinhash | E g message\n")
    for sk in used:
        fh.write("S\t%s\t%s\n" % (sk, seqmap[sk]))
    ids = {}
    for i, key in enumerate(sorted(merged), 1):
        ids[key] = "c%d" % i
        cfg, c = key
        fh.write("C\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n" % ((ids[key],) + cfg[:5] + (c,) + cfg[5:]))
    for key in sorted(merged):
        for sk in sorted(merged[key]):
            res = merged[key][sk]
            fh.write("R\t%s\t%s\t%s\n" % (ids[key], sk, res if isinstance(res, str) else "%s\t%s" % res))
    for (g,), msg in sorted(errs.items()):
        fh.write("E\t%d\t%s\n" % (g, msg))
    return len(used), len(merged), sum(len(d) for d in merged.values())


def main():
    argv = sys.argv[1:]
    outp = DEFAULT_OUT
    zp, dp = DEFAULT_ZIP, None
    if "--out" in argv:
        outp = argv[argv.index("--out") + 1]
    if "--zip" in argv:
        zp = argv[argv.index("--zip") + 1]
    if "--dir" in argv:
        dp, zp = argv[argv.index("--dir") + 1], None
    out = Out(zip_path=zp, directory=dp)
    cfgs, rows, errs, seqs, nruns, checks = convert(out, PKG)
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        ns, nc, nr = write(cfgs, rows, errs, seqs, fh)
    print("runs=%d sequences=%d configs=%d rows=%d genetic-code errors=%s checks=%s -> %s (%d bytes)" % (
        nruns, ns, nc, nr, sorted(errs), dict(checks), outp, os.path.getsize(outp)), file=sys.stderr)


if __name__ == "__main__":
    main()
