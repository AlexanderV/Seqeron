#!/usr/bin/env python3
"""DRAFT converter of the returned ./out of this package into oracle rows (finalised by the integrating session).

STATUS: the INSDC feature-table parser follows the public INSDC feature-table specification (`>Feature <id>`, then
`start<TAB>stop<TAB>CDS`; '<' / '>' mark partial ends; start > stop = complement strand). The ORF-FASTA header parser is a
tolerant GUESS: the exact header layout of ORFfinder 0.4.3 `-outfmt 0/1` was NOT observable in the authoring sandbox (the NCBI
binary could not be run there), so `--probe` prints the first headers of a real run and the integrating session adjusts
HEADER_PATTERNS once, then asserts that fmt 0, fmt 1 and fmt 3 agree on every sequence (the converter refuses to emit a row
otherwise).

  python3 parse_orffinder_out.py --selftest
  python3 parse_orffinder_out.py --probe  <out dir or out zip>
  python3 parse_orffinder_out.py --run    <out dir or out zip> > orffinder_cli_oracle.tsv     (after HEADER_PATTERNS is confirmed)

Output TSV (one row per (job, sequence)), tab separated:
  job_id  args  seq_id  sequence  orfs
where orfs = ';'-joined `strand,start,end,partial5,partial3` with 1-based inclusive forward-strand input coordinates exactly as
ORFfinder prints them (start <= end for both strands), partial5/partial3 in {0,1} in the ORF's own 5'->3' sense; '-' = none.
"""
import io
import os
import re
import sys
import zipfile

# ----- INSDC feature table ------------------------------------------------------------------------------------------
LOC = re.compile(r"^([<>]?)(\d+)\t([<>]?)(\d+)\t(\w+)")


def parse_feature_table(text):
    """{seq id: [(strand, start, end, partial5, partial3)]} with 1-based inclusive forward coordinates."""
    out, cur = {}, None
    for line in text.splitlines():
        if line.startswith(">Feature"):
            cur = line.split(None, 1)[1].strip() if " " in line else ""
            out.setdefault(cur, [])
            continue
        m = LOC.match(line)
        if m and cur is not None and m.group(5) == "CDS":
            p1, a, p2, b = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
            if a <= b:
                out[cur].append(("+", a, b, int(p1 == "<"), int(p2 == ">")))
            else:  # complement: first coordinate is the 5' end of the ORF = the larger forward coordinate
                out[cur].append(("-", b, a, int(p2 == ">"), int(p1 == "<")))
    return out


# ----- ORF FASTA headers (GUESS, see module docstring) ----------------------------------------------------------------
HEADER_PATTERNS = [
    # >lcl|ORF1_<seqid>:<from>:<to>  (forward from < to, reverse from > to) -- unverified guess
    re.compile(r"^>lcl\|ORF(?P<n>\d+)_(?P<id>[^:\s]+)[:\s](?P<a>\d+)[:\-.]+(?P<b>\d+)"),
]


def parse_orf_fasta(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith(">"):
            continue
        for pat in HEADER_PATTERNS:
            m = pat.match(line)
            if m:
                a, b = int(m.group("a")), int(m.group("b"))
                out.setdefault(m.group("id"), []).append(("+", a, b, 0, 0) if a <= b else ("-", b, a, 0, 0))
                break
        else:
            raise ValueError("unrecognised ORF FASTA header (adjust HEADER_PATTERNS): " + line[:120])
    return out


def read_fasta(text):
    recs, cur = [], None
    for line in text.splitlines():
        if line.startswith(">"):
            cur = [line[1:].strip().split()[0] if line[1:].strip() else "", []]
            recs.append(cur)
        elif cur is not None:
            cur[1].append(line.strip())
    return {r[0]: "".join(r[1]) for r in recs}


# ----- run directory access -----------------------------------------------------------------------------------------------
class Out:
    def __init__(self, path):
        self.zf = zipfile.ZipFile(path) if path.endswith(".zip") else None
        self.root = path
        self.names = self.zf.namelist() if self.zf else None

    def read(self, rel):
        if self.zf:
            cand = [n for n in self.names if n.endswith("/" + rel) or n == rel]
            return self.zf.read(cand[0]).decode("utf-8", "replace") if cand else None
        p = os.path.join(self.root, rel)
        return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else None

    def jobs(self):
        return [l.split("\t") for l in self.read("jobs.tsv").splitlines()]


def selftest():
    t = ">Feature lcl|ORF1_x\n1\t99\tCDS\n\t\t\tcodon_start\t1\n<5\t>120\tCDS\n200\t150\tCDS\n"
    r = parse_feature_table(t)["lcl|ORF1_x"]
    assert r == [("+", 1, 99, 0, 0), ("+", 5, 120, 1, 1), ("-", 150, 200, 0, 0)], r
    print("selftest ok")


def probe(path):
    o = Out(path)
    for jid, group, infile, tmo, useout, args in o.jobs():
        if jid in ("default__designed__fmt0", "default__designed__fmt3", "default__designed__fmt1", "default__designed__fmt2"):
            t = o.read("runs/%s/orfs.txt" % jid) or ""
            print("=== %s\n%s\n" % (jid, t[:1500]))


def run(path):
    o = Out(path)
    print("job_id\targs\tseq_id\tsequence\torfs")
    seqs = {}
    for f in ("designed.fa", "random.fa", "circular.fa", "region.fa"):
        seqs[f] = read_fasta(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs", f)).read())
    for jid, group, infile, tmo, useout, args in o.jobs():
        if not jid.endswith("fmt3") or group in ("misc", "edge", "circedge", "single"):
            continue
        t3 = o.read("runs/%s/orfs.txt" % jid)
        if t3 is None:
            continue
        ft = parse_feature_table(t3)
        t0 = o.read("runs/%s0/orfs.txt" % jid[:-1])
        if t0 is not None:
            fa = parse_orf_fasta(t0)
            for k, v in fa.items():  # coordinates (and count) of fmt 0 and fmt 3 must agree
                key = [x for x in ft if x.endswith("_" + k) or x == k]
                got = sorted((s, a, b) for s, a, b, _, _ in ft[key[0]]) if key else None
                assert got == sorted((s, a, b) for s, a, b, _, _ in v), (jid, k, got, v)
        base = os.path.basename(infile)
        for sid, sq in seqs.get(base, {}).items():
            key = [x for x in ft if x.endswith("_" + sid) or x == sid]
            orfs = ft[key[0]] if key else []
            print("%s\t%s\t%s\t%s\t%s" % (jid, args, sid, sq, ";".join(",".join(map(str, x)) for x in orfs) or "-"))


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    elif "--probe" in sys.argv:
        probe(sys.argv[-1])
    elif "--run" in sys.argv:
        run(sys.argv[-1])
    else:
        print(__doc__)
