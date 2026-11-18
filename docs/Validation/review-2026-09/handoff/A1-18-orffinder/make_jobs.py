#!/usr/bin/env python3
"""Run matrix of the NCBI ORFfinder oracle package: prints one TSV line per ORFfinder invocation.

columns: id  group  infile(relative to the package dir, '-' = none)  timeout_s  useout(1 = append '-out <dir>/orfs.txt')  args

  python3 make_jobs.py            full matrix
  QUICK=1 python3 make_jobs.py    reduced matrix (smoke test, a few minutes)
  SINGLE=0 python3 make_jobs.py   without the one-sequence-per-process runs
  python3 make_jobs.py --summary  counts per group

Every option of ORFfinder 0.4.3 USAGE.txt is exercised: -h -help -xmlhelp -logfile -conffile -version -version-full -dryrun
-in -id -b -e -c -g -s -ml -n -strand -out -outfmt (0 ORF list FASTA, 1 CDS FASTA, 2 text ASN.1, 3 feature table).
NOTE: USAGE.txt has no `-partial` option; partial ORFs (running off either end of the sequence) are covered by the
`partial*` designs of the inputs, `-s 2`, `-c` and the -b/-e grid (a region end cuts ORFs).
"""
import os
import sys

QUICK = os.environ.get("QUICK") == "1"
SINGLE = os.environ.get("SINGLE", "1") != "0"
jobs = []


def sel(seq):
    seq = list(seq)
    return seq[::4] if QUICK else seq


def job(jid, group, infile, args, timeout=300, useout=True):
    jobs.append((jid, group, infile, timeout, 1 if useout else 0, " ".join(str(a) for a in args)))


def fmts(*f):
    return f


MAIN = ["inputs/designed.fa", "inputs/random.fa"]
short = lambda p: os.path.basename(p)[:-3]  # noqa: E731
B = {"false": "0", "true": "1"}

# ---- 1. no options at all (defaults: -g 1 -s 1 -ml 75 -n false -strand both -c false -outfmt 0) -----------------
for f in MAIN:
    for fmt in (0, 1, 2, 3):
        job("default__%s__fmt%d" % (short(f), fmt), "default", f, ["-outfmt", fmt])
for fmt in (0, 3):
    job("default_stdout__designed__fmt%d" % fmt, "default", "inputs/designed.fa", ["-outfmt", fmt], useout=False)

# ---- 2. factorial: genetic code x start mode x nested x strand, -ml 30 --------------------------------------------
for g in sel((1, 2, 4, 11)) if QUICK else (1, 2, 4, 11):
    for s in (0, 1, 2):
        for n in ("false", "true"):
            for strand in ("both", "plus", "minus"):
                if strand != "both" and g not in (1, 11):
                    continue
                if QUICK and strand != "both" and g != 1:
                    continue
                for f in MAIN:
                    base = ["-g", g, "-s", s, "-ml", 30, "-n", n, "-strand", strand]
                    key = "g%d_s%d_ml30_n%s_%s" % (g, s, B[n], strand)
                    for fmt in (0, 3):
                        job("fact__%s__%s__fmt%d" % (short(f), key, fmt), "factorial", f, base + ["-outfmt", fmt])
                    if strand == "both" and n == "false":
                        job("fact__%s__%s__fmt1" % (short(f), key), "factorial", f, base + ["-outfmt", 1])
                    if strand == "both" and g in (1, 11) and s in (0, 1, 2) and (n == "true" or s == 1):
                        job("fact__%s__%s__fmt2" % (short(f), key), "factorial", f, base + ["-outfmt", 2])

# ---- 3. minimum-length sweep (clamp: values below 30 become 30) ---------------------------------------------------
ML = (-5, 0, 1, 10, 29, 30, 31, 45, 60, 74, 75, 76, 90, 150, 300, 1000)
for ml in sel(ML) if QUICK else ML:
    for n in ("false", "true"):
        for f in MAIN:
            for fmt in (0, 3):
                job("ml__%s__ml%d_n%s__fmt%d" % (short(f), ml, B[n], fmt), "minlength", f,
                    ["-ml", ml, "-n", n, "-outfmt", fmt])
for ml in (0, 29, 30, 31, 45, 75):
    for s in (0, 2):
        for fmt in (0, 3):
            job("ml__designed__ml%d_s%d__fmt%d" % (ml, s, fmt), "minlength", "inputs/designed.fa", ["-ml", ml, "-s", s, "-outfmt", fmt])

# ---- 4. every genetic code 1..31 plus out-of-range values (designed set only) -----------------------------------
for g in range(0, 34):
    if QUICK and g not in (0, 3, 5, 27, 28, 31, 32):
        continue
    for s in (0, 1):
        for fmt in (0,) if s == 0 else (0, 3):
            job("gcode__designed__g%d_s%d__fmt%d" % (g, s, fmt), "gcode", "inputs/designed.fa", ["-g", g, "-s", s, "-ml", 30, "-outfmt", fmt])

# ---- 5. circular (-c) ------------------------------------------------------------------------------------------------
for c in ("true", "false"):
    for s in (0, 1, 2):
        for n in ("false", "true"):
            for strand in ("both", "plus", "minus"):
                if QUICK and strand != "both":
                    continue
                for fmt in (0, 3):
                    job("circ__circular__c%s_s%d_n%s_%s__fmt%d" % (B[c], s, B[n], strand, fmt), "circular", "inputs/circular.fa",
                        ["-c", c, "-s", s, "-ml", 30, "-n", n, "-strand", strand, "-outfmt", fmt])
for s in (0, 1, 2):
    for n in ("false", "true"):
        for fmt in (0, 3):
            job("circ__designed__ctrue_s%d_n%s__fmt%d" % (s, B[n], fmt), "circular", "inputs/designed.fa", ["-c", "true", "-s", s, "-ml", 30, "-n", n, "-outfmt", fmt])
for fmt in (0, 1, 2, 3):
    job("circ__circular__ctrue_g11_s1_ml75__fmt%d" % fmt, "circular", "inputs/circular.fa", ["-c", "true", "-g", 11, "-outfmt", fmt])
for c in ("t", "f", "yes", "no", "TRUE", "1", "maybe"):
    job("circ_boolspelling__c%s" % c, "circular", "inputs/circular.fa", ["-c", c, "-ml", 30, "-outfmt", 3])
for n in ("t", "f", "yes", "no", "TRUE", "1", "maybe"):
    job("nested_boolspelling__n%s" % n, "minlength", "inputs/designed.fa", ["-n", n, "-ml", 30, "-outfmt", 3])

# ---- 6. region -b / -e (1-based, -e 0 = to the end) ----------------------------------------------------------------
PAIRS = [(1, 0), (1, 300), (1, 150), (151, 0), (151, 300), (100, 200), (2, 299), (31, 270), (1, 100), (101, 200), (201, 300),
         (50, 250), (2, 0), (1, 299), (60, 70), (300, 300), (1, 301), (1, 1000), (0, 0), (0, 100), (250, 100), (100, 100), (5, 3), (-5, 100)]
CORE = [(1, 150), (151, 0), (100, 200), (2, 299), (31, 270), (50, 250), (1, 100), (201, 300), (60, 70), (1, 301)]
for b, e in sel(PAIRS) if QUICK else PAIRS:
    for fmt in (0, 1, 3):
        job("region__b%d_e%d_s1_nf_both__fmt%d" % (b, e, fmt), "region", "inputs/region.fa", ["-b", b, "-e", e, "-ml", 30, "-outfmt", fmt])
    for fmt in (0, 3):
        job("region__b%d_e%d_s1_nt_both__fmt%d" % (b, e, fmt), "region", "inputs/region.fa", ["-b", b, "-e", e, "-ml", 30, "-n", "true", "-outfmt", fmt])
for b, e in sel(CORE) if QUICK else CORE:
    for s in (0, 2):
        for n in ("false", "true"):
            job("region__b%d_e%d_s%d_n%s_both__fmt0" % (b, e, s, B[n]), "region", "inputs/region.fa", ["-b", b, "-e", e, "-s", s, "-ml", 30, "-n", n, "-outfmt", 0])
            job("region__b%d_e%d_s%d_n%s_both__fmt3" % (b, e, s, B[n]), "region", "inputs/region.fa", ["-b", b, "-e", e, "-s", s, "-ml", 30, "-n", n, "-outfmt", 3])
    for strand in ("plus", "minus"):
        for fmt in (0, 3):
            job("region__b%d_e%d_s1_nf_%s__fmt%d" % (b, e, strand, fmt), "region", "inputs/region.fa", ["-b", b, "-e", e, "-ml", 30, "-strand", strand, "-outfmt", fmt])
    for fmt in (0, 3):
        job("region__b%d_e%d_ctrue__fmt%d" % (b, e, fmt), "region", "inputs/region.fa", ["-b", b, "-e", e, "-c", "true", "-ml", 30, "-outfmt", fmt])
for b, e in ((1, 150), (151, 0), (100, 200)):
    job("region__designed__b%d_e%d__fmt3" % (b, e), "region", "inputs/designed.fa", ["-b", b, "-e", e, "-ml", 30, "-outfmt", 3])

# ---- 7. one sequence per process (designed categories where -n / partial semantics matter) --------------------------
SINGLE_CATS = ("xnest", "xoverlap", "strandnest", "strandoverlap", "strandequal", "partial", "inframe")
if SINGLE:
    ids = []
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs", "manifest.tsv")).read().splitlines()[1:]:
        cid, f, cat, _, _ = line.split("\t")
        if f == "designed.fa" and cat.replace("_rc", "") in SINGLE_CATS:
            ids.append(cid)
    if QUICK:
        ids = ids[::6]
    for cid in ids:
        for key, args, fl in (("s1_nf", ["-s", 1, "-n", "false"], (0, 3)), ("s1_nt", ["-s", 1, "-n", "true"], (0, 3)),
                              ("s2_nt", ["-s", 2, "-n", "true"], (0,)), ("s0_nf", ["-s", 0, "-n", "false"], (0,))):
            for fmt in fl:
                job("single__%s__%s_ml30__fmt%d" % (cid, key, fmt), "single", "work/single/%s.fa" % cid, args + ["-ml", 30, "-outfmt", fmt], timeout=60)

# ---- 8. edge-case files, one sequence per process ---------------------------------------------------------------------
EDGE = sorted(os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs", "edge")))
for name in EDGE:
    for key, args in (("default", []), ("ml30_s2_nt", ["-ml", 30, "-s", 2, "-n", "true"]), ("ml30_ctrue", ["-ml", 30, "-c", "true"])):
        for fmt in (0, 3):
            job("edge__%s__%s__fmt%d" % (name[:-3], key, fmt), "edge", "inputs/edge/" + name, args + ["-outfmt", fmt], timeout=60)

# ---- 9. circular genomes without a stop codon (timeout: a circular ORF scan could loop) ---------------------------------
CEDGE = sorted(os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs", "circular_edge")))
for name in CEDGE:
    for c in ("true", "false"):
        for s in (0, 2):
            job("circedge__%s__c%s_s%d__fmt3" % (name[:-3], B[c], s), "circedge", "inputs/circular_edge/" + name, ["-c", c, "-s", s, "-ml", 30, "-outfmt", 3], timeout=60)

# ---- 10. other options / error paths ---------------------------------------------------------------------------------------
D = "inputs/designed.fa"
job("misc__h", "misc", "-", ["-h"], 60, False)
job("misc__help", "misc", "-", ["-help"], 60, False)
job("misc__xmlhelp", "misc", "-", ["-xmlhelp"], 60, False)
job("misc__version", "misc", "-", ["-version"], 60, False)
job("misc__version_full", "misc", "-", ["-version-full"], 60, False)
job("misc__dryrun_in", "misc", D, ["-dryrun"], 60, False)
job("misc__dryrun_noinput", "misc", "-", ["-dryrun"], 60, False)
job("misc__noinput", "misc", "-", [], 60, False)
job("misc__missing_file", "misc", "inputs/does_not_exist.fa", ["-outfmt", 0], 60)
job("misc__logfile", "misc", D, ["-ml", 30, "-logfile", "orffinder.log", "-outfmt", 0], 120)
job("misc__conffile_missing", "misc", D, ["-conffile", "no_such.ini", "-ml", 30, "-outfmt", 0], 60)
for v in (4, 5, -1):
    job("misc__outfmt%d" % v, "misc", D, ["-ml", 30, "-outfmt", v], 60)
for v in ("x", "left", "BOTH", "Plus"):
    job("misc__strand_%s" % v, "misc", D, ["-ml", 30, "-strand", v, "-outfmt", 0], 60)
for v in (-1, 3, 99):
    job("misc__s%d" % v, "misc", D, ["-ml", 30, "-s", v, "-outfmt", 0], 60)
job("misc__ml_text", "misc", D, ["-ml", "abc", "-outfmt", 0], 60)
job("misc__id_lambda", "misc", "-", ["-id", "NC_001416.1", "-ml", 100, "-outfmt", 3], 120)
job("misc__id_and_in", "misc", D, ["-id", "NC_001416.1", "-ml", 30, "-outfmt", 0], 120)
job("misc__in_stdin_dash", "misc", "-", ["-in", "-", "-ml", 30, "-outfmt", 0], 30)

if __name__ == "__main__":
    if "--summary" in sys.argv:
        from collections import Counter
        c = Counter(j[1] for j in jobs)
        for k, v in c.items():
            print("%-12s %5d" % (k, v))
        print("%-12s %5d" % ("TOTAL", len(jobs)))
    else:
        ids = [j[0] for j in jobs]
        assert len(ids) == len(set(ids)), "duplicate job ids"
        for j in jobs:
            print("\t".join(str(x) for x in j))
