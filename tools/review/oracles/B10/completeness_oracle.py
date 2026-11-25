#!/usr/bin/env python3
"""B10 ASSEMBLY-COMPLETE-001 oracle (run with `python3 -I`).

Independent re-implementations of the sourced definitions used by
GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers and AssessCompleteness:

* Merqury (marbl/merqury master):
  - build/filt.sh -> kmerHistToPloidyDepth.jar (HistToPloidyDepth.java, getPloidyDepthTable):
    the solid threshold `filt` = effectiveMinDepth (transcribed below, full loop incl. the parts that do
    not influence effectiveMinDepth are omitted because they run only after it is fixed).
  - eval/spectra-cn.sh: completeness = distinct(asm ∩ read.gt$filt) / distinct(read.gt$filt) * 100.
  - eval/qv.sh / spectra-cn.sh: ERROR = 1-(1-ASM_ONLY/TOTAL)^(1/k), QV = -10*log10(ERROR)
    (awk transcribed literally), ASM_ONLY/TOTAL = `meryl statistics` "present" (instances).
  - meryl canonical counting re-implemented here (min of k-mer and reverse complement, ACGT order,
    windows containing a non-ACGT base skipped).
* GenomeScope 1.0 (schatzlab/genomescope genomescope.R): error_rate = uniroot of
  f1(x) = sum_{i=1..k} (1-x)^(k-i) x^i choose(k,i) * total_kmers - total_error_kmers  (solved with
  scipy brentq, independent of the closed form used in C#); total_len = (total - error) / peak, R round.
* BUSCO (gitlab ezlab/busco src/busco/busco_tools/hmmer.py): percentages over n and the one-line
  summary "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}" with round(x, 1).
* compleasm v0.2.6 (huangnengCSU/compleasm compleasm.py): pass filter
  (Protein_mapped_rate >= min_length_percent) | (Identity >= min_identity), record_1st_gene_label
  (by_length=False) and get_region_clusters, transcribed.
  Hits: Biopython 1.88 PairwiseAligner local alignment (BLAST+ blastn scoring +2/-3, gap open 5,
  extend 2 -> open_gap_score -7, extend_gap_score -2) of the marker and its reverse complement against
  each whole sequence; further hits are found by masking the aligned region with N and re-aligning while
  the best score is >= 2*min(28, m) (= the score of the exact seed word the C# search requires); a hit is kept
  only when it shares an exact ACGT word of min(28, marker length) with the marker (megablast word_size).
"""
import json
import math
import random
import sys

from Bio.Align import PairwiseAligner
from scipy.optimize import brentq

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


# ---------------------------------------------------------------- Merqury / GenomeScope
def merqury_filt(hist):
    """HistToPloidyDepth.getPloidyDepthTable -> effectiveMinDepth. hist = [(depth, count)] ascending."""
    prev_count = sys.float_info.max
    prev_d = -1.0
    line_count = 0
    count_min = 0
    effective_min_depth = 0
    for depth, count in hist:
        count = float(count)
        if line_count == 0:
            prev_d = count - prev_count
            prev_count = count
            line_count += 1
            continue
        d = count - prev_count
        if d * prev_d < 0:
            if d > prev_d:
                count_min += 1
                if count_min == 1:
                    effective_min_depth = depth - 1
        # (ploidy peak bookkeeping omitted: it never changes effectiveMinDepth)
        if count_min >= 1:
            break
        prev_count = count
        prev_d = d
    return effective_min_depth


def meryl_count(seqs, k, canonical=True):
    counts = {}
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - k + 1):
            w = s[i:i + k]
            if any(c not in "ACGT" for c in w):
                continue
            if canonical:
                w = min(w, rc(w))
            counts[w] = counts.get(w, 0) + 1
    return counts


def genomescope_error_rate(error_kmers, total_kmers, k):
    if error_kmers == 0:
        return 0.0

    def f1(x):
        return sum((1 - x) ** (k - i) * x ** i * math.comb(k, i) for i in range(1, k + 1)) * total_kmers - error_kmers

    return brentq(f1, 0.0, 1.0, xtol=1e-15, rtol=1e-15)


def merqury_qv(asm_only, total, k):
    """awk: ERROR = 1-(1-$1/$2)^(1/k); QV = -10*log(ERROR)/log(10)."""
    if total == 0:
        return None, None
    err = 1 - (1 - asm_only / total) ** (1 / k)
    qv = math.inf if err == 0 else -10 * math.log(err) / math.log(10)
    return err, qv


def kmer_estimate(spectrum, assembly=None, expected_coverage=0, canonical=True):
    k = len(next(iter(spectrum)))
    hist = {}
    for c in spectrum.values():
        hist[c] = hist.get(c, 0) + 1
    hist = sorted(hist.items())
    filt = merqury_filt(hist)
    total = sum(c * h for c, h in hist)
    error = sum(c * h for c, h in hist if c <= filt)
    above = [(c, h) for c, h in hist if c > filt]
    peak = expected_coverage or max(above, key=lambda t: (t[1], -t[0]))[0]  # which(y==max(y))[1]
    out = {
        "k": k, "filt": filt, "peak": peak, "total": total, "error": error,
        "errorRate": genomescope_error_rate(error, total, k),
        "errorRateClosedForm": 1 - (1 - error / total) ** (1 / k),
        "size": round((total - error) / peak),
        "solid": sum(h for c, h in above),
    }
    if assembly is not None:
        asm = meryl_count(assembly, k, canonical)
        solid = {km for km, c in spectrum.items() if c > filt}
        asm_solid = len(set(asm) & solid)
        asm_total = sum(asm.values())
        asm_only = sum(c for km, c in asm.items() if km not in spectrum)
        err, qv = merqury_qv(asm_only, asm_total, k)
        out.update({"asmSolid": asm_solid, "completeness": asm_solid / len(solid),
                    "completenessPercent": 100 * asm_solid / len(solid),
                    "asmTotal": asm_total, "asmOnly": asm_only, "asmErrorRate": err, "qv": qv})
    return out


# ---------------------------------------------------------------- BUSCO / compleasm
ALIGNER = PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                          open_gap_score=-7, extend_gap_score=-2)
MIN_SCORE = 56


def best_local(query, target):
    if not target:
        return None
    alns = ALIGNER.align(target, query)
    if alns.score < 2 * min(MIN_SCORE // 2, len(query)):  # score of one exact seed word
        return None
    a = alns[0]
    (ts, te), (qs, qe) = (a.coordinates[0][0], a.coordinates[0][-1]), (a.coordinates[1][0], a.coordinates[1][-1])
    t_aln, q_aln = a[0], a[1]
    cols = len(t_aln)
    ident = sum(1 for x, y in zip(t_aln, q_aln) if x == y and x != "-")
    return int(ts), int(te), int(qs), int(qe), ident / cols


def find_hits(marker, assembly):
    marker = marker.upper()
    m = len(marker)
    hits = []
    for strand, query in (("+", marker), ("-", rc(marker))):
        for sid, seq in assembly:
            seq = seq.upper()
            while True:
                r = best_local(query, seq)
                if r is None:
                    break
                ts, te, qs, qe, ident = r
                # megablast seed requirement: the hit must share an exact ACGT word of w = min(28, m)
                w = min(28, m)
                words = {query[i:i + w] for i in range(m - w + 1) if set(query[i:i + w]) <= set("ACGT")}
                region = seq[ts:te]
                if any(region[i:i + w] in words for i in range(len(region) - w + 1)):
                    hits.append({"id": sid, "strand": strand, "start": ts, "end": te - 1,
                                 "identity": ident, "coverage": (qe - qs) / m})
                seq = seq[:ts] + "N" * (te - ts) + seq[te:]
    return hits


def get_region_clusters(regions):  # compleasm, verbatim
    sorted_regions = sorted(regions, key=lambda x: x[0], reverse=False)
    clusters = []
    for (start, stop) in sorted_regions:
        if not clusters:
            clusters.append([start, stop])
        else:
            last_cluster = clusters[-1]
            if last_cluster[0] <= start <= last_cluster[1]:
                clusters[-1][0] = min(last_cluster[0], start)
                clusters[-1][1] = max(last_cluster[1], stop)
            else:
                clusters.append([start, stop])
    return clusters


def classify(hits, min_identity=0.9, min_complete=0.9, min_length_percent=0.6):
    valid = [h for h in hits if h["coverage"] >= min_length_percent or h["identity"] >= min_identity]
    if not valid:
        return "Missing"
    if len(valid) == 1:
        return "Single" if valid[0]["coverage"] >= min_complete else "Fragmented"
    complete = [h for h in valid if h["coverage"] >= min_complete]
    if not complete:
        return "Fragmented"
    if len(complete) == 1:
        return "Single"
    if len({h["id"] for h in complete}) > 1:
        return "Duplicated"
    clusters = get_region_clusters([(h["start"], h["end"] + 1) for h in complete])  # PAF half-open
    return "Single" if len(clusters) == 1 else "Duplicated"


def busco_summary(labels):
    n = len(labels)
    s = labels.count("Single")
    d = labels.count("Duplicated")
    f = labels.count("Fragmented")
    mi = n - s - d - f
    pct = lambda x: abs(round(x / n * 100, 1))
    return {"n": n, "S": s, "D": d, "F": f, "M": mi,
            "C%": (s + d) / n * 100, "D%": d / n * 100,
            "notation": "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}".format(pct(s + d), pct(s), pct(d), pct(f), pct(mi), n)}


# ---------------------------------------------------------------- scenarios
def rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate, rng):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in s)


def scenarios():
    rng = random.Random(41)
    g = rnd(300, rng)
    flank = lambda n: rnd(n, rng)
    cases = {}
    L, R = flank(100), flank(100)
    cases["exact_forward"] = (g, [("c1", L + g + R)])
    cases["reverse_strand"] = (g, [("c1", L + rc(g) + R)])
    cases["two_contigs"] = (g, [("c1", L + g + R), ("c2", flank(50) + g + flank(50))])
    cases["tandem_same_contig"] = (g, [("c1", L + g + g + R)])
    cases["plus_and_minus_same_contig"] = (g, [("c1", L + g + flank(80) + rc(g) + R)])
    cases["fragment_70pct"] = (g, [("c1", L + g[:210] + R)])
    cases["fragment_40pct"] = (g, [("c1", L + g[:120] + R)])
    cases["absent"] = (g, [("c1", flank(500))])
    cases["one_bp_deletion"] = (g, [("c1", L + g[:150] + g[151:] + R)])
    cases["complete_plus_fragment"] = (g, [("c1", L + g + R), ("c2", flank(40) + g[:120] + flank(40))])
    sub = list(g)
    for i in range(10, 300, 20):  # 15 substitutions, every 20 bp -> identity 0.95
        sub[i] = COMP[sub[i]]
    cases["substitutions_5pct"] = (g, [("c1", L + "".join(sub) + R)])
    out = {}
    for name, (marker, asm) in cases.items():
        hits = find_hits(marker, asm)
        out[name] = {"marker": marker, "assembly": asm, "hits": hits, "label": classify(hits)}
    return out


def random_busco(count, seed):
    rng = random.Random(seed)
    res = []
    for _ in range(count):
        g = rnd(rng.randint(60, 240), rng)
        asm = []
        for c in range(rng.randint(1, 3)):
            parts = [rnd(rng.randint(0, 60), rng)]
            for _ in range(rng.randint(0, 2)):
                kind = rng.random()
                piece = g
                if kind < 0.3:
                    a = rng.randint(0, len(g) // 2)
                    piece = g[a:a + rng.randint(len(g) // 4, len(g))]
                elif kind < 0.5:
                    piece = mutate(g, rng.choice([0.01, 0.03]), rng)
                if rng.random() < 0.4:
                    piece = rc(piece)
                parts += [piece, rnd(rng.randint(0, 80), rng)]
            asm.append(("c%d" % c, "".join(parts)))
        hits = find_hits(g, asm)
        res.append({"marker": g, "assembly": asm, "label": classify(hits)})
    return res


def random_kmers(count, seed):
    rng = random.Random(seed)
    res = []
    for _ in range(count):
        k = rng.choice([3, 4, 5, 7, 11])
        genome = rnd(rng.randint(30, 400), rng)
        reads = []
        cov = rng.randint(4, 25)
        for _ in range(cov * len(genome) // 50 + 1):
            a = rng.randint(0, max(0, len(genome) - 50))
            r = genome[a:a + 50]
            if rng.random() < 0.5:
                r = mutate(r, 0.01, rng)
            reads.append(r if rng.random() < 0.5 else rc(r))
        spectrum = meryl_count(reads, k, True)
        assembly = [genome[: len(genome) // 2], mutate(genome[len(genome) // 2:], 0.02, rng)]
        expected = rng.choice([0, 0, 0, cov])
        result = kmer_estimate(spectrum, assembly, expected)
        if result["qv"] == math.inf:
            result["qv"] = "inf"  # JSON has no infinity
        res.append({"spectrum": spectrum, "assembly": assembly, "expected": expected, "result": result})
    return res


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "scenarios":
        sc = scenarios()
        for name, v in sc.items():
            print(name, v["label"], [(h["id"], h["strand"], h["start"], h["end"], round(h["identity"], 6), round(h["coverage"], 6)) for h in v["hits"]])
        print("summary", busco_summary([v["label"] for v in sc.values()]))
        if len(sys.argv) > 2:
            json.dump(sc, open(sys.argv[2], "w"))
    elif mode == "kmers":
        cases = {
            # (spectrum, assembly, expected coverage)
            "valley": ({"AAA": 1, "AAC": 1, "AAG": 1, "ACA": 2, "ACC": 5, "ACG": 5, "ACT": 5, "AGA": 4, "AGC": 6}, ["AAAACCAGCAATT"], 0),
            "monotone": ({"AAA": 4, "AAC": 4, "AAG": 2, "AAT": 1}, None, 0),
            "single_count": ({"AAA": 30, "AAC": 30}, ["AAACC"], 0),
            "rising_start": ({"AAA": 2, "AAC": 3, "AAG": 3}, None, 0),
            "expected_cov": ({"AAA": 4, "AAC": 4, "AAG": 2, "AAT": 1}, None, 2),
        }
        for name, (sp, asm, e) in cases.items():
            print(name, kmer_estimate(sp, asm, e))
    elif mode == "random-busco":
        json.dump(random_busco(int(sys.argv[2]), int(sys.argv[3])), open(sys.argv[4], "w"))
    elif mode == "random-kmers":
        json.dump(random_kmers(int(sys.argv[2]), int(sys.argv[3])), open(sys.argv[4], "w"), allow_nan=True)
    elif mode == "notation":
        for labels in (["Single"] * 2 + ["Duplicated"] + ["Missing"], ["Single"] * 1 + ["Fragmented"] * 2,
                       ["Single"] * 199 + ["Missing"] * 1801, ["Duplicated"] + ["Missing"] * 1999):
            print(busco_summary(labels))
