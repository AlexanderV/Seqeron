#!/usr/bin/env python3
"""Deterministic test-input generator for the NCBI ORFfinder oracle package (A1-18 / A1-19).

Pure standard library; own splitmix64 PRNG so the output is byte-identical on every Python 3.x / platform.
Writes ./inputs/*.fa, ./inputs/edge/*.fa, ./inputs/circular_edge/*.fa, ./inputs/manifest.tsv, ./inputs/inputs.sha256.

  python3 make_inputs.py            # (re)generate into ./inputs
  python3 make_inputs.py --check    # regenerate in memory and compare with the committed ./inputs (exit 1 on any difference)

Files (see manifest.tsv for one row per sequence: id, file, category, length, note):
  designed.fa  hand-designed sequences (nesting in frame and across frames/strands, strand overlaps, ORFs running off
               both ends, every start codon of tables 1/2/4/11, every stop of those tables, ORF lengths 6..90 nt
               for the -ml clamp) plus the reverse complement of every forward design
  random.fa    random sequences: iid bases at several GC levels, gene-dense mosaics of implanted ORFs on both
               strands, and point-mutated copies of those
  circular.fa  ORFs (both strands, all codon phases) cut by the sequence origin, for -c t
  region.fa    60 sequences of exactly 300 nt with ORFs straddling the positions used by the -b/-e grid
  edge/*.fa    one sequence per file: length 0..3, 6, 29, 30, lower case, N/IUPAC, U, gaps, CRLF, long header, ...
  circular_edge/*.fa  circular genomes without any stop codon (a circular ORFfinder could in principle loop): run with a timeout
"""
import hashlib
import os
import sys

M64 = (1 << 64) - 1
BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# Union of the stop / start codons of NCBI tables 1, 2, 4, 11 (the -g values of the run matrix).
STOPS_ALL = {"TAA", "TAG", "TGA", "AGA", "AGG"}
STARTS_ALL = ["ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA", "TTA"]
ALL_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
NEUTRAL = [c for c in ALL_CODONS if c not in STOPS_ALL and c not in STARTS_ALL]  # 51 codons: neither start nor stop


class Rng:
    def __init__(self, seed):
        self.s = seed & M64

    def next(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & M64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
        return z ^ (z >> 31)

    def below(self, n):
        return self.next() % n

    def between(self, lo, hi):  # inclusive
        return lo + self.below(hi - lo + 1)

    def rand(self):
        return (self.next() >> 11) / float(1 << 53)

    def choice(self, seq):
        return seq[self.below(len(seq))]


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def neutral(rng, n):
    return "".join(rng.choice(NEUTRAL) for _ in range(n))


def noise(rng, n, gc=0.5):
    out = []
    for _ in range(n):
        if rng.rand() < gc:
            out.append("G" if rng.below(2) else "C")
        else:
            out.append("A" if rng.below(2) else "T")
    return "".join(out)


# ---------------------------------------------------------------------------------------------------------
# Tiny definition-based ORF scan, used only to verify that a designed sequence contains the designed ORFs.
# ---------------------------------------------------------------------------------------------------------
def scan_strand(s, stops, is_start, nested, include_open):
    out = []
    n = len(s)
    for f in range(3):
        pend = []
        i = f
        while i + 3 <= n:
            c = s[i:i + 3]
            if c in stops:
                out.extend((p, i + 3, True) for p in pend)
                pend = []
            elif is_start(c) and (nested or not pend):
                pend.append(i)
            i += 3
        if include_open:
            out.extend((p, i, False) for p in pend)
    return out


def orfs_forward_coords(seq, is_start, nested=False, include_open=False, stops=("TAA", "TAG", "TGA")):
    """[(strand, start, end_exclusive, terminated)] in forward coordinates of seq."""
    n = len(seq)
    res = [("+", a, b, t) for a, b, t in scan_strand(seq, set(stops), is_start, nested, include_open)]
    r = rc(seq)
    res += [("-", n - b, n - a, t) for a, b, t in scan_strand(r, set(stops), is_start, nested, include_open)]
    return res


def atg_only(c):
    return c == "ATG"


def any_sense(c):
    return c not in ("TAA", "TAG", "TGA")


class Case:
    def __init__(self, cid, cat, seq, note=""):
        assert all(ch in BASES for ch in seq), cid
        self.id, self.cat, self.seq, self.note = cid, cat, seq, note


# ---------------------------------------------------------------------------------------------------------
# designed.fa
# ---------------------------------------------------------------------------------------------------------
def orf_text(rng, body_codons, start="ATG", stop="TAA"):
    return start + neutral(rng, body_codons) + stop


def place(cv, a, strand, text):
    """Write an ORF text (5'->3' on its own strand) at forward positions [a, a+len)."""
    s = text if strand == "+" else rc(text)
    cv[a:a + len(s)] = list(s)


def has_orf(seq, strand, a, b, is_start=atg_only):
    return (strand, a, b, True) in orfs_forward_coords(seq, is_start, nested=False)


def build_with_orfs(rng, length, specs, is_start=atg_only, tries=4000, flank_gc=None):
    """specs = [(strand, a, body_codons, start, stop)], written in order on a 'CCC' background (later ones overwrite).
    Rejection-sampled until every ORF is present with exactly its designed forward span."""
    for _ in range(tries):
        cv = list(("CCC" * (length // 3 + 1))[:length])
        spans = []
        for strand, a, nb, start, stop in specs:
            t = orf_text(rng, nb, start, stop)
            place(cv, a, strand, t)
            spans.append((strand, a, a + len(t)))
        seq = "".join(cv)
        if all(has_orf(seq, st, a, b, is_start) for st, a, b in spans):
            return seq, spans
    raise RuntimeError("could not realise design %r" % (specs,))


def designed(rng):
    cases = []

    def add(cid, cat, seq, note=""):
        cases.append(Case(cid, cat, seq, note))

    pre = "CCC" * 4
    post = "CCC" * 4

    # ---- A. in-frame structure ------------------------------------------------------------------------
    for k in (2, 3, 4):
        s = pre + "".join("ATG" + neutral(rng, rng.between(3, 9)) for _ in range(k)) + neutral(rng, 5) + "TAA" + post
        add("inframe_atg_x%d" % k, "inframe", s, "%d in-frame ATGs sharing one TAA" % k)
    add("inframe_alt_then_atg", "inframe", pre + "CTG" + neutral(rng, 6) + "ATG" + neutral(rng, 6) + "TAA" + post,
        "CTG start upstream of ATG: -s 0 reports the ATG ORF, -s 1 the CTG ORF (longest)")
    add("inframe_ttg_ctg_atg", "inframe", pre + "TTG" + neutral(rng, 4) + "CTG" + neutral(rng, 4) + "ATG" + neutral(rng, 4) + "TAG" + post,
        "chain of three starts")
    for st in STARTS_ALL:
        add("start_%s" % st, "start", pre + st + neutral(rng, 14) + "TGA" + post,
            "only start codon %s; reported by -s 1 depending on the -g table, by -s 2 as a stop-to-stop region" % st)
    for st in STARTS_ALL:
        add("start_%s_then_atg" % st, "start", pre + st + neutral(rng, 8) + "ATG" + neutral(rng, 8) + "TAA" + post,
            "alt start %s upstream of an in-frame ATG" % st)
    for stop in sorted(STOPS_ALL):
        add("stop_%s" % stop, "stop", pre + "ATG" + neutral(rng, 15) + stop + post, "ORF ended by %s (stop in tables: see -g)" % stop)
    add("tga_readthrough", "stop", pre + "ATG" + neutral(rng, 6) + "TGA" + neutral(rng, 6) + "TAA" + post, "TGA = Trp in tables 2/4: one ORF there, two in table 1/11")
    add("aga_agg_pair", "stop", pre + "ATG" + neutral(rng, 5) + "AGA" + neutral(rng, 4) + "ATG" + neutral(rng, 4) + "AGG" + post, "AGA/AGG stops in table 2 only")
    add("every_stop_in_a_row", "stop", pre + "ATG" + "TAATAGTGAAGAAGG" + post, "ORF of 6 nt then four more stop codons")

    # ---- B. cross-frame nesting on the same strand (Ignore nested ORFs, -n) --------------------------------
    for shift in (1, 2):
        for outer, inner in ((30, 8), (45, 12), (24, 5)):
            L = 3 * (outer + 2)
            spec = [("+", 0, outer, "ATG", "TAA"), ("+", 3 * 4 + shift, inner, "ATG", "TAA")]
            seq, spans = build_with_orfs(rng, L + 6, spec)
            add("xnest_s%d_%d_%d" % (shift, outer, inner), "xnest", seq,
                "inner ORF fully inside outer, other frame (shift %d), spans %s" % (shift, spans))
    for gap_start, gap_end in ((1, 1), (2, 2), (1, 2), (2, 1), (4, 4), (5, 5), (30, 30)):
        outer = 36
        a_in = gap_start
        # inner: ATG..stop; its end is `gap_end` nt before the end of the outer ORF
        end_outer = 3 * (outer + 2)
        inner_len = end_outer - gap_end - a_in
        if inner_len % 3 != 0 or inner_len < 9:
            # adjust the inner end so that its length is a multiple of three (keeps the requested start gap)
            inner_len -= inner_len % 3
        inner_body = inner_len // 3 - 2
        spec = [("+", 0, outer, "ATG", "TAA"), ("+", a_in, inner_body, "ATG", "TAA")]
        try:
            seq, spans = build_with_orfs(rng, end_outer + 6, spec, tries=12000)
        except RuntimeError:
            continue
        add("xnest_edge_a%d_e%d" % (gap_start, gap_end), "xnest", seq,
            "inner starts %d nt after / ends ~%d nt before the outer ORF; spans %s" % (gap_start, gap_end, spans))
    # an inner ORF starting 1 / 2 nt after the outer start is only possible with alternative starts (ATG cannot overlap ATG at offset 1/2)
    def any_alt_start(c):
        return c in STARTS_ALL

    for gap_start, outer_start, inner_start in ((1, "ATT", "TTG"), (1, "CTG", "TGG"), (2, "ATT", "TTG"), (1, "GTG", "TGT")):
        if inner_start not in STARTS_ALL:
            continue
        spec = [("+", 0, 30, outer_start, "TAA"), ("+", gap_start, 8, inner_start, "TAA")]
        try:
            seq, spans = build_with_orfs(rng, 3 * 32 + 6, spec, is_start=any_alt_start, tries=20000)
        except RuntimeError:
            continue
        add("xnest_altstart_a%d_%s_%s" % (gap_start, outer_start, inner_start), "xnest", seq,
            "inner %s ORF starts %d nt after outer %s ORF (alt-start modes only); spans %s" % (inner_start, gap_start, outer_start, spans))
    for overhang in (1, 2, 4, 5, 9):
        outer = 20
        end_outer = 3 * (outer + 2)
        a_in = 3 * 6 + 1 + (overhang % 2)
        inner_len = end_outer + overhang - a_in
        inner_len -= inner_len % 3
        spec = [("+", 0, outer, "ATG", "TAA"), ("+", a_in, inner_len // 3 - 2, "ATG", "TAA")]
        try:
            seq, spans = build_with_orfs(rng, end_outer + overhang + 12, spec, tries=12000)
        except RuntimeError:
            continue
        add("xoverlap_right_%d" % overhang, "xoverlap", seq, "inner ORF starts inside and ends beyond the outer ORF; spans %s" % spans)
        spec = [("+", 12 + 1 + (overhang % 2), 20, "ATG", "TAA"), ("+", 0, 14 + overhang // 3, "ATG", "TAA")]
        try:
            seq, spans = build_with_orfs(rng, 12 + 70 + 12, spec, tries=12000)
        except RuntimeError:
            continue
        add("xoverlap_left_%d" % overhang, "xoverlap", seq, "ORFs overlapping in two frames, neither contains the other; spans %s" % spans)

    # ---- C. opposite strands -----------------------------------------------------------------------------
    for rel in ("inner_rev", "inner_fwd"):
        for off, nb in ((12, 10), (13, 6), (14, 14)):
            if rel == "inner_rev":
                spec = [("+", 0, 40, "ATG", "TAA"), ("-", off, nb, "ATG", "TAA")]
            else:
                spec = [("-", 0, 40, "ATG", "TAA"), ("+", off, nb, "ATG", "TAA")]
            seq, spans = build_with_orfs(rng, 3 * 42 + 12, spec, tries=12000)
            add("strand_nest_%s_%d_%d" % (rel, off, nb), "strandnest", seq, "ORF of the other strand fully inside; spans %s" % spans)
    for off in (60, 63, 64, 65, 90):
        spec = [("+", 0, 24, "ATG", "TAA"), ("-", off, 20, "ATG", "TAA")]
        seq, spans = build_with_orfs(rng, off + 3 * 22 + 6, spec, tries=12000)
        add("strand_overlap_%d" % off, "strandoverlap", seq, "forward and reverse ORF overlapping partially; spans %s" % spans)
    spec = [("+", 12, 10, "ATG", "TAA"), ("-", 0, 40, "ATG", "TAA")]
    seq, spans = build_with_orfs(rng, 3 * 42 + 12, spec, tries=12000)
    add("strand_envelop_fwd_in_rev", "strandnest", seq, "forward ORF inside a longer reverse ORF; spans %s" % spans)

    # equal spans / shared boundaries across strands are only possible with s=2 (stop-to-stop regions): see the proof in README
    for body in (12, 20, 31):
        for _ in range(4000):
            x = neutral(rng, body)
            seq = "TAA" + "TTA" + x + "TAA" + "TTA"
            if has_orf(seq, "+", 3, 3 + 3 + 3 * body + 3, any_sense) and has_orf(seq, "-", 3, 3 + 3 + 3 * body + 3, any_sense):
                add("strand_equal_span_s2_%d" % body, "strandequal", seq,
                    "same forward span [3,%d) as a forward and as a reverse stop-to-stop region (-s 2 only)" % (3 + 3 + 3 * body + 3))
                break
    for body in (8, 16):
        for _ in range(4000):
            x = neutral(rng, body)
            seq = "TAA" + "TTA" + x + "TAA" + "CCCCCC"
            if has_orf(seq, "+", 3, 3 + 3 + 3 * body + 3, any_sense):
                add("strand_share_boundary_s2_%d" % body, "strandequal", seq, "forward stop-to-stop region whose first codon TTA is the reverse of a TAA stop")
                break

    # ---- D. partial ORFs (running off one or both ends) --------------------------------------------------
    for lead in (0, 1, 2):
        add("partial5_lead%d" % lead, "partial", "C" * lead + neutral(rng, 18) + "TAA" + post, "no start codon before the stop: 5'-partial stop-to-stop region from the sequence start")
        add("partial5_atg_lead%d" % lead, "partial", "C" * lead + neutral(rng, 6) + "ATG" + neutral(rng, 6) + "TAA" + post, "ATG inside, upstream part open")
    for trail in (0, 1, 2):
        add("partial3_trail%d" % trail, "partial", pre + "ATG" + neutral(rng, 18) + "C" * trail, "ATG ... no stop, ends %d nt after the last codon" % trail)
    add("partial3_halfstop_TA", "partial", pre + "ATG" + neutral(rng, 12) + "TA", "ends with the first two bases of a stop")
    add("partial3_halfstop_T", "partial", pre + "ATG" + neutral(rng, 12) + "T", "ends with the first base of a stop")
    add("partial3_atg_last_codon", "partial", "CCC" * 10 + "ATG", "ATG is the last complete codon (open ORF of one codon)")
    add("partial3_atg_plus1", "partial", "CCC" * 10 + "ATGC", "ATG as last complete codon + 1 trailing base")
    add("partial_both_neutral", "partial", neutral(rng, 30), "no start, no stop in the leading frame: runs off both ends")
    add("partial_both_atg", "partial", neutral(rng, 5) + "ATG" + neutral(rng, 24), "ATG inside, no stop: open at the 3' end, start missing at the 5' end is irrelevant")
    add("partial_atg_at_0_nostop", "partial", "ATG" + neutral(rng, 25), "ATG at position 0 and no stop")
    add("partial_whole_orf", "partial", "ATG" + neutral(rng, 20) + "TAA", "ORF covers the whole sequence exactly")
    add("partial_stop_at_0", "partial", "TAA" + neutral(rng, 20), "stop at position 0")
    add("partial_stop_last", "partial", neutral(rng, 20) + "TAA", "stop is the last codon, no start before it")
    add("partial_atg_pos1_whole", "partial", "C" + "ATG" + neutral(rng, 20) + "TAA", "ORF in frame 2 covering the sequence from base 2 to the end")

    # ---- E. -ml boundaries: exactly one forward ORF of L nt (stop included) in a CCC background -----------
    for L in range(6, 93, 3):
        cv = list("CCC" * 8) + list("ATG" + neutral(rng, L // 3 - 2) + "TAA") + list("CCC" * 8) if L > 6 else list("CCC" * 8) + list("ATGTAA") + list("CCC" * 8)
        add("mllen_%03d" % L, "mllen", "".join(cv), "single ATG..stop ORF of %d nt incl. stop (%d nt without)" % (L, L - 3))

    # ---- F. repeats / low complexity -----------------------------------------------------------------------
    add("polyA_60", "lowcomplex", "A" * 60, "no codon is a start or a stop")
    add("atg_x40", "lowcomplex", "ATG" * 40, "stop-free ATG repeat: open ORF in frame 1")
    add("atgtaa_x10", "lowcomplex", "ATGTAA" * 10, "ten 6-nt ORFs")
    add("taa_x20", "lowcomplex", "TAA" * 20, "only stops")
    add("atgatgtaa_pad", "lowcomplex", pre + "ATGATGTAA" + post, "")
    add("cag_x30", "lowcomplex", "CAG" * 30, "CAG repeat: reverse strand CTG start / stop-free")
    add("tttaa_x16", "lowcomplex", "TTTAA" * 16, "period-5 repeat: all frames mix of stops")
    add("gc_rich_stopfree", "lowcomplex", "GCC" * 3 + "ATG" + "GCCGCT" * 12 + "TGA" + "GCC" * 3, "GC-rich ORF")

    # ---- G. sequence-length boundaries ---------------------------------------------------------------------
    for n in (30, 31, 32, 33, 34, 35):
        add("len_%d_orf_from_0" % n, "lenedge", ("ATG" + neutral(rng, 8) + "TAA" + "CCCCC")[:n], "length %d, 30-nt ORF at 0" % n)

    return cases


def reverse_complements(cases):
    out = []
    for c in cases:
        if c.cat in ("lenedge",):
            continue
        out.append(Case(c.id + "_rc", c.cat + "_rc", rc(c.seq), "reverse complement of " + c.id))
    return out


# ---------------------------------------------------------------------------------------------------------
# random.fa
# ---------------------------------------------------------------------------------------------------------
def gene_dense(rng, n_orfs, max_body, allow_overlap=True):
    parts = []
    for _ in range(n_orfs):
        gap = rng.between(0, 40)
        parts.append(noise(rng, gap, 0.5))
        body = rng.between(4, max_body)
        starts = rng.choice(["ATG"] * 4 + ["TTG", "CTG", "GTG", "ATT", "ATC", "ATA", "TTA"])
        sense = [c for c in ALL_CODONS if c not in STOPS_ALL]
        mid = "".join(rng.choice(sense) for _ in range(body))
        stop = rng.choice(["TAA", "TAG", "TGA", "AGA", "AGG"])
        o = starts + mid + stop
        parts.append(o if rng.below(2) == 0 else rc(o))
    s = "".join(parts)
    if allow_overlap and len(s) > 90 and rng.below(2):
        # splice a second copy shifted by 1-2 nt over the first half to force cross-frame overlaps
        k = rng.between(1, 2)
        a = rng.below(len(s) // 2)
        seg = s[a:a + len(s) // 3]
        s = s[:a + k] + seg + s[a + k + len(seg):]
    return s


def mutate(rng, s, rate):
    out = list(s)
    for i in range(len(out)):
        if rng.rand() < rate:
            out[i] = rng.choice(BASES)
    return "".join(out)


def random_set(rng):
    cases = []
    for i in range(60):
        gc = (0.30, 0.40, 0.50, 0.60, 0.70)[i % 5]
        n = rng.between(60, 900)
        cases.append(Case("rnd_iid_%03d" % i, "random_iid", noise(rng, n, gc), "iid bases GC=%.2f length %d" % (gc, n)))
    for i in range(60):
        s = gene_dense(rng, rng.between(1, 8), rng.between(10, 120))
        cases.append(Case("rnd_dense_%03d" % i, "random_dense", s, "mosaic of implanted ORFs on both strands"))
    for i in range(40):
        base = gene_dense(rng, rng.between(2, 6), 80)
        s = mutate(rng, base, rng.choice([0.01, 0.02, 0.05]))
        cases.append(Case("rnd_mut_%03d" % i, "random_mut", s, "point-mutated gene-dense mosaic"))
    for i in range(10):
        n = rng.between(900, 2400)
        cases.append(Case("rnd_long_%03d" % i, "random_long", noise(rng, n, 0.5 + 0.05 * (i - 5)), "long iid sequence"))
    return cases


# ---------------------------------------------------------------------------------------------------------
# circular.fa
# ---------------------------------------------------------------------------------------------------------
def rotate(s, r):
    r %= len(s)
    return s[r:] + s[:r]


def circular_set(rng):
    cases = []
    # one ORF (ATG + body + stop) in a CCC background, rotated so the origin falls at several places inside / next to it
    for strand in ("+", "-"):
        for body in (12, 30):
            for pad_extra in (0, 1, 2):
                orf = "ATG" + neutral(rng, body) + "TAA"
                g = "CCC" * 5 + orf + "CCC" * 5 + "C" * pad_extra
                if strand == "-":
                    g = rc(g)
                a = (15 if strand == "+" else len(g) - 15 - len(orf))  # forward start of the span
                n = len(g)
                cuts = {
                    "start1": a + 1, "start2": a + 2, "start3": a + 3, "mid": a + len(orf) // 2,
                    "stop1": a + len(orf) - 2, "stop2": a + len(orf) - 1, "after": a + len(orf), "before": a,
                }
                for name, r in cuts.items():
                    add_id = "circ_%s_b%d_p%d_%s" % ("f" if strand == "+" else "r", body, pad_extra, name)
                    cases.append(Case(add_id, "circular", rotate(g, r),
                                      "%s-strand ORF of %d nt, genome %d nt, origin cut at %s (rotation %d)" % ("forward" if strand == "+" else "reverse", len(orf), n, name, r)))
    # nested in frame across the origin, and cross-frame nesting across the origin
    for k in range(6):
        orf = "ATG" + neutral(rng, 4) + "ATG" + neutral(rng, 14) + "TAG"
        g = "CCC" * 4 + orf + "CCC" * 4
        cases.append(Case("circ_nest_%d" % k, "circular", rotate(g, rng.between(10, len(g) - 10)), "two in-frame ATGs, origin somewhere inside"))
    for k in range(8):
        spec = [("+", 12, 24, "ATG", "TAA"), ("+", 12 + 3 * 3 + 1 + (k % 2), 7, "ATG", "TAA")] if k < 4 else [("+", 12, 24, "ATG", "TAA"), ("-", 12 + 14 + (k % 3), 8, "ATG", "TAA")]
        seq, spans = build_with_orfs(rng, 12 + 78 + 12, spec, tries=12000)
        cases.append(Case("circ_cross_%d" % k, "circular", rotate(seq, rng.between(20, len(seq) - 20)), "overlapping ORFs of two frames/strands, origin cut inside; designed spans before rotation %s" % spans))
    for k in range(6):
        s = gene_dense(rng, 2, 25)
        cases.append(Case("circ_dense_%d" % k, "circular", s, "short gene-dense circular genome %d nt" % len(s)))
    for k in range(6):
        s = gene_dense(rng, 3, 60)
        cases.append(Case("circ_dense_long_%d" % k, "circular", rotate(s, rng.below(len(s))), "gene-dense, random rotation"))
    return cases


# ---------------------------------------------------------------------------------------------------------
# region.fa : 60 sequences of exactly 300 nt
# ---------------------------------------------------------------------------------------------------------
def region_set(rng):
    cases = []
    # ORFs that start/end near the grid positions 50, 100, 150, 200, 250 (1-based), forward and reverse
    anchors = [(0, 40), (30, 30), (60, 20), (90, 33), (100, 30), (120, 25), (145, 20), (150, 31), (175, 22), (190, 20), (210, 28), (240, 18)]
    for i in range(30):
        specs = []
        k = rng.between(2, 5)
        for _ in range(k):
            a0, nb = rng.choice(anchors)
            a = a0 + rng.between(0, 2)
            nb = max(4, min(nb, (300 - a) // 3 - 3))
            specs.append((rng.choice("+-"), a, nb, "ATG", "TAA"))
        try:
            seq, spans = build_with_orfs(rng, 300, specs, tries=3000)
        except RuntimeError:
            seq, spans = build_with_orfs(rng, 300, specs[:1], tries=3000)
        cases.append(Case("reg_design_%03d" % i, "region", seq, "designed ORFs, forward spans %s" % spans))
    for i in range(15):
        cases.append(Case("reg_dense_%03d" % i, "region", (gene_dense(rng, 6, 40) + noise(rng, 300, 0.5))[:300], "gene-dense, cut to 300 nt"))
    for i in range(15):
        cases.append(Case("reg_iid_%03d" % i, "region", noise(rng, 300, 0.35 + 0.02 * i), "iid 300 nt"))
    for c in cases:
        assert len(c.seq) == 300, c.id
    return cases


# ---------------------------------------------------------------------------------------------------------
# edge/*.fa  (raw FASTA text, one record per file) and circular_edge/*.fa
# ---------------------------------------------------------------------------------------------------------
def edge_files(rng):
    orf30 = "ATG" + neutral(rng, 8) + "TAA"
    f = {}
    f["e01_len0.fa"] = ">e01_len0\n"
    f["e02_len1.fa"] = ">e02_len1\nA\n"
    f["e03_len2.fa"] = ">e03_len2\nAT\n"
    f["e04_len3_atg.fa"] = ">e04_len3_atg\nATG\n"
    f["e05_len6_atgtaa.fa"] = ">e05_len6_atgtaa\nATGTAA\n"
    f["e06_len29.fa"] = ">e06_len29\n" + (orf30)[:29] + "\n"
    f["e07_len30_orf.fa"] = ">e07_len30_orf\n" + orf30 + "\n"
    f["e08_lowercase.fa"] = ">e08_lowercase\n" + ("ccc" * 4 + "atg" + neutral(rng, 30).lower() + "taa" + "ccc" * 4) + "\n"
    f["e09_mixedcase.fa"] = ">e09_mixedcase\n" + ("CCC" * 4 + "ATG" + neutral(rng, 10) + "taa" + "aTg" + neutral(rng, 10).lower() + "TGA" + "CCC") + "\n"
    f["e10_with_N.fa"] = ">e10_with_N\n" + ("CCC" * 4 + "ATG" + neutral(rng, 10) + "NNN" + neutral(rng, 10) + "TAA" + "CCC" * 4) + "\n"
    f["e11_N_in_start.fa"] = ">e11_N_in_start\n" + ("CCC" * 4 + "NTG" + neutral(rng, 20) + "TAA" + "CCC" * 4) + "\n"
    f["e12_N_in_stop.fa"] = ">e12_N_in_stop\n" + ("CCC" * 4 + "ATG" + neutral(rng, 20) + "TAN" + neutral(rng, 6) + "TAA" + "CCC" * 4) + "\n"
    f["e13_iupac_R_Y.fa"] = ">e13_iupac_R_Y\n" + ("CCC" * 4 + "ATG" + neutral(rng, 10) + "TAR" + neutral(rng, 6) + "YAA" + "CCC" * 4) + "\n"
    f["e14_U_rna.fa"] = ">e14_U_rna\n" + ("CCC" * 4 + "AUG" + neutral(rng, 20).replace("T", "U") + "UAA" + "CCC" * 4) + "\n"
    f["e15_gap_dash.fa"] = ">e15_gap_dash\n" + ("CCC" * 4 + "ATG" + neutral(rng, 10) + "---" + neutral(rng, 10) + "TAA" + "CCC" * 4) + "\n"
    s = "CCC" * 4 + "ATG" + neutral(rng, 40) + "TAA" + "CCC" * 4
    f["e16_crlf.fa"] = ">e16_crlf\r\n" + s[:60] + "\r\n" + s[60:] + "\r\n"
    f["e17_wrapped_60_blank_lines.fa"] = ">e17_wrapped\n\n" + s[:60] + "\n" + s[60:] + "\n\n"
    f["e18_long_header.fa"] = ">e18_long_header " + "x" * 300 + " |gi|123|\n" + s + "\n"
    f["e19_header_with_spaces.fa"] = ">e19 header with spaces and | pipes | lcl|x\n" + s + "\n"
    f["e20_no_trailing_newline.fa"] = ">e20_no_trailing_newline\n" + s
    f["e21_two_records_same_id.fa"] = ">dup\n" + s + "\n>dup\n" + rc(s) + "\n"
    f["e22_protein_input.fa"] = ">e22_protein_input\nMKVLAAGIVGLLLAQPSTWQ\n"
    f["e23_all_N.fa"] = ">e23_all_N\n" + "N" * 90 + "\n"
    f["e24_plus_empty_record_between.fa"] = ">e24_a\n" + s + "\n>e24_empty\n>e24_b\n" + s + "\n"
    return f


def circular_edge_files(rng):
    f = {}
    f["ce1_ccc_only.fa"] = ">ce1_ccc_only\n" + "CCC" * 30 + "\n"
    f["ce2_atg_nostop_frame.fa"] = ">ce2_atg_nostop_frame\n" + "ATG" + "CCC" * 29 + "\n"
    f["ce3_atg_nostop_len_not_mult3.fa"] = ">ce3_atg_nostop_len_not_mult3\n" + "ATG" + "CCC" * 29 + "C" + "\n"
    f["ce4_stop_only_one.fa"] = ">ce4_stop_only_one\n" + "CCC" * 10 + "TAA" + "CCC" * 10 + "\n"
    f["ce5_atg_only_one_stop_elsewhere.fa"] = ">ce5_atg_only_one_stop_elsewhere\n" + "ATG" + "CCC" * 10 + "TAA" + "CCC" * 10 + "\n"
    f["ce6_no_stop_any_strand_atg_both.fa"] = ">ce6_no_stop_any_strand_atg_both\n" + "ATG" + "CCC" * 10 + "CAT" + "CCC" * 10 + "\n"
    return f


# ---------------------------------------------------------------------------------------------------------
def fasta(records, width=70):
    out = []
    for cid, seq in records:
        out.append(">" + cid + "\n")
        for i in range(0, len(seq), width):
            out.append(seq[i:i + width] + "\n")
    return "".join(out)


def build():
    """Returns {relative path: text}."""
    files = {}
    manifest = ["id\tfile\tcategory\tlength\tnote"]

    d = designed(Rng(20261118))
    d += reverse_complements(d)
    r = random_set(Rng(20261119))
    c = circular_set(Rng(20261120))
    g = region_set(Rng(20261121))
    for fname, cases in (("designed.fa", d), ("random.fa", r), ("circular.fa", c), ("region.fa", g)):
        ids = [x.id for x in cases]
        assert len(ids) == len(set(ids)), "duplicate ids in " + fname
        files["inputs/" + fname] = fasta([(x.id, x.seq) for x in cases])
        for x in cases:
            manifest.append("%s\t%s\t%s\t%d\t%s" % (x.id, fname, x.cat, len(x.seq), x.note.replace("\t", " ")))
    for name, text in sorted(edge_files(Rng(20261122)).items()):
        files["inputs/edge/" + name] = text
        manifest.append("%s\tedge/%s\tedge\t-\tmalformed / edge-case input" % (name[:-3], name))
    for name, text in sorted(circular_edge_files(Rng(20261123)).items()):
        files["inputs/circular_edge/" + name] = text
        manifest.append("%s\tcircular_edge/%s\tcircular_edge\t-\tcircular genome without a terminating stop (run with a timeout)" % (name[:-3], name))
    files["inputs/manifest.tsv"] = "\n".join(manifest) + "\n"
    sums = ["%s  %s" % (hashlib.sha256(files[k].encode("utf-8", "surrogateescape")).hexdigest(), k[len("inputs/"):])
            for k in sorted(files) if k != "inputs/inputs.sha256"]
    files["inputs/inputs.sha256"] = "\n".join(sums) + "\n"
    return files


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    files = build()
    if "--check" in sys.argv:
        bad = 0
        for k, text in files.items():
            p = os.path.join(here, k)
            try:
                with open(p, "rb") as fh:
                    cur = fh.read()
            except OSError:
                print("MISSING", k)
                bad += 1
                continue
            if cur != text.encode("utf-8"):
                print("DIFFERS", k)
                bad += 1
        extra = set()
        for root, _, names in os.walk(os.path.join(here, "inputs")):
            for n in names:
                rel = os.path.relpath(os.path.join(root, n), here).replace(os.sep, "/")
                if rel not in files:
                    extra.add(rel)
        for e in sorted(extra):
            print("EXTRA", e)
            bad += 1
        print("inputs OK (%d files)" % len(files) if not bad else "inputs check FAILED (%d problems)" % bad)
        sys.exit(1 if bad else 0)
    for k, text in files.items():
        p = os.path.join(here, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "wb") as fh:
            fh.write(text.encode("utf-8"))
    n_seq = sum(t.count(">") for k, t in files.items() if k.endswith(".fa"))
    print("wrote %d files, %d FASTA records" % (len(files), n_seq))


if __name__ == "__main__":
    main()
