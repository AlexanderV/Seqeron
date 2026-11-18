#!/usr/bin/env python3
"""GENOMIC-ORF-001 oracle: independent six-frame ORF finders vs NCBI COrf port vs EMBOSS getorf 6.6.0.

Three independent implementations are computed for every (sequence, options) row and asserted equal where
their definitions coincide; the TSV then locks the library against them.

  1. `brute`  - definition-based Python finder written from the Rosalind ORF / NCBI ORFfinder definitions:
                Biopython 1.88 `Seq.translate(table=t)` decides STOP ('*'); starts from Biopython
                `CodonTable.start_codons` (alt) / 'ATG' / any sense codon; reverse strand via Biopython
                `reverse_complement`. Every START is paired with the first in-frame STOP (nested), optionally
                reduced to the most upstream START per STOP (longest), optionally with open-ended ORFs.
  2. `ncbi`   - line-by-line port of NCBI C++ Toolkit `COrf::FindOrfs` (src/algo/sequence/orf.cpp,
                `FindForwardOrfs`/`s_FindOrfs`, raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public master),
                the engine behind NCBI ORFfinder: `longest_orfs = False` == nested, `True` == longest per stop;
                allowable_starts ["ATG"] / table starts / ["STOP"] (any sense codon = stop-to-stop); its min length
                excludes the stop codon; open-ended ORFs are always reported (3'-partial). Not applicable to
                nested + any-sense (COrf has no such mode).
  3. `getorf` - EMBOSS getorf 6.6.0 `-find 3 -table T -minsize M` for the longest-ORF-per-stop, alt-start,
                open-ended, stop-excluded-length configuration (the one-codon open ORF that getorf never writes,
                a START that is the last complete codon of its frame, is dropped from the comparison).

Output TSV (embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/orf_oracle.tsv), tab separated:
  table  startMode(atg|alt|any)  nested(0|1)  open(0|1)  excludeStop(0|1)  minLength  bothStrands(0|1)  sequence  orfs
where orfs = ';'-joined `strand,frame,start,end,terminated` in library enumeration order (forward frames 1..3, then
reverse frames 1..3, ascending scanned-strand start within a frame); start/end are 0-based INCLUSIVE
forward-strand coordinates of the span incl. the stop codon, frame is the frame within the scanned strand,
terminated is 1/0. Empty = no ORF; '-' placeholder used for empty.
Usage: python3 generate_orf_oracle.py > orf_oracle.tsv      (needs biopython and the `getorf` binary, Debian emboss)
"""
import os
import random
import subprocess
import sys
import tempfile

from Bio.Data import CodonTable
from Bio.Seq import Seq

rnd = random.Random(20261118)
DUAL = {27, 28, 31}
TABLES = [t for t in sorted(CodonTable.unambiguous_dna_by_id) if t not in DUAL]
BASES = "ACGT"
# EMBOSS 6.6.0 EGC files whose start/stop sets equal NCBI gc.prt 4.6 (= Biopython); EGC.3 and EGC.13 carry older start sets
GETORF_TABLES = [1, 2, 4, 5, 6, 9, 10, 11, 12, 14, 15, 16, 21, 22, 23]


def rc(s):
    return str(Seq(s).reverse_complement())


def codon_sets(t):
    ct = CodonTable.unambiguous_dna_by_id[t]
    return set(ct.start_codons), set(ct.stop_codons)


def is_stop(codon, t):
    return str(Seq(codon).translate(table=t)) == "*"


# --------------------------------------------------------------------------------------------------
# 1. brute: definition-based
# --------------------------------------------------------------------------------------------------
def brute_strand(s, t, start_mode, nested, include_open):
    starts_alt, _ = codon_sets(t)
    n = len(s)
    out = []  # (frame0, start, end_excl, terminated)
    for f in range(3):
        pos = list(range(f, n - 2, 3))
        stop_at = {p: is_stop(s[p:p + 3], t) for p in pos}

        def is_start(p):
            c = s[p:p + 3]
            if stop_at[p]:
                return False
            if start_mode == "atg":
                return c == "ATG"
            if start_mode == "alt":
                return c in starts_alt
            return True  # any sense codon

        groups = {}  # end_excl (or None for open) -> [starts]
        for idx, p in enumerate(pos):
            if not is_start(p):
                continue
            q = next((x for x in pos[idx + 1:] if stop_at[x]), None)
            if q is not None:
                groups.setdefault(q + 3, []).append(p)
            elif include_open:
                last_end = pos[-1] + 3 if pos else f
                groups.setdefault(None, []).append(p)
        for key in sorted(k for k in groups if k is not None):
            ps = sorted(groups[key])
            for p in (ps if nested else ps[:1]):
                out.append((f, p, key, True))
        if None in groups:
            ps = sorted(groups[None])
            last_end = pos[-1] + 3
            for p in (ps if nested else ps[:1]):
                out.append((f, p, last_end, False))
    # library order: frame, then ascending start
    out.sort(key=lambda r: (r[0], r[1]))
    return out


def brute(seq, t, start_mode, nested, include_open, both):
    n = len(seq)
    res = []
    for f, p, e, term in brute_strand(seq, t, start_mode, nested, include_open):
        res.append(("+", f + 1, p, e - 1, term))
    if both:
        for f, p, e, term in brute_strand(rc(seq), t, start_mode, nested, include_open):
            res.append(("-", f + 1, n - e, n - 1 - p, term))
    return res


# --------------------------------------------------------------------------------------------------
# 2. ncbi: port of COrf::FindForwardOrfs / s_FindOrfs (ACGT only, so the gap/N machinery is not exercised)
# --------------------------------------------------------------------------------------------------
def ncbi_forward(seq, min_bp, t, allowable_starts, longest_orfs, stop_to_stop):
    n = len(seq)
    stops = [[], [], []]
    for i in range(n - 2):
        if is_stop(seq[i:i + 3], t):
            stops[i % 3].append(i)
    ranges = []  # (from, to_inclusive_incl_stop_or_open, to_fuzz)

    def add(frm, to, to_fuzz):
        if not to_fuzz:
            to += 3
        ranges.append((frm, to, to_fuzz))

    for frame in range(3):
        st = stops[frame] + [n, n]
        frm = frame
        for i in range(len(st) - 1):
            frm0 = frm
            stop = st[i]
            gap_after = stop >= n
            if stop >= min_bp + frm:
                to = ((stop - frm) // 3) * 3 + frm - 1
                if to + 1 >= min_bp + frm:
                    starts = None
                    if allowable_starts:
                        inframe_to_open = to + 1
                        starts = {inframe_to_open}
                        pos = inframe_to_open - 3
                        while pos >= frm:
                            if seq[pos:pos + 3] in allowable_starts:
                                starts.add(pos)
                            pos -= 3
                        frm = min(starts)
                    if to + 1 >= min_bp + frm:
                        if frm != frm0 and stop_to_stop:
                            add(frm0, to, gap_after)
                        if not (stop_to_stop and frm != frm0 and longest_orfs):
                            add(frm, to, gap_after)
                        if not longest_orfs and allowable_starts:
                            starts.discard(min(starts))
                            for s in sorted(starts):
                                if to + 1 < min_bp + s:
                                    break
                                add(s, to, gap_after)
                    elif stop_to_stop:
                        frm = frm0
                        add(frm, to, gap_after)
            if gap_after:
                break
            frm = stop + 3
    return ranges


def ncbi(seq, min_bp, t, start_mode, longest, both):
    n = len(seq)
    starts_alt, _ = codon_sets(t)
    if start_mode == "atg":
        allowed, s2s = ["ATG"], False
    elif start_mode == "alt":
        allowed, s2s = sorted(starts_alt), False
    else:
        allowed, s2s = [], True  # ["STOP"] alone
    min_bp = max(min_bp, 3)
    res = []

    def frame_order(rs, comp):
        # library order: by frame (scanned-strand coordinate % 3), ascending start within frame
        return sorted(rs, key=lambda r: (r[0] % 3, r[0]))

    for frm, to, fuzz in frame_order(ncbi_forward(seq, min_bp, t, set(allowed), longest, s2s), False):
        res.append(("+", frm % 3 + 1, frm, to, not fuzz))
    if both:
        comp = rc(seq)
        for frm, to, fuzz in frame_order(ncbi_forward(comp, min_bp, t, set(allowed), longest, s2s), True):
            e = to
            res.append(("-", frm % 3 + 1, n - 1 - e, n - 1 - frm, not fuzz))
    return res


# --------------------------------------------------------------------------------------------------
# 3. getorf
# --------------------------------------------------------------------------------------------------
def getorf_orfs(seq, t, minsize, both):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "s.fa")
        open(fa, "w").write(">s\n" + seq + "\n")
        cmd = ["getorf", "-sequence", fa, "-outseq", "stdout", "-find", "3", "-table", str(t),
               "-minsize", str(minsize), "-auto"]
        if not both:
            cmd.append("-noreverse")
        r = subprocess.run(cmd, capture_output=True, text=True)
        txt = r.stdout
    recs = []
    for line in txt.splitlines():
        if line.startswith(">"):
            hdr = line[line.index("[") + 1:line.index("]")]
            s, e = [int(x.strip()) for x in hdr.split("-")]
            recs.append((s, e, "REVERSE" in line))
    out = set()
    for s, e, rev in recs:
        if not rev:
            out.add(("+", s - 1, e - 1))  # excl-stop forward span (0-based inclusive)
        else:
            out.add(("-", e - 1, s - 1))
    return out


def to_getorf_form(orfs):
    """library ORFs -> (strand, lo, hi) excl-stop spans, dropping one-codon open ORFs (getorf never writes them)."""
    out = set()
    for strand, frame, a, b, term in orfs:
        if not term and (b - a + 1) == 3:
            continue
        if strand == "+":
            out.add((strand, a, b - 3 if term else b))
        else:
            out.add((strand, a + 3 if term else a, b))
    return out


# --------------------------------------------------------------------------------------------------
def gen_seq():
    kind = rnd.random()
    if kind < 0.15:
        n = rnd.randint(0, 24)
        return "".join(rnd.choice(BASES) for _ in range(n))
    parts = []
    target = rnd.choice([40, 90, 150, 240, 360])
    bias = rnd.choice(["ACGT", "ACGT", "AT", "ATTG", "ACCG"])
    codons = ["ATG", "TTG", "CTG", "GTG", "ATA", "ATT", "ATC", "TAA", "TAG", "TGA", "AGA", "AGG"]
    while sum(map(len, parts)) < target:
        r = rnd.random()
        if r < 0.45:
            parts.append(rnd.choice(codons))
        elif r < 0.55:
            parts.append(rnd.choice(["ATG", "TGA"]) * rnd.randint(1, 4))
        else:
            parts.append("".join(rnd.choice(bias) for _ in range(rnd.randint(1, 9))))
    return "".join(parts)


ROSALIND = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"


def fmt(orfs):
    return ";".join(f"{s},{f},{a},{b},{1 if t else 0}" for s, f, a, b, t in orfs) or "-"


def main():
    rows = []
    seqs = [ROSALIND, "", "ATG", "ATGTAA", "ATGATGATGATGTGA", "TTAGGGGGGCAT", "ATGGGGATGCCCTAA",
            "ATG" * 12, "TAA" * 5, "A" * 20, "ATGAAAAAATAA", "CATGCCCTAA"]
    seqs += [gen_seq() for _ in range(330)]
    n_getorf = n_ncbi = 0
    for idx, seq in enumerate(seqs):
        combos = []
        if idx < 12:
            combos.append((1, "atg", 1, 0, 0, 0, 1))  # defaults: ATG only, nested, no open, incl. stop
            combos.append((1, "atg", 1, 0, 0, 1, 1))
            combos.append((1, "atg", 0, 1, 1, 3, 1))
        for _ in range(4):
            combos.append((rnd.choice(TABLES), rnd.choice(["atg", "alt", "any"]), rnd.randint(0, 1),
                           rnd.randint(0, 1), rnd.randint(0, 1), rnd.choice([0, 3, 6, 9, 12, 18, 30, 60]),
                           rnd.choice([1, 1, 1, 0])))
        # one getorf-comparable combination per sequence
        combos.append((rnd.choice(GETORF_TABLES), "alt", 0, 1, 1, rnd.choice([3, 6, 9, 12, 18, 30]), rnd.choice([1, 1, 0])))
        for t, mode, nested, op, excl, minlen, both in combos:
            orfs = brute(seq, t, mode, bool(nested), bool(op), bool(both))
            # apply the library's length rule: measured = length - 3 if (excl and terminated) else length
            kept = []
            for s_, f_, a_, b_, term in orfs:
                length = b_ - a_ + 1
                measured = length - 3 if (excl and term) else length
                if measured >= minlen:
                    kept.append((s_, f_, a_, b_, term))
            # NCBI cross-check: min excludes stop for terminated, as-is for open; run COrf with the post-filter
            # picking (floor 3) and also with its own min when the convention coincides.
            if not (mode == "any" and nested):
                nc = ncbi(seq, 3, t, mode, longest=not nested, both=bool(both))
                nc = [o for o in nc if op or o[4]]
                nc_kept = [o for o in nc if ((o[3] - o[2] + 1) - (3 if (excl and o[4]) else 0)) >= minlen]
                key = lambda o: (o[0], o[1], o[2])
                assert sorted(nc_kept, key=key) == sorted(kept, key=key), (seq, t, mode, nested, op, excl, minlen, both,
                                                                         nc_kept, kept)
                if excl and minlen >= 3:
                    # COrf's own min_length_bp (stop excluded); open ORFs are always reported by COrf
                    nc2 = ncbi(seq, minlen, t, mode, longest=not nested, both=bool(both))
                    nc2 = [o for o in nc2 if op or o[4]]
                    assert sorted(nc2, key=key) == sorted(kept, key=key), ("own-min", seq, t, mode, nested, op, minlen)
                n_ncbi += 1
            if mode == "alt" and not nested and op and excl and minlen >= 3 and seq and t in GETORF_TABLES:
                g = getorf_orfs(seq, t, minlen, bool(both))
                mine = to_getorf_form(kept)
                assert g == mine, ("getorf", seq, t, minlen, both, sorted(g), sorted(mine))
                n_getorf += 1
            rows.append((t, mode, nested, op, excl, minlen, both, seq, kept))
    for t, mode, nested, op, excl, minlen, both, seq, kept in rows:
        print("\t".join(map(str, (t, mode, nested, op, excl, minlen, both, seq if seq else "-", fmt(kept)))))
    print(f"rows={len(rows)} ncbi-checked={n_ncbi} getorf-checked={n_getorf}", file=sys.stderr)


if __name__ == "__main__":
    main()
