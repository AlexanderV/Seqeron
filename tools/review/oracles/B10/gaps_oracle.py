"""Independent oracle for ASSEMBLY-GAPS-001 (GenomeAssemblyAnalyzer FindGaps / AnalyzeGapDistribution /
AnalyzeScaffolds / ExtractContigs).

References transcribed (sources fetched 2026-11-25):
- QUAST quast_libs/qutils.py `split_by_ns` (verbatim below) and `correct_fasta` (`len(seq) >= min_contig`
  filter applied to the broken-scaffold file); qconfig.Ns_break_threshold = 10.
  https://raw.githubusercontent.com/ablab/quast/master/quast_libs/qutils.py
- seqtk `stk_gap` (seqtk.c): a gap is a maximal run of non-ACGT characters of length >= -l, printed as
  0-based half-open BED [i-l, i).  https://raw.githubusercontent.com/lh3/seqtk/master/seqtk.c
- numpy.median (sample median).

Seqeron conventions checked against them: only N/n is a gap symbol (QUAST `correct_seq` upper-cases and
maps IUPAC codes to N; seqtk treats every non-ACGT as gap) -> the random alphabet below has no IUPAC
codes when comparing to seqtk; QUAST is fed `seq.upper()` and the original-case slice is compared.

Usage: python3 -I gaps_oracle.py OUT.json [n_cases] [seed]
Writes cases {seq, minGap, minContig, gaps:[[start,endInclusive,len]], contigs:[[start,endExclusive]],
extracted:[seq...], median}.
"""
import json
import random
import sys

import numpy as np


# --- QUAST qutils.split_by_ns, verbatim (ablab/quast master) -------------------------------------
def split_by_ns(seq, name, splitted_fasta, Ns_break_threshold=1, min_contig=1, total_contigs=0):
    cur_contig_start = 0
    cumul_contig_length = 0
    while (cumul_contig_length < len(seq)) and (seq.find('N', cumul_contig_length) != -1):
        start = seq.find("N", cumul_contig_length)
        end = start + 1
        while (end != len(seq)) and (seq[end] == 'N'):
            end += 1

        cumul_contig_length = end + 1
        if end - start >= Ns_break_threshold:
            splitted_fasta.append(
                (name.split()[0] + "_" +
                 str(total_contigs + 1),
                 seq[cur_contig_start: start]))
            total_contigs += 1
            cur_contig_start = end
    if len(seq) - cur_contig_start >= min_contig:
        splitted_fasta.append(
            (name.split()[0] + "_" +
             str(total_contigs + 1),
             seq[cur_contig_start: ]))
        total_contigs += 1
    return total_contigs


def quast_broken_pieces_exact(seq, threshold):
    """Coordinates of the split_by_ns pieces: the verbatim loop above, recording (start, endExclusive)
    instead of the slice; asserted equal to the verbatim function's slices in main()."""
    up = seq.upper()
    cur, cumul, pieces = 0, 0, []
    while cumul < len(up) and up.find('N', cumul) != -1:
        start = up.find('N', cumul)
        end = start + 1
        while end != len(up) and up[end] == 'N':
            end += 1
        cumul = end + 1
        if end - start >= threshold:
            pieces.append((cur, start))
            cur = end
    pieces.append((cur, len(up)))
    return [(a, b) for a, b in pieces if b > a]


def seqtk_gap(seq, min_size):
    """seqtk stk_gap transcription: BED [i-l, i) for maximal non-ACGT runs with l >= min_size."""
    out, l = [], 0
    n = len(seq)
    for i in range(n + 1):
        is_base = i < n and seq[i] in "ACGTacgt"
        if i == n or is_base:
            if l > 0 and l >= min_size:
                out.append((i - l, i))
            l = 0
        else:
            l += 1
    return out


def main():
    out_path = sys.argv[1]
    n_cases = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
    cases = []
    for _ in range(n_cases):
        length = rng.randint(0, 60)
        seq, i = [], 0
        while i < length:
            if rng.random() < 0.3:
                run = rng.choice([1, 1, 2, 3, 5, 9, 10, 11, 15, 25])
                seq.extend(rng.choice("Nn") for _ in range(run))
                i += run
            else:
                seq.append(rng.choice("ACGTacgt"))
                i += 1
        seq = "".join(seq)
        min_gap = rng.choice([1, 2, 3, 5, 10, 11, 20])
        min_contig = rng.choice([0, 1, 2, 3, 5, 8, 20])
        gaps = [(a, b - 1, b - a) for a, b in seqtk_gap(seq, min_gap)]
        contigs = quast_broken_pieces_exact(seq, min_gap)
        # cross-check the coordinate recovery against the verbatim function's piece strings
        verbatim = []
        split_by_ns(seq.upper(), "s", verbatim, min_gap, 1)
        assert [p for _, p in verbatim if p] == [seq.upper()[a:b] for a, b in contigs], seq
        extracted = [seq[a:b] for a, b in contigs if b - a >= max(min_contig, 1)]
        lens = [g[2] for g in gaps]
        cases.append({
            "seq": seq, "minGap": min_gap, "minContig": min_contig,
            "gaps": gaps, "contigs": contigs, "extracted": extracted,
            "median": float(np.median(lens)) if lens else 0.0,
            "mean": float(np.mean(lens)) if lens else 0.0,
        })
    with open(out_path, "w") as f:
        json.dump(cases, f)
    print(f"{len(cases)} cases written; "
          f"{sum(1 for c in cases if any(0 < g < c['minGap'] for g in [r[1]-r[0] for r in seqtk_gap(c['seq'],1)]))} "
          "with an N-run shorter than minGap")


if __name__ == "__main__":
    main()
