# Evidence Artifact: REP-DIRECT-001

**Test Unit ID:** REP-DIRECT-001
**Algorithm:** Exact direct-repeat detection reported as maximal repeated pairs (left- and right-maximal), forward strand, A/C/G/T only
**Date Collected:** 2026-09-30 (campaign 2026-09, batch B04, finding F11)

---

## Sources opened

### MUMmer `repeat-match` — source code (`src/tigr/repeat-match.cc`, A. Delcher, rev. 2002-12-12)

**URL:** https://raw.githubusercontent.com/mummer4/mummer/master/src/tigr/repeat-match.cc (fetched, 46 kB) +
`src/tigr/tigrinc.cc`, `include/mummer/tigrinc.hh` (fetched); compiled locally with `g++ -O2` and run.
**Authority rank:** 2 (original tool source code; MUMmer 3/4, Kurtz et al. 2004 Genome Biol 5:R12)

Key points (read from the code):

1. File header: "This program identifies maximal exact repeat regions (longer than MIN_MATCH_LEN threshold) in
   input genome." Usage: "Find all maximal exact matches"; `-f` "Forward strand only, don't use reverse
   complement"; `-n #` minimum exact match length.
2. Exhaustive mode (`-E`, the reference definition in the code): for `i < j`, `if (Data[i-1] == Data[j-1]) continue;`
   (left-maximal; `Data[0]` is a unique start sentinel) and `match = Longest_Prefix_Match(Data+i, Data+j)`
   (right-maximal = full common-prefix length); report `(i, j, match)` when `match >= Min_Match_Len`.
   Overlapping copies are reported (no spacing constraint).
3. Suffix-tree mode `List_Matches`: skips pairs with `Data[i-1] == Data[j-1] || Data[i+n] == Data[j+n]` — the
   same left/right-maximality test. Output is 1-based `Start1 Start2 Length`.
4. `Read_String` lower-cases the input (case-insensitive); non-alphabetic characters become `x`.

### MUMmer `mummer` — `src/essaMEM/mummer.cpp`

**URL:** https://raw.githubusercontent.com/mummer4/mummer/master/src/essaMEM/mummer.cpp (fetched).
Option `-n`: "match only the characters a, c, g, or t" (`nucleotides_only`). Adopted as the non-ACGT convention:
N / IUPAC / gap symbols never match (consistent with REP-STR-001 / REP-INV-001 ACGT-only rules in this batch).

### Maximal pair definition (Gusfield 1997; Abouelhoda, Kurtz & Ohlebusch 2004; Brodal et al. 2000)

WebSearch snippets (full PDFs blocked by the proxy): "A repeated pair ((i1,j1),(i2,j2)) is called left maximal if
S[i1−1] ≠ S[i2−1] and right maximal if S[j1+1] ≠ S[j2+1]. A repeated pair is called maximal if it is left and right
maximal." "The algorithm of Gusfield computes maximal repeated pairs of a sequence S of length n in O(|Σ|n + z)
time, where z is the number of maximal repeated pairs." (Abouelhoda et al. 2004, J Discrete Algorithms 2:53–86;
Brodal, Lyngsø, Pedersen & Stoye 2000.) REPuter (Kurtz & Schleiermacher 1999, Bioinformatics 15:426) — title and
abstract: "fast computation of maximal repeats in complete genomes" (forward and palindromic modes).

Wikipedia "Direct repeat" / "Repeated sequence (DNA)" (definitions: same-orientation copies, spacer ≥ 0) remain the
biological definition; they do not define a reporting convention.

---

## Convention adopted

A pair of 0-based positions `i < j` with length `L` is reported iff

- `S[i..i+L) = S[j..j+L)` over A/C/G/T (case-insensitive; any other symbol never matches);
- right-maximal: `L` is the full common-prefix length of the suffixes at `i` and `j`;
- left-maximal: `i = 0` or `S[i−1] ≠ S[j−1]` (a non-ACGT left neighbour never matches);
- `minLength ≤ L ≤ maxLength` — a longer maximal repeat is **not** reported (not truncated);
- `Spacing = j − i − L ≥ minSpacing` (negative `minSpacing` admits overlapping copies; `int.MinValue` returns every
  maximal pair = `repeat-match -f`).

Output sorted by (FirstPosition, SecondPosition); a position pair appears at most once.

---

## Reference cross-checks (all 0 mismatches)

| Comparison | Cases | Result |
|---|---|---|
| Independent Python brute force of the definition vs compiled `repeat-match -f -n L` (random 1–120 bp, 5 alphabets, planted repeats, L 1–8) | 3 000 (2 473 non-empty) | 0 mismatches |
| C# `FindDirectRepeats` vs Python brute force (0–150 bp, N / IUPAC / `-` / lowercase, random min/max/spacing incl. ±int limits) | 8 000 (4 844 non-empty, 373 960 pairs) | 0 mismatches |
| C# (maxLength = int.MaxValue, minSpacing = int.MinValue) vs `repeat-match -f` (1–200 bp) | 2 000 (213 882 pairs) | 0 mismatches |
| same, 1–5 000 bp (low-complexity alphabets) | 200 (12 788 609 pairs) | 0 mismatches |
| same, 50 kb / 200 kb / 1 Mb random + 50 planted repeats (L = 10/14/16) | 3 (959 / 112 / 136 pairs) | identical; C# 0.2 / 0.3 / 1.1 s vs repeat-match 0.0 / 0.1 / 0.9 s |

Worked values (repeat-match -f, converted to 0-based):

| Sequence | -n | Maximal pairs (i, j, L) |
|---|---|---|
| ACGTACGTTTTTTTTTACGTACGT | 4 | (0,4,4) (0,16,8) (0,20,4) (3,15,5) (7,8,8) (7,9,7) (7,10,6) (7,11,5) (7,12,4) (15,19,5) |
| ACGTACGTACGT | 2 | (0,4,8) (0,8,4) |
| AAAAAATTTTAAAAAA | 5 | (0,1,5) (0,10,6) (0,11,5) (1,10,5) (10,11,5) |
| ACGTATTACGTATTACGTA | 4 | (0,7,12) (0,14,5) |
| A×40 | 2 | (0, j, 40 − j) for j = 1..38 |

---

## Defects found in the previous implementation (fixed)

1. Every `(i, j, len)` window for every `len ∈ [minLength, maxLength]` was reported: one repeat of length `L`
   yielded `O((L − minLength)²)` nested hits (`AAAAAATTTTAAAAAA`, 4–6, spacing ≥ 1 → 14 hits; maximal pairs: 5;
   snapshot `ACGTACGTTTTTTTTTACGTACGT` → 20 hits; maximal pairs: 4).
2. Negative `minSpacing` produced self-pairs `FirstPosition == SecondPosition` (`ACGTACGTACGT`, 4, 4, −4 → 15 hits
   including (0,0), (1,1), …).
3. N / non-ACGT runs reported as repeats (`N×20`, 5, 5, 1 → NNNNN pairs).
4. Cost `O(r · n · (m + k))` with a `Substring` + suffix-tree lookup per window and an `(i, j, len)` hash set.
