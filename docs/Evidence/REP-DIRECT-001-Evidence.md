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

---

## Enumeration variants (B04 completeness audit WP2, 2026-09-30)

### Sources opened
- mummer4 `src/tigr/repeat-match.cc` (already compiled for F11): without `-f` the reverse complement is added to the
  tree (`Data = % S $ revcomp(S) $`); `List_Matches` skips pairs with both leaves in revcomp and keeps a reverse pair
  only when `k ≥ i` (`k = Genome_Len − (j − String_Separator) − n + 2`), printing `L = i`, `R = k + n − 1` and `r`;
  `Verify_Match` checks `Data[a+t] = Complement(Data[b−t])`, i.e. Start2 is the **last** base of copy 2 (1-based).
  Both strings end in the same `$`, so `Add_String` merges a revcomp suffix equal to a suffix of S into one leaf
  ("Suffix can't appear twice"), and an exact whole-genome palindrome aborts ("Genome is exact palindrome").
- Vmatch 2.3.1 (Kurtz; ISC licence) — Ubuntu `vmatch` binary package + `vmatch_2.3.1+dfsg.orig.tar.xz` from
  archive.ubuntu.com: manual `src/doc/virtman.tex` Appendix A "Basic Notions" (palindromic match
  `u_i…u_{i+l−1} ≈ wcc(v_{j+r−1}…v_j)`, `i ≤ j` for self-comparison; k-mismatch match `d_H(x,y) ≤ k`;
  "A k-mismatch match is maximal if it is not contained in another k-mismatch match of the same kind";
  "A supermaximal repeat is a maximal repeat that never occurs as a substring of any other maximal repeat";
  wildcards "always lead to a pair of mismatching characters"); options `-p`, `-h`, `-allmax` ("compatibility with
  REPuter"), `-seedlength` (seed = max(⌊ℓ/(k+1)⌋, m)), `-supermax`; source `kurtz/extendHD.c` (REPuter maximum-error
  extension: per seed, left/right tables of the first k+1 mismatches, all splits). REPuter itself is not open source;
  Kurtz et al. 2001 (NAR 29:4633) and the Vmatch successor by the same author were used.
- Gusfield 1997 §7.12.1 / Theorem 7.12.4 (supermaximal repeat ⇔ internal node, all children leaves, left-diverse).

### Reference cross-checks (C# harness `scratchpad/rc/h` vs binaries / brute force; 0 mismatches unless stated)

| Comparison | Cases | Result |
|---|---|---|
| `FindReverseComplementRepeats(s, L, ∞, int.MinValue)` vs `repeat-match -n L` `r` lines (1–200 bp, 6 alphabets, planted revcomp copies) | 3 000 (318 713 pairs) | 0 mismatches after classifying repeat-match's shared-`$` leaf loss: 31 cases / 106 pairs present in C# and missing in repeat-match, every one a pair whose revcomp suffix equals a suffix of S; 14 self-reverse-complement inputs aborted by repeat-match |
| same, repeat-match run on `S + N` (unique sentinel removes both quirks; pairs touching the sentinel dropped) | 3 000 (279 490) + 100 cases of 1–5 kb (6 607 616 pairs) | 0 mismatches |
| vs `vmatch -p -l L` (`mkvtree -dna -pl -allout`) | 2 000 (186 444 pairs) | 0 mismatches |
| vs Python brute force of the definition (N / IUPAC / `-` / U / lowercase, random maxLength and minSpacing) | 3 000 (7 470 pairs) | 0 mismatches |
| `FindApproximateDirectRepeats` (both modes, random maxLength/minSpacing, N/IUPAC/`-`) vs Python brute force (per-diagonal; literal containment) | 3 000 (24 637 repeats) | 0 mismatches |
| `excludeContained = true` vs `vmatch -l m -h k -allmax` (k = 1–4; k = 0 vs `vmatch -l m`), ACGT/ACGTN/lowercase, 1–150 bp | 3 000 (432 579 repeats) | 0 mismatches |
| same, 1–1 500 bp | 40 (659 247 repeats) | 0 mismatches |
| Python check: `vmatch -h k -allmax` = per-diagonal maximal set minus repeats contained in a repeat on another diagonal | 300 (9 363) | 0 mismatches (hypothesis confirmed) |
| `FindSupermaximalRepeats` vs `vmatch -supermax -l m` (all position pairs) | 3 000 × ≤ 200 bp (24 767 pairs) + 60 × ≤ 3 kb (13 667) | 0 mismatches |
| vs Python brute force (maximal repeat strings not contained in another; N/IUPAC/U) | 3 000 (3 808 repeats) | 0 mismatches |
| `FindDirectRepeats` regression after the engine refactor vs `repeat-match -f` | 2 000 (223 155 pairs) | 0 mismatches |

Vmatch's **default** `-h k` output (no `-allmax`) keeps one E-value-best extension per seed and is not a set
definition (it can report a shorter-than-maximal window, e.g. (26,81,12,−1) where (26,81,13,−2) exists); the
library therefore reproduces the definitional set (`-allmax`) by default; since WP15 the default output is available as
`DegenerateRepeatReporting.BestPerSeed` (§WP15).

Timing (1 Mb random DNA, Release): direct (min 20) 0.63 s; reverse-complement (min 20) 1.2 s; k-mismatch
(min 30, k = 2, both modes) 0.76 s; supermaximal (min 15) 0.53 s.

### Worked values locked in tests
| Call | Result (0-based) | Reference |
|---|---|---|
| RC `AAAAAAAACGTTGCAACGTAAAA`, 3 | (6,6,6) (7,7,12) (15,15,4) | repeat-match `7 12r 6`, `8 19r 12`, `16 19r 4`; vmatch -p |
| RC `TTGCATGCAAAAAATTTTTTTGCATGCAA`, 4 | 14 pairs incl. (0,15,14) | repeat-match = vmatch -p |
| RC `GGATCCTTTTTTTGGATCCAAAAAAA`, 4 | (6,19,4) (6,19,5) (6,19,6) … | vmatch -p |
| k-mm `ACGTTGCAAGCTTACGGGGGGACGATGCAAGCATACGG`, 10, 2 | (0,21,17, 2 mm) | vmatch -h 2 |
| k-mm `AAAAAAAACGTTGCAACGTAAAA`, 5, 1, excludeContained | 7 repeats | vmatch -h 1 -allmax |
| k-mm `ACGTACGTACTTTTACGTNCGTAC`, 8, 1 | (0,4,8,1) (0,14,10,1) | vmatch -h 1 (N = wildcard) |
| supermax `CAGCAGCAGTTTCAGCAG`, 3 | CAGCAG @ 0,3,12 | vmatch -supermax |

## WP8 — degenerate repeats: k-differences (edit distance) and palindromic (B04 completeness audit, 2026-10-01)

### Sources opened
- Vmatch 2.3.1 source, Ubuntu `vmatch_2.3.1+dfsg.orig.tar.xz` (archive.ubuntu.com pool) + Debian patches
  `vmatch_2.3.1+dfsg-9.debian.tar.xz`: `src/doc/virtman.tex` (options `-e`, `-h`, `-p`, `-allmax` "for compatibility with
  REPuter", `-seedlength` = max(⌊ℓ/(k+1)⌋, m); App. A "Basic Notions": edit distance d_E, direct match `i < j`,
  palindromic match `i ≤ j`, k-mismatch / k-differences match, containment, maximality "not contained in another match
  of the same kind"); `kurtz/extendED.c` (`editextend`: left/right greedy fronts, every combination with
  `lookindex + (dist − lookindex) = dist ≤ k`, both instances ≥ ℓ, swap to `pos1 ≤ pos2`, `acceptmatch`); `kurtz/frontSEP.c`
  (Ukkonen/Myers furthest-reaching fronts; `evalentrybackward` stops a left extension that scans an exact run ≥ the seed
  length — "seed … detected while scanning"); `kurtz/mcontain.c` (`matchcontainer`); `Vmengine/fself.c` (per-seed then
  global container for `-allmax`); `Vmengine/extendgen.c`.
- Kurtz et al. 2001 NAR 29:4633 (REPuter degenerate repeats, edit distance); Ukkonen 1985; Myers 1986.

### Reference build
`apt install vmatch` (2.3.1+dfsg-9) = stock binary. Same release built from source (Debian patches, `Makedef-debian`,
genometools for prototype generation); `vmatch.x` reproduces the stock output byte for byte. One switch added: env
`VM_NOPRUNE` passes `reachlength = UINT_MAX` to `extendedleftSEP` (disables the left-extension seed shortcut only).
Brute force (`scratchpad/wp8/oracle.c`): DP over all (i, j, l, r), dominance closure for containment, `acceptmatch`
for direct edit matches, palindromic containment over both orientations, output `l, r ≥ ℓ`.

### Findings
1. `vmatch -e k -allmax` (direct) = definition + `acceptmatch` (direct edit only; without it the brute force reports
   trivial whole-sequence shifts in 297/300 cases) **except** for the seed shortcut: with `VM_NOPRUNE` 0 differences.
2. Palindromic: Vmatch runs `-p` self-comparisons as a query of `revcomp(S)` against `S`; containment is decided in that
   space (both orientations), output filtered to `i ≤ j`. A brute force restricting containers to `i ≤ j` differs on
   264 matches with `i = j`, `l ≠ r` (e.g. `(7,6,6,6,1)` kept); with both orientations 0 differences.
3. Hamming (`-h`, `-p -h`): stock Vmatch = definition (no shortcut in `extendHD.c`).

### Cross-checks (`FindDegenerateRepeats`, `minSpacing = int.MinValue`; harness `scratchpad/wp8/cmp8.py` + `xc8`)
| Mode | Cases | Repeats | vs brute force | vs Vmatch (`VM_NOPRUNE`) | vs stock Vmatch |
|---|---|---|---|---|---|
| `-e k` direct, 8–50 bp | 1 500 | 21 155 | 0 | 0 | 35 cases / 59 repeats differ (shortcut) |
| `-p -e k`, 8–50 bp | 1 500 | 17 740 | 0 | 0 | 23 / 39 (shortcut) |
| `-h k` direct (`-allmax`), 8–50 bp | 1 500 | 21 895 | 0 | 0 | 0 |
| `-p -h k`, 8–50 bp | 1 500 | 17 717 | 0 | 0 | 0 |
| `-e k` direct, 100–1 500 bp | 200 | 170 527 | — | 0 | 9 / 234 (shortcut) |
| `-p -e k`, 100–1 500 bp | 200 | 124 752 | — | 0 | 4 / 316 (shortcut) |
| `-h k`, 100–1 500 bp | 200 | 138 127 | — | 0 | 0 |
| `-p -h k`, 100–1 500 bp | 200 | 101 181 | — | 0 | 0 |
| `-e 2 -l 30`, direct + `-p`, 1 Mb random + 500 planted 100-bp copies (5 % noise) | 1 | 2 333 | — | 0 | 0 |

k = 1–3 (edit) / 1–4 (Hamming), ℓ = k+1 … k+1+max(8, n/6); inputs: uniform ACGT, AC-only, ACGT with 5 % N, motif
copies with substitutions/indels and spacers, single N. Timing 1 Mb (m 30, k 2): direct edit 0.83 s, palindromic edit
1.78 s, Hamming 0.88 / 1.69 s.

### Worked values locked in tests (`RepeatFinder_DegenerateRepeats_Tests`; tuples (i, j, l, r, d))
| Call | Result | Reference |
|---|---|---|
| `ACGTTGCATGCAAACGTAGCATGCAGGGTTTACGTTGCTTGCAAACG` -l 8 -e 1 | (0,13,12,12,1) (0,31,16,16,1) (5,18,8,8,1) (8,39,9,8,1) | stock = VM_NOPRUNE = brute force |
| `ATCTGGTGTACTCTGCCCACGACTATCGGTGTACTCTGC` -l 15 -e 3 | (0,21,16,18,3) (0,23,17,16,3) (0,24,18,15,3) | VM_NOPRUNE = brute force; stock (0,22,16,17,3) |
| `ACACACACACACACACAC` -l 4 -e 1 | (0,1,16,17,1) (0,2,17,16,1) | all three |
| `GATTACANGATTACA` -l 6 -e 1 | (0,7,7,8,1) (0,8,8,7,1) | all three |
| `GGACCATGAAGG` -p -l 5 -e 3 | 9 matches incl. (0,0,5,7,3) (0,0,7,5,3) (1,4,7,7,3) | VM_NOPRUNE = brute force; stock (1,4,7,6,3) |
| `TTGACCGTAACCCCCGTTACGGTCAACC` -p -l 8 -e 1 / -h 1 | 3 / 2 matches | all three |

## WP15 — Vmatch default output (best per seed) and stock compatibility (B04 audit, 2026-10-01)

### Sources opened
- Vmatch 2.3.1 source (same tree as §WP8):
  - `src/doc/virtman.tex`. Option `-allmax`: "Report all maximal matches … mainly used for compatibility with REPuter". The paragraph after Table "Exclude": "If the options -best, -allmax, and -complete are not used, then for each seed a best match, i.e. one with a minimum E-value is output … There is no limit on the number of matches reported".
  - `Vmengine/extendgen.c`. Without `-allmax`, `processallmax = NULL` and `processfinal` is called once per seed.
  - `kurtz/extendHD.c`:
    - `hammingextend`, with the lookleft/lookright tables and the boundary special cases;
    - `extendmismatchesleft`, whose `searchlength` shortcut returns h − 1 when the run between two mismatches is longer than the seed length;
    - the enumeration order: distance, then left errors.
  - `kurtz/extendED.c`:
    - `editextend`: the `maxleftextend`/`maxrightextend` early exit, the combination order (distance → lookindex → left diagonal → right diagonal), the swap to `pos1 ≤ pos2`, and `acceptmatch`;
    - best selection with `cmpmatches`.
  - `include/extcmp.c` (`cmpmatches`), in priority order:
    1. E-value via `incgetEvalue(evalues, 1.0, d, len)`;
    2. identity via `EVALIDENTITY`, `100(1 − |d|/len)`;
    3. length;
    4. a full tie returns 1 (replace).
  - `kurtz/evalues.c` (`incprecomputehammingEvalues`, `incgetEvalue`, the `averagequot` table, `SMALLESTEVALUE` 1e-300). `Vmatch/procmatch.c` sets the match probability to 1/(mapsize − 1) = 1/4 (checked: the printed E-value 9.42e-17 for (36, d 1) equals 36 · 0.25³⁵ · 0.75³ · n²/2).
  - `kurtz/frontSEP.c` and `kurtz/front.gen`:
    - fronts: `evalentryforward/backward`, `frontspecparms` band, `accessfront`;
    - minus infinity = −max(ulen, vlen);
    - the same-position rule `uptr == vptr`.
  - `kurtz/mcontain.c` (`matchcontainer`: a new match identical to a stored one is dropped, so the first seed's distance stays).
  - `Vmengine/vmatfind.c` (`processleafedge` / `processbranch` pair emission order).
  - `Vmengine/fquery.c`, `procexqu.c` and `Vmatch/procfinal.c`:
    - `-p` is a query of revcomp(S) streamed by query start;
    - `fetchpositions` drops rows with relpos1 > relpos2.

### Reference build
- Stock: `apt` `vmatch` 2.3.1+dfsg-9.
- Source build (§WP8). Besides `VM_NOPRUNE` (`extendedleftSEP` reach = UINT_MAX), WP15 adds `VM_NOPRUNE_H`, which passes `searchlength = UINT_MAX` to `extendmismatchesleft` in `extendHD.c`.
- Without the env switches, the rebuilt `vmatch.x` output is byte-identical to stock.
- Harness `scratchpad/wp15/cmp15.py` (+ `xc15`, `big15.py`) compares rows as multisets of (l, i, r, j, |d|), because the default output may repeat a row.

### Findings
1. **Default output.** Without `-allmax`, Vmatch prints one row per seed. The same match is printed once for every seed whose best extension it is.
   - Stock Vmatch rarely prints such duplicates: the left-extension shortcut stops seeds that cross another seed.
   - With the shortcut off, the duplicates are frequent.
   - The E-value criterion prefers a short exact or low-distance window to a longer maximal one. Example: `ATCTGGTGTACTCTGCCCACGACTATCGGTGTACTCTGC -l 15 -e 3` prints only (0,24,16,15,1).
2. **Hamming shortcut.** `extendHD.c` has its own shortcut (`searchlength`). It never changes the `-allmax` set (§WP8: stock = definition), but it does change the best-per-seed rows (2 982 / 6 000 direct cases differ between stock and shortcut-off).
3. **Distance label in `-allmax`.** `matchcontainer` keeps the first stored copy of a match, so Vmatch prints the distance of the **first seed** (in its seed order) that reaches the match. That can exceed the edit distance:
   - `CAAATTTT -p -l 8 -e 3 -allmax` prints (0,0,8,8,3) although d_E(CAAATTTT, AAAATTTG) = 2;
   - `ACAAAAAAAACAC -l 7 -e 3 -allmax` (stock) prints (0,2,11,11,3), where d_E = 2.

   `FindDegenerateRepeats` keeps reporting the edit distance by default. The 6 000-case WP8 harness had not hit such a case, but the WP15 run did: 1 case, 2 rows, vs the shortcut-off build. `vmatchCompatible` reproduces the label by extending seeds in Vmatch's order: `vmatfind.c` emission order for direct seeds, query start for `-p`.

### Implementation (WP15)
- `DegenerateRepeatReporting { AllMaximal, BestPerSeed }`.
- New overloads:
  - `FindDegenerateRepeats(…, reporting, vmatchCompatible = false)`;
  - `FindApproximateDirectRepeats(string, …, excludeContained, reporting, vmatchCompatible = false)`.
- MCP `find_degenerate_repeats` / `find_approximate_direct_repeats` get the optional parameters `reporting` (`allMaximal` | `bestPerSeed`) and `vmatchCompatible`.
- **Literal port, `VmatchExtension`.** Hamming tables, edit fronts with every convention above, and candidates produced in Vmatch's order. `VmatchEvalues` builds the table exactly as `incprecomputehammingEvalues` (same double arithmetic order), and `CompareVmatchMatches` = `cmpmatches`.
- **`vmatchCompatible`.** The reach is the seed length (stock shortcut). For edit `-allmax`, candidates go through the per-seed and global containers with first-seed-wins distances.

### Cross-checks (0 = no differing row; reference: stock `vmatch` for compatible, `VM_NOPRUNE`+`VM_NOPRUNE_H` build for complete)
The edit `-allmax` complete row is the existing F47 default; it reports the edit distance, not the first-seed label.

| Inputs | Configurations | Rows (Vmatch) | Differing rows |
|---|---|---|---|
| 6 000 × 8–50 bp (k 1–3; ℓ = k+1 … k+1+max(8, n/6); uniform / AC / 5 % N / motif copies with indels / single N) | 16 = {h, e} × {direct, -p} × {allmax, best} × {stock, complete} | 1 423 142 | 0 in all 15 configurations except edit -p `-allmax` complete (default): 1 case / 2 rows, the first-seed distance label (finding 3) |
| 2 000 × 8–50 bp, k 1–4 | 8 Hamming configurations | 233 692 | 0 |
| 300 × 100–1 500 bp | 16 | 2 338 723 | 0 |
| 1 Mb random + 500 planted 100-bp copies (half reverse-complemented, 5 % noise), ℓ 30, k 2 | e/h × d/p best stock, e d/p allmax stock, e/h d best complete | 6 237 | 0 |
| 200 kb, 100 planted copies, ℓ 20, k 3 | 10 (all stock + complete best) | 12 929 | 0 |

The shortcut changes the output in:

| Configuration | Cases changed (of 6 000) |
|---|---|
| direct edit `-allmax` | 109 |
| `-p` edit `-allmax` | 97 |
| best-per-seed edit | 3 624 (direct) / 2 961 (`-p`) |
| best-per-seed Hamming | 2 982 (direct) / 2 773 (`-p`) |

Every one of those cases matches Vmatch in both settings.

Timing (Release):
- 1 Mb (ℓ 30, k 2): best-per-seed 1.0–2.7 s, against stock Vmatch 0.1–0.5 s.
- 200 kb (ℓ 20, k 3, seed length 5): 1.5–27 s, against 3–16 s.

### Locked values (`RepeatFinder_VmatchReporting_Tests`, MCP tests)
Tuples (i, j, l, r, d). The 121-bp input is `CCGGCCCCTGAGTCCGAGGAGGATCACAGTCTACACTGCTCACTCCAACCGAGGATCACAGTTTACACTGCTCACTCCAACCGAGGGTGCTTGGATCACAGTCTACATGCTCACTCCAACC`.

| Call | Stock | Shortcut off |
|---|---|---|
| 121-bp copies `-l 12 -h 2` | 18,50,36,36,1; 20,92,15,15,0; 31,63,23,23,0; 36,107,14,14,0; 52,92,15,15,1; 68,107,14,14,0 | 18,50,36,36,1 ×2; 20,92,15,15,0; 36,107,14,14,0; 52,92,15,15,1 ×2; 68,107,14,14,0 |
| same `-l 12 -e 2` | 18,50,36,36,1; 20,92,30,29,1; 31,63,23,23,0; 36,107,14,14,0; 52,92,30,29,2; 63,103,19,18,1; 68,107,14,14,0 | 18,50,36,36,1 ×2; 20,92,30,29,1 ×2; 52,92,30,29,2 ×3 |
| same `-p -l 12 -h 2` | 26,26,12,12,2; 26,58,12,12,2; 58,58,12,12,2 | each ×2 |
| `ATCTGG…GTGTACTCTGC -l 15 -e 3` / `-h 3` | 0,24,16,15,1 / 1,24,15,15,3 | 0,24,16,15,1 ×2 / 1,24,15,15,3 |
| `GGACCATGAAGG -p -l 5 -e 3` | 9 rows | 19 rows |
| `-allmax` `CAAATTTT -p -l 8 -e 3` | 0,0,8,8,3 | default 0,0,8,8,2 (edit distance) |

