# Evidence Artifact: ASSEMBLY-REPEAT-001

**Test Unit ID:** ASSEMBLY-REPEAT-001
**Algorithm:** Assembly repeat analysis (`GenomeAssemblyAnalyzer.FindRepetitiveRegions`, `FindTandemRepeats`, `CalculateRepeatContent`, `RepeatAnnotation`; MCP `find_repetitive_regions`, `find_tandem_repeats`, `repeat_content`) and the linguistic-complexity routing of `CalculateLocalQuality`
**Date Collected:** 2026-11-25 (B10 sweep U3)

---

## Online Sources

### NCBI C++ Toolkit — WindowMasker (`src/algo/winmask`, ncbi/ncbi-cxx-toolkit-public master)

**URL:** https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/winmask/{seq_masker.cpp, seq_masker_window.cpp, seq_masker_score_mean.cpp, win_mask_gen_counts.cpp} and `include/algo/winmask/seq_masker_window.hpp`
**Accessed:** 2026-11-25 (curl, neutral User-Agent; files read in full)
**Authority rank:** 3 (reference implementation of Morgulis et al. 2006, Bioinformatics 22:134; the paper itself was not opened)

**Key Extracted Points:**

1. `CSeqMaskerWindow::FillWindow`: a window of `window_size` letters must consist of A/C/G/T only (`LOOKUP` 1..4; any other letter, including lower case and N, resets `iter`); `Advance` refills when it meets such a letter → only all-ACGT windows are scored.
2. `CSeqMaskerScoreMean`: score = `sum / num` over the window's `num = window_size − unit_size + 1` units (unsigned integer division).
3. `CSeqMasker::DoMask` (first pass; trigger = mean, so `ts = s`): `s < textend` closes an open interval when `window.Start() > cend`; `s ≥ textend` and `ts < cutoff` keeps the interval connected (`Start() > cend + 1` closes it, else `cend = window.End()`); `ts ≥ cutoff` opens or extends it (`start = Start()`, `cend = end = End()`). Intervals are `[start, end]` inclusive, 0-based.
4. `win_mask_gen_counts.cpp`: counting walks each sequence, resets on an ambiguous letter, and increments `counts[unit]` when `unit ≤ runit` and `counts[runit]` when `runit ≤ unit` → canonical (both-strand) counts; a self-complementary unit satisfies both tests and is incremented twice per occurrence.
5. `T_threshold` / `T_extend` are taken from percentiles of the unit-count distribution; the merge pass (`merge_pass`) and discontiguous units are optional extras.

### RepeatMasker — `ProcessRepeats` (rmhubley/RepeatMasker master)

**URL:** https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/ProcessRepeats
**Accessed:** 2026-11-25 (curl; `generateTableOutput` and the aggregate-statistics loop read)
**Authority rank:** 3

1. `.tbl` header: "bases masked: n bp (p %)", `p = maskedlength / totseqlen · 100`; `maskedlength` is the length of the masked sequence (`maskSequence`) or the aggregated `annotatedlength`.
2. Aggregation loop comment: "Define length of coverage based on possibly overlapping fragments": annotations are walked in query order; an annotation overlapping the previous one on the same sequence contributes only the bases beyond the previous end (`length = End − Start + 1 − (pastEnd − Start + 1)`), a contained one contributes no bases ("Allow it to contribute to annotation count but don't count it's bases"); the length is added to the annotation's class (`totlength{class}`) and to `annotatedlength` → every base counted once, class lengths add up to the total.

### Tandem Repeats Finder — Benson 1999 (NAR 27:573) / TRF 4.10.0

Opened earlier in this project: docs/Evidence/REP-APPROX-001-Evidence.md (TRF 4.10.0 README + compiled TRF parity: 99.8–100 % of rows identical on 700 random sequences). `RepeatFinder.FindApproximateTandemRepeats` is that implementation; **the TRF source was not rebuilt in this session** (the sandbox refused to build the cloned AGPL source), so TRF-level parity is inherited from REP-APPROX-001 and complemented by an independent scorer below.
TRF's definitions used: consensus pattern, copy number (aligned consensus columns / consensus size), "percent matches between adjacent copies", alignment score (match 2, mismatch −7, indel −7), Minscore 50, MaxPeriod ≤ 2000.

### Linguistic complexity — canonical `SequenceComplexity` (Troyanskaya et al. 2002; Orlov & Potapov 2004; Rosalind LING)

Defined in docs/Evidence of the complexity unit; LC = Σ V_i / Σ min(a^i, N − i + 1), m = 10, a = 4 extended by any other symbol (e.g. N). The private routine removed here (distinct 4-mers / min(valid 4-mers, 256), N-words skipped) is not a published measure.

---

## Reference numbers (oracles in `tools/review/oracles/B10/`, `python3 -I`)

### WindowMasker port (`repeat_windowmasker_oracle.py`)

Literal port of `FillWindow`/`Advance`, the mean score, `DoMask` first pass and the canonical counting of `win_mask_gen_counts.cpp` (double increment of palindromes). 3000 random cases (seed 1; k ∈ {3..8, 15}; W = k..k+40; T_threshold 1–6, T_extend 1..T; 1–3 sequences with planted repeat motifs and reverse-complement copies, N/R/Y/n symbols, lower case; forward-only and both-strand counting) → 9391 regions, C# `FindRepetitiveRegions` identical on **3000/3000** cases. 14 of them are locked as `WindowMaskerOracle` fixtures.
Hand values: (CAG)×200, k 6, T 3, W 60 → [0, 599], copies 199; ACGTACGT×100, k 8, W 50 → [0, 799], copies 398 (self-complementary unit: 199 × 2 strands); ACGTG, ACGTG, CACGT with k = W = 5 → no region forward-only, three regions with copies 3 on both strands.

### TRF alignment score (`repeat_trf_wdp_oracle.py`)

Independent cyclic wraparound DP (local alignment of the sequence against unlimited copies of the reported consensus, 2/−7/−7). 400 random flanked arrays (period 2–25, 4–25 copies, substitution 0–15 %, indel 0–5 %): 348 reported by TRF; the best repeat's alignment score equals the independent score in **346/348**; the two others have consensus size 21 and 24 (> 20: TRF's narrow-band WDP, expected to be ≤ the full DP: 265 vs 269, 229 vs 233).
Hand values (all reproduced): (AT)×30 → [0, 59], AT, 30 copies, 100 %, score 120; (CAG)×30 → score 180; (ATG)×7·ATC·(ATG)×7 → [0, 44], 15 copies, 40/42 = 95.238 % matches, score 2·44 − 7 = 81; ACGTGCATTG·(CAG)×9·TGCATC → [9, 36] GCA, 28/3 copies, score 56; (CAG)×15 + 12 N + (CAG)×15 bridged into one repeat (−84 < +90), with 40 N split into [0, 44] and [85, 129].

### Linguistic complexity

300 random windows (A, AC, ACGT, ACGTN alphabets, length 1–200): independent Python Σ V_i / Σ min(a^i, N−i+1) equals `SequenceComplexity.CalculateLinguisticComplexity` on **300/300** (tolerance 1e-12). Hand: "ACGTACGT" → 26/32 = 0.8125 (old private 0.8); poly(A) 500 → 10/3301 (old 1/256); "ACNGT" → 1.

### Repeat content

Per-base marking oracle (annotation first in (start, input order) owns the base): 2000 random annotation sets (1–3 sequences, 3 classes, overlaps and containment) → total, percentage, class lengths and class sum identical (in-test oracle `CalculateRepeatContent_RandomAnnotations_MatchPerBaseOracle`).

---

## Findings (B10 F48–F53)

- **F48** `FindRepetitiveRegions`: unsourced window heuristic (window step `windowSize/10`, "more than `windowSize/(2k)` high-copy k-mers", region end = start of the first non-repetitive window − 1, so ends overshot/undershot the evidence) replaced by the WindowMasker masking model (source 1–3). `bothStrands` (default, WindowMasker) and `extendCopies` options.
- **F49** `FindRepetitiveRegions` N rule: private loops skipped only windows containing `N` — a window with R/Y/other IUPAC or `-` was counted as a distinct k-mer. Now `KmerAnalyzer.CountKmers` with `AcgtOnly` (Jellyfish / WindowMasker rule): any non-ACGT symbol excludes the unit and the window (B06 request, DUP_MAP).
- **F50** `FindTandemRepeats`: unsourced model (first copy as the unit, ≥ 80 % identity per following copy, substitutions only, no scoring, integer copies, 8-bp (AT)×4 reported) replaced by Tandem Repeats Finder via `RepeatFinder.FindApproximateTandemRepeats`; B04's "legitimately different" is not borne out — the 80 % rule has no source (TRF reports by alignment score ≥ 50 with indels allowed). `Copies` is now `double` (TRF copy number), `Unit` the TRF consensus, `Purity` TRF percent matches / 100, new `minScore`.
- **F51** `CalculateRepeatContent`: overlapping/contained annotations were summed twice (repeat content could exceed 100 %); now RepeatMasker's once-per-base accounting; negative/inverted intervals, negative genome length and more masked bases than genome bases rejected; null throws `ArgumentNullException`.
- **F52** `CalculateLinguisticComplexity` (private, backing `CalculateLocalQuality.Complexity`): unsourced single-k ratio replaced by canonical `SequenceComplexity.CalculateLinguisticComplexity` (m = 10, N as an extra symbol). Behaviour change: window complexity values (and thus `FindSuspiciousRegions` low-complexity scores, sweep U5) change; locked tests updated with oracle values.
- **F53** `FindRepetitiveRegions` with `windowSize` < 10 never advanced (`i += windowSize / 10` = 0): infinite loop. Gone with the rewrite; `windowSize` ≥ `kmerSize` is validated eagerly.

## Deviations kept (documented)

WindowMasker: T_threshold / T_extend are counts (percentile derivation needs a genome-wide count distribution), merge pass and discontiguous units omitted, T_low/T_high count-database cut-offs omitted. TRF: `minCopies` is an additional library filter (TRF itself reports ≥ 1.9 copies). Repeat content: when two overlapping annotations belong to different classes the overlap is credited to the one first in (start, input order); RepeatMasker credits by its resolved annotation order.
