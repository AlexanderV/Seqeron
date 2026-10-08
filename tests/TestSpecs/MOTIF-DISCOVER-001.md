# Test Specification: MOTIF-DISCOVER-001

**Test Unit ID:** MOTIF-DISCOVER-001
**Area:** Matching
**Algorithm:** Motif Discovery via Overrepresented k-mers (observed/expected enrichment)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-29

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Compeau & Pevzner, *Bioinformatics Algorithms*, expected k-mer occurrences | 1 | https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string | 2026-06-14 |
| 2 | monaLisa `getKmerFreq` (O/E ratio, log2 enrichment) | 3 | https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html | 2026-06-14 |
| 3 | RSAT oligo-analysis (van Helden et al. 1998, J Mol Biol 281:827) — Bernoulli expected frequency, ratio | 1 | https://raw.githubusercontent.com/rsa-tools/rsat-code/master/perl-scripts/oligo-analysis | 2026-09-29 |

### 1.2 Key Evidence Points

1. Expected occurrences of a specific k-mer in a length-N string under the i.i.d. uniform (each base p=1/4) background is `E = (N − k + 1) / 4^k` — Source 1.
2. Overrepresentation is the observed/expected (O/E) ratio; value > 1 means overrepresented — Sources 1, 2.
3. `N − k + 1` is the number of length-k windows; `4^k` is the number of distinct DNA k-mers — Source 1.

### 1.3 Documented Corner Cases

- `k > N`: zero length-k windows, no k-mer can be counted (no motifs returned) — Source 1.
- The self-overlap caveat affects only the *probability* statistic, not the deterministic observed count or the O/E denominator — Source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Using a floor/clamp on the expected count (e.g. `max(E, 0.1)`) is not part of the published statistic and distorts the O/E ratio — Source 1 (formula has no clamp).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DiscoverMotifs(DnaSequence sequence, int k = 6, int minCount = 2)` | MotifFinder | Canonical | Returns `DiscoveredMotif` records (Sequence, Count, Positions, Enrichment) |
| `DiscoverMotifs(DnaSequence, int k, int minCount, IReadOnlyList<double> background)` | MotifFinder | Overload | Bernoulli background (RSAT); uniform background ≡ default |
| MCP `discover_motifs` | AnalysisTools | Wrapper | delegates to the default overload |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Count equals the number of occurrences of the k-mer; Positions are the 0-based window starts of those occurrences | Yes | Source 1 (window enumeration) |
| INV-2 | Enrichment = Count / E where E = (N − k + 1) / 4^k | Yes | Source 1 |
| INV-3 | Every returned motif has Count ≥ minCount | Yes | Method contract |
| INV-4 | Enrichment > 0 for every returned motif (E > 0 when any k-mer exists) | Yes | Source 1 (N−k+1 ≥ 1) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Repeated k-mer found | "ATGCATGCATGC", k=4, minCount=2 | result contains "ATGC" with Count=3 | Source 1 (window count) |
| M2 | Positions reported | "ATGCATGCATGC", k=4, minCount=2 | "ATGC" Positions = {0,4,8} exactly | INV-1 |
| M3 | Exact O/E enrichment (tandem) | "ATGCATGCATGC", k=4 | "ATGC" Enrichment = 768/9 = 85.3333… (3 / (9/256)) | Source 1 formula |
| M4 | Exact O/E enrichment (homopolymer) | "AAAAAAAAAA", k=3 | "AAA" Count=8, Enrichment = 64.0 (8 / (8/64)) | Source 1 formula |
| M5 | minCount filter | "ATGCAAAA", k=4, minCount=2 | every returned motif has Count ≥ 2 | INV-3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null sequence | sequence = null | ArgumentNullException | Validation contract |
| S2 | k < 1 | k = 0 | ArgumentOutOfRangeException | Validation contract |
| S3 | k > N | "AAA", k=5 | empty result | Corner case (no windows) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | No floor on expected count | very long sequence so E > 0.1, plus a long-k case where the old `max(E,0.1)` clamp would have changed the value | Enrichment = Count / ((N−k+1)/4^k) exactly | Guards against re-introducing the clamp defect |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinderTests.cs` — `#region Motif Discovery Tests` contains `DiscoverMotifs_FindsRepeatedKmer`, `DiscoverMotifs_ReturnsPositions`, `DiscoverMotifs_CalculatesEnrichment`, `DiscoverMotifs_FiltersByMinCount`, plus `DiscoverMotifs_NullSequence_ThrowsException`, `DiscoverMotifs_ZeroK_ThrowsException`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (repeated k-mer found) | ⚠ Weak | Existing `DiscoverMotifs_FindsRepeatedKmer` uses `Any(...)`, no Count assertion |
| M2 (positions) | ⚠ Weak | `DiscoverMotifs_ReturnsPositions` uses `Does.Contain` (no exact set / boundary) |
| M3 (exact O/E tandem) | ❌ Missing | No exact enrichment value asserted anywhere |
| M4 (exact O/E homopolymer) | ⚠ Weak | `DiscoverMotifs_CalculatesEnrichment` only asserts `> 1` (permissive) |
| M5 (minCount filter) | ✅ Covered | `DiscoverMotifs_FiltersByMinCount` asserts Count ≥ 2 (kept; mirrored in canonical file) |
| S1 (null) | ✅ Covered | `DiscoverMotifs_NullSequence_ThrowsException` (mirrored) |
| S2 (k < 1) | ✅ Covered | `DiscoverMotifs_ZeroK_ThrowsException` (mirrored) |
| S3 (k > N) | ❌ Missing | Not tested |
| C1 (no floor) | ❌ Missing | Not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_DiscoverMotifs_Tests.cs` — all MUST/SHOULD/COULD cases with exact evidence values.
- **Remove:** the six weak/duplicate `DiscoverMotifs_*` tests from `MotifFinderTests.cs` (their coverage is superseded by the canonical file with exact values).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_DiscoverMotifs_Tests.cs` | Canonical | 9 |
| `MotifFinderTests.cs` (`Motif Discovery` region) | removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote with exact Count=3 | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote with exact Positions {0,4,8} | ✅ Done |
| 3 | M3 | ❌ Missing | Added exact Enrichment 768/9 | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewrote with exact Enrichment 64.0 + Count 8 | ✅ Done |
| 5 | M5 | ✅ Covered | Mirrored into canonical file | ✅ Done |
| 6 | S1 | ✅ Covered | Mirrored | ✅ Done |
| 7 | S2 | ✅ Covered | Mirrored | ✅ Done |
| 8 | S3 | ❌ Missing | Added k>N empty-result test | ✅ Done |
| 9 | C1 | ❌ Missing | Added no-floor exact-ratio test | ✅ Done |

**Total items:** 9
**✅ Done:** 9 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact Count=3 in canonical file |
| M2 | ✅ | Exact Positions {0,4,8} |
| M3 | ✅ | Enrichment = 768/9 within 1e-10 |
| M4 | ✅ | Enrichment = 64.0, Count=8 |
| M5 | ✅ | minCount filter asserted |
| S1 | ✅ | ArgumentNullException |
| S2 | ✅ | ArgumentOutOfRangeException |
| S3 | ✅ | k>N empty |
| C1 | ✅ | No-floor exact O/E |

**Total in-scope cases:** 9 | **✅:** 9

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `minCount` is a presentation filter, not part of the published statistic (not correctness-affecting per record) | M5 |

---

## 7. Open Questions / Decisions

1. The checklist signature reads `DiscoverMotifs(sequences, k)`, but the registered Type is "Overrepresented k-mers" and the implemented canonical method operates on a single `DnaSequence` (per-sequence k-mer overrepresentation). The single-sequence overrepresentation method is the canonical one; the cross-sequence variant is the separate unit MOTIF-SHARED-001 (`FindSharedMotifs`). External evidence (Compeau & Pevzner expected-count formula) defines the single-sequence statistic, so testing targets `DiscoverMotifs(DnaSequence, k, minCount)`.

## 8. Review 2026-09 additions (B05)

| ID | Test | Expected | Evidence |
|----|------|----------|----------|
| R1 | Bernoulli background A.3/C.2/G.2/T.3, "ATGCATGCATGC", k=4 | ATGC 92.5925925925926, TGCA/GCAT/CATG 61.72839506172839 | RSAT exp_occ = ∏q · (N−k+1); Python Fractions |
| R2 | Unnormalised background (3,2,2,3); "A"×10, k=3 | ATGC 92.5926; AAA 37.03703703703704 | normalisation as CreatePwm; RSAT |
| R3 | Background (1,1,1,1) overload | identical (Sequence, Count, Enrichment) to default, 50 random cases | RSAT equiprobable ≡ Bernoulli ¼ |
| R4 | Invalid background (null, 3 values, zero, NaN) | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |
| R5 | X+X, X = LCG 512-mer, k=512 | count 2, positions {0,512}, enrichment 7.008550233381348e+305 (finite; was +∞); bg (26,24,24,26) → 3.2383751592973165e+306 | exact Fractions |
| R6 | Counts vs `KmerAnalyzer.CountKmers`; positions ascending; first-occurrence order | equal; ascending; ACGT,CGTA,GTAC,TACG,GTAA | canonical counter |

## 9. Review 2026-09 follow-up — RSAT oligo-analysis significance, Markov backgrounds, both strands

Method: `MotifFinder.DiscoverMotifs(DnaSequence, int k, int minCount, OligoBackgroundModel, OligoStrandMode, bool countOverlapping)` → `OligoAnalysisResult`.
Reference: RSAT `oligo-analysis` 1.169 (rsa-tools/rsat-code 10043f2) run with perl 5.38, full-precision debug print after `MultiTestCorrections`; independent Python port (scipy `binom.sf`) of the same source agreeing with RSAT to ≤ 2e-12 on 68 configurations (plus 18 `-bgfile` runs to RSAT's %5g) and with the C# code to ≤ 3e-13 on 400 random configurations (all backgrounds, strands, `-noov`).
Test file: `Unit/Analysis/MotifFinder_OligoAnalysis_Tests.cs`; binomial tail: `Unit/Core/StatisticsHelper_BinomialUpperTail_Tests.cs`.

| ID | Test | Expected | Evidence |
|----|------|----------|----------|
| S1 | t2 (63 nt), k=4, `-1str -bg equi -lth occ 2` | n 60, tested 10, NPO 256; ATGC occ 6, occ_P 1.4844942103072834e-07, occ_E 1.4844942103072835e-06, occ_sig 5.8284214918614703; TGCA occ 4, occ_P 9.5343622068379445e-05 | RSAT run; scipy binom.sf(5,60,1/256) |
| S2 | t2, k=4, `-2str` input Bernoulli, `-lth occ 2` | tested 15, NPO 136; ATGC|GCAT occ 9, exp_freq 0.0077771093320170084, occ_P 1.0760647330191343e-09, occ_sig 7.7920703428970466; CATG palindrome occ 5, occ_P 4.0637121408222015e-06 | RSAT run |
| S3 | t2, k=4, `-markov 1` / `-markov 2 -lth occ 3` | ATGC exp_freq 0.032915191520534237, occ_P 0.013954892377198637, occ_E 0.55819569508794542 (40 tested); GCAT exp_freq 0.055540625279942669, occ_P 0.65459731021733769, occ_sig −0.41803420775951439 (4 tested) | RSAT run |
| S4 | Markov table f(AA..TT) = 1..16 | ψ=0: exp_freq(ATGC) = 600/331296 (hand); ψ=0.01: 0.0018478880336134456 (RSAT prints %5g 0.00184789), occ_P 1.8299363778610906e-09, occ_sig 8.436534013635193 | RSAT `-bgfile`; MarkovModel.pm port |
| S5 | t3 (35 nt), k=2, `-noov` | AA occ 5 at {1,3,5,27,29}, occ_P 0.058469994275128084 (10 tested); `-2str -noov` AA|TT occ 10, occ_P 0.0070505538526756187; `-2str -ovlp` occ 17, occ_P 1.2344994901541932e-07 | RSAT run |
| S6 | t2, `-2str -markov 1 -lth occ 3` | ATGC|GCAT exp_freq 0.06289152665530649 = p(ATGC)+p(GCAT), occ_P 0.012300624704130635, occ_sig 1.3080128404304179 | RSAT manual formula (RSAT code: 2·p(min), a quirk); Python port |
| S7 | Equiprobable single strand vs legacy overload | Ratio = Enrichment (1e-12), same words/positions/order | definition |
| S8 | `-2str` count | occ(W|W') = occ(W)+occ(W'), pairs partition the N−k+1 windows | RSAT `SumReverseComplements` |
| S9 | invalid arguments | null → ArgumentNullException; k<1, Markov order > k−2, table order ≥ k, ψ∉[0,1] → ArgumentOutOfRangeException; empty/mixed/non-ACGT table, ψ=0 null transition → ArgumentException | RSAT FatalError conditions |
| S10 | X+X, X 600-mer, k=600 | occ_P underflows to 0 but occ_sig = −(ln C(601,2) − 1200 ln 4 + ln T)/ln 10 (finite); NPO = +∞ | log-space tail |
| B1 | `StatisticsHelper.BinomialUpperTail` / `LogBinomialUpperTail` | 7 cases = mpmath exact sums (scipy ≤ 5e-14); p = e^−800: ln tail −1586.8786371229292547 (mpmath); bounds; pmf differences | Loader 2000; RSAT sum_of_binomials |

Heavy tier: `Properties/MotifOligoAnalysisProperties.cs` (definitions consistency, monotone in occ, Markov-0 = input Bernoulli, uniform composition Markov-0 = 4^−k, `-2str` = forward + revcomp, ms consistency), `Metamorphic/MotifOligoAnalysisMetamorphicTests.cs` (reverse-complement input invariance under `-2str`), `Fuzzing/MotifOligoAnalysisFuzzTests.cs`.

## 10. Review 2026-09 audit group E — RSAT oligo-analysis options and dyad-analysis

Methods: `MotifFinder.AnalyzeOligos(sequences, k, OligoAnalysisOptions)` → `OligoAnalysisReport`; `OligoBackgroundModel.Lexicon`;
`OligoCalibration.Parse/FromEntries`; `MotifFinder.FindSharedMotifs(…, double pseudoFrequency)`; `MotifFinder.AnalyzeDyads(sequences, DyadAnalysisOptions)`;
`StatisticsHelper.LogPoissonRangeProbability` / `LogNegativeBinomialRangeProbability`.
Reference: RSAT `oligo-analysis` 1.169 / `dyad-analysis` 1.78 (rsa-tools/rsat-code 10043f2) run with perl 5.38 and a %.17g dump of every
pattern field before `PrintResult`; oracle fixes only where RSAT code is defective (stated per row); mpmath / scipy for Poisson / negative binomial.
Random cross-checks (C# vs RSAT): z-score 100 runs / 3,390 patterns ≤ 1.5e-14; `-pseudo` 100 / 3,128 ≤ 1.5e-14; `-lexicon` 100 / 3,620 ≤ 1.6e-14;
`-oneN`/`-onedeg` (fixed oracle) 100 / 26,168 ≤ 2.9e-14 (z ≤ 4e-13); calibration (exact-sum oracle) 100 / 3,822 ≤ 1.7e-14; mseq `-pseudo` 40 runs ≤ 4.7e-14;
dyad single strand 200 runs / 14,746 dyads ≤ 1.3e-14, `-expfreq` 40/40, both strands vs an independent Python port 150/150.
Test files: `Unit/Analysis/MotifFinder_OligoAnalysisOptions_Tests.cs`, `Unit/Analysis/MotifFinder_DyadAnalysis_Tests.cs`, `Unit/Core/StatisticsHelper_PoissonNegBinRange_Tests.cs`; MCP `OligoAnalysisOptionsAndDyadTests`.

| ID | Test | Expected | Evidence |
|----|------|----------|----------|
| E1 | Default options vs `DiscoverMotifs` (equi, input `-2str` ±`-noov`, Markov 1/2, lexicon) | bit-identical words, counts, positions, exp_freq, exp_occ, ratio, occ_P/E/sig, tested | definition |
| E2 | t2 k=4 `-1str -bg equi -return occ,proba,zscore -lth occ 2` | ATGC exp_var 0.22613525390625, z 12.124455829017547; TGCT ovlp 1.015625, exp_var 0.23345947265625, z 3.6542034172140836 | RSAT run |
| E3 | t2 k=4 `-2str -noov` / `-ovlp` zscore | noov: ATGC occ 6, overlaps 3, exp_var = exp_occ 0.46662655992102053, z 8.1003774086907097; ovlp: exp_var 0.43396550795746169, z 12.953679437171939 | RSAT run |
| E4 | `-markov 1 -pseudo 0.1`; `-2str -bg equi -pseudo 0.2` | ATGC exp_freq 0.030014297368480814, occ_P 0.0091641466020121916, z 3.6626693266117103; pair exp_freq 0.0091911764705882356 = 2(0.8/256 + 0.2/136) | RSAT run |
| E5 | `-lexicon` k=4; k=5 with `-pseudo 0.01` | ATGC exp_freq 0.025, segments a·tgc (0.25, 0.1), occ_P 0.003853224500294343; ATGCA 0.025605858701702098, CATGC occ_E 0.10730819082725597 | RSAT run |
| E6 | `-onedeg` k=3 `-lth occ 6`; `-2str -oneN` k=4 | NPO 528, 71 tested, AYG occ 8 exp_freq 0.03071422572556359 occ_P 0.00054574388255023234; ATGN|NCAT occ 11 occ_P 3.5186386968671786e-06 | RSAT with `Degenerate` fixed (RSAT 1.169: 0 tested, no output) |
| E7 | `-calibN` (cal3.tab) k=3 | ATG negbin occ_P 0.0067142352257039137 (RSAT 5-digit terms 0.0067141894090253264); TGC Poisson 0.0044559807752478468 (RSAT `$prev_value` defect: 0.003529988861723204); GCA 0.009079857800153978, z 4 | RSAT run + exact-sum oracle; mpmath |
| E8 | `-calib1` two sequences `-2str -noov` | ATG|CAT occ 6, overlaps 5, exp_occ 1.8, exp_var 2.88, exp_freq 0.092307692307692313, occ_P 0.03602153062820438 | exact-sum oracle |
| E9 | mseq `-markov 1 -pseudo 0.1` (ms.fa) | ACGT exp_ms 1.1451050486269749, ms_P 0.073696629017800608, ms_E 18.866337028556956 | RSAT run |
| E10 | Calibration parsing, rc inference, errors | `aac` → GTT inherits (1.5, 1); FormatException / ArgumentException | RSAT `ReadCalibration` |
| D1 | dyad t2 `-l 3 -sp 0-2 -1str -lth occ 2` | 13 tested, spacings (58,54,4),(57,54,3),(56,54,2); ATGn{0}CAT occ 2, overlaps 1, exp_freq 0.008062348830959418, occ_P 0.079911679827561213, z 2.41 | RSAT run |
| D2 | dyad `-l 2 -sp 0-6 -type dr -1str -ovlp -lth occ 2` | NPD 112, 6 tested, ATn{2}AT occ_P 0.015985324783370104 | RSAT run |
| D3 | dyad RSAT defaults (`-2str -noov`) | spacings (58,50,8)…; ATGn{0}CAT occ_P 0.079911679827561213; 11 tested (RSAT 9: pairs seen only as the larger member lost) | RSAT run + pair rule |
| D4 | `-type rep` palindromic monads; `MinCount = 0`; dyad / monad tables | counted once (RSAT twice); 3267 tested; table value / monad fallback | brute force; RSAT run |
| P1 | Poisson / negative binomial range probabilities | 8 + 5 cases = mpmath (≤ 1e-12) | mpmath; scipy |

## 11. Review 2026-09 audit round 2 (G2 + G4) — `-seqtype` and degenerate matching sequences

Tests: `Unit/Analysis/MotifFinder_OligoSequenceType_Tests.cs` (10), `Unit/Analysis/MotifFinder_SharedMotifsDegenerate_Tests.cs` (6),
MCP `OligoSequenceTypeAndDegenerateSharedTests` (5). Oracle: RSAT oligo-analysis 1.169 copies `oligo-analysis-g2` / `-g4`
(Evidence, audit round 2).

| ID | Case | Expected | Source |
|----|------|----------|--------|
| S1 | `-seqtype prot -l 2` MKLLVAAGLLKLMKXLLA* / mkllvqqKLLAA | n 26 (3 windows discarded), NPO 400, 14 tested; LL occ 5 exp_freq 0.14387633769322233 ovlp 1.3793103448275863 z 0.63864701205606378 occ_P 0.31619335865366638; KL occ_P 0.0865394727065604 | RSAT run |
| S2 | `-seqtype other -bg equi -l 2` "Hello, World! hello world. AbC abc" | alphabet 13, NPO 169, n 28, 18 tested; ll ovlp 1.0769230769230769 z 4.250165852579884 occ_P 0.01194994022319412; ab occ 2 (case folded) | RSAT run |
| S3 | `-seqtype prot -markov 1 -l 3` LLLLAAAAMKLLLAA | AAA exp_freq 0.2040816326530612, LLL 0.27332361516034986; ovlp 1.0525 (RSAT fallback 1.3125 / 1); z NaN (exp_var < 0) | RSAT run, ovlp fix |
| S4 | DNA with N/R `-2str -l 4` ACGTNACGTACGRTTACGTAC / ttacgtnnACGTAC | n 16, NPO 136, 4 tested; ACGT occ 5 occ_P 3.4909208986462239e-09 z 20.294579804259762; GTAC occ_P 3.0388746348848393e-05 | RSAT run |
| S5 | DNA `-markov 2 -l 5` with Y / N | uncounted sub-word → exp_freq 0, untested (RSAT: division by zero) | RSAT run, fix |
| S6 | Pure-ACGT strings (white space, lower case) vs `DnaSequence` path | bit-identical (5 option sets) | identity |
| S7 | Protein + both strands / degeneracy / calibration / DNA-only background; unknown type; null | ArgumentException / ArgumentOutOfRangeException / ArgumentNullException | contract |
| M1 | `-l 3 -1str -onedeg -return mseq,proba` ACGTACGGATCC / ATGCATGAAC / ACGATGTT | NPO 528; AYG mseq 3 (RSAT sum 4 > 3) exp_freq 0.032666666666666663 exp_ms 0.69998598816389546 ms_P 0.01270294085234148 ms_E 6.7071527700363012; ACN ms_P 0.41105031540287185; ATN 0.085427115865445352 | oracle `-g4` |
| M2 | `-l 3 -2str -oneN -lth mseq 2` | NPO 24; ANG|CNT exp_freq 0.12444444444444444 exp_ms 1.9639256569777699 ms_P 0.28055070694211726; ACN|NGT 0.72456658149428699; ATN 0.35326100605480787 | oracle `-g4` |
| M3 | mseq = brute-force union of matching sequences; `None` = plain overload; `-pseudo` uses the degenerate NPO | identity | definition |

