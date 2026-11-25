# Evidence Artifact: ASSEMBLY-COMPLETE-001

**Test Unit ID:** ASSEMBLY-COMPLETE-001
**Algorithm:** Assembly completeness — BUSCO-like marker search (`GenomeAssemblyAnalyzer.AssessCompleteness`, `FindMarkerHits`) and k-mer spectrum evaluation (`EstimateCompletenessFromKmers`: Merqury completeness / QV, GenomeScope error rate / genome size); MCP `assess_completeness`, `estimate_completeness_from_kmers`
**Date Collected:** 2026-11-25 (B10 sweep U2)

---

## Online Sources

All fetched with a neutral User-Agent (`Seqeron-review/1.0`) and read in full where cited.

### BUSCO — `src/busco/busco_tools/hmmer.py` (ezlab/busco, GitLab master)

**URL:** https://gitlab.com/ezlab/busco/-/raw/master/src/busco/busco_tools/hmmer.py
**Authority rank:** 3 (reference implementation of Simão et al. 2015; Manni et al. 2021)

1. `_sort_matches`: matches above the profile score cut-off are complete (odb10: −2 ≤ (mean length − size)/σ ≤ 2; odb12: lower ≤ size ≤ 1.5·median), very large (counted complete) or fragment (shorter).
2. `_get_busco_percentages`: `s_percent = round(single_copy / total_buscos * 100, 1)`, `d_percent = round(multi_copy / total_buscos * 100, 1)`, `f_percent`, `complete_percent = round((single_copy + multi_copy) / total_buscos * 100, 1)`, `missing_percent` — **every class over n**; missing = n − S − D − F.
3. `record_results`: one-line summary `"C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}"`.

### compleasm — `compleasm.py` (huangnengCSU/compleasm, tag v0.2.6 and main)

**URL:** https://raw.githubusercontent.com/huangnengCSU/compleasm/v0.2.6/compleasm.py (and `main`)
**Authority rank:** 3 (Huang & Li 2023, Bioinformatics 39:btad595 — BUSCO-compatible classifier over genome alignments)

1. Records: `Protein_mapped_rate = (Protein_End − Protein_Start) / Protein_length` (PAF 0-based half-open); `Identity` from miniprot.
2. Valid hits (lite mode): `(Protein_mapped_rate >= min_length_percent) | (Identity >= min_identity)`; defaults `--min_identity 0.4`, `--min_length_percent 0.6`, `--min_complete 0.9` ("The length threshold for complete gene").
3. `record_1st_gene_label(by_length=False)`: none → Missing; one record → Single if rate ≥ min_complete else Fragmented; several → complete regions = rate ≥ min_complete; 0 → Fragmented, 1 → Single, > 1 → Duplicated if on different contigs, else `get_region_clusters` (sorted by start, merge when `last.start <= start <= last.stop`) → 1 cluster Single, else Duplicated.
4. Summary percentages over the total (`len(single_genes) / total_genes * 100`).
5. main (2025+) deprecates `--min_complete` in favour of lineage length cut-offs (`_is_complete`: without a cut-off `size >= 0.8 * qlen`); Seqeron keeps the aligned-fraction mode because nucleotide markers carry no lineage length table.

### BLAST+ defaults (word size / scoring)

`blastn` default task is megablast with `word_size 28`; the blastn task uses reward 2, penalty −3, gap open 5, gap extend 2 — already encoded and sourced in `SequenceAligner.BlastDna` (BLAST+ User Manual, application options).

### Merqury — `eval/qv.sh`, `eval/spectra-cn.sh`, `build/filt.sh`, `eval/kmerHistToPloidyDepth.jar` (marbl/merqury master)

**URL:** https://raw.githubusercontent.com/marbl/merqury/master/… ; jar source `javax/arang/kmer/HistToPloidyDepth.java` read from the jar in a shallow clone.
**Authority rank:** 3 (Rhie et al. 2020, Genome Biol 21:245)

1. `qv.sh` / `spectra-cn.sh`: `meryl difference output $asm.0.meryl $asm.meryl $read` (full read db), `ASM_ONLY` and `TOTAL` = `meryl statistics … head -n4 | tail -n1` (the "present" = instance count), `ERROR = 1-(1-ASM_ONLY/TOTAL)^(1/k)`, `QV = -10*log(ERROR)/log(10)`.
2. `spectra-cn.sh`: `meryl intersect output $asm.solid.meryl $asm.meryl $read_solid`; completeness = `100 * distinct(asm.solid) / distinct(read_solid)` (`head -n3` = distinct).
3. `filt.sh`: `filt` = line 2 last field of `kmerHistToPloidyDepth.jar` output = `effectiveMinDepth`; `meryl greater-than $filt` → solid = count > filt.
4. `HistToPloidyDepth.getPloidyDepthTable`: first line: `prevD = count − Double.MAX_VALUE`; each next line `d = count − prevCount`; `if (d * prevD < 0 && d > prevD)` the first time → `effectiveMinDepth = depth − 1`; default 0.
5. meryl counts canonical k-mers; windows with a non-ACGT base are not counted.

### GenomeScope 1.0 — `genomescope.R` (schatzlab/genomescope master)

**URL:** https://raw.githubusercontent.com/schatzlab/genomescope/master/genomescope.R
**Authority rank:** 3 (Vurture et al. 2017, Bioinformatics 33:2202)

1. `total_kmers = sum(as.numeric(x*y))`, `total_error_kmers = sum(error_kmers * x[1:error_xcutoff_ind])`.
2. Error rate: `uniroot` of `f1(x) = sum_{i=1..k} (1−x)^(k−i) x^i choose(k,i) * total_kmers − total_error_kmers` — whose closed form is `1 − (1 − E/T)^(1/k)` (the same form as Merqury's ERROR).
3. `total_len = (total_kmers − total_error_kmers) / (2*kcov)`; `bp_format` uses R `round`.
4. The error k-mers and `kcov` come from a negative-binomial mixture fit (`nlsLM`, 4-peak model, up to `NUM_ROUNDS` coverage trims) — not reproduced (see Deviations).

---

## Reference numbers (oracle `tools/review/oracles/B10/completeness_oracle.py`, `python3 -I`)

### Marker hits / classification (Biopython 1.88 `PairwiseAligner` local, +2/−3, open −7, extend −2; masking for further hits; 28-bp exact-word requirement)

300-bp marker G, 100-bp flanks L/R (oracle seed 41):

| Scenario | Oracle hits (id strand start–end identity coverage) | Class |
|---|---|---|
| L+G+R | c1 + 100–399 1.0 1.0 | Single |
| L+RC(G)+R | c1 − 100–399 1.0 1.0 | Single |
| two contigs | c1 + 100–399; c2 + 50–349 | Duplicated |
| L+G+G+R (abutting) | c1 + 100–399; c1 + 400–699 | **Single** (compleasm merges [100,400) and [400,700)) |
| L+G+80bp+G+R | c1 + 100–399; c1 + 480–779 | Duplicated |
| G and RC(G) on one contig | c1 + 100–399; c1 − 480–779 | Duplicated |
| 70 % fragment | c1 + 100–309 1.0 0.7 | Fragmented |
| 40 % fragment | c1 + 100–223 0.991935 0.413333 | Fragmented (identity ≥ 0.9) |
| absent | — | Missing |
| 1-bp deletion | c1 + 100–398 0.996667 1.0 | Single |
| complete + 40 % fragment | c1 + 100–399; c2 + 40–159 1.0 0.4 | Single |
| 3 substitutions | c1 + 100–399 0.99 1.0 | Single |
| substitution every 20 bp | — (no 28-bp exact word) | Missing |
| lower case, one N | c1 + 100–399 0.996667 1.0 | Single |
| RC(G[100..199]) | c1 − 100–202 102/103 103/300 | Fragmented |

Thresholds on the 40 % fragment: identity 1.0 → Missing; identity 1.0 + min_length 0.4 → Fragmented. 70 % fragment with coverageThreshold 0.7 → Single (inclusive).
Multi-marker (g1 single, g2 two contigs, g3 = reversed G absent, g4 48.4 % fragment): S1 D1 F1 M1, C% 50, D% 25 (old D/C = 50), notation `C:50.0%[S:25.0%,D:25.0%],F:25.0%,M:25.0%,n:4`.
BUSCO notation (Python `round(x, 1)` of `count / n * 100`): (n4 S2 D1 M1) `C:75.0%[S:50.0%,D:25.0%],F:0.0%,M:25.0%,n:4`; (n3 S1 F2) `C:33.3%[S:33.3%,D:0.0%],F:66.7%,M:0.0%,n:3`; (n2000 S199) `C:10.0%…` (199/2000·100 = 9.950000000000001); (n2000 D1) `C:0.1%[S:0.0%,D:0.1%]…` (round(0.05, 1) = 0.1).

**Random cross-check:** 400 random marker/assembly sets (seed 2026; markers 60–240 bp; exact / partial / 1–3 % mutated / reverse-complement copies, spacers 0–80 bp, 1–3 contigs; labels S111 D184 F29 M76): C# classification = oracle **400/400**. Planted-copy fuzz test (`AssemblyCompletenessFuzzTests`, 150 cases) locked.

### k-mer evaluation

| Spectrum (+ assembly) | filt | peak | T | E | error rate | size | solid | asm∩solid | completeness | TOTAL | ASM_ONLY | QV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 canonical 3-mers 1:5,2:2,3:1,4:3,5:6,6:3 + `ACGTTGCAAC`,`GGNNCATTAC` (canonical) | 3 | 5 | 72 | 12 | 0.05896397111897156 | 12 | 12 | 2 | 2/12 | 12 | 4 | 8.981858105339732 |
| same, canonical = false | 3 | 5 | 72 | 12 | same | 12 | 12 | 1 | 1/12 | 12 | 9 | 4.317519438166913 |
| same, expectedCoverage 4, no assembly | 3 | 4 | 72 | 12 | same | 15 | 12 | — | NaN | | | |
| same + `AAACC` | | | | | | | | 0 | 0 | 3 | 0 | +∞ |
| plateau 1:3,2:1,3:1,4:3 | 0 (zero slope) | 1 | 20 | 0 | 0 | 20 | 8 | | | | | |
| rising 2:1,3:2 | 2 | 3 | 8 | 2 | 0.09143970358393017 | 2 | 2 | | | | | |
| AAA4,CCC4,GGG2,TTT1 | 0 | 4 | 11 | 0 | 0 | 3 (old 2) | 4 | | | | | |
| same, expectedCoverage 2 | | 2 | | | | 6 = round(5.5) (old 5) | | | | | | |
| 100 6-mers ×1 + 1000 ×30 | 29 | 30 | 30100 | 100 | 5.544779023926213e-4 | 1000 | 1000 | | | | | |

GenomeScope `uniroot` (scipy `brentq`) vs closed form: max |Δ| 2.3e-16 over the random set.

**Random cross-check:** 400 random read sets (seed 2027; k ∈ {3,4,5,7,11}; genome 30–400 bp; 50-bp reads at 4–25×, half with a 1 % error, both strands; meryl-style canonical counting re-implemented in Python; assembly = genome halves with 2 % mutations; 333 with a valley filt > 0): every field (k, filt, peak, T, E, error rate, size, solid, asm∩solid, completeness, TOTAL, ASM_ONLY, Merqury E, QV) C# = oracle **400/400**. A naive in-test reference is locked as `AssemblyCompletenessFuzzTests` (400 cases).

---

## Findings (B10 F41–F47)

- **F41** `FindMarkerHits` was not an alignment: forward strand only (a marker on the reverse strand → Missing), coverage = assembly k-mer occurrences / distinct marker k-mers (can exceed 1; a contig with two copies gave one hit), `identity = min(1, coverage·1.1)` fabricated, hits kept at `coverage ≥ 0.5·threshold` / `identity ≥ 0.8·threshold` (unsourced). → seeded (28-bp megablast word) Smith-Waterman-Gotoh local alignment on both strands with BLAST blastn scoring via canonical `SequenceAligner.LocalAlignAffine` / `CalculateStatistics`; new public `FindMarkerHits` + `MarkerHit` (0-based inclusive coordinates).
- **F42** Classification: complete copies on one contig were never Duplicated (one hit per contig) and fragments used the relaxed magic thresholds → compleasm `pass_tids` filter + `record_1st_gene_label` incl. region clustering; new `minLengthFraction` (compleasm `min_length_percent`, 0.6).
- **F43** `DuplicationPercent` was D/C·100; BUSCO D% = D/n·100. Added `SingleCopyPercent`, `FragmentedPercent`, `MissingPercent`, `BuscoNotation` (BUSCO evaluation order and Python rounding).
- **F44** `EstimateCompletenessFromKmers` completeness = (#k-mers ≥ peak/2)/(distinct − singletons) — unsourced and computable without an assembly. → Merqury k-mer completeness with an assembly (new overload + `KmerCompletenessEstimate`, Merqury QV); without an assembly NaN (MCP `null`).
- **F45** Error rate was singletons / distinct k-mers (not a per-base rate; e.g. 0.25 for {4,4,2,1}). → GenomeScope/Merqury `1 − (1 − E/T)^(1/k)` with E = k-mers at or below Merqury's solid threshold.
- **F46** Genome size = all k-mers (errors included) / mode-excluding-count-1, integer division; peak ignored the error valley. → (T − E)/peak with the peak = mode above Merqury `filt`, R `round`.
- **F47** Validation: null inputs threw NullReferenceException; mixed k-mer lengths, duplicates and counts ≤ 0 were accepted; thresholds unchecked in the library; empty markers silently Missing → ArgumentNullException / ArgumentOutOfRangeException / ArgumentException. MCP docs/schemas rewritten (old "k-mer based" / "clean spectrum completeness 1").

## Deviations / limitations kept (documented)

- **Not BUSCO's search:** BUSCO/compleasm align lineage protein HMM profiles (HMMER score and length cut-offs, gene predictors / miniprot); nucleotide markers given by the caller are aligned directly — no introns/frameshifts, no profile score cut-off; the 28-bp exact seed is the significance filter, so markers diverged enough to lose every exact 28-mer are Missing (as megablast). The lineage datasets and HMM search are out of scope (multi-GB data, separate tools) — **LIMITATIONS proposal**.
- **Genome size without GenomeScope's model fit:** the NB-mixture (heterozygosity, duplication, 4 peaks, `nlsLM` rounds) is not fitted; the observed mode above the Merqury valley stands for 2·kcov. For heterozygous genomes whose het peak is highest, the size is ~2× too large unless `expectedCoverage` is the homozygous peak — **LIMITATIONS proposal**.
- Merqury quirks reproduced faithfully: a zero slope never marks the valley (plateau → filt 0); a histogram rising from its first line gives filt = second depth − 1.
- compleasm quirk reproduced: abutting complete copies (half-open stop = next start) merge into one cluster → Single.
- `identityThreshold` default kept at 0.9 (compleasm's protein default is 0.4).

## R1-WP1 addendum (2026-11-25): ASCII fold (A15-7, A11-6)

Sources opened: `Infrastructure/SequenceCase.cs` remarks (ASCII-only fold; B01 F21), `SequenceExtensions.IsValidDna` (IUPAC-IUB unambiguous DNA = Biopython `unambiguous_dna_letters`). `IsAcgt(span)` vs `IsValidDna(span)` differ only for lower-case input, which `FindRepetitiveRegions` has already upper-cased. Non-ACGT (incl. non-ASCII) symbols make a WindowMasker unit invalid exactly as `N` does (tests compare the result with the `N` variant). Note: `CalculateLocalQuality` passes the window to `SequenceComplexity.CalculateLinguisticComplexity`, which folds with `ToUpperInvariant` itself (B04 file, left as is; cross-batch note in the WP report).

