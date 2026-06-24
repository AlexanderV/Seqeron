# Validation Report: CODON-RARE-001 — Rare Codon Detection

- **Validated:** 2026-06-24   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.FindRareCodons(string codingSequence, CodonUsageTable table, double threshold = 0.15)`
  (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs:609)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia / codon usage bias** and **Kazusa Codon Usage Database** (E. coli K12,
  species 316407 / 83333) — the standard rare-codon metric is the per-amino-acid relative
  **fraction** (Kazusa "fraction" column; each synonymous family sums to ~1.0). Rare codons
  correspond to low-abundance tRNAs and slow translation.
- **Shu et al. (2006)** — CUA (Leu) is a rare E. coli codon (~0.04 fraction); consecutive
  rare Leu codons cause ~3-fold translation inhibition.
- **Sharp & Li (1987)** — directional synonymous codon usage bias; the rare criterion is
  applied **per synonymous family** (relative to synonyms), not globally.
- **Web search (2026-06)** independently corroborates the canonical "RIL" rare E. coli set —
  **AGA, AGG (Arg), AUA (Ile), CUA (Leu)**, plus CGA (Arg) — as the recognised rare codons
  supplemented by tRNAs in BL21-CodonPlus(RIL) strains.

### Metric & threshold confirmed
- **Metric:** per-amino-acid relative synonymous fraction f / Σf_family (the stored table
  values are already per-family fractions in [0,1]). This is a relative-within-family
  criterion, matching the spec text "less than 15% relative frequency within their synonymous
  codon family." Not an absolute per-1000 frequency.
- **Threshold:** default **0.15**; comparison is **strict `<`** (codon exactly at threshold
  is NOT flagged).

### Edge-case semantics check
Empty/null → empty; not divisible by 3 → trailing incomplete codon ignored; no rare → empty;
all rare → all positions; DNA T → RNA U internally; case-insensitive; unknown codon →
fraction 0 (always flagged); single-codon AAs (Met AUG=1.00, Trp UGG=1.00) never flagged
unless threshold > 1.0. All defined and sourced.

### Independent cross-check (numbers)
Embedded EColiK12 table fractions vs. literature / Kazusa:

| Codon | AA | Table | Kazusa/literature | Rare (<0.10)? |
|-------|----|-------|-------------------|---------------|
| AGA | R | 0.04 | ~0.04 | yes |
| AGG | R | 0.02 | ~0.02 | yes |
| CGA | R | 0.06 | ~0.06–0.07 | yes |
| CUA | L | 0.04 | ~0.04 | yes |
| AUA | I | 0.07 | ~0.07 | yes |
| CUG | L | 0.50 | ~0.50 | no (common) |
| AGA (Yeast) | R | 0.48 | ~0.48 | no (common) |

Hand-trace `AUGAGA` @ threshold 0.10: AUG (M, 1.00 — not rare); AGA (R, 0.04 < 0.10) → one
hit (position 3, "AGA", "R", 0.04). Matches code.

### Findings / divergences (PASS-WITH-NOTES)
- The validation prompt and the TestSpec narrative mention **rare-codon clusters/runs**. The
  implemented method performs **per-codon** detection only — it reports each individual rare
  codon position and does NOT detect or group consecutive runs. This is the method's actual
  documented contract (return tuple is per-codon); cluster detection is simply out of scope,
  not a defect. The Evidence doc already lists "does not consider codon context (consecutive
  rare codons)" as a known limitation. Flagged as a description-overreach note, not a code bug.
- The metric is the relative **fraction** f/Σf, not the relative-adaptiveness w = f/f_max.
  Both are legitimate per-family rare criteria; spec and code both use the fraction form, so
  they agree (w is used separately for CAI in CalculateCAI).
- Minor rounding differences between the embedded table and live Kazusa snapshots (e.g. CGA
  0.06 vs ~0.07) do not change any rare/common classification at tested thresholds.

## Stage B — Implementation

### Code path reviewed
`CodonOptimizer.FindRareCodons` (CodonOptimizer.cs:609-629):
1. Empty/null → `yield break`.
2. `ToUpperInvariant()`, `T`→`U`.
3. `SplitIntoCodons` (line 687) — only complete codons (`i + 2 < length`).
4. Per codon: `freq = table.CodonFrequencies.GetValueOrDefault(codon, 0)`; if
   `freq < threshold` → `yield (i*3, codon, TranslateCodon(codon), freq)`.

### Formula realised correctly?
Yes. Stored frequencies are per-family relative fractions, so `freq < threshold` is the
validated per-synonymous-family rare criterion with strict `<`. Position = `i*3` (0-based
nucleotide index of codon start). AA via `StandardGeneticCode` (unknown → "X", fraction 0).

### Cross-verification table recomputed vs code (via tests, 21/21 pass)
| Input | Threshold | Expected | Verified |
|-------|-----------|----------|----------|
| "" | 0.15 | empty | ✓ M01 |
| AUGAGA | 0.10 | [(3,AGA,R,0.04)] | ✓ M02/M09 |
| AUGAGAAGGCGA | 0.10 | AGA,AGG,CGA | ✓ M03 |
| AUGUAG | 0.10 | UAG (0.07) detected | ✓ M05 |
| AUGCUA | 0.04 | CUA NOT detected (strict <) | ✓ M06 |
| AUGCUGUGG | 0.15 | empty | ✓ M07 |
| ATGAGA (DNA) | 0.10 | AGA after T→U | ✓ M10 |
| AUGAUA | default 0.15 | AUA (0.07) detected | ✓ S01 |
| AUGAGA Yeast | 0.10 | AGA NOT rare (0.48) | ✓ C03 |

### Variant/delegate consistency
Single canonical method; the 0.15 default overload and explicit-threshold form are the same
method. No `*Fast` variant. Consistent.

### Numerical robustness
Frequencies are dictionary lookups in [0,1]; no division, overflow, or NaN risk. Threshold 0
and 1.0 behave as defined (C01/C02).

### Test quality audit
21 tests assert exact sourced values (positions, codons, AAs, frequencies to 0.001), exact
counts, strict-`<` boundary, organism-specific classification, and invariants
(positions multiple of 3, all freq < threshold, determinism). Deterministic, real assertions
— not no-throw tautologies. Covers all Stage-A edge cases.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — metric (per-family fraction), threshold (0.15, strict `<`),
  and rare E. coli codon set independently confirmed against Kazusa and the rare-codon
  literature.
- **Stage B: PASS** — code faithfully realises the per-family `freq < threshold` criterion;
  positions and frequencies match; 21/21 per-codon unit tests green.
- **State: CLEAN** (for the per-codon contract).

## 2026-06-24 — Limitation closed: rare-codon cluster / run detection added

The PASS-WITH-NOTES item "cluster/run detection mentioned in the spec narrative is NOT
implemented" has been **resolved** by adding two opt-in, source-backed sliding-window methods
to `CodonOptimizer` (per-codon `FindRareCodons` default behaviour unchanged):

- `CalculateMinMaxProfile(string, CodonUsageTable, int windowSize = 18)` — the **%MinMax**
  profile of **Clarke & Clark (2008), "Rare Codons Cluster", PLoS ONE 3(10):e3412**
  (doi:10.1371/journal.pone.0003412). Per amino acid: `Xij`, `Xmax,i`, `Xmin,i`,
  `Xavg,i = (1/n)ΣXij`; per window `%Max = Σ(Xij−Xavg)/Σ(Xmax−Xavg)×100` (positive) or
  `%Min = Σ(Xavg−Xij)/Σ(Xavg−Xmin)×100` (returned negative). Default window 18 codons.
- `FindRareCodonClusters(string, CodonUsageTable, double rareThreshold = 0.15, int windowSize = 7,
  int minRareCodons = 4)` — the **Sherlocc** rare-codon-cluster rule of **Chartier, Gaudreault &
  Najmanovich (2012), Bioinformatics 28(11):1438–1445** (doi:10.1093/bioinformatics/bts149):
  "a seven position-wide window … containing at least four pause positions out of seven";
  a "pause"/"slow" position is a codon below the rare threshold. Overlapping qualifying windows
  merge into maximal clusters.

Both windows/thresholds are copied verbatim from the retrieved sources (18 codons; 7-codon
window, ≥4 rare) — no invented parameters. New tests:
`tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_RareCodonClusters_Tests.cs` (24 cases,
exact %MinMax values and exact cluster boundaries). Evidence/algorithm-doc/TestSpec extended.

**State: re-validation pending.** The ROOT `ALGORITHMS_CHECKLIST_V2.md` Status for
CODON-RARE-001 is reset ☑ → ☐ pending an independent re-validation of the new capability.
