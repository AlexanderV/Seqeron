# Validation Report: CODON-RARE-001 — Rare Codon Detection

- **Validated:** 2026-06-25 (fresh re-validation; see bottom section)   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.FindRareCodons` (per-codon, :663),
  `CodonOptimizer.CalculateMinMaxProfile` (%MinMax, :720),
  `CodonOptimizer.FindRareCodonClusters` (Sherlocc RCC, :825)
  (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **State:** ✅ CLEAN

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

## 2026-06-25 — FRESH RE-VALIDATION (campaign E2: +%MinMax +Sherlocc clusters)

Re-validated from authoritative EXTERNAL first sources retrieved this session, ignoring the
repo's own TestSpec/Evidence. The unit now spans **three** public methods:
`FindRareCodons` (per-codon), `CalculateMinMaxProfile` (%MinMax), `FindRareCodonClusters`
(Sherlocc RCC). Source `CodonOptimizer.cs:663` / `:720` / `:825`. Tests:
`CodonOptimizer_FindRareCodons_Tests.cs` (21) + `CodonOptimizer_RareCodonClusters_Tests.cs` (24).

- **Stage A: PASS-WITH-NOTES**
- **Stage B: PASS**
- **State: ✅ CLEAN**

### Sources opened this session (verbatim confirmation)
- **Clarke TF, Clark PL (2008) "Rare Codons Cluster", PLoS ONE 3(10):e3412** (PLOS full text).
  Confirmed verbatim: 18-codon sliding window; `Xij` = actual codon usage frequency, `Xmax,i`
  / `Xmin,i` = most/least common synonymous codon frequency, `Xavg,i` = sum of synonymous
  frequencies ÷ number of synonymous codons; the metric is "the difference between Xij and
  Xavg,i divided by the difference between Xmax,i (or Xmin,i) and Xavg,i". Sign convention
  verbatim: "If the codon usage frequencies for a given window are greater than the average, a
  value will be returned for %Max; if less than the average, … %Min. … %Min values are plotted
  and reported as negative numbers." Range −100 … +100; 0 = average usage.
- **Chartier M, Gaudreault F, Najmanovich R (2012) Bioinformatics 28(11):1438–1445,
  doi:10.1093/bioinformatics/bts149** (PMC3465090 full text). Confirmed verbatim cluster rule:
  "a … seven position-wide window … searching for windows with at least four pause positions
  out of seven." A position is "slow"/"pause" when its codon-usage average falls below a
  threshold (paper fits an extreme-value distribution; thresholds 13–18 tested).

### Formula check — matches source exactly (CodonOptimizer.cs:771-799)
`sumXij`, `sumXavg`, `sumMaxDelta = Σ(Xmax−Xavg)`, `sumMinDelta = Σ(Xavg−Xmin)`;
`%Max = (ΣXij−ΣXavg)/sumMaxDelta×100` when ΣXij>ΣXavg, `%Min = −(ΣXavg−ΣXij)/sumMinDelta×100`
when ΣXij<ΣXavg, else 0. `Xavg` is the per-family arithmetic mean (sum÷count). Single-codon
AAs / unknown / stop contribute 0 to both numerator and denominator (no NaN). Window slides
one codon; produces `codonCount − windowSize + 1` windows. ✓ faithful to Clarke & Clark.
Cluster code (`:825-901`): 7-codon window, `windowRare ≥ minRareCodons` (default 4); a codon
is rare when `freq < rareThreshold` (default 0.15, strict `<`, same criterion as
`FindRareCodons`); overlapping/adjacent qualifying windows merge into maximal clusters. ✓ the
verbatim Sherlocc window/count rule.

### Independent cross-check (hand-computed THIS session, then confirmed vs LIVE library)
E. coli K12 Arg family (Kazusa 316407): CGU 0.38, CGC 0.40, CGA 0.06, CGG 0.10, AGA 0.04,
AGG 0.02 → Σ=1.00, Xavg=0.166667, Xmax=0.40, Xmin=0.02. Leu family: Xavg=0.166667, Xmax=0.50.

| Window | Hand-computed (Clarke&Clark formula) | Live `CalculateMinMaxProfile` |
|--------|--------------------------------------|-------------------------------|
| 3× CGC (most common Arg), w=3 — **%Max>0** | (0.40−0.16667)/(0.40−0.16667)×100 = **+100** | +100.0000 ✓ |
| 3× AGA (rarest Arg), w=3 — **%Min<0** | −(0.16667−0.04)/(0.16667−0.02)×100 = **−86.36364** | −86.36364 ✓ |
| CUG(0.50)+AGA(0.04), w=2 — **%Max>0** | (0.54−0.33333)/((0.50−0.16667)+(0.40−0.16667))×100 = **+36.47059** | +36.47059 ✓ |

Cluster hand-check: sequence = 10× common (CGC) + 7× rare (AGA). The 7 rare codons sit at
codon indices 10–16. `FindRareCodons` flagged all 7 at nt positions 30,33,36,39,42,45,48 ✓.
`FindRareCodonClusters` returned ONE cluster `codons 7..16, RareCount=7`: a 7-wide window at
start 7 (covers 7–13) already holds codons 10–13 = 4 rare ≥ 4, so the qualifying-window union
legitimately begins at codon 7; only the 7 true rare codons are counted. ✓ Sherlocc-correct,
and isolated rare codons (C4) correctly yield NO cluster while per-codon detection still flags
them — the exact capability the campaign added.

### Edge-case semantics (all defined & sourced)
empty/null → empty; non-multiple-of-3 → trailing partial codon dropped (`SplitIntoCodons`,
`i+2 < length`); sequence shorter than window → empty; window<1 / minRare<1 → throws; %Max and
%Min sign branches both exercised; threshold boundary (3 rare = no cluster, 4 rare = cluster,
strict `<` for %MinMax rare); single-codon AA → 0 not NaN; DNA T→U normalised; organism-specific
(AGA rare in E. coli, common 0.48 in yeast).

### Test quality audit (HARD gate — PASS)
45 tests total across the two fixtures. Expected values trace to Clarke&Clark/Chartier and the
hand-math above (e.g. −86.36363636363637, +100, +36.470588235294116, cluster boundaries 0..6 /
0..9 / 7..16), NOT code echoes. Real assertions (exact positions, codons, AAs, frequencies to
1e-10/1e-12, counts, bounds invariant [−100,100], determinism) — no no-throw tautologies, no
green-washing. Every public method/overload and every Stage-A path covered.

### Findings / divergences (PASS-WITH-NOTES, not defects)
1. **Sherlocc "slow/pause" criterion is simplified.** The paper derives the per-position rare
   cutoff statistically (extreme-value fit over a 7-codon usage average); the implementation
   uses a simpler per-codon fraction threshold (default 0.15) — the SAME criterion as
   `FindRareCodons`. The published **window=7, ≥4-of-7** rule is reproduced verbatim and the
   per-position criterion is documented and fully parameterizable (`rareThreshold`). Documented
   simplification of an intentionally lighter-weight, table-driven variant; not a code defect.
2. The %MinMax metric uses the per-family relative **fraction** (Kazusa "fraction" column),
   consistent with Clarke & Clark's `Xij` and with `FindRareCodons`; `w = f/f_max` is used
   separately for CAI. Agree.

### Verdict
Stage A **PASS-WITH-NOTES** (formulas, windows, thresholds, sign convention all confirmed
verbatim against the two primary papers; one documented simplification of the Sherlocc
per-position criterion). Stage B **PASS** (code realises both formulas exactly; live library
reproduces every hand-computed number; 45/45 tests green; full `dotnet test Seqeron.sln`
= 0 failed, Genomics 18787 passed). **State: ✅ CLEAN.**
