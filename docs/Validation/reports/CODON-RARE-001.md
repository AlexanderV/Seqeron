# Validation Report: CODON-RARE-001 — Rare Codon Detection

- **Validated:** 2026-06-12   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.FindRareCodons(string codingSequence, CodonUsageTable table, double threshold = 0.15)`
  (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs:609)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia / literature on codon usage bias** — rare codons correspond to low-abundance
  tRNAs; slow translation elongation; consecutive rare codons inhibit translation. Confirmed.
- **Kazusa Codon Usage Database** (E. coli, species 316407 / 83333) — the per-amino-acid
  relative **fraction** column is the standard rare-codon metric. For E. coli the arginine
  codons AGA, AGG, CGA, leucine CUA and isoleucine AUA have the lowest fractions in their
  synonymous families. Confirmed against the live Kazusa table.
- **Shu et al. (2006)** — CUA (Leu) is a rare E. coli codon (~0.04 fraction); consecutive rare
  Leu codons cause ~3-fold translation inhibition. Confirmed.
- **Sharp & Li (1987)** — directional synonymous codon usage bias; rare codons have low
  relative adaptiveness. Confirms that the rare-codon criterion is applied **per synonymous
  family** (relative to synonyms), not globally.
- **Microbial Cell Factories (2009) / RIL-codon literature** — confirms the canonical "classic"
  rare E. coli codons: **AGA, AGG (Arg), AUA/ATA (Ile), CUA/CTA (Leu)** — the "RIL" codons
  recognised by the argU/ileY/leuW tRNAs supplied in BL21-CodonPlus(RIL).

### Metric & threshold confirmed
- **Metric:** per-amino-acid **relative synonymous fraction** f / Σf_family (the Kazusa
  "fraction" column; each synonymous family sums to ~1.0). This is a relative-within-family
  criterion, **not** an absolute per-1000 frequency. The spec text ("less than 15% relative
  frequency within their synonymous codon family") matches this exactly.
- **Threshold:** default **0.15**; comparison is **strict `<`** (a codon at the threshold is
  NOT flagged). Per-family application is achieved because the stored table values are already
  per-family fractions.

### Edge-case semantics check
Empty → empty; sequence not divisible by 3 → trailing incomplete codon ignored; no rare
codons → empty; all rare → all positions; DNA T → RNA U internally; case-insensitive;
unknown codon → fraction 0 (always flagged); single-codon AAs (Met AUG=1.00, Trp UGG=1.00)
never flagged unless threshold > 1.0. All defined and sourced.

### Independent cross-check (numbers)
Kazusa E. coli fractions used in the table vs. live Kazusa lookup:
| Codon | AA | Spec/table | Live Kazusa (316407) | Rare (<0.10)? |
|-------|----|-----------|----------------------|---------------|
| AGA | R | 0.04 | ~0.04 | yes |
| AGG | R | 0.02 | ~0.02 | yes |
| CGA | R | 0.06 | ~0.07 | yes |
| CUA | L | 0.04 | ~0.04 | yes |
| AUA | I | 0.07 | ~0.07 | yes |
| CUG | L | 0.50 | ~0.50 | no (common) |

Known rare codon AGA (Arg) is correctly identified as rare. Worked example `AUGAGA`:
AUG (M, fraction 1.00 — not rare) and AGA (R, fraction 0.04 < 0.10) → one hit at nucleotide
position 3, codon "AGA", AA "R", frequency 0.04.

### Findings / divergences (PASS-WITH-NOTES)
- Minor rounding divergences between the embedded table and the current live Kazusa snapshot
  (CGA 0.06 vs ~0.07; UAG 0.07 vs ~0.05). These do not change any rare/common classification
  at the tested thresholds and are within normal Kazusa-version variation. Documented, no defect.
- The metric is the relative **fraction** f/Σf, not the relative-adaptiveness w = f/f_max. Both
  are legitimate per-family rare-codon criteria; the spec describes the fraction form and the
  code implements the fraction form, so they agree. (w is used elsewhere for CAI.)

## Stage B — Implementation

### Code path reviewed
`CodonOptimizer.FindRareCodons` (CodonOptimizer.cs:609-629):
1. Empty/null → `yield break`.
2. Uppercase, `T`→`U`.
3. `SplitIntoCodons` (CodonOptimizer.cs:687) — only complete codons (`i + 2 < length`).
4. For each codon: `freq = table.CodonFrequencies.GetValueOrDefault(codon, 0)`; if
   `freq < threshold` → `yield (i*3, codon, TranslateCodon(codon), freq)`.

### Formula realised correctly?
Yes. The stored frequencies are per-family relative fractions, so `freq < threshold` is the
validated per-synonymous-family rare criterion with strict `<`. Position = `i*3` (0-based
nucleotide index of the codon start). AA via the standard genetic code; unknown codon → "X"
and fraction 0 (always flagged), matching documented behaviour.

### Cross-verification table recomputed vs code (via tests)
- M02 `AUGAGA`, thr 0.10 → {AGA @ pos 3}. ✓
- M03 `AUGAGAAGGCGA`, thr 0.10 → {AGA, AGG, CGA}. ✓
- M06 `AUGCUA`, thr 0.04 → CUA NOT flagged (strict `<`). ✓
- S01 default thr 0.15: `AUGAUA` → AUA (0.07) flagged. ✓
- C03 `AUGAGA`: E. coli flags AGA (0.04); Yeast does NOT (0.48). ✓ (per-organism, per-family.)

### Variant/delegate consistency
Single canonical method; the default-threshold overload is the same method (0.15 default).
Same relative-fraction tables feed `OptimizeSequence`/`AvoidRareCodeons` and CAI consistently.

### Test quality audit
`CodonOptimizer_FindRareCodons_Tests.cs` — 32 tests covering MUST/SHOULD/COULD + invariants:
empty, single/multiple rare, nucleotide position, strict-`<` boundary, default 0.15, incomplete
codon, case-insensitivity, all-rare, mixed, threshold 0/1, cross-organism, position-multiple-of-3,
frequency-below-threshold, determinism. Assertions check exact sourced codons/positions/
frequencies (e.g. AGA freq 0.04). Real, deterministic, edge-covering.

### Findings / defects
None. Implementation faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — biology/metric/threshold confirmed against Kazusa, Shu et al.,
  Sharp & Li, and RIL-codon literature; only minor Kazusa-version rounding notes.
- **Stage B: PASS** — code matches the validated per-family relative-fraction criterion with
  strict `<`, correct positions/AA/frequency; worked examples recomputed and confirmed.
- **State: CLEAN** — no defect; no code changes required.
- **Tests:** `~RareCodon` filter = 32 passed; full `Seqeron.Genomics.Tests` = 4484 passed,
  0 failed (baseline preserved).
