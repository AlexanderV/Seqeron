---
type: concept
title: "Rare codon analysis (FindRareCodons + %MinMax / Sherlocc clusters)"
tags: [annotation, algorithm]
mcp_tools:
  - find_rare_codons
sources:
  - docs/Evidence/CODON-RARE-001-Evidence.md
  - docs/algorithms/Codon_Optimization/Rare_Codon_Detection.md
  - docs/Validation/reports/CODON-RARE-001.md
source_commit: 8ce0af79a29f9fbddc217026508abc6f2c572e61
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-rare-001-evidence
      evidence: "ID: CODON-RARE-001 ... Method Under Test: CodonOptimizer.FindRareCodons(string, CodonUsageTable, double)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:codon-optimization
      source: codon-rare-001-evidence
      evidence: "AvoidRareCodons replaces only codons whose frequency falls below a threshold — the same rare-codon criterion FindRareCodons defines (default 0.15)"
      confidence: high
      status: current
---

# Rare codon analysis (FindRareCodons + %MinMax / Sherlocc clusters)

Identifies the **rare codons** in a protein-coding sequence — codons whose usage frequency in the
target organism falls below a threshold — and, in the 2026-06-24 addendum, the **clusters/runs** of
rare codons that matter most biologically. Where [[relative-synonymous-codon-usage|RSCU]],
[[codon-adaptation-index|CAI]] and [[effective-number-of-codons|ENC/Nc]] reduce a gene to a
bias *summary*, rare-codon analysis **localizes** the problematic codons and stretches;
[[codon-optimization]]'s `AvoidRareCodons` strategy is the actuator that removes exactly what this
unit detects. Validated as [[codon-rare-001-evidence|CODON-RARE-001]]; see [[test-unit-registry]]
for how the unit is tracked. The independent two-stage re-validation verdict — Stage A
PASS-WITH-NOTES / Stage B PASS / CLEAN, 45/45 tests, every %MinMax and Sherlocc oracle reproduced
against the live library — is recorded in the report [[codon-rare-001-report]].

## Why rare codons matter

Rare codons correspond to **low-abundance tRNAs**, so translation elongation slows at them. A single
rare codon is usually tolerable, but **consecutive rare codons** can stall ribosomes and inhibit
translation — Shu et al. (2006) measured a **~3-fold** inhibition from five consecutive rare CUA
(Leu) codons in E. coli, with 5′-end runs worse than internal ones. Because elongation rate shapes
**cotranslational folding**, rare-codon placement also has regulatory and folding significance. This
is why detecting **clusters**, not just isolated rare codons, is the biologically important capability.

## Base detection: `FindRareCodons`

A codon is **rare** when its frequency within its synonymous family (from the organism's
`CodonUsageTable`) is **strictly below a threshold** — default **0.15** (< 15% relative frequency in
its family). The method returns, per rare codon, its **nucleotide position** (0-indexed, a multiple
of 3), the translated **amino acid**, and the **actual frequency**.

Standard E. coli K12 rare codons (Kazusa MG1655, freq < 0.10): **AGA** 0.04, **AGG** 0.02, **CGA**
0.06 (Arg); **CUA** 0.04 (Leu); AUA 0.07 (Ile, borderline); UAG 0.07 (stop).

**Invariants.** Positions are multiples of 3 in `[0, len−3]`; codon length always 3; frequency
∈ [0, 1]; every reported codon has freq < threshold; amino acid matches the standard genetic code;
deterministic. **Edge behaviour.** Empty → empty; length not ÷3 → trailing partial codons ignored;
DNA `T` converted to `U` internally; threshold 0 → all rare except freq-0 codons; threshold 1 → all
rare; **unknown codon → frequency 0, always reported**.

**Documented base limitations** (all addressed or scoped): no codon **context**, no positional
weighting (5′ vs internal), single threshold across all amino acids, reports frequency not tRNA
availability. The addendum closes the context limitation.

## Cluster detection (addendum): two published methods

Both are **opt-in**; the per-codon default (`FindRareCodons`) is unchanged.

### %MinMax — Clarke & Clark (2008)

A sliding window scored against **per-amino-acid synonymous** statistics. For each codon of amino
acid *i* (with *n* synonymous codons): Xij = its actual frequency, Xmax,i / Xmin,i = the most / least
common synonymous codon's frequency, Xavg,i = the family mean (Σfreq / n). Per window:

    %Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i) × 100     if ΣXij > ΣXavg,i
    %Min = Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i) × 100     if ΣXij < ΣXavg,i

Each window yields **one signed value**: −100% = a window of only the rarest synonymous codons,
+100% = only the most common, 0% = at the mean. **Rare-codon clusters appear as negative (%Min)
peaks.** Default window **18 codons** (Clarke & Clark 2008; the 2018 %MinMax revision cites z = 17).
These are the same **per-family relative fractions** that
[[relative-synonymous-codon-usage|RSCU]] normalizes.

### Sherlocc — Chartier, Gaudreault & Najmanovich (2012)

A **7-codon-wide** window. A **'slow'/pause position** is any codon with usage frequency **≤
threshold** — the *same* per-codon rare criterion as `FindRareCodons` (default 0.15). A window with
**≥ 4 of 7** slow positions is reported as a **rare-codon cluster (RCC)**. So an isolated rare codon
is flagged by `FindRareCodons` but is a Sherlocc cluster **only** when ≥ 4 fall inside a 7-codon
window — the exact capability the addendum adds.

### Cluster corner cases

- **Single-codon amino acids (Met AUG, Trp UGG)** — Xmax = Xmin = Xavg = Xij, contributing 0 to both
  the %MinMax numerator and denominator, so the window value stays defined (**no divide-by-zero /
  NaN**). Same single-codon-family fact that fixes their `w ≡ 1` in [[codon-adaptation-index|CAI]].
- **Window longer than the sequence** → no windows / no clusters.
- **Overlapping qualifying windows** are merged into **one maximal cluster** so a long rare run is
  reported once (implementation choice — Sherlocc reports regions, not raw windows).

## Worked oracles

E. coli K12 Arg family (CGU 0.38, CGC 0.40, CGA 0.06, CGG 0.10, AGA 0.04, AGG 0.02; Xavg = 0.16667,
Xmax = 0.40, Xmin = 0.02):

- `AGA·AGA·AGA` (%MinMax, w 3) → **−86.36%**; `CGC·CGC·CGC` → **+100%**; `CUG·AGA` (w 2) → **+36.47%**.
- `7× AGA` (Sherlocc, ≥4) → **1 cluster** (codons 0–6, 7 rare); `4 AGA + 3 CGC` → **1 cluster**
  (exactly 4); `3 AGA + 4 CGC` → **no cluster** (3 < 4).

Base method: `AUGAGA` → AGA@3; `AUGAGAAGG` → AGA@3, AGG@6; `AUGCUG` → none; `CUGCUA` → CUA@3.

## Place in the codon-usage family

Rare-codon analysis is the family's **localization / diagnostic** view: it names *which* codons and
*which* stretches are problematic, rather than scoring the whole gene like
[[relative-synonymous-codon-usage|RSCU]] / [[codon-adaptation-index|CAI]] /
[[effective-number-of-codons|ENC/Nc]], and it feeds directly into
[[codon-optimization|codon optimization]]'s `AvoidRareCodons` strategy (which replaces exactly the
sub-threshold codons this unit flags). See [[algorithm-validation-evidence]] for the shared
evidence-artifact pattern.
