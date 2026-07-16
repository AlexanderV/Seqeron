---
type: concept
title: "MoRF prediction (dip-in-disorder heuristic)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/DISORDER-MORF-001-Evidence.md
  - docs/Validation/reports/DISORDER-MORF-001.md
  - docs/algorithms/ProteinPred/MoRF_Prediction.md
source_commit: 384cf12d39a790299dd1ba7a83c1528c59e4def3
created: 2026-07-09
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: disorder-morf-001-evidence
      evidence: "Test Unit ID: DISORDER-MORF-001 ... Algorithm: MoRF (Molecular Recognition Feature) Prediction"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:protein-low-complexity-seg
      source: disorder-morf-001-evidence
      evidence: "MoRF prediction is the second unit of the protein disorder / features family after SEG low-complexity (DISORDER-LC-001); both operate on protein compositional/disorder signals but detect different features"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:intrinsic-disorder-prediction-top-idp
      source: disorder-morf-001-evidence
      evidence: "The per-residue disorder profile comes from the repository's PredictDisorder, which uses the TOP-IDP amino-acid scale normalized to [0,1]; MoRF prediction reads that profile to find an ordered dip inside disorder"
      confidence: high
      status: current
---

# MoRF prediction (dip-in-disorder heuristic)

A **Molecular Recognition Feature (MoRF)** is a short, loosely structured protein segment embedded
**within a longer intrinsically disordered region** that **undergoes a disorder-to-order transition
when it binds a partner** — disordered in isolation, folded in the bound complex. Seqeron predicts
them with the **"dip within disorder" heuristic** (Oldfield et al. 2005 / Cheng, Oldfield et al.),
validated under test unit **DISORDER-MORF-001**; the pre-implementation evidence record is
[[disorder-morf-001-evidence]] and the independent two-stage re-validation verdict (Stage A
PASS-WITH-NOTES / Stage B PASS / CLEAN, suite 6609/0) is [[disorder-morf-001-report]] — one row of
the [[validation-ledger]]. [[test-unit-registry]] tracks the unit. See
[[algorithm-validation-evidence]] for the artifact pattern.

This is the **second ingested unit of the protein disorder / features family** (DISORDER-LC / MORF /
PRED / PROPENSITY / REGION), after the [[protein-low-complexity-seg|SEG low-complexity]] anchor. It
is a **distinct algorithm**: SEG partitions a protein by *compositional* complexity, whereas MoRF
prediction reads a *per-residue disorder profile* and finds an ordered dip inside disorder. The two
are siblings, not the same operation — low-complexity regions overlap with but are not identical to
intrinsically disordered regions.

## The dip-in-disorder criterion

The heuristic "identifies short regions of order within longer regions of disorder – or **'dips'** –
in disorder prediction profiles" (Cheng/Oldfield, PMC2570644). Seqeron reports a MoRF where **all**
of the following hold:

1. an **ordered run** — a contiguous stretch of residues whose per-residue disorder score is
   **below 0.5** (the PMC2570644 disorder/order threshold),
2. of **total length within the 10–70 residue band** (Mohan et al. 2006 MoRF length range; α-MoRFs
   are "around 20 residues", candidates ≤ 30), and
3. **flanked on both sides** by ≥1 disordered residue (score ≥ 0.5) inside a predicted disordered
   region — the "within a longer region of disorder" requirement.

A dip at the very **start or end** of the sequence is *not* flanked by disorder on both sides and is
therefore **not** a MoRF.

## Underlying disorder score (TOP-IDP)

The per-residue disorder profile comes from the repository's `PredictDisorder` —
[[intrinsic-disorder-prediction-top-idp|intrinsic-disorder prediction]] (DISORDER-PRED-001), the
shared anchor of the disorder family — which uses the **TOP-IDP amino-acid scale** (Campen et al.
2008) normalized to `[0,1]`:

```
score(residue) = (TOP-IDP_raw + 0.884) / 1.871      # Campen Table 2; higher = more disordered
```

so a residue is **ordered** when `score < 0.5` and **disordered** when `score ≥ 0.5`. Representative
normalized values:

| Residue | TOP-IDP raw | Normalized | Class at 0.5 |
|---------|-------------|------------|--------------|
| P | 0.987 | 1.000 | disordered |
| E | 0.736 | 0.866 | disordered |
| L | −0.326 | 0.298 | ordered |
| I | −0.486 | 0.213 | ordered |
| W | −0.884 | 0.000 | ordered |

`PredictDisorder` applies **window averaging**, which smooths the profile near segment boundaries;
tests use flanks long enough that interior residues reach the pure per-residue score.

## Structural sub-types (Mohan 2006)

MoRFs are classified by the secondary structure they adopt **in the bound state**:

- **α-MoRF** — forms an α-helix (the class the Cheng/Oldfield "dip" mining targeted).
- **β-MoRF** — forms a β-strand.
- **ι-MoRF** — irregular secondary structure.

## MoRF score

Each reported MoRF carries a **score in `[0,1]` that increases with dip depth** (how far the ordered
run drops below the 0.5 threshold). The bounded normalization is a documented derivation (the
algorithm doc), not a source-prescribed value; the 0.5 threshold it is measured against is
source-backed (PMC2570644).

## Implementation

`DisorderPredictor.PredictMoRFs(string sequence, int minLength = 10, int maxLength = 70)` (in
`DisorderPredictor.cs`) scans the profile that `PredictDisorder` **already produced** — it does not
re-window — for maximal ordered runs, applies the length band and the both-sided flank test, and
yields each survivor as `(Start, End, Score)` with **0-based inclusive** coordinates,
**non-overlapping and ordered by `Start`**. Cost is **O(n·w)** time / **O(n)** space (`n` = length,
`w` = 21-residue disorder window); the dip scan itself is O(n). No substring/pattern search is
involved, so the repository suffix tree is **not applicable** to this unit. The score is the mean dip
depth below the threshold, normalized by the maximum possible depth:

```
Score = clamp01( (0.5 − mean_{i∈[s,e]} d(i)) / 0.5 )
```

Because every `d(i) < 0.5` inside the dip, `Score ∈ (0, 1]` and rises monotonically with dip depth.

**Two distinct thresholds — do not conflate.** The MoRF dip is defined against the **0.5**
order/disorder threshold from the MoRF literature (PMC2570644), which is **independent of the 0.542
TOP-IDP decision cutoff** that `PredictDisorder` uses for general IDR calling; the two serve
different purposes.

**Worked oracle.** `new string('P',25) + new string('L',30) + new string('P',25)` → exactly one MoRF
at `Start 29, End 50` (length 22, in band); mean disorder over the dip ≈ 0.362033 → Score ≈ 0.2759.
Replacing the ordered `L` core with the more order-promoting `I` (TOP-IDP raw −0.486) deepens the dip
(mean ≈ 0.300196, Score ≈ 0.399608), pinning score monotonicity (INV-05).

## Canonical oracle and corner cases

- **Dip-in-disorder** — a short ordered homopolymer window (e.g. `L`, score 0.298) inside long P/E
  disordered flanks → **exactly one MoRF** at the dip's coordinates.
- **Fully ordered protein** → no surrounding disorder → **no MoRF**.
- **Fully disordered protein** → no ordered dip → **no MoRF**.
- **Ordered run shorter than 10 or longer than 70 residues** → not a MoRF (Mohan length band).
- **Terminal dip** (not flanked by disorder on both sides) → not a MoRF.
- **Two separate dips** → two non-overlapping MoRFs. Case-insensitive input; null/empty → empty.

## Deviation / assumption

**One assumption**, scoped to the flank-length detail only. Oldfield et al. 2005 defines the exact
numeric dip parameters (precise flank lengths, the ordered-run window) but its **Methods section is
paywalled and could not be retrieved**. Seqeron therefore implements the fully-retrievable
**qualitative** criterion above. The **0.5 threshold**, the **10–70 residue band**, and the
**order-within-disorder shape** are all source-traceable and are **not** assumptions — only the
flank-length detail is a correctness-affecting modeling choice.

## References

Mohan A. et al. (2006) *J Mol Biol* 362(5):1043–1059 (PMID 16935303); Oldfield C.J. et al. (2005)
*Biochemistry* 44(37):12454–12470 (PMID 16156658); Cheng Y., Oldfield C.J. et al. *Biochemistry*
(PMC2570644, α-MoRF "dip" mining); Campen A. et al. (2008) *Protein Pept Lett* 15(9):956–963
(PMC2676888, TOP-IDP scale used by `PredictDisorder`); Wikipedia "Molecular recognition feature".
Full citations in [[disorder-morf-001-evidence]] (do not duplicate here). A
[[research-grade-limitations|research-grade]] implementation of the standard dip-in-disorder MoRF
heuristic.
