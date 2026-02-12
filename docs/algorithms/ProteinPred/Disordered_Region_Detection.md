# Disordered Region Detection

## Documented Theory

### Purpose

Identifies contiguous intrinsically disordered regions (IDRs) within a protein sequence from per-residue disorder predictions, and classifies each region by amino acid composition bias.

IDRs are protein segments that lack fixed three-dimensional structure under physiological conditions. They are enriched in polar, charged amino acids and depleted in bulky hydrophobic residues (Dunker et al. 2001). Identifying and classifying these regions is essential for understanding protein function, as IDRs mediate molecular recognition, flexible linking, and signaling (van der Lee et al. 2014).

### Core Mechanism

**Region Identification:**

Given a list of per-residue disorder predictions (each marked as disordered or ordered based on a threshold), the algorithm performs a linear scan to identify maximal contiguous runs of disordered residues. A run is emitted as a region only if its length ≥ a minimum length parameter.

```
for each residue i:
    if residue[i].IsDisordered:
        extend current region
    else:
        if current region length ≥ minLength:
            emit region(start, end, meanScore, classification)
        reset current region
// handle trailing region at sequence end
```

For each emitted region:
- `MeanScore` = average of constituent residue DisorderScores
- `Confidence` = (meanScore − cutoff) / (1.0 − cutoff) — normalized distance from TOP-IDP decision boundary (Campen 2008)
- `RegionType` = classification by amino acid composition

**Region Classification:**

Each disordered region is classified by the dominant amino acid composition:
1. Proline-rich: fraction of P > 0.25 — P has highest TOP-IDP propensity (Campen 2008); 0.25 = 5× random frequency
2. Acidic: fraction of (E + D) > 0.25 — Das & Pappu (2013) f− boundary in diagram-of-states
3. Basic: fraction of (K + R) > 0.25 — Das & Pappu (2013) f+ boundary in diagram-of-states
4. Ser/Thr-rich: fraction of (S + T) > 0.25 — van der Lee (2014) Table 1, P/T/S tandem repeats
5. Long IDR: length > 30 and no composition bias above 0.25 — Ward (2004), van der Lee (2014)
6. Standard IDR: fallback (short, no composition bias)

Priority order follows most-specific single-AA bias first (Proline, Campen 2008), then
charge-based classes (Das & Pappu 2013), then S/T repeats (van der Lee 2014 Table 1),
then length-based classification (Ward 2004).

### Properties

- **Deterministic:** Same input always produces same output.
- **Linear scan:** Single pass over the prediction list.
- **Non-overlapping:** Emitted regions are non-overlapping and ordered by position.
- **Coverage invariant:** Every disordered residue is either in an emitted region or in a below-minLength run that was filtered out.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n) | Single pass over n residue predictions |
| Space | O(n) | Output regions proportional to input length |

---

## Implementation Notes

**Implementation location:** [DisorderPredictor.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `IdentifyDisorderedRegions(predictions, threshold, minLength)`: Private static method. Scans the per-residue prediction list to find contiguous disordered residues. Emits `DisorderedRegion` records for runs ≥ minLength. Called internally by `PredictDisorder()`.
- `ClassifyDisorderedRegion(region)`: Private static method. Examines amino acid composition of a region and returns a classification string. Priority: Proline-rich > Acidic > Basic > Ser/Thr-rich > Long IDR > Standard IDR.
- `CalculateConfidence(meanScore, length)`: Private static method. Computes confidence as normalized distance from TOP-IDP cutoff: (meanScore − 0.542) / (1.0 − 0.542), clamped to [0, 1]. Derived from Campen et al. (2008) scale.

Both canonical methods are private and tested indirectly through `PredictDisorder()`, which is the public API.

---

## Deviations and Assumptions

None. All parameters are traceable to peer-reviewed sources:

| # | Parameter | Value | Source |
|---|-----------|-------|--------|
| 1 | Classification threshold | 0.25 | ≥5× random single-AA frequency (1/20); for charged pairs matches f+/f− boundary in Das & Pappu (2013) |
| 2 | Classification priority | Pro > Acidic > Basic > S/T > Long > Standard | Most-specific bias first: Campen (2008), Das & Pappu (2013), van der Lee (2014) Table 1, Ward (2004) |
| 3 | Long IDR threshold | >30 residues | Ward et al. (2004); van der Lee et al. (2014) |
| 4 | Confidence formula | (meanScore − cutoff) / (1.0 − cutoff) | Normalized TOP-IDP distance, Campen et al. (2008) |

---

## Sources

- Campen A, et al. (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett. 15(9):956-963. https://doi.org/10.2174/092986608785849164
- Das RK, Pappu RV (2013). "Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues." Proc Natl Acad Sci USA. 110(33):13392-13397. https://doi.org/10.1073/pnas.1304749110
- Dunker AK, et al. (2001). "Intrinsically disordered protein." J Mol Graph Model. 19(1):26-59. https://doi.org/10.1016/s1093-3263(00)00138-8
- van der Lee R, et al. (2014). "Classification of intrinsically disordered regions and proteins." Chemical Reviews. 114(13):6589-6631. https://doi.org/10.1021/cr400525m
- Ward JJ, et al. (2004). "Prediction and functional analysis of native disorder in proteins from the three kingdoms of life." J Mol Biol. 337(3):635-645. https://doi.org/10.1016/j.jmb.2004.02.002
