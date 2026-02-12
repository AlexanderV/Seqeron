# Disorder Prediction

## Documented Theory

### Purpose

Intrinsically Disordered Protein (IDP) prediction identifies regions of protein sequences that lack stable 3D structure under physiological conditions. These regions exist as dynamic conformational ensembles and play critical roles in cell signaling, regulation, and molecular recognition (Dunker et al. 2001, Wikipedia "Intrinsically disordered proteins").

### Core Mechanism

The prediction uses the TOP-IDP scale — an optimized amino acid disorder propensity scale from Campen et al. (2008), PMC2676888. The TOP-IDP scale was derived by surveying 517 amino acid scales and applying simulated annealing to maximize discrimination between ordered and disordered regions (ARV = 0.761).

For each residue, the score is the average of normalized TOP-IDP values over a sliding window (default 21 residues):

$$S_i = \frac{1}{|W_i|} \sum_{c \in W_i} \frac{\text{TOP-IDP}(c) - \text{TOP-IDP}_{\min}}{\text{TOP-IDP}_{\max} - \text{TOP-IDP}_{\min}}$$

where $\text{TOP-IDP}_{\min} = -0.884$ (Trp) and $\text{TOP-IDP}_{\max} = 0.987$ (Pro).

Residues with score ≥ 0.542 (maximum-likelihood cutoff from Campen et al. 2008) are classified as disordered.

### Properties

- **Deterministic:** Yes — same input produces same output.
- **Score range:** All scores are in [0, 1] (inherent from normalized TOP-IDP averaging).
- **Case-insensitive:** Input is uppercased before processing.
- **Unknown residues:** Amino acids not in the standard 20 are silently ignored (zero contribution to that component).

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n × w) where w = window size | Implementation analysis |
| Space | O(n) | One ResiduePrediction per residue |

---

## Implementation Notes

**Implementation location:** [DisorderPredictor.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `PredictDisorder(sequence, windowSize=21, disorderThreshold=0.542, minRegionLength=5)`: Main entry point. Returns `DisorderPredictionResult` containing per-residue predictions, identified disordered regions, overall disorder content, and mean score.
- `CalculateDisorderScore(window)`: Private method computing the average normalized TOP-IDP value for a sequence window.
- `CalculatePerResidueScores(sequence, windowSize, threshold)`: Private helper iterating sliding window across sequence.
- `GetDisorderPropensity(char)`: Public utility returning TOP-IDP propensity for a single amino acid.
- `IsDisorderPromoting(char)`: Public utility checking Dunker et al. (2001) classification. Returns true for {A, R, G, Q, S, P, E, K}.
- `DisorderPromotingAminoAcids`: Public property returning 8 disorder-promoting amino acids — Dunker et al. (2001).
- `OrderPromotingAminoAcids`: Public property returning 8 order-promoting amino acids — Dunker et al. (2001).
- `AmbiguousAminoAcids`: Public property returning 4 ambiguous amino acids {D, H, M, T} — Dunker et al. (2001).
- `CalculateHydropathy(string)`: Public utility returning mean Kyte-Doolittle hydropathy for a sequence.

---

## Deviations and Assumptions

None. All parameters sourced from peer-reviewed publications:
- Scoring: TOP-IDP scale averaged over window — Campen et al. (2008) Table 2.
- Threshold: 0.542 — Campen et al. (2008) maximum-likelihood cutoff.
- Classification: Dunker et al. (2001) disorder/order/ambiguous sets.
- Hydropathy: Kyte & Doolittle (1982).
- Window size: 21 residues — Campen et al. (2008).

---

## Sources

- Kyte J, Doolittle RF (1982). J Mol Biol. 157(1):105-132.
- Dunker AK et al. (2001). J Mol Graph Model. 19(1):26-59.
- Uversky VN et al. (2000). Proteins. 41(3):415-427.
- Campen A et al. (2008). Protein Pept Lett. 15(9):956-963.
- Shannon CE (1948). Bell System Technical Journal. 27(3):379-423.
- Wikipedia: "Intrinsically disordered proteins", "Hydrophilicity plot", "Amino acid"
