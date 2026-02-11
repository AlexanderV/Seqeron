# Acceptor Site Detection

## Documented Theory

### Purpose

The 3' splice site (acceptor site) marks the end of an intron in pre-mRNA. Detection of acceptor sites is essential for predicting intron boundaries and gene structure. The acceptor site is characterized by an almost invariant AG dinucleotide preceded by a polypyrimidine tract (PPT) and the consensus motif (Y)nNCAG|G (Shapiro & Senapathy 1987; Burge et al. 1999).

### Core Mechanism

1. **AG dinucleotide recognition:** The 3' end of the intron terminates with AG (RNA) or AG (DNA). This is nearly invariant across eukaryotic introns (>99%).

2. **Polypyrimidine tract (PPT):** A pyrimidine-rich region (C and U) typically 15–20 nt long, located 5–40 nt upstream of the 3' splice site. PPT quality (fraction of pyrimidines) is a major determinant of splice site strength. U2AF65 binds the PPT during spliceosome assembly.

3. **Position Weight Matrix (PWM) scoring:** Nucleotide frequencies at each position around the acceptor site are compiled into a weight matrix. Log-odds scoring quantifies how well a candidate sequence matches the consensus versus background. Key conserved positions:
   - Position -2: A (100% conserved)
   - Position -1: G (100% conserved)
   - Position -3: C (~70%)
   - Positions -5 to -15: pyrimidine-enriched (PPT)
   - Position 0 (first exonic): G (~50%)

4. **U12-type acceptor:** Minor spliceosome introns use AC instead of AG at the 3' splice site (~0.4% of human introns; Patel & Steitz 2003).

### Properties

- Deterministic for a given sequence and parameters.
- Score is a composite of PPT quality and PWM match.
- Higher PPT pyrimidine fraction → higher acceptor score.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n) | Linear scan over sequence |
| Space | O(1) per site | Constant working memory per candidate |

---

## Implementation Notes

**Implementation location:** [SpliceSitePredictor.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `FindAcceptorSites(sequence, minScore, includeNonCanonical)`: Scans sequence for AG (canonical) and AC (U12) dinucleotides starting at position 15. Scores each candidate via `ScoreAcceptorSite` or `ScoreU12AcceptorSite`. Returns sites with score ≥ minScore.
- `ScoreAcceptorSite(sequence, position)`: Private. Computes PPT score (count of C/U in positions [i-15, i-3), normalized by 12, scaled ×2) plus AcceptorPwm log-odds score. Normalizes to [0, 1] via `(score/(count+1) + 2) / 4`.
- `ScoreU12AcceptorSite(sequence, position)`: Private. Returns fixed 0.6 for AC dinucleotide.
- `CalculateConfidence(score, minExpected, maxExpected)`: Utility. Linear interpolation clamped to [0, 1].

---

## Deviations and Assumptions

1. **Sparse AcceptorPwm:** The implementation uses only 8 PWM positions (-15, -10, -5, -4, -3, -2, -1, 0) rather than a full position-by-position matrix for all upstream positions. This is a simplification compared to MaxEntScan's 23-position model.

2. **PPT scoring separate from PWM:** The implementation scores PPT quality separately (pyrimidine count in [-15, -3)) and adds it to the PWM score, rather than integrating PPT into the PWM matrix. This is an implementation-specific design choice.

3. **Normalization heuristic:** The formula `(score/(count+1) + 2) / 4` is not derived from any published method. It is an empirical mapping of raw log-odds to [0, 1].

4. **U12 acceptor fixed score:** The U12 acceptor scoring returns a fixed 0.6 for any AC dinucleotide without context-dependent scoring. Real U12 acceptor models use dedicated PWMs.

5. **Position reported as i+1:** `FindAcceptorSites` reports the position *after* the AG dinucleotide (the first exonic nucleotide), not the AG itself.

---

## Sources

- Shapiro MB, Senapathy P (1987). Nucleic Acids Research 15(17):7155–7174.
- Burge CB, Tuschl T, Sharp PA (1999). The RNA World, 2nd ed. CSHL Press, pp. 525–560.
- Yeo G, Burge CB (2004). J Comput Biol 11(2–3):377–394.
- Patel AA, Steitz JA (2003). Nat Rev Mol Cell Biol 4(12):960–970.
