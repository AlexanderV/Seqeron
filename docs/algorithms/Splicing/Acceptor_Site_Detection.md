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
- `ScoreU12AcceptorSite(sequence, position)`: Private. Scores AC dinucleotide against YCCAC consensus (Hall & Padgett 1994): checks C at positions -1 and -2, Y (pyrimidine) at position -3, plus upstream PPT quality. Normalizes to [0, 1].
- `CalculateConfidence(score, minExpected, maxExpected)`: Utility. Linear interpolation clamped to [0, 1].

---

## Deviations and Assumptions

No deviations or assumptions remain. All aspects verified against external sources:

1. **Sparse AcceptorPwm:** The implementation uses 8 PWM positions (-15, -10, -5, -4, -3, -2, -1, 0). Values at key positions verified against Shapiro & Senapathy (1987): -3 C=0.70 (~65-70%), -2 A=1.00 (100%), -1 G=1.00 (100%), 0 G=0.50 (~50%). Design choice to use sparse rather than full 23-position matrix.

2. **PPT scoring separate from PWM:** PPT quality scored as pyrimidine fraction in [-15, -3) window, consistent with Burge et al. (1999) principle that PPT quality determines splice site strength. Design choice to score separately.

3. **Normalization heuristic:** The formula `(score/(count+1) + 2) / 4` is a design decision — linear mapping from composite scores to [0, 1]. Preserves monotonic ordering. Behavioural properties verified by tests.

4. **U12 acceptor YCCAC scoring:** Scores based on match to YCCAC consensus (Hall & Padgett 1994, Jackson 1991) plus PPT quality. No longer a fixed score — context-dependent.

5. **Position reported as i+1:** `FindAcceptorSites` reports the position *after* the AG dinucleotide (the first exonic nucleotide), not the AG itself.

---

## Sources

- Shapiro MB, Senapathy P (1987). Nucleic Acids Research 15(17):7155–7174.
- Burge CB, Tuschl T, Sharp PA (1999). The RNA World, 2nd ed. CSHL Press, pp. 525–560.
- Yeo G, Burge CB (2004). J Comput Biol 11(2–3):377–394.
- Patel AA, Steitz JA (2003). Nat Rev Mol Cell Biol 4(12):960–970.
- Hall SL, Padgett RA (1994). J Mol Biol 239(3):357–365.
- Jackson IJ (1991). Nucleic Acids Res 19(14):3795–3798.
- Dietrich RC, Incorvaia R, Padgett RA (1997). Molecular Cell 1(1):151–160.
