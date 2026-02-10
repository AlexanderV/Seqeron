# Pre-miRNA Hairpin Detection

## Documented Theory

### Purpose

Pre-miRNA detection identifies potential precursor microRNA (pre-miRNA) hairpin structures in RNA sequences. Pre-miRNAs are ~60–120 nt stem-loop structures that serve as intermediates in miRNA biogenesis, produced by Drosha/DGCR8 cleavage of pri-miRNA and subsequently processed by Dicer into mature ~22 nt miRNA duplexes (Bartel 2004, Cell 116:281-297).

### Core Mechanism

Computational pre-miRNA detection relies on structural features:

1. **Stem identification:** Scan for self-complementary regions that can form a double-stranded stem. Watson-Crick (A-U, G-C) and wobble (G-U) base pairs are counted.
2. **Loop identification:** The unpaired region between the two stem arms forms a terminal loop.
3. **Structural validation:** Check that stem length (≥18 bp), loop size (3–25 nt), and overall length (~55–120 nt) match known pre-miRNA parameters.
4. **Mature/star extraction:** The mature miRNA (~22 nt) is typically extracted from the 5' arm, with the star sequence from the 3' arm.
5. **Energy estimation:** Approximate folding free energy indicates thermodynamic stability.

Formal criteria (Ambros et al. 2003, RNA 9:277-279):
- Must fold into a hairpin precursor structure
- ~22 nt mature miRNA resides in one arm of the hairpin

### Properties

- **Deterministic:** Same input always produces same output.
- **Scanning approach:** Examines all windows of valid length.
- **Simplified model:** Uses consecutive base pairing from ends (no bulge/mismatch tolerance in stem).

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n² × L) where n = sequence length, L = max hairpin length | Implementation analysis |
| Space | O(L) per candidate | Implementation analysis |

---

## Implementation Notes

**Implementation location:** [MiRnaAnalyzer.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength)`: Public method. Scans input sequence for candidate subsequences, converts T→U, calls `AnalyzeHairpin` on each candidate. Returns `IEnumerable<PreMiRna>`.
- `AnalyzeHairpin(sequence, matureLength)`: Private helper. Validates hairpin by checking consecutive complementary bases from ends inward. Requires stem ≥ 18 bp, loop 3–25 nt. Extracts mature/star sequences, generates dot-bracket notation, estimates free energy.

### PreMiRna Record Fields

| Field | Description |
|-------|-------------|
| Start | 0-based start position in input sequence |
| End | 0-based end position (inclusive) |
| Sequence | The candidate hairpin subsequence (RNA, uppercase) |
| MatureSequence | Extracted mature miRNA from 5' arm |
| StarSequence | Extracted star sequence from 3' arm |
| Structure | Dot-bracket notation string |
| FreeEnergy | Estimated folding energy (kcal/mol approximation) |

---

## Deviations and Assumptions

1. **Consecutive stem requirement:** The implementation requires consecutive base pairs from position 0 inward. Real pre-miRNAs contain internal bulges and mismatches. This is a simplification.
2. **Simplified energy model:** Uses `-stemLength * 1.5 + loopSize * 0.5` instead of Turner nearest-neighbor parameters. Values are approximations, not thermodynamically precise.
3. **No bulge tolerance:** A single non-pairing base breaks the stem count. Genuine pre-miRNAs often have 1–3 nt internal loops/bulges.
4. **5' arm default:** Mature miRNA is always extracted from 5' arm. In biology, either arm can produce the dominant strand.

---

## Sources

- Bartel DP (2004). MicroRNAs: Genomics, Biogenesis, Mechanism, and Function. Cell 116:281-297.
- Ambros V et al. (2003). A Uniform System for microRNA Annotation. RNA 9:277-279.
- Bartel DP (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell 136:215-233.
- Krol J et al. (2004). Structural Features of miRNA Precursors. J Biol Chem 279:42230-42239.
- Wikipedia: MicroRNA. https://en.wikipedia.org/wiki/MicroRNA
- miRBase. https://mirbase.org/
