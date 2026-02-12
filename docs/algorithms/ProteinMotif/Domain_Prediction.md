# Domain Prediction & Signal Peptide Prediction

## Documented Theory

### Purpose

**Domain prediction** identifies known structural/functional protein domains within a sequence by scanning for conserved signature patterns. Protein domains are compact, independently folding units typically 50–250 amino acids long (Wetlaufer, 1973; Xu & Nussinov, 1998). Detection of domains provides functional annotation of proteins.

**Signal peptide prediction** identifies N-terminal signal sequences that target proteins for secretion or membrane insertion. Signal peptides are 16–30 amino acids long and are cleaved after translocation (von Heijne, 1986; Owji et al., 2018).

### Core Mechanism

#### Domain Finding

The implementation uses regex pattern matching against known domain signature consensus sequences derived from PROSITE (Hulo et al., 2006). For each candidate domain type, a regex pattern encoding the conserved residues and spacing is matched against the input sequence.

Detected domains:
- **Zinc Finger C2H2** (PF00096/PS00028): `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`
- **WD40 Repeat** (PF00400): `[LIVMFYWC]-x(5,12)-[WF]-D`
- **SH3 Domain** (PF00018): `[LIVMF]-x(2)-[GA]-W-[FYW]-x(5,8)-[LIVMF]`
- **PDZ Domain** (PF00595): `[LIVMF]-[ST]-[LIVMF]-x(2)-G-[LIVMF]-x(3,4)-[LIVMF]-x(2)-[DEN]`
- **Protein Kinase ATP-binding** (PF00069/PS00017): `[AG]-x(4)-G-K-[ST]` (Walker A motif)

This is a simplified approach compared to full Pfam HMM profiles.

#### Signal Peptide Prediction

The implementation uses a rule-based approach based on von Heijne's tripartite model:

1. **Scanning:** Candidate cleavage sites are tested at positions 15–35.
2. **(-1, -3) rule:** Positions -1 and -3 from the cleavage site must be small amino acids {A, G, S} — von Heijne (1983) Eur J Biochem 133:17–21.
3. **Region scoring:** For each candidate site, three regions are scored:
   - **N-region** (1–5 residues): Positive charge density (K, R). Normalized so 2 charges = 1.0, based on mean n-region charge ≈ +2.0 — von Heijne (1986) Nucl Acids Res 14:4683–4690.
   - **H-region** (minimum 7 residues): Hydrophobic amino acid fraction {A, I, L, M, F, V, W}. Direct ratio measurement — von Heijne (1985) J Mol Biol 184:99–105.
   - **C-region** (5 residues): Small/polar amino acid fraction {A, G, S, T, N}. Direct ratio measurement — von Heijne (1984) J Mol Biol 173:243–251.
4. **Total score:** Weighted mean with 1:2:1 ratio (n:h:c). The h-region receives double weight reflecting its role as "both necessary and sufficient for membrane targeting" — von Heijne (1985).
5. **Threshold:** Score ≥ 0.4 required for prediction.
6. **Probability:** Equals Score (direct quality measure).
7. **Output:** Cleavage position, three regions, score, and probability.

### Properties

- **Deterministic:** Both methods produce identical output for identical input.
- **No optimality guarantee:** Pattern-based domain detection is heuristic; it may miss true domains (false negatives) or match spurious patterns (false positives). Domain patterns are PROSITE consensus patterns expressed as regex.
- **Signal peptide detection:** Rule-based predictor using the tripartite model with evidence-based parameters from von Heijne (1983, 1984, 1985, 1986).

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time (FindDomains) | O(n × d) where n = sequence length, d = number of domain patterns | Implementation analysis |
| Time (PredictSignalPeptide) | O(1) (scans fixed-length N-terminal region ≤70 aa) | Implementation analysis |
| Space | O(1) per call (output proportional to matches) | Implementation analysis |

---

## Implementation Notes

**Implementation location:** [ProteinMotifFinder.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `FindDomains(proteinSequence)`: Scans for 5 known domain signature patterns using regex matching. Returns `IEnumerable<ProteinDomain>` with name, Pfam accession, start/end positions, score, and description.
- `PredictSignalPeptide(proteinSequence, maxLength)`: Predicts signal peptide by scanning N-terminal region for cleavage sites satisfying the -1/-3 rule and scoring the tripartite structure. Returns `SignalPeptide?` (nullable).

Supporting types:
- `ProteinDomain` record: Name, Accession, Start, End, Score, Description.
- `SignalPeptide` record: CleavagePosition, NRegion, HRegion, CRegion, Score, Probability.

---

## Sources

- von Heijne G (1986). Signal sequences. The limits of variation. J Mol Biol 184(1):99-105. https://doi.org/10.1016/0022-2836(85)90046-4
- von Heijne G (1986). A new method for predicting signal sequence cleavage sites. Nucl Acids Res 14(11):4683–4690. https://doi.org/10.1093/nar/14.11.4683
- von Heijne G (1984). How signal sequences maintain cleavage specificity. J Mol Biol 173(2):243–251. https://doi.org/10.1016/0022-2836(84)90192-x
- von Heijne G (1983). Patterns of amino acids near signal-sequence cleavage sites. Eur J Biochem 133(1):17-21. https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
- Walker JE et al. (1982). Distantly related sequences in α- and β-subunits of ATP synthase. EMBO J 1(8):945-951. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x
- Krishna SS, Majumdar I, Grishin NV (2003). Structural classification of zinc fingers. Nucleic Acids Res 31(2):532-550. https://doi.org/10.1093/nar/gkg161
- Hulo N et al. (2006). The PROSITE database. Nucleic Acids Res 34:D227-D230. https://doi.org/10.1093/nar/gkj063
- El-Gebali S et al. (2019). The Pfam protein families database in 2019. Nucleic Acids Res 47(D1):D427-D432. https://doi.org/10.1093/nar/gky995
- PROSITE PS00028 — Zinc finger C2H2 type. https://prosite.expasy.org/PS00028
- PROSITE PS00017 — ATP/GTP-binding site motif A (P-loop). https://prosite.expasy.org/PS00017
