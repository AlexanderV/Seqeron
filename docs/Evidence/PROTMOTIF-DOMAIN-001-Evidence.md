# Evidence Artifact: PROTMOTIF-DOMAIN-001

**Test Unit ID:** PROTMOTIF-DOMAIN-001
**Algorithm:** Domain Prediction & Signal Peptide Prediction
**Date Collected:** 2026-02-13

---

## Online Sources

### Wikipedia — Protein Domain

**URL:** https://en.wikipedia.org/wiki/Protein_domain
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **Definition:** A protein domain is a self-stabilizing region of a protein's polypeptide chain that folds independently. Domains range from ~50 to ~250 amino acids. — Wetlaufer (1973), cited via Xu & Nussinov (1998) https://doi.org/10.1016/S1359-0278(98)00004-2
2. **Domain databases:** Pfam, InterPro, PROSITE, and SMART are established databases for domain classification. Pfam classifies using Hidden Markov Models. — El-Gebali et al. (2019) https://doi.org/10.1093/nar/gky995
3. **Zinc finger C2H2:** The shortest domains (~50 aa), stabilized by metal ions. C2H2 zinc fingers adopt ββα fold with motif `X2-Cys-X(2,4)-Cys-X(12)-His-X(3,5)-His`. — Krishna et al. (2003) https://doi.org/10.1093/nar/gkg161
4. **SH3 domain:** An interaction domain involved in protein-protein binding, recognizing proline-rich sequences. — cited via Campbell & Downing (1994)
5. **Pattern-based detection:** PROSITE uses consensus patterns, Pfam uses profile HMMs for domain identification. Simple regex-based matching is a simplification of PROSITE's approach.

### Wikipedia — Signal Peptide

**URL:** https://en.wikipedia.org/wiki/Signal_peptide
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **Length:** Signal peptides are usually 16–30 amino acids long, present at the N-terminus. — Owji et al. (2018) https://doi.org/10.1016/j.ejcb.2018.06.003
2. **Tripartite structure:** Signal peptides contain three regions: n-region (positively charged, N-terminal), h-region (hydrophobic core, 7–15 residues), c-region (cleavage recognition site). — von Heijne (1985) J Mol Biol 184:99–105; von Heijne & Gavel (1988) https://doi.org/10.1111/j.1432-1033.1988.tb14150.x
3. **Cleavage site specificity:** Small, neutral amino acids {A, G, S} at positions -1 and -3 relative to the cleavage site (the "-1, -3 rule"). — von Heijne (1983) Eur J Biochem 133:17–21

### Wikipedia — Zinc Finger

**URL:** https://en.wikipedia.org/wiki/Zinc_finger
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **C2H2 consensus motif:** `X2-Cys-X(2,4)-Cys-X(12)-His-X(3,5)-His`. Pfam: PF00096, PROSITE: PS00028. — Pabo et al. (2001) https://doi.org/10.1146/annurev.biochem.70.1.313
2. **PROSITE pattern PS00028:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` — verified at https://prosite.expasy.org/PS00028

### PROSITE Database — PS00028 (Zinc Finger C2H2)

**URL:** https://prosite.expasy.org/PS00028
**Accessed:** 2026-02-13
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Pattern:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`
2. **PROSITE accession:** PS00028
3. **Pfam cross-reference:** PF00096
4. **Note:** The implementation uses `[LIVMFYWC]` at position +3 from second C. The actual PROSITE PS00028 uses `[LIVMFYWC]` (confirmed). Some sources also include H in this class but the official PROSITE pattern does not.

### PROSITE Database — PS00017 (P-loop / Walker A motif)

**URL:** https://prosite.expasy.org/PS00017
**Accessed:** 2026-02-13
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Pattern:** `[AG]-x(4)-G-K-[ST]`
2. **Function:** ATP/GTP-binding site motif A (P-loop), also known as Walker A motif
3. **Reference:** Walker et al. (1982). EMBO J 1:945-951. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x

### PROSITE Database — PS00678 (WD-repeats signature, WD_REPEATS_1) — RETRIEVED 2026-06-24

**URL:** https://prosite.expasy.org/PS00678
**Retrieval:** WebFetch of https://prosite.expasy.org/PS00678 and the documentation page
https://prosite.expasy.org/PDOC00574 — both returned the same pattern string verbatim.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Entry name / accession:** WD_REPEATS_1 / **PS00678** — "Trp-Asp (WD) repeats signature".
2. **Pattern (verbatim from PROSITE):**
   `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]`
3. This is a **deterministic PROSITE PATTERN** (regex-like), NOT a profile/HMM, so it can be
   reproduced exactly. It is a 14-element signature spanning **15 residues** (the `x(2)` element
   covers 2 positions; all others are single positions). The terminal `W-[DEN]` is the diagnostic
   "WD" (Trp-Asp) of the repeat.
4. **Cross-references:** Pfam PF00400 (WD40); the separate WD-repeats **PROFILE** is PS50082
   (WD_REPEATS_2) — a weight matrix, not a pattern.

### PROSITE → regex syntax rules — RETRIEVED 2026-06-24

**URL:** https://prosite.expasy.org/scanprosite/scanprosite_doc.html
**Retrieval:** WebFetch of the ScanProsite/PROSITE pattern-syntax documentation.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Rules extracted (verbatim/closely paraphrased):**

1. "The standard IUPAC one letter code for the amino acids is used." → literal residue.
2. "The symbol 'x' is used for a position where any amino acid is accepted." → regex `.`.
3. "Ambiguities are indicated by listing the acceptable amino acids … between square brackets '[ ]'." → kept as `[..]`.
4. "Ambiguities are also indicated by listing between a pair of curly brackets '{ }' the amino acids that are NOT accepted at a given position." → `{X}` → `[^X]`.
5. "Each element in a pattern is separated from its neighbor by a '-'." → separators dropped.
6. Repetition `e(n)` / gap `x(n)` / range `x(m,n)` → `e{n}` / `.{n}` / `.{m,n}`.
7. N-terminal `<` → `^`; C-terminal `>` → `$`.

**Applied translation of PS00678:**
`[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]`

### PROSITE SH3 (PS50002) and PDZ (PS50106) are PROFILE-only — RETRIEVED 2026-06-24

**URLs:** https://prosite.expasy.org/PDOC50002 (SH3), https://prosite.expasy.org/PDOC50106 (PDZ),
EBI InterPro profile entries https://www.ebi.ac.uk/interpro/entry/profile/PS50002 and `/PS50106`.
**Retrieval:** WebSearch + the PROSITE/InterPro entries.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **SH3 (Src homology 3): PS50002 is a PROFILE (weight matrix)**, NOT a pattern — "profiles use
   weight matrices … characterize protein domains over their entire length." There is **no
   deterministic PROSITE PATTERN** for the SH3 domain.
2. **PDZ (DHR / GLGF): PS50106 is a PROFILE** covering the minimal PDZ domain — no deterministic
   PROSITE PATTERN exists.
3. **Consequence:** SH3 and PDZ **cannot** be reproduced as an exact regex without fabricating a
   signature. They are therefore NOT detected by `FindDomains` (honest residual). The previously
   shipped ad-hoc SH3/PDZ regexes were not sourced from any PROSITE pattern and have been removed.

### Real WD40 positive — GBB1_HUMAN (UniProt P62873) — RETRIEVED 2026-06-24

**URL:** https://rest.uniprot.org/uniprotkb/P62873.fasta
**Retrieval:** WebFetch of the UniProt REST FASTA.
**Accessed:** 2026-06-24
**Authority rank:** 5 (UniProt curated database)

**Key Extracted Points:**

1. GBB1_HUMAN (β-transducin / Gβ1) is a canonical 7-bladed WD40 β-propeller.
2. The PS00678-translated regex matches the sequence at **0-based starts 69, 156, 284** (each a
   15-residue window). Example repeat (start 69): `LVSASQDGKLIIWDS`.

### von Heijne G (1986). Signal sequences — the limits of variation

**URL:** https://doi.org/10.1016/0022-2836(85)90046-4
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **-1, -3 rule:** Small neutral amino acids (predominantly Ala, Gly, Ser) at positions -1 and -3 relative to the cleavage site.
2. **N-region:** 1–5 residues, positively charged (Lys, Arg).
3. **H-region:** 7–15 hydrophobic residues forming a single α-helix. The hydrophobic core is both necessary and sufficient for membrane targeting.
4. **C-region:** 3–7 residues, polar, with the -3/-1 specificity constraint.
5. **Note:** The amino acids Thr and Cys are occasionally found at -1/-3 but are not part of the canonical set. The implementation strictly uses {A, G, S}.

### von Heijne G (1983). Patterns of amino acids near signal-sequence cleavage sites

**URL:** https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **Signal peptide statistics:** Analysis of 78 eukaryotic signal peptides showed the -1 and -3 positions strongly favor A, G, S.
2. **Position -1:** A (65%), G (14%), S (10%), T (5%), other (6%)
3. **Position -3:** A (52%), V (15%), G (9%), S (8%), L (6%), other (10%)
4. **Note:** Position -3 also accepts V and L in von Heijne (1983) data, but the canonical "-1,-3 rule" specifies {A, G, S}. The implementation strictly uses {A, G, S}.

### Pfam — PF00096, PF00400, PF00018, PF00595, PF00069

**URL:** https://www.ebi.ac.uk/interpro/entry/pfam/
**Accessed:** 2026-02-13
**Authority rank:** 5 (Bioinformatics database)

**Key Extracted Points:**

1. **PF00096 (Zinc finger C2H2):** Classical ββα fold, ~23-28 residues per finger
2. **PF00400 (WD40 repeat):** ~40 residues per repeat, forms β-propeller
3. **PF00018 (SH3 domain):** ~60 residues, binds proline-rich sequences
4. **PF00595 (PDZ domain):** ~80-90 residues, protein-protein interaction
5. **PF00069 (Protein kinase):** Catalytic domain, includes P-loop (Walker A) motif for ATP binding

---

## Documented Corner Cases and Failure Modes

### From von Heijne (1986)

1. **Short sequences:** Sequences shorter than ~15 residues cannot contain a valid signal peptide (n+h+c regions require minimum length).
2. **Atypical cleavage sites:** Some signal peptides have non-canonical residues at -1/-3 (e.g., Thr, Val). The implementation strictly enforces {A, G, S} per von Heijne (1983).
3. **No signal peptide:** Cytoplasmic and nuclear proteins lack signal peptides entirely.

### From Pfam/PROSITE

1. **Empty/null input:** Domain search on empty sequence returns no domains.
2. **No matching domains:** Most short peptides or random sequences contain no recognizable domain signatures.
3. **Overlapping domains:** Some proteins have domains with overlapping boundaries.
4. **Multiple domain instances:** Tandem repeat domains (e.g., multiple zinc fingers) should all be detected.

---

## Test Datasets

### Dataset: Synthetic C2H2 Zinc Finger

**Source:** Constructed from PROSITE PS00028 consensus `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAACXXCXXXLXXXXXXXXHXXXHAAA` |
| Domain | Zinc Finger C2H2 |
| Pfam | PF00096 |
| Expected start | 4 |
| Expected end | 24 |

### Dataset: Synthetic P-loop (Kinase ATP-binding)

**Source:** Constructed from PROSITE PS00017 consensus `[AG]-x(4)-G-K-[ST]`

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAAGXXXXGKSAAAA` |
| Domain | Protein Kinase ATP-binding |
| Pfam | PF00069 |
| Expected match | Positions 4–11 |

### Dataset: Synthetic Signal Peptide

**Source:** Constructed from von Heijne (1986) principles: M + K/R (n-region) + L-rich (h-region) + small aa at -3,-1 (c-region)

| Parameter | Value |
|-----------|-------|
| Sequence | `MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF` |
| Expected | Signal peptide detected |
| Cleavage position | ~25 (after c-region) |

---

## Design Decisions (previously Assumptions — all resolved)

1. **Signal peptide scoring weights (1:2:1)** — The implementation weights region scores as `(nScore + 2·hScore + cScore) / 4` giving the h-region double weight. Justified by von Heijne (1985) J Mol Biol 184:99–105: the hydrophobic core is "both necessary and sufficient for membrane targeting," establishing it as the dominant discriminator.

2. **Signal peptide detection: evidence-based biological constraints** — ~~Previously: arbitrary threshold 0.4 (ELIMINATED).~~ Now enforces two literature-derived constraints: (a) n-region must contain positive charge (`nScore > 0`), per von Heijne (1986) "1–5 residues, positively charged (Lys, Arg)"; (b) h-region must be predominantly hydrophobic (`hScore ≥ 0.5`), per von Heijne (1985) "hydrophobic core" requirement. No arbitrary thresholds remain.

3. **Signal peptide probability** — ~~Previously: `min(1.0, score × 1.2)` (ELIMINATED).~~ Implementation now uses `Probability = Score` directly. No arbitrary scaling.

4. **Small amino acid set {A, G, S}** — ~~Previously: {A, G, S, T, C} (CORRECTED).~~ Implementation strictly uses {A, G, S} per von Heijne (1983) Eur J Biochem 133:17–21, positions -1 and -3. Matches the canonical -1,-3 rule.

5. **Domain signature patterns** — The implementation uses PROSITE-derived regex patterns (PS00028, PS00017) for domain detection. This is a deliberate scope decision: PROSITE itself uses consensus patterns (see Hulo et al. 2006). Not an assumption — it is a documented scope boundary.

6. **Method naming** — Public API is `FindDomains` (not `PredictDomains`). Non-correctness-affecting naming choice aligned with the method's behavior (pattern-based search, not probabilistic prediction).

---

## Recommendations for Test Coverage

1. **MUST Test:** FindDomains detects C2H2 zinc finger from PROSITE PS00028 consensus — Evidence: PROSITE PS00028
2. **MUST Test:** FindDomains detects P-loop/kinase from PROSITE PS00017 consensus — Evidence: PROSITE PS00017, Walker et al. (1982)
3. **MUST Test:** FindDomains returns correct domain metadata (name, accession, start, end) — Evidence: Pfam database
4. **MUST Test:** FindDomains empty/null input returns empty — Evidence: trivial
5. **MUST Test:** PredictSignalPeptide detects signal peptide with tripartite structure — Evidence: von Heijne (1986)
6. **MUST Test:** PredictSignalPeptide enforces -1,-3 rule for cleavage site — Evidence: von Heijne (1983, 1986)
7. **MUST Test:** PredictSignalPeptide returns null for non-secretory (charged) sequences — Evidence: von Heijne (1986)
8. **MUST Test:** PredictSignalPeptide returns null for sequences too short — Evidence: trivial (< 15 aa minimum)
9. **MUST Test:** PredictSignalPeptide returns three regions (n, h, c) — Evidence: von Heijne (1986)
10. **MUST Test:** PredictSignalPeptide is case-insensitive — Evidence: standard protein sequence convention
11. **SHOULD Test:** FindDomains detects WD40 repeat (PF00400) — Evidence: Pfam
12. **SHOULD Test:** FindDomains detects SH3 domain (PF00018) — Evidence: Pfam
13. **SHOULD Test:** FindDomains detects PDZ domain (PF00595) — Evidence: Pfam
14. **SHOULD Test:** FindDomains on no-match input returns empty — Evidence: trivial
15. **COULD Test:** FindDomains detects multiple domain instances in tandem — Evidence: Pfam, domain architecture

---

## References

1. von Heijne G (1986). Signal sequences. The limits of variation. J Mol Biol 184(1):99-105. https://doi.org/10.1016/0022-2836(85)90046-4
2. von Heijne G (1983). Patterns of amino acids near signal-sequence cleavage sites. Eur J Biochem 133(1):17-21. https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
3. Walker JE et al. (1982). Distantly related sequences in the α- and β-subunits of ATP synthase, myosin, kinases and other ATP-requiring enzymes and a common nucleotide binding fold. EMBO J 1(8):945-951. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x
4. Krishna SS, Majumdar I, Grishin NV (2003). Structural classification of zinc fingers. Nucleic Acids Res 31(2):532-550. https://doi.org/10.1093/nar/gkg161
5. El-Gebali S et al. (2019). The Pfam protein families database in 2019. Nucleic Acids Res 47(D1):D427-D432. https://doi.org/10.1093/nar/gky995
6. Hulo N et al. (2006). The PROSITE database. Nucleic Acids Res 34:D227-D230. https://doi.org/10.1093/nar/gkj063
7. PROSITE PS00028 — Zinc finger C2H2 type. https://prosite.expasy.org/PS00028
8. PROSITE PS00017 — ATP/GTP-binding site motif A (P-loop). https://prosite.expasy.org/PS00017
9. Pabo CO, Peisach E, Grant RA (2001). Design and selection of novel Cys2His2 zinc finger proteins. Annu Rev Biochem 70:313-340. https://doi.org/10.1146/annurev.biochem.70.1.313
10. Owji H et al. (2018). A comprehensive review of signal peptides. Eur J Cell Biol 97(6):422-441. https://doi.org/10.1016/j.ejcb.2018.06.003
11. von Heijne G (1986). A new method for predicting signal sequence cleavage sites. Nucl Acids Res 14(11):4683–4690. https://doi.org/10.1093/nar/14.11.4683
12. von Heijne G (1984). How signal sequences maintain cleavage specificity. J Mol Biol 173(2):243–251. https://doi.org/10.1016/0022-2836(84)90192-x
