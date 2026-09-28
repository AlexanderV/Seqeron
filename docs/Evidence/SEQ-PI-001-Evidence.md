# Evidence Artifact: SEQ-PI-001

**Test Unit ID:** SEQ-PI-001
**Algorithm:** Isoelectric Point (pI) Calculation
**Date Collected:** 2026-06-13

---

## Online Sources

### EMBOSS `iep` application documentation

**URL:** https://emboss.sourceforge.net/emboss/apps/iep.html
**Accessed:** 2026-06-13 (fetched via WebFetch of the URL above)
**Authority rank:** 3 (reference implementation — EMBOSS)

**Key Extracted Points:**

1. **Purpose:** "Calculate the isoelectric point of a protein from its amino acid composition assuming that no electrostatic interactions change the propensity for ionization."
2. **pI definition / method:** The program reports, at a series of pH values, the number of bound electrons and the net charge; the isoelectric point is the pH at which net charge = 0 (where positive and negative charges balance).
3. **pKa values (Epk.dat), quoted verbatim from the fetched page:** Amino (N-terminus) = 8.6; Carboxyl (C-terminus) = 3.6; C (Cysteine) = 8.5; D (Aspartic acid) = 3.9; E (Glutamic acid) = 4.1; H (Histidine) = 6.5; K (Lysine) = 10.8; R (Arginine) = 12.5; Y (Tyrosine) = 10.1.

### Peptides R package — `charge_pI.cpp` (net charge formula)

**URL:** https://raw.githubusercontent.com/cran/Peptides/master/src/charge_pI.cpp
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation — CRAN Peptides, Osorio et al. 2015)

**Key Extracted Points:**

1. **Net charge formula (Henderson–Hasselbalch, Moore 1985):** positive (basic) groups contribute `+1 / (1 + 10^(pH − pKa))`; negative (acidic) groups contribute `−1 / (1 + 10^(pKa − pH))` (written in the source as `-1 / (1 + 10^(-(pH - pKa)))`). Verbatim from `calculateCharge()`: positive term `1.0 / (1.0 + pow(10, (1 * (pH - pKvalue(scale,'n')))))`; negative term `(-1 / (1 + pow(10, (-1 * (pH - pKvalue(scale, seq[i]))))))`.
2. **Termini:** N-terminus is treated as a positive group (key `'n'`), C-terminus as a negative group (key `'c'`); both are added once per chain.
3. **Acidic vs basic grouping:** acidic (negative) = D, E, C, Y, C-term; basic (positive) = R, K, H, N-term.

### Peptides R package — `charge()` documentation (worked example, EMBOSS scale)

**URL:** https://rdrr.io/cran/Peptides/man/charge.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation — CRAN Peptides)

**Key Extracted Points:**

1. **Available pKa scales:** "Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon", "Stryer".
2. **Worked example (EMBOSS scale), sequence `FLPVLAGLTPSIVPKLVCLLTKKC`:** net charge = 3.037398 at pH 5; 2.914112 at pH 7; 0.7184524 at pH 9.
3. **Charge basis:** "the net charge of a protein sequence based on the Henderson-Hasselbalch equation described by Moore, D. S. (1985)."

### seqinr R package — `computePI` documentation (worked example, Bjellqvist scale)

**URL:** https://rdrr.io/cran/seqinr/man/computePI.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation — CRAN seqinr; replicates ExPASy Compute pI)

**Key Extracted Points:**

1. **Definition:** "Isoelectric point is the pH at which the protein has a neutral charge."
2. **Method:** uses pK values of Bjellqvist et al., same algorithm as ExPASy Compute pI.
3. **Worked example (Bjellqvist scale):** pI of `ACDEFGHIKLMNPQRSTVWY` = 6.78454. (Used as a cross-scale reference only — this unit targets the EMBOSS scale, which gives a different pI for the same sequence; see Assumptions.)

### ExPASy Compute pI/Mw documentation

**URL:** https://web.expasy.org/compute_pi/pi_tool-doc.html (summarized from WebSearch result snippet of the same page)
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference tool documentation)

**Key Extracted Points:**

1. **Method:** pI is computed from pK values of Bjellqvist et al. (1993), determined from polypeptide migration in immobilised pH gradient gels.
2. **Limitation:** predictions for highly basic proteins and small proteins can be problematic; buffer capacity affects accuracy.

## 2026-09 review — reference-implementation evidence (supersedes the EMBOSS values above)

### EMBOSS 6.6.0 source + `iep` binary
**Opened:** `https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat`, `.../nucleus/embiep.c`; Debian `emboss` 6.6.0 `iep` binary run locally; EMBOSS iep 6.6 doc page (`https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/iep.html`).

1. The shipped `Epk.dat` has **`Amino 7.5`**, `Carboxyl 3.6` (side chains as above); `embIepPkReadFile` defaults to amino 7.50 / carboxyl 3.55 when absent. The doc page's printed Epk.dat (`Amino 8.6`, also used by the Peptides R "EMBOSS" scale) is stale: the doc's own worked outputs are reproduced only with 7.5 — LACI_ECOLI `iep` = **6.8385** (8.6 → 6.8820); IFNA2_HUMAN `-disulphides 2 -lysinemodified 2` = **5.7240** (8.6 → 5.7322).
2. `embIepGetCharge`: basic K,R,H,N-term fraction protonated `H/(H+K)`; acidic Y,C,D,E,C-term `c − pro` — identical to Henderson–Hasselbalch.
3. `embIepCompC`: B → D `(int)(0.5 + nB·5.5/9.8)` (rest N); Z → E `(int)(0.5 + nZ·6.0/9.9)` (rest Q).
4. `embIepPhConverge`: bisection over pH [1, 14] to 1e-4; no sign change → no pI.

| Sequence | `iep` pI | Sequence | `iep` pI |
|---|---|---|---|
| A, AG | 5.5500 | DDDD | 3.2279 |
| D | 3.7498 | KKKK | 11.2772 |
| E | 3.8497 | DDDDDDDD | 2.9549 |
| K | 9.1500 | RRRRRRRR | 13.3450 |
| R | 10.0000 | ACDEFGHIKLMNPQRSTVWY | 6.9681 |
| H | 7.0004 | FLPVLAGLTPSIVPKLVCLLTKKC | 9.5678 |
| C | 5.5289 | DKDK = KDDK | 5.9023 |
| Y | 5.5494 | PETER | 4.2577 |
| B / AB | 3.7498 | MAEGEITTFT | 3.6135 |
| Z / AZ | 3.8497 | BBBB / ZZZZ | 3.4918 / 3.6135 |

Charge table (`iep`, 2 dp): FLPVLAGLTPSIVPKLVCLLTKKC +3.03 / +2.70 / +0.46 at pH 5/7/9; LACI_ECOLI +7.75 / −0.63 / −5.99.

### Biopython 1.88 `Bio.SeqUtils.IsoelectricPoint` (Bjellqvist 1993/1994)
**Opened:** installed source `Bio/SeqUtils/IsoelectricPoint.py`; Peptides `charge_pI.cpp` Bjellqvist side-chain table (identical); WebSearch snippet of Biopython docs for the terminal table.

pKs: Nterm 7.5, K 10.0, R 12.0, H 5.98; Cterm 3.55, D 4.05, E 4.45, C 9.0, Y 10.0; N-terminal residue A 7.59, M 7.0, S 6.93, P 8.36, T 6.82, V 7.44, E 7.7; C-terminal residue D 4.55, E 4.75. `pi()` bisects over [4.05, 12] to 1e-4 (clamps outside).

| Sequence | Biopython `pi()` | exact root (brentq on `charge_at_pH`) |
|---|---|---|
| ACDEFGHIKLMNPQRSTVWY | 6.784552 (seqinr doc 6.78454) | 6.784543 |
| INGAR / PETER | 9.7500 / 4.5321 | 9.750025 / 4.532088 |
| FLPVLAGLTPSIVPKLVCLLTKKC | 9.3902 | 9.390178 |
| A / D / K / E | 5.5700 / 4.2994 / 8.7501 / 4.5993 | same |
| AKD / SKE / PKD | 6.1315 / 5.9377 / 6.5108 | same |
| KKKK / LACI_ECOLI | 10.4777 / 6.3901 | same |
| MAEGEITTFT / DDDD / RRRRRRRR | 4.05 / 4.05 / 12.0 (clamped) | 3.794634 / 3.521692 / 12.845100 |

`charge_at_pH`: FLPV… 3.030883 / 2.737303 / 0.757930 (pH 5/7/9); INGAR 0.760092 (pH 7); PETER −1.035860 (pH 7).

---

## Documented Corner Cases and Failure Modes

### From ExPASy Compute pI/Mw documentation

1. **Small proteins / highly basic proteins:** predicted pI may be inaccurate; poor buffer capacity increases error. (Accuracy caveat, not a correctness rule for the charge model.)

### From EMBOSS `iep`

1. **No electrostatic interactions:** the model assumes each ionizable group titrates independently (no coupling), so pI is a function of amino-acid composition only, not sequence order.

---

## Test Datasets

### Dataset: Peptides EMBOSS-scale net-charge reference

**Source:** Peptides R package `charge()` documentation, EMBOSS scale (https://rdrr.io/cran/Peptides/man/charge.html)

| Parameter | Value |
|-----------|-------|
| Sequence | FLPVLAGLTPSIVPKLVCLLTKKC |
| Net charge @ pH 5 | 3.037398 |
| Net charge @ pH 7 | 2.914112 |
| Net charge @ pH 9 | 0.7184524 |
| pKa scale | EMBOSS |

### Dataset: EMBOSS pKa table (Epk.dat)

**Source:** EMBOSS `iep` documentation (https://emboss.sourceforge.net/emboss/apps/iep.html)

| Group | pKa | Sign |
|-------|-----|------|
| N-terminus | 8.6 | + |
| C-terminus | 3.6 | − |
| C (Cys) | 8.5 | − |
| D (Asp) | 3.9 | − |
| E (Glu) | 4.1 | − |
| H (His) | 6.5 | + |
| K (Lys) | 10.8 | + |
| R (Arg) | 12.5 | + |
| Y (Tyr) | 10.1 | − |

### Dataset: Derived pI values (SUPERSEDED 2026-09 — used the stale Amino 8.6; see the 2026-09 review section)

**Source:** Derived in this session from the EMBOSS pKa table + the Peptides charge formula. The charge function used was independently confirmed against the Peptides EMBOSS worked example above (reproduces 3.037398 / 2.914112 / 0.718452 exactly to 6 dp), so these pI values are traceable to the retrieved sources.

| Sequence | Derived pI (±0.01) | Notes |
|----------|--------------------|-------|
| FLPVLAGLTPSIVPKLVCLLTKKC | 9.67 | charge crosses 0 above pH 9 (consistent with +0.72 @ pH 9) |
| A | 6.10 | termini only: midpoint of N(8.6)/C(3.6) = 6.10 |
| AG | 6.10 | termini only (no ionizable side chains) |
| D | 3.75 | one acidic side chain + termini |
| K | 9.70 | one basic side chain + termini |
| DDDD | 3.23 | acidic-dominated |
| KKKK | 11.27 | basic-dominated |
| ACDEFGHIKLMNPQRSTVWY | 7.36 | one of each residue (EMBOSS scale) |

---

## Assumptions

1. **ASSUMPTION: Empty / null sequence returns neutral 7.0.** No authoritative source defines pI for a zero-length protein (a real protein always has both termini). The repository convention (sibling statistics methods return a neutral/zero sentinel for empty input) is followed: empty or null → 7.0. This is a non-correctness-affecting input-guard convention, not an algorithm output; documented and tested as such.
2. **ASSUMPTION: pKa scale selection = EMBOSS.** Multiple published scales exist (Bjellqvist, EMBOSS, Lehninger, …) giving slightly different pI. The single-pKa-per-residue model in this repository matches the EMBOSS scale (not the position-dependent Bjellqvist model), so the EMBOSS Epk.dat values are adopted as the authoritative constant set. The seqinr Bjellqvist worked value (6.78454) is therefore NOT used as an expected value for this implementation; it is recorded only to document scale dependence.

---

## Recommendations for Test Coverage

1. **MUST Test:** Net charge of `FLPVLAGLTPSIVPKLVCLLTKKC` equals 3.037398 / 2.914112 / 0.7184524 at pH 5/7/9 (validates the charge formula + EMBOSS pKa). — Evidence: Peptides `charge()` doc.
2. **MUST Test:** pI of `FLPVLAGLTPSIVPKLVCLLTKKC` ≈ 9.67; basic peptide. — Evidence: derived from Peptides charge example (positive charge persists past pH 9).
3. **MUST Test:** pI bounds 0 ≤ pI ≤ 14 for any input (invariant). — Evidence: bisection interval [0,14], EMBOSS iep.
4. **MUST Test:** Termini-only sequence (`A`, `AG`) → pI = midpoint of N/C-term pKa = 6.10. — Evidence: EMBOSS pKa (8.6, 3.6).
5. **SHOULD Test:** Acidic-only (`DDDD`) low pI ≈ 3.23; basic-only (`KKKK`) high pI ≈ 11.27. — Rationale: monotonic response to charge composition.
6. **SHOULD Test:** Empty and null → 7.0. — Rationale: documented input-guard convention.
7. **COULD Test:** Order-independence (pI is composition-only): permutation of a sequence yields identical pI. — Rationale: EMBOSS "no electrostatic interactions" assumption.

---

## References

1. EMBOSS. iep — Calculate the isoelectric point of proteins. EMBOSS application documentation. https://emboss.sourceforge.net/emboss/apps/iep.html (accessed 2026-06-13).
2. Osorio D, Rondón-Villarreal P, Torres R. (2015). Peptides: A Package for Data Mining of Antimicrobial Peptides. The R Journal 7(1):4–14. Source: `src/charge_pI.cpp`, `charge()` doc. https://github.com/cran/Peptides ; https://rdrr.io/cran/Peptides/man/charge.html (accessed 2026-06-13).
3. Charif D, Lobry JR. (seqinr). computePI — theoretical isoelectric point. CRAN seqinr documentation. https://rdrr.io/cran/seqinr/man/computePI.html (accessed 2026-06-13).
4. Bjellqvist B, Hughes GJ, Pasquali C, Paquet N, Ravier F, Sanchez JC, Frutiger S, Hochstrasser D. (1993). The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. Electrophoresis 14:1023–1031. https://doi.org/10.1002/elps.11501401163 (DOI confirmed via WebSearch; PubMed PMID 8125050).
5. ExPASy. Compute pI/Mw documentation. https://web.expasy.org/compute_pi/pi_tool-doc.html (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation (SEQ-PI-001).
- **2026-09-28**: Review 2026-09 — EMBOSS source/binary evidence (Amino 7.5; the 8.6 listing and the "Derived pI" table above are superseded), Biopython Bjellqvist evidence.
