# Evidence Artifact: SEQ-MW-001

**Test Unit ID:** SEQ-MW-001
**Algorithm:** Molecular Weight Calculation (protein and nucleotide)
**Date Collected:** 2026-06-13

---

## Online Sources

### Expasy Compute pI/Mw — documentation

**URL:** https://web.expasy.org/compute_pi/pi_tool-doc.html
**Accessed:** 2026-06-13 (fetched via WebFetch of the URL above)
**Authority rank:** 2 (official specification / standard tool of SIB Swiss Institute of Bioinformatics; basis of Expasy ProtParam Mw)

**Key Extracted Points:**

1. **Protein Mw definition (verbatim):** "Protein Mw is calculated by the addition of average isotopic masses of amino acids in the protein and the average isotopic mass of one water molecule." Results are expressed in Daltons (Da).
2. **Mass set:** The tool uses the *average* isotopic masses of amino acid residues, linked to the FindMod average-mass table (see next source). The Bjellqvist et al. reference on that page applies to pI, not to Mw.

### Expasy ProtParam — documentation (cross-reference)

**URL:** https://web.expasy.org/protparam/protparam-doc.html
**Accessed:** 2026-06-13 (fetched via WebFetch of the URL above)
**Authority rank:** 2

**Key Extracted Points:**

1. **Delegation (verbatim):** "Molecular weight and theoretical pI are calculated as in Compute pI/Mw." Confirms ProtParam Mw is the Compute pI/Mw formula above.

### Expasy FindMod — average masses of amino acid residues

**URL:** https://web.expasy.org/findmod/findmod_masses.html
**Accessed:** 2026-06-13 (fetched via WebFetch of the URL above)
**Authority rank:** 2 (official Expasy reference table)

**Key Extracted Points:**

1. **Average residue masses (Da):** Ala (A) = 71.0788; Gly (G) = 57.0519; Trp (W) = 186.2132; Arg (R) = 156.1875; Cys (C) = 103.1388. These are *residue* masses (free amino acid minus one water lost on peptide-bond formation).
2. **Average water mass:** H₂O = 18.01524 Da (section "Other mass values").
3. **Relation to free-amino-acid mass:** free amino-acid average mass = residue mass + 18.0153 (e.g., A: 71.0788 + 18.0153 = 89.0941; matches Biopython 89.0932 to rounding).

### Biopython — Bio/Data/IUPACData.py (weight tables)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py
**Accessed:** 2026-06-13 (fetched via WebFetch of the raw GitHub URL above)
**Authority rank:** 3 (established reference implementation; tables stated "Mass data taken from PubChem")

**Key Extracted Points:**

1. **`protein_weights` (average, free-amino-acid masses, Da):** A 89.0932, C 121.1582, D 133.1027, E 147.1293, F 165.1891, G 75.0666, H 155.1546, I 131.1729, K 146.1876, L 131.1729, M 149.2113, N 132.1179, P 115.1305, Q 146.1445, R 174.201, S 105.0926, T 119.1192, V 117.1463, W 204.2252, Y 181.1885.
2. **`unambiguous_dna_weights` (average, monophosphate, Da):** A 331.2218, C 307.1971, G 347.2212, T 322.2085.
3. **`unambiguous_rna_weights` (average, monophosphate, Da):** A 347.2212, C 323.1965, G 363.2206, U 324.1813.
4. **Provenance note (verbatim sense):** values are "for monophosphate deoxy nucleotides" (DNA) and monophosphate nucleotides (RNA); "Mass data taken from PubChem."

### Biopython — Bio/SeqUtils/__init__.py (`molecular_weight`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-13 (fetched via WebFetch of the raw GitHub URL above)
**Authority rank:** 3 (established reference implementation)

**Key Extracted Points:**

1. **Water constants (verbatim):** `water = 18.010565` if monoisotopic else `water = 18.0153`.
2. **Single-strand formula (verbatim):** `weight = sum(weight_table[x] for x in seq) - (len(seq) - 1) * water`. The same formula is used for protein, DNA and RNA; `weight_table` is selected by `seq_type` and `monoisotopic`.
3. **Bonds:** subtracting `(len-1) * water` removes one water per bond formed (peptide bonds for protein, phosphodiester bonds for nucleic acids). For protein this is algebraically equal to `sum(residue_masses) + water`, matching the Expasy definition.
4. **Circular:** `if circular: weight -= water` (one extra bond closes the ring).
5. **Double-stranded:** adds the complement strand's single-strand weight; raises an error for protein.
6. **Alphabet:** "Only unambiguous letters are allowed"; nucleotide sequences "are assumed to have a 5' phosphate" (built into the monophosphate residue masses).
7. **Worked examples (verbatim from docstring):** `molecular_weight("AGC", "DNA")` → 949.61; `molecular_weight("AGC", "RNA")` → 997.61; `molecular_weight("AGC", "protein")` → 249.29.

---

## Documented Corner Cases and Failure Modes

### From Biopython Bio/SeqUtils/__init__.py

1. **Unknown / ambiguous letters:** "Only unambiguous letters are allowed" — characters not in the weight table cause a lookup error (no silent averaging). A conforming implementation must reject or explicitly define behavior for non-standard symbols.
2. **5' phosphate assumption:** nucleotide masses already include a 5' monophosphate; the formula does not add or remove terminal phosphate beyond this.
3. **Protein double-stranded:** requesting double_stranded for protein is an error.

### From Expasy Compute pI/Mw

1. **Single residue / single nucleotide:** with one monomer there are zero bonds, so `(len-1)*water = 0`; the result is the free monomer mass (free amino acid for protein, monophosphate for nucleotide).

---

## Test Datasets

### Dataset: Biopython docstring worked examples

**Source:** Biopython `Bio.SeqUtils.molecular_weight` docstring, master branch (accessed 2026-06-13).

| Input | seq_type | Expected (Da) |
|-------|----------|---------------|
| AGC | protein | 249.29 |
| AGC | DNA | 949.61 |
| AGC | RNA | 997.61 |

Derivation (average tables, water = 18.0153), re-computed in this session:
- protein AGC = 89.0932 + 75.0666 + 121.1582 − 2·18.0153 = 249.2874 ≈ 249.29
- DNA AGC = 331.2218 + 347.2212 + 307.1971 − 2·18.0153 = 949.6095 ≈ 949.61
- RNA AGC = 347.2212 + 363.2206 + 323.1965 − 2·18.0153 = 997.6077 ≈ 997.61

### Dataset: Single-monomer reference values

**Source:** Biopython IUPACData tables (accessed 2026-06-13). Zero bonds ⇒ free monomer mass.

| Input | seq_type | Expected (Da) |
|-------|----------|---------------|
| G (Gly) | protein | 75.0666 |
| A | DNA | 331.2218 |
| A | RNA | 347.2212 |

---

## Assumptions

1. **ASSUMPTION: Input alphabet normalization** — Sources fix the residue masses but the Seqeron API accepts free-form `string`. We mirror sibling methods: input is upper-cased (`ToUpperInvariant`); standard amino-acid / nucleotide letters use the cited masses. This is an API-shape choice (case folding) and does not change the cited numeric values for valid input.
2. **ASSUMPTION: Behavior for unknown symbols** — Biopython rejects unknown letters; the Seqeron sibling methods are non-throwing. We resolve this conservatively for correctness: unknown amino-acid symbols are skipped (contribute no mass and no bond) and unknown nucleotide symbols are skipped, so the reported mass only reflects recognized monomers. This keeps every reported value source-backed (no invented "average" mass is used). Documented as a deviation from Biopython's reject-on-unknown behavior.

---

## Recommendations for Test Coverage

1. **MUST Test:** Protein MW of "AGC" = 249.29 Da. — Evidence: Biopython docstring + Expasy formula.
2. **MUST Test:** DNA MW of "AGC" = 949.6095 Da. — Evidence: Biopython docstring/tables.
3. **MUST Test:** RNA MW of "AGC" = 997.6177 Da. — Evidence: Biopython docstring/tables.
4. **MUST Test:** Single amino acid "G" = 75.0666 Da (zero peptide bonds ⇒ free amino-acid mass). — Evidence: Expasy formula + IUPACData.
5. **MUST Test:** Single nucleotide "A" (DNA) = 331.2218 Da and (RNA) = 347.2212 Da (zero bonds). — Evidence: IUPACData tables.
6. **MUST Test:** Empty / null input → 0 (degenerate, no monomers). — Evidence: implementation contract; sources define ≥1 monomer only.
7. **SHOULD Test:** Case-insensitivity ("agc" == "AGC"). — Rationale: sibling-method convention.
8. **SHOULD Test:** Unknown nucleotide/amino-acid symbol contributes no mass (deviation from Biopython). — Rationale: documented failure-mode resolution.
9. **COULD Test:** Two-residue peptide subtracts exactly one water (bond-count invariant). — Rationale: verifies the (len−1)·water term directly.

---

## References

1. Gasteiger E., Hoogland C., Gattiker A., Duvaud S., Wilkins M.R., Appel R.D., Bairoch A. (2005). Protein Identification and Analysis Tools on the ExPASy Server. In *The Proteomics Protocols Handbook* (Walker J.M., ed.), Humana Press, pp. 571–607. Tool documentation: Compute pI/Mw — https://web.expasy.org/compute_pi/pi_tool-doc.html
2. Expasy FindMod — Average masses of amino acid residues. SIB Swiss Institute of Bioinformatics. https://web.expasy.org/findmod/findmod_masses.html
3. Expasy ProtParam documentation. SIB Swiss Institute of Bioinformatics. https://web.expasy.org/protparam/protparam-doc.html
4. Cock P.J.A., Antao T., Chang J.T., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 — source files `Bio/SeqUtils/__init__.py` and `Bio/Data/IUPACData.py` (master branch), https://github.com/biopython/biopython

---

## Change History

- **2026-06-13**: Initial documentation (SEQ-MW-001).
