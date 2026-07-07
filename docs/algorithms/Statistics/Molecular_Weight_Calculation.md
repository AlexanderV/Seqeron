# Molecular Weight Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-MW-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Computes the average-isotopic molecular weight (in daltons) of a protein sequence or a
DNA/RNA sequence. For proteins the result is the sum of the average isotopic masses of the
constituent amino acids plus the average isotopic mass of one water molecule [1]. For nucleic
acids it is the sum of average monophosphate mononucleotide masses minus one water per
phosphodiester bond [4]. The calculation is exact (specification-driven) for the standard
alphabets; it is the same quantity reported by Expasy Compute pI/Mw / ProtParam and by
Biopython `Bio.SeqUtils.molecular_weight`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A biopolymer is a chain of monomers joined by condensation bonds, each of which releases one
water molecule. A protein of *n* amino acids has *n − 1* peptide bonds; a single-stranded
nucleic acid of *n* nucleotides has *n − 1* phosphodiester bonds. The molecular weight is the
sum of the monomer masses corrected for the water lost during polymerization [1][4].

### 2.2 Core Model

Protein (average isotopic mass) [1]:

> "Protein Mw is calculated by the addition of average isotopic masses of amino acids in the
> protein and the average isotopic mass of one water molecule." [1]

Using free-amino-acid average masses `m_aa` and water mass `W` this is equivalent to the
Biopython single-strand formula [4]:

`MW = Σ m_aa(residue_i) − (n − 1) · W`

Nucleic acid (single strand) [4]:

`MW = Σ m_nt(monomer_i) − (n − 1) · W`

where `m_nt` are average 5'-monophosphate mononucleotide masses [5] and `W = 18.0153 Da`
(average isotopic mass of water) [2][4].

Numeric tables (average masses, Da):

- Amino acids (free): A 89.0932, C 121.1582, D 133.1027, E 147.1293, F 165.1891, G 75.0666,
  H 155.1546, I 131.1729, K 146.1876, L 131.1729, M 149.2113, N 132.1179, P 115.1305,
  Q 146.1445, R 174.201, S 105.0926, T 119.1192, V 117.1463, W 204.2252, Y 181.1885 [5].
- DNA monophosphates: A 331.2218, C 307.1971, G 347.2212, T 322.2085 [5].
- RNA monophosphates: A 347.2212, C 323.1965, G 363.2206, U 324.1813 [5].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | MW(empty) = MW(null) = 0 | No monomers; the formula is defined for n ≥ 1 [1][4] |
| INV-02 | MW > 0 for any non-empty recognized sequence | Every monomer mass exceeds W [2][5] |
| INV-03 | Exactly one water removed per bond: MW(2 monomers) = m₁ + m₂ − W | (n − 1) = 1 for n = 2 [4] |
| INV-04 | Case-insensitive: MW(s) = MW(uppercase(s)) | Input is upper-cased before lookup |

### 2.5 Comparison with Related Methods

| Aspect | Average isotopic Mw (this method) | Monoisotopic Mw |
|--------|-----------------------------------|------------------|
| Mass set | Averaged over all stable isotopes | Most abundant isotope only |
| Water constant | 18.0153 Da [4] | 18.010565 Da [4] |
| Typical use | Bulk biochemistry, gel/size estimates | High-resolution mass spectrometry |

Only the average-isotopic variant is implemented; monoisotopic is not in scope (§5.3).

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | Protein, one-letter codes | Case-insensitive; standard 20 AAs recognized |
| sequence | string | required | DNA or RNA sequence | Case-insensitive; A/C/G/T (DNA) or A/C/G/U (RNA) |
| isDna | bool | true | Selects DNA vs RNA mass table | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Molecular weight in daltons (Da); 0 for null/empty |

### 3.3 Preconditions and Validation

Null or empty input returns 0 (no exception). Input is upper-cased
(`ToUpperInvariant`). Only recognized monomers contribute mass and bonds; unknown symbols are
skipped (no mass, no bond) — this deviates from Biopython, which rejects unknown letters
(§5.4). A sequence containing no recognized monomers returns 0. No coordinate system applies.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return 0 if input is null/empty.
2. Upper-case the sequence; for each character, if it is in the relevant mass table, add its
   mass and increment the monomer count.
3. If no recognized monomers, return 0.
4. Return `accumulatedMass − (monomerCount − 1) · 18.0153`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The amino-acid, DNA, and RNA average-mass tables and the water constant are the lookup
constants that define output; their values and origins are given in §2.2 and cited to
Biopython IUPACData [5] and Expasy [1][2]. The DNA/RNA table is selected by the `isDna` flag.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateMolecularWeight / CalculateNucleotideMolecularWeight | O(n) | O(1) | One pass; fixed-size lookup tables |

This is not a search/matching unit, so the repository suffix tree does not apply (no occurrence
enumeration or pattern lookup is performed).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateMolecularWeight(string)`: average-isotopic protein Mw (Da).
- `SequenceStatistics.CalculateNucleotideMolecularWeight(string, bool)`: average-isotopic DNA/RNA Mw (Da).

### 5.2 Current Behavior

Single forward pass; the water correction `(n − 1) · W` is applied once after accumulation.
Constants are named (`AverageWaterMass`, `AminoAcidWeights`, `DnaNucleotideWeights`,
`RnaNucleotideWeights`) and source-cited inline. Unknown symbols are silently skipped so every
reported mass derives only from cited monomer masses.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Protein Mw = Σ average amino-acid masses + one water [1], realized as Σ free-aa − (n−1)·W [4].
- Nucleic-acid Mw = Σ monophosphate masses − (n−1)·W [4].
- Average mass tables and water = 18.0153 Da copied from Biopython IUPACData / SeqUtils [4][5].

**Intentionally simplified:**

- Non-standard / ambiguous symbols: skipped rather than rejected; **consequence:** an
  ambiguous input yields the mass of its recognized monomers only, instead of an error.

**Not implemented:**

- Monoisotopic masses; double-stranded and circular nucleic-acid corrections; **users should
  rely on:** Biopython `Bio.SeqUtils.molecular_weight` (monoisotopic/double_stranded/circular flags) — no current in-repo alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Unknown symbols skipped | Deviation | Differs from Biopython reject-on-unknown; no invented mass | accepted | Reported mass uses only cited monomers |
| 2 | Average-only mass set | Assumption | Result is average, not monoisotopic | accepted | Matches Expasy/ProtParam default |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | 0 | No monomers (INV-01) |
| single monomer | free monomer mass (zero bonds) | (n−1)·W = 0 [1][4] |
| lowercase input | same as uppercase | Input upper-cased (INV-04) |
| unknown symbol | contributes no mass / no bond | §5.4 deviation |

### 6.2 Limitations

Average masses only (not monoisotopic); single-stranded only (no double_stranded/circular
handling); ambiguous/modified residues and non-standard nucleotides are ignored, not modeled.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double protein = SequenceStatistics.CalculateMolecularWeight("AGC");          // 249.2874 Da
double dna     = SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: true);  // 949.6095 Da
double rna     = SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: false); // 997.6177 Da
```

**Numerical walk-through (DNA "AGC"):**
331.2218 (A) + 347.2212 (G) + 307.1971 (C) − 2 × 18.0153 = 949.6095 Da, rounding to the
Biopython docstring value 949.61 [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateMolecularWeight_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateMolecularWeight_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-MW-001-Evidence.md](../../../docs/Evidence/SEQ-MW-001-Evidence.md)

## 8. References

1. Gasteiger E., Hoogland C., Gattiker A., et al. 2005. Protein Identification and Analysis Tools on the ExPASy Server. *The Proteomics Protocols Handbook*, Humana Press, 571–607. Compute pI/Mw documentation: https://web.expasy.org/compute_pi/pi_tool-doc.html
2. Expasy FindMod. Average masses of amino acid residues. SIB Swiss Institute of Bioinformatics. https://web.expasy.org/findmod/findmod_masses.html
3. Expasy ProtParam documentation. SIB Swiss Institute of Bioinformatics. https://web.expasy.org/protparam/protparam-doc.html
4. Cock P.J.A., Antao T., Chang J.T., et al. 2009. Biopython. *Bioinformatics* 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 — `Bio/SeqUtils/__init__.py` (`molecular_weight`): https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py
5. Biopython `Bio/Data/IUPACData.py` (`protein_weights`, `unambiguous_dna_weights`, `unambiguous_rna_weights`), mass data from PubChem. https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py
