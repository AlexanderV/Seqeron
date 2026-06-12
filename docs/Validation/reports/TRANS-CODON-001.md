# Validation Report: TRANS-CODON-001 — Genetic code / codon tables

- **Validated:** 2026-06-12   **Area:** Translation
- **Canonical method(s):** `GeneticCode.Translate(codon)`, `IsStartCodon`, `IsStopCodon`, `GetCodonsForAminoAcid`, `GetByTableNumber`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect found)

## Stage A — Description

### Sources opened & what they confirm
- **NCBI Genetic Codes** (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — the definitive
  source. Fetched the verbatim `AAs`, `Starts`, `Base1`, `Base2`, `Base3` strings for tables 1, 2, 3, 11:

  | Tbl | AAs |
  |-----|-----|
  | 1  | `FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG` |
  | 2  | `FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG` |
  | 3  | `FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG` |
  | 11 | `FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG` |

  `Starts`:
  - T1:  `---M------**--*----M---------------M----------------------------`
  - T2:  `----------**--------------------MMMM----------**---M------------`
  - T3:  `----------**----------------------MM---------------M------------`
  - T11: `---M------**--*----M------------MMMM---------------M------------`

  Base1/2/3 are identical across all tables (the codon order is `T,C,A,G` in each position).

### Decode check (codon[i] = Base1[i]Base2[i]Base3[i], T→U)
Decoding the NCBI strings programmatically yields, for **Table 1**:
- 64 codons → 20 amino acids + stop.
- Stops: **UAA, UAG, UGA** (positions where AAs=`*`).
- Starts: **AUG, UUG, CUG** (positions where Starts=`M`).

Alternative-table differences from Standard, decoded from NCBI directly:
- **Table 2 (Vert. Mito):** AUA→M (was I), UGA→W (was *), AGA→*, AGG→* (were R). Stops {UAA,UAG,AGA,AGG}; Starts {AUU,AUC,AUA,AUG,GUG}.
- **Table 3 (Yeast Mito):** CUU/CUC/CUA/CUG→T (were L), AUA→M (was I), UGA→W (was *). Stops {UAA,UAG}; Starts {AUA,AUG,GUG}.
- **Table 11 (Bacterial/Plastid):** codon table identical to Standard; Starts {UUG,CUG,AUU,AUC,AUA,AUG,GUG}; Stops {UAA,UAG,UGA}.

### Edge-case semantics
TestSpec / Evidence define: length≠3 → ArgumentException; null/empty → ArgumentException;
non-IUPAC nucleotide → ArgumentException; IUPAC ambiguity codon (N/R/Y…) → 'X' (Biopython/EMBOSS
convention); DNA T ↔ RNA U normalized; case-insensitive; stop → '*'. All standard and sourced.

### Findings / divergences
None. The Evidence doc's `AAs`/`Starts` strings and per-table difference tables reproduce the
NCBI strings exactly. NCBI official name for Table 11 ("Bacterial, Archaeal and Plant Plastid")
matches the code's `Name`.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Core/GeneticCode.cs`
- `CreateStandardCodonTable()` (lines 230–320): the 64-entry literal dictionary.
- `CreateVertebrateMitochondrial/YeastMitochondrial/BacterialPlastid` (174–225): start from the
  standard table and apply NCBI overrides + start/stop sets.
- `Translate` (63–79): length guard → ToUpper + T→U → table lookup → IUPAC→'X' → throw.

### Codon tables match NCBI EXACTLY (programmatic, all 64 codons)
Parsed the 64 literal `["XXX"] = 'Y'` entries from the source and compared against the
NCBI-decoded Table 1 map: **identical, 0 diffs**. Recomputed the override-applied maps for
tables 2/3/11 vs the NCBI-decoded maps: **all 3 identical, 0 diffs**. Start/stop sets in the
code compared to the NCBI `Starts`/`AAs`-derived sets: **all 4 tables match exactly** (counts
3/3, 5/4, 3/2, 7/3 for starts/stops of T1/T2/T3/T11).

### Test quality audit (`GeneticCodeTests.cs`, 64 tests)
- `Translate_CompleteStandardCodonTable_MatchesNcbi` hard-codes all 64 expected AAs; I verified
  the hard-coded expectations equal the NCBI-decoded values (not merely self-consistent with the code).
- Alternative-table differences each asserted (UGA=W, AGA/AGG=*, AUA=M; CUx=T; etc.) with exact
  start/stop set counts and members.
- Edge cases covered with real assertions: length≠3, empty, null, invalid nucleotide (XYZ,12G),
  ambiguous→X (ANN,NNN,NAA,NCC), mixed/lower case, DNA-format starts/stops, degeneracy
  (M=1, W=1, L/S/R=6, I=3, stop=3), `GetByTableNumber` valid/invalid.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.** No defect found; no code or test changes.
- Build: `Seqeron.Genomics.Tests` builds with 0 warnings/0 errors.
- Tests: GeneticCode filter = 64 passed / 0 failed; full suite = **4484 passed, 0 failed** (baseline held).
- Scope note: library implements 4 of NCBI's 33 tables (1, 2, 3, 11); `GetByTableNumber` correctly
  throws for unsupported numbers. This is a documented, bounded scope — not a defect.
