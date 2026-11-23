# Named substitution matrices (Biopython 1.88 parity), fractional scores, text parser

> Audit R14-8 (B09 F301). Code: `SubstitutionMatrix` (`AlignmentTypes.cs`, `.Registry.cs`, `.Parse.cs`, generated `.BiopythonData.cs`), `CodonSubstitutionMatrix`; the pairwise aligners (`SequenceAligner.cs` `PairScoring.Build`).

## The names

`SubstitutionMatrix.Get(name)` / `TryGet` resolve (case-insensitive) the EMBOSS 6.6.0 data matrices (EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, EDNAFULL, ..., with the NCBI names BLOSUM62, PAM250, NUC.4.4, ...) and, since R14-8, **every one of the 30 matrices of Biopython 1.88 `Bio.Align.substitution_matrices.load()`**:

| Kind | Matrices |
|------|----------|
| integer, from the EMBOSS / NCBI set | BLOSUM45, 50, 62, 80, 90, PAM30, 70, 250, NUC.4.4 |
| integer, new (`BiopythonNames`) | BLASTN, BLASTP, FENG, HOXD70, LEVIN, MCLACHLAN, MEGABLAST, RAO, STR, TRANS |
| fractional, new | BENNER6, BENNER22, BENNER74, DAYHOFF, GENETIC, GONNET1992, JOHNSON (2 decimals), JONES, MDM78 (2 decimals), RISLER |
| codon matrix | SCHNEIDER (`CodonSubstitutionMatrix.Schneider`: 64 x 64 over the codons, Schneider, Cannarozzi & Gonnet 2005) |

The cells are the decimal text of Biopython's data files (`scripts/alignment/generate_biopython_matrices.py`): **every cell of the 29 residue matrices equals Biopython 1.88**, and **812 optimal global / local scores** with integer and fractional affine gap scores equal `PairwiseAligner` (`scripts/alignment/generate_matrix_oracle.py`, `biopython_matrices_oracle.tsv.gz`).

## Fractional scores

A matrix constructed from `decimal[,]` stores every cell exactly as an integer in units of 1 / `Denominator` (the smallest power of ten, at most 10^6, that makes every cell integral: 10 for GONNET1992, 100 for JOHNSON). The real score is `GetValue` / `TryGetValue` / `ValueAt`; `IsPositive(a, b)` is the "similar pair" test of the statistics and the report marks. The integer accessors (`GetScore`, `TryGetScore`, `ScoreAt`, `MaxAbsScore`) belong to integer matrices and throw `InvalidOperationException` for a fractional one instead of returning a wrong number. The pairwise aligners (`GlobalAlign` / `LocalAlign` / `SemiGlobalAlign` and the affine, linear-space, general-gap, co-optimal and top-N ones with a `MatrixScoring`) score a fractional matrix exactly: the cells join the gap scores in the search for the smallest common integer unit (the same exact scaling that already served fractional gap scores), `ExactScore` is the real value. The EMBOSS program ports (needle, water, stretcher, matcher, ...) and ClustalW read integer matrices only, as the programs do, and refuse a fractional one.

## Residues U, O, J (parity decided by the oracle)

Measured with the installed `blastp 2.12.0 -ungapped -matrix BLOSUM62` on pairs that differ in one residue: U (selenocysteine) scores as C (U/U 9, U/C 9, U/A 0; `MKTAUYIAK` / `MKTAAYIAK` = 39), O (pyrrolysine) as X (-1 with every residue, -4 with `*`), J (leucine / isoleucine) as the NCBI J column (J/I 3, J/L 3, J/M 2, J/A -1). The **BLASTP** matrix of Biopython reproduces exactly these values (27 symbols `ABCDEFGHIJKLMNOPQRSTUVWXYZ*`). The named BLOSUM / PAM matrices stay as Biopython and EMBOSS have them (no U, O, J: a sequence with them is refused, as `Bio.Align` raises `KeyError`): choose `BLASTP` to score like blastp.

## Parse(text)

`SubstitutionMatrix.Parse(text, name)` reads the NCBI / EMBOSS / Biopython data-file format: `#` comment lines, a header line of single-character symbols, one row per symbol (symbol, then one score per column; integers or decimals; rows in any order; lower case folded). It refuses an empty text, a multi-character symbol (a codon matrix), a duplicated symbol, a missing or short row and a non-number. `ToText()` writes the same format (round trip equal on all 29 residue matrices). MCP: the `matrix` parameter of the pairwise tools takes either a name or the text of a matrix (any value with a line break is parsed).
