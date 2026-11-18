# Pairwise Alignment Statistics and Formatting

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-STATS-001 |
| Related Projects | Seqeron.Genomics.Alignment, Seqeron.Genomics.Infrastructure |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 (B09 re-validation: EMBOSS 6.6.0 binary + source, Biopython 1.88) |

## 1. Overview

Given a completed pairwise alignment (two equal-length strings over the alphabet plus the gap symbol `-`), this unit computes the standard summary statistics reported by global/local aligners — **Identity**, **Similarity** (BLAST *positives*), **Gaps** (count and percentage of the alignment length) and the number of **gap openings** — and renders a human-readable three-line alignment with the EMBOSS `needle`/`water` markup line. Residue-pair similarity is decided either by a scalar match/mismatch model (`ScoringMatrix`, DNA) or by a residue-pair `SubstitutionMatrix` (built-in **BLOSUM62**, **PAM250**, **EDNAFULL**, cell-for-cell equal to the EMBOSS 6.6 data files and Biopython 1.88 `substitution_matrices`). The statistics and markup reproduce the EMBOSS `needle`/`water` 6.6.0 binary and `ajalign.c` source on a 301-alignment oracle (12 439 aligned columns) and Biopython 1.88 `Alignment.counts` for mismatches / gap openings. Exact, deterministic, O(L).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

After a dynamic-programming aligner produces two gapped sequences of equal length L, biologists summarise the relationship with three percentages over the L aligned columns. The denominator L **includes gap columns**, which is the EMBOSS `needle` convention [2] and the convention corroborated for percent identity generally [secondary corroboration]. These summaries let users judge relatedness independently of the raw score.

### 2.2 Core Model

For an alignment of two rows `a`, `b` of equal length L (residues compared **case-insensitively**, as EMBOSS and the aligners of this class upper-case their input), each column `i` is classified:

- **gap column** if `a[i] = '-'` or `b[i] = '-'` (EMBOSS "one or more sequences with a gap"; a both-gap column, impossible in a pairwise alignment, would count the same way);
- **identical** if not a gap and `a[i] = b[i]`;
- **mismatch** otherwise.

`Matches + Mismatches + Gaps = L`. A non-gap column is **similar** (a *positive*) when its substitution score is positive. EMBOSS `alignConsStats` (ajalign.c, EMBOSS 6.6.0; plurality 51 % ⇒ threshold 1.02 for two sequences) makes this exact: with `M` the matrix and `self(x) = [M(x,x) > 0]`, a column `(x, y)` is similar iff `M(x,y) > 0 ∧ (self(x) ∨ self(y))`; for BLOSUM62, PAM250 and EDNAFULL this is the same as `M(x,y) > 0`, the BLAST/Biopython `positives` definition. In particular an **identical** column is similar only if its own score is positive (X/X in BLOSUM62 = −1, N/N in EDNAFULL = −1 are identities but not positives). Then [2]:

- Identity% = `Matches / L × 100`, Similarity% = `Positives / L × 100`, Gaps% = `Gaps / L × 100` — the denominator is always the alignment **Length including gap columns** (EMBOSS needle; **not** the shorter sequence, not the number of non-gap columns, not BLAST's per-HSP rounding);
- **GapOpens** = number of maximal `-` runs, counted per row (BLAST tabular `gapopen`, Biopython `open_gaps`); independent of any gap scoring. A run in one row directly followed by a run in the other counts as two.

The published EMBOSS worked example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62, gap open 10 / extend 0.5) gives L = 149, Identity 65 (43.6 %), Similarity 90 (60.4 %), Gaps 9 (6.0 %) [2]; it is reproduced exactly through the real BLOSUM62 table (and its 4 gap runs equal Biopython `open_gaps`).

**Markup line** (EMBOSS `needle`/`water` srspair, `alignConsStats` + `alignSim`): with `h = self(x) + [M(x,y) > 0]` (a symbol outside the matrix alphabet contributes 0): `h = 2` → `|` if `x = y` else `:`; `h = 1` → `.` (aligned, not similar); `h = 0` → space; gap column → space. Consequences, all confirmed with the needle 6.6.0 binary: L/I, A/S, K/R in BLOSUM62 → `:`; A/T (0) → `.`; X/X → space (identity without a bar). EMBOSS quirk reproduced verbatim: when the **first** row's residue has a non-positive self score (an IUPAC ambiguity code, or X) the column is marked `.` (pair positive) or space (pair non-positive) although the header may count it as similar (R/A → `.` but A/R → `:`; X/P → space but P/X → `.`). For the scalar `ScoringMatrix` model (no residue-pair information): identical → `|` iff `Match > 0`, non-identical → `:` iff `Mismatch > 0`, else `.`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Matches + Mismatches + Gaps = AlignmentLength` | The three classes partition every column [2][4] |
| INV-02 | `Identity% ≤ Similarity% ≤ 100%` for a model with positive identical scores | Positive-score set ⊇ identical set [2][5] |
| INV-03 | `Identity% + mismatch% + Gap% = 100%` | Length-denominator partition [2] |
| INV-04 | Markup line length equals the displayed sequence-block length; characters ∈ {`|`, `:`, `.`, space} | srspair markup (ajalign.c `alignSim`) |
| INV-05 | `GapOpens ≤ Gaps`; `GapOpens = 0` iff `Gaps = 0` | A run has ≥ 1 column |

### 2.5 Comparison with Related Methods

| Aspect | This unit (EMBOSS needle) | BLAST HSP | Biopython `Alignment.counts` |
|--------|---------------------------|-----------|------------------------------|
| Percent denominator | alignment length incl. gap columns | alignment length incl. gaps (identities/length) | — (counts only) |
| Similar column | M(x,y) > 0 (and a positive self score) | positive score | `positives` = M(x,y) > 0 |
| Gap count | gap *columns* (= total gap length) | total gap length; `gapopen` = runs | `gaps` total, `open_gaps` runs |
| Mismatches | non-gap non-identical columns | `mismatch` | `mismatches` |

Agreement with Biopython on `mismatches`, `open_gaps` and `identities`, and with EMBOSS on Length/Identity/Similarity/Gaps, is checked on the oracle table.

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `alignment` | `AlignmentResult` | required | Completed alignment (two equal-length gapped rows) | non-null, rows equal length |
| `scoring` | `ScoringMatrix?` | `SimpleDna` | Scalar model (overloads without a matrix): similar iff `Match > 0` (identical) / `Mismatch > 0` (non-identical) | — |
| `matrix` | `SubstitutionMatrix` | — | Residue-pair model (overloads with a matrix): `SubstitutionMatrix.Blosum62` / `Pam250` / `EdnaFull` or a user table | non-null |
| `lineWidth` | `int` | 60 | Residues per display block (FormatAlignment only) | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Matches` | int | Identical aligned residue pairs (EMBOSS identity count) |
| `Mismatches` | int | Non-gap, non-identical columns |
| `Gaps` | int | Columns with a gap on either row |
| `AlignmentLength` | int | L = total columns (includes gaps); denominator of all percentages |
| `Identity` | double | `Matches/L × 100` |
| `Similarity` | double | `Positives/L × 100` |
| `GapPercent` | double | `Gaps/L × 100` |
| `Positives` | int | Similar columns (identical or not) |
| `GapOpens` | int | Gap runs (BLAST `gapopen`, Biopython `open_gaps`) |
| (FormatAlignment) | string | Three-line blocks: row1, markup (`|` `:` `.` space), row2, blank line; `\n` newlines |

`SubstitutionMatrix` tables are symmetric, lookups case-insensitive; BLOSUM62/PAM250: `ARNDCQEGHILKMFPSTWYVBZX*`, EDNAFULL: `ATGCSWRYKMBVHDNU`. Every cell equals EMBOSS `EBLOSUM62` / `EPAM250` / `EDNAFULL` (data files of EMBOSS 6.6.0) and Biopython 1.88 `BLOSUM62` / `PAM250` / `NUC.4.4` (+ U row/column).

### 3.3 Preconditions and Validation

`alignment` must be non-null (`ArgumentNullException`). Rows of different length throw `ArgumentException` (previously an `IndexOutOfRangeException`/silently ignored tail). An empty alignment returns `AlignmentStatistics.Empty` / `""`. `lineWidth` must be positive (`ArgumentOutOfRangeException`). The gap symbol is `-`; residues are compared case-insensitively.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; return Empty/`""` for an empty alignment.
2. Walk the L columns once: track gap runs per row (`GapOpens`); classify each column (gap / identical / mismatch) on upper-cased residues; ask the model (scalar or `SubstitutionMatrix`) whether the column is similar.
3. Percentages = count / L × 100 (gap-inclusive L).
4. (FormatAlignment) emit the markup character of the model's column classification (§2.2), wrapping every `lineWidth` columns.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

`SubstitutionMatrix` stores a dense `int[n,n]` with a 128-entry index (O(1) lookup). BLOSUM62: Henikoff & Henikoff 1992; PAM250: Dayhoff, Schwartz & Orcutt 1978 (as distributed by NCBI/EMBOSS); EDNAFULL: Todd Lowe's IUPAC nucleotide matrix (NCBI NUC.4.4, match 5 / mismatch −4). The gap-convention note on `ScoringMatrix` records that the affine aligners charge `GapOpen + k·GapExtend` (BLAST+/Biopython) whereas EMBOSS charges `gapopen + (k−1)·gapextend`; statistics do not depend on gap scoring.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateStatistics | O(L) | O(1) | single pass over L columns |
| FormatAlignment | O(L) | O(L) | builds the output string |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs); types in [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs).

- `SequenceAligner.CalculateStatistics(AlignmentResult, ScoringMatrix?)` and `CalculateStatistics(AlignmentResult, SubstitutionMatrix)`.
- `SequenceAligner.FormatAlignment(AlignmentResult, int, ScoringMatrix?)` and `FormatAlignment(AlignmentResult, int, SubstitutionMatrix)`.
- `SubstitutionMatrix.Blosum62 / Pam250 / EdnaFull`, `TryGetScore`, `GetScore`.
- MCP: `alignment_statistics` and `format_alignment` accept `matrix` (BLOSUM62 | PAM250 | EDNAFULL) or the scalar `match` / `mismatch` and delegate to the methods above; the result carries `positives` and `gapOpens`.

### 5.2 Current Behavior

Both overload pairs share one column classifier core (`CalculateStatisticsCore`, `FormatAlignmentCore`), so statistics and markup cannot disagree. This unit performs no search/matching, so the repository suffix tree is **not applicable** — a single linear pass over an already-computed alignment.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited sources):**

- Identity / Similarity / Gaps with gap-inclusive denominator (needle docs [2]; confirmed by the needle/water 6.6.0 binary on 301 alignments).
- Similarity per EMBOSS `alignConsStats` (ajalign.c [6]) incl. the identical-with-non-positive-self-score rule and the markup quirk above.
- `mismatches`, `gap openings`, identities agree with Biopython 1.88 `Alignment.counts` on all 301 alignments.
- BLOSUM62 / PAM250 / EDNAFULL cell-for-cell.

**Intentionally simplified:** none for the counted statistics. The scalar `ScoringMatrix` overloads cannot express graded residue-pair scores (they decide similarity per class, identical vs non-identical); for protein or IUPAC DNA use the `SubstitutionMatrix` overloads.

**Not implemented:** residue-pair matrices other than the three built-ins (any other table can be supplied via the public `SubstitutionMatrix` constructor); BLAST's integer-rounded percent identity; EMBOSS's multi-sequence (>2) consensus statistics.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Prior Similarity = `(matches+mismatches)/L` | Deviation | Over-reported similarity | fixed | ALIGN-STATS-001 (2026-06) |
| 2 | Residues compared case-sensitively | Deviation | `a` vs `A` counted as mismatch | fixed | B09 F11 |
| 3 | Unequal rows not rejected | Deviation | exception / wrong counts | fixed | B09 F12 |
| 4 | No gap-opening count, no positives | Missing | BLAST/Biopython parity | added | B09 F13 |
| 5 | Similarity only from a scalar model; markup `' '` for mismatches, `.` declared unreachable; identical X/X counted similar | Deviation | protein/IUPAC similarity unavailable, markup differed from needle | fixed | B09 F14 (SubstitutionMatrix, EMBOSS classifier) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty alignment | `AlignmentStatistics.Empty` / `""` | Denominator undefined |
| Null alignment | `ArgumentNullException` | Validation |
| Unequal rows | `ArgumentException` | Not an alignment |
| `lineWidth ≤ 0` | `ArgumentOutOfRangeException` | Validation |
| All-gap column run | counted as gaps; Identity 0 % if no identical column | Gap classification [2][4] |
| Lower-case residues | same as upper case | EMBOSS case-insensitive |
| Symbol outside the matrix alphabet | counted identical/mismatch by character, never similar (EMBOSS "Bad character" → 0) | ajalign.c |
| Identical pair with non-positive score (X/X, N/N) | identity, not similar, markup space | ajalign.c [6] |
| Perfect identity | Identity = Similarity = 100 %, Gap 0 % (positive Match) | Partition invariant |

### 6.2 Limitations

Statistics assume `-` as the sole gap symbol. A gap run in one row immediately followed by a run in the other row counts as two openings (Biopython `open_gaps` convention). The markup reproduces EMBOSS for two sequences only.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var aln = SequenceAligner.GlobalAlign("ACGTACGT", "ACCTAACGT", SequenceAligner.SimpleDna);
var stats = SequenceAligner.CalculateStatistics(aln);            // SimpleDna => Similarity == Identity
string view = SequenceAligner.FormatAlignment(aln, lineWidth: 60);

// Protein: EMBOSS needle HBA_HUMAN vs HBB_HUMAN rows => Length 149, Matches 65, Positives 90, Gaps 9, GapOpens 4
var prot = SequenceAligner.CalculateStatistics(proteinAlignment, SubstitutionMatrix.Blosum62);
```

**Numerical / biological walk-through (optional):**

EMBOSS HBA vs HBB [2] (reproduced with the real BLOSUM62 table): L = 149, M = 65 ⇒ Identity 65/149 = 43.624…% ≈ 43.6%; similar columns = 90 ⇒ Similarity 90/149 = 60.402…% ≈ 60.4%; G = 9 ⇒ Gaps 9/149 = 6.040…% ≈ 6.0%.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAligner_CalculateStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_CalculateStatistics_Tests.cs) — covers `INV-01`–`INV-05`, the 301-alignment oracle `TestData/Alignment/needle_water_stats.tsv` (generator `scripts/alignment/generate_alignment_stats_oracle.py`)
- Evidence: [ALIGN-STATS-001-Evidence.md](../../../docs/Evidence/ALIGN-STATS-001-Evidence.md)
- Related algorithms: [Global_Alignment_Needleman_Wunsch](../Alignment/Global_Alignment_Needleman_Wunsch.md)

## 8. References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
2. EMBOSS needle application documentation, release 6.6. EMBOSS. https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
3. EMBOSS Alignment Formats (markup legend; its `:`/`.` thresholds differ from what the 6.6.0 binary and source do — the binary and ajalign.c [6] are followed). EMBOSS. https://emboss.sourceforge.net/docs/themes/AlignFormats.html
4. NCBI. BLAST QuickStart — Comparative Genomics. NCBI Bookshelf NBK1734. https://www.ncbi.nlm.nih.gov/books/NBK1734/
5. Pérez-Mejías A (amaurypm). pseqsid — pairwise sequence identity/similarity. GitHub. https://github.com/amaurypm/pseqsid
6. EMBOSS 6.6 source, `ajax/core/ajalign.c` (`alignConsStats`, `alignSim`). https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajalign.c (opened 2026-11-18).
7. Henikoff S, Henikoff JG. 1992. Amino acid substitution matrices from protein blocks. *PNAS* 89:10915–10919. Dayhoff MO, Schwartz RM, Orcutt BC. 1978. A model of evolutionary change in proteins. *Atlas of Protein Sequence and Structure* 5(S3):345–352.
8. Cock PJA et al. 2009. Biopython. *Bioinformatics* 25:1422–1423 (Biopython 1.88 `Bio.Align.Alignment.counts`, `Bio.Align.substitution_matrices`).
