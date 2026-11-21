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

**Also implemented (audit round 1, A1-9, B09 F42/F43):** BLAST's integer-rounded percentages and BLAST positives (§5.5) and statistics of an alignment of more than two sequences (§5.6: EMBOSS infoalign / distmat, ClustalW conservation line).

**Named built-in matrices (audit A2-15, B09 F90).** `SubstitutionMatrix.Get(name)` / `TryGet` / `BuiltInNames` (`SubstitutionMatrix.Registry.cs`, data `SubstitutionMatrix.EmbossData.cs` generated by `scripts/alignment/generate_emboss_matrices.py`): every matrix of the EMBOSS 6.6.0 data directory, 70 of them - `EBLOSUM30` ... `EBLOSUM90` (13), `EBLOSUM62-12`, `EBLOSUMN`, `EPAM10` ... `EPAM500` (50), `EDNAFULL`, `EDNAMAT`, `EDNASIMPLE`, `ENUC.4.4` - by EMBOSS file name or NCBI name without the `E` (`BLOSUM45`, `BLOSUM50`, `BLOSUM80`, `PAM30`, `PAM70`, `PAM250`, `NUC.4.4` = EDNAFULL as before), cell for cell equal to the files (70 / 70 matrices, `SubstitutionMatrix_Registry_Tests`; the NCBI aliases equal Biopython 1.88 `substitution_matrices` on every common symbol, 0 differing cells for BLOSUM45 / 50 / 62 / 80 / 90, PAM30 / 70 / 250, NUC.4.4). `SubstitutionMatrix.IsNucleotide` separates the nucleotide matrices (the MCP tools `alignment_statistics`, `format_alignment`, the matrix-scored aligners and `multiple_align_clustalw`'s `scoreMatrix` accept every name). `ENUC.4.2` is a rectangular 15 x 4 table (IUPAC rows x A T G C) and not a pairwise matrix; `EGONNET` is not part of the EMBOSS 6.6.0 data directory of the installed program (the Gonnet tables of ClustalW are `SubstitutionMatrix.ClustalWGonnet250` / `ClustalWMatrixFile` series). Any other table can still be supplied through the public constructor.

**EMBOSS `distmat`: all corrections and options (audit A2-15, B09 F90).** `MultipleAlignmentStatistics.CalculateDistanceMatrix(rows, MsaDistanceOptions)` (`MultipleAlignmentStatistics.Distmat.cs`) ports `emboss/distmat.c` function by function with its single-precision arithmetic: `-nucmethod` 0-5 (uncorrected, Jukes-Cantor, Kimura 2-parameter, Tamura, Tajima-Nei, Jin-Nei Gamma with `-calculatea` / `-parametera`), `-protmethod` 0-2 (Kimura Protein), `-gapweight`, `-ambiguous` (`distmat_checkambigProt`; nucleotides: `ajBaseAlphacharCompare` with its own defects) and `-position` 123 / 12 / 1 / 2 / 3; rows are typed one by one as the EMBOSS reader does (a nucleotide-looking row has X turned into N). 1 100 random alignments (DNA with IUPAC codes, protein with X/B/Z/J; every method and option) equal the printed values of the installed `distmat` string for string. **Reproduced defects of the program:** the Jin-Nei `slen` is not reset for a pair without scored positions; `-ambiguous` on nucleotides never matches an unambiguous base (it examines all but the last alternative of the first code) and does not terminate when the first code is not an IUB letter, e.g. a gap (`NotSupportedException`); `-position` other than 1, 2, 3, 12, 123 and 1 / 2 / 3 for proteins are fatal (`ArgumentException`). **Not reproduced:** Tajima-Nei and the protein Kimura method with a position other than 123 loop to the alignment length instead of the position-selected length and read heap memory beyond their buffers (the program prints values that depend on the allocator; about half of 100 sampled runs equal the NUL-padded model): `NotSupportedException`. The earlier overload `CalculateDistanceMatrix(rows, method, residueType, gapWeight)` delegates to it, so it is now single precision as well (a hand value 10.0 is 10.000002); the printed two decimals are the sourced values. MCP `alignment_distance_matrix` takes `method` (`uncorrected`, `jukes-cantor`, `kimura`, `tamura`, `tajima-nei`, `jin-nei`), `ambiguous`, `position`, `calculateA`, `parameterA`.

**ClustalW's `-type=auto` guess** (the residue type of `clustalw` when `-TYPE` is omitted) is `SequenceAligner.GuessClustalWSequenceType` (audit A2-14, F92: the SECOND sequence decides, DNA at >= 85 % A / C / G / T / U / N; equal to the program on 500 oracle sets; the note above that 4 of 151 random DNA alignments were classified protein is this rule at work). The column statistics of this page keep an explicit residue type or the EMBOSS first-row rule.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Prior Similarity = `(matches+mismatches)/L` | Deviation | Over-reported similarity | fixed | ALIGN-STATS-001 (2026-06) |
| 2 | Residues compared case-sensitively | Deviation | `a` vs `A` counted as mismatch | fixed | B09 F11 |
| 3 | Unequal rows not rejected | Deviation | exception / wrong counts | fixed | B09 F12 |
| 4 | No gap-opening count, no positives | Missing | BLAST/Biopython parity | added | B09 F13 |
| 5 | Similarity only from a scalar model; markup `' '` for mismatches, `.` declared unreachable; identical X/X counted similar | Deviation | protein/IUPAC similarity unavailable, markup differed from needle | fixed | B09 F14 (SubstitutionMatrix, EMBOSS classifier) |

### 5.5 BLAST-style statistics (A1-9, B09 F42)

`AlignmentStatistics` carries the BLAST+ report values (NCBI C++ toolkit `src/objtools/align_format/align_format_util.cpp` and `showalign.cpp`, opened 2026-11-18):

- `BlastPercent(n, d)` = `CAlignFormatUtil::GetPercentMatch`: an **integer**, `(int)(0.5 + 100·n/d)` (round half up) but **never 100 unless n = d** (an alignment with mismatches is at most 99%); used for `Identities = a/L (p%)`, `Positives = b/L (q%)` and `Gaps = g/L (r%)`. Properties `BlastPercentIdentity`, `BlastPercentPositives`, `BlastPercentGaps`; `ToBlastString(protein)` prints the pairwise-report line; `BlastTabularPidentText` is the `-outfmt 6` `pident` (100·Matches/L, three decimals in BLAST+ 2.12; "100.000" only when all columns match).
- `BlastPositives` = identities + non-identical positive-score columns (BLAST+ `x_FillIdentityInfo`: `match + '+'`); it differs from the EMBOSS `Positives` (§2.2) only for an identical pair whose own score is not positive (X/X in BLOSUM62), which BLAST counts as an identity.
- Oracle (`scripts/alignment/generate_blast_stats_oracle.py`, `TestData/Alignment/blast_stats_oracle.tsv`): BLAST+ 2.12.0 `blastn` (megablast and blastn tasks, `-dust no`) and `blastp` (`-comp_based_stats 0 -seg no`) on 320 random homologous pairs, best HSP both as `-outfmt "6 qseq sseq length nident mismatch gapopen gaps positive pident"` and as pairwise text (`-outfmt 0`). Seqeron's `Matches`, `Mismatches`, `GapOpens`, `Gaps`, `AlignmentLength`, `BlastPositives`, the three integer percentages, `pident` text and the whole report line equal BLAST on all 320 (24 rows exercise the "99, not 100" cap: pident 99.5-99.99).

### 5.6 Multiple-alignment statistics (A1-9, B09 F43)

`MultipleAlignmentStatistics` ([MultipleAlignmentStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs), project Analysis because it calls the EMBOSS-`cons` consensus `MotifFinder.GenerateEmbossConsensus` — not re-implemented) works on gapped rows (upper-cased, `.`/`~` read as `-`; unequal rows rejected):

- `CalculateSequenceStatistics(rows, referenceIndex?, residueType)` = EMBOSS `infoalign` (`infoalign_Compare`, `emboss/infoalign.c`, opened): per row against the reference row or the consensus (plurality 50 % of the total weight, identity 0): leading/trailing gaps of the row are ignored; `SeqLen = Ident + Similar + Differ`; `AlignLen` = columns between the row's first and last residue; `Gaps` = internal gap runs, `GapLen` = internal gap positions; a residue is Identical (case-insensitive equal), else Similar when the EDNAFULL/EBLOSUM62 score with the reference symbol is > 0, else Different (also against a reference gap); `% Change = 100·(AlignLen − Ident)/AlignLen`. Oracle: infoalign 6.6.0 `-refseq k` and `-refseq 0` on 150 random DNA/IUPAC and protein MSAs (3-8 rows, end/internal gaps) = 3 268 per-row records (600 runs, 219 against the consensus), all equal (percent within 1e-4: EMBOSS prints a float).
- `CalculateDistanceMatrix(rows, method, residueType, gapWeight)` = EMBOSS `distmat` (opened `emboss/distmat.c`): matches = columns with the same non-gap symbol, g = columns where either row has a gap; Uncorrected `D = 1 − m/(L − g + g·gapWeight)`; Jukes-Cantor `−b·ln(1 − D/b)`, b = 3/4 (nucleotide) or 19/20 (protein); saturated/undefined entries are NaN/∞ as in distmat. Oracle: distmat 6.6.0 (`-nucmethod/-protmethod 0|1`, `-gapweight 0|0.5`; two decimals ×100) on the same 150 MSAs = 600 matrices, all equal (the corrections and options of A2-15: 1 100 more, see above) (two decimals; EMBOSS computes in single precision, so near-saturated Jukes-Cantor values agree to 1e-4 relative). `CalculatePairwiseIdentityMatrix(rows, gapWeight)` = 100·(1 − D) (default gap weight 0: identity over the columns where neither row has a gap). Biopython 1.88 `DistanceCalculator('identity')` compares characters with `skip_letters = ()` and divides by the alignment length, i.e. it equals the gap-weight-1 identity plus the both-gap columns (checked on 77 matrices).
- Per-column statistics (`Calculate`): residue counts (gap `-` included), gap count/fraction, most common residue, and the ClustalW conservation symbol (`*` identical without gap, `:` all residues in one strong group STA/NEQK/NHQK/NDEQ/QHRK/MILV/MILF/HY/FYW, `.` in one weak group CSA/ATV/SAG/STNK/STPA/SGND/SNDEQK/NDEQHK/NEQHRK/FVLIM/HFY; nucleotide alignments only `*`; any gap column blank). Oracle: ClustalW 2.1 `-convert -output=clustal -type=dna|protein` conservation line on 150 MSAs. Summary: fully conserved and gap-free column counts, gap percent, mean pairwise identity.
- MCP (Analysis server): `multiple_alignment_statistics` and `alignment_distance_matrix`; the Alignment-server `alignment_statistics` reports the BLAST fields of §5.5.

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

### 7.x EMBOSS `showalign` (A3-13c, B09 F148)

`MultipleAlignmentDisplay.ShowAlign(names, rows, ShowAlignOptions)` / MCP `showalign` renders an alignment as the program does: blocks of `width` columns with the position-number line (every 10th column; relative to the reference residues with `-n` / `+n` flanks when a row is the reference), the tick ruler (`|` every 10th, `:` every 5th, `.` at reference gaps, `V` / `v` flank ticks), the reference row, the other rows converted against the reference (`show` N / I / S / D / A with the case-encodes-similarity rule, `.` for suppressed residues) in `order` I / A / S, and the consensus (the `cons` algorithm with the gap option, case kept) at the bottom or top. 600 / 600 random cases equal `showalign` 6.6.0. `-sbegin` / `-send` (A4-9, `ShowAlignOptions.Begin` / `End`, MCP `sbegin` / `send`) cut the displayed columns (1-based, negative = from the end, beyond the alignment clamps, an end before the begin shows one column) while the consensus, numbers and ticks stay those of the whole alignment: 400 / 400 random ranges equal the program. `-html` / `-highlight` (A3-13c2, F191; `ShowAlignOptions.Html` / `Highlight`, ranges parsed by `EmbossRangeSpec.Parse` = the program's `range` syntax): the display is wrapped in `<pre>` ... `</pre>` and the ranges of every sequence / consensus line are wrapped in `<font color=COLOUR>` ... `</font>` (empty colour = `red`; names, numbers and ticks are not coloured): 300 / 300 texts equal the program. (`consambig` ignores `-sbegin` / `-send`: measured on 300 rows.)

## 8. References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
2. EMBOSS needle application documentation, release 6.6. EMBOSS. https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
3. EMBOSS Alignment Formats (markup legend; its `:`/`.` thresholds differ from what the 6.6.0 binary and source do — the binary and ajalign.c [6] are followed). EMBOSS. https://emboss.sourceforge.net/docs/themes/AlignFormats.html
4. NCBI. BLAST QuickStart — Comparative Genomics. NCBI Bookshelf NBK1734. https://www.ncbi.nlm.nih.gov/books/NBK1734/
5. Pérez-Mejías A (amaurypm). pseqsid — pairwise sequence identity/similarity. GitHub. https://github.com/amaurypm/pseqsid
6. EMBOSS 6.6 source, `ajax/core/ajalign.c` (`alignConsStats`, `alignSim`). https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajalign.c (opened 2026-11-18).
7. Henikoff S, Henikoff JG. 1992. Amino acid substitution matrices from protein blocks. *PNAS* 89:10915–10919. Dayhoff MO, Schwartz RM, Orcutt BC. 1978. A model of evolutionary change in proteins. *Atlas of Protein Sequence and Structure* 5(S3):345–352.
8. Cock PJA et al. 2009. Biopython. *Bioinformatics* 25:1422–1423 (Biopython 1.88 `Bio.Align.Alignment.counts`, `Bio.Align.substitution_matrices`).
9. NCBI C++ toolkit, `src/objtools/align_format/align_format_util.cpp` (`GetPercentMatch`), `showalign.cpp` (`s_DisplayIdentityInfo`, `x_FillIdentityInfo`) and BLAST+ 2.12.0 `blastn`/`blastp` binaries. https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/objtools/align_format/align_format_util.cpp (opened 2026-11-18).
10. EMBOSS 6.6 source `emboss/infoalign.c`, `emboss/distmat.c`; binaries `infoalign`, `distmat` 6.6.0. https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/infoalign.c (opened 2026-11-18).
11. ClustalW 2.1 (`clustalw -convert`, conservation line `*` `:` `.`); Larkin MA et al. 2007. Clustal W and Clustal X version 2.0. *Bioinformatics* 23:2947-2948.
12. Biopython 1.88 `Bio.Phylo.TreeConstruction.DistanceCalculator('identity')` (source inspected locally).
