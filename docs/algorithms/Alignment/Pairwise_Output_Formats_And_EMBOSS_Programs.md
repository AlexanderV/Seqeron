# Pairwise alignment output formats and the EMBOSS programs needle, water and stretcher

**Units:** ALIGN-GLOBAL-001, ALIGN-LOCAL-001, ALIGN-MULTI-001 (audit A2-27 / A2-28, B09 F120 - F122).
**Code:** `EmbossAlignmentWriter`, `EmbossAlignFormat`, `EmbossAlignOptions`, `EmbossPairAlignment`, `EmbossAlignedSequence` (`EmbossAlignmentTypes.cs`), `SequenceAligner.StretcherAlign`, `SequenceAligner.NeedleAlign`, `SequenceAligner.WaterAlign`.
**MCP:** `write_alignment`, `stretcher_align`.

## 1. EMBOSS `-aformat` writers (F120)

`EmbossAlignmentWriter` reproduces the 29 alignment formats of EMBOSS 6.6.0 (`ajax/core/ajalign.c` `alignFormat[]`) byte for byte:

* pairwise layouts (exactly two rows): `pair`, `simple`, `unknown`, `multiple`, `srspair`, `srs`, `markx0`, `markx1`, `markx2`, `markx3`, `markx10`, `match`, `score`, `sam`, `tcoffee`. They print the program header (`Program`, `Rundate`, `Commandline`, `Align_format`, `Report_file`), the sequence list (with `-aaccshow` / `-adesshow` / `-ausashow`), `Matrix`, `Gap_penalty`, `Extend_penalty`, the statistics block (`Length`, `Identity`, `Similarity`, `Gaps`, `Score`, computed in single precision as `ajAlignSetSubStandard` does) and, once per file, the tail;
* sequence-set layouts (any number of rows, the `ajseqwrite.c` writers): `fasta`, `a2m`, `msf` (GCG checksums, `CompCheck`), `clustal` (CLUSTAL W 1.83 layout), `mega`, `meganon`, `nexus` / `paup`, `nexusnon` / `paupnon`, `phylip`, `phylipnon`, `selex`, `treecon`.

Facts established from the source and confirmed by the binaries: the number columns of `markx0` use the `%*d` running widths of `alignWriteMark` (the `num[]` flag persists over blocks); `markx2` shows the second row with `.` where it equals the first; the T-Coffee percent identity is 0 unless the program stored statistics (stretcher does); stretcher prints integer scores and drops the descriptions; EMBOSS prints dates as `Thu 19 Nov 2026 21:43:53` (`%2d` day) in headers and `dd/mm/yy` in MSF / MEGA / NEXUS; the type that MSF / NEXUS / MEGA print comes from the output set (`seqret`) or from the residues (alignment programs: NEXUS always `datatype=DNA`, MEGA `Protein` for every protein-compatible text); `ajStrAssignSubS` clamps a start behind the end to the last character (MEGA blocks of rows of different length).

Verification: `emboss_aformat_oracle.tsv` (`scripts/alignment/generate_emboss_aformat_oracle.py`): needle / water / stretcher, 29 formats x 114 alignment sets (all matrices, gap penalties, `-awidth`, `-aglobal`, `-adesshow`, `-aaccshow`, `-ausashow`, one sequence against up to three, faked run dates through `scripts/alignment/faketime.c`), 3 303 files identical; `emboss_seqout_oracle.tsv` (`generate_emboss_seqout_oracle.py`): `seqret -osformat` of 11 formats on gapped sets of 1 - 9 rows, 1 309 files identical. Not implemented (EMBOSS needs the `cons` consensus): the pair-markup formats for three or more sequences; `debug` (memory trace); `mase` is not an `-aformat`.

## 2. stretcher (F121)

`StretcherAlign` ports `emboss/stretcher.c`: the Myers-Miller (1988) linear-space divide-and-conquer with the program's integer arithmetic, its midpoint rules (`CC + RR >= midc` with the gap-closed test, then the type-2 scan from the right), the edit-script macros and `stretcher_Calcons`. A gap of length k costs `gapopen + (k - 1) * gapextend`; defaults 16 / 4 (nucleotide, EDNAFULL) and 12 / 2 (protein, EBLOSUM62). Oracle `stretcher_oracle.tsv`: 689 / 689 pairs equal (rows and score).

## 3. needle and water (F122)

`NeedleAlign` ports `embAlignPathCalcWithEndGapPenalties` + `embAlignWalkNWMatrixUsingCompass` in single precision (float matrices, gap penalties rounded by `ajRoundFloat` to 16 mantissa bits, so 0.1 becomes 0.10000038), including `-endweight`, `-endopen`, `-endextend` and the traceback tie rules (`U_FEPS` 1.192e-6). `WaterAlign` ports `embAlignPathCalcSW` + `embAlignWalkSWMatrix` (float path, double maxima, the last cell within `U_FEPS` of the maximum, gap walk-back within 0.01). Oracles `needle_oracle.tsv` (887 / 887 pairs, 100+ with end weights) and `water_oracle.tsv` (860 / 860 pairs, rows, starts and scores).

Both programs are not always optimal: their recurrences forbid some gap-state transitions and mis-initialise the boundary of a one-residue sequence (needle prints score 0), so the plain affine optima (`GlobalAlignAffine`, `SemiGlobalAlignAffine(EndGaps.All)`, `LocalAlignAffine`) are upper bounds that the programs reach on more than 90 % of the committed pairs.

## 4. Remaining (registered in B09.md)

Biopython `Bio.Align` writers / readers (A2-27b), `needleall`, `supermatcher` seeds (A2-28b / c), est2genome / wordmatch / dotmatcher (A2-28d), the ksw2 variants (A2-29).
