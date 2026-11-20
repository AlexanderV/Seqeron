# Evidence Artifact: ALIGN-STATS-001

**Test Unit ID:** ALIGN-STATS-001
**Algorithm:** Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Alignment Formatting
**Date Collected:** 2026-06-13 (re-validated 2026-11-18, B09)

---

## Online Sources

### EMBOSS needle — application documentation (release 6.6)

**URL:** https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — EMBOSS, the standard global aligner)

**Key Extracted Points:**

1. **Identity definition:** "the percentage of identical matches between the two sequences over the reported aligned region (including any gaps in the length)." Identity counts only exact character matches; the denominator is the full alignment length **including gap columns**.
2. **Similarity definition:** "the percentage of matches between the two sequences over the reported aligned region (including any gaps in the length)." Similarity is broader than Identity: it includes positions that are not identical but score positively under the substitution matrix. In the worked example 65 positions are identical but 90 are "matching" — i.e. 25 non-identical positions score positively and are counted as similar.
3. **Gaps definition:** the count and proportion of gap positions within the aligned sequences; denominator is the alignment length.
4. **Score definition:** "the sum of the matches taken from the scoring matrix, minus penalties arising from opening and extending gaps."
5. **Worked example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62):** verbatim output block —
   ```
   # Length: 149
   # Identity:      65/149 (43.6%)
   # Similarity:    90/149 (60.4%)
   # Gaps:           9/149 ( 6.0%)
   # Score: 292.5
   ```
   Confirms: percentage = count / Length × 100, denominator = Length (149) which includes the 9 gap columns. 65/149 = 43.6%, 90/149 = 60.4%, 9/149 = 6.0%.

### EMBOSS — Alignment Formats (AlignFormats / srspair markup legend)

**URL:** https://emboss.sourceforge.net/docs/themes/AlignFormats.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (EMBOSS reference documentation)

**Key Extracted Points:**

1. **Standard (pair / srspair) markup line legend (as printed on the page):** `|` identical; `:` similar (score > 1.0); `.` small positive score; space mismatch or gap. **Superseded 2026-11-18:** the needle/water 6.6.0 binary and `ajalign.c` write `.` for every aligned non-similar column and space only for gap columns (and for identities / first-row ambiguity residues without a positive self score) — see the next section.
2. **Three-line layout:** the markup/consensus line is placed between the two aligned sequence lines.

### NCBI BLAST — pairwise alignment metric definitions

**URL:** (search) "BLAST pairwise alignment identities positives gaps definition NCBI" — corroborating definitions from NCBI/UC-Berkeley BLAST guide and NCBI BLAST QuickStart (https://www.ncbi.nlm.nih.gov/books/NBK1734/)
**Accessed:** 2026-06-13
**Authority rank:** 3 (NCBI BLAST reference output)

**Key Extracted Points:**

1. **Identities:** percentage of residues with a direct (exact) match in the alignment.
2. **Positives (= Similarity):** number of residues that are either identical or have similar chemical properties — operationally, residues that **score positively in the substitution matrix**. Confirms the "positive substitution score ⇒ similar" rule independently of EMBOSS.
3. **Gaps:** dashes (`-`) indicate gaps — positions where one sequence has no counterpart.

### Percent-identity denominator convention (corroboration)

**URL:** (search) "percent identity alignment number identical positions divided by alignment length including gaps definition"
**Accessed:** 2026-06-13
**Authority rank:** 4 (secondary, corroborates the EMBOSS primary definition)

**Key Extracted Points:**

1. **Denominator includes gaps:** "the percent identity is calculated by counting the number of identical aligning residues, dividing by the total length of the aligned region, including gaps in both sequences, and multiplying by 100." Independently confirms that the alignment length (with gap columns) is the denominator.

### pseqsid — reference implementation of sequence identity/similarity

**URL:** https://github.com/amaurypm/pseqsid
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Identity:** "percentage of identical residues in equivalent positions with respect to the sequence length." With the `alignment` length option the denominator is "total alignment length including gap-only columns."
2. **Similarity:** "percentage of identical or similar residues in equivalent positions." A position counts as similar when the residues score positively under the substitution matrix; gap-only columns are ignored from the similar count.


### B09 re-validation (2026-11-18): EMBOSS 6.6 source, needle/water 6.6.0 binaries, Biopython 1.88

**Opened:** `https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajalign.c` (User-Agent `Seqeron-review/1.0`); local Debian `emboss` 6.6.0 (`needle`, `water`, data files `EBLOSUM62`, `EPAM250`, `EDNAFULL`); PyPI `biopython==1.88`. Authority rank 1 (reference implementation source + binary).

**Key Extracted Points (alignConsStats / alignSim, two sequences, weights 1, plurality 51 % ⇒ threshold 1.02):**

1. A column's per-residue "matching" weight is `[M(x,x)>0] + [M(x,y)>0]` (the inner loop includes the sequence itself); `issim` (counted in the header Similarity) iff some residue reaches ≥ 1.02, i.e. `M(x,y) > 0` and a positive self score of x or y; `isident` iff both residues equal. Identity count = `isident` columns (the header Identity counts identical residues regardless of score).
2. Markup is built from the consensus residue of the FIRST row: weight 2 → upper-case (`|` if identical, else lower-case `:`), weight 1 → `#` → `.`, weight 0 → gap char → space. Hence X/X (self −1) → space; IUPAC/X residue in row 1 with non-positive self score → `.` or space even when the pair score is positive (needle: R/A `.`, A/R `:`; header similarity counts both).
3. Gap column: any row `-` → "Gaps" count and a space in the markup.
4. Header: Length = alignment columns incl. gaps; Identity/Similarity/Gaps counts over Length.

**Binary experiments (needle 6.6.0, `-aformat srspair`):** every single-residue pair of BLOSUM62 with W flanks and every IUPAC DNA code pair with A flanks, both orders (EDNAFULL): all markup characters equal the rule in point 2.

**Oracle table** `tests/.../TestData/Alignment/needle_water_stats.tsv` (generator `scripts/alignment/generate_alignment_stats_oracle.py`, seed 20261104, regenerates byte-identical): 301 alignments (needle global + water local; random DNA / IUPAC-DNA with EDNAFULL / protein with EBLOSUM62 and EPAM250 incl. B, Z, X; the published HBA/HBB example as #0) = 12 439 aligned columns. Per row: EMBOSS Length/Identity/Similarity/Gaps and markup line, Biopython 1.88 `Alignment.counts(substitution_matrix)` mismatches and `open_gaps`. Seqeron equals every field.

**Matrix tables:** all 576 (BLOSUM62), 576 (PAM250) and 256 (EDNAFULL) cells equal the EMBOSS 6.6 data files; BLOSUM62 / PAM250 / NUC.4.4 equal Biopython 1.88 `substitution_matrices` on all 576 / 576 / 225 shared cells.

**BLAST / Biopython conventions recorded:** identities/positives/gaps/`gapopen` (BLAST tabular: gap *openings*) and Biopython `identities`, `mismatches`, `positives`, `gaps`, `open_gaps` — same counts as the oracle columns; percent denominators are the alignment length including gaps (EMBOSS and BLAST HSP header), never the shorter sequence.

**Gap-scoring note (not a statistic):** EMBOSS needle/water charge a gap of length k as `gapopen + (k-1)·gapextend`; BLAST+ and this library's affine aligners charge `GapOpen + k·GapExtend`.

---

## Documented Corner Cases and Failure Modes

### From EMBOSS needle / AlignFormats

1. **Gap columns counted in the denominator:** the percentage denominator is always the full alignment Length, so identity + (similar-only) + (mismatch-only) + gap fractions partition 100%.
2. **Similarity ≥ Identity always:** the similar set is a superset of the identical set (identical columns always score positively), so Similarity% ≥ Identity% for any scoring model where identical residues score positively.
3. **DNA simple match/mismatch model:** when the substitution score for a non-identical pair is non-positive (Mismatch ≤ 0, as in all DNA models exposed by this class), no non-identical column is similar, so Similarity equals Identity.

### From BLAST

1. **Empty / no-counterpart positions:** dashes represent gaps; gap columns are not identities and not positives.

---

## Test Datasets

### Dataset: EMBOSS needle HBA_HUMAN vs HBB_HUMAN (published worked example)

**Source:** EMBOSS needle documentation, release 6.6 (Rice, Longden & Bleasby 2000).

| Parameter | Value |
|-----------|-------|
| Alignment Length | 149 |
| Identity (count / %) | 65 / 149 = 43.6% |
| Similarity (count / %) | 90 / 149 = 60.4% |
| Gaps (count / %) | 9 / 149 = 6.0% |
| Score | 292.5 |

> Since 2026-11-18 the library ships BLOSUM62 (`SubstitutionMatrix.Blosum62`) and this protein example is reproduced literally (needle 6.6.0 rows, 65/90/9, 4 gap runs).

### Dataset: Formula reconstruction from the EMBOSS counts

**Source:** derived directly from the EMBOSS counts above; arithmetic only.

| Quantity | Computation | Expected |
|----------|-------------|----------|
| Identity % | 65 / 149 × 100 | 43.624161… |
| Similarity % | 90 / 149 × 100 | 60.402684… |
| Gap % | 9 / 149 × 100 | 6.040268… |

### Dataset: Hand-constructed DNA alignment (deterministic, exact)

**Source:** constructed for this unit; values derived by direct application of the EMBOSS formulas to a known column layout (no implementation output used).

Aligned columns (SimpleDna model: Match +1, Mismatch −1):

```
seq1: A C G T - A C G T
seq2: A C C T A A C - T
```

| Column | c1 | c2 | Class |
|--------|----|----|-------|
| 1 | A | A | identity |
| 2 | C | C | identity |
| 3 | G | C | mismatch (Mismatch −1 ⇒ not similar) |
| 4 | T | T | identity |
| 5 | - | A | gap |
| 6 | A | A | identity |
| 7 | C | C | identity |
| 8 | G | - | gap |
| 9 | T | T | identity |

Length = 9, Matches = 6, Mismatches = 1, Gaps = 2.
Identity = 6/9 × 100 = 66.6̅%; Similarity = 6/9 × 100 = 66.6̅% (no positive non-identical column under SimpleDna); Gap% = 2/9 × 100 = 22.2̅%.
Markup line (needle 6.6.0 style, SimpleDna): `||.` `|` ` ` `||` ` ` `|` → `||.| || |` (mismatch `.`, gap columns space).

### Dataset: Similarity > Identity case (positive-scoring mismatch)

**Source:** constructed; uses a scoring model with Mismatch = +1 so that every non-identical column scores positively and is counted as similar (direct application of the "positive substitution score ⇒ similar" rule from EMBOSS/BLAST/pseqsid).

```
seq1: A C G T
seq2: A C C T   (column 3 non-identical, Mismatch = +1 ⇒ similar)
```

With `ScoringMatrix(Match: 1, Mismatch: 1, GapOpen: 0, GapExtend: -1)`:
Length = 4, Matches = 3, Mismatches = 1, Gaps = 0.
Identity = 3/4 × 100 = 75%; Similarity = (3+1)/4 × 100 = 100%; Gap% = 0%.
Markup line: `||:|` ( `:` marks the similar column).

---

## Assumptions

None. (The 2026-06 assumption that the srspair `.` tier is unreachable and collapses to `:` was refuted by the needle/water 6.6.0 binary and `ajalign.c` on 2026-11-18; the rule is now implemented exactly, see above.)

---

## Recommendations for Test Coverage

1. **MUST Test:** Identity/Similarity/Gaps formula and denominator (count / Length × 100, Length includes gaps) reproduced from the EMBOSS 149-column counts (43.6% / 60.4% / 6.0%). — Evidence: EMBOSS needle worked example.
2. **MUST Test:** DNA SimpleDna alignment ⇒ Similarity equals Identity (no positive non-identical columns). — Evidence: EMBOSS/BLAST "positive substitution score ⇒ similar"; DNA Mismatch < 0.
3. **MUST Test:** Positive-scoring mismatch (Mismatch = +1) ⇒ Similarity > Identity, similar column counted. — Evidence: EMBOSS/BLAST/pseqsid similarity rule.
4. **MUST Test:** Gap counting — column with either side `-` counts as gap, not identity/mismatch. — Evidence: BLAST gap definition; EMBOSS gap count.
5. **MUST Test:** FormatAlignment markup `|`/`:`/`.`/space per needle 6.6.0 (identity, similar, mismatch, gap; X/X and first-row ambiguity quirk). — Evidence: ajalign.c + binary.
6. **SHOULD Test:** Empty alignment ⇒ AlignmentStatistics.Empty / "" — Rationale: documented empty handling, denominator undefined.
7. **SHOULD Test:** Null alignment ⇒ ArgumentNullException; non-positive lineWidth ⇒ ArgumentOutOfRangeException — Rationale: input validation contract.
8. **COULD Test:** FormatAlignment line wrapping at lineWidth (multi-block output) — Rationale: formatting completeness.
9. **COULD Test:** Invariant Identity ≤ Similarity ≤ 100 and fractions partition to 100% (property test) — Rationale: O(n) but cheap structural invariant.

---

## References

1. Rice P, Longden I, Bleasby A (2000). EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
2. EMBOSS needle application documentation (release 6.6). https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html (accessed 2026-06-13)
3. EMBOSS Alignment Formats (markup legend). https://emboss.sourceforge.net/docs/themes/AlignFormats.html (accessed 2026-06-13)
4. NCBI. BLAST QuickStart — Comparative Genomics. NCBI Bookshelf NBK1734. https://www.ncbi.nlm.nih.gov/books/NBK1734/ (accessed 2026-06-13)
5. pseqsid — pairwise sequence identity/similarity reference implementation. https://github.com/amaurypm/pseqsid (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-11-18**: B09 re-validation against the EMBOSS 6.6 source, needle/water 6.6.0 binaries and Biopython 1.88; oracle table, matrix cross-check, corrected markup and similarity rules.

---

## BLAST integer percentages and multiple-alignment statistics (audit A1-9, B09 F42/F43)

- **Opened:** NCBI C++ toolkit `src/objtools/align_format/align_format_util.cpp` (`GetPercentMatch`: `(int)(0.5 + 100n/d)`, 100 only if n == d, else min 99), `showalign.cpp` (`s_DisplayIdentityInfo`: Identities/Positives/Gaps lines; `Positives = positive + match`), `blast_tabular`/`tabular.cpp` (tabular, `pident` 3 decimals checked with the local binary); EMBOSS 6.6 `emboss/infoalign.c` (`infoalign_Compare`, reference = row or `embConsCalc` consensus with plurality 50 %, setcase 0, identity 0) and `emboss/distmat.c` (`distmat_calc_match`, `distmat_uncorrected`, `distmat_JukesCantor`, output x100); ClustalW 2.1 `-convert -output=clustal` (conservation line); Biopython 1.88 `DistanceCalculator._pairwise` (source inspected: identity model `skip_letters = ()`, divides by `len(seq1)`, so both-gap columns count as matches).
- **Numbers:** BLAST+ 2.12.0 blastn + blastp, 320 best HSPs (24 with 99.5 <= pident < 100 shown as 99%): all counts, positives, integer percentages, pident text and report lines equal. infoalign 3 268 per-row records, distmat 600 matrices, ClustalW 151 conservation lines (DNA only `*`; the protein groups STA/NEQK/... and CSA/ATV/... reproduced exactly), Biopython identity 77 matrices.
- **Findings:** EMBOSS `infoalign` ignores a row's leading/trailing gaps, counts every residue facing a reference gap as Different, and `% change = (AlignLen - Ident)/AlignLen` includes similar residues and gap positions; `distmat` divides by `L - g + g*gapweight` (g = columns where either row has a gap) and computes in single precision (near-saturated Jukes-Cantor values differ from double precision by about 1e-5 relative); ClustalW's own `-type=auto` guess is not a simple ACGTUN fraction (4 of 151 random DNA MSAs were classified protein), so the residue type is explicit or the EMBOSS first-row rule.

## Audit A2-15 (B09 F90): distmat corrections / options and the full EMBOSS matrix set
- **Opened:** EMBOSS 6.6 `emboss/distmat.c` (all of it: `distmat_Tamura`, `_Kimura`, `_KimuraProt`, `_calc_match`, `_uncorrected`, `_TajimaNei`, `_JinNei`, `_JukesCantor`, `_checkRY`, `_checkSubs`, `_checkambigProt`, `_getSeq`, `_outputDist`) and `ajax/core/ajbase.c` (`ajBaseAlphacharCompare`, `baseInit`) from raw.githubusercontent.com/kimrutherford/EMBOSS (neutral User-Agent), `/usr/share/EMBOSS/data` (`Matrices.protein`, `Matrices.nucleotide`, `Ebases.iub`, the matrix files), Biopython 1.88 `substitution_matrices` (NCBI aliases).
- **Numbers:** 1 100 distmat 6.6.0 runs (every method, position, `-ambiguous`, `-calculatea`, gap weights 0 / 0.5 / 1 / 2.5, 50 runs the program rejects or that do not terminate): the first comparison equal on 1 053 of 1 110 test cases; the 57 differences were (a) 1 case of the EMBOSS reader: every ROW is typed on its own and a nucleotide-looking row has X turned into N (`ajSeqSetNuc`; a protein alignment whose row `KCX...V...` is read as nucleotide stops matching the X of its neighbours) - added; (b) 56 runs of Tajima-Nei / protein Kimura with a position other than 123 whose C code loops over `mlen` instead of the selected length, i.e. reads heap memory past the buffers (about half of them equal the NUL-padded model, the rest do not) - refused (`NotSupportedException`) and asserted. After that 1 100 / 1 100 equal. 70 EMBOSS matrices cell for cell (the C# table is generated from the files; the test reads an independent TSV parsed by the Python generator); BLOSUM45 / 50 / 62 / 80 / 90, PAM30 / 70 / 250, NUC.4.4 = Biopython on every common symbol (0 differing cells).
- **Findings:** `distmat -ambiguous` on nucleotides is broken in 6.6.0 (`ajBaseAlphacharCompare` uses `len = |list| - 1`: an unambiguous base examines nothing, `1/0 = inf`, and a non-IUB first character gives `(unsigned)-1` iterations = no termination); `-position` for proteins accepts only 123 and 12 and for 12 the buffer copy reads one character past the end when the length is not a multiple of 3 (NUL, which C's `strchr` treats as a member of every set); the printed value is `100 D` with two decimals of a `float`, `-0.00` appears for a negative zero.

## Audit A3-13c (B09 F148): showalign

Sources actually opened: EMBOSS 6.6.0 `emboss/showalign.c` in full (`showalign_Getrefseq`, `_NiceMargin`, `_Convert` with `_MakeAll` / `_MakeIdentity` / `_MakeNonidentity` / `_MakeSimilar` / `_MakeDissimilar`, `_Order`, `_Output`, `_OutputNums`, `_OutputTicks`, `_OutputSeq`), `nucleus/embcons.c` (`embConsCalc`, with the gap option), `nucleus/embshow.c` (`embShowUpperRange`), `showalign.acd`, installed `showalign` 6.6.0 (raw.githubusercontent.com/kimrutherford/EMBOSS, User-Agent `Seqeron-review/1.0`).

`MultipleAlignmentDisplay.ShowAlign` ports the program except `-html` / `-highlight` (not part of this item; registered). Measured on the binary: the alignment is NOT upper-cased (the consensus letter is the letter of the winning row, lower case kept; the Make* functions compare case-insensitively); the consensus is `embConsCalc` with `-gaps Y` (default: a gap column can be the consensus; `identical[0]` and `matching[0]` count gaps), shared with `cons` (`EmbossConsensusCore`, extended by `gaps`, `keepCase`, a matrix override and the `ajSeqSetNuc` x -> n switch, which showalign does not apply); `ajStrAssignSubS` clips BOTH ends of a window into the string (a number line shorter than the margin yields its last character: the line is one blank longer than the margin); the code of a character that is not a matrix label scores 0 against everything; `qsort` of equal names / equal similarities is taken as stable.

Oracle `showalign_oracle.jsonl` (`generate_showalign_oracle.py`): 600 random alignments (1-8 rows, 4-120 columns, nucleotide with ambiguity letters / protein with B Z X, upper and lower case, gaps `- . ~`, ragged rows, unique names) with random `-refseq` (0, number, name), `-[no]bottom`, `-show` A I N S D, `-order` I A S, `-[no]similarcase`, `-[no]consensus`, `-[no]number`, `-[no]ruler`, `-width` 1-80, `-margin` -1 0 3 12, `-uppercase`, `-matrix`, `-plurality`, `-setcase`, `-identity`, `-[no]gaps`: **600 / 600 texts equal**. Tests: `MultipleAlignmentDisplay_ShowAlign_Tests` (3), MCP `showalign`.

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

Sources actually opened: EMBOSS `ajax/core/ajseqread.c` (`seqAppend`: `ajStrKeepSetAlphaC(.., "*.~?#+-")`), `ajseqtype.c` (`ajSeqTypeCheckIn`, `seqTypeFix`: an ungapped type calls `ajStrRemoveGap`; `seqType[]` table: every ACD type of the programs is ungapped, `gapany` for consambig / showalign), `ajstr.c` (`ajStrRemoveGap` keeps `A-Za-z*`), `ajseq.c` (`ajSeqIsNuc`, `ajSeqGapStandard`) from raw.githubusercontent.com/kimrutherford/EMBOSS (User-Agent `Seqeron-review/1.0`), the installed ACD files (`/usr/share/EMBOSS/acd`: `type:` of every sequence input) and the installed programs. Measured with `scripts/alignment/measure_emboss_reader_chars.py` (each of the characters 1..255 except newline / carriage return / `>` inserted into the middle of a sequence, both sequences, nucleotide and protein): for showalign the output equals the output without the character for exactly 199 characters (all except the 52 ASCII letters and `*`); `*` is kept (a protein stop; an error in a nucleotide-typed input).

Oracle `emboss_gap_oracle.tsv` (`generate_emboss_gap_oracle.py`, seed 20261126): 18 programs x 240 rows, showalign (240 alignments): **equal on every row**. Audit example: `needle AC-GTTAGC.TT* ACGTTAGCTT` prints `# Length: 11`, `# Score: 58.0` (the port scored 43 on the uncleaned text).

## Audit A4-1 (B09 F157): cancellation

Mechanism: `AlignmentCancellation` (SequenceAligner.Cancellation.cs, A2-40); tests `SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `MultipleAlignmentDisplay.ShowAlign` poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

## Audit A4-9 (B09 F164): showalign -sbegin / -send

Sources actually opened: EMBOSS 6.6.0 `emboss/showalign.c` (`begin = ajSeqsetGetBegin(seqset) - 1; end = ajSeqsetGetEnd(seqset) - 1; for(pos = begin; pos <= end; pos += width)`), `ajax/core/ajseq.c` (`ajSeqsetGetBegin` / `ajSeqsetGetEnd`: a positive value is 1-based, a negative one counts from the end, `ajCvtSposToPosStart` clamps into [imin, len - 1]), installed `showalign` 6.6.0.

`ShowAlignOptions.Begin` / `End` (MCP `showalign` `sbegin` / `send`) cut the DISPLAYED columns; the consensus is that of the whole alignment and the numbers / ticks keep the positions of the whole alignment (measured: `-sbegin 3 -send 12` shows the columns 3 - 12 numbered from the full alignment; `-sbegin 30` = the last column; `-sbegin 8 -send 4` = column 8). `scripts/analysis/generate_showalign_range_oracle.py` -> `showalign_range_oracle.jsonl`, **400 rows** (0, 1-based, negative, beyond the alignment, end before begin): equal on 400 / 400.
