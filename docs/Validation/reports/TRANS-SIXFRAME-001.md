# Validation Report: TRANS-SIXFRAME-001 — Six-Frame Translation and ORF finding

- **Validated:** 2026-06-15   **Area:** Translation
- **Canonical method(s):** `Translator.TranslateSixFrames(DnaSequence, GeneticCode?)`,
  `Translator.FindOrfs(DnaSequence, GeneticCode?, int minLength, bool searchBothStrands)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)
- **NCBI The Genetic Codes, table 1** — reconstructed the standard code locally from the verbatim
  `AAs/Starts/Base1/Base2/Base3` lines recorded in the Evidence doc. Independent decode gives
  start codons **{ATG, CTG, TTG}** and stop codons **{TAA, TAG, TGA}** — matches the description and
  `GeneticCode.CreateStandard()`. (Decode script `/tmp/sixframe_check.py`.)
- **Biopython 1.85** (installed locally; `Bio.SeqUtils.six_frame_translations`, `Seq.translate`,
  `Seq.reverse_complement`) — reference implementation. Forward frames = translate(seq) at offsets
  0/1/2; reverse frames = translate(reverse_complement(seq)) at offsets 0/1/2 ("independent offset"
  convention). The `six_frame_translations` text output labels these exactly as the repo does
  (−1 = rc offset 0, etc.).
- **EMBOSS transeq / getorf docs** (as cited in Evidence) — six frames total (`-frame 6`); two
  reverse-numbering conventions exist (phase-locked default vs. the "alternative" independent-offset
  used here); getorf `-find 1` = START→STOP region, inclusive of the STOP, both strands by default.
- **Wikipedia — Reading frame** — six reading frames (3 forward + 3 reverse-complement read 5'→3').

### Formula check
- Six frames keyed {+1,+2,+3,−1,−2,−3}; forward = offset f−1 of input, reverse = offset f−1 of
  reverse complement; only complete codons consumed (trailing 1–2 nt dropped). All match Biopython.
- ORF = START→STOP; EndPosition inclusive of stop's last base; protein includes start residue,
  excludes the stop. Matches getorf `-find 1`.

### Edge-case semantics check
- Empty → 6 empty frames; null → `ArgumentNullException`; partial trailing codon dropped; no START →
  no ORF; ORF below minLength filtered; ORF running off the end with no STOP → open ORF emitted.
  All have defined, sourced behaviour.

### Independent cross-check (numbers)
Independent reimplementation built from NCBI table-1 raw lines, and Biopython 1.85, both reproduce the
39-nt dataset (`ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`) exactly:

| Frame | Protein (independent + Biopython) |
|-------|-----------------------------------|
| +1 | `MAIVMGR*KGAR*` |
| +2 | `WPL*WAAERVPD` |
| +3 | `GHCNGPLKGCPI` |
| −1 | `LSGTLSAAHYNGH` |
| −2 | `YRAPFQRPITMA` |
| −3 | `IGHPFSGPLQWP` |

ORF datasets reproduced independently: `GGGATGAAACCCTAAGGG` → (Start 3, End 14, +1, `MKP`);
`GGTTGAAAGGGTAACC` → (Start 2, End 13, +3, `LKG`); `CCCTTAGGGTTTCATCCC` revcomp = `GGGATGAAACCCTAAGGG`
→ reverse ORF (−1, `MKP`); `ATGAAATAA` minLength 3 → empty; `ATGAAACCCGGG` → open ORF (0, 11, +1, `MKPG`).

### Findings / divergences (Stage A)
1. **Reverse-frame numbering convention** — repo follows Biopython "independent offset" (−1 = rc
   offset 0), the EMBOSS-documented *alternative* (not the phase-locked default). Documented as an
   accepted convention; confirmed identical to Biopython 1.85. Note, not a defect.
2. **Initiator residue at alternative start codons** — NCBI states the initiator codon is "by default
   translated as methionine" even when TTG/CTG. The repo translates the start codon by its standard
   residue (TTG→L), not forced to M (EMBOSS getorf offers `-methionine` to force M; default keeps
   the codon's own residue in some pipelines). Documented assumption; both behaviours exist in tools.
   Note, not a defect.
3. **minLength unit** — getorf `-minsize` is nucleotides; repo `minLength` is amino acids. Documented
   API-shape deviation. Note.

→ Description is biologically/mathematically correct; three documented convention/API notes.

## Stage B — Implementation

### Code path reviewed
- `Translator.cs:115-136` `TranslateSixFrames`; `:138-162` `TranslateSequence`; `:78-97` `FindOrfs`;
  `:164-229` `FindOrfsInSequence`; `:235-250` `OrfResult`. `GeneticCode.cs` standard table + start/stop.

### Formula realised correctly?
Yes. Forward loop fills +1/+2/+3 from offsets 0/1/2 of the input; reverse fills −1/−2/−3 from offsets
0/1/2 of one shared reverse complement. Loop bound `i + 3 <= length` implements the Biopython
trailing-codon truncation. ORF scanner enters on START, accumulates the actual start residue, emits on
in-frame STOP with inclusive end `i+2`, and emits an open ORF at end-of-strand. All verified against
the independent reimplementation (byte-identical on every dataset).

### Cross-verification table recomputed vs code
Every test value (M1–M13, S1–S2, C1–C2, plus the new open-ORF case) recomputed from external sources
and re-derived through the independent reimplementation — all match the values the tests assert.

### Variant/delegate consistency
M4/M5 confirm `TranslateSixFrames` forward/reverse frames equal `Translate(dna/revComp, frame:0/1/2)`.
Consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every expected protein/position traces to NCBI table 1 + Biopython
  1.85 (retrieved this session), not to repo output. Exact-value assertions (`Is.EqualTo`), no
  `Greater`/`AtLeast`/ranges where exact values are known.
- **No green-washing:** no weakened assertions, widened tolerances, or skipped tests introduced.
- **Coverage gap found & fixed:** the documented Stage-A edge case "ORF runs to sequence end with no
  STOP → open ORF emitted" (INV-05 / doc §6.1, getorf incomplete-ORF) had **no test**. Added
  `FindOrfs_OrfRunsToSequenceEndWithoutStop_EmitsOpenOrf` (hand-verified: `ATGAAACCCGGG` → Start 0,
  End 11, +1, `MKPG`), exercising the `Translator.cs:218-227` branch. All other public methods/overloads
  and Stage-A branches are covered.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6529` (one unrelated pre-existing
  benchmark skipped: `MFE_Benchmark_AllScenarios`); `dotnet build` 0 errors. Changed test file builds
  warning-free.

### Findings / defects
No code defect. One test-coverage gap (open-ended ORF) fixed in this session by adding a sourced test.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (3 documented convention/API notes — reverse-frame numbering,
  initiator residue, minLength unit — none affecting correctness).
- **Stage B:** PASS-WITH-NOTES (code faithful to validated description; added one missing edge-case
  test for the open-ended-ORF branch).
- **End-state:** ✅ CLEAN — no defect; coverage gap completely fixed; build + full suite green.
- **Test-quality gate:** PASS.
