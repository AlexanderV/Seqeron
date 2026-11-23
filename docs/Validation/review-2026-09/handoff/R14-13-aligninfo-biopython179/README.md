**STATUS: INTEGRATED (2026-11-23, B09.md F308, `returned/e6f5e012-R14-13-aligninfo-biopython179-out.zip`: cases.jsonl 64 cases, results.jsonl 64 results, Biopython 1.79 on Python 3.9.25).** The 6 199 calls are embedded as `r14_13_aligninfo_{cases,results}.jsonl.gz` and `AlignInfo_Returned_Oracle_Tests` requires every value and every exception (kind and message) of `SummaryInfo` / `Pssm` to equal the returned one: 6 199 / 6 199, no divergence of the port (the 128 calls with `chars_to_ignore` given as a string are Python's `isinstance` check, counted). No further run.

# Hand-off R14-13-aligninfo-biopython179 (temporary Docker package for the user's Mac)

**Why.** The audit item R14-13 asks for the members of `Bio.Align.AlignInfo.SummaryInfo`. Biopython 1.88 (the version the rest of the work is checked against) has no such members any more (the class is deprecated and keeps only `get_column`); the last release that has
them is **1.79**. The sandbox that builds the port is not allowed to run downloaded code, so the reference class cannot be run there: this package runs it on your Mac.

**What it does.** Builds a `linux/amd64` image (Python 3.9, `pip install biopython==1.79`), writes `cases.jsonl` (64 small alignments: the three of the Biopython tests, one with sequence weights, 60 random nucleotide / protein alignments with gaps, point differences and weights; per alignment about 100 calls:
`dumb_consensus` / `gap_consensus` with 5 thresholds x 2 ambiguity letters and `require_multiple`, `replacement_dictionary` with 5 letter sets and the deprecated `skip_chars`, `pos_specific_score_matrix` with 6 `chars_to_ignore` values and axis sequences, `information_content` with expected-frequency tables (complete, with an extra
letter, with a letter missing, empty), 3 logarithm bases, 3 pseudo counts, ranges, bad arguments, `ic_vector`, `print_info_content`, `get_column`) and runs every call on a fresh `SummaryInfo`. Every result and every exception (type and message) is recorded in `out/results.jsonl`.

**Commands (the 3 you need).**
```
cd docs/Validation/review-2026-09/handoff/R14-13-aligninfo-biopython179
./run.sh                      # builds the image, runs the calls, writes ./out and ./returned/R14-13-aligninfo-biopython179-out.zip
ls -l returned                # a minute or two (emulated on Apple Silicon)
```
**Send back.** The archive `returned/R14-13-aligninfo-biopython179-out.zip` (or the `out` directory) placed in `docs/Validation/review-2026-09/handoff/R14-13-aligninfo-biopython179/returned/`.

**How it is integrated.** The session that integrated it unzipped it into an empty scratch directory (data only), copies `cases.jsonl` + `results.jsonl` to `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (gzip), and added the test `AlignInfo_Returned_Oracle_Tests` that builds each alignment as a `MultipleSeqAlignment`, runs the same calls on `SummaryInfo`
and requires every value and every error message to be equal (floats to the bit: both sides compute `count / total * log(...)` with the same `libm`). A call on which the port differs is a defect of the port (fixed), unless Biopython itself misbehaves (then the row is refused with the evidence, as for the other oracles).
This README has `STATUS: INTEGRATED`.
