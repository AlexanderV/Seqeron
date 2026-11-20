# A2-29b: the real ksw2 / minimap2 functions as the oracle (Docker package)

**Why.** `Ksw2Aligner` (`src/.../Seqeron.Genomics.Alignment/Ksw2Aligner*.cs`) ports `ksw_extz`, `ksw_extz2_sse`, `ksw_extd`, `ksw_extd2_sse` (lh3/ksw2, minimap2 v2.26). The exact cases (no Z-drop, no approximate flags) are
verified against an independent dynamic program (`ksw2_exact_oracle.tsv`, 1 320 cases). Z-drop, the end bonus on cut extensions, the approximate-maximum / approximate-drop flags and the byte-level behaviour of the SIMD code
can only be fixed by the real C program, which cannot be compiled and run in the review sandbox (compiling downloaded source is refused). This package does it on your machine.

**Run (needs Docker; x86_64 or Apple Silicon with emulation).**

```bash
cd docs/Validation/review-2026-09/handoff/A2-29b-ksw2
./run.sh                      # builds the image (pinned minimap2 v2.26 + lh3/ksw2 sources, sha256 checked), runs cases.txt, writes ./out/
zip -r returned.zip out && mkdir -p returned && cp -r out returned/ && rm returned.zip
```

`cases.txt` (2 000 cases: 4 variants x 500, 170 per variant with Z-drop / approximate flags) is produced by `python3 scripts/alignment/generate_ksw2_oracle.py cases cases.txt`.

**Send back.** The directory `docs/Validation/review-2026-09/handoff/A2-29b-ksw2/returned/out/` (at least `ksw2_oracle.txt`, `versions.txt`, `checksums.txt`; the `uname -m` and the CPU flags are recorded in `versions.txt`).

**Integration (next session).** `returned/out/ksw2_oracle.txt` (`id variant max zdropped max_q max_t mqe mqe_t mte mte_q score n_cigar reach_end cigar`) is converted to `tests/.../TestData/Alignment/ksw2_oracle.tsv` and a test
`Ksw2Aligner_RealProgram_Tests` compares every field and the CIGAR of all 2 000 cases (fields of `Ksw2Result`). Differences decide where the port differs from the SSE4.1 code (the SSE2-only build `ksw2_sse2.txt` is a cross-check of the
two code paths of the C file). The item A2-29b stays DOABLE until this is integrated.
