**STATUS: RETURNED (2026-11-23, run 2: returned/out/, returned/A2-29b-ksw2-out.run2.zip); NO FURTHER RUN NEEDED. INTEGRATION = B09 WP68.** Findings of the run: (1) SSE4.1 and SSE2 harnesses: 1000 rows each, byte-identical files (ksw2_sse41.txt = ksw2_sse2.txt in size), ASan/UBSan clean on all 2000 cases. (2) The scalar `ksw_extz` / `ksw_extd` of lh3/ksw2 has a heap-buffer-overflow: ksw2_extz.c:116 `eh[j].h = h1, eh[j].e = KSW_NEG_INF;` writes at j = st > qlen when a row lies completely right of the query (tlen > qlen + w + 1); proven by ASan (asan_scalar.txt, first case 26: qlen 30-ish, w 8) and valgrind (valgrind_scalar.txt: invalid write 0 bytes after a block of size 72); 33 cases abort (`ksw2_scalar.txt.crashes`: id and exit code 134/139), the other 967 scalar rows are valid unless the same overflow touched them silently (the integration decides per case with the port's own band bounds). The crashed / overflow cases are REFUSED as undefined behaviour of the reference program (not bent). (3) mm2/ksw2.h differs from the recorded checksum (the .c files match); recorded in the integration.

**STATUS: RERUN-NEEDED (2026-11-23).** Run 1: the SSE4.1 harness wrote its 1000 result lines, the scalar harness (`ksw_extz` / `ksw_extd` of lh3/ksw2) aborted with `malloc(): invalid size (unsorted)` (returned/run_scalar.err.run1) and, because stdout was buffered, lost every scalar line. The harness now prints each line at once, names the running case (`CASE id` on stderr), `inside.sh` resumes after a crash and writes the crashed ids to `<out>.crashes`, and an ASan/UBSan build of both harnesses plus valgrind produce the reports `asan_scalar.txt`, `asan_sse41.txt`, `valgrind_scalar.txt`: they prove whether the abort is in the reference program or in the harness. Run `./run.sh` again and return the whole `out/` folder (zip). **Earlier status:** AWAITING-USER (since 2026-11-19). Integration target: `TestData/Alignment/ksw2_oracle.tsv` (audit item A2-29b: Z-drop, end bonus, approximate flags of the SIMD code); nothing has been returned yet (no `returned/` directory).

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
