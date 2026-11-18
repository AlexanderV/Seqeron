# TEMPORARY hand-off package A1-10c-muscle - the MUSCLE 3.8.1551 oracle (Docker, for the user's Mac)

**Why.** B09 item A1-10c part (a): `SequenceAligner.MultipleAlignMuscle` is a step-by-step port of MUSCLE **3.8.1551** (Edgar 2004, NAR 32:1792; the final v3.8 release; source =
Debian package `muscle3` 3.8.1551-2, `sources.debian.org`) for `-maxiters 1` (progressive draft) and `-maxiters 2` (plus the tree-dependent refinement of iteration 2).
The installed `muscle` of the audit machine is **5.1** - a different program (PPP / probcons-like pipeline) - so it cannot be the oracle, and the sandbox refuses to run a
binary compiled from downloaded source. This package lets **you** build and run the original v3.8 once. The item stays **DOABLE (not ticked) for its verification** until the
returned output is integrated (see below); the port itself, its tests of the published definitions and an independent re-implementation are already committed.

## What it does

1. `Dockerfile` (base `debian:12.7-slim`, linux/amd64 manifest pinned by digest) fetches the 183 files of muscle3 3.8.1551-2 from `sources.debian.org` (`fetch_src.py`, neutral
   User-Agent `Seqeron-review/1.0`, sha256 of every file recorded), compiles them with the upstream Makefile flags (`-O3 -funroll-loops -DNDEBUG=1`, no patch) and records the
   binary's sha256. If the plain build fails with the image's g++, the same files are rebuilt with `-fpermissive -std=gnu++98` (no change of floating-point semantics).
2. `inputs/sets.tsv` (committed, 420 sets; `make_inputs.py` re-creates it deterministically): 270 protein, 110 DNA, 40 RNA families with the option variants of the port:
   `-maxiters 1|2`, `-le/-sp/-sv`, `-gapopen/-gapextend/-center`, `-cluster1/-cluster2`, `-sueff`, `-distance1`, `-maxtrees 0`; leading M (methionine hack), wildcards, lower case,
   duplicates, sequences shorter than 5 residues, 2-sequence sets.
3. `run_sets.py` runs every set with `muscle3 -in x.fa -out o.fa -quiet -seqtype ... <options>` and writes `out/muscle3_oracle.tsv` (rows in input order + the order MUSCLE wrote them).

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
cd docs/Validation/review-2026-09/handoff/A1-10c-muscle
./run.sh                  # builds the image (a few minutes: compiles MUSCLE), runs 420 muscle3 calls, zips ./out
ls -l returned/           # A1-10c-muscle-out.zip
```

## What to send back

The archive `returned/A1-10c-muscle-out.zip` (or `.tar.gz`) - place it in this repository under
`docs/Validation/review-2026-09/handoff/A1-10c-muscle/returned/` (it contains `out/ENV.txt` with the date, `uname -m`, the binary checksum and `muscle3 -version`,
`out/fetch_src.sha256`, `out/muscle3_oracle.tsv`, `out/run.log`). A run that prints "failed: N" with N > 0 is still useful (the failing sets are listed with MUSCLE's message).

## How the next session integrates it

`out/muscle3_oracle.tsv` becomes `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/muscle3_oracle.tsv` (embedded resource) and a test
`SequenceAligner_MuscleOracle_Tests` asserts, for every set, that `MultipleAlignMuscle(seqs, options)` returns the oracle rows (input order) and `TreeOrder` returns the oracle order.
Every set that differs is a defect of the port (or of this README's reading of the source); nothing is "accepted" - the port is corrected against the oracle, the independent
re-implementation (`scripts/alignment/muscle3_reference.py`) is compared again, and F65 of `B09.md` is amended. After that A1-10c (a) is ticked and the MCP / catalog text
`multiple_align_muscle` loses its "awaiting oracle" qualifier. If `ENV.txt` shows a `muscle3 -version` other than 3.8.1551 nothing is integrated until that is explained.
