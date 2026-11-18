# TEMPORARY hand-off package A1-10d - ClustalW 2.1 end-gap options (Docker, for the user's Mac)

**Why.** B09 item A1-10d (end-gap handling of ClustalW's profile stage): `UserParameters::useEndGaps` ("use end gaps" - gaps at the ends of the
sequences are counted as gaps in `ProfileBase::calcGapCoeff`) and `UserParameters::endGapPenalties` ("end gap penalties" - the terminal gaps of
a profile alignment are penalised: `profile[0][GAPCOL] = gapCoef`, `openPenalty*/extPenalty*/gapPenalty*` no longer return 0 at the profile
ends) are ported in `SequenceAligner.ClustalWMultiple.cs` (internal `CwConfig.UseEndGaps / EndGapPenalties`, default false) but **cannot be
reached from the `clustalw` command line**: `-ENDGAPS` only executes `setUseEndGaps(false)` (verified no-op on the 54 protein oracle sets that carry it) and nothing
sets `endGapPenalties`. They exist in the program's interactive menu only. The sandbox the audit ran in refuses to execute a binary that was
compiled from the downloaded source, so this package lets **you** run it once. The item stays **DOABLE (not ticked)** until the returned output is
integrated (see below).

## What it does

1. `Dockerfile` (base `debian:12.7-slim`, linux/amd64 manifest pinned by digest) fetches `src/` of Debian's clustalw 2.1+lgpl-7 from
   `sources.debian.org` (`fetch_src.py`, neutral User-Agent `Seqeron-review/1.0`), applies **one** edit (the `-ENDGAPS` handler of
   `interface/CommandLineParser.cpp` reads the environment variables `CW_USEENDGAPS` / `CW_ENDGAPPEN` instead of the constant `false`), compiles with
   g++ -O2 and records the binary's sha256.
2. `inputs/sets.tsv` (committed; `make_inputs.py` re-creates it from the committed oracle TSVs): 220 protein sets (no `-NOVGAP`) of
   `clustalw_protein_msa_oracle.tsv` and 120 DNA sets of `clustalw_msa_full_oracle.tsv`, each with its options and the rows of the *installed* clustalw 2.1.
3. `run_sets.py` runs every set four times: control (no variable), `CW_USEENDGAPS=1`, `CW_ENDGAPPEN=1`, both, always with `-ENDGAPS`
   and the set's own options, and writes `out/end_gaps_oracle.tsv` (id, type, mode, rows) and `out/control_check.txt` (how many control runs equal the
   committed rows of the installed binary - this validates the self-built binary; it should be 340 of 340).

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
cd docs/Validation/review-2026-09/handoff/A1-10d
./run.sh                  # builds the image (a few minutes: compiles ClustalW), runs ~1 360 clustalw calls, zips ./out
ls -l returned/           # A1-10d-out.zip
```

## What to send back

The archive `returned/A1-10d-out.zip` (or `.tar.gz`) - place it in this repository under
`docs/Validation/review-2026-09/handoff/A1-10d/returned/` (it contains `out/ENV.txt` with the date, `uname -m` and the binary checksum,
`out/control_check.txt`, `out/end_gaps_oracle.tsv`, `out/CommandLineParser.patched.cpp`).

## How the next session integrates it

`out/end_gaps_oracle.tsv` becomes `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/clustalw_end_gaps_oracle.tsv`; the options
`UseEndGaps` / `EndGapPenalties` are added to `ClustalWOptions` and `ClustalWProteinOptions` (the ported code paths exist), the MCP tool gets
`useEndGaps` / `endGapPenalties`, and a test asserts every row of the three non-control modes (and that each mode changes the rows of many sets).
If `control_check.txt` is not 340 of 340 the self-built binary differs from the installed one and nothing is integrated until that is explained.
