# INTEGRATED (2026-11-19; not needed from the user: run locally with the Debian `muscle3` package, audit A2-1 / B09 F89) - hand-off package A1-10c-muscle-r3 - ROUND 3 of the MUSCLE 3.8.1551 oracle (Docker, for the user's Mac)

> **STATUS (2026-11-19): INTEGRATED / NOT NEEDED.** `apt-get install muscle3` (Debian 3.8.1551-3) works in the sandbox; `python3 scripts/alignment/generate_muscle3_r3_oracle.py` ran the 60 sets with this package's harness `run_sets.py` and wrote `tests/.../TestData/Alignment/muscle3_r3_oracle.tsv` (`returned/ENV-local.txt`); `SequenceAligner_MuscleOracleRound3_Tests` asserts 60 / 60. The Docker run below is no longer required; if a returned `returned/A1-10c-muscle-r3-out.zip` is ever added, the generator prefers it. This package stays as provenance (inputs + harness).

**Why.** Round 2 (`../A1-10c-muscle-r2`, 638 sets) is returned and integrated: `SequenceAligner.Muscle*` equals the unmodified MUSCLE v3.8.1551 on **630 / 630** recorded sets (after the fixes it found: the k-mer
bit count, `-brenner` in iteration 2, the leaf-index / sequence-id weight lookup of `-usetree`, negative zero in the files, `-maketree` ignoring `-distance2`, the refusals). Eight sets could not be recorded
(`-matrix blosum62nostar`: the program's own warning prints the byte 0xff, which the round-2 harness could not decode), and a few corners were never exercised. This package runs 60 sets:

| group | sets | what is decided |
|---|---|---|
| `-matrix blosum62nostar` (R330 ... R350) | 8 | the 8 unrecorded sets of round 2, same ids and data; the port reads such a file the way `readmx.cpp` does (last column not read, "Matrix is not symmetrical" warning) |
| `-matrix` with DNA / RNA (R638 ... R645) | 8 | the port refuses `-matrix` for nucleotides; the source reads the matrix with the nucleotide alphabet and then switches `g_Alpha` to amino acid and the score to SP - the result (or crash) decides what the port implements |
| one- and two-sequence profiles / trees / scores / refinements (R646 ... R668) | 23 | `-root1/-root2 midlongestspan|minavgleafdist` on a one-sequence profile crashed in round 2 (segmentation fault, modelled as `NotSupportedException`); all other one- / two-sequence cases are recorded |
| `align` with one, two and three sequences (R669 ... R681) | 13 | tree files, `-usetree`, roots, clusterings, `-computeweights` on tiny inputs |
| score flags `-le -sp -sv -spn` on nucleotides and proteins (R682 ... R689) | 8 | the program replaces a protein score by SPN for nucleotides and SPN by LE for proteins |

The harness is the round-2 harness with one change: console output is decoded with `errors='backslashreplace'` (the program prints non-UTF-8 bytes), so nothing is lost.

## Run it (Apple Silicon: the image is x86-64 and is emulated; docker must be running)

```bash
cd docs/Validation/review-2026-09/handoff/A1-10c-muscle-r3
./run.sh                  # builds the image (a few minutes: compiles MUSCLE), runs 60 muscle3 calls, zips ./out
ls -l returned/           # A1-10c-muscle-r3-out.zip
```

`run.sh` needs only `docker` (and `zip` or `tar`). The image is built from the unmodified Debian source `muscle3 3.8.1551-2` with the upstream flags; no account data, the download uses the User-Agent
`Seqeron-review/1.0`. Nothing aborts: a set on which MUSCLE fails (or crashes) is recorded with its exit status and message - those rows are part of the oracle.

## What to send back

`returned/A1-10c-muscle-r3-out.zip` (or `.tar.gz`), kept in this repository under `docs/Validation/review-2026-09/handoff/A1-10c-muscle-r3/returned/` (`out/ENV.txt`, `out/build.log`, `out/preflight.txt`,
`out/fetch_src.sha256`, `out/muscle3_r3_oracle.tsv`, `out/run.log`).

## How the next session integrates it

`python3 scripts/alignment/generate_muscle3_r2_oracle.py` is extended to merge the round-3 rows (the 8 sets of round 2 keep their ids) into `TestData/Alignment/muscle3_r2_oracle.tsv`; the test
`SequenceAligner_MuscleOracleRound2_Tests` then asserts all 690 recorded sets (630 + 60). Differences are fixed in the port stage by stage; `-matrix` for nucleotides is implemented from the source
when the oracle shows the program runs. If `ENV.txt` shows a `muscle3 -version` other than 3.8.1551, nothing is integrated until that is explained.
