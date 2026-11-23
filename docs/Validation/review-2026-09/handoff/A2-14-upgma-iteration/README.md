**STATUS: INTEGRATED (2026-11-23, B09.md F297)** - `returned/A2-14-upgma-iteration-out.zip` (800 sets with rows, 0 without) is `TestData/Alignment/clustalw_upgma_iteration_oracle.tsv`, locked by `SequenceAligner_ClustalWUpgmaIteration_Tests`; the port equals the program on 800 / 800 sets.

# TEMPORARY hand-off package A2-14-upgma-iteration - oracle for ClustalW 2.1 `-CLUSTERING=UPGMA` with `-ITERATION` and the profile commands (Docker, for the user's Mac)

**Why.** `handoff/A2-14-clustalw-upgma` (returned 2026-11-19) verified `-CLUSTERING=UPGMA -ALIGN` on 600 / 600 sets with the unmodified Debian source built at `-O0` (the installed binary and an `-O2` build of the same source
segmentation-fault on every UPGMA run: undefined behaviour of the program). Three more code paths use UPGMA and are ported from the source (`SequenceAligner.ClustalWIteration.cs`: `GetWeightsForProfileAlignUpgma`,
`UpgmaProfileWeights`) but have no binary to compare with: `TreeInterface::getWeightsForProfileAlignUPGMA` is reached by `-CLUSTERING=UPGMA` together with `-ITERATION=alignment|tree`, `-PROFILE` and `-SEQUENCES`.
This package builds the same unmodified source at `-O0` (only the O0 build; no patch) and runs 800 sets through exactly those commands.

| | |
|---|---|
| `inputs/sets.tsv` | 800 sets (400 DNA, 400 protein): 250 `align-alignment` (`-ITERATION=alignment [-NUMITER]`), 200 `align-tree` (`-ITERATION=tree`), 175 `profile`, 175 `sequences`; options `-NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -NUMITER -QUICKTREE`; written by `make_inputs.py` (deterministic: the committed file is reproduced exactly) |
| `run_sets.py` | for `profile` / `sequences` first aligns the two parts with the program itself (NJ, no UPGMA) to make the profiles, then runs the command with `-CLUSTERING=UPGMA`; records profiles, rows, exit status, the binary's output |

## Run it (Apple Silicon: x86-64 image, emulated; docker must be running)

```bash
cd docs/Validation/review-2026-09/handoff/A2-14-upgma-iteration
./run.sh                  # builds the image (a few minutes: one build of ClustalW at -O0), runs 800 sets, zips ./out
ls -l returned/           # A2-14-upgma-iteration-out.zip
```

`run.sh` needs only `docker` (and `zip` or `tar`). Nothing aborts: a set on which the binary fails is recorded with its exit status. No account data; the source download uses the User-Agent `Seqeron-review/1.0`.

## What to send back

`returned/A2-14-upgma-iteration-out.zip` (or `.tar.gz`), kept in this repository under `docs/Validation/review-2026-09/handoff/A2-14-upgma-iteration/returned/` (`out/ENV.txt`, `out/build_O0.log`, `out/preflight.txt`,
`out/fetch_src.sha256`, `out/upgma_iteration_oracle.tsv`, `out/run.log`).

## How the next session integrates it

`out/upgma_iteration_oracle.tsv` becomes `tests/.../TestData/Alignment/clustalw_upgma_iteration_oracle.tsv` (converted by a script in `scripts/alignment/`, as `generate_clustalw_upgma_oracle.py` did) and a test
`SequenceAligner_ClustalWUpgmaIteration_Tests`: `MultipleAlignClustalW` with `Clustering = Upgma` and `Iteration = Alignment | Tree`, `AlignClustalWProfiles` / `AlignClustalWSequencesToProfile` with `Clustering = Upgma`
must equal the rows on every set the program finished. Differences are fixed stage by stage (weights of `getWeightsForProfileAlignUPGMA` before `calcPairwiseForProfileAlign`, `distanceTree` staying false, profile 1 weight 0
when it has fewer than two sequences). The `-ITERATION=tree` sets also exercise the UPGMA tree weights inside `iterationOnTreeNode`.
