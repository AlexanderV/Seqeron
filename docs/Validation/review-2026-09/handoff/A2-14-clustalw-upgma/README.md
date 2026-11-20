**STATUS: INTEGRATED (2026-11-19, B09 F101, `returned/A2-14-clustalw-upgma-out.zip` -> `clustalw_upgma_oracle.tsv`).**

# TEMPORARY hand-off package A2-14-clustalw-upgma - oracle for ClustalW 2.1 `-CLUSTERING=UPGMA -ALIGN` (Docker, for the user's Mac)

**Why.** B09 audit item A2-14 asks for ClustalW's `-CLUSTERING=UPGMA`. The port is written (`SequenceAligner.ClustalWUpgma.cs`: `UPGMAAlgorithm` / `Node` / `RootedGuideTree` of `src/tree/UPGMA`, wired into
`MultipleAlignClustalW` as `ClustalWOptions.Clustering = ClustalWClustering.Upgma`, DNA and PROTEIN), but its OUTPUT cannot be compared with a program here: the installed Debian
`/usr/bin/clustalw` (2.1+lgpl-7) **segmentation-faults on every `-CLUSTERING=UPGMA` run** (`-ALIGN` and `-TREE`, 2, 3 and 4 sequences alike; only the two-sequence tree file is written before the crash). The same upstream
source built at `-O0` ran the UPGMA alignment to the end (AddressSanitizer build in the audit sandbox: "Group 1: Sequences: 2 Score:228 / Group 2: Sequences: 3 Score:199 / Alignment Score 216" on three 13-mers).
The sandbox's safety classifier refused to run further self-compiled copies, so this package builds ClustalW from the unmodified Debian source (O0 and O2) and runs the 600 sets for you.

| | |
|---|---|
| `inputs/sets.tsv` | 600 sets (300 DNA, 300 protein; 2-9 sequences; options `-NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -PWMATRIX -QUICKTREE`), written by `make_inputs.py` (deterministic: the committed file is reproduced exactly) |
| `run_sets.py` | runs `clustalw -TYPE=<t> -ALIGN -CLUSTERING=UPGMA -OUTPUT=FASTA -OUTORDER=INPUT -NEWTREE=<f>` with the O0 build (O2 when O0 fails), records rows, the `.dnd` guide tree, exit status, the binary used |

## Run it (Apple Silicon: x86-64 image, emulated; docker must be running)

```bash
cd docs/Validation/review-2026-09/handoff/A2-14-clustalw-upgma
./run.sh                  # builds the image (a few minutes: two builds of ClustalW), runs 600 sets, zips ./out
ls -l returned/           # A2-14-clustalw-upgma-out.zip
```

`run.sh` needs only `docker` (and `zip` or `tar`). Nothing aborts: a set on which a binary fails is recorded with its exit status. No account data; the source download uses the User-Agent `Seqeron-review/1.0`.

## What to send back

`returned/A2-14-clustalw-upgma-out.zip` (or `.tar.gz`), kept in this repository under `docs/Validation/review-2026-09/handoff/A2-14-clustalw-upgma/returned/` (`out/ENV.txt`, `out/build_O0.log`,
`out/build_O2.log`, `out/preflight.txt`, `out/fetch_src.sha256`, `out/clustalw_upgma_oracle.tsv`, `out/run.log`).

## How the next session integrates it

`out/clustalw_upgma_oracle.tsv` (rows + `.dnd` text per set) becomes `tests/.../TestData/Alignment/clustalw_upgma_oracle.tsv` and a test `SequenceAligner_ClustalWUpgma_Tests` (rows equal for every set the
program finished; the two-sequence case is expected to come out UNALIGNED: `RootedClusterTree::treeFromDistMatrix` saves no alignment step for `nSeqs == 2`; the `.dnd` text may be compared with a port of
`RootedTreeOutput::printPhylipTree`). Differences are fixed in the port stage by stage (join order and tie rules of `getNodeWithMinDist`, `Node::merge` heights, `RootedGuideTree::calcSeqWeights`).
If the O0 and O2 builds disagree, nothing is integrated until that is explained (undefined behaviour in the source).
