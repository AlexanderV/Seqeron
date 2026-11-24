**STATUS: INTEGRATED (B09 WP71, F323 / F324, 2026-11-24): `scripts/alignment/convert_fasta34_battery_returned.py` -> `fasta34_battery_oracle.tsv.gz` (340 raw outputs: the port equals all 4 931 printed lines) and `mafft_fastaparttree_full_oracle.tsv.gz` (200 sets: `MafftProgram.Run` equal, rc / rows / names / tree / standard error); tests `SequenceAligner_MafftFasta34Battery_Tests`, `SequenceAligner_MafftFastaPartTreeFull_Tests`. `--fastapair` / `--fastswpair` need alignment blocks this battery does not contain: package `../A2-19c-fastapair`. (Was: RETURNED (2026-11-24, returned/A2-19c-fasta34-battery-out.zip, 9.2 MB: 360 raw `fasta34` outputs (20 probes + 340 calls) in battery/, 200 sets of `mafft --fastaparttree` without --quiet in results.tsv + sets/*.stderr.gz; fasta34 sha256 7491909b...). NO FURTHER RUN NEEDED. Integration = B09 WP71 (compare `MafftFasta34` score by score, embed the sets).)**


```
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-fasta34-battery
./run.sh
```
then commit `returned/A2-19c-fasta34-battery-out.zip` (or `.tar.gz`) to the repository and tell the next session.

# Hand-off A2-19c-fasta34-battery: the raw `fasta34` output (FASTA 3.4, W. R. Pearson) and `mafft --fastaparttree` without --quiet

**Why.** B09 F312 ports `splittbfast -S` (`mafft --fastaparttree`) and the part of FASTA 3.4 that it uses (the `opt` score of `fasta34 -z3 -m10 -p|-n -Q -H -b<n> -E<n> -d0 query library 2|6`). The earlier run
(`../A2-19c-fastaparttree/returned/`) gave 60 sets run with `--quiet` (rows only). This package returns (1) the raw output of `fasta34` on a battery of 340 calls exactly as `splittbfast.c` makes them (the 20 probes of
`getfastascoremtx` + 320 random queries / libraries, protein and DNA, 1 - 60 library sequences of 6 - 400 letters, unknown letters, low complexity, duplicates), so that the port is compared with the program score by score,
and (2) 200 sets (the 60 of the earlier run + 140 more, up to 260 sequences) run with `mafft --fastaparttree ...` WITHOUT `--quiet` (rows, tree file, the full standard error of each run gzip'd).

**What the package contains.** Dockerfile (Debian 12 `mafft` 7.505-1 + `fasta34` built from `../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z`, parsed not executed by `fetch_fasta34.py`; no download is needed), `make_inputs.py`
(deterministic; writes `inputs/battery.tsv`, `inputs/sets.tsv`; already run), `run.sh` (builds, runs, zips `out/`), `in_container.sh`. No account data, no secrets; downloads (only if the shar is missing) use the User-Agent `Seqeron-review/1.0`.

**Output.** `out/versions.txt`, `out/battery.tsv` (id, kind, return code, SHA-256 of the raw output, first 200 bytes of standard error), `out/battery/<id>.out.gz` (the raw `-m10` output of every call), `out/results.tsv`
(id, options, sequences, return code, rows, names, SHA-256 of the standard error, tree file), `out/sets/<id>.stderr.gz`.

**What the next session does with it.** Compares `MafftFasta34` (B09 port) with every raw output (the `opt` of each library sequence), embeds the sets as `mafft_fastaparttree_oracle.tsv.gz` (rows, tree, names; the standard error
contains the process number in its last lines, so it is compared with that line removed) and adds the equality tests.
