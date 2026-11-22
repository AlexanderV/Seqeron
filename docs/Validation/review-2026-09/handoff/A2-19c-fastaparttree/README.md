# Hand-off A2-19c-fastaparttree: `mafft --fastaparttree` (MAFFT 7.505, `splittbfast` with the pairwise scores of `fasta34`)

**Why.** The port of MAFFT 7.505 (B09 F283 - F287) covers the progressive, iterative, pair-library and `--auto` methods. The tree-splitting methods `--parttree`, `--dpparttree` (`splittbfast.c`) are
registered as DOABLE; `--fastaparttree` needs the program `fasta34` (FASTA 3.4, W. R. Pearson) in the PATH, which the audit machine does not have (it has `fasta36`, a different program: not an oracle).

**What to run (one command).** `./run.sh` (Docker Desktop; the image is `linux/amd64`): builds Debian 12 `mafft` 7.505-1 and `fasta34` (sources: see `fetch_fasta34.py`; if the download fails put
`fasta34*.tar.gz` into `inputs/` and run again), runs the 60 sets of `inputs/sets.tsv` with `mafft --quiet --fastaparttree ...`, writes `out/results.tsv` (id, options, sequences, return code, aligned rows,
SHA-256 of the standard error, its first 200 bytes) and `out/versions.txt`, and zips `out/` into `returned/`.

**Send back.** `docs/Validation/review-2026-09/handoff/A2-19c-fastaparttree/returned/` (the archive). No account data, no secrets: downloads use the neutral User-Agent `Seqeron-review/1.0`.

**What the lead does with it.** Port `splittbfast.c` (`--parttree`, `--dpparttree`: oracle = the installed mafft, locally) and the `fasta34` branch of `dndfast7.c` (`--fastaparttree`: oracle = this package's results), regenerate
`mafft_parttree_oracle.tsv.gz` from `inputs/sets.tsv` with the installed program for the non-fasta variants.
