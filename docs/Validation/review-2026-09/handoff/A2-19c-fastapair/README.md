**STATUS: INTEGRATED in WP75 (B09 F346 - F349; archive stored in `returned/A2-19c-fastapair-out.zip`: 180 sets, `results.tsv`, per set `direct.json.gz`, `fastalog.tar.gz`, `stderr.gz`; converter `scripts/alignment/convert_fastapair_returned.py` -> `fasta34_align_oracle.tsv.gz` (886 raw fasta34 calls) and `mafft_fastapair_oracle.tsv.gz` (180 sets); tests `SequenceAligner_MafftFastaPair_Tests`: 886 / 886 raw outputs equal byte for byte, dndfast7 179 / 180 and the complete `mafft` run 179 / 180 equal, the 180th (set 161) is refused: the program reads an unassigned variable). The success path of the nucleotide reader has one example in this run: more sets are in the package `A2-19c-fastapair2`. Was: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP71; result to be returned as `returned/A2-19c-fastapair-out.zip`). The returned battery of `A2-19c-fasta34-battery` (340 raw `fasta34 -m10` outputs) does NOT suffice for `mafft --fastapair` / `--fastswpair`: these options run the program `dndfast7`, which calls `fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup` (no `-H`, `-d` = ALL alignments, `-A` = the Smith-Waterman variant of `--fastswpair`, ktup 1 for protein) and reads the ALIGNMENT BLOCKS of that output (`>>+==========+n`, `; fa_opt:` / `; sw_s-w opt:`, `_overlap`, `_start`, `_stop`, `_display_start`, the aligned sequences) with `ReadFasta34m10` / `ReadFasta34m10_nuc` (`io.c` L3552 - L3825) into local homologies and the file `hat3`; the battery calls used `-H -d0` (list of scores only) and no `-A`. The package gives the raw output of exactly those calls and the files `dndfast7` writes (`hat2`, `hat3`).**

# TEMPORARY hand-off package A2-19c-fastapair (for the user's Mac)

**Why.** B09 item A2-19c-s5: `mafft --fastapair` (FASTA 3.4 optimised scores + local alignments as pair library, strategy "F-") and `--fastswpair` ("H-", fasta34 `-A`, Smith-Waterman). They need fasta34 in the PATH
(the script's own check "== Install FASTA ===" is what the port prints today, equal to an installation without fasta34). To port `dndfast7` and the alignment part of FASTA 3.4 (`do_walign`, `calcons` / the `-m10` alignment writer, `-A`) and verify them, the raw
text of every call and the resulting `hat3` are needed.

**What it does.** For each of the 180 sets of `inputs/sets.tsv` (2 - 40 protein / DNA sequences of 6 - 300 residues; duplicates, unknown letters, low complexity, lower case; options `--fastapair`, `--fastswpair`, `--retree 1`, `--maxiterate 2/3`, `--treeout --distout`, `--reorder`,
`--anysymbol`, `--op --ep`, `--nofft`): (1) `mafft <options> in.fa` WITHOUT `--quiet`, with `FASTA_4_MAFFT` pointing to a logging wrapper (`fasta34-log`) that keeps, for every call of fasta34, the arguments, the query file, the library file and the RAW output
(`out/sets/<id>.fastalog.tar.gz`); rows, names, return code, standard error (`out/sets/<id>.stderr.gz`), `in.fa.tree` and `in.fa.hat2` (`out/results.tsv`); (2) `dndfast7 [-A] < in.fa` run directly in an empty directory: its `hat2`, `hat3` and standard error (`out/sets/<id>.direct.json.gz`).
Debian 12 `mafft` 7.505-1, `fasta34` built from `../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z` (same as the battery; the shar is copied by `run.sh`, parsed and not executed by `fetch_fasta34.py`). No account data, no secrets.

## Run it (Apple Silicon: the image is x86-64 and is emulated; 20 - 40 minutes)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-fastapair
./run.sh
ls -l returned/          # A2-19c-fastapair-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-19c-fastapair/returned/
```

## How the next session integrates it

1. `scripts/alignment/convert_fastapair_returned.py` (to be written then) embeds `fastapair_oracle.tsv.gz` (sets, rows, names, tree / hat2, standard error without the process-number line, `hat2` / `hat3` of the direct runs, every raw fasta34 call normalised).
2. Port: the alignment half of FASTA 3.4 (`do_walign` with `FLOCAL_ALIGN` / banded alignment, the `-m10` alignment block writer of `pro_m10` / `calcons`, `-A` = `smith_waterman`), `ReadFasta34m10[_nuc]`, `cutal`, `putlocalhom2`, `dndfast7` (`hat2` from `mtx2`, `hat3`), the `tbfast` call with `hat3`; strategies "F-" / "H-" of the script. Test: `hat2`, `hat3`, rows, standard error of every set equal; the raw fasta34 text of the port's alignment writer equal to the logged output (the `MafftFasta34` list part is already equal on 340 calls: `SequenceAligner_MafftFasta34Battery_Tests`).
3. `MafftProgram.Run` then runs `--fastapair` / `--fastswpair` (done in WP75: `SequenceAligner.MafftFasta34Align.cs` / `MafftFasta34Program.cs` / `MafftDndfast7.cs`; the script's "Install FASTA" message is gone: the installation of the oracle has fasta34).
