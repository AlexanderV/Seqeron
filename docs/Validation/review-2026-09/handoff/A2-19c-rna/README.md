**STATUS (B09 WP86, F377): RUN 3 INTEGRATED (150 sets, all deterministic: `--contrafold` 16 sets rc 0 equal the port byte for byte after the port of `contrafoldwrap` + CONTRAfold 2.02 + the RNA step of `tbfast`; `--dafs` 8 sets rc 1 "ERROR IN mccaskill_mea" and `--dafspair` 8 sets rc 1 "Cannot open hat4" equal the port; `--scarnapair` 18 sets rc 1, `--mccaskill` protein 1 set rc 1 equal the port; the 73 `--mccaskill` + 25 `--rnaalifold` + 1 protein-abort sets are the open port, DOABLE: `mccaskillwrap` and the RNA step of `tbfast` are ported and verified with a sandbox oracle, what is missing is `mxscarnamod` itself (licence: owner decision) and the RNAalifold partition function). Run 3 stays as `returned/A2-19c-rna-out.run3.zip`. RUN 4 NEEDED (DOABLE-awaiting-user): the `dafs` built from the public release 0.0.3 has no `--mafft-out` / `--mafft-in` (the sets stop in `mccaskillwrap -G`; the stderr of `dafs` was not kept: the package now records `dafs --help`, the failing call and `strings dafs | grep -i mafft` in `out/build_rnaprogs.log`), and the port of the CONTRAfold numerics (single precision, polynomial `Fast_Exp` / `Fast_LogExpPlusOne`) and of MXSCARNA needs the programs' own output: `out/probes.json` (111 sequences of 1 - 340 letters through `contrafold predict IN --posteriors 0.01 OUT`, `contrafold predict IN`, `mxscarnamod -m -writebpp IN`). User command (Apple Silicon: `--platform linux/amd64`, emulated; about 15 minutes):**
```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-rna
./run.sh                  # builds the image, runs the 150 sets twice and the 111 probes, zips ./out into returned/A2-19c-rna-out.zip
ls -l returned/           # commit returned/A2-19c-rna-out.zip (runs 1 - 3 stay as .run1 / .run2 / .run3.zip)
```

**STATUS: RUN 3 RETURNED (`returned/A2-19c-rna-out.zip`: CONTRAfold 2.02 and DAFS built; --contrafold 16 sets rc=0, --dafs/--dafspair 16 sets still rc=1, rest as run 2; all deterministic); integration pending (next B09 WP). Runs 1 and 2 stay as `.run1.zip` / `.run2.zip`.**

**STATUS: RUN 2 INTEGRATED (B09 WP83, F374: 150 sets; the 51 stops with return code 1 equal the port, 98 computed runs refused by the port = DOABLE port); RUN 3 NEEDED: the image now also builds CONTRAfold 2.02 and DAFS (`build_rnaprogs.sh`), so that the 32 `--contrafold` / `--dafs` / `--dafspair` sets (stopped at "Install CONTRAfold" / "Install DAFS" in run 2: neither program is part of the MAFFT archive) really run. User command (Apple Silicon: the scripts use `--platform linux/amd64`, emulated; about 10 minutes with the two extra builds):**
```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-rna
./run.sh                  # builds the image, runs the 150 sets twice, zips ./out into returned/A2-19c-rna-out.zip (also after an error or Ctrl-C)
ls -l returned/           # commit returned/A2-19c-rna-out.zip (run 2 stays as returned/A2-19c-rna-out.run2.zip)
```
**STATUS: RUN 1 RETURNED 2026-11-25 (`returned/A2-19c-rna-out.run1.zip`): the image build failed - `E: Package 'vienna-rna' has no installation candidate` on debian:12 (the package is in `non-free`), which aborted the single `apt-get install` (no curl, no mafft either), so nothing ran. PACKAGE FIXED IN B09 WP81 (F370): the base packages are installed one by one (a missing package is logged, not fatal), ViennaRNA 2.5.1 (the version of the Debian package; fallbacks 2.4.18, 2.6.4, then `non-free`) is built from the upstream release with `--without-perl --without-python --without-swig --disable-openmp --disable-lto` and `-O2` (no `-march=native`) by `build_vienna.sh`, every step logged, the build never fails, the zip is always written; the container script was run end to end here in `RNA_QUICK=1` mode with the installed mafft 7.505 and vienna-rna 2.5.1 (6 sets, `--rnaalifold --localpair` rc 0). RUN 2 NEEDED (DOABLE-AWAITING-USER): `cd docs/Validation/review-2026-09/handoff/A2-19c-rna && RNA_QUICK=1 ./run.sh` (optional check), then `./run.sh`, commit `returned/A2-19c-rna-out.zip`.**

**STATUS: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP78; result to be returned as `returned/A2-19c-rna-out.zip`).** The RNA options of the script `mafft` (`--mccaskill`, `--scarnapair`, `--rnaalifold`, `--dafs` / `--dafspair`, `--contrafold`, `--larapair`, `--foldalignpair`) run programs that the Debian binary package `mafft` 7.505-1 does not have: the installed script stops with "Install MXSCARNA" / "Install DAFS" / "Install CONTRAfold" (reproduced by the port, equal on every case of the command oracle), and `--rnaalifold` runs the progressive stage up to the first RNA profile step and ends with `sh: 1: RNAalifold: not found` / `Cannot open /tmp/_rnaalifoldin` (RNAalifold of ViennaRNA is not installed here). Their real behaviour (the base-pairing probabilities of McCaskill's algorithm inside tbfast / dvtditr, `mxscarnamod`, the RNA profile scores `-o -c -r -R -e`) can only be observed with the extended build and ViennaRNA.

# TEMPORARY hand-off package A2-19c-rna (for the user's Mac)

**What the package does.** The image (`Dockerfile`, Debian 12, `--platform linux/amd64`) installs ViennaRNA (`vienna-rna`: `RNAalifold`) and builds MAFFT 7.505 *with its extensions* from the official source archive `mafft-7.505-with-extensions-src.tgz` (`build_ext.sh`: `core` and `extensions`, `make` / `make install` into `/opt/mafft-ext`; every step is logged and never fails the image build; the programs that exist are listed in `out/versions.txt`). `make_inputs.py` (fixed seed) writes `inputs/sets.tsv`: 150 sets, 144 RNA families with stems (2 - 8 sequences) for 18 option lines (`--mccaskill` with `--globalpair` / `--localpair` / `--genafpair` / `--scarnapair` / `--maxiterate`, `--scarnapair`, `--rnaalifold`, `--dafs`, `--dafspair`, `--contrafold`, output options), and DNA / protein sets for the refusals. `in_container.sh` runs each set twice with `mafft <options> in.fa` (the extended build when it exists) and records return code, standard output, standard error, files and whether both runs are equal.

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-rna
RNA_QUICK=1 ./run.sh     # optional: 6 sets, a check of the package
./run.sh                 # the 150 sets
ls -l returned/          # A2-19c-rna-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-19c-rna/returned/
```
`run.sh` never stops on an error inside the container (`set +e` around `docker build` / `docker run`): whatever happens, `out/` is zipped into `returned/A2-19c-rna-out.zip` (steps in `out/steps.log`, the build log in `out/build_ext.log`, the programs and their checksums in `out/versions.txt`). `RUN_TIMEOUT=300` (seconds per run) can be raised.

## How the next session integrates it

1. A converter (`scripts/alignment/convert_rna_returned.py`, to be written with the first result) turns the deterministic runs into `TestData/Alignment/mafft_rna_oracle.tsv.gz` (the format of the other mafft oracles).
2. `tbfast` / `dvtditr` / `disttbfast` get the RNA branch the options switch on: `rnaprediction` ('m' = McCaskill inside the program: `rna.c`, `mccaskillwrap.c` of the extensions; 'r' = `RNAalifold -p`), the scoring of RNA profiles (`-o -c -r -R -e -s`), the base-pairing probabilities as the `rnakozo` tables of `Salignmm` / `MSalignmm` / `partSalignmm`; `--scarnapair` is `mxscarnamod` through `tbfast _ ... -s`.
3. Tests: `SequenceAligner_MafftRna_Tests` (`Oracle_HasEveryGroup`, `Group_EqualsMafft`), then `MafftProgram.Run` / `mafft_command` stop answering `--mccaskill` etc. with the installed refusal (the refusal stays for an installation without the programs: the reference installation).
4. If a program of the extensions cannot be built (see `out/build_ext.log`), the option stays registered as DOABLE-awaiting-user with that log as the proof; `--dafs` needs DAFS (http://www.ncrna.org/, not part of the MAFFT archive), `--contrafold` needs CONTRAfold, `--larapair` LaRA, `--foldalignpair` FOLDALIGN: if the archive does not contain them the package records the script's own "Install ..." text, which is what the port already reproduces.
