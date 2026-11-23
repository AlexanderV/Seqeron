**STATUS: INTEGRATED (2026-11-23, B09.md F296)** - `returned/A2-12-out.zip`: the 130 `sap_real` rows are in `tcoffee_external_oracle.tsv` and locked by `SequenceAligner_TCoffeeSap_Tests` (glue `sap_pair` ported); the DaliLite v5 probe failed on a missing `dali.conf` and asked the wrong question (the glue calls DaliLite v3's `DaliLite -pairwise`): continued in `handoff/A2-12b-dali-sap-stderr` (AWAITING-USER).

# TEMPORARY hand-off package A2-12 - the real SAP 1.1.3 and DaliLite v5 next to t_coffee 13.41 (Docker, for the user's Mac)

**STATUS: AWAITING THE USER'S RUN (2026-11-19).** Item A2-12 of `B09.md` (`sap_pair`, `fugue_pair`, `dali_pair`).

**What the audit machine does, and what is already done.** SAP, DaliLite and `fugueali` are not installed there (not in apt, PyPI or the reachable indexes; the sandbox also refuses to run
a downloaded / compiled binary - this is why SAP and DaliLite need this package). T-Coffee 13.41 reacts to the absence in three different ways; all three are reproduced and tested
(`SequenceAligner_TCoffeeExternal_Tests`, oracle `tcoffee_external_oracle.tsv`, features `sap_pair`, `dali_pair`, `fugue_pair`):

| method | installed program | port |
|--------|-------------------|------|
| `sap_pair`, `-mode 3dcoffee` (= `-in Msap_pair -template_file SELF_P_`) | "SAP is not installed. TMalign will be used instead": the method becomes `pdb_pair` + TMalign | `TCoffeePairMethod.SapPair` = `TMalignPair` glue (rows equal) |
| `dali_pair` | stops when it loads the method ("ERROR: UNSPECIFIED"), exit status 1 | `TCoffeeAbortedException` |
| `fugue_pair` | "The Program fugueali Needed by T-COFFEE Could not be found [FATAL]", exit status 1 | `TCoffeeAbortedException` |

FUGUE (`fugueali`, Shi, Blundell & Mizuguchi 2001) is licence-gated: `install.pl` of T-Coffee itself lists its source as "empty" and its address
`http://www-cryst.bioc.cam.ac.uk/fugue/download.html` is not reachable; neither PyPI (`fugue` is the unrelated fugue-project), conda (bioconda has no fugue), the Debian archive nor GitHub
(search `fugueali`, `FUGUE homstrad`) has it: **BLOCKED with that proof** (the committed behaviour is the program's refusal).

**What is still open (this package).** The behaviour with the REAL programs: (1) real SAP: `sap_pair` of `util_dp_drivers.c` (L2049: `normalize_pdb_file`, the `<template1>_<template2>.sap_results`
cache, `is_sap_file`, the weight `WEIGHT 100`) is a different glue from `pdb_pair`; (2) real DaliLite: the glue of `tc_generic_method.pl` (`DaliLite -pairwise pdb1.pdb pdb2.pdb`, the `Query` / `Sbjct`
lines of the output) expects the DaliLite v3 command line, DaliLite v5 (the only source that is reachable: https://github.com/liisaholm/DaliLite.v5) has another one (`import.pl`, `dali.pl`).

## What it does

1. `Dockerfile` (Ubuntu 24.04, linux/amd64) installs the same `t-coffee` 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 package as the audit machine plus `tm-align`, `mustang`, `ncbi-blast+`, builds **SAP 1.1.3**
   from the GitHub archive `https://github.com/mathbio-nimr-mrc-ac-uk/SAP/archive/v.1.1.3.tar.gz` (sha256 `1ee5025f8a900cd9d9c490f7038b98d80a619e3015f2dc97b869ea3033c459e0`, the source of the bioconda
   recipe; `./bootstrap && ./configure && make && make install`) and **DaliLite v5** (`git clone` of the GitHub repository, `cd bin && make`, the commit is recorded). Builds never abort the image: logs and
   `BUILD_FAILED` markers are returned.
2. `in_container.sh` (a) records the environment (`out/ENV.txt`: date, `uname -m`, t_coffee version, the SAP binary's sha256, the DaliLite commit); (b) runs the feature `sap_real` of
   `scripts/alignment/generate_tcoffee_external_oracle.py` (the SAME synthetic CA-only structures as the committed TMalign / mustang oracle: random-walk backbones, indels, substitutions, rigid motions) with
   `-method sap_pair`, `-mode 3dcoffee` and `-method sap_pair,slow_pair` and the REAL `sap` in t_coffee's plug-in directory; every row records the alignment t_coffee wrote and the raw `*.sap_results` files of its
   cache (SAP's own output); (c) `dali_probe.py` runs DaliLite v5's own `import.pl` / `dali.pl` on 24 pairs of the same kind of structures and records every output (a refusal of CA-only models is a result too).

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
cd docs/Validation/review-2026-09/handoff/A2-12-sap-dali
./run.sh                  # builds the image (SAP and DaliLite are compiled: a few minutes), runs the sets, zips ./out
ls -l returned/           # A2-12-out.zip
```

## What to send back

The archive `returned/A2-12-out.zip` (or `.tar.gz`) placed in this repository under `docs/Validation/review-2026-09/handoff/A2-12-sap-dali/returned/`
(`out/ENV.txt`, `out/sap_build.log`, `out/dali_build.log`, `out/tcoffee_external_oracle.tsv`, `out/dali_probe.tsv`, `out/gen_sap.log`, `out/dali_probe.log`, `out/run_console.log`). A run that prints a failed step is still useful.

## How the next session integrates it

* `out/tcoffee_external_oracle.tsv` (feature `sap_real`: rows of `t_coffee` with the real SAP + `extra.sap` = SAP's raw output per pair of templates) becomes the oracle of a port of `sap_pair` (util_dp_drivers.c L2049: the
  parser of the `.sap_results` files and the library weights); the committed `SapPair` (= TMalign fallback) then only applies when SAP is absent, and the port is verified row for row against the returned rows.
* `out/dali_probe.tsv` decides whether a DaliLite v5 shim for `DaliLite -pairwise` can be built (and from which output the `Query` / `Sbjct` parser of `tc_generic_method.pl` is fed); a refusal of the synthetic models
  means the backbone-complete variant of the generator is needed first.
* If `ENV.txt` shows another t_coffee than 13.41.0.28bdc39 nothing is integrated until that is explained.
