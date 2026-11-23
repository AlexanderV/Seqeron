**STATUS: RUN 4 RETURNED (2026-11-23, returned/A2-12b-out.run4.zip): DaliLite v3.3 (`DaliLite_3.3`, tarball sha256 85f6dd8d..., 1 367 393 bytes, resumed with curl from the Internet Archive) was downloaded and BUILT (`BUILD_OK`); s1/s2/s3 re-run. s4 (`dali_pair` with the real DaliLite v3.3) crashed in the generator, not in DaliLite: the jobs run in parallel threads and shared one `/tmp/dali_log`, so a log file was removed between `listdir` and `open` (FileNotFoundError). Fixed (the logging wrapper writes into `$DALI_LOG_DIR`, one directory per job; generator reads it). ONE MORE RUN for s4 only: `git pull`, then `docker build --no-cache -t seqeron-a212b . && ./run.sh` (the image is rebuilt from cache layers except the fetch step; ~25 min); return `returned/A2-12b-out.zip`.**

# TEMPORARY hand-off package A2-12b - what is left of item A2-12 after the user's first run (Docker, for the user's Mac)

**What the first run (`handoff/A2-12-sap-dali`, integrated 2026-11-23, F296) decided.** The glue `sap_pair` of t_coffee 13.41 with the REAL SAP 1.1.3 is ported and equals the program on 130 / 130 sets
(`SequenceAligner_TCoffeeSap_Tests`). Two things stay open and need a second run:

1. **Why the program stops.** In every oracle set in which a pair of two usable structures has no SAP result file (5 sets) t_coffee ended with exit status 1 and no alignment; by the source (`util_dp_drivers.c`
   L2144) a missing result only produces "SAP failed to align" and the job returns. The port reproduces the observed stop (`TCoffeeAbortedException`), but the statement that comes out of the program's
   own message is missing. This package re-runs the same 130 sets (same seeds) with t_coffee's stderr / stdout recorded (`extra.stderr`, `extra.stdout`).
2. **DaliLite.** The 24-pair probe of the first run failed for a trivial reason (`dali.conf` missing: DaliLite v5 `dali.pl` stopped with "Config file not found"), and it was the wrong question anyway:
   the glue of t_coffee (`tc_generic_method.pl`, `seq2pdb_pair`) calls **`DaliLite -pairwise pdb1.pdb pdb2.pdb`** and parses the `Query` / `Sbjct` lines - the command line of DaliLite **v3**. The official
   DaliLite v5 tarball (`http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz`, listing read by the audit: `bin/` holds `import.pl`, `dali.pl`, `dccp2dalicon.pl`, ... and NO `DaliLite`) has no such command.
   This package (a) runs the 60 committed `dali_pair` sets with v5's `bin/` on the PATH (expectation: exit status 1 exactly like the committed rows, i.e. nothing changes) and (b) tries to fetch DaliLite v3
   from the Internet Archive snapshot of its download page (`web.archive.org` is blocked on the audit machine; the live v3 download page is 404) and, if it builds, runs the same sets with the real v3 behind a
   logging wrapper (`extra.dali` = arguments, return code, stdout and stderr of every `DaliLite` call). If v3 cannot be fetched the item is closed as BLOCKED-FOR-EVERYONE with the fetch log as proof.
3. **`-mode 3dcoffee`** (= `-in Msap_pair -template_file SELF_P_`): the 8 committed sets could not be completed offline (name look-ups of the profile templates); the user's machine has a network.

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
cd docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr
./run.sh                  # builds the image (SAP and DaliLite v5 are compiled, DaliLite v3 is fetched and compiled if the archive has it), runs the sets, zips ./out
ls -l returned/           # A2-12b-out.zip
```

## What to send back

`returned/A2-12b-out.zip` placed in this repository under `docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr/returned/`: `out/ENV.txt`, `out/dalilite3_fetch.log`, `out/dali5_*`, `out/s1/` (sap_real with stderr),
`out/s2/` (3dcoffee, real SAP), `out/s3/` (dali_pair, DaliLite v5 on the PATH), `out/s4/` (dali_pair, DaliLite v3; only when it was built), `out/run_console.log`. A run that prints a failed step is still useful.

## How the next session integrates it

* `out/s1`: the rows must equal the committed `sap_real` rows (same seeds; a difference means the SAP run is not deterministic and is analysed first); `extra.stderr` of the 5 stopped sets names the stop and replaces the
  observed-rule comment of `TcSapPairEntries`.
* `out/s3`: appended to `tcoffee_external_oracle.tsv` as feature `dali_v5_path` (rows must equal the committed `dali_pair` rows: still "stops").
* `out/s4` (if present): the port of `DaliPair` for the DaliLite v3 glue (`Query` / `Sbjct` parser of `seq2pdb_pair`, `.` -> `-`) is verified row for row against it; otherwise the fetch log is the BLOCKED proof.
* `out/s2`: `-mode 3dcoffee` with the real SAP becomes a supported mode if the rows are reproducible by the port.
* If `ENV.txt` shows another t_coffee than 13.41.0.28bdc39 nothing is integrated until that is explained.
