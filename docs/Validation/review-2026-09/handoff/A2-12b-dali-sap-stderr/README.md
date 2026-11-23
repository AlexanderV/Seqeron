**STATUS: RUN 9 RETURNED, RUN 10 NEEDED (2026-11-24, returned/A2-12b-out.run9.zip). Run 9: the source rebuild still fails in the Fortran sources (the Makefile's `FC = g77`; `util.f` lines are cut at column 72, `-ffixed-line-length-none` had no effect: the Makefile rule must use its own variable), but the `libg2c.so.0` shim worked: `ldd -r` finds all symbols of `puutos`, and of `puu` / `dalicon` all but `G77_fdate_0` (g77 intrinsic FDATE; now added to the shim). Yet every call (137) still ends `sh: 1: .../Bin/puu: not found` (also `puutos`, which has no unresolved symbol): the cause is not visible from the logs. The package now writes `out/diag_puu.txt` first (file / interpreter / NEEDED / ldd -r / direct run / run through sh / strace of the exec of `puu`, `puutos`, `dalicon`, `dssp`). RUN 10 = two commands: `git pull`, `cd docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr`, then `./diag.sh` (builds the image, ~15 min, 1 minute of diagnostics, writes `returned/A2-12b-diag.zip`) and afterwards `./run.sh` (the image is cached now: only the run, ~50 min, writes `returned/A2-12b-out.run10.zip`). Return BOTH archives.**







```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr
docker build --no-cache --platform linux/amd64 -t seqeron-a212b . && ./run.sh        # ~15 min; Apple Silicon: emulated amd64
ls -l returned/        # A2-12b-out.run7.zip : commit it under docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr/returned/
```

(`./run.sh` builds the image itself, so the explicit `docker build` only forces the DaliLite fetch / build layer to be rebuilt with the new wrapper.) Integration: B09.md item **A2-12b-s5** (the `.saved` files and the `call_*.txt` of `out/s4`, the t_coffee rows of the 60 sets: lock the `Query` / `Sbjct` parser and the library entries of a successful DaliLite against the real rows; `TCoffeeDaliLiteAlignment` and the `DaliInstalled` path of `TcPairJobEntries` are ready).

**(History) STATUS: RETURNED, FINAL (2026-11-23, returned/A2-12b-out.run6.zip): the s4 run with DaliLite v3.3 linked at its compiled-in install path worked: forwarder.pl found, 60 sets, 163 logged DaliLite calls (158 return code 0, 5 return code 2), t_coffee rows 26 OK / 34 FAIL. DaliLite 3.3 itself accepts the CA-only structures of the test sets without an alignment (`suspect pdb1.pdb is CA-only`, `Analyse () vs. ()`, no result file), a few calls print `getcwd() failed` / `cannot create dali.lock` (the working directory is removed by t_coffee while the program runs) - this is the real program's behaviour and is the oracle. NO FURTHER RUN NEEDED. Integration = WP68 (s1-s4). Earlier archives (run 1, 3, 4, 5) are kept as returned data.**

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
