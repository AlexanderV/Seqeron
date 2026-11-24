**STATUS: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP71). Small package: re-runs only the 15 `dali_pair` sets that ended with exit status 1 in run 10 of `handoff/A2-12b-dali-sap-stderr` (the real DaliLite 3.3 never prints `Query` / `Sbjct` to standard output, so the library of those sets is empty), with t_coffee's standard error recorded. Not needed for the byte-for-byte lock of the 60 rows (that is integrated: `SequenceAligner_TCoffeeSapDaliRuns_Tests`, section s4f); it closes the one thing run 10 did not record: the TEXT of the 15 stops in the all-pair-jobs-ran condition (runs 5 / 6 show the same message "You have not provided any method or enough Sequences [FATAL]" for the 12 sets whose pair jobs all ran in them too; sets 18, 29 and 38 had pair jobs missing in runs 5 / 6, so their message is derived from the source, t_coffee.c L5005, not observed).**

# TEMPORARY hand-off package A2-12b-s5-stderr (for the user's Mac)

Run 10 of `A2-12b-dali-sap-stderr` ran the 60 `dali_pair` sets with the real DaliLite 3.3 and t_coffee 13.41 (45 alignments, 15 exit status 1) but `TC_STDERR` was not set, so the standard error of the 15 stops was not kept.
This package runs ONLY those 15 sets (positions 5 8 10 15 16 17 18 20 21 23 24 25 29 32 38 of the job list; same seeds, `FULL_BACKBONE=1`, one worker) with `TC_STDERR=1`.

## Run it (Apple Silicon: the image is x86-64 and is emulated)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-12b-s5-stderr
./run.sh                  # builds the image (DaliLite 3.3 from the authors' server, built and repaired as in runs 7 - 10), runs the 15 sets, zips ./out; about 10 minutes
ls -l returned/           # A2-12b-s5-stderr-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-12b-s5-stderr/returned/
```

## What to send back

`returned/A2-12b-s5-stderr-out.zip` (`out/ENV.txt`, `out/dalilite3_fetch.log`, `out/s4/tcoffee_external_oracle.tsv` with `extra.stderr` / `extra.stdout` / `extra.dali`, `out/run_console.log`).

## How the next session integrates it

`scripts/alignment/convert_tcoffee_sap_dali_returned.py` takes the 15 rows as section `s4g` of run `11`; `SequenceAligner_TCoffeeSapDaliRuns_Tests` asserts for each: exit status 1 (`FAIL:1`), `extra.stderr` contains "You have not provided any method or enough Sequences", the port stops with the same message. A row that is not `FAIL:1` is analysed first (the lock files between the processes of the glue are timing dependent; runs 5 and 6 already showed sets that end differently).
