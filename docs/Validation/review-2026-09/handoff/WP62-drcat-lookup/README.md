# Hand-off WP62-drcat-lookup (temporary Docker package for the user's Mac)

**Why.** A USA `db:id` of a database that the configuration does not define is looked up by `queryDbfind` (`ajquery.c`) in the DRCAT resource catalogue (`DRCAT.dat` of the EMBOSS data) and fetched over the network from the URL / method of
the entry. The audit machine has no network access to the databases, so the lookup (the messages, the fetched entries, the failure modes) can only be observed from a machine that has.

**What it does.** Builds a `linux/amd64` image (Ubuntu 22.04, package `emboss` 6.6.0+dfsg, the version the port is built against; no account data) and runs the 20 commands of `rows.txt` (seqret / entret / showdb / drfinddata / dbtell on public databases,
unknown entries and unknown databases, wildcards) with a 60 s limit each. Every command writes standard output, error stream, exit status and written files to `out/rows/NNN/`; `out/environment.txt` holds `uname -m`, the program version, the checksum
of `DRCAT.dat` and the address of the server (to show that the network was there).

**Commands.**
```
cd docs/Validation/review-2026-09/handoff/WP62-drcat-lookup
./run.sh                      # builds the image, runs the rows, writes ./out and ./returned/WP62-drcat-lookup-out.zip
ls -l returned
```
**Send back.** The archive `returned/WP62-drcat-lookup-out.zip` (or the `out` directory) placed in `docs/Validation/review-2026-09/handoff/WP62-drcat-lookup/returned/`.

**How it is integrated.** The next session reads `out/rows` into oracle rows of `seqret` / `entret` for the USAs of DRCAT databases (the fetched entries become test files, the messages of failures decide the messages of the port) and implements
`queryDbfind` / the URL access methods of `ajquery.c` against them; the oracle of `seqret` (`emboss_prog_seqret...`) gets the rows that run without network once the catalogue is a data file of the port.
