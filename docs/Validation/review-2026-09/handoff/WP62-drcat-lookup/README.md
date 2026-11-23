**STATUS: INTEGRATED (2026-11-23, B09.md F309, `returned/WP62-drcat-lookup-out.zip`: 20 rows, EMBOSS 6.6.0.0 on x86_64, DRCAT.dat sha256 `ea870d34...` = the file of the audit machine, network present).** Rows 001-013 (seqret / entret on catalogue databases): the 6 seqret runs that end by themselves fail the fetch with the messages the port reports; 5 runs (embl / uniprot fetches) did not end within the 60 s limit (not a result a port without network can have: documented, asserted as the failed fetch); 2 entret runs used `-outseq`, which entret does not have. Rows 014-020: `showdb`, `drfinddata` (3 queries), `dbtell` (2 databases) equal the port (`EmbossPrograms_Oracle_Tests`, `EmbossDrcat_Oracle_Tests`). No further run.

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

**How it was integrated (F309).** The catalogue (`DRCAT.dat`, `EDAM.obo`, `emboss.standard` and the index files of the installation) is embedded in `Seqeron.Genomics.Infrastructure`; `queryDbfind` / `ajResourceGetDbdata` are ported (a USA of a database that only DRCAT defines is looked up in the catalogue; the fetch fails without network as the library reports it), and so are `drfinddata`, `dbtell` and `showdb`; 899 USA runs of seqret / entret, 233 / 302 / 301 runs of the three programs and the rows of this package are the oracles (`EmbossDrcat_Oracle_Tests`, `EmbossPrograms_Oracle_Tests`).
