**STATUS: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP76, follow-up of `../A2-12b-s6-dalilite-old`, whose result is stored in that package's `returned/`). Question: is there an older DaliLite release (1.0 ... 3.2, or an undated `DaliLite.tar.gz` / `DaliLite.zip`) whose `DaliLite -pairwise pdb1.pdb pdb2.pdb` prints the alignment as `Query` / `Sbjct` lines to STANDARD OUTPUT (the only place the glue `tc_generic_method.pl` reads them from; DaliLite 3.3 never does: 0 of 272 logged calls of run 10)? Why a new package: the first package asked the Internet Archive's CDX index and got 0 captured URLs for all 12 prefixes (CDX is unreliable: it also gave 0 when the 3.3 snapshots existed), and `www.ebi.ac.uk/dali/DaliLite/*` serves a 26 344-byte HTML soft 404 for every name. This package does not use the CDX: it probes SNAPSHOTS directly (`https://web.archive.org/web/<YYYY>id_/<url>` for every candidate release name x every candidate base URL x the years 2003 - 2013) and the availability API. From the sandbox of the author of this package `web.archive.org` answers "Connection reset by peer" to every request (the proxy), so the probe has to run on the user's machine. NOT integrated.**

# TEMPORARY hand-off package A2-12b-s6b-dalilite-snapshots (for the user's Mac)

Fetch-and-inspect only (a small Python / curl image: nothing that is downloaded is extracted to disk, compiled or run). What it does:

1. **connectivity** (`out/connectivity.log`): HEAD of `web.archive.org`, the availability API, `www.ebi.ac.uk/dali/DaliLite/`, `ekhidna2.biocenter.helsinki.fi/dali/` (so that a blocked network is visible).
2. **availability API** (`https://archive.org/wayback/available?url=<url>` and the endpoint named in the task, `https://archive.org/wait/available?url=<url>`) for every (base URL, release name): lines `AVAILABLE` in `out/requests.log`.
3. **directory listings** of every base URL for every year 2003 - 2013 through `web/<year>id_/` (`out/listings/*.html`, read as text): every `href` that names an archive with `dali` in it is added to the candidates.
4. **snapshot probes**, release names `DaliLite_<v>.tar.gz`, `DaliLite.v<v>.tar.gz`, `DaliLite-<v>.tar.gz`, `DaliLite_<v>.tgz` for v = 1.0, 1.1, 1.2, 2.0 ... 2.4.5, 3.0 ... 3.3, plus `DaliLite.tar.gz`, `.zip`, `.tgz`, `.tar.Z`, `dalilite.tar.gz`, `DaliLite.v2/v3/v4.tar.gz`, x 19 base URLs (`www.ebi.ac.uk/dali/DaliLite/`, `www.ebi.ac.uk/DaliLite/`, `ekhidna.biocenter.helsinki.fi/dali_lite/(downloads/)`, `ekhidna.../dali/`, `ekhidna2.../dali/(downloads/)`, `ekhidna2.../dali_lite/downloads/`, `www.bioinformatics.org/dali/`, `.../project/`, `ftp.ebi.ac.uk/pub/databases/dali/`, `.../pub/software/dali/` ...): stage A asks three years (2003, 2008, 2013; the archive answers with the capture nearest to the year), stage B asks every year 2003 - 2013 for each URL that got any answer but a plain 404 / 403 / failure. A body is fetched only after a 2xx answer that is not HTML, and kept only if its magic bytes are gzip / zip / compress / bzip2 / tar.
5. **inspection as text**: for every kept archive `out/inspect/<name>/` holds `members.txt`, CHANGES, README, INSTALL, Makefile, the `DaliLite` script, the `*.pl` scripts, and `pairwise_print_lines.txt` (the lines of the script that print `Query` / `Sbjct` / `print` / `echo` / `cat`); `out/report.txt` has one line per archive: sha256, size, members, the newest release named in CHANGES, `script_has_Query`, `script_has_Sbjct`.
6. Every request is logged in `out/requests.log` (id, method, URL, curl return code, HTTP status, final URL, content type, bytes); a failing request never stops the run; the zip is always written.

## Run it (Apple Silicon: the image is x86-64 and is emulated; it only runs Python and curl)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-12b-s6b-dalilite-snapshots
./run.sh                  # builds the image, probes (several thousand small requests, 8 in parallel; it may take 30 - 90 minutes), zips ./out
ls -l returned/           # A2-12b-s6b-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-12b-s6b-dalilite-snapshots/returned/
```

## What to send back

`returned/A2-12b-s6b-out.zip` (`out/requests.log`, `out/fetch.log`, `out/report.txt`, `out/inspect/`, `out/listings/`, `out/avail/`, `out/connectivity.log`, `out/ENV.txt`, `out/docker_build.log`). An empty `report.txt` with a `requests.log` that shows archive answers (HTTP 200 / 302 / 404 from `web.archive.org`, not connection failures) is a result too: it is the proof that the archive holds no older release.

## How the next session integrates it

* No archive older than 3.3, or none whose `DaliLite` script prints `Query` / `Sbjct` to standard output: the item is closed as BLOCKED-FOR-EVERYONE with `out/requests.log` as the proof (the `Query` / `Sbjct` branch of the glue then only sees the EBI web service's HTML, `daliweb_pair`, which needs the network).
* A release whose `DaliLite` script prints `Query` / `Sbjct` to standard output: the next session writes a build package (same repairs as `../A2-12b-dali-sap-stderr/fetch_dali3.py`: install path, libg2c shim, logging wrapper) and runs the 60 `dali_pair` sets with it; `StructureAlignerOutputs[DaliPair]` then carries lines for `TCoffeeDaliLiteAlignment` and the 60 rows lock the parse branch.
