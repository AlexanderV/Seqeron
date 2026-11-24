**STATUS: RETURNED 2026-11-24 (archive stored in `returned/A2-12b-s6-out.zip`). RESULT: only DaliLite 3.3 was obtained (sha256 85f6dd8d...); every `www.ebi.ac.uk/dali/DaliLite/*` URL answers HTTP 200 with the same 26 344-byte HTML page (a soft 404, sha256 fb1872b5...), the Internet Archive CDX returned 0 captured URLs for all 12 prefixes (one request failed) - the CDX is known to be unreliable (it also returned 0 lines when the 3.3 snapshots existed), so this is NOT yet a proof that no older release exists; follow-up package with direct snapshot probes = `../A2-12b-s6b-dalilite-snapshots` (created in WP76, F351, DOABLE-awaiting-user). Was: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP71). Question: is there an older DaliLite release (v2.x / v1) whose `DaliLite -pairwise pdb1.pdb pdb2.pdb` prints the alignment as `Query` / `Sbjct` lines to STANDARD OUTPUT (the only place the glue `tc_generic_method.pl` (`seq2pdb_pair`, `-method=DaliLite`) reads them from; DaliLite 3.3 never does: 0 of 272 logged calls of run 10 of `A2-12b-dali-sap-stderr`, its `DaliLite` script writes `index.html` / `molXA.html` only)? What this machine could establish, with the sandbox's network: the live server of the authors serves ONLY `DaliLite_3.3.tar.gz` (sha256 85f6dd8d... = the Internet Archive copy used by runs 4 - 10), `DaliLite.v4.tar.gz` and `DaliLite.v5.tar.gz`; HEAD requests for DaliLite_{1, 1.0, 2.0 ... 2.4.5, 3.0 - 3.2}.tar.gz, `.tgz`, `DaliLite.v{1..4}.tar.gz` under `http://ekhidna2.biocenter.helsinki.fi/dali/` and `.../dali_lite/downloads/{,v1,v2,v3}/` all answered 404; `ftp.ebi.ac.uk/pub/databases/dali/` holds a README only, `ftp.ebi.ac.uk/pub/software/` has no DaliLite; `web.archive.org` answers 403 to this machine (the sandbox's proxy), so the Internet Archive's index of older releases could not be read here. The CHANGES file of DaliLite 3.3 lists releases back to v2.3 (30 Aug 2004) and says of v2.4 "renamed links in -pairwise HTML output page": every release back to 2.3 wrote HTML files too. This package asks the archive's index (CDX) and downloads what it lists. NOT yet integrated.**

# TEMPORARY hand-off package A2-12b-s6-dalilite-old (for the user's Mac)

Fetch-and-inspect only (a small Python image: nothing that is downloaded is compiled or run). It (1) asks the live servers for the release names `DaliLite_<v>.tar.gz`, (2) asks the Internet Archive's CDX index for every captured archive under the distribution places of DaliLite
(`ekhidna.biocenter.helsinki.fi/dali_lite|dali|dali_server|downloads`, `www.ebi.ac.uk/DaliLite`, `ftp.ebi.ac.uk/pub/databases/dali`, `www.bioinformatics.org/project`, ...), (3) downloads each candidate (resumed curl, neutral User-Agent `Seqeron-review/1.0`) and reads
`CHANGES`, `README`, `INSTALL`, `Makefile`, the script `DaliLite` and the Perl scripts as TEXT, (4) writes `out/report.txt` with one line per archive (sha256, size, newest release named in CHANGES, whether the `DaliLite` script contains the strings `Query` / `Sbjct`).

## Run it (Apple Silicon: the image is x86-64 and is emulated; it only runs Python and curl)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-12b-s6-dalilite-old
./run.sh                  # a few minutes (downloads may be large: the archive's copies of the data archives are skipped by name: only names that contain "dali" and end in .tar.gz / .tgz / .tar / .zip / .Z / .gz are fetched)
ls -l returned/           # A2-12b-s6-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-12b-s6-dalilite-old/returned/
```

## What to send back

`returned/A2-12b-s6-out.zip`: `out/fetch.log` (every URL asked and the answer), `out/report.txt`, `out/inspect/NN/` (the text files read out of each archive), `out/ENV.txt`. An empty candidate list is a result too (it is the proof that the archive has no older release).

## How the next session integrates it

* No candidate / no archive older than 3.3 / no script with `Sbjct`: item closed as BLOCKED-FOR-EVERYONE with `out/fetch.log` as the proof (the `Query` / `Sbjct` branch of the glue then only ever sees the EBI web service's HTML, `daliweb_pair`, which needs the network: DOABLE separately).
* A release whose `DaliLite` script prints `Query` / `Sbjct` to standard output: the next session writes the build package (same repairs as `A2-12b-dali-sap-stderr/fetch_dali3.py`: install path, libg2c shim, logging wrapper) and runs the 60 `dali_pair` sets with it; `StructureAlignerOutputs[DaliPair]` then carries lines for `TCoffeeDaliLiteAlignment` and the 60 rows lock the parse branch.
