# TEMPORARY hand-off package A1-18 / A1-19 - NCBI ORFfinder 0.4.3 oracle (Docker, for the user's Mac)

> **STATUS: DONE (2026-11-18, B09 WP10, F57).** The user ran the package (1750 runs; `ORFfinder -version` 0.4.3, x86_64, Debian 12, ORFfinder.gz sha256
> `f30a0806...09f2` as pinned); the archive is `returned/A1-18-orffinder-out.zip`. It was integrated as `scripts/analysis/generate_orffinder_oracle.py` ->
> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/orffinder_cli_oracle.tsv` and `GenomicAnalyzer_FindOpenReadingFrames_OrfFinderCli_Tests`
> (54 820 results, 232 162 ORFs, all equal). `parse_orffinder_out.py` (draft with guessed header formats) is superseded by that converter - the real
> formats are 0-based `>lcl|ORF<n>_<id>:<from>:<to>` headers (see Evidence GENOMIC-ORF-001). Result: `-c` is a no-op in the binary; `-n`, `-b/-e`, `-ml`
> clamp, `-strand`, partial flags and tables 27/28/31 implemented (INTEGRATION.md section 2 table is the plan; the measured rules are in
> `docs/algorithms/Analysis/Open_Reading_Frame_Detection.md` section 5.5). The package files stay for reproducibility; nothing here needs to be re-run.

**Why.** B09 audit items A1-18 (ORFfinder "Ignore nested ORFs" `-n`) and A1-19 (`-b/-e` region, `-c` circular, `-ml` clamp,
`-strand`, `-s`, `-g`, partial ORFs) can only be decided by *measuring* the official NCBI standalone program: its cross-frame
`-n` filter, the `-b/-e` and `-c` semantics and the `-ml < 30 -> 30` clamp are defined nowhere in the open-source `COrf` library
(see B09 F23/F25). The sandbox the audit ran in refuses to execute a downloaded binary, so this package lets **you** run it once.
The items stay **DOABLE (not ticked)** until the returned output is integrated (see `INTEGRATION.md`).

## What it does

1. `Dockerfile` (base `debian:12.7-slim`, linux/amd64 manifest pinned by digest) downloads the **official** files from
   `https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/` (`linux-i64/ORFfinder.gz`, `USAGE.txt`, `CHANGELOG.txt`), verifies their
   sha256 (ORFfinder.gz `f30a0806...09f2`, 8 386 303 bytes, Last-Modified 2021-05-04), unpacks, runs `ldd`.
2. `make_inputs.py` (committed output in `inputs/`, regenerated and compared inside the container) = 629 FASTA records (244 designed, 170 random, 122 circular, 60 region, 33 edge):
   * `designed.fa` - nesting in frame; nesting / overlap across frames and across strands (so `-n` has something to act on, incl.
     1-nt and 2-nt boundary offsets, alternative-start nesting, equal spans across strands under `-s 2`); ORFs running off the 5' end,
     3' end or both ends; half stop codons at the end; every start codon of tables 1/2/4/11 and every stop of those tables;
     ORFs of 6, 9, ... 90 nt for the `-ml` clamp; low-complexity repeats; every forward design also as its reverse complement;
   * `random.fa` - 170 random sequences (iid at five GC levels, gene-dense mosaics of ORFs on both strands, point-mutated copies, 900-2400 nt);
   * `circular.fa` - 122 genomes whose ORFs (both strands, every codon phase, nested / cross-frame) are cut by the origin, for `-c true`;
   * `region.fa` - 60 sequences of exactly 300 nt for the `-b/-e` grid (ORFs straddling 100/150/200/250);
   * `edge/*.fa` (24 one-record files: length 0..3, 6, 29, 30, lower case, N, IUPAC, U, gaps, CRLF, long header, protein, ...),
     `circular_edge/*.fa` (6 circular genomes without a stop; run with a timeout).
   `inputs/manifest.tsv` lists id, file, category, length and design note of every record.
3. `make_jobs.py` prints the run matrix: **1750 ORFfinder invocations** (`python3 make_jobs.py --summary`): defaults; the factorial
   `-g {1,2,4,11} x -s {0,1,2} x -n {false,true} x -strand {both,plus,minus}` (`-ml 30`) in `-outfmt` 0, 3 (1 and 2 on subsets);
   a `-ml` sweep (-5 ... 1000) with and without `-n`; every genetic code 0..33; `-c true/false` on the circular set and on the
   designed set; a `-b/-e` grid (24 pairs incl. invalid ones) x strand x `-n` x `-s` x `-c`; 116 designed records one-per-process;
   edge files; circular genomes without stop; `-h -help -xmlhelp -version -version-full -dryrun -logfile -conffile -id` and invalid
   option values. (USAGE.txt has no `-partial` option: partial ORFs are covered by the `partial*` designs, `-s 2`, `-c` and the region ends.)
4. `in_container.sh` runs them in parallel (all CPUs) with a per-run timeout and stores, per run, `cmd.txt`, `stdout.txt`,
   `stderr.txt`, `rc.txt`, `meta.tsv` and the `-out` file `orfs.txt` under `out/runs/<job id>/`; plus `out/ENV.txt` (date, `uname -m`,
   `/etc/os-release`, sha256 of every file, `ldd`, `ORFfinder -version`, `-version-full`), `out/USAGE.txt`, `out/CHANGELOG.txt`,
   `out/runs.tsv` (id, group, exit code, milliseconds, byte counts), `out/jobs.tsv`, `out/SHA256SUMS`.

## The commands (Docker Desktop must be running; on Apple Silicon the x86-64 binary runs through emulation, which is why every
`docker` call uses `--platform linux/amd64`)

```bash
cd docs/Validation/review-2026-09/handoff/A1-18-orffinder
QUICK=1 ./run.sh          # optional smoke test (~600 runs, a few minutes); look at out/ENV.txt and out/rc_histogram.txt
./run.sh                  # the full run (1750 runs; `JOBS=8 ./run.sh` to change the parallelism); builds the image, runs, zips
```

`run.sh` ends by writing **`returned/A1-18-orffinder-out.zip`** (zip of `./out`; `tar.gz` if `zip` is missing). It is resumable
inside one container run only; re-running starts from scratch (`rm -rf out` is done for you).

## What to send back

Put the archive where the next session can read it, in the repository:
`docs/Validation/review-2026-09/handoff/A1-18-orffinder/returned/A1-18-orffinder-out.zip`
(`git add returned/*.zip`, commit, push to the campaign branch `claude/busy-sagan-8ztcyb`; or attach the file to the session).
If the archive exceeds ~25 MB, send the result of `QUICK=1 ./run.sh` plus `out/runs.tsv` + `out/ENV.txt` and say so.
Also tell the next session if the build stopped at the sha256 check (NCBI replaced the file): then rebuild with
`docker build --platform linux/amd64 --build-arg ORFFINDER_GZ_SHA256=<new> .` and send the new value along.

## What the next session does with it

`INTEGRATION.md` (same directory): the exact C# changes planned per option (`OrfScanner`, `GenomicAnalyzer.FindOpenReadingFrames`,
`OrfSearchOptions`, MCP `find_open_reading_frames`), the oracle TSV generator that turns `out/` into
`tests/.../TestData/Analysis/orffinder_cli_oracle.tsv`, and the tests that must pass. Nothing here changes any library code.

## Privacy

Nothing in the package contains an e-mail, name or account id; the only network use is the three NCBI FTP downloads inside the
Docker build (User-Agent `Seqeron-review/1.0`) and, in one `misc` run (`-id NC_001416.1`, 120 s timeout), whatever the program
itself contacts (that run is allowed to fail).

## Authoring-time validation (the sandbox has the Docker CLI but no Docker daemon, and refuses to run the NCBI binary)

* **Not done here:** `docker build`, any execution of ORFfinder. The Dockerfile has only been checked by reading.
* Done: the binary was downloaded **but never executed**, only inspected: `file` = ELF x86-64 dynamically linked; `readelf -d` NEEDED =
  libz, libbz2, libuv, libnghttp2, libresolv, libdl, libdw, libstdc++, libm, libgcc_s, libc (-> the Debian packages installed in
  the Dockerfile; GLIBC <= 2.14); sha256 of `.gz` / unpacked binary / USAGE.txt / CHANGELOG.txt recorded in the Dockerfile.
* `python3 -I make_inputs.py --check` (byte-identical regeneration), `python3 -I -m py_compile`, `bash -n run.sh in_container.sh`;
  `python3 make_jobs.py` = 1750 unique job ids; `in_container.sh` was exercised end to end with a **stub** `ORFfinder` script
  (615 QUICK runs, 4 parallel, all outputs/`runs.tsv`/`SHA256SUMS` produced) to validate the runner mechanics.
