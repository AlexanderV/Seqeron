# S1B hand-off: the production chain at 3x10^9 symbols and the arm64 tests on your Mac

**Status of this package.** Written and run on a linux/amd64 machine (4 vCPU, 15 GB) with a Docker daemon that was started by hand: the image built there, and the chain ran inside it at 10^7 symbols (three kinds of text). **Never built or run on linux/arm64, on Docker Desktop or on macOS.** The scripts also ran natively (10^6 to 10^8 symbols, the self-test of 65 checks, a build killed at 2.5 s and resumed). Section 9 lists what is verified and what is not. If something here does not behave as written, that is a finding, not your mistake: stop and send the console text.

## 1. What this answers

Two things the orchestrator could not measure (it has no Mac, 15 GB of RAM and about 9 GB of disk):

1. **The production chain at 3x10^9 symbols**: `SesaBuilder.Build` (the build that the product uses; the Bench command `sesa resume build`), then the checks on its container. The questions: is **RssAnon at most 320 MB at every size** with the default parameters and four threads (owner decision of 2026-11-18, A1 relaxed); how long does each phase take; how much disk does it need at its peak; is the container correct (the verifier `e8-verify --container`, eleven checks); do count / exists / locate answers equal an independent oracle.
2. **arm64**: the same container built on Apple silicon must be **byte-identical** to the one built on x86-64 (the output depends on the text and on PSI only); the CRC-32C **instruction** path, the window decode and the platform code run; and the test sets of the CI policy run on arm64.

Both are Linux questions: Docker Desktop runs a Linux VM, and inside it the program is a Linux arm64 process. Native macOS is not supported (the pre-flight refuses it, section 5).

## 2. Your machine, from the E8 run of 2026-11-18

Its pre-flight (`analysis/e8/results/mac-20261118/preflight.last.txt`) showed: Apple silicon (aarch64), 18 cores visible, a Docker Desktop VM with **8.3 GB of memory**, 1.1 GB of swap, cgroup v2, a disk image of 1.96 TB with **1.86 TB free**, ext4, .NET 10.0.15. If that is still so, everything below fits (the 3x10^9 run needs about 86 GB, the VM's 8.3 GB is plenty). **Check again**: Docker Desktop settings change, and `preflight` says what it sees.

Docker Desktop, Settings, Resources: **Memory 8 GB or more** (16 GB if you want to run the test sets next to something else), **CPUs 4 or more**, **Swap at the minimum** (with swap on, RssAnon can read lower than what the build asked for: a run with swap in use proves nothing about the gate), **Disk image size** with at least 100 GB free. Nothing to set for file sharing: the data lives in a **named Docker volume** inside the VM, never in a Mac folder (a bind mount goes through virtiofs / FUSE, which is slow and unreliable for mapped files; the pre-flight refuses it).

**Apple silicon or Intel.** The wrapper asks the Docker daemon for its architecture and builds and runs the image for **that** architecture only. On an Intel Mac that is linux/amd64 (native, but then the arm64 question stays open); an arm64 image on an Intel Mac would run under QEMU emulation, and its timings mean nothing: the wrapper and the pre-flight refuse it.

**Keep the Mac awake**: `caffeinate -i <command>` in front of every long command, or run detached (`--detach`, then `docker logs -f seqeron-s1b-run`). A sleeping Mac suspends the VM; run the same command again afterwards (section 4).

## 3. Step by step

`W=docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff/s1b-docker.sh`, run from the repository checkout. **Times are ESTIMATES unless marked measured** (basis in section 6).

| step | command | what it does | expected |
|---|---|---|---|
| 1 | `$W image` | builds the image for your Docker architecture: the chain programs, the scripts, and (default `full` target) the .NET SDK, the sources and warm builds of the three test projects. `--chain-only` skips the SDK and the tests (0.4 GB). | 10-20 min the first time (downloads the base images and the NuGet packages of the test projects) |
| 2 | `$W preflight --n 3g` | what the machine is; refuses with numbers (section 5) | 1 min |
| 3 | `$W chain --n 100m` | **the arm64 chain at 10^8**: text, build, hash, verify, Structure level, queries, form, checker. The checker compares the container's SHA-256 with the x86-64 reference. | 3-8 min |
| 4 | `caffeinate -i $W chain --n 3g` | **the production chain at 3x10^9** (genome-like text, Resume on): see 3.1 | 25-60 min |
| 5 | `caffeinate -i $W chain --n 2g --structure yes --tag reader` | the same at 2x10^9, **where the reader works**: Structure level and the engine's answers (3.2) | 40-120 min |
| 6 | `caffeinate -i $W tests all` | the six test sets of the CI policy (section 7) | 1-2.5 h |
| 7 | `$W check` | the checker again over everything that was done | seconds |
| 8 | `$W export <a folder on the Mac>` | copies `form.json`, `form.md` and the small result files into the folder and makes `s1b-return.tar.gz`. **Send `form.json` back by hand**; nothing in the package sends anything anywhere. | seconds |

Run step 3 first: it costs minutes and tells whether the image, the volume and the checker work, **before** step 4 spends an hour. Step 4 and step 5 can be repeated with `--tag b`: a second container of the same text must have the same SHA-256 (the checker says so).

### 3.1 The production chain at 3x10^9 (step 4)

`chain --n 3g` does, in order: **gen** (the genome-like text of E3-E8 with `e8-gen`, 7 ns per symbol measured, a 3 GB file: the build reads its input twice, so it needs a file and cannot be fed from a stream) → **build** (`sesa resume build`, four threads, PSI on, fsync on, **Resume on**) under the independent sampler → **hash** (SHA-256 of the 40 GB container) → **verify** (`e8-verify --container`: eleven checks; no in-RAM reference; reads the sections in place) → **Structure level** (see below) → **queries** (3.3) → **form** and **checker**. Every 30 s it prints a line with the last phase marker, the RssAnon now and the size of the work directory. Where it ends is written in the data volume; the pre-flight and every step print what they refuse or warn about.

* **Which text.** `--kind genome` (default): the genome-like family of E3-E8 (random background, interspersed repeats with 10-15 % divergence, satellite arrays, N runs). It is the harder of the measured texts for the 320 MB gate (296 MB at 10^8 in 1b.9d, against 262 MB for random DNA) and the one whose generator is tested byte for byte. `--kind dna` (random ACGT: the stream of the benchmark, `dna_10m` and `dna_100m` are its prefixes), `--kind rep` (a 100 kbp block copied with 1 % substitutions) and `--kind prot` (protein-like, 20 letters) are switches for a second, smaller run (`--n 300m --kind prot`).
* **Resume.** `--resume auto` (default) = **on from 10^9 symbols**, off below. **Use it for the 3x10^9 run**: a kill after 40 minutes (the Mac sleeps, Docker restarts, you press Ctrl-C) then costs the stage in flight, not the run. What it costs: one `fsync` and one CRC-32C pass over the output of each of six stages and a state file (MEASURED at 10^7: 0.07 s of 3.0 s; at 10^8: see the ADR table; ESTIMATE at 3x10^9: 1-3 minutes). The output is identical with Resume on or off (a test and the kill rehearsal below). Below 10^9 the default is off, as in the build itself (no state file, no lock, no extra `fsync`).
* **Options of the build** (reasons in the ADR): `SesaBuilder.Build` with the default parameters of the product except `Resume`; threads 4 (the gate is at four); `ReservationMode` automatic (the container is reserved on disk up front with `fallocate`, so a full disk is an `IOException`, never SIGBUS); fsync on (the atomic publication of the file); `ScaleN` 0 (the real geometry, not the scale-down emulation of E3-E7); no `DOTNET_GCHeapHardLimit` (the gate is the RssAnon of the whole process as the product runs it).
* **The Structure level.** The reader of the container (`SuffixIndexFactory.Open` at the level `Structure`: the whole structure re-derived in one pass) **takes at most 2,147,483,645 symbols** (a 32-bit surface, stage S2). **At 3x10^9 it refuses**, and `chain` records that as `not applicable` after the checksum pass over every section (a real call: it also proves the sections' CRC-32C at that size). The checker prints `N/A S1`, not a pass. At 3x10^9 the check of the container's content is `e8-verify --container` (which does not use the reader). So the Structure level and the engine are run at **2x10^9 (step 5)**, the largest round size under the limit. Cost of the level: MEASURED 21 s at 10^8; ESTIMATE 20-60 min at 2x10^9 (an O(n) pass that scans the suffix array once per 4 Mi positions: 477 scans of 8 GB; it is disk-bound when the suffix array does not fit the page cache). `--structure auto` runs it up to 3x10^8 symbols and above the limit; between, it is **off** and you start it on its own: `$W structure --n 2g --kind genome --tag reader`.
* **Disk** (MEASURED at 10^8, ESTIMATE above): the container (13.4-13.7 bytes per symbol, reserved at the start) + the scratch (at most 4.2 bytes per symbol during PH1-PH3, **11.9 at the start of PH4**, when the suffix array is complete, falling to 0) + the text (1 byte per symbol). **3x10^9: about 3 GB + 40.5 GB + 36 GB = 80 GB at the peak** (the pre-flight asks for 86 GB). The text and the container stay after the run (`clean --n 3g` removes them; the rows and hashes stay).

### 3.2 The reader run at 2x10^9 (step 5)

Exactly the same chain at `--n 2g`: the Structure level (a bitmap of n/8 bytes: about 0.3 GB of RssAnon, not the build's gate), and the queries answered by the production engine as well as by the search over the suffix array. About 2/3 of the time and disk of step 4 plus the Structure level.

### 3.3 The query set and its oracle

1,500 queries, fixed by a seed (20261122) and the text: for each of the lengths 12, 20 and 32, 300 windows **cut from the text** at random positions and 200 random strings over the alphabet of the text (almost all absent). They are answered **(a)** by a binary search over the container's suffix array and the input text (two reads per step; works at any size the builder makes, above 2^31 as well) and **(b)** by the production engine (`SuffixIndexFactory.Open`: `Contains`, `CountOccurrences`, `FindAllOccurrences`; only where the reader works). The **oracle** is one pass over the input text with a rolling hash for each of the three lengths (about 8 ns per symbol for all three, so about 25 s at 3x10^9): every window whose hash matches a query is compared byte by byte, so it counts **exactly** every occurrence of every query and records the first 1,000 positions. **What it proves:** that on these 1,500 queries the container's suffix array is complete and ordered as far as a search can see, and the engine's count, exists and locate answers equal the text's own. **What it does not prove:** anything about other queries or other lengths; damage that no query reaches (one wrong suffix-array entry is found by `e8-verify`, not by a query set: 300 random queries touch a small share of the ranks); positions beyond the first 1,000 of a query with more occurrences (the count is exact, the positions are checked for what was reported); that the container is the container of this text beyond the length check (`e8-verify --input` compares the text symbol by symbol and its SHA-256). A hash cannot cause a false answer: every candidate is verified against the bytes.

## 4. Stopping, resuming, and what to do when it fails

* **Stop:** Ctrl-C in the foreground, or `docker stop seqeron-s1b-run` when detached. The script kills the build and tells you to run the same command again.
* **Resume:** run the **same command** again (`$W chain --n 3g`). With Resume on the build finds `state.json`, verifies the CRC-32C of what it keeps, **zeroes every byte that no kept stage vouches for** (1b.9e.1) and continues with the first stage that is not done (stages: packed text, rank file, suffix array, LCP and child table, BWT, PSI). `gen` and the hash reuse what is there. With Resume off a killed build starts again from the beginning (its scratch is removed first). After a Mac sleep, a Docker restart or a crash of the VM: the same.
* **Rehearse it first** (10^8 symbols, 3 minutes): `$W build --n 100m --resume on --kill-after 9 --tag kr` (SIGKILL after 9 s; exit code 4, the message says what to do), then `$W build --n 100m --resume on --tag kr`, then `$W hash --n 100m --tag kr`: the container must equal the one of the uninterrupted run (the checker compares the two when both are in the form). MEASURED here on x86-64: see the ADR.
* **What each failure means**
  * `REFUSED` by the pre-flight: the line says what and how much (disk, memory, a bind mount, emulation, a missing program). Fix it; nothing was started.
  * Exit code **137** from the build: SIGKILL. In Docker Desktop that is almost always the VM's out-of-memory killer: raise the VM memory (Settings, Resources) and run the same command again. (`docker stop` and `--kill-after` give the same code.)
  * Exit code **3** from the build: an `IOException`, nearly always a full disk (the message is in `runs/<label>/build.N.err`). Nothing is corrupted and the container is not truncated: free space, run the same command again (with Resume on it continues).
  * **R1 or R2 FAIL** in the checker: RssAnon above 320 MB. That is the finding, not an error of yours. Send `form.json` **and** `runs/<label>/build.N.series.tsv` (the RssAnon every 100 ms and the phase markers in `build.N.out`): they show which phase. The container may still be correct (the verifier does not look at memory).
  * **I1 FAIL** (container hash differs from the reference or from your other run): the most serious result, a nondeterminism or an architecture dependence of the builder. Keep both containers and send the form.
  * **V1 FAIL**: a check of `e8-verify` failed; `runs/<label>/verify.out` names it.
  * **Q1 FAIL**: the first mismatches are in `runs/<label>/queries.out` (rows `q_mismatch`).
  * The Docker daemon is not reachable / the image is for another architecture: the wrapper says so and exits with 2.
* **Disk:** `docker system df` shows what Docker holds. `$W clean --n 3g` removes a run's text, container and scratch; `$W clean-volume` deletes the whole volume (asks first).

## 5. What the pre-flight refuses or warns about

`$W preflight --n 3g` prints `[ ok ]`, `[warn]` or `[REFUSE]` for: the architecture and **emulation** (kernel against CPU format and image architecture, an emulator mapped into the process); cores (warn below 4: the build runs four threads regardless); **memory** visible to the container and its cgroup limit (warn below 4 GB, refuse below 1 GB) and **swap** (warn); **disk**: the free space of the data volume against the need at the peak for 10^8, 10^9, 3x10^9 and your `--n` (refuses when the run does not fit; warns when it fits barely); the **file system** of the data directory (refuses fuse, virtiofs, 9p, nfs, osxfs, tmpfs, which is what a Mac bind mount looks like; warns on overlay: the data directory is then in the container's layer, not in a volume); a 64 MB write, `mmap`, `msync` and `fallocate` test with its speed; .NET and the program; the **CRC-32C path** (the instruction, checked against the table path; a wrong value refuses); whether `RssAnon` is readable. **On macOS natively it refuses**: the gate is a Linux `RssAnon` gate (`/proc/<pid>/status`, anonymous resident pages without the mapped file pages of the container); macOS has only `ps` and `footprint`, which count other things, and the build's disk reservation (`fallocate`) is Linux code. Inside Docker everything is Linux. If you do run the sampler where `/proc` has no RssAnon, it asks `ps` for the RSS, writes `source: ps-rss`, and the checker prints `N/A R2` (not judged).

## 6. Estimates and their basis

All times are for the build with four threads, PSI on, fsync on. **MEASURED** = on the orchestrator's 4-vCPU x86-64 Linux box (virtio disk, warm cache), by this package's scripts. **ESTIMATE** = scaled, not measured. The measured rows and the scaling are in ADR-004, section 5a, "1b.9f".

| | 10^8 | 10^9 | 3x10^9 | basis |
|---|---|---|---|---|
| disk at the peak, genome-like | 4.7 GB asked, 2.5 GB used + text 0.1 | 30 GB asked, 27 GB used | **86 GB asked, 80 GB used** | MEASURED at 10^8 (container 13.49 B/char, scratch peak 11.9, text 1), linear above |
| RssAnon (the gate: 320 MB) | 287 MB (genome), 251 MB (random DNA) | ESTIMATE 290-320 | **unknown: this run decides** | MEASURED at 10^8; the E8 prototype went 170 -> 226 -> 310 MB at 1, 2.2, 3 Gbp |
| build, this box | 17-25 s | ESTIMATE 3-6 min | ESTIMATE 9-25 min | MEASURED 169-249 ns per symbol at 10^8; x30 linear, up to x1.5 for what has never run at this geometry |
| build, your Mac | ESTIMATE 10-25 s | ESTIMATE 2-6 min | ESTIMATE **5-25 min** | the E8 prototype ran 123-137 ns per symbol on that Mac against 205-300 on this box (0.5-0.6x) |
| `e8-verify --container` | 5 s | ESTIMATE 1-3 min | ESTIMATE 3-8 min | MEASURED 4-6 s at 10^8 (this box, cached); 125 s at 3 Gbp for the prototype's files on your Mac |
| SHA-256 of the container | 1 s | 10-30 s | 30-90 s | disk speed |
| queries (set of 1,500) | under 1 s | ESTIMATE 10 s | ESTIMATE 1-2 min | MEASURED oracle 0.7 s at 10^8; the search is about 200,000 small reads |
| Structure level (reader) | 21 s | ESTIMATE 3-8 min | **not applicable (n > 2,147,483,645)**; at 2x10^9 ESTIMATE 20-60 min | MEASURED 19-22 s at 10^8 |
| text generation | 1 s | 10-30 s | 1-2 min | MEASURED 7 ns per symbol (`e8-gen`), 19 ns (random DNA, Python) |

The chain at 3x10^9 (without the Structure level) is therefore ESTIMATE **25-60 minutes** on your Mac, dominated by the build; with step 5 about two to three hours of machine time in all, before the tests.

**Not measured, and why it may differ**: the external levels of the difference-cover recursion of PH1 run only above about 5x10^8 symbols at the real geometry (the reduced string has 3.2 M symbols at 10^8 and fits the base case of 16 Mi; at 3x10^9 it has about 100 M): they were run only in the scaled geometry and in tests. The page cache (8.3 GB for the VM against a text of 1.1 GB packed, a rank file and 12 GB of suffix array) is smaller than the working set at 3x10^9: the E8 sweep and 1b.9d say the build slows when the cache is below text + rank file (0.5 GB at 3x10^9), not when it is below the scratch.

## 7. The test sets (step 6)

`$W tests all` runs, in the image, with `--network none`, the sets of the CI policy: **esa-default** (`SuffixTree.Esa.Tests`, Debug), **esa-exhaustive** and **esa-fuzz** (Release, `Category=Exhaustive` / `Fuzz`), **esa-checked** (`ESA_EXPECT_BOUNDS=true dotnet test -p:EsaBounds=true`: every pointer access bounds-checked), **emproto** (`SuffixTree.EmProto.Tests`, Release), **persistent** (`SuffixTree.Persistent.Tests`). Each writes `tests/<set>.log` and `tests/<set>.json` (passed, failed, skipped, seconds); the checker reads the json. Test count of the reference head (x86-64): esa-default 4,846 + 19 skipped, esa-exhaustive 82, esa-fuzz 39, esa-checked 4,847 + 18 skipped, emproto 535, persistent 672. **A different count is a warning, a failed test is a failure**: on arm64 a platform may skip tests that need AVX2 or add none.

| set | the orchestrator's box (4 vCPU x86-64, bare Linux) | ESTIMATE on your Mac under Docker | realistic? |
|---|---|---|---|
| persistent | 0.5 min | 0.3-0.7 min | yes |
| esa-default | about 20 min | 12-25 min | yes |
| esa-exhaustive | about 17 min | 10-20 min | yes |
| esa-fuzz | about 15 min | 9-18 min | yes |
| esa-checked | about 20 min | 12-25 min | yes |
| emproto | about 16 min | 10-20 min | yes |
| **all six** | about 1 h 30 | **1-2 h** | yes: run it over a meal, `caffeinate -i`; or one set at a time |

Basis: the box's times are the orchestrator's measurements at head `522f8c5f`; the Mac factor 0.6-1.2 is the E8 observation (0.5-0.6x per symbol for the prototype, cores 18 against 4, tests mostly single-threaded) widened by the Docker overlay file system for the tests that write temporary files. The six sets need under 2 GB of memory each as far as the box showed; if one is killed with exit code 137, raise the VM memory. The image has all packages and warm builds: no network at run time.

## 8. The form and the checker

`form.json` (and `form.md`, the same numbers for a human) is **filled by the run**: after every step `chain` regenerates it from the files in the volume (the program's rows, the verifier's rows, the sampler's summary and series, the SHA-256 files, the test totals). `$W check` prints, per gate, **PASS / FAIL / WARN / N/A / NOT RUN** and ends with `RESULT: PASS`, `FAIL` or `INCOMPLETE` (exit codes 0, 1, 3; 2 = not a readable form). The gates are fixed **in the checker**, not read from the form:

| gate | what | PASS when |
|---|---|---|
| R1 | RssAnon of the build, its own sampler (10 ms) | every attempt <= 320 MB |
| R2 | RssAnon, the independent sampler of `/proc/<pid>/status` (100 ms) | every attempt <= 320 MB (a figure from `ps` is not judged) |
| B1 | the build finished; container bytes equal the build's own `file_size` | exit 0, sizes equal |
| I1 | container SHA-256 | equals the reference of `reference_hashes.tsv` for the same text, size and PSI (built and verified on x86-64) and every other run of the same text in the form; **N/A** when there is no reference (a new size) |
| I2 | the input text's SHA-256 | equals the reference's |
| V1 | `e8-verify --container` | exit 0, all eleven checks `pass` (a skipped check is a FAIL) |
| S1 | the Structure level of the reader | `pass`; `N/A` above 2,147,483,645 symbols; `NOT RUN` is allowed (optional) |
| Q1 | the query set | 0 mismatches of the search and of the engine; at least 500 queries |
| D1 | peak scratch | at most 13 bytes per symbol (a WARN, not a FAIL) |
| P1 / X1 | no emulation; CRC-32C check value and table = hardware | PASS (X1 WARN when the CPU has no instruction) |
| T1 | the test sets | no failed test (a count that differs from the reference head's is a WARN) |

Every one of these faults was tried in `test_s1b_check.py` (30 tests): a doctored RssAnon (own sampler, independent sampler, a killed attempt), a hash mismatch against the reference and against the other run, a missing section or field, a verifier check that is not a pass, a query mismatch, a build that did not finish, an emulated process, a failed test set.

## 9. What is verified here, and what is not

**Verified on the orchestrator's machine (linux/amd64, 4 vCPU, 15 GB):** the context (648 files, 6.3 MB) builds both programs offline with an empty package cache and the three test projects from a pristine copy of the context restoring from a local feed; the image (`chain` target) built with a Docker 29 daemon started by hand, and the chain ran inside it at 10^7 on the three kinds of text (PASS in the checker, container hashes equal to the native build's); the scripts natively at 10^6, 10^7, 2x10^7 and 10^8 on four kinds of text (self-test of 65 checks: a second run gives the same container, doctored forms, a missing program, too little disk, a bind-mount-like file system, a build killed and resumed with Resume on and off); the unit tests of the checker (30), the sampler (6), the generators (8: byte-identical to the benchmark's in-memory ones, `dna_10m` reproduced from the file of the benchmark) and the query pieces (20 tests in `SuffixTree.Esa.Tests`).

**NOT verified:** the image on **linux/arm64**, on **Docker Desktop**, on **macOS**; cgroup v2 as Docker Desktop presents it (the cgroup readers of E8 ran against synthetic trees); `drop_caches`; **anything at 3x10^9** (the largest build ever run is 10^8 in this package and 5x10^8 for the E8 prototype); the arm64 CRC-32C instruction (tested through a software stand-in of its steps, 1b.7e; it never ran on hardware); the `full` image's test sets inside a container beyond the Persistent set; the behaviour of the reader and the engine on a container above 2^31 symbols (the code refuses by design; no such container was built here); an `fsync` on a real device and a power loss.

## 10. Files of the package (`docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff/`)

`Dockerfile`, `Dockerfile.dockerignore` (image; allow-list of 648 files), `s1b-docker.sh` (host wrapper: image, preflight, plan, chain, gen, build, hash, verify, structure, queries, tests, form, check, status, clean, export, shell, clean-volume; `--dry-run`), and inside the image `s1b.sh` (entry point), `preflight.sh`, `s1blib.sh` (helpers; the cgroup, memory and file-system helpers are those of `analysis/e8/e8lib.sh`), `s1b_sampler.py` (independent RssAnon and disk sampler), `s1b_gen.py` (streaming generators of the random DNA, repetitive and protein-like texts), `s1b_form.py` (fills the form), `s1b_check.py` (the checker), `reference_hashes.tsv` (the reference containers), `s1b_refs.py` (maintainer: collects verified containers into the table); tests: `selftest_s1b.sh`, `test_s1b_check.py`, `test_s1b_sampler.py`, `test_s1b_gen.py`; the programs: `apps/SuffixTree.Esa.Bench` (`sesa resume build`, `sesa handoff env | open | queries`) and `apps/SuffixTree.EmProto` (`e8-gen`, `e8-verify`); `analysis/e8/check_docker_context.py` checks the context without a daemon (`--package-dir`, `--publish`, `--allow-tests`, `--tests-source`).

**Privacy.** Every container runs with `--network none`; the scripts open no connection and write nothing outside the volume. The image build downloads public base images, Ubuntu packages and (target `full`) NuGet packages; it uploads nothing, and the labels carry only a title, a description and the git hash. `form.json` holds architecture, kernel release, core count, memory, file-system type, .NET version, the figures and the hashes: no host name, user name, folder or address. The sandbox of the orchestrator is not named anywhere in the package.
