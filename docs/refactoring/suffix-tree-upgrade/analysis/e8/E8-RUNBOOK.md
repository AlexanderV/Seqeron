# E8 runbook: the B-DCS pipeline at 1, 2.2 and 3 Gbp on your Mac

**Status of this package: the image is UNBUILT and UNTESTED.** It was written in a session that had no Docker daemon. Nothing in it was ever run inside a container, on macOS, or on arm64. What was run, on a Linux x86-64 machine
(4 vCPU, cgroup v1, no Docker): the program and the verifier up to 200 Mbp, the shell and Python scripts at that scale, and the build of the program from exactly the files the image build would see (offline, empty package cache). Section 9 lists
what is verified and what is not. If something in this runbook does not behave as written, that is a finding, not your mistake: stop and send the console text.

## 1. What E8 is, what you run, what it answers

E8 is the whole B-DCS construction (command `e6`: PH0 prepare, PH1 difference-cover ranks, PH2/PH3 suffix array, PH4 LCP fused with PH5 child table; settings as in E6: layout A, codec Packed, PH4 base mode, no chain encoding,
4 threads, 96 MiB GC limit, the scale rule of E6 / E7 with `--scale-n 3000000000`, i.e. at 3 Gbp the reference block sizes) on a deterministic genome-like text of **1, 2.2 and 3 Gbp**. The largest run so far is 500 Mbp. Nobody has run
the 2.147 Gbp boundary (2^31): positions above it are unverified until your run. The protocol, cheap first:

| step | what | ESTIMATE on the prototype machine (4 vCPU x86, virtio disk; E7: 205-300 ns per symbol) | stop rule |
|---|---|---|---|
| a | 1 Gbp, generous memory, then the full verification | build about 4-5 min, verification 1-3 min | none: the baseline of step d |
| b | 2.2 Gbp, generous (crosses 2^31) | about 9-11 min + 2-6 min | none |
| c | 3 Gbp, generous | about 11-15 min + 3-10 min | none |
| e | 1 Gbp with 1 and 2 threads (the outputs must be identical to step a's) | about 8 + 5 min | none |
| d | the limited-cache sweep: the same 1 Gbp text at decreasing container memory | 1x to 6x the time of step a per cell, 8 cells | a cell is stopped after **6 times** the generous time ("does not finish in practice at this memory"); the sweep stops after two such cells in a row |

**The question E8 answers:** how much page cache PH4's 12 B/char of scratch needs beyond the 1x (text + rank file) that ADR-004 6b fixed as the floor (E7: the whole pipeline in one constant limit did not finish below an effective 0.7 of text + rank
at 100 Mbp). Everything else (times, I/O, disk) is reported on the way.

## 2. What the numbers mean, and what they do not

* Docker Desktop is a **Linux virtual machine on a virtual disk file**. The file may be cached by macOS (a third layer of caching you cannot see or drop), and Docker's `--memory` limits one container inside the VM. Your numbers are "this Mac under
  Docker Desktop", not "the algorithm". They are **not comparable** with the prototype machine (x86, virtio disk, 15 GB RAM), and Apple-silicon (arm64) timings do not transfer to x86 or the other way round.
* On arm64 the suffix comparison uses the scalar 64-bit path instead of the AVX2 vector path (the results are identical, a unit test pins both); it is slower on tie-rich text. This is expected and shows up in `env_lce_path` of the preflight.
* What does carry over: whether the outputs are **correct** (the verifier), the **memory** of the process (RssAnon, G2), the **disk** footprint (G3), the shape of the cache table (where the knee is relative to text + rank).
* Page cache dropping needs a privileged container. Without `--cold` every run says `cold_start=no` ("warm start") in its rows.

## 3. Before you start

* **A Mac with Docker Desktop** (chip, RAM and disk unknown to the author). Plan for: at least 8 GB of RAM for the Docker VM (16 GB if the Mac has 32 GB), 4 CPUs, and **100 GB or more of free disk for the Docker disk image**
  (the 3 Gbp run needs about 70 GB at its peak: ASCII text 3 GB + scratch 36 GB at the peak of PH4 + outputs 28.5 GB, ESTIMATES; the 1 and 2.2 Gbp texts stay on the volume, 3.2 GB more).
* **Docker Desktop settings** (Settings, Resources): Memory 8 GB or more (the sweep limits single containers below that, so the VM must be larger than the largest limit; cells above the VM's memory are skipped); CPUs 4 or more;
  Swap at the minimum; **Disk image size at least 100 GB** (the default is smaller; the preflight prints the actual size of the file system). Nothing to set for file sharing: the data lives in a **named Docker volume** inside the VM,
  never in a macOS folder (a bind mount goes through virtiofs / gRPC FUSE, which is slow and unreliable for the memory-mapped files of this program).
* **Apple Silicon or Intel:** the wrapper asks the Docker daemon for its architecture and builds and runs the image for **that** architecture only. Never run an amd64 image on Apple Silicon: it would run under emulation (Rosetta / QEMU) and the
  timings would be meaningless. The wrapper refuses a platform that differs from the daemon's, and the preflight inside the container refuses when it detects emulation (kernel versus CPU format, an emulator mapped into the process,
  the .NET process architecture).
* **Keep the Mac awake** for the long steps: `caffeinate -i <command>` (macOS), or run detached (`--detach`, below). A sleeping Mac suspends the Docker VM; the run in progress is then purged and repeated by the next call.
* **Nothing about you goes into the image.** The build context is your local checkout, filtered by `Dockerfile.dockerignore` (an allow-list: the program, the two libraries it needs, the scripts; no `.git`, no tests, no wiki, no `bin`/`obj`).
  Labels carry only the title, a description and the git hash you build from. No credentials, no tokens, no `ARG` for either. The build pulls the two public .NET base images from mcr.microsoft.com and a few Debian packages
  (bash, python3, procps, coreutils, bc, util-linux, findutils) with `apt-get`; it fetches **no NuGet package** (the Dockerfile replaces the repository's root `Directory.Build.props`, which adds an analyzer package, by a file in
  `docker/` without it) and no source from the network.

## 4. Step by step

All commands are run from anywhere; the wrapper finds the repository root. `W=docs/refactoring/suffix-tree-upgrade/analysis/e8/e8-docker.sh` below.

1. **Build** (once; 3-8 min the first time, the base images are 0.7 GB): `$W build`. It prints the architecture it builds for. If it fails, send the last 40 lines of the output.
2. **Preflight** (1 min): `$W preflight`. Read it top to bottom. `[ ok ]` is fine, `[warn]` is something to know, `[REFUSE]` means the run will not start until you change it; the end of the report says in plain words what to change in Docker Desktop.
   It prints: architecture and whether the process runs under emulation, cores, memory seen by the VM and by the container (cgroup v1 and v2 are both read), swap, free disk of the data volume against the need of the largest run
   (n x (1 + 12 + 9.5) B plus a margin; **refuses** when it does not fit), the file-system type of the data volume (**refuses** fuse, 9p, virtiofs, nfs, tmpfs), a 64 MB write + mmap + msync + fallocate test with its speed, whether `drop_caches`
   is possible, the .NET version and which vector path the program will use. Run it again after every Docker Desktop change.
3. **Plan**: `$W plan` prints the table of section 1 with the state of every step.
4. **Step a** (1 Gbp, about 10-15 min with the text generation and the verification): `caffeinate -i $W run a`. The first call generates the text (about 10-30 s plus the write) and checks the free disk.
   Console, per run: a header with the estimate and the abandon time; `verifying the outputs ...`; one summary line with `rc=0`, `killed=no`, `total_s`, `ext_rssanon_peak` (MB), `verify_all=pass`.
   **If `verify_all` is not `pass`, or `rc` is not 0: stop and send the rows file (section 6). Do not continue to the big steps.**
5. **Step b** (2.2 Gbp, 15-20 min): `caffeinate -i $W run b`. **Step c** (3 Gbp, 15-30 min): `caffeinate -i $W run c`. These are the runs the gates G1-G3 are judged on. Generous means: no memory limit on the container (the VM's memory is the limit).
6. **Step e** (identity across thread counts, about 15 min): `$W run e` (1 thread, then 2; hashes only).
7. **Step d, the sweep** (the main question; 1 to 3 hours, up to about 5 if several cells are abandoned): `caffeinate -i $W sweep`. It runs step a first if it is not done, then one container per cell with
   `--memory M --memory-swap M` (no swap) for M = 16, 6, 4, 3, 2.2, 1.6, 1.25 and 1.0 times (text + rank file) + 160 MiB, largest first (at 1 Gbp: 7792, 3024, 2080, 1600, 1216, 928, 768, 640 MiB; cells larger than the VM are skipped).
   Each cell prints its effective cache fraction, its time and its slowdown against step a when it ends. Cells that need more than 6 times the generous time are stopped; after two in a row the sweep stops. Use `--cold` (a privileged container, so that
   every run can drop the page cache) when you want cold starts: `$W --cold sweep`. To run chosen cells only: `$W sweep --cells "3024 1216 768"`. Replicates of a cell: `$W --memory 1216 run d --size 1g --memory-mb 1216 --runs 3`.
   The 3 Gbp text can be swept too if there is time (`--size 3g`, cells from `e8.sh cells 3000000000`); it is optional.
8. **Report and return** (section 6): `$W report`, then `$W export <folder>`.

Run in the background (survives a closed terminal, not a sleeping Mac): `$W --detach run c`, then `docker logs -f seqeron-e8-run`.

## 5. Reading the rows and the tables

`$W report` prints the markdown tables of RESULTS-E8 and judges the gates (also into the volume: `rows/RESULTS-E8.tables.md`, `rows/gates.json`). The raw rows are tab-separated:
`R  E8  <case>  <run>  <metric>  <value>  <unit>  <notes>`. Cases: `n1000m.gen` (generous), `n1000m.gen.t1` / `.t2` (threads), `n1000m.m1216` (limited, MiB), `n2200m.gen`, `n3000m.gen`. Metrics from the program (`total_s`, `ph*_s`,
`rssanon_peak`, `io_*`, `scratch_peak`, `output_bytes`) and from the **independent sampler** (`ext_*`: RssAnon every 20 ms; `/proc/<pid>/io`; `eff_resident_ph3` = the time-weighted fraction of the packed text and the DCS rank file that was
resident during PH3 (mincore); `cg_*` = cgroup v1 or v2 counters: usage, anonymous, file, mapped, dirty, writeback, limit hits (`failcnt` / events max) and OOM kills; `ext_disk_*` = scratch + outputs + input and the free space of the volume;
`marker_*` = when the scratch sub-directories ph1, ranks, ph3, ph4, ph5 appeared; `psi_memory_*` = memory stall time of the VM), and from the verifier (`verify_<check>`, `verify_fact_*`, `sha256_*`).

## 6. Stopping, resuming, and what to send back

* **Stop:** Ctrl-C (or `docker stop seqeron-e8-run`). The run in progress leaves a partial scratch directory; **it does not matter**: every call purges the unfinished run (its rows are removed because they lack `run_done`) and repeats it.
  Finished runs are skipped. After a Mac sleep, a Docker restart or a crash: run the same command again.
* **Free disk:** the runner removes the outputs of every run after it is verified or hashed (`--keep` keeps them) and always removes the scratch. The texts stay (`$W shell`, then `rm /data/text/genome_n3000m.txt`; they are regenerated by `gen`
  and checked by their SHA-256). `$W clean-volume` deletes the whole volume (asks first).
* **Send back:** `$W report` (prints the tables and writes `return/e8-return-<time>.tar.gz` in the volume), then `$W export <a folder on the Mac>` copies `rows/`, `sums/`, `logs/` and `return/` out. Send **the archive (it contains the rows file,
  the verifier output of every generous run, the hashes, the console logs of every e6 run, the sampler series and the preflight output)** and the **text of the preflight report** from your terminal. The rows file is small (tens of KB). Nothing in
  it names a host, a user, a folder or a network: the notes carry the architecture, the kernel release, the cgroup version, the limits and the commit hash of the checkout.

## 7. The gates (fixed by the orchestrator before the run; verbatim; `tables_e8.py` judges them as PASS / FAIL / NOT RUN)

* **G1** the verifier passes at 1, 2.2 and 3 Gbp and the output SHA-256 are identical between the generous and every limited-cache run of the same text and between thread counts.
* **G2** independently sampled RssAnon of the dotnet process at most 256 MB at every size (E7: 124 MB at 500 Mbp).
* **G3** peak scratch at most 12.5 B/char and peak total disk at most 70 GB at 3 Gbp.
* **G4** reported, not gated: time per phase and total (ESTIMATE from E7 for 3 Gbp: 11-15 min on the prototype machine, linear), I/O in B/char, and the cache table of (d) with the answer to the open question.
* **Abandon rule:** a run that takes more than 6 times the generous run of the same text is stopped and reported as "does not finish in practice at this memory".

"Peak total disk" is the peak of scratch + outputs + the input text of the run, sampled every second (apparent sizes). The verdict of the sweep is "the smallest container memory with a slowdown of at most 3x" (the E7 criterion), with the effective cache
fraction of each cell next to it.

## 8. What the verifier checks, and what it does not (`e8-verify`, run after every generous run)

Eleven checks, each reported as pass / FAIL / skipped (a skipped check is never a pass): `files` (all seven files exist with the lengths implied by n, bits, the LCP width and the overflow count), `meta` (the metadata is internally consistent),
`text` (symbol counts of the packed text equal the metadata; padding and unused bits zero), `input` (the packed text equals the ASCII input symbol by symbol; its SHA-256 equals the metadata's), `sa_permutation` (SA[0] = n; every entry at most n
and distinct: a bitmap in a temporary file, not in RAM), `sa_order` (every adjacent pair of suffixes strictly ordered), `lcp_values` (every stored LCP equals the directly computed common prefix of its adjacent pair, the overflow list matches
the escapes), `lcp_width` (the base width is the one the ADR-002 rule gives), `stats` (NODE_COUNT, internal nodes, LRS depth and rank, max LCP recomputed from the LCP stream by a stack scan: the file must be byte-identical),
`child` (the child table cell by cell against the lcp-interval definition, next links, up and down fields, and no other cell set), `sha256` (of all seven files).
**The limits, stated plainly:** (1) the common prefix of an adjacent pair is computed directly **up to 65,536 symbols**; a pair with a longer stored LCP ("long pair": N runs, satellite arrays, a^n) is checked up to the cap, at its mismatch
position (both symbols, hence the order) and at 16 further offsets (the last 8 symbols of the claimed match and 8 spread over it): a stored value too large by up to 8 is always caught, by more it is missed with probability about 4^-8 per
fault, and a corruption in the middle of a long match is found only by chance; the number of long pairs is reported (`verify_fact_long_pairs`). (2) The child scan needs a stack of 16 bytes per open interval (an N run of 15 Mbp: 240 MB); above
64 Mi entries it reports `skipped`, never `pass`. (3) The input check needs the ASCII text (the runner passes it). (4) `name=` in the metadata is not checked. (5) A fault that keeps every checked relation true is not seen: the verifier is not a proof.
Its memory is small (24 MB RssAnon at 200 Mbp) and its time was 59 ns per symbol at 200 Mbp with everything cached (12 s), i.e. ESTIMATE 3 minutes at 3 Gbp cached, 5-15 minutes when the files are not cached.

## 9. What is verified here, and what is not

Verified on the author's machine (Linux x86-64, 4 vCPU, 15 GB RAM, cgroup **v1**): the generator (byte-equal to the in-RAM family of E3-E7; at 200 Mbp byte-identical to the E7 input); the verifier (77 unit tests, a mutation audit, 33 injected faults on real outputs at
100 and 200 Mbp, all detected; outputs at 100 and 200 Mbp pass and have the SHA-256 of the references that E6 / E7 made with the in-RAM kernels); the boundary tests with synthetic values above 2^31; the scripts at 20 and 200 Mbp (generous run, thread runs,
a cgroup-limited cell, the report); the cgroup readers against synthetic v1 and v2 trees (`selftest_e8.sh`); the program built from the files of the image build context, offline, with an empty package cache (`check_docker_context.py`), giving outputs
byte-identical to the normal build.
**NOT verified** (unverifiable there, or beyond the budget): that the Dockerfile builds (the apt step, BuildKit's reading of `Dockerfile.dockerignore`, the Debian release that the `10.0` tag resolves to); that anything runs inside the image, on Docker
Desktop, on macOS, on arm64 (the program has x86 intrinsics guarded by `IsSupported` and a tested scalar fallback; its P/Invoke list is empty: every native call is a .NET framework API); **cgroup v2 on a real kernel** (this machine has only v1; v2 was exercised
against synthetic files only); `drop_caches`, `--privileged` and the `--memory` limit as Docker Desktop's VM presents them; the cache cells at the real block sizes and at 1-3 Gbp; everything above 500 Mbp; the verifier above 200 Mbp; positions above 2^31
end to end (the boundary tests reach the packed text, the rank table, the geometry and the mapped writers with synthetic values only). Section 10 of the PROGRESS entry and the ADR subsection list the boundary-audit findings.

**Architecture neutrality of the program (static audit; the publish is framework-dependent IL, built once on the build platform and used for both architectures):** P/Invoke and `LibraryImport`: **none** in `apps/SuffixTree.EmProto`, in the two libraries it references
(`src/SuffixTree/Algorithms/SuffixTree.Esa` and `.Core`) and in the linked sources of the benchmark project; every native call is a .NET framework API (`NativeMemory`, `MemoryMappedFile`, `RandomAccess`, `File`), so there is no native library to be missing on arm64.
`StructLayout(Sequential)`: `Ph4/Spool.cs:9,18,26` (three scratch-record structs of `uint` / `ulong` fields), `Ph5/ChildScan.cs:91` (`IntervalEntry`, four `uint`), `Ph5/NsvScan.cs:23`; they are written to scratch files and read back by the same process on the same
machine, little endian on both architectures: safe. Output files are little endian by construction (`BitConverter` / `MemoryMarshal` on a little-endian machine; arm64 Linux and macOS are little endian; `env_little_endian` in the preflight). x86-only intrinsics:
`Text/PackedText.cs:60,154,210-212` (`Avx2`, behind `UseVectors = Avx2.IsSupported`; scalar 64-bit fallback, identical results, pinned by a test), `Ph4/PlcpScan.cs:33`, `Ph4/LcpBuild.cs:383` and `src/.../Esa/Memory/Buf.cs:152` (`Sse.Prefetch0` behind
`Sse.IsSupported`, a JIT-time constant: no prefetch on arm64). The cross-platform `Vector128/256` APIs (`LcpKernel.Direct.cs`, the static shuffle mask of `PackedText`) have software fallbacks. Consequence: arm64 is **correct but slower** on tie-rich text.
`dotnet publish` with `UseAppHost=false` needs no runtime-specific pack, so the build stage never downloads one.

## 10. Choices made without asking (change them if you disagree)

* The child table is verified by a **full scan of every cell** (a stack scan of the LCP stream that predicts each cell and compares it, then counts the non-zero cells), not on a sample of intervals: it is the stronger check, it covers small and large intervals alike,
  and it costs one pass over the LCP and the expected cells; its price is the stack (above).

* `--scale-n 3000000000` for all three sizes: at 3 Gbp it gives the reference block sizes, at 1 and 2.2 Gbp the E7 series-B ratios (n / K = 3000, P = 375 chunks of PH4). `E8_SCALE_N` overrides; `E8_EXTRA` appends arguments to the `e6` command.
* The memory of one run is the **container's** limit (the cell list above); the program's own settings are those of E6 / E7 (96 MiB GC limit, 4 threads or fewer if the VM has fewer cores).
* Limited cells are **hashed** (SHA-256 of the seven files, compared with step a's); only the generous runs get the full verifier (a cell with different hashes gets the verifier, too).
* The abandon time of a cell is 6 x the process wall time of run 1 of the generous run of the same text; a cell needs that run first (the runner refuses otherwise).
* One run per cell by default (E7 used 3 in its main cells). `--runs 3` where a cell is cheap.
* The text is the genome-like family of E3-E7 (seed 20240002; random background, interspersed repeats with 10-15 % divergence, 171-bp satellite arrays, six N runs of 0.85 % in all, the largest n/200 = 15 Mbp at 3 Gbp). **The family's own scale rules apply to
  n**: the satellite arrays and N runs grow with the text.
* The sampler's disk figures are sampled once a second; the program's own `scratch_peak` row is exact.

## 11. If something goes wrong

* `REFUSED: this cell needs a container started with --memory ...`: the cell was started without the wrapper; use `$W sweep` or `docker run --memory 1216m --memory-swap 1216m ...`.
* Exit code 137 / `killed=oom-or-signal`: the cgroup OOM killer; for a cell this is a result ("does not finish at this memory"), for a generous step it means the VM has too little memory: raise it.
* The disk fills up: `docker system df`, remove other images / volumes, or raise the disk image size; the preflight says what the next step needs. A full disk during a run leaves a partial scratch that the next call purges.
* The preflight says emulation: you are running an image for the wrong architecture. `$W build` again (it rebuilds for the daemon's architecture) and remove the old image (`docker image rm seqeron-e8`).
* A step looks stuck: the runner prints only at the start and end of a run. `docker exec` into the container and `ls /data/work/scratch.*` shows the phase (`ph1`, `ranks.bin`, `ph3`, `ph4`, `ph5`); `top` shows the dotnet process. Abandon a generous step after about 3 times the ESTIMATE and report it.

## 12. Files of the package (`docs/refactoring/suffix-tree-upgrade/analysis/e8/`)

`Dockerfile`, `Dockerfile.dockerignore`, `docker/Directory.Build.props` (image build); `e8-docker.sh` (host wrapper: build, preflight, plan, gen, run, sweep, report, export, shell, clean-volume; `--dry-run` prints the docker commands); inside the image `e8.sh`
(entry point: preflight, plan, gen, run, verify, report, clean, shell), `preflight.sh`, `run_e8.sh` (the runner), `e8lib.sh` (shared helpers), `sample_e8.py` (the independent sampler), `tables_e8.py` (rows to the tables and the gates), `inject_faults.py`
(fault injection on real outputs); `selftest_e8.sh` (cgroup readers), `check_docker_context.py` (the offline closure check); `mutation/` (the verifier audit). The program: `apps/SuffixTree.EmProto` commands `e8-gen`, `e8-verify`, `e8-env`.
