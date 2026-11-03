# Hand-off package: competitor benchmark scripts

Everything needed to re-measure the competitors on another machine and append rows to `bench.tsv`.
All harnesses print the same line protocol and are driven by the same Python driver (`lib/benchlib.py`), so
rows are comparable with the Seqeron and Go rows already in `bench.tsv`.

Status on the machine that produced `bench.tsv` (4 vCPU Intel Xeon under KVM, 15 GB RAM, Ubuntu 24.04, gcc 13.3, cmake 3.28, go 1.24.7, .NET 10.0.403, rustc 1.97):
**every script below was run there and its rows are already in `bench.tsv`** (the owner lifted the "do not build competitors" rule for open-source tools mid-task; nothing was left unrun except the items in section 3).
They are shipped so the owner can reproduce, re-run on different hardware, or extend the size range; rows from another machine must NOT be mixed into the same `bench.tsv` without a different `tool`/`version` tag (e.g. add the host name to `--version`).

Package: `handoff_package.tar.gz` (generator, driver, harness sources, all `competitors/run_*.sh`, hashes, provenance; ~120 KB, no data files, no third-party sources).
On the target machine: unpack to a directory `bench/`, regenerate the inputs (section 0), then `competitors/run_all_handoff.sh` (or the scripts one by one), then `python3 lib/assemble.py bench`.

## 0. Inputs (deterministic, byte-identical for every subject)

```
python3 gen/gen_inputs.py <bench-dir>            # ~5 s, writes data/*.txt, data/*.fa, queries/*.txt, MANIFEST.tsv/.sha256
python3 gen/gen_extra_100m.py <bench-dir>        # optional dna_100m (used only for the Seqeron stop-point probe)
(cd <bench-dir> && sha256sum -c MANIFEST.sha256) # must print OK for every file -> inputs identical to the ones measured here
```
`MANIFEST.sha256` shipped with this package lists the hashes of the files measured here. Generator uses only the Python 3 stdlib
(`random.Random(str)` seeding is SHA-512 based and stable across CPython 3.x).

Texts: `dna_1m dna_10m dna_20m dna_50m` (nested prefixes of one seeded ACGT stream), `rep_10m` (100 kbp block, 100 copies, 1 % substitutions each),
`prot_1m` (20-letter uniform). Per text: `present20` (100 000 20-mers sampled from the text), `absent20` (100 000 random 20-mers),
`locate8` (1 000 random 8-mers; `locate4` for protein), `memq100k` (100 kbp window of the text with 1 % substitutions).

## 1. Driver and row format

* `lib/benchlib.py measure --raw RAW --tool T --version V --input NAME [--repeats 3] [--drop-caches] [--phase P] [--op-prefix X] [--timeout S] -- CMD...`
  runs CMD `repeats` times, parses its `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes` lines, adds wall time, CPU and peak RSS (`wait4` ru_maxrss), appends raw rows to RAW.
* `lib/assemble.py <bench-dir>` (also writes explicit `NA` rows for operations a tool's API lacks) merges every `runs/raw*.tsv`, takes the median over repeats and writes `bench.tsv`
  (columns `tool version input n operation metric value unit repeats notes`). To merge into an existing `bench.tsv`, copy your `raw_competitors.tsv` into
  `<bench-dir>/runs/` and re-run `assemble.py`; it is idempotent.
* `lib/crosscheck.py bench.tsv` verifies that every subject returns identical count / locate checksums (all tools measured here agree; a disagreement is a bug in a harness or a tool).
* Environment knobs for every script: `BENCH_DIR`, `WORK_DIR`, `RAW_FILE`, `REPEATS`, `INPUTS_DNA`, `TMO` (per-run timeout seconds).
* Run the scripts **one at a time and do not run anything else on the machine while they measure** (they time 1/4-thread throughput).
* Privacy: no script sends any identity anywhere; they only `git clone` public repositories, `pip install` public packages and `cargo` public crates.

## 2. Scripts (all under `competitors/`)

| script | tool(s) | what it does | needs | approx run time (4 vCPU) |
|---|---|---|---|---|
| `run_libsais_libdivsufsort.sh` | libsais v2.9.1 (Apache-2.0), libdivsufsort 2.0.1 (MIT) | clone @ pinned SHA, **cmake** build (libsais with OpenMP), link `sa_c/sa_bench.c`, SA build time/memory + binary-search count/locate, 1/4 threads, libsais_omp 2/4 thread build | git cmake gcc | ~7 min of measurement + ~1 min build (inputs up to 50 Mbp) |
| `run_sdsl.sh` | sdsl-lite v2.1.1 (GPL-3.0) | build sdsl-lite (install.sh), compile `sdsl/sdsl_bench.cpp`; FM-index `csa_wt<wt_huff<>>`, `csa_wt<wt_huff<rrr_vector<127>>,512,1024>`, `csa_sada<>`, CST `cst_sct3<>` | git cmake g++ make | ~50 min measurement (RRR FM-index alone 32 min) + ~1.5 min sdsl build |
| `run_rustbio.sh` | rust-bio (crate `bio`, MIT) | `cargo build --release`, SA-IS + BWT + FM-index, full-SA locate | cargo, crates.io | ~2 min measurement + ~1 min cargo build |
| `run_python.sh` | pydivsufsort (MIT), `suffix-trees` (pure-Python Ukkonen) | venv + pip; SA build + `sa_search` queries; suffix-trees at 10 k ... 1 M chars until it gives up | python3 venv, PyPI | pydivsufsort ~2.5 min; suffix-trees 1 Mbp ~13 s/run, 10 Mbp ~5 min/run (7.5 GB RAM), 20 Mbp needs > 12 GB and is stopped by the `ulimit -v` cap after ~13 min |
| `run_aligners.sh` | BWA `index -a bwtsw` (GPL-3.0), `bowtie2-build` (GPL-3.0) | FM-index build time, peak RSS, index bytes only (no count/locate API) | `bwa` + `bowtie2` on PATH, or `BWA_BIN` / `BT2_BUILD_BIN` | ~2.5 min measurement + ~4 min to build both from source |

Each script prints a `[BLOCKED/FAILED]` line instead of aborting when one (tool, input) fails; failed repeats are recorded in the raw file as `failed_rc` rows with the error text.
Pinned source URLs / commit SHAs / licences are written to `competitors/_work/SOURCES.txt` on every run.

## 3. Things worth doing on a bigger machine (not done here)

* libsais/libdivsufsort at 1 Gbp (`libsais64`) -- needs >= 16 GB RAM; add a 1 Gbp input with the generator's `rand_dna`.
* sdsl-lite `cst_sct3` at >= 20 Mbp, and `csa_wt` with `bit_vector`-based sampling variations (the RRR + 512 sampling variant is extremely slow at locate, by design).
* BWA/Bowtie2 `-a is` vs `bwtsw`, `--threads 4` build.
* ropebwt3 / r-index / grlBWT (pangenome-scale, run-length FM-index) -- scripts not written because they need a repeat-rich multi-genome input (the `rep_10m` text is the closest proxy here).
