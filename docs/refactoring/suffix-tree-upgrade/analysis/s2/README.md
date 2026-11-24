# analysis/s2: the measurements and prototypes of work package S2-0 (ADR-005)

Everything here is design-time evidence for [ADR-005](../../adr/ADR-005-genome-scale.md). Nothing here is shipped code and nothing is part of `Seqeron.sln`: the one project (`widthbench/`) is built by hand, managed code only (BCL `NativeMemory`
and `MemoryMappedFile`; no P/Invoke, no native binary), and the only file outside this directory that it reads is the product project `SuffixTree.Esa` (as the oracle of one check). Machine of every row: 4 vCPU Intel Xeon (KVM), 15 GB RAM, Linux,
.NET 10.0.403, Release, tiered compilation and dynamic PGO at their defaults, `DOTNET_GCHeapHardLimit` unset, nothing else running (the scripts print the load and the busiest processes before they time anything). Scratch is `$S2_OUT`
(default `/tmp/stbench/s20`), never inside the repository.

## Reproduce

```
export S2_OUT=/tmp/stbench/s20                      # scratch; D=/tmp/stbench holds the S1 bench inputs (data/, queries/) and b9/sa_bench_libsais (analysis/bench/run_1b9_gate.sh builds it)
S2=docs/refactoring/suffix-tree-upgrade/analysis/s2
dotnet build $S2/widthbench/WidthBench.csproj -c Release -o $S2_OUT/bin/wb
dotnet build apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj -c Release -o $S2_OUT/bin/bench      # the production builder (internal) through its command line
dotnet build apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj -c Release -o $S2_OUT/bin/emproto        # only for the genome-like generator e8-gen
dotnet $S2_OUT/bin/wb/WidthBench.dll selftest                                                           # every layout reads back every value it stored (n = 1 .. 10^6)
```

| row file | what | command | time on this machine |
|---|---|---|---|
| `results/width_synthetic.tsv` kinds `thr`, `scan`, `null` | independent random reads, sequential slices, the decode floor, 5 layouts x n = 10^7 .. 10^9 x 1 and 4 threads | `scripts/run_width_synthetic.sh` | 12 min |
| `results/width_synthetic.tsv` kind `chase` | dependent random reads (the latency), the chain restarts every 1,024 steps | `scripts/run_width_chase.sh` | 5 min |
| `results/width_descent.tsv`, labels `dna_10m`, `dna_100m`, `genome_100m`, `rep_10m` | the count and locate descent of a real container with the stored width swapped | `scripts/run_width_descent.sh <container> <label> <variants> [reps] [queries]` | 3-6 min per container |
| `results/width_descent.tsv`, label `dna_500m group ...` | the same at n = 5x10^8 in groups that fit the RAM | `scripts/run_width_descent_500m.sh <container>` | about 25 min |
| `results/text_bench.tsv` | the read cost of 2-bit and 3-bit packed DNA against the byte text (random read, dependent read, scan, longest common extension) | `scripts/run_text_bench.sh` | 3 min |
| `results/child_stats.txt` | how far the child cells and PSI increments are from their own index (the delta forms) on random DNA, a genome-like text and a repeat text | `widthbench stats --container <file>` | seconds |
| `results/bchar_table.tsv` | bytes per character per section for each width at n = 10^8 .. 6.2x10^9 (arithmetic from the section layouts) | `python3 scripts/bchar.py --tsv results/bchar_table.tsv` | instant |
| `results/int_inventory.tsv` | every place where a 32-bit type limits a size or a position (file, line, category), by pattern | `python3 scripts/int_inventory.py > results/int_inventory.tsv` | instant |
| `results/drift_s2.tsv` | the libsais `dna_10m` single-thread build (S1 reference 0.328 s) at the start and end of every session; more than 10 % off = "drift" | `scripts/drift.sh <label>` | 1 min |
| `results/build_chain_rows.txt` | the production chain (Bench `sesa build`) on 10^8 and 5x10^8 random DNA, Core and Full: phase times, RssAnon, container sizes | see below | 18 s / 80 s / 105 s |
| `results/libsais_check_500m.txt` | SHA-256 of the libsais suffix array of the 5x10^8 text against SA[1..n] of the container (the checksum method of the S2 gate; equal) | `scripts/check_sa_vs_libsais.sh <text> <container> <libsais source dir>` | 40 s, 2 GB of disk |

Containers used by the descent rows (all built by `sesa build` of the Bench, 4 threads, default parameters):

```
python3 $S2/scripts/gen_dna.py 500000000 $S2_OUT/data/dna_500m.txt           # the first 10^8 bytes are dna_100m (SHA-256 32bca5fd...), the first 10^7 are dna_10m
dotnet $S2_OUT/bin/bench/SuffixTree.Esa.Bench.dll sesa build --input <text> --work <dir> --label <name> --threads 4 --keep true --rss true [--psi false]
dotnet $S2_OUT/bin/emproto/SuffixTree.EmProto.dll e8-gen --n 100000000 --out $S2_OUT/data/genome_100m.txt --force --case gen.genome
```

## Reading the rows

* `kind layout n threads metric median min max unit reps notes`. `median`, `min`, `max` are over the repetitions (5-11); the repetitions of the descent rows are **interleaved** across the variants (the order rotates), so that a slow drift of the machine
  hits every variant alike; the first two passes of every variant are warm-up (tiered compilation with PGO reaches tier 1 as in the product). The ratio columns in the ADR are medians of those rows.
* `chase` rows: ns per read of one thread (the wall time of the slowest thread divided by its reads); `thr` and `scan`: wall time divided by the reads of all threads (the aggregate throughput as ns per read); `descent`: ns per query per thread
  (wall time of a pass / (patterns x threads)).
* Layouts: `w4` 32-bit entries, `w5` 40-bit byte-packed entries, `w8` 64-bit, `hyb` the hybrid (32-bit slots, escape value, 16-byte record per 64 slots, 40-bit long values; the share of escapes is that of n = 6.2x10^9: 30.7 %), `p34` 34 bits per entry,
  15 per 64-byte line. Descent variants: `v0` int arithmetic (the code of `MappedSearch.Run`), `v1` the same with `long` (w4 tagged cells), `v1b` a second copy of `v1` (**the null experiment of the instrument**), `v2` w4 packed cells,
  `v3` w5 tagged, `v4` w5 packed, `v5` w8 tagged, `v6` hybrid suffix array with w5 cells, `v7` p34, `v8` / `v9` delta child cells (2 bytes tagged / 1 byte packed) with a binary-search escape list, `v8h` / `v9h` the same with a hash table for the escapes.
* The first version of the `chase` loop followed a random function and fell into its cycle of about sqrt(n) entries, which stays in the cache; those rows were removed and re-measured with a chain that restarts every 1,024 steps (the file
  `results/width_synthetic.tsv` has only the valid rows; the script `run_width_chase.sh` says so).
* The synthetic arrays are anonymous memory (`NativeMemory.AlignedAlloc`, 4 KiB pages: `/sys/kernel/mm/transparent_hugepage/enabled` is `madvise` and .NET does not call `madvise`, so a mapped file and these arrays have the same page size).
  The descent rows read the container's text and LCP in place and copy the suffix array and the child table into anonymous memory in the layout of the variant.

## Limits of these measurements (stated, not hidden)

* The arrays are at most 10^9 entries (8 GB for the 8-byte layout) because of the 15 GB of RAM; the page cache is not involved (anonymous memory). What 6.2x10^9 entries (25-50 GB per array) would add - page-table walks that miss
  the cache (the PTEs of a 50 GB array are 100 MB), arrays larger than RAM, NUMA - cannot be measured here (ADR-005 section 5).
* The descent is measured on real containers of 10^7, 10^8 and 5x10^8 symbols only; above that the cost model of the ADR is arithmetic over the synthetic rows.
* A spread of 25-40 % between the minimum and the maximum of a single-thread row is normal on this VM (steal time); the four-thread rows and the interleaved descent rows have a spread of 3-10 % and the null experiment (`v1` against `v1b`)
  shows what the instrument resolves (0.95-1.03).
