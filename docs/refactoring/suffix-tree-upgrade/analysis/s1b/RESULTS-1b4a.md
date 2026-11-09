# 1b.4a results: prototypes P1 (suffix-link kernel) and P2 (end-to-end MEM)

Work package 1b.4a of S1b (ADR-001 section 5). Decision and discussion: PROGRESS-S1.md, Log "1b.4a". This file is the data guide and the extended tables.

## Provenance

* Code: `apps/SuffixTree.Esa.Bench` (commands `links`, `mem`; prototype library `Links/`), measured with commit `71c5daba` plus the two-line fix of `LinksCommands.cs` that is committed in `cd6f63d2` (the `# process` headers therefore say `71c5dabaf3-dirty`; the fix only makes the cross-check skip links that were not computed). Tests: `tests/SuffixTree/SuffixTree.Esa.Tests/Links` (33 default + 5 explicit `Exhaustive`).
* Machine: 4 vCPU `Intel(R) Xeon(R) Processor` VM, 16 GB, `lscpu` reports L2 2 MiB per core and an L3 of 480 MiB (shared with the host, so every structure of the 10 Mbp index, about 180 MB, is cache-resident here: random reads cost L3 latency, not DRAM latency; `dna_50m` below has about 900 MB). Release build, defaults of the runtime (tiered compilation, dynamic PGO on, workstation non-concurrent GC), the kernels are `AggressiveOptimization`. The machine was quiet during every measurement (`uptime`, `ps`; build servers shut down before).
* Inputs: `/tmp/stbench` (`gen_inputs.py` files, `MANIFEST.sha256`), `a^n`, Fibonacci, Thue-Morse from the bench app's `--synthetic` generators (1,000,000 symbols; 125k / 250k / 500k for the scaling rows). The P2 query is `queries/<name>.memq100k.txt` of the same generator (100,000 bp window with 1 % substitutions), minLength 20, the query of `BENCHMARK.md` section 4.
* Files: `p1_links_raw.tsv` (30 headline processes: 6 inputs x 5; 9 scaling processes; 3 x `dna_50m`), `p1_fanout_raw.tsv` (pyramid fanout 8 / 16 / 32 / 64, `dna_10m` and `rep_10m`, 3 processes each), `p1_summary.tsv` (medians over the processes), `p2_mem_raw.tsv` (60 processes), `aggregate.py` (all tables below come from `aggregate.py p1 / p2`), `run_p1.sh` / `run_p2.sh` (drivers, one fresh process per run, the kernel order rotated by the run number), `mutation_p1_runner.py` and `mutation_p1_result.txt` (hand-written mutants of the prototype code).
* Row format: `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes` as the other bench output; a `# process input=... run=... commit=... options: ...` line precedes the rows of each process.
* ns per link are **net**: the median of back-to-back `Stopwatch.GetTimestamp` pairs (22 ns here) is subtracted; `max ns` is the raw largest single link of any process and is dominated by scheduler noise, so the boundedness argument uses steps (deterministic), not max ns.

## P1: ns per link, steps (sample of 200,000 internal nodes, seed 1; 145,561 for `prot_1m`; medians over 5 fresh processes)

(raw per-process ranges are in `aggregate.py p1` output; the `a^n` descent row counts 259 links because the total budget of 4 x 10^8 steps was used: budget exceeded after 259 of 200,000)

| input | kernel | median ns ISA / Psi | p99 ns ISA / Psi | max ns raw ISA / Psi | steps median / p99 / max | links counted (budget exceeded) |
|---|---|---|---|---|---|---|
| dna_10m | (a) scan | 296 / 315 | 547 / 567 | 366410 / 4902670 | 4 / 18 / 46 | 200000 (0) |
| dna_10m | (b) descent | 1278 / 1165 | 2363 / 2236 | 768859 / 150538 | 37 / 48 / 57 | 200000 (0) |
| dna_10m | (c) RMQ + PSV/NSV | 410 / 415 | 919 / 869 | 173924 / 252871 | 17 / 120 / 304 | 200000 (0) |
| dna_10m | (c2) PSV/NSV | 324 / 311 | 619 / 575 | 211505 / 213749 | 4 / 19 / 45 | 200000 (0) |
| rep_10m | (a) scan | 408 / 354 | 763 / 689 | 253609 / 135208 | 2 / 96 / 658 | 200000 (0) |
| rep_10m | (b) descent | 1608 / 1553 | 3010 / 2665 | 5368370 / 163202 | 117 / 227 / 258 | 200000 (0) |
| rep_10m | (c) RMQ + PSV/NSV | 580 / 550 | 1030 / 1010 | 115733 / 2688950 | 43 / 112 / 310 | 200000 (0) |
| rep_10m | (c2) PSV/NSV | 444 / 361 | 847 / 744 | 235671 / 2276760 | 2 / 27 / 79 | 200000 (0) |
| prot_1m | (a) scan | 153 / 153 | 362 / 380 | 100319 / 202266 | 12 / 104 / 204 | 145561 (0) |
| prot_1m | (b) descent | 372 / 385 | 1024 / 1049 | 130768 / 604985 | 40 / 70 / 90 | 145561 (0) |
| prot_1m | (c) RMQ + PSV/NSV | 268 / 274 | 570 / 590 | 362836 / 235193 | 53 / 143 / 277 | 145561 (0) |
| prot_1m | (c2) PSV/NSV | 163 / 168 | 353 / 368 | 864600 / 3258180 | 12 / 50 / 67 | 145561 (0) |
| a^n | (a) scan | 161 / 157 | 354 / 329 | 58848 / 167808 | 2 / 2 / 2 | 200000 (0) |
| a^n | (b) descent | 2037610 / 2037980 | 4024330 / 4056480 | 5622990 / 4479440 | 1522360 / 2994670 / 2999390 | 259 (0) |
| a^n | (c) RMQ + PSV/NSV | 358 / 366 | 615 / 649 | 137611 / 4782850 | 108 / 150 / 176 | 200000 (0) |
| a^n | (c2) PSV/NSV | 162 / 161 | 300 / 337 | 79847 / 293748 | 2 / 2 / 2 | 200000 (0) |
| fibonacci | (a) scan | 141 / 152 | 256 / 279 | 185480 / 225005 | 2 / 2 / 90168 | 200000 (0) |
| fibonacci | (b) descent | 507 / 489 | 904 / 911 | 735145 / 224351 | 79 / 100 / 114 | 200000 (0) |
| fibonacci | (c) RMQ + PSV/NSV | 193 / 187 | 540 / 522 | 100062 / 107989 | 9 / 105 / 279 | 200000 (0) |
| fibonacci | (c2) PSV/NSV | 127 / 141 | 256 / 282 | 857012 / 544136 | 2 / 2 / 79 | 200000 (0) |
| thue-morse | (a) scan | 117 / 132 | 234 / 283 | 144557 / 261477 | 2 / 2 / 2609 | 200000 (0) |
| thue-morse | (b) descent | 418 / 445 | 769 / 843 | 346567 / 135189 | 54 / 66 / 73 | 200000 (0) |
| thue-morse | (c) RMQ + PSV/NSV | 182 / 187 | 518 / 551 | 297638 / 108284 | 7 / 94 / 307 | 200000 (0) |
| thue-morse | (c2) PSV/NSV | 125 / 137 | 250 / 282 | 349736 / 1371390 | 2 / 2 / 63 | 200000 (0) |

## P1: cost of the navigation structures (medians [min..max] over 5 processes)

| input | nodes | ISA build s (scattered pass) | Psi build s (from SA + ISA) | pyramid build s (F = 16) | ISA / Psi / pyramid bytes per char |
|---|---|---|---|---|---|
| dna_10m | 6229745 | 0.0708 [0.0705..0.0781] | 0.0766 [0.0732..0.0858] | 0.0273 [0.0269..0.0277] | 4.000 / 4.000 / 0.267 |
| rep_10m | 8568247 | 0.0866 [0.0780..0.0920] | 0.0906 [0.0893..0.0949] | 0.0332 [0.0323..0.0437] | 4.000 / 4.000 / 0.267 |
| prot_1m | 295123 | 0.0063 [0.0059..0.0079] | 0.0056 [0.0053..0.0066] | 0.0054 [0.0050..0.0105] | 4.000 / 4.000 / 0.267 |
| a^n | 999999 | 0.0025 [0.0023..0.0026] | 0.0031 [0.0027..0.0032] | 0.0050 [0.0049..0.0054] | 4.000 / 4.000 / 0.267 |
| fibonacci | 999995 | 0.0087 [0.0077..0.0092] | 0.0062 [0.0058..0.0064] | 0.0052 [0.0051..0.0058] | 4.000 / 4.000 / 0.267 |
| thue-morse | 999902 | 0.0133 [0.0111..0.0134] | 0.0105 [0.0087..0.0108] | 0.0054 [0.0052..0.0060] | 4.000 / 4.000 / 0.267 |

## P1: ISA against Psi (per fresh process, both measured in the same process)

| input | (a) scan | (b) descent | (c) rmq | (c2) psv |
|---|---|---|---|---|
| dna_10m | 1.02 [0.99..1.13] | 1.00 [0.77..1.21] | 0.95 [0.92..1.05] | 1.04 [0.90..1.08] |
| rep_10m | 0.86 [0.71..0.89] | 0.99 [0.84..0.99] | 0.94 [0.92..0.95] | 0.87 [0.81..0.91] |
| prot_1m | 1.01 [0.95..1.03] | 1.03 [0.96..1.06] | 1.01 [0.98..1.04] | 1.00 [0.99..1.04] |
| a^n:125000  [--sample 50000 --worst 300] | 1.06 [1.06..1.06] | 0.96 [0.96..0.96] | 1.00 [1.00..1.00] | 1.09 [1.09..1.09] |
| a^n:250000  [--sample 50000 --worst 300] | 0.90 [0.90..0.90] | 0.96 [0.96..0.96] | 1.03 [1.03..1.03] | 1.05 [1.05..1.05] |
| a^n:500000  [--sample 50000 --worst 300] | 0.99 [0.99..0.99] | 0.99 [0.99..0.99] | 1.03 [1.03..1.03] | 1.01 [1.01..1.01] |
| a^n:1000000 | 1.01 [0.86..1.01] | 1.01 [0.95..1.04] | 1.02 [0.68..1.42] | 1.00 [0.98..1.02] |
| fibonacci:125000  [--sample 50000 --worst 300] | 0.95 [0.95..0.95] | 1.04 [1.04..1.04] | 0.98 [0.98..0.98] | 0.96 [0.96..0.96] |
| fibonacci:250000  [--sample 50000 --worst 300] | 1.06 [1.06..1.06] | 1.06 [1.06..1.06] | 1.14 [1.14..1.14] | 1.04 [1.04..1.04] |
| fibonacci:500000  [--sample 50000 --worst 300] | 1.12 [1.12..1.12] | 1.05 [1.05..1.05] | 0.99 [0.99..0.99] | 1.11 [1.11..1.11] |
| fibonacci:1000000 | 1.10 [1.08..1.12] | 1.00 [0.94..1.05] | 0.97 [0.97..1.00] | 1.10 [1.08..1.11] |
| thue-morse:125000  [--sample 50000 --worst 300] | 0.90 [0.90..0.90] | 1.05 [1.05..1.05] | 1.00 [1.00..1.00] | 1.00 [1.00..1.00] |
| thue-morse:250000  [--sample 50000 --worst 300] | 0.75 [0.75..0.75] | 1.02 [1.02..1.02] | 0.86 [0.86..0.86] | 0.98 [0.98..0.98] |
| thue-morse:500000  [--sample 50000 --worst 300] | 1.20 [1.20..1.20] | 1.04 [1.04..1.04] | 1.02 [1.02..1.02] | 1.19 [1.19..1.19] |
| thue-morse:1000000 | 1.11 [1.07..1.19] | 1.04 [1.01..1.10] | 1.02 [0.98..1.03] | 1.09 [0.87..1.11] |
| dna_50m  [--sample 200000 --worst 300] | 1.04 [1.04..1.05] | 1.00 [0.86..1.01] | 0.99 [0.98..1.01] | 1.08 [1.04..1.12] |

## P1: boundedness by text size

| input | max scan cost over ALL nodes | max node depth | (a) steps max sample+worst | (b) steps max sample+worst | (c) steps max sample+worst | (c2) steps max sample+worst |
|---|---|---|---|---|---|---|
| dna_10m | 50 | 15 | 52 | 57 | 405 | 46 |
| rep_10m | 780 | 97 | 782 | 273 | 381 | 87 |
| prot_1m | 273 | 7 | 275 | 90 | 281 | 68 |
| a^n:125000  [--sample 50000 --worst 300] | 1 | 124999 | 2 | 374994 | 157 | 2 |
| a^n:250000  [--sample 50000 --worst 300] | 1 | 249999 | 2 | 749994 | 179 | 2 |
| a^n:500000  [--sample 50000 --worst 300] | 1 | 499999 | 2 | 1499994 | 163 | 2 |
| a^n:1000000 | 1 | 999999 | 2 | 2999994 | 176 | 2 |
| fibonacci:125000  [--sample 50000 --worst 300] | 29504 | 36 | 29506 | 96 | 292 | 71 |
| fibonacci:250000  [--sample 50000 --worst 300] | 59013 | 39 | 59015 | 105 | 332 | 84 |
| fibonacci:500000  [--sample 50000 --worst 300] | 118027 | 40 | 118029 | 108 | 300 | 78 |
| fibonacci:1000000 | 236064 | 42 | 236066 | 114 | 325 | 79 |
| thue-morse:125000  [--sample 50000 --worst 300] | 10418 | 22 | 10420 | 63 | 266 | 58 |
| thue-morse:250000  [--sample 50000 --worst 300] | 20834 | 23 | 20836 | 67 | 259 | 68 |
| thue-morse:500000  [--sample 50000 --worst 300] | 41669 | 24 | 41671 | 70 | 354 | 65 |
| thue-morse:1000000 | 83334 | 25 | 83336 | 74 | 352 | 74 |
| dna_50m  [--sample 200000 --worst 300] | 67 | 17 | 69 | 64 | 473 | 53 |

(a) grows linearly with n on Fibonacci (exact worst case over all nodes: 29,504 / 59,013 / 118,027 / 236,064 steps at n = 125k / 250k / 500k / 1M, about 0.236 n) and on Thue-Morse (about n / 12); (b) grows as 3 n on `a^n` (node depth n - 1); (c) stays at 157-405 steps and (c2) at 2-123 steps on every input and size.

## P1: pyramid fanout (kernels (c) and (c2), Psi, 3 processes per cell)

| input | fanout F | pyramid B/char | pyramid build s | (c) median ns | (c) steps p99 / max | (c2) median ns | (c2) steps p99 / max |
|---|---|---|---|---|---|---|---|
| dna_10m | 8 | 0.571 | 0.036 | 504 [498..534] | 88 / 219 | 389 [383..409] | 17 / 30 |
| dna_10m | 16 | 0.267 | 0.027 | 472 [447..494] | 120 / 304 | 368 [361..376] | 19 / 45 |
| dna_10m | 32 | 0.129 | 0.024 | 452 [412..465] | 180 / 548 | 370 [356..372] | 18 / 48 |
| dna_10m | 64 | 0.063 | 0.021 | 403 [399..464] | 284 / 861 | 338 [336..374] | 18 / 46 |
| rep_10m | 8 | 0.571 | 0.035 | 611 [551..822] | 82 / 246 | 445 [390..449] | 22 / 54 |
| rep_10m | 16 | 0.267 | 0.032 | 617 [607..660] | 112 / 310 | 371 [309..414] | 27 / 79 |
| rep_10m | 32 | 0.129 | 0.032 | 671 [662..693] | 178 / 458 | 406 [367..411] | 38 / 123 |
| rep_10m | 64 | 0.063 | 0.029 | 620 [606..716] | 284 / 936 | 406 [345..414] | 67 / 222 |

## P1: `dna_50m` (31,181,032 internal nodes, 3 processes; the structures total about 900 MB)

| kernel | median ns ISA / Psi | p99 ns ISA / Psi | steps median / p99 / max |
|---|---|---|---|
| (a) scan | 579 / 607 | 874 / 899 | 4 / 17 / 45 |
| (b) descent | 2188 / 1985 | 3752 / 3450 | 41 / 53 / 64 |
| (c) RMQ+PSV/NSV | 681 / 667 | 1446 / 1330 | 18 / 120 / 362 |
| (c2) PSV/NSV | 574 / 633 | 981 / 1029 | 4 / 18 / 44 |

ISA build 0.625 s, Psi build 0.723 s, pyramid 0.135 s; 31181032 internal nodes.

## P2: end-to-end MEM, 100 kbp query with 1 % substitutions, minLength 20

dna_10m: 817 MEMs, rep_10m: 160,979 MEMs (the counts of `BENCHMARK.md` section 4). Counters of the driver (all variants of an input equal except the loci): 100,000 query starts; 2 loci: 200,000 links followed, 200,000 extension steps, dna_10m 21,707 child lookups (extension from a node) and 342,090 rescan lookups, rep_10m 70,165 and 186,061; 1 locus: 100,000 links, 100,000 extension steps, dna_10m 10,852 / 171,047, rep_10m 21,537 / 114,312. The MEM set of the driver equals the set of `FindMaximalExactMatches` of the in-memory Ukkonen tree of the same process in all 30 comparisons (3 headline variants x 5 processes x 2 inputs); the checksum (sum of r * 31 + q * 17 + len) is identical in all 60 processes of an input (dna_10m 77915707478, rep_10m 25158196373282).

**dna_10m** (100,000 starts)

| variant | fresh processes | first call s (cold JIT) | steady median s (last 5 of 10 calls) | share of the time in the link kernel |
|---|---|---|---|---|
| (c) RMQ+PSV/NSV, Psi, 2 loci (ADR kernel, generic-algorithm cost) | 5 | 0.130 [0.125..0.150] | 0.097 [0.090..0.106] | 54 % |
| (c), ISA, 2 loci | 5 | 0.119 [0.114..0.127] | 0.090 [0.080..0.101] | 51 % |
| (c2) PSV/NSV, Psi, 2 loci | 5 | 0.103 [0.094..0.106] | 0.073 [0.065..0.077] | 49 % |
| (c2), ISA, 2 loci | 3 | 0.090 [0.083..0.105] | 0.066 [0.055..0.071] | 42 % |
| (c2), Psi, 1 locus (ESA matching statistics) | 3 | 0.086 [0.079..0.090] | 0.053 [0.050..0.068] | 50 % |
| (c), Psi, 1 locus | 3 | 0.108 [0.104..0.108] | 0.064 [0.058..0.069] | 53 % |
| (a) scan, Psi, 2 loci | 3 | 0.092 [0.089..0.140] | 0.065 [0.061..0.071] | 48 % |
| (b) descent, Psi, 2 loci | 3 | 0.208 [0.190..0.214] | 0.160 [0.156..0.185] | 74 % |
| Ukkonen in-memory tree, same process (15 processes) | 15 | 0.071 [0.065..0.083] | 0.048 [0.041..0.053] | - |

**rep_10m** (100,000 starts)

| variant | fresh processes | first call s (cold JIT) | steady median s (last 5 of 10 calls) | share of the time in the link kernel |
|---|---|---|---|---|
| (c) RMQ+PSV/NSV, Psi, 2 loci (ADR kernel, generic-algorithm cost) | 5 | 0.211 [0.209..0.269] | 0.160 [0.154..0.170] | 50 % |
| (c), ISA, 2 loci | 5 | 0.202 [0.185..0.218] | 0.151 [0.135..0.158] | 44 % |
| (c2) PSV/NSV, Psi, 2 loci | 5 | 0.177 [0.147..0.199] | 0.113 [0.107..0.117] | 40 % |
| (c2), ISA, 2 loci | 3 | 0.116 [0.114..0.118] | 0.085 [0.085..0.101] | 32 % |
| (c2), Psi, 1 locus (ESA matching statistics) | 3 | 0.113 [0.113..0.122] | 0.090 [0.081..0.100] | 38 % |
| (c), Psi, 1 locus | 3 | 0.163 [0.160..0.168] | 0.125 [0.119..0.128] | 48 % |
| (a) scan, Psi, 2 loci | 3 | 0.131 [0.128..0.144] | 0.103 [0.099..0.108] | 39 % |
| (b) descent, Psi, 2 loci | 3 | 0.263 [0.244..0.287] | 0.232 [0.231..0.270] | 67 % |
| Ukkonen in-memory tree, same process (15 processes) | 15 | 0.494 [0.450..0.550] | 0.386 [0.363..0.456] | - |

Paired ratios of the steady median per repetition (`aggregate.py p2`): (c2) / (c) with Psi, 2 loci: dna_10m 0.74 [0.68..0.80], rep_10m 0.69 [0.66..0.73]; ISA / Psi with (c): 0.85 [0.80..1.06] and 0.94 [0.79..1.00], with (c2): 0.88 and 0.76; 1 locus / 2 loci with (c2) Psi: 0.71 and 0.80; (a) scan / (c): 0.63 and 0.64; (b) descent / (c): 1.62 and 1.45.

## Mutation sanity of the prototype code

`mutation_p1_runner.py`, 22 hand-written mutants of `Links/*.cs` (helper, kernels, driver, node set, tables): 20 killed (H5 by a hang), 2 survive (K2: scan starts at `b` instead of `b + 1`; K5: PSV/NSV starts at `b` instead of `b + 1`) and are equivalent: `LCP[b] >= depth - 1` holds for every endpoint `b` inside the target, so the first position read is skipped either way and only the step count differs by one.
