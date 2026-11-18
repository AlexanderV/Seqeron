### Environment (preflight rows)

| item | value |
|---|---|
| machine | aarch64 |
| image built for | arm64 |
| /proc/cpuinfo format | arm |
| emulation | no |
| cores | 18 |
| CPU quota | max |
| memory visible (bytes) | 8317202432 |
| cgroup version | 2 |
| container memory limit | max |
| container swap allowance | max |
| disk size (bytes) | 1962694664192 |
| disk free at preflight (bytes) | 1860086378496 |
| file system | ext2/ext3/ext4 |
| write / mmap / msync test | fallocate=ok mmap_write_msync_s=0.021 readback=ok write_mb_s=3196 |
| drop_caches possible | no |
| .NET | .NET 10.0.15 |
| process architecture | Arm64 |
| suffix comparison path | scalar 64-bit words |
| page size | 4096 |
| preflight result | OK |

### Generous runs: time, memory, I/O and disk (steps a, b, c; G4 reported, G2 and G3 judged here)

| case | run | finished | total s | process s | PH0 | PH1 | PH2/3 | PH4+5 | ns/sym | RssAnon prog MB | RssAnon sampler MB | I/O counted B/char | rchar+wchar B/char | device read B/char | scratch peak B/char (prog) | scratch peak B/char (sampler) | outputs B/char | peak disk GB | cold start |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n1000m.gen | 1 | yes | 123.3 | 123.7 | 2.8 | 9.9 | 37.0 | 73.6 | 123 | 169.8 | 169.8 | 75.5 | 67.5 | 20.02 | 12.00 | 12.00 | 9.50 | 18.2 | no |
| n1000m.gen | 2 | yes | 124.5 | 124.9 | 2.7 | 9.8 | 36.5 | 75.5 | 125 | 172.0 | 172.0 | 75.5 | 67.5 | 20.32 | 12.00 | 11.88 | 9.50 | 18.4 | yes |
| n2200m.gen | 1 | yes | 291.3 | 291.7 | 5.9 | 21.9 | 86.8 | 176.6 | 132 | 226.3 | 226.3 | 75.6 | 67.7 | 28.17 | 12.00 | 12.00 | 9.51 | 40.7 | yes |
| n3000m.gen | 1 | NO (no) | 406.8 | 407.2 | 8.3 | 30.8 | 124.5 | 243.1 | 136 | 309.7 | 309.7 | 75.4 | 67.5 | 29.88 | 12.00 | 11.94 | 9.54 | 55.6 | yes |
| n3000m.gen | 2 | yes | 411.5 | 411.9 | 8.6 | 30.5 | 122.5 | 249.9 | 137 | 309.2 | 309.2 | 75.4 | 67.6 | 28.79 | 12.00 | 11.94 | 9.54 | 55.6 | yes |

Phase seconds are the program's own; "process s" includes the runtime start. ESTIMATE from E7 for comparison: 250 ns per symbol (205-300), i.e. 4.2 / 9.2 / 12.5 min at 1 / 2.2 / 3 Gbp on the prototype machine (4 vCPU x86, virtio disk).

### Verifier (e8-verify) after the generous runs (G1)

| case | run | all | verify s | files | meta | text | input | permutation | order | lcp values | lcp width | stats | child | sha256 | long pairs | LCE cap | stack peak |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n1000m.gen | 1 | pass | 43.4 | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | 8106796 | 65536 | 5000000 |
| n1000m.gen | 2 | pass | 43.2 | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | 8106796 | 65536 | 5000000 |
| n2200m.gen | 1 | pass | 102.3 | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | 18306793 | 65536 | 11000000 |
| n3000m.gen | 1 | not run | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| n3000m.gen | 2 | pass | 125.4 | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass | 25106794 | 65536 | 3000002 |

### Identity of the outputs against the generous run of the same text (G1)

| case | run | exit code | seven files SHA-256 equal |
|---|---|---|---|
| n1000m.gen.t1 | 1 | 0 | yes |
| n1000m.gen.t2 | 1 | 0 | yes |
| n1000m.m1216 | 1 | 0 | yes |
| n1000m.m1600 | 1 | 0 | yes |
| n1000m.m2080 | 1 | 0 | yes |
| n1000m.m3024 | 1 | 0 | yes |
| n1000m.m640 | 1 | 0 | yes |
| n1000m.m768 | 1 | 0 | yes |
| n1000m.m928 | 1 | 0 | yes |

### The limited-cache sweep (step d; G4 reported): the container limit, the effective cache, the time against the generous run

| text | container memory MiB | result | effective cache of text + rank in PH3 | total s | slowdown | PH0 s | PH1 s | PH2/3 s | PH4+5 s | device read B/char | cgroup failcnt / events | RssAnon sampler MB | residency sampled |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n1000m | 3024 | finished (1 run) | 1.000 | 129.8 | 1.05x | 2.8 | 10.2 | 38.2 | 78.6 | 27.25 | 34425 | 162.8 | yes |
| n1000m | 2080 | finished (1 run) | 1.000 | 128.5 | 1.04x | 2.8 | 10.1 | 38.5 | 77.1 | 28.87 | 46641 | 156.2 | yes |
| n1000m | 1600 | finished (1 run) | 1.000 | 130.5 | 1.06x | 2.8 | 10.2 | 38.6 | 79.0 | 29.67 | 51118 | 163.4 | yes |
| n1000m | 1216 | finished (1 run) | 0.999 | 131.6 | 1.07x | 2.8 | 10.2 | 39.4 | 79.2 | 31.17 | 62397 | 160.8 | yes |
| n1000m | 928 | finished (1 run) | 0.999 | 132.0 | 1.07x | 2.9 | 10.4 | 39.8 | 78.9 | 34.48 | 103982 | 158.4 | yes |
| n1000m | 768 | finished (1 run) | 0.970 | 133.2 | 1.08x | 2.8 | 10.3 | 41.7 | 78.4 | 38.23 | 117855 | 146.0 | yes |
| n1000m | 640 | finished (1 run) | 0.858 | 146.6 | 1.19x | 2.9 | 10.8 | 45.6 | 87.3 | 52.93 | 193164 | 170.7 | yes |

Slowdown = total_s of the cell (median of its runs) / total_s of run 1 of the generous run of the same text. 'Effective cache' = time-weighted resident fraction of the packed text and the DCS rank file over PH3 (mincore, 200 ms), as in E7. Abandon rule: a run that takes more than 6 times the generous run of the same text is stopped and reported as 'does not finish in practice at this memory'.

#### Verdict: the smallest memory with a slowdown of at most 3x (the E7 criterion)

| text | smallest container memory with slowdown <= 3x | slowdown | effective cache | as a multiple of text + rank file | per symbol |
|---|---|---|---|---|---|
| n1000m | 640 MiB | 1.19x | 0.858 | 1.33 x (text + rank file = 481 MiB) | 0.67 B/char |

The open question of ADR-004 6c (3): how much page cache PH4's 12 B/char of scratch needs beyond the 1x (text + rank file). Read it from the table above: the container limit includes the anonymous memory (about 0.1-0.15 GB), the page cache of the text, the rank file and the scratch pages in flight; the first cell whose PH4+5 time grows against the generous run is the one where the scratch no longer fits next to the text. The interpretation is for the orchestrator; the numbers are above.

### Gates

| gate | text (verbatim) | result | numbers |
|---|---|---|---|
| G1 | the verifier passes at 1, 2.2 and 3 Gbp and the output SHA-256 are identical between the generous and every limited-cache run of the same text and between thread counts | PASS | verifier by size: 1000m pass, 2200m pass, 3000m pass; identity comparisons 9, differences 0; thread comparisons: n1000m ['yes', 'yes'] |
| G2 | independently sampled RssAnon of the dotnet process at most 256 MB at every size (E7: 124 MB at 500 Mbp) | FAIL | n1000m 172.0 MB, n2200m 226.3 MB, n3000m 309.7 MB |
| G3 | peak scratch at most 12.5 B/char and peak total disk at most 70 GB at 3 Gbp | PASS | scratch peak 11.94/12.00 B/char, peak total disk 55.6 GB (3 Gbp generous run) |
| G4 | reported, not gated: time per phase and total (ESTIMATE from E7 for 3 Gbp: 11-15 min on the prototype machine, linear), I/O in B/char, and the cache table of (d) with the answer to the open question | reported | see the tables above |

Abandon rule: a run that takes more than 6 times the generous run of the same text is stopped and reported as 'does not finish in practice at this memory'. Cells stopped by it: none.

