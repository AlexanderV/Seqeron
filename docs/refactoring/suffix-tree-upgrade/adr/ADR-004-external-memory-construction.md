# ADR-004: External-memory construction of the persistent ESA (streaming straight into the SESA files)

Status: **proposed** (research and design step, 2026-11-04; no code was written, nothing was built or timed). Trigger: the categorical owner requirement of 2026-11-04 in [PROGRESS-S1.md](../PROGRESS-S1.md): the persistent index is built by streaming into the memory-mapped files, "build in RAM, then write" is rejected, working memory must not grow with n. Companions: [ADR-001](ADR-001-index-contract.md) (contract), [ADR-002](ADR-002-sa-construction.md) (in-RAM SA-IS; stays the in-memory builder), [ADR-003](ADR-003-esa-container-format.md) (SESA container). Paths are repo-relative.

Evidence marks, used on every number: **V** = read in the source text (PDF text or fetched page); **U** = secondhand (search snippet, citation inside another paper, page that could not be opened); **D** = derived here by arithmetic from V numbers (formula given); **E** = engineering estimate or hypothesis of this ADR, to be falsified by the prototype (section 5). Sources actually opened are listed in section 7. Published timings come from 2005-2023 hardware, most of them from 2005-2017.

## 0. Decision in brief

1. **Primary: "B-DCS"** = Karkkainen's blockwise suffix sorting with a difference-cover sample (lexicographic-order block generation, TCS 2007), adapted to external memory: (a) the difference-cover sample (DCS) ranks live in a scratch file produced by recursion that ends in the existing in-RAM SA-IS kernel; (b) suffix positions are distributed into lexicographic buckets of fixed capacity, every bucket is sorted in RAM and written at its final SA offset (so SA, then LCP and the child table, are produced in final order by independent, parallel tasks); (c) **LCP is not computed per bucket**: it is computed afterwards from the SA stream by a chunked external Phi/PLCP pass (two radix-distribution passes, at most 2n character comparisons whatever the text), which makes bucket boundaries irrelevant; (d) the child table comes from one scan of the LCP stream with an **external (spilling) stack**; (e) the navigation pack stores **Psi, not ISA** (Psi is built in one pass with sigma output streams). Published basis: Karkkainen 2007 (algorithm, worst-case O(n log n + v n), measured 6.4x slowdown on periodic DNA), Haag et al. 2025 (sketch of the same scan-per-bucket scheme for external memory). **The EM combination itself is not a published, measured algorithm; it is this ADR's hypothesis (E).**
2. **Fallback: SAscan family** (block SA by our SA-IS kernel + gap-array merge by streaming the text; pSAscan, MIT licence): robust and simple, but O(n^2/M) I/O, the final order exists only after the last merge, LCP needs a separate external pass. Switch to an fSAIS-class induced sorter only when n/M > ~75 (published crossover).
3. Biggest risk: random reads into the DCS rank file and the text depend on the OS page cache (rank file 0.2-0.3 B/char + 2-bit text 0.25 B/char); the prototype gate E3 must show the build does not degrade by more than 3x when the page cache is limited to half of that.
4. The persistent builder never needs a 64-bit SA-IS kernel: n appears only in record widths and file offsets; the int32 kernel is used only for a base case of at most ~16 M symbols.
5. Ranking for this requirement set: (1) B-DCS (C5): final-order bucket output, published worst-case bound, parallel, moderate effort; (2) SAscan/pSAscan (C4): robust and simple but quadratic I/O, no streamed final order or LCP; (3) fSAIS / eSAIS (C3 / C2): the only near-linear robust EM algorithms, but an external priority queue, 230-450 B/char (eSAIS) and, for eSAIS, GPL. Not recommended: pDC3 (superseded by eSAIS, 30n+ working space), BWTDisk (n^2 CPU, >= 10 h at this M), prefix-partitioned trees (O(n^2), tree not SA/LCP), CaPS-SA and PSAC (RAM / cluster), BCR (string collections).

## 1. Requirement and assumptions (the owner may correct any of them)

| id | assumption |
|---|---|
| A1 | Working set M <= 256 MB of **anonymous** memory (managed heap + native buffers + thread stacks + runtime baseline, 30-60 MB U from ADR-002 section 3) for every n. Page cache of mapped files is not counted, but access must be sequential or cache-friendly. |
| A2 | n up to 3x10^9 (hg19 file 3,070,128,194 V DKMS05 Table 1; GRCh38 2,945,849,068 and T2T-CHM13 3,117,292,071 V CaPS-SA Table 1) with 64-bit capable positions; no hard upper limit by design (both strands = 6.2x10^9 must work). |
| A3 | Alphabets: DNA sigma 4-5 (2-bit or 3-bit packed), bytes, UTF-16 and integer alphabets (dense order-preserving remap, ADR-002 section 2). |
| A4 | Output in SA order; SA, LCP, child table and pack sections are appended in rank order or at precomputed offsets by sequential writers. |
| A5 | Prototype machine: 4 vCPU, 15 GB RAM, ~22 GB free disk (ext4 on a virtio disk; NVMe not verified). Production assumption: NVMe 2-4 GB/s sequential (E). |
| A6 | The text is read through a read-only mapping and may be read randomly (page-cache resident: 0.25 B/char 2-bit, 1 B/char bytes). |
| A7 | **Scratch files outside the final index are allowed** (Q1). Even the DCS ranks are O(n) data (0.04-0.14 n entries) that cannot be RAM-resident at M = 256 MB for n = 3x10^9, so a build with no scratch at all is impossible for a robust algorithm. |
| A8 | The output bytes do not depend on thread count, M, bucket capacity or splitters (the SA and LCP are unique); BUILD_INFO may record the parameters. |
| A9 | No third-party code is run or copied. GPL sources (eSAIS, Bowtie2) were read only as documentation; the implementation is clean-room from the papers (PROCESS.md rule 1). |

## 2. Candidates

Table 2.1 cost model and memory. M = RAM, B = block, sort(x) = external sort of x words, scan(x) = x words read.

| id | algorithm (authors) | I/O and time | RAM / scratch disk |
|---|---|---|---|
| C1 | pDC3, pDC7 (Dementiev, Karkkainen, Mehnert, Sanders, ALENEX 2005 / JEA 2008) | DC3 sort(30n)+scan(6n) V (DKMS05 recurrence); DC7 sort(24.75n)+scan(3.5n), DC13 sort(30.1n)+scan(2.9n) V (thesis section 6.2) | M free (STXXL); peak scratch 46n (32-bit) / 58n (40-bit) V (thesis Table 6.1); "more than 30 n" working space V (FGM Fig. 2) |
| C2 | eSAIS (Bingmann, Fischer, Osipov, ALENEX 2013; LCP variant) | sort(17n)+scan(9n) V (thesis Thm 6.5) | M free; scratch 25n/28n (SA), 44n/54n (SA+LCP, in-RAM RMQ), 52n/63n (external RMQ) V (Table 6.1); the in-RAM RMQ grows to 300 MiB on the human genome V (thesis section 6.3.2) |
| C3 | fSAIS (Karkkainen, Kempa, Puglisi, Zhukova, ALENEX 2017); siblings DSA-IS, EM-SA-DS, SAIS-PQ | O(scan(n) log_{M/B}(n/B)) V (thesis Table 5.1) | scratch 8.2n (8.1n worst) V; EM-SA-DS needs M = Omega(sqrt(nB)), DSA-IS M = Omega(sqrt(nW)) V; SAIS-PQ 23n, SAIS-PQ+ 15n V |
| C4 | SAscan / pSAscan (Karkkainen, Kempa, Puglisi, SEA 2014 / CPM 2015); LCP by LCPscan, EM-SPhi, EM-SI (ESA 2016) | "O(n^2/M) before logarithmic factors" V (pSAscan README); worst case O(scan(n)(1 + (n/M) log_sigma n)) V (thesis section 6.4.1) | M free (tested 3.5-120 GiB); scratch 7.5n incl. input and 5n output V; 40-bit output V |
| C5 | **Karkkainen blockwise suffix sorting with difference cover (TCS 2007)**; implemented in Bowtie / Bowtie 2 as KarkkainenBlockwiseSA | worst case O(n log n + v n), average O(n log n + sqrt(v) n) on random text, extra space O(n log n / sqrt(v)) bits, v in [3, n^(2/3)] V (abstract, Thm 12); measured space 20n bits (v=16), 5n (v=256), 1.8n (v=2048) V | no scratch (BWT written directly, V); in the EM adaptation the DCS moves to a file: RAM = O(b_max + v + n/b_max) (E) |
| C6 | BWTDisk / bwte (Ferragina, Gagie, Manzini, LATIN 2010 / Algorithmica 2012) | BWT: O(n/M) scan passes, I/O O(n^2/(M B log n)), CPU O(n^2/M); SA: I/O O(n^2/(M B)); Psi: same class V (Thm 3.1, 4.2, 4.3) | n bits of disk working space V |
| C7 | Prefix-partitioned suffix trees: TDD, ST-Merge, Trellis (semi-disk), B2ST, DiGeST, WaveFront, ERa (Mansour et al., VLDB 2011), PCF | O(n^2) worst case for all disk-based ones V (ERa Table 2); B2ST O(cn), c = 2n/M V | Trellis needs the string in RAM V; B2ST temporary results 343 GB for 2.6x10^9 symbols = 132 B/char D (ERa section 3) |
| C8 | Sample-splitter parallel sorts: CaPS-SA (Khan, Rubel, Molloy, Dhulipala, Patro, WABI 2023 / AMB 2024); PSAC (Flick, Aluru, SC 2015); Fischer, Kurpicz (ALENEX 2019) | CaPS-SA O(n log n + sum L_i) comparisons w.h.p., worst case O(n^2) V; PSAC O(T_sort(n,p) log n) U | CaPS-SA working space 4 w n bytes, w = 4 or 8 (16 or 32 B/char), external memory "left to future work" V |
| C9 | Distributed DCX with bucketing (Haag, Kurpicz, Sanders, Schimek, arXiv 2024 / ESA 2025) | DCX with q buckets so that only n/(pq) X-tuples are materialised V; section 5 sketches the distributed external-memory variant (scan the text once per bucket, splitters from a sample) V | 14-15x input (large), 20-26x (medium) V; no EM implementation |
| C10 | BCR, extLCP, eGSA (Bauer, Cox, Rosone 2013; Cox, Garofalo, Rosone, Sciortino 2016) | scan-only, disk <= 2x output (extLCP) V (abstract) | RAM independent of the number of strings V (BCR abstract, via search record U); **string collections, not one long text** |
| C11 | Prefix-free parsing / Big-BWT (Boucher et al. 2018/2019) | U (search record only) | grows with dictionary + parse; designed for repetitive collections, not one human genome |

Newer than 2019 that I found: CaPS-SA (C8, 2023/24), Haag et al. (C9, 2024/25), the journal version of the EM LCP algorithms (JEA 2019, same content as ESA 2016, U), libsais 2.8-2.9 (64-bit, OpenMP, generalized SA + LCP; in-RAM only; Apache-2.0; U search record). I found **no 2019-2026 external-memory single-text SA+LCP construction paper that supersedes pSAscan / fSAIS / eSAIS**; two 2020 papers (nSAIS "Succinct suffix sorting in external memory", "An Improved Algorithm for Building Suffix Array in External Memory") are known to me only as listings that I could not open (U).

Table 2.2 numbers reported by the authors.

| id | result | where |
|---|---|---|
| C1 | "inputs of up to 4 GByte in hours" and 4x10^9 characters "over night" on a 2.0 GHz Xeon pair, 1 GB RAM, four 80 GB disks V; hg19-era Genome n = 3,070,128,194, maxlcp 21,999,999 V; DC3 6.88 us per input byte on Genome (2.86 GB) V; I/O 200-300 B per input byte U (FGM citing DKMS) | DKMS05 abstract, Table 1, conclusion; FGM Fig. 2 |
| C2 | 80 GiB Wikipedia, 4 GiB RAM: SA 2.4 days, 18 TiB I/O = 2.41 us/B and 230 B/char D; SA+LCP 5.0 days, 35 TiB = 5.03 us/B and 448 B/char D (all V); about 2x faster and 60 % of the I/O of DC3 on real inputs V | thesis section 6.3.1-6.3.2 |
| C3 | 2x faster than eSAIS, 1/3 of the disk, fastest EM SA algorithm "for inputs about 75 times larger than the available RAM", pSAscan is faster for smaller inputs; no LCP V | thesis section 6.4.1 (the ALENEX 2017 paper itself could not be opened) |
| C4 | 200 GiB with 3.5 GiB RAM about 170 h = 2.9 us/B D; 1 TiB with 120 GiB RAM a little over a week, 7.2 TiB disk V (README); **512 GiB of 1000 Genomes reads, 120 GiB RAM: pSAscan 2.51 days, 16.63 TiB I/O (33 B/char D), 3.75 TiB disk; LCP given SA: EM-SI 2.16 d / 25.89 TiB (52 B/char D), EM-SPhi 2.69 d / 20.77 TiB (42 B/char D), LCPscan 5.04 d / 62.22 TiB** V | README; ESA 2016 Table 4 (data are reads, not an assembled genome) |
| C5 | 1 GB random DNA, 2.6 GHz Pentium 4, 0.90 GB memory: dnabwt (2-bit text) 8,250 s = 8.25 us/char CPU; **repeat-DNA (random block of 1024 repeated) 52,668 s = 6.4x**; deep-shallow (Manzini-Ferragina) on repeat-DNA 256 MB: 70,125 s vs 237 s random = 296x V | TCS 2007 Table 1 |
| C6 | Genome 2.86 GB, block m = 400 MB, one core of a 2.5 GHz Phenom: bwt-disk 3.72 us/B total (step 3 = 2.32 us/B), DC3 6.88 us/B, bwt-disk I/O < 6x input (14.76x on Random2); "designed to be lightweight, not to be fast" V | arXiv 0909.4341 Fig. 2 |
| C7 | ERa: human genome (2.6 Gbp) in 19 min on an 8-core desktop with 16 GB RAM; PWaveFront 15 min on 1024 CPUs of an IBM BlueGene/L (the authors' claim) V; PCF (parallel B2ST): 3 GB in 425 s on 172 processes V; B2ST 12 GB DNA in 26 h on one machine with 2 GB RAM U | ERa abstract; PMC3962650; search record |
| C8 | CaPS-SA, 32 threads: GRCh38 93 s / 45 GB, T2T 106 s / 48 GB (divsufsort 556 s / 25 GB); 1 thread GRCh38 35 min vs divsufsort 9 min; pathological Bacteria(1K) (mean LCP 10,897): 259 s vs parallel-dc3 164 s, 1 thread 5,294 min vs divsufsort 726 min V. PSAC: human genome in < 8 s on 1024 cores U | Table 2, Table 4 |
| C9 | up to 5x faster than prior distributed algorithms on up to 12,288 cores; on a repetitive input (first MB of DNA repeated 245,760 times) prefix doubling "did not terminate or crashed" while DCX (X = 39) was insensitive V | ESA 2025 html |

Table 2.3 properties that decide the choice.

| id | LCP? | output in final order, streamable? | repetitive text (a^n, satellites, N runs) | alphabet | implementation, licence | effort in managed C# (mmap, unsafe) |
|---|---|---|---|---|---|---|
| C1 | DC3-LCP exists (disk 88n/109n), "justifiable time only for very small instances" V | last merge yes; 5-10 sort passes before | "almost completely insensitive to the input" V | any integer | Dementiev esuffix on STXXL (licence not checked U) | very high (pipelined external sorters) |
| C2 | yes, first fully external LCP V | final PQ extraction yes | recursion depth worst case ("Skyline") is slower but bounded V | 8/32/40/64-bit templates V | GPL v3, STXXL V | very high (external priority queue) |
| C3 | no V | yes | bounded (worst-case I/O bound) | configurable char and offset type V | MIT/X11 V | high (radix heap, pre-inducing) |
| C4 | no; separate EM-SPhi / EM-SI pass | **only the final merge**; partial SA is rewritten once per block | data-independent cost model; the LCP phase depends on irreducible LCPs (de Bruijn is its worst case) V | bytes only, byte 255 reserved V | MIT/X11, needs libdivsufsort 2.0.1 V | moderate |
| C5 | no (BWT only); LCP not produced | **yes: blocks come out in lexicographic order** V | DC bound O(n log n + v n); 6.4x measured on periodic DNA V | alphabet independent; 2-bit program exists V | Bowtie 2 class KarkkainenBlockwiseSA, GPL v3 or later, threaded, block temp files V (Bowtie 1 has the same file, U) | moderate (+ DCS recursion) |
| C6 | no (BWT; SA/Psi variants) | no (partial BWT rewritten per pass) | data-independent passes | bytes | bwt-disk page not opened (licence U) | low |
| C7 | tree, not SA/LCP | per S-prefix sub-tree in order | O(n^2) V | small | ERa code not found (U); WaveFront patent US8914415B2 U | high |
| C8 | yes (LCP-aware merge) | yes in RAM; EM not designed | worst case O(n^2) V | bytes | CaPS-SA: README states no licence, the paper says "open-source license" (LICENSE file not opened) V | low in RAM, EM undefined |
| C9 | no | yes (sketched) | robust (DCX) V | any | licence not stated in the paper (repository not opened) V | very high (distributed) |
| C10 | yes, for collections | yes | n/a | DNA reads | various | not applicable to one text: as I understand BCR (one symbol of every string per round) it would need n rounds for a single string of length n (E, the full paper was not opened) |

## 3. Requirement-by-requirement scoring

Scale: ++ meets, + partly or with caveat, - fails, ? unknown. The B-DCS column is a **design assessment before any measurement** (E).

| req | C1 pDC3 | C2 eSAIS | C3 fSAIS | C4 pSAscan | **C5 B-DCS (rec.)** | C6 BWTDisk | C7 ERa class | C8 CaPS-SA class |
|---|---|---|---|---|---|---|---|---|
| R1 RAM independent of n | ++ | + (RMQ 300 MiB) | ++ | ++ | + (RAM fixed; DCS and text via page cache) | ++ | + (sub-tree and string windows) | - (16-32 B/char) |
| R2 streamed sequential output in final order | + | + | + | - (final pass only, k rewrites) | ++ | - | + (tree-shaped) | + (RAM only) |
| R3 worst-case robustness | ++ | + | ++ | ++ (SA) / + (LCP) | + (published bound; EM adaptation untested) | ++ | - | - |
| R4 alphabets incl. 2-bit DNA | ++ | ++ | + | - (bytes) | ++ (key = 64 bits of packed symbols) | + | + | + |
| R5 n >= 3x10^9, 64-bit | + | ++ (40-bit V) | ++ | ++ (<= 1 TiB V) | ++ | + | + | + (LARGE_IDX) |
| R6 managed C# | -- | -- | - | + | + | ++ | - | + |
| R7 parallel | + | - (1 core V) | + | ++ | ++ (independent buckets) | - (1 CPU V) | ++ | ++ |
| R8 I/O and 3 Gbp wall time | table below | table below | table below | table below | table below | table below | - | n/a |
| R9 LCP and child table in one stream | + (very slow) | + | - | - | ++ (design: fused Phi pass + scan) | - | + | + (LCP only) |
| R10 navigation pack buildable in EM | ISA via sort | ISA via sort | ISA via sort | ISA via sort | ++ (Psi in one pass; ISA by 2-pass distribution) | + (Psi, quadratic V) | - | - |

### 3.1 R8: I/O volume per input character and expected wall time for 3x10^9 characters

| id | I/O bytes per input char | wall time for 3 Gbp | basis |
|---|---|---|---|
| C1 pDC3 | 200-300 U; sort(30n)+scan(6n) with 4-5 B words D | 5.7 h on the FGM machine (6.88 us/B V, D); >= 2 h on a current CPU (E) | FGM Fig. 2 |
| C2 eSAIS | 230 (SA), 448 (SA+LCP) D | 2.0 h (SA), 4.2 h (SA+LCP) D at 2.41 / 5.03 us/B (4 GiB RAM, 2012 hardware); more with M = 256 MB (E) | thesis p. 226 |
| C3 fSAIS | about half of eSAIS (E from "twice as fast") | about 1 h SA only (E); LCP extra | thesis section 6.4.1 |
| C4 pSAscan + EM-SPhi | SA 33 at n/M = 4.3 D, scaled by n/M to about 90 at n/M = 12 (E); LCP 42 D | SA 0.39 us/B at 12 cores D, about 2.8x more per byte at n/M = 12 => 1.5-3 h on 4 fast cores; LCP about the same again => 3-6 h (E) | ESA 2016 Table 4 |
| C5 B-DCS | **about 90 B/char** D (stage sums in 4.2), 270 GB for 3 Gbp: 135 s at 2 GB/s, 540 s at 500 MB/s, so I/O is not the bottleneck | **25-90 min** on 4 vCPU (E); compare dnabwt 8.25 us/char CPU on a 2007 core V, about 20 min if a current core is 6x faster and 4 threads scale (E) | section 4.2 |
| C6 BWTDisk | < 6x input V (BWT only) | 3.72 us/B at n/m = 7 V => >= 10 h at n/m = 30-100 (E) | FGM Fig. 2 |
| C7 ERa | tree files 26x input (Trellis, V ERa section 6) | 19 min on 8 cores and 16 GB RAM V (RAM 60x ours); not comparable | ERa |

### 3.2 R9: the lcp-interval scan, its stack and a^n

Bottom-up scan (Abouelhoda-Kurtz-Ohlebusch style): push when `LCP[i] > top.l`, pop while `LCP[i] < top.l`.
* **Claim K1.** The stack holds the open intervals containing the cursor with strictly increasing lcp values, so its height is at most `1 + max LCP`, and at most n+1. A random text over sigma has `max LCP` about 2 log_sigma n (about 30 for DNA at 3x10^9, E); the lcp-interval tree depth is the stack height.
* **a^n**: the augmented SA is `$, a, aa, ...` and `LCP = 0, 0, 1, 2, ..., n-1`: every entry pushes, nothing pops before the end, the height is **about n** (the ESA of a^n is a chain of n-1 internal nodes below the root, ADR-001 section 4). Real inputs reach the same order: the hg19-era Genome file has maxlcp 21,999,999 (V DKMS05 Table 1; long N runs), i.e. 22 M stack entries = 264-352 MB at 12-16 B per entry, **more than M**. So a RAM stack violates R1 and must spill.
* **Design.** Stack entry = `(l, lb, lastChildIndex)` (12-16 B). The top two blocks of 256 K entries are RAM-resident; when the top block fills, the older block is appended to `stack.spill`; when a pop leaves the window, the last spilled block is read back (backwards, sequentially). At least one full block stays hot, so push/pop alternation at a block boundary does not thrash. Amortised I/O is O(1/block) per operation, the file is accessed only at its end (LIFO). Worst-case spill volume is 12-16 B x (max LCP + 1) = 36-48 GB for a^n at 3x10^9, known before the scan from the PLCP maximum, so the build can fail early with a disk-space error. Optional later optimisation: run-length encode arithmetic progressions `(l0, lb0, count)` (a^n and N runs become O(1)); not needed for correctness.
* **Child-table writes behind the cursor.** The Abouelhoda layout writes `up/down/nextlIndex` at the index of an older interval (`down[top]`, `nextlIndex[last]`). All targets are on the stack, pops happen in LIFO order, each old interval is popped exactly once, so the dirty pages form at most one forward hot region (cursor) and one backward sweep; the number of such writes is at most 2n. This is cache-friendly on a mapped file but **not proven strictly bounded**, so E5 measures dirty pages and write amplification on a^n, Fibonacci and an N run. **Fallback layout F**: store `PSV[i]` (previous smaller value) written forward and `NSV[i]` (next smaller value) written by a **backward** scan of the LCP stream with the same stack; then every write is sequential (forward or backward) at the price of reading the LCP stream twice (1.1 B/char). ADR-001 prototype P3 already evaluates an RMQ/PSV/NSV-only variant; ADR-004 adds the requirement "all CHILD writes must be sequential in one direction or LIFO-local".
* The scan is single-threaded and cheap (it reads 1.1 B/char of LCP and writes 4 B/char of CHILD; E: 3-10 ns per entry, i.e. 10-30 s for 3x10^9 entries). It runs fused with PH4-B3, so it adds no I/O pass of its own for the LCP.

### 3.3 R10: the navigation pack and the inverse permutation

* Need: suffix link `LCA(Psi(lb), Psi(rb))` with `Psi(i) = ISA[SA[i]+1]` (ADR-001 section 5). ISA is the inverse permutation of SA.
* **External inverse permutation (ISA)**: it is an external sort of n pairs `(SA[i], i)` by the first field, but because the keys are a dense permutation of [0, n) a **one-level radix distribution** suffices: pass 1 reads SA and writes the tuples into P = n / R chunk files by `SA[i] / R` (R positions fit in RAM, P <= ~400 so that P buffers of 64-256 KB fit); pass 2 reads each chunk file, scatters it into an in-RAM array of R entries and appends it to the ISA section. I/O = read SA (w n) + write and read the tuples (2 x 2w n) + write ISA (w n) = 6w n = 24 B/char for w = 4, scratch 8n (D). This "distribute, then scatter" primitive (called DIST below) is also the building block of the chunked Phi/PLCP LCP and of the sample-rank inversion in the DCS recursion.
* **Psi without ISA** (cheaper and fully streaming): with `BWT[j] = T[SA[j]-1]` and `C[c]` = number of suffixes with first symbol smaller than c, `Psi[C[c] + k]` = rank of the k-th (0-based) occurrence of c in BWT (Psi is the inverse of LF and is increasing inside every first-symbol bucket; checked here, D). One scan of BWT in rank order with sigma output cursors writes Psi sequentially in sigma places: DNA sigma <= 5, bytes sigma <= 256 (256 windows of 64 KB), larger alphabets in ceil(sigma/64) passes. I/O = 1 B/char (BWT, emitted for free when the SA slice is produced) + w B/char (Psi) = 5 B/char instead of 24. The pack therefore holds **PSI + RMQ/PSV/NSV** (ADR-003 kind 8 `ISA` is replaced by a new kind `PSI`; `ISA` stays optional for `RankOf`/`ExtendLeft`, built on request by DIST). The FGM paper also computes Psi in EM by scans (Thm 4.3 V) but at O(n^2/(M B log n)) I/O.
* **Can the pack be avoided?** Not for the matching-statistics family (MEM, MUM, LCS, anchors): the alternatives are (a) a stored suffix-link table, 4-8 B per internal node = 2.5-5 B/char on DNA (0.62 internal nodes per char, BENCHMARK.md via ADR-001) and needs the same Psi lookups at build time (about 1.9x10^9 nodes x 2 random reads x 100 ns = 380 s single-thread, E): it is not smaller than Psi+RMQ (about 4.4 B/char) and costs an extra pass; (b) O3 re-descent, quadratic on a^n. **Consequence for ADR-001 section 5 R-SL(1)**: the "lazy in-heap pack on first use" cannot exist for out-of-core sizes (4.4 B/char = 13 GB of heap at 3x10^9). Rule: above a `PackHeapLimit` the family is `Unsupported` unless the pack section was persisted at creation (profile "Full").
* Related: the semi-external LCP algorithm (SE-SI) needs `PLCPsucc[0..2n)` in RAM, i.e. n/4 bytes = 750 MB at 3x10^9 V (ESA 2016 section 6), hence the fully external route.

## 4. Recommendation

### 4.1 Why B-DCS

It is the only family whose natural output is the one the owner asked for: **independent tasks produce consecutive slices of the final SA**, so the SA, the LCP and the child table are appended in rank order with no merge phase and linear scratch I/O, and the buckets are embarrassingly parallel. It has a published worst-case bound that does not depend on repetitiveness (DCS: two suffixes with a common prefix of v-1 symbols are ordered in constant time, V TCS 2007 Lemma 7; in my reading this is two sample-rank lookups at an offset l < v, D), and a measured benign behaviour (6.4x on periodic DNA while the deep-shallow sorter is 296x slower, both V). It is alphabet independent and has a 2-bit program (V). Its weakness is that the published implementation keeps the DCS in RAM (O(n / sqrt(v)) bits); in the EM version the DCS is a file. The alternatives are worse on the requirements that are categorical: pSAscan rewrites partial SAs and has no LCP stream, eSAIS/fSAIS need an external priority queue (eSAIS measured 230-450 B/char of I/O) and a much larger engineering effort, CaPS-SA is RAM-only and quadratic on repeats, ERa-style trees do not produce SA/LCP.

### 4.2 Data flow (stages, files, buffers)

Notation: `w` = bytes of a stored position (4 if n+1 <= 2^32-1, 5 in scratch or 8 above, section 4.3); `K` = records per bucket task; `X` = difference-cover period (power of two, from a precomputed cover table with |D| about sqrt(1.5 X); the table is not verified here); sample fraction `f = |D|/X`, 0.04-0.05 for X = 1024 (E). Numbers below are for n = 3x10^9, DNA, w = 4, 4 threads.

| stage | what it does | reads | writes | RAM | parallel |
|---|---|---|---|---|---|
| PH0 prepare | stream-copy the input into the TEXT section (codec 4, 2-bit or b-bit packed; ADR-003 reserves codec 4 for S2, pull it forward), alphabet scan, symbol counts C[], SHA-256 of the text | input (1 B/char) | TEXT 0.25 B/char | 4 MB | 1 |
| PH1 DCS | sample positions `i mod X in D`; sort the sample by its X-prefix (buckets of the sample: per bucket one scan of the sample positions, keep records whose X-prefix is in range, sort, name equal prefixes); names go to `names.L0` (position order, via DIST or a q-way merge by sample index); if names are not unique recurse on `names.L0` as a text over a 32-bit alphabet with its own X_l; the base case (<= N_base = 16 M symbols) runs the existing `SaisKernel` on the int32 reduced string and inverts the SA in RAM; result `rank.L0[s]` = rank of sample suffix s | text, sample | names 0.17 B/char, rank 0.17 B/char (scratch), levels >= 1 are < 5 % | 150-200 MB | by sample bucket |
| PH2 partition | oversampled random positions (seeded) sorted by the **DCS comparator** give 64 super-bucket splitters; one scan of all n positions binary-searches the splitters (key of 32 symbols first, DCS comparator only on key ties) and appends the position (5 B) to one of 64 super-bucket files; each super-bucket file is split again in a second pass into buckets of <= K using splitters from a reservoir sample of the file; oversized buckets are split once more (Bowtie rebuilds samples up to 5 times V; Karkkainen section 7 has a deterministic procedure with block sizes between b_max/2 and b_max V) | text, 5n | 2 x 5n (scratch, 5n at a time) | 64 x 128 KB x 4 threads = 32 MB | by text range |
| PH3 bucket sort | per bucket: read positions (increasing text order, so the key reads are monotone), key = first 32 symbols as a 64-bit word, LSD radix sort of (key, position), runs of equal keys are sorted by the **DCS comparator** (compare l < X symbols, then `rank[i+l]` vs `rank[j+l]`); write the SA slice at offset `prefix(bucket sizes)` through a mapped window of <= 4 MB; emit BWT bytes (scratch) and the first-k-mer boundaries of the BUCKET section from the sorted keys | position file, text, rank file (random) | SA 4 B/char at final offsets; BWT 1 B/char | 4 tasks x K x 32 B (K = 1 M) = 128 MB + 4 x (1 MB in + 4 MB out) | **4 independent tasks**, completion watermark for downstream |
| PH4 LCP (chunked Phi/PLCP) | B1: scan SA in rank order, tuple `(SA[i], SA[i-1], i)` (12 B) distributed by `SA[i]/R` into P = n/R files (R = 8 M positions, P = 375); B2 per chunk in text order: scatter Phi, run the PLCP loop of `LcpKernel.ComputePlcp` (carry l, vector compare, <= 2n compares in total), count escapes (>= 255, >= 65535) -> **choose U8/U16/U32 now, before any LCP entry is written** (ADR-002 section 4 rule); emit `(i, PLCP[j])` distributed by `i/R` (8 B); B3 per rank chunk: scatter into the compact LCP, append to the LCP section and the overflow list, and **feed the same chunk to PH5** | SA, text, tuples | tuples 12n + 8n (scratch, bounded mode: text-range rounds with S_max), PLCP compact, LCP 1.1 B/char | 375 x 64 KB + Phi chunk 32 MB + 16 MB = 100 MB | by chunk (carry restart costs <= max PLCP per chunk, 3x10^9 SIMD bytes for a^n) |
| PH5 scan | one pass over the LCP stream in rank order: child table with the external stack (3.2), NODE_COUNT (`1 + I + n + 1`), LRS_DEPTH/LRS_RANK, NEXTSYM (`text[SA[i]+LCP[i]]`, random text reads) | LCP, SA, text | CHILD 4 B/char, NEXTSYM 1 B/char | 2 x 4 MB stack window + 3 writer windows | 1 (sequential by definition) |
| PH6 pack (optional) | PSI by sigma-way distribution of the BWT scratch; RMQ: block minima + sparse table over blocks from one LCP scan; PSV (forward) and NSV (backward) if layout F | BWT, LCP | PSI 4 B/char, RMQ about 0.4 B/char | sigma windows x 1 MB | PSI by symbol group |
| PH7 seal | header, section table, CRC-32C of sections (one extra read pass at ~10 GB/s), `fsync`, magic written last, `rename` over the target, delete scratch | all | header | 4 MB | - |

I/O sum (D): PH0 1.25 + PH1 0.6 + PH2 20 (position files written and read twice) + PH3 5 (SA 4 + BWT 1) + PH4 47 (SA read 4 + Phi tuples 24 + PLCP 2.2 + rank tuples 16 + LCP 1.1) + PH5 6 (LCP 1.1 + CHILD 4 + NEXTSYM 1) + PH6 8.5 (BWT 1 + PSI 4 + RMQ 0.4 + LCP 1.1 + PSV/NSV 2) = **about 90 B/char**, i.e. 270 GB for 3x10^9. Final file (Core, 2-bit text): TEXT 0.25 + SA 4 + LCP 1.1 + CHILD 4 + BUCKET 0.1 + NEXTSYM 1 = 10.5 B/char = 31 GB (+ pack 4.4 B/char = 44 GB Full). Peak scratch: max(5n positions, 12n Phi tuples) = 36 GB at 3x10^9. **Bounded-scratch mode** repeats B1/B2 in text-range rounds of rho positions (each round re-reads the SA sequentially), keeps the compact PLCP in text order (1.1 B/char), and replaces the second distribution by a gather `PLCP[SA[i]]` from that file (random reads of a 1.1 B/char file through the page cache, as in pass 3 of ADR-002 section 4); scratch = 12 rho + 1.1 n, so S_max = 4 GB at 3x10^9 gives about 9 rounds.

**RAM budget (M = 256 MB), worst stage PH3**: runtime baseline 50 MB (U) + 4 tasks x 32 MB bucket buffers + 4 x 5 MB I/O windows + splitter table 64 KB + stacks 4 x 1 MB = about 200 MB; margin 56 MB. `M(n) = c1 K + c2 F B_buf + c3 (n/K) s_split` with fan-out F <= 400 and a splitter table of 16 B per bucket (48 KB at 3x10^9); the only n-dependent term is the splitter table. Larger n changes the number of buckets and scans, not the buffers; a second-level fan-out keeps F bounded (a sqrt(n) dependence of M is also the published requirement of EM-SA-DS, V).

### 4.3 Answers to the specific design questions

* **How the SESA sections are appended.** All sizes except the LCP representation are known after PH0 (n, sigma, w): TEXT, ALPHABET, SA, CHILD, BUCKET, NEXTSYM, PSI are placed first at offsets computed up front; the section table is written last (ADR-003 section 6 already writes the magic last). `LCP` and `LCP_OVERFLOW` are placed after them because their sizes are known only after PH4-B2 (before the first LCP entry is written), then RMQ, BUILD_INFO. The section order of ADR-003 section 3 is a writer convention, not a format rule (offsets are in the table), so readers are unaffected. SA slices are written by up to four tasks at disjoint offsets, each sequentially; LCP, CHILD, NEXTSYM, PSI by one sequential writer each. Writers use a **sliding mapped window** (<= 4-64 MB `MemoryMappedViewAccessor`, `Flush` and dispose when full), so the process never holds more than a few windows of dirty file pages; if dirty-page behaviour of mapped windows is poor on a platform, the same writer can use `RandomAccess.Write` (BCL pwrite) on the same file with identical bytes.
* **LCP between partition boundaries.** Not needed in the default path: the Phi/PLCP stage reads the SA as one stream (`SA[i-1]` is simply the previous entry, also across bucket boundaries). An optional in-bucket fast path would need `LCE(last of bucket b-1, first of bucket b)` per boundary, which is O(max LCP) per boundary (q x n for a^n), and the DCS gives order but not LCP values; therefore it is **not** the default.
* **How repeats are handled.** Order: DCS comparator, at most X-1 symbol compares and two rank lookups for any pair, independent of the LCP; N runs, satellites and a^n only make many pairs reach the comparator (cost model: ties cost 2 random reads in the rank file, about 160 ns, E). LCP: Phi/PLCP, at most 2n symbol compares for any text (PLCP[j+1] >= PLCP[j]-1). Stack: external. Bucket capacity: splitters are real suffixes ordered by the comparator, so even a^n is split into equal buckets (the suffixes `a^k` are distinct).
* **How int64 positions are stored.** Records and files use a width chosen from n at PH0: `w = 4` (unsigned, n+1 <= 4,294,967,295 covers the human genome of 3.1x10^9, whereas ADR-002's int32 kernels stop at 2,147,483,591), `w = 5` (40-bit) in scratch for 2^32 <= n < 2^40 (the EM tools use 40-bit V), `w = 8` above. SESA `POS_WIDTH` stays 4 or 8 (ADR-003 section 3); an optional width 5 saves 37 % of SA/CHILD/PSI for both strands of a mammal (6.2x10^9). In code positions are `ulong`, chunk-local indices are `int`; `Buf<T>` (int-indexed) is used only on chunks of <= 2^31 elements, file offsets are `long`. The int32 kernel is called only on reduced strings of <= N_base symbols, so **no 64-bit SA-IS kernel and no generated `long` clone is needed for the persistent builder** (ADR-002 section 2 item 6 and the S2 "64-bit path" stay relevant only to the in-memory backend, which cannot exceed ~1.07-2.1 Gsymbols anyway).
* **Crash safety and atomicity.** Output is created as `<path>.tmp` (first 8 bytes zero) and scratch under `<path>.tmp.d/`, on the same volume as the target. A small `state.json` records the input fingerprint (length, SHA-256 from PH0), the parameters (X, K, w, S_max) and the completed stages; a stage is marked done only after its outputs are `fsync`ed; every stage is idempotent (SA slices go to fixed offsets, scratch is recreated), so a resume skips completed stages (`Resume = true`) and a mismatching fingerprint discards the scratch. Publication: header with CRCs, `fsync`, magic last, `fsync`, `rename`. A crash leaves either no target or a `.tmp` file without magic, which the dispatcher rejects (ADR-003 section 6). Pre-flight: free space >= final size + S_max + stack-spill bound (computed after PH4-B2, before PH5); failure is early and loses no data.
* **Verification of the result.** The existing `Structure` validation of ADR-003 section 6 is O(n) with random reads and is run by the prototype on every test output; production uses `Header` plus per-section CRCs.

### 4.4 Reuse of what exists, and what must be added

| existing | use in the EM builder |
|---|---|
| `SaisKernel` (int32), `PreparedText`, `SymbolReaders` | base case of the DCS recursion (names string as int symbols), the n <= N_base path, and the SA-IS block sorter of the fallback |
| `LcpKernel.ComputePlcp<TSym>` (PLCP loop, vector compare, prefetch), `CompactLcp` selection rule | PH4-B2 PLCP loop on a Phi chunk (the loop must be factored out of the function that builds Phi from an in-RAM SA: new entry `PlcpScan(text, phiChunk, firstIndex, carry)`); the selection rule of ADR-002 section 4 decides U8/U16/U32 after B2 |
| `NativeArena`, `Buf<T>` (int-indexed, `NativeArena` not thread-safe) | bucket buffers, Phi chunks: **one arena per task**; add an `ArenaBudget` counter that throws when the sum of live blocks exceeds the cap (enforces R1 from inside) |
| generic algorithms (`SuffixTreeAlgorithms.*`), `EsaNavigator` | unchanged: they read the mapped SESA |
| new | `MappedWindowWriter`, `MappedTextReader` (long-indexed, 2-bit/b-bit packed, `ulong` key extraction), `Dist` (distribute-and-scatter primitive), `SpillStack`, `DcsTable` (difference-cover table and `l(a,b)` offsets), comparator, radix sorter for (key, position), stage runner with journal, I/O and RSS counters |

### 4.5 Which parts of ADR-002 / ADR-003 / ADR-001 / PROGRESS must change

| document | section | change |
|---|---|---|
| ADR-002 | header, section 3 (workspace plan), last sentence "the SA is written section by section from these buffers" | applies to the **in-memory builder only**; the persistent builder is ADR-004 with M <= 256 MB. The 12-17 B/char table stays valid for the in-memory index |
| ADR-002 | section 2 item 7 and "Deferred" (S2: 64-bit SA path, external-memory build) | the 64-bit SA-IS path is **not needed** for the persistent index (4.3); external-memory build is promoted from S2 to an S1b blocker |
| ADR-002 | section 4 (LCP) | Phi/PLCP on an in-RAM `int[]` stays for the in-memory builder; the persistent builder uses the chunked DIST variant; the U8/U16/U32 rule is unchanged but evaluated after B2 |
| ADR-003 | section 6 last sentence ("builds sections directly from the build buffers; no second in-memory copy") and section 3 (order of sections written) | replaced by the staged streaming writer of 4.3; the writer assigns offsets up front, the table is written last; CRC-32C per section computed in a final pass or incrementally by the sequential writers |
| ADR-003 | section 4 kind 8 `ISA` | add a kind `PSI` (id from the reserved range 12-4095, proposal 12; `PSI[r] = ISA[SA[r]+1]`, `PSI[0]` = 0 by convention); the pack is PSI + RMQ; `ISA` becomes optional on request; ADR-003 section 5 size table unchanged (PSI is 4 B/char like ISA) |
| ADR-003 | section 4 kind 1 TEXT codec 4 | 2-bit (b-bit) packed text must exist in S1b, not S2 |
| ADR-003 | section 3 `POS_WIDTH` | optional width 5 (40-bit) for 2^32 <= n < 2^40; width 4 is unsigned and already covers 3.1x10^9 |
| ADR-003 | section 9 (mapping and lifetime) | add the write side: sliding windows, no `madvise` (still no P/Invoke) |
| ADR-001 | section 5 R-SL(1) | the lazy in-heap pack cannot exist above `PackHeapLimit`; the pack is persisted at creation (profile "Full") or the matching-statistics family is `Unsupported` there |
| PROGRESS-S1 | S1b work packages | 1b.7 (container writer) is replaced by: 1b.7a `Dist` + `SpillStack` + counters, 1b.7b `MappedWindowWriter`/`MappedTextReader` + codec 4, 1b.7c PH4-PH5 (LCP, child) on a given SA file, 1b.7d PH1-PH3 (DCS, partition, buckets), 1b.7e PH6-PH7 (pack, seal, resume); the prototype gates of section 5 precede 1b.7c-d |

### 4.6 Fallback and switch rules

* **F1 = SAscan family** (C4), clean-room from the paper (the pSAscan code is MIT, read-only reference): blocks of ~25 M symbols (about 10 B/char with the SA-IS kernel, E), gap array by streaming the text through the block's BWT, merge with the partial SA; reuses `SaisKernel` for the block SAs; then LCP by the PH4 stage (reads the SA file). Costs: I/O about 90 B/char at n/M = 12 (E), 1.5-3 h for 3 Gbp (E). Choose F1 when a gate of section 5 fails.
* **F2 = fSAIS class** (C3) only when n/M > ~75 (published crossover V): at M = 256 MB that is n > ~20 GB; not planned for S1.
* **F0 = n <= N_base**: one bucket, the existing in-RAM SA-IS + LCP + the PH5 writer; working memory is bounded by N_base whatever n is (Q6).

## 5. Prototype plan with measurable gates

Rules for all experiments: no third-party binary is run; oracles are `SaisReference`/`LcpReference` and the committed libsais hashes (`analysis/sa_vectors.tsv`: `dna_1m`, `dna_10m`, `rep_10m`, `prot_1m`, adversarial families at 10^5 and 10^6); inputs from `analysis/bench/gen/gen_inputs.py` (present locally in `/tmp/stbench/data`: `dna_1m`, `dna_10m`, `dna_20m`, `dna_50m`, `rep_10m`, `prot_1m`; `gen_extra_100m.py` makes `dna_100m`, which is not present yet); larger and genome-like inputs come from a C# generator in the prototype (the Python generator builds the text in memory and is slow above 100 M; planned families: random + interspersed repeats with 10-15 % divergence + 171-bp satellite arrays + N runs). Disk: 22 GB free means **500 Mbp is the largest local run** (output about 5 GB + scratch about 12n = 6 GB); a 1 Gbp Core build fits only in bounded-scratch mode (S_max = 4 GB: output 10.5 GB + scratch about 5 GB); 3 Gbp needs about 31 GB output + bounded scratch + stack spill, i.e. about 50-60 GB (Q8).

| id | experiment (smallest first) | inputs | gate | abandon or change the design if |
|---|---|---|---|---|
| E0 | harness only: RssAnon sampler (`/proc/self/status` every 10 ms; BENCHMARK.md already samples RssAnon every 20 ms V), I/O counting wrappers, `ArenaBudget` | no-op stages | sampler overhead < 2 %, counters match `/proc/self/io` within 1 % | - |
| **E1** | **in-RAM simulation of PH2-PH3 without DCS** (direct comparison with a vectorised LCE), bucket capacity scaled to n/1000 (the 3x10^9 / 3 M ratio), SA only; a run that exceeds 10 minutes counts as a failure | `dna_10m`, `rep_10m` (moderate repeats only: 100 copies of a 100 kbp block with 1 % substitutions), `prot_1m`, a^n (1 M, 10 M), Fibonacci, Thue-Morse, 171-bp tandem x 10^4, N run 1 M (the adversarial families come from the 1a.3/1a.7 test generators) | SA equal to the libsais hashes; dna_10m <= 3 s on 1 thread; report the share of positions in key ties and the sum of LCE work per family | **if every family is within 10x of dna_10m, DCS is not needed and PH1 is dropped** (large simplification); if not, E2 |
| E2 | add DCS (in-RAM ranks, X = 512/1024/2048) to E1 | same | every adversarial family <= 10x dna_10m; dna_10m slowdown from DCS < 25 % | a^n or tandem > 25x after tuning X and K: use the F1 fallback |
| E3 | scale-down emulation of the memory ratio: n = 100 M with K = n / 3000, rank file and text read through mappings, then repeated with the page cache limited to 0.5x (rank + text). Limiting method: a cgroup v1 memory limit if the prototype user may create a sub-cgroup (the memory controller is mounted here; `touch` in the cgroup root was denied, `mkdir` is **untested**), otherwise a helper process that pins anonymous memory so that only 0.5x (rank + text) is left for the cache | `dna_100m` | RssAnon <= 256 MB; slowdown under the limited cache <= 3x | > 5x: the design needs a minimum page cache; document it or drop to F1 |
| E4 | PH4 (chunked Phi/PLCP via DIST) on a given SA file | `dna_10m`, `rep_10m`, 100 M, a^n 50 M | LCP equal to `LcpKernel`; I/O <= 50 B/char; throughput >= 5 MB/s on 4 threads | LCP phase > 3x the SA phase |
| E5 | PH5 with the external stack on a^n (50 M), N run 30 M, Fibonacci; both child layouts | as listed | child table equal to the reference; spill volume equals 12-16 B x (max LCP + 1); mapped dirty pages <= 64 MB sampled from `/proc/meminfo` Dirty; write amplification < 2 | layout A thrashes: use layout F |
| E6 | end to end SA+LCP+child on 100 Mbp (dna and genome-like), 4 threads | `dna_100m`, genome-like 100 M | **total <= 90 s (>= 1.1 MB/s) and SA stage <= 30 s**; I/O <= 100 B/char; RssAnon <= 256 MB incl. runtime; output identical for 1/2/4 threads and for two K values | total > 4.5 min (extrapolates to > 2.3 h for 3 Gbp) |
| E7 | scaling: 50 M, 100 M, 200 M, 500 M | same family | fitted exponent b in `t = a n^b` <= 1.15; RssAnon slope vs n <= 5 MB per doubling | b > 1.3 or RSS grows with n |
| E8 | optional, needs disk: 1 Gbp bounded-scratch, then 3 Gbp | genome-like | 3 Gbp extrapolated <= 3 h (the pSAscan class V/E); measured if disk allows | extrapolated > 6 h |
| E9 | PH6: Psi by sigma-way distribution vs ISA by DIST; suffix-link smoke test through the existing P1 prototype of ADR-001 | `dna_10m`, `prot_1m` | Psi <= 6 B/char of I/O; link kernel results identical with Psi and with ISA | - |

**RAM-cap enforcement and verification**: (1) from inside: every big buffer comes from a `NativeArena` that carries an `ArenaBudget` (sum of live blocks <= 160 MB); `DOTNET_GCHeapHardLimit` set to 96 MB caps the managed heap; windows are bounded; (2) from outside: the RssAnon sampler (anonymous resident pages only; `VmHWM`/`ru_maxrss` also count mapped file pages and are therefore reported but not gated) with the gate `max RssAnon <= 256 MB`; (3) `ulimit -v` is **not** used (the .NET runtime reserves tens of GB of address space); (4) a cgroup memory limit is used only as the page-cache stress test of E3, not as the RAM gate, because mapped file pages count in cgroup v1 and are reclaimable.

**Oracles per stage**: PH1/PH3 against `SaisReference` (n <= 2000 exhaustive and random, adversarial families) and the libsais hashes up to 10 M; PH4 against `LcpReference`/Kasai; PH5 against a reference child table built in RAM and, for tiny texts, against the Ukkonen tree (`CalculateLogicalHash`, ADR-001 R-HASH); a differential fuzz harness with scaled-down K (K = 3..50, F = 2..4, X = 7..64) so that every multi-bucket, multi-level and spill path is exercised on texts of 10^3-10^5 symbols.

**Acceptance for adopting the recommendation**: E1/E2, E3, E4-E7 all pass. **Abandon (switch to F1)**: an adversarial family stays > 25x after tuning; RssAnon cannot be held <= 256 MB without O(n) state; total at 100 Mbp > 4.5 min; scaling exponent > 1.3; the build needs a page cache larger than 2x the text + rank file to stay within 3x.

## 5a. Prototype results (measured)

Added 2026-11-09 after the owner approved the prototypes E0-E3 (the recommendation text of this ADR is unchanged: the owner decides B-DCS versus F1). Code: `apps/SuffixTree.EmProto`; raw rows and the full tables: [analysis/em/RESULTS.md](../analysis/em/RESULTS.md). Everything here is **measured** on the prototype machine (A5) unless marked otherwise.

### E0 harness - gates met

| gate | measured | result |
|---|---|---|
| sampler overhead < 2 % | 0.64 % (15 interleaved pairs, `SaisKernel` on 10 M random DNA; A/A noise band 0.9917-1.0083, i.e. inside the noise) | PASS |
| `/proc/self/io` agreement within 1 % | writes 0 %, reads 0.047 % (302 MB written, 336 MB read through `CountedFile`) | PASS |
| `ArenaBudget` | enforced and released, thread-safe (tests); a wrapper class around `NativeArena`, which is not modified | PASS |

### E1 in-RAM simulation without the difference-cover sample - gate met, DCS needed

Measured with `e1` of `apps/SuffixTree.EmProto` (1 thread, K = n/1000, fan-out 64, direct comparison with an AVX2 LCE over the b-bit packed text; tables and raw rows in [RESULTS.md](../analysis/em/RESULTS.md)). Every suffix array equals the libsais hash or the `SaisKernel` result.

| family | n | `sa_s` median [min-max], fresh processes | ratio to `dna_10m` | key-tie share | LCE work |
|---|---|---|---|---|---|
| `dna_10m` | 10 M | 1.006 s [0.973-1.041], 5 runs (gate <= 3 s: **met**; 0.556 s with 4 threads) | 1 | 0 % | 1.4 x 10^6 symbols |
| `rep_10m` | 10 M | 2.58 s [2.43-2.65], 3 | 2.6 | 77 % | 3.4 x 10^9 |
| `prot_1m` | 1 M | 0.258 s, 3 | 0.26 | 0 % | 5.7 x 10^5 |
| Thue-Morse | 1 M | 4.14 s [4.03-4.27], 3 | 4.1 | 99.99 % | 3.0 x 10^11 |
| Fibonacci | 1 M | 9.41 s [9.28-9.61], 3 | 9.4 | 99.99 % | 7.3 x 10^11 |
| genome-like (own generator) | 10 M | 2.22 s [2.11-2.33], 3 | 2.2 | 8.6 % | 2.6 x 10^10 |
| a^n | 1 M | 118 s [117.9-119.1], 3 | **117** | 100 % | 1.0 x 10^13 |
| a^n | 10 M | **> 600 s, aborted** (failure by the rule) | > 600 | - | - |
| 171-bp tandem x 10^4 | 1.71 M | 505 s [502-511], 3 | **502** | 100 % | 2.1 x 10^13 |
| N run 1 M in 10 M random DNA | 10 M | 412 s [410-415], 3 | **410** | 10 % | 1.2 x 10^13 |

With the merge sort of the final code (fewer comparator calls) the numbers are lower by 11-39 % but the picture is the same: a^n 1 M 88x, tandem 301x, N run 240x, a^n 10 M still above 600 s. **Rule of E1: not every family is within 10x, so the difference-cover sample (PH1) is needed and E2 was run.** The cost of a failing family is not the number of comparator calls (1.4-2.6 x 10^7 for all of them) but their length: 10^6 matched symbols per text symbol.

### E2 with the difference-cover sample (in-RAM ranks) - gates met

The difference covers are constructed, not tabulated: Singer perfect difference sets for `X = q^2 + q + 1` (q prime; nominal 512 / 1024 / 2048 become 553 / 993 / 1893 with |D| = 24 / 32 / 44) and a pruned square-root construction for any other X (exact 512 / 1024 / 2048 give |D| = 33 / 47 / 67, i.e. 6.5 / 4.6 / 3.3 % of the positions). **Proof:** the default test covers every X from 2 to 2300 (all residue differences occur) and checks the Singer sets for perfection. The comparator is Karkkainen's (compare l < X symbols, then two sample ranks).

Round 1 (Singer 993, K = n/1000, 1 thread, ranks built in RAM and charged to the time): `dna_10m` 1.052 s; every adversarial family within 10x: a^n 1 M 0.8x, tandem 1.1x, N run 1.8x, Thue-Morse 1 M 1.0x, Fibonacci 1 M 0.9x, rep_10m 2.8x, genome-like 1.5x, and at 10 M a^n 6.0x, Thue-Morse 6.7x, **Fibonacci 7.9x (the largest)**. The cost of the DCS on `dna_10m`: 7.4 % (paired A/B, 9 pairs; gate 25 %). All outputs equal the libsais hashes or `SaisKernel`. **No tuning round was needed.** Sweep over X (553, 993, 1893 Singer; 512, 1024, 2048 exact), K = n/300 and n/3000 and 4 threads: every configuration passes (worst ratio 8.5x); smaller X is faster on repetitive text (a^n 10 M 4.1 s at 512, 6.3 s at 2048) and costs a denser rank file (6.5 % against 2.3 % of the positions). Details: [RESULTS.md](../analysis/em/RESULTS.md).

## 6. Honest limitations and open questions

What I could **not** verify:
* The primary evidence for pSAscan (CPM 2015) and fSAIS (ALENEX 2017) could not be opened (author-site PDFs returned 404); their per-dataset tables are unread, numbers come from the README, the ESA 2016 Table 4 and the thesis summary. B2ST, DiGeST, TDD, Trellis, WaveFront primary papers, Flick-Aluru, Fischer-Kurpicz 2019, the JEA 2008 pDC3 journal version (I read the ALENEX 2005 version), the ALENEX 2013 eSAIS paper (I read the thesis chapter), BCR experiments, nSAIS and the 2020 EM paper were not read; PFP only as a search record.
* Nothing was run (the machine was busy); every number about B-DCS is an estimate (E). The literature numbers come from 2005-2017 hardware, 8-bit symbols, mostly one input family each.
* **B-DCS in external memory is an extrapolation.** Karkkainen's paper measured 1 GB on one core with the DCS in RAM; the EM bucketing is only a sketch in Haag et al. section 5; the recursion that builds the DCS and the file-resident rank lookups are untested. The cost of the DCS comparator on satellite-rich genomes (T2T has about 60 Mbp of alpha satellites) is unknown; E1/E2 exist to measure it.
* The difference-cover table for X = 512..4096 (size, offsets `l(a,b)`) is not verified; known constructions give |D| about sqrt(1.5 X).
* The effect of the page cache is the main hardware dependence: the text and the rank file are read randomly (about 0.45-0.55 B/char for 2-bit DNA, 1.2-1.3 B/char for bytes).
* Whether sliding mapped windows keep dirty pages low on this kernel, and whether the cgroup memory controller can be used by the prototype user, are untested.
* A WaveFront patent exists (US8914415B2, U); not relevant to the recommendation but relevant if the ERa/WaveFront family is ever considered.
* The ESA 2016 numbers are for reads (1000 Genomes), not an assembled genome, and for 512 GiB with 120 GiB RAM (n/M = 4.3); my scaling to n/M = 12 is linear in n/M (E).

Open questions for the owner:
* **Q1** Are scratch files allowed (about 5-12 B/char peak, bounded by `S_max`, on the output volume) and is the DCS rank file (0.04-0.14 entries per char, i.e. 0.2-0.6 B/char) acceptable as "index-derived scratch"? A strictly scratch-free build is only possible with a non-robust algorithm.
* **Q2** Is the 256 MB cap on anonymous memory including the .NET runtime baseline and thread stacks (assumed yes)?
* **Q3** DNA with N and IUPAC symbols: sigma is 5-16, not 4. Pack as b = 3-4 bits (key = 16-21 symbols), or 2 bits plus an exception list for N runs? Are both strands (n = 6.2x10^9, `POS_WIDTH` 8 or 5) in scope?
* **Q4** Is a design that reads the text and the rank file randomly through the OS page cache acceptable, and what is the minimum RAM of the target machine?
* **Q5** Must the output be byte-identical for any thread count and M (the design provides it)?
* **Q6** For n <= N_base (16 M symbols) is "in-RAM SA-IS, then stream the sections" acceptable (working memory is bounded by N_base)?
* **Q7** Priority: SA+LCP+child only (Core), or the pack (Psi + RMQ) in the same build? ISA for `RankOf`/`ExtendLeft`: needed in Core?
* **Q8** Disk for the 1-3 Gbp gates: 40 GB free allows 1 Gbp unbounded, 60 GB allows 3 Gbp; otherwise the largest local run is 500 Mbp and 3 Gbp stays an extrapolation.
* **Q9** Licence policy: GPL sources (eSAIS, Bowtie 2) are read as documentation only, MIT sources (pSAscan, fSAIS) are used as read-only references and nothing is copied (PROCESS.md rule 1 extended to all third-party code).

## 6a. Owner answers to Q1-Q4 (2026-11-04) and the constraints that follow

* **Q1 scratch files: allowed, "without fanaticism"** - the build must not be slowed down badly and must not fill the disk with temporary files. Constraints: scratch lives on the output volume, is bounded by `S_max`, is released stage by stage (a scratch file is deleted as soon as its consumer finished), and the prototype gates E1-E3 report **peak scratch (B/char) and total I/O (B/char)** next to wall time. Proposed acceptance (owner to confirm): peak scratch <= 8 B/char and total I/O <= 100 B/char; if a design needs more, report it instead of silently exceeding.
* **Q2: the 256 MB cap includes the .NET runtime baseline and thread stacks.** Measured as peak anonymous memory (`RssAnon`, plus a cgroup limit when the prototype user may create one). Mapped-file pages are not counted.
* **Q3: text codec = b bits per symbol, 4 bits for DNA with N/IUPAC, with a 2-bit fast path - and it must be extensible without rewriting.** Requirements: the TEXT section header records `codec` and `bits_per_symbol` (the ADR-003 codec field is used, not replaced); every stage takes `b` as a parameter (the 64-bit bucket key holds `64/b` symbols; splitter, bucket and merge logic never assume sigma = 4); no file layout depends on `b`; unit tests are parameterised over b in {2, 3, 4, 8} from the first prototype (small inputs); widening to 8-bit bytes or 16-bit UTF-16 is a **new codec value, not a format change**. The 2-bit path is chosen when the alphabet scan finds sigma <= 4, b = 4 when sigma <= 16 (A, C, G, T, N and IUPAC). An exception-list codec for long N runs is a possible later codec id, not built now. Both strands are out of scope here (S3, query side).
* **Q4: dependence on the OS page cache for the random reads is acceptable.** No minimum RAM for the target machine is fixed; gate E3 states how throughput degrades with the available page cache.
* Accepted by default (the owner did not object): **Q5** byte-identical output for any thread count and `M`; **Q6** in-RAM SA-IS as the base case for n <= N_base (16 M symbols, memory bounded by N_base, independent of n); **Q7** Core first (SA + LCP + child table), the navigation pack (Psi + RMQ) afterwards; **Q8** local runs up to what the free disk allows (about 22 GB free now, so roughly 100-200 Mbp with scratch), 3 Gbp as a hand-off package for a machine with about 60 GB; **Q9** GPL sources (eSAIS, Bowtie 2) are read as documentation only, MIT sources (pSAscan, fSAIS) as read-only references, nothing is copied.

## 7. Sources actually opened (2026-11-04)

Local: PROGRESS-S1.md (top), ADR-001, ADR-002, ADR-003, BENCHMARK.md (grep), PROCESS.md rules, `src/SuffixTree/Algorithms/SuffixTree.Esa` (file list, `Buf`, `NativeArena`, `LcpKernel` signatures), `analysis/bench/gen/gen_inputs.py`.
* Karkkainen, Fast BWT in small space by blockwise suffix sorting, TCS 387 (2007): https://www.cs.helsinki.fi/u/tpkarkka/publications/tcs07.pdf (abstract, sections 1-9, Table 1)
* Dementiev, Karkkainen, Mehnert, Sanders, ALENEX 2005 version: https://ae.iti.kit.edu/documents/people/sanders/papers/DKMS05.pdf (abstract, Theorem for DC3, Table 1, experiments)
* Bingmann, Scalable String and Suffix Sorting (dissertation, includes eSAIS, chapter 6, Table 5.1): https://arxiv.org/pdf/1808.00963
* Karkkainen, Kempa, Faster External Memory LCP Array Construction, ESA 2016: https://drops.dagstuhl.de/storage/00lipics/lipics-vol057-esa2016/LIPIcs.ESA.2016.61/LIPIcs.ESA.2016.61.pdf (Tables 1-4)
* Ferragina, Gagie, Manzini, Lightweight data indexing and compression in external memory: https://arxiv.org/pdf/0909.4341 (Thm 3.1, 4.2, 4.3, Fig. 2); abstract view https://ar5iv.arxiv.org/html/0909.4341
* Mansour, Allam, Skiadopoulos, Kalnis, ERa: https://arxiv.org/pdf/1109.6884 (abstract, Table 2, section 3); Comin, Farreras (PCF): https://pmc.ncbi.nlm.nih.gov/articles/PMC3962650
* Khan et al., CaPS-SA: https://pmc.ncbi.nlm.nih.gov/articles/PMC11056320/ (Tables 1, 2, 4, memory, limitations); README https://raw.githubusercontent.com/jamshed/CaPS-SA/main/README.md
* Haag, Kurpicz, Sanders, Schimek: https://arxiv.org/pdf/2412.10160 (section 5), https://arxiv.org/abs/2412.10160, https://drops.dagstuhl.de/storage/00lipics/lipics-vol351-esa2025/html/LIPIcs.ESA.2025.47/LIPIcs.ESA.2025.47.html
* pSAscan page and README: https://www.cs.helsinki.fi/group/pads/pSAscan.html, https://raw.githubusercontent.com/dominikkempa/psascan/master/README.md; fSAIS README: https://raw.githubusercontent.com/dominikkempa/fsais/master/README.md; eSAIS page: https://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-LCP-Arrays-in-External-Memory/
* Bowtie 2 blockwise suffix array source (read as documentation): https://raw.githubusercontent.com/BenLangmead/bowtie2/master/blockwise_sa.h
* extLCP abstract: https://arxiv.org/abs/1605.04098
* Search records only (U): Flick-Aluru SC15 (7.3 s, 1024 cores), Fischer-Kurpicz ALENEX 2019 (KIT record), libsais (docs.rs, bazel registry), Boucher et al. prefix-free parsing, DiGeST thesis (UVic), WaveFront patent, nSAIS listing, BCR TCS 2013 record.
