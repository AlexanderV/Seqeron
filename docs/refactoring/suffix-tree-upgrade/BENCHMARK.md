# Suffix-tree benchmark: Seqeron.SuffixTree vs Go stdlib and open-source competitors

Measured on one machine, 3 repeats (median) unless stated.  Seqeron = commit `242318429a5e` (working tree clean; library sources were copied to
`bench/seqeron/snapshot` and built from there, nothing under `/home/user/Seqeron` was touched).  All numbers below are generated from `bench.tsv`
by `lib/make_report.py`; raw per-repeat rows are in `runs/raw*.tsv`.  Environment: `ENV.txt`.  Provenance and licences of every third-party source: `SOURCES.txt`.

## 0. Read this first (what is and is not comparable)

* **Different data structures.** Seqeron and `sdsl cst_sct3` / Python `suffix-trees` are *suffix trees*; Go `index/suffixarray`, libsais, libdivsufsort and pydivsufsort
  are plain *suffix arrays* (4 B/char, binary-search queries); rust-bio and `sdsl csa_*` are *FM-indexes* (compressed suffix arrays);
  BWA and Bowtie2 are *FM-indexes for read alignment* (index build only was measured).  A suffix array answers count/locate/contains but has no LRS, MEM/MUM,
  matching statistics, suffix links or tree navigation; those operations are exactly what a suffix tree is bought for.  Numbers therefore compare **cost for the shared
  functionality (build, size, count, locate)**; they do not say the other tools replace the tree's extra functions.
* libsais / libdivsufsort / pydivsufsort are *construction* libraries: their count/locate rows are *our* binary-search code (`competitors/sa_c/sa_bench.c`) over their SA.
* Every subject returned **identical count and locate checksums** on every input it ran (`crosscheck.md`; 0 disagreements) so the query sets really are the same workload.
* The VM reports a 480 MiB L3 and 4 vCPU (KVM); random-access structures up to a few hundred MB may be served largely from that cache, so the 1-10 Mbp rows are cache-friendly relative to a bare-metal desktop.

## 1. Headline: 10 Mbp random DNA (`dna_10m`)

| metric | Seqeron in-memory | Seqeron persistent (mmap) | Go `index/suffixarray` | best other tool (that metric) | Seqeron in-mem vs best | Seqeron persistent vs best |
|---|---|---|---|---|---|---|
| build wall time (s) | 14.6 | 14.9 | 0.73 | 0.328 (libsais) | 44.6x slower | 45.4x slower |
| index bytes per input char | 169 (GC heap; +2 B/char text string) | 56.2 (file; text included) | 4 (+1 B/char text) | 0.27 (sdsl-csa_wt_rrr) | 624.4x larger | 207.8x larger |
| peak RSS during build (MB) | 1,761 | 1,179 (of which anonymous 46; rest = mmap'd file pages) | 53 | 52 (min over 9 tools) | - | - |
| count, 20-mers present (queries/s, 1 thread) | 589,810 | 388,371 | 573,327 | 747,415 (libsais) | 1.3x slower | 1.9x slower |
| count, 4 threads, one shared index (queries/s) | 2,561,620 | 1,716,118 | 2,408,974 | 2,858,062 (libsais) | 1.1x slower | 1.7x slower |
| locate, 1 000 random 8-mers, ~152 hits each (patterns/s) | 49,653 | 39,943 | 463,542 | 1,140,278 (rust-bio) | 23.0x slower | 28.5x slower |
| locate throughput (positions returned/s) | 7,565,278 | 6,085,774 | 70,626,253 | 173,735,030 | - | - |
| count latency p50 / p99 (ns, 1 thread) | 1,969 / 2,832 | 2,580 / 3,961 | 1,654 / 2,649 | libsais 1,349 / 2,175 | - | - |

"Best other tool" ranges over Go, libsais, libdivsufsort, rust-bio, sdsl (4 configs) and pydivsufsort -- i.e. over suffix arrays *and* FM-indexes, so its bytes/char is an FM-index and its
build time a suffix-array construction.  Same-kind comparators (suffix trees) are in section 3.

## 2. All tools, all inputs (medians of 3 repeats; see section 6 for exceptions)

Columns: build wall time; index bytes per input char (Seqeron in-memory = GC heap delta, persistent = file size, Go/libsais/divsufsort = SA only, FM-indexes/CST = `size_in_bytes`, BWA/Bowtie2 = all index files); peak RSS during build; single-thread 20-mer count (present / absent), locate of 1 000 8-mers, 4-thread shared-index count.  `-` = not applicable / not run.  For `prot_1m` the 8-mer locate set has almost no hits (20^8 >> 1e6): use the `locate4` column.


### dna_1m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-inmem | 0.647 | 169 | 206 | 1,147,184 | 1,396,889 | 493,411 | - | 4,720,062 | 3 |
| seqeron-st-persistent | 0.921 | 56.2 | 149 | 943,168 | 1,173,994 | 378,616 | - | 4,040,552 | 3 |
| go-suffixarray | 0.0548 | 4.01 | 8 | 1,451,218 | 1,497,450 | 1,381,132 | - | 5,681,200 | 3 |
| libsais | 0.0261 | 4 | 7 | 1,766,037 | 1,766,237 | 1,891,311 | - | 6,810,965 | 3 |
| libdivsufsort | 0.0479 | 4 | 7 | 1,777,806 | 1,691,060 | 1,907,876 | - | 6,810,468 | 3 |
| rust-bio | 0.101 | 22.3 | 35 | 947,831 | 2,028,829 | 5,334,494 | - | 3,669,292 | 3 |
| sdsl-csa_wt_huff | 0.131 | 0.56 | 31 | 2,095,314 | 3,759,384 | 70,082 | - | 7,708,584 | 3 |
| sdsl-csa_wt_rrr | 0.135 | 0.271 | 31 | 39,582 | 74,607 | 66 | - | 152,848 | 3 |
| sdsl-csa_sada | 0.173 | 0.667 | 37 | 23,785 | 44,945 | 11,934 | - | invalid (race) | 3 |
| sdsl-cst_sct3 | 0.219 | 1.65 | 43 | 34,839 | 48,268 | 33,267 | - | 135,552 | 3 |
| pydivsufsort | 0.0511 | 4 | 34 | 124,961 | 128,103 | 106,096 | - | - | 3 |
| bwa-index | 0.317 | 1.75 | 13 | - | - | - | - | - | 3 |
| bowtie2-build | 0.567 | 9.81 | 98 | - | - | - | - | - | 3 |

### dna_10m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-inmem | 14.6 | 169 | 1,761 | 589,810 | 681,703 | 49,653 | - | 2,561,620 | 3 |
| seqeron-st-persistent | 14.9 | 56.2 | 1,179 | 388,371 | 494,770 | 39,943 | - | 1,716,118 | 3 |
| go-suffixarray | 0.73 | 4 | 53 | 573,327 | 599,754 | 463,542 | - | 2,408,974 | 3 |
| libsais | 0.328 | 4 | 52 | 747,415 | 736,424 | 633,517 | - | 2,858,062 | 3 |
| libdivsufsort | 0.662 | 4 | 52 | 709,676 | 724,083 | 616,734 | - | 2,797,430 | 3 |
| rust-bio | 1.87 | 22.3 | 290 | 337,999 | 728,647 | 1,140,278 | - | 1,528,175 | 3 |
| sdsl-csa_wt_huff | 1.64 | 0.561 | 154 | 659,606 | 1,161,632 | 2,147 | - | 2,593,714 | 3 |
| sdsl-csa_wt_rrr | 1.72 | 0.27 | 154 | 35,825 | 56,173 | 7 | - | 138,395 | 1 |
| sdsl-csa_sada | 2.06 | 0.698 | 201 | 22,212 | 36,301 | 1,018 | - | invalid (race) | 3 |
| sdsl-cst_sct3 | 2.69 | 1.65 | 192 | 18,156 | 23,742 | 1,912 | - | 72,546 | 1 |
| pydivsufsort | 0.537 | 4 | 79 | 112,284 | 115,321 | 96,167 | - | - | 3 |
| bwa-index | 5.01 | 1.75 | 57 | - | - | - | - | - | 3 |
| bowtie2-build | 5.31 | 2.26 | 113 | - | - | - | - | - | 3 |

### rep_10m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-inmem | 9.15 | 193 | 2,001 | 406,748 | 1,060,786 | 31,417 | - | 1,822,682 | 3 |
| seqeron-st-persistent | 6.94 | 64.1 | 1,338 | 307,818 | 854,737 | 21,766 | - | 1,323,536 | 3 |
| go-suffixarray | 0.609 | 4 | 53 | 509,021 | 664,608 | 456,911 | - | 1,892,410 | 3 |
| libsais | 0.279 | 4 | 52 | 670,690 | 755,105 | 644,427 | - | 2,555,993 | 3 |
| libdivsufsort | 0.55 | 4 | 52 | 671,145 | 740,119 | 648,866 | - | 2,589,932 | 3 |
| rust-bio | 1.73 | 22.3 | 273 | 420,622 | 1,127,131 | 1,328,564 | - | 1,910,824 | 3 |
| sdsl-csa_wt_huff | 1.43 | 0.56 | 154 | 647,631 | 1,550,976 | 5,382 | - | 2,544,974 | 3 |
| sdsl-csa_wt_rrr | 1.46 | 0.164 | 154 | 64,666 | 123,065 | 17 | - | 235,830 | 1 |
| sdsl-csa_sada | 1.88 | 0.375 | 201 | 46,556 | 86,014 | 3,110 | - | invalid (race) | 3 |
| sdsl-cst_sct3 | 2.77 | 2.22 | 192 | 16,187 | 26,637 | 4,401 | - | 65,886 | 1 |
| pydivsufsort | 0.437 | 4 | 79 | 105,737 | 112,690 | 95,824 | - | - | 3 |
| bwa-index | 6.13 | 1.75 | 56 | - | - | - | - | - | 3 |
| bowtie2-build | 10 | 2.26 | 113 | - | - | - | - | - | 3 |

### prot_1m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-inmem | 1.12 | 190 | 232 | 1,072,479 | 1,485,293 | 2,050,086 | 730,733 | 4,527,560 | 3 |
| seqeron-st-persistent | 1.34 | 45 | 128 | 1,575,900 | 2,229,731 | 5,573,478 | 694,964 | 6,680,381 | 3 |
| go-suffixarray | 0.0436 | 4.01 | 8 | 1,489,924 | 1,532,184 | 1,676,161 | 1,412,532 | 5,467,715 | 3 |
| libsais | 0.03 | 4 | 7 | 1,712,898 | 1,749,408 | 1,929,280 | 1,852,218 | 6,664,857 | 3 |
| libdivsufsort | 0.0423 | 4 | 7 | 1,747,483 | 1,765,609 | 2,026,658 | 1,960,601 | 6,780,176 | 3 |
| rust-bio | 0.139 | 65 | 99 | 316,522 | 2,209,171 | 9,621,890 | 10,983,654 | 1,676,953 | 3 |
| sdsl-csa_wt_huff | 0.139 | 0.97 | 31 | 935,309 | 3,838,941 | 3,902,355 | 77,900 | 3,644,648 | 3 |
| sdsl-csa_wt_rrr | 0.15 | 0.571 | 31 | 18,115 | 70,634 | 71,573 | 76 | 71,231 | 3 |
| sdsl-csa_sada | 0.168 | 1.09 | 37 | 21,745 | 95,251 | 94,161 | 15,038 | invalid (race) | 3 |
| sdsl-cst_sct3 | 0.215 | 1.93 | 43 | 25,583 | 58,385 | 58,735 | 37,573 | 96,744 | 3 |
| pydivsufsort | 0.0383 | 4 | 34 | 119,587 | 123,280 | 127,070 | 105,997 | - | 3 |

### dna_20m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-persistent | 39.1 | 56.2 | 2,322 | 338,069 | 410,291 | 19,602 | - | 1,469,098 | 3 |
| go-suffixarray | 1.55 | 4 | 104 | 490,346 | 505,359 | 324,338 | - | 2,011,314 | 3 |
| libsais | 0.739 | 4 | 102 | 606,052 | 604,171 | 532,525 | - | 2,410,347 | 3 |
| libdivsufsort | 1.44 | 4 | 102 | 608,049 | 593,739 | 505,240 | - | 2,501,848 | 3 |
| bwa-index | 11.8 | 1.75 | 62 | - | - | - | - | - | 3 |
| bowtie2-build | 10.1 | 1.84 | 134 | - | - | - | - | - | 3 |

### dna_50m

| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |
|---|---|---|---|---|---|---|---|---|---|
| seqeron-st-persistent | 389 | 56.2 | 5,904 | 248,884 | 285,753 | 5,481 | - | 1,235,441 | 3 |
| go-suffixarray | 4.24 | 4 | 257 | 374,504 | 381,718 | 208,977 | - | 1,502,073 | 3 |
| libsais | 2.09 | 4 | 252 | 428,094 | 429,592 | 374,705 | - | 1,943,297 | 3 |
| libdivsufsort | 4.11 | 4 | 252 | 420,786 | 432,399 | 366,693 | - | 1,880,427 | 3 |

## 3. Suffix trees only: scaling and same-kind comparison

Seqeron persistent (mmap) tree vs the other suffix-array / compressed-tree builds at growing n (random DNA):

| n (Mbp) | Seqeron persistent build s | Seqeron chars/s | file GB | peak RSS GB (file-backed + anon) | anon MB | count q/s | locate8 pat/s | MEM 100 kbp query s | Go build s | libsais build s | libsais_omp4 build s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.921 | 1,086,166 | 0.06 | 0.15 | 10 | 943,168 | 378,616 | 0.124 | 0.0548 | 0.0261 | 0.0185 |
| 10 | 14.9 | 670,877 | 0.56 | 1.18 | 46 | 388,371 | 39,943 | 0.163 | 0.73 | 0.328 | 0.179 |
| 20 | 39.1 | 512,086 | 1.12 | 2.32 | 86 | 338,069 | 19,602 | 0.168 | 1.55 | 0.739 | 0.362 |
| 50 | 389 | 128,558 | 2.81 | 5.90 | 358 | 248,884 | 5,481 | 0.2 | 4.24 | 2.09 | 0.997 |

* **Largest Seqeron size that finished: 50 Mbp** (persistent tree, 389 s, 3 of 3 repeats, 2.81 GB file, peak RSS 5.9 GB of which only 358 MB anonymous; the rest is the mmap'd file).
  Build rate falls from 670,877 chars/s at 10 Mbp to 128,558 chars/s at 50 Mbp (at 50 Mbp the Ukkonen phase is 92 % of the time, i.e. random access to the node file, not the finalize passes).
* **100 Mbp did not finish**: one attempt with a 600 s cap was killed in the Ukkonen phase (`runs/raw.tsv` row `dna_100m ... BLOCKED`).  At the kill the main storage file held 1 531 MB of
  24-byte nodes = ~64 M of the ~162 M nodes a 100 Mbp tree needs, i.e. ~40 % of the text after 600 s (estimate), RSS 4.8 GB (4.3 GB file-backed, 0.45 GB anonymous); the sample log is `runs/mon_100m.keep.txt`.
  Extrapolating the 10 -> 20 -> 50 Mbp trend, a 100 Mbp build would need well over 30 min.  The text length is also an `int` in the API (< 2^31 chars) -- 3 Gbp is out of reach by construction, not only by time.
* Construction is single-threaded: CPU time = wall time (380 s CPU for 389 s wall at 50 Mbp).  libsais's OpenMP build gives 2.1x at 4 threads on 50 Mbp.
* Build stage breakdown at 10 Mbp (persistent): Ukkonen 10.7 s, finalize (child arrays + leaf counts + write-back + trim) 4.22 s.

Same-kind comparators (index = suffix *tree*), 1 Mbp -> 10 Mbp:

| tool | n (Mbp) | build s | B/char | peak RSS MB | count q/s | locate8 pat/s |
|---|---|---|---|---|---|---|
| seqeron-st-inmem | 1 | 0.647 | 169 | 206 | 1,147,184 | 493,411 |
| seqeron-st-inmem | 10 | 14.6 | 169 | 1,761 | 589,810 | 49,653 |
| seqeron-st-persistent | 1 | 0.921 | 56.2 | 149 | 943,168 | 378,616 |
| seqeron-st-persistent | 10 | 14.9 | 56.2 | 1,179 | 388,371 | 39,943 |
| sdsl-cst_sct3 | 1 | 0.219 | 1.65 | 43 | 34,839 | 33,267 |
| sdsl-cst_sct3 | 10 | 2.69 | 1.65 | 192 | 18,156 | 1,912 |
| suffix-trees-py (pure Python) | 0.01 | 0.0211 | 2,832 (RSS/n) | 28 | 270,596 (contains) | 307,691 |
| suffix-trees-py (pure Python) | 0.1 | 0.6 | 914 (RSS/n) | 91 | 72,277 (contains) | 112,595 |
| suffix-trees-py (pure Python) | 1 | 13 | 766 (RSS/n) | 766 | 7,110 (contains) | 30,564 |
| suffix-trees-py (pure Python) | 3 | 49.7 | 757 (RSS/n) | 2,271 | 2,026 (contains) | 11,187 |
| suffix-trees-py (pure Python) | 10 | 196 | 753 (RSS/n) | 7,533 | 592 (contains) | 3,520 |

Python `suffix-trees` 0.4.0 **gives up between 10 and 20 Mbp**: 10 Mbp builds in 196 s using 7.5 GB; the 20 Mbp attempt hit `MemoryError` at the 12 GB address-space cap after ~800 s (one repeat observed, recorded as BLOCKED).
Seqeron's in-memory tree is 13x faster to build than that Python tree at 10 Mbp and uses 4.3x less memory; `sdsl cst_sct3` builds the same 10 Mbp in 2.69 s with 1.65 B/char (compressed tree, navigation via `sdsl` API, count by descent -- 18,156 q/s).

## 4. Seqeron-specific measurements

| metric | in-memory 10 Mbp | persistent 10 Mbp | persistent 50 Mbp |
|---|---|---|---|
| LongestRepeatedSubstring first call (s) | 0.00019 | 0.0015 | 0.0017 |
| MEM, 100 kbp query w/ 1 % substitutions, minLength 20: first call / warm (s) | 0.0623 / 0.0459 | 0.163 / 0.143 | 0.2 / 0.155 |
| MEM matches returned | 817 | 817 | 802 |
| MUM (`Both`) time (s) / matches | 0.0391 / 816 | 0.105 / 816 | 0.103 / 802 |
| FindExactMatchAnchors time (s) | 0.032 | 0.0877 | 0.0852 |
| MEM on the repetitive 10 Mbp text (`rep_10m`): time (s) / matches | 0.399 / 160,979 | 1.18 / 160,979 | - |
| cold `Load` time, page cache dropped (ms) | - | 4.84 | 7.11 |
| warm `Load` time (ms) | - | 3.93 | 3.92 |
| first 100 000 count queries right after cold / warm load (s) | - | 0.488 / 0.346 | 1.81 / 0.535 |
| 4-thread / 1-thread count throughput (one shared tree) | 4.60x | 4.43x | 4.42x |
| 4-thread / 1-thread locate8 throughput | 3.87x | 3.88x | 4.21x |

Go reference for the same two rows: cold `Index.Read` of the serialized index 111 ms at 10 Mbp (47.9 MB serialized, 4.79 B/char; no mmap, whole index copied into RAM);
4T/1T count scaling 4.01x.

Notes on the numbers above
* 4T/1T ratios slightly above 4.0 on a 4-vCPU VM are measurement noise in the 1-thread baseline (turbo/steal), read them as "scales ~linearly to 4 threads".  Concurrent reads of one shared tree returned correct results in every run (`checksum_ok = 1`).
* Persistent query throughput is 1.24-1.52x lower than in-memory on 10 Mbp (count / absent / locate8 / 4T count); the persistent tree is queried through the page cache here (file just built / loaded, 15 GB RAM) -- a true out-of-core run (index larger than RAM) was **not** done.
* In-memory index = 169 B/char of managed heap (16,229,747 nodes for 10,000,000 leaves); process peak RSS 1.76 GB at 10 Mbp.
* Locate is the weak spot: 7.6 M positions/s (in-memory) and 6.1 M/s (persistent) against 71 M/s for Go and 97 M/s for a libsais SA slice (a suffix array returns a contiguous run of the SA; the tree has to walk the subtree).
  At 50 Mbp the 8-mer locate drops to 5,481 patterns/s (762,118 hits per 1 000 patterns).
* Build memory is the in-memory tree's problem, not the persistent tree's: persistent build anonymous memory stays at 358 MB for 50 Mbp.

## 5. What each tool provides (functionality actually exercised here)

| tool | index kind | exercised operations | not provided / not measured |
|---|---|---|---|
| Seqeron.SuffixTree (in-memory) | suffix tree (Ukkonen), UTF-16 chars | Build, Contains, CountOccurrences, FindAllOccurrences, LongestRepeatedSubstring, FindMaximalExactMatches, FindMaximalUniqueMatches, FindExactMatchAnchors, 4-thread shared reads | - |
| Seqeron.SuffixTree.Persistent | suffix tree in a memory-mapped file (v6 format) | same + Create(file) / Load(file) | out-of-core (index > RAM) not tested |
| Go `index/suffixarray` | suffix array (SA-IS), bytes | New, Lookup(n=1 / n=-1), Write, Read | no count API (count = len(Lookup)); no LRS / MEM / MUM / anchors; whole index must be in RAM |
| libsais 2.9.1, libdivsufsort 2.0.1 | suffix-array construction (int32) | construction (libsais also OpenMP 2/4 threads) | no query API: count/locate = our binary search |
| pydivsufsort 0.0.20 | libdivsufsort bindings | divsufsort, `sa_search` from a Python loop | queries are Python-call-bound |
| rust-bio 2.3.0 | SA + BWT + FM-index (Occ k=3), locate with the full `usize` SA | suffix_array, backward_search, Interval::occ | no LRS/MEM; sampled SA not used |
| sdsl-lite 2.1.1 `csa_wt<wt_huff>` / `csa_wt<wt_huff<rrr>,512,1024>` / `csa_sada` | compressed suffix arrays (FM-index / Psi) | construct_im, count, locate | no LRS/MEM API |
| sdsl-lite `cst_sct3<>` | compressed suffix tree | construct_im, descent with child()/edge(), size(v), locate via csa[lb..rb], LRS (user DFS, 1 Mbp only) | no MEM/MUM API |
| Python `suffix-trees` 0.4.0 | pure-Python Ukkonen suffix tree | STree build, find, find_all | no LRS/MEM measured |
| BWA 0.7.19 `index -a bwtsw`, Bowtie 2 2.5.5 `bowtie2-build` | FM-index for read alignment | index build time / peak RSS / on-disk bytes | no count/locate harness (use the aligners' own `mem`/`bowtie2`) |

## 6. Methodology and caveats

* **Procedure.** One process = one measurement; each (tool, input, phase) repeated 3 times (median reported, min/max in the `notes` column of `bench.tsv`).  Exceptions: `sdsl csa_wt_rrr` and `cst_sct3` on 10 Mbp inputs ran once
  (`csa_wt_rrr` took 6-16 min per run at 10 Mbp because its locate costs ~256 LF-steps per hit on average at sample rate 512; `cst_sct3` took ~75 s); `suffix-trees` 20 Mbp and Seqeron 100 Mbp are single *blocked* attempts.  `bench.tsv` column `repeats` states the real count.
* **Query timing.** untimed warm-up pass over the full query set (JIT/tier-1, caches), then full passes repeated until >= 1 s; throughput = queries / mean pass time.  Query sets are fixed files (`queries/`), 100 000 present 20-mers, 100 000 random 20-mers, 1 000 random 8-mers (`locate8`).
  Locate sums every returned position (so lists are really materialised); count for Go is `len(Lookup(q, -1))` (Go has no count call).  Per-query latency uses one clock read per query (adds ~25-40 ns).
* **Threads.** each thread runs the *whole* query set `reps` times against one shared index (same per-thread working set at 1 and 4 threads); throughput = threads x queries x reps / wall.
* **Memory.** `peak RSS` = `/proc/self/status VmHWM` (cross-checked with `wait4` `ru_maxrss`, column `ext_peak_rss`).  RSS of the persistent tree counts file-backed mmap pages that the kernel may drop, so the anonymous part (`RssAnon`, sampled every 20 ms) is reported separately.
  Seqeron in-memory bytes/char = `GC.GetTotalMemory(true)` delta around `Build` (text string excluded); Go bytes/char = after-GC `HeapAlloc` delta of `suffixarray.New` (the Index references the caller's text); rust-bio bytes/char is analytic (8n SA + n BWT + (n/3+1)*sigma*8 Occ), the RssAnon growth is in `rss_anon_growth`.
* **Cache state.** page cache dropped (`echo 3 > /proc/sys/vm/drop_caches`) only for the `cold_load` rows; all other query rows run on a warm page cache.
* **Alphabets.** Seqeron indexes UTF-16 `char`s (persistent file stores ASCII at 1 B/char); the other tools index bytes.  `prot_1m` (20 letters) was run for every tool that accepts bytes.
* **Seqeron settings.** Release build, .NET 10.0.403 (runtime 10.0.14), workstation non-concurrent GC, TieredPGO on, `ISuffixTree` string overloads (`CountOccurrences(string)` etc.).  Library code was built from a snapshot of commit 242318429a5e.
* **sdsl `csa_sada` is not safe for concurrent reads in our harness**: 4 threads sharing one instance returned a wrong count checksum (`checksum_ok = 0`) and then crashed with `std::bad_alloc` in the locate test (first run crashed the process; the harness was changed to stop after the first failed concurrent test).  Its 4-thread count number is therefore **not valid**; the other sdsl structures passed.  We did not dig into the cause.
* **libsais** rows use the serial `libsais()`; the `build_omp2` / `build_omp4` rows use `libsais_omp`.  libdivsufsort 2.0.1 serial.  Both with `-O3`, no `-march=native`.
* **Go count/locate** includes the cost of allocating the returned slice (`Lookup` copies); suffix-array C harnesses return slice indices without allocation, which is part of their locate advantage.
* **Randomness / generality.** inputs are random or lightly mutated; real genomes have different repeat structure.  `rep_10m` (100 copies of a 100 kbp block, 1 % substitutions) is the only repeat-rich input.
* **Not measured:** out-of-core (index > RAM) behaviour; Seqeron `FindAllOccurrences` span overloads; any Seqeron operation not listed above (approximate search, reverse complement, k-common substrings...); libsais64 / > 2^31 characters; sdsl `cst` MEM/MUM (no API).
* **Anything marked `NA` in `bench.tsv`** (value `NA`, repeats 0) is an operation the tool's own API does not offer; it is deliberately visible instead of absent.

## 7. Status per requested subject

| subject | status | where |
|---|---|---|
| Seqeron in-memory (1, 10 Mbp, rep 10 Mbp, protein 1 Mbp) | MEASURED, 3 repeats | `bench.tsv` tool `seqeron-st-inmem` |
| Seqeron persistent (1, 10, 20, 50 Mbp, rep, protein) | MEASURED, 3 repeats; 100 Mbp BLOCKED (600 s cap) | `seqeron-st-persistent` |
| Go `index/suffixarray` (1, 10, 20, 50 Mbp, rep, protein) | MEASURED, 3 repeats | `go-suffixarray` |
| libsais, libdivsufsort | MEASURED (built from source with cmake) | `libsais`, `libdivsufsort` |
| sdsl-lite CSA/FM-index + CST | MEASURED (1 repeat for RRR/CST at 10 Mbp) | `sdsl-*` |
| rust-bio | MEASURED (crate fetched from crates.io) | `rust-bio` |
| pydivsufsort, Python `suffix-trees` | MEASURED; suffix-trees BLOCKED at 20 Mbp (MemoryError at 12 GB cap) | `pydivsufsort`, `suffix-trees-py` |
| BWA, Bowtie2 | MEASURED: index build time / size only | `bwa-index`, `bowtie2-build` |

Hand-off scripts (to re-run or to extend on a bigger machine): `HANDOFF.md`.

## 8. Files

`bench.tsv` (final rows), `runs/raw*.tsv` (every repeat), `crosscheck.md`, `ENV.txt`, `SOURCES.txt`, `MANIFEST.sha256` (input hashes), `gen/` (input generator), `lib/` (driver, assembler, report generator),
`seqeron/SeqeronBench` (C# harness) and `seqeron/snapshot` (library snapshot), `go/` (Go harness), `competitors/` (C, C++, Rust, Python harnesses + run scripts), `run_main.sh run_all.sh run_stretch.sh run_100m.sh` (exact Seqeron/Go command lines).

Exact commands: `python3 gen/gen_inputs.py bench`; `dotnet build seqeron/SeqeronBench -c Release -o seqeron/out`; `(cd go && go build -o gobench .)` (GOTOOLCHAIN=local GOPROXY=off GOFLAGS=-mod=mod);
`./run_all.sh` (Go 1-50 Mbp, Seqeron in-memory, Seqeron persistent 1-10 Mbp), `./run_stretch.sh` (persistent 20/50 Mbp), `./run_100m.sh`; `competitors/run_libsais_libdivsufsort.sh`, `run_sdsl.sh`, `run_rustbio.sh`, `run_python.sh`, `run_aligners.sh`;
`python3 lib/assemble.py bench` -> `bench.tsv`; `python3 lib/crosscheck.py bench.tsv`; `python3 lib/make_report.py bench` -> this file.
