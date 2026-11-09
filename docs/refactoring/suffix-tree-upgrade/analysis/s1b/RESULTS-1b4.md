# Work package 1b.4: results (production navigator behind the generic algorithms)

Raw rows: `nav_raw.tsv` (fresh process per run, run scripts `run_nav.sh` / `run_nav_all.sh`, aggregation `aggregate_nav.py`); the first set of MEM / MUM / anchors / LCS rows, taken with the Ukkonen tree built BEFORE the ESA calls, is kept in `nav_raw_discarded_ukkonen_first.tsv` and not used (see below).

Conditions: this machine (4 cores, 480 MiB L3), `Release`, defaults of the runtime (workstation non-concurrent GC, tiered compilation and dynamic PGO on), build servers shut down, nothing else running, load average about 1 (the harness itself); medians over 5 fresh processes with [min .. max]; the query is `<name>.memq100k` of the benchmark generator (100 kbp, 1 % substitutions), minLength 20; the first call of a process is the cold call (JIT included, the navigation pack already built), `warm` is the median of calls 2..10, `steady` of calls 6..10. The ESA calls run BEFORE the Ukkonen tree is built (probable cause, not checked by counting collections: its 1.9 GB heap of small objects makes the gen-2 collections that the large result lists provoke expensive; with the tree built first the ESA took 0.182 s instead of 0.153 s on `rep_10m` MEM).

## dna_10m: MEM / MUM / anchors of the 100 kbp query, seconds

| operation | ESA cold | ESA warm | ESA steady | Ukkonen cold | Ukkonen warm | Ukkonen steady | ESA / Ukkonen (cold, steady) | results |
|---|---|---|---|---|---|---|---|---|
| mem | 0.109 [0.099 .. 0.111] | 0.067 [0.054 .. 0.081] | 0.066 [0.053 .. 0.078] | 0.068 [0.066 .. 0.073] | 0.040 [0.038 .. 0.042] | 0.037 [0.034 .. 0.040] | 1.60x, 1.77x | 817 |
| mum | 0.075 [0.073 .. 0.097] | 0.046 [0.044 .. 0.048] | 0.045 [0.043 .. 0.048] | 0.053 [0.051 .. 0.056] | 0.035 [0.032 .. 0.043] | 0.034 [0.032 .. 0.037] | 1.41x, 1.33x | 816 |
| anchors | 0.076 [0.071 .. 0.076] | 0.044 [0.043 .. 0.047] | 0.044 [0.043 .. 0.048] | 0.050 [0.049 .. 0.051] | 0.031 [0.030 .. 0.035] | 0.031 [0.030 .. 0.034] | 1.52x, 1.43x | 816 |

Answers equal to the Ukkonen tree's: MEM 5/5, MUM 5/5, anchors 5/5 processes (MEM and MUM lists exactly, anchors as (query position, length) pairs). ESA build 0.999 [0.898 .. 1.028] s, navigation pack 0.0975 [0.0940 .. 0.0994] s (pyramid alone, warm: 0.0183 [0.0162 .. 0.0184] s), 4.129 B/char (ISA 4.000 + pyramid 0.129, 5 levels), Ukkonen build 12.9 [12.5 .. 13.8] s, peak RSS 1.90 GB with both structures in the process.

## rep_10m: MEM / MUM / anchors of the 100 kbp query, seconds

| operation | ESA cold | ESA warm | ESA steady | Ukkonen cold | Ukkonen warm | Ukkonen steady | ESA / Ukkonen (cold, steady) | results |
|---|---|---|---|---|---|---|---|---|
| mem | 0.195 [0.188 .. 0.250] | 0.156 [0.150 .. 0.157] | 0.153 [0.151 .. 0.157] | 0.535 [0.512 .. 0.653] | 0.382 [0.371 .. 0.404] | 0.387 [0.378 .. 0.448] | 0.36x, 0.40x | 160979 |
| mum | 0.064 [0.061 .. 0.065] | 0.038 [0.034 .. 0.041] | 0.038 [0.032 .. 0.041] | 0.034 [0.033 .. 0.036] | 0.020 [0.014 .. 0.022] | 0.019 [0.013 .. 0.021] | 1.90x, 1.94x | 682 |
| anchors | 0.059 [0.056 .. 0.062] | 0.034 [0.030 .. 0.042] | 0.034 [0.029 .. 0.039] | 0.033 [0.032 .. 0.038] | 0.021 [0.018 .. 0.024] | 0.019 [0.017 .. 0.024] | 1.80x, 1.74x | 715 |

Answers equal to the Ukkonen tree's: MEM 5/5, MUM 5/5, anchors 5/5 processes (MEM and MUM lists exactly, anchors as (query position, length) pairs). ESA build 1.126 [1.049 .. 1.157] s, navigation pack 0.1156 [0.1104 .. 0.1188] s (pyramid alone, warm: 0.0238 [0.0228 .. 0.0290] s), 4.129 B/char (ISA 4.000 + pyramid 0.129, 5 levels), Ukkonen build 8.8 [8.4 .. 9.4] s, peak RSS 2.29 GB with both structures in the process.

## LCS family, two 1 Mbp texts (dna_1m against a deterministic second text of 1 Mbp), seconds

| second text | operation | ESA cold | ESA steady | Ukkonen cold | Ukkonen steady | length |
|---|---|---|---|---|---|---|
| mutated copy (1 % substitutions, 0.1 % indels) | lcs_info | 0.261 [0.249 .. 0.263] | 0.223 [0.209 .. 0.254] | 0.266 [0.254 .. 0.273] | 0.254 [0.238 .. 0.258] | 1078 chars |
| mutated copy (1 % substitutions, 0.1 % indels) | lcs_all | 0.223 [0.212 .. 0.256] | 0.230 [0.220 .. 0.246] | 0.243 [0.230 .. 0.255] | 0.255 [0.249 .. 0.265] | 1078 chars |
| mutated copy (1 % substitutions, 0.1 % indels) | lcs_distinct | 0.230 [0.220 .. 0.246] | 0.228 [0.212 .. 0.250] | 0.261 [0.226 .. 0.265] | 0.250 [0.238 .. 0.265] | 1 ties |
| independent random text | lcs_info | 0.326 [0.321 .. 0.340] | 0.314 [0.297 .. 0.318] | 0.349 [0.331 .. 0.361] | 0.347 [0.334 .. 0.365] | 21 chars |
| independent random text | lcs_all | 0.310 [0.298 .. 0.326] | 0.306 [0.297 .. 0.318] | 0.351 [0.316 .. 0.355] | 0.348 [0.329 .. 0.354] | 21 chars |
| independent random text | lcs_distinct | 0.304 [0.299 .. 0.351] | 0.299 [0.294 .. 0.305] | 0.362 [0.344 .. 0.370] | 0.341 [0.330 .. 0.351] | 1 ties |

## Production link kernel against the prototype (c2), ns per link, same nodes, paired halves

| input | production | prototype (c2, ISA, fanout 32, plus the SA[Lb] read that completes the handle) | ratio | differences |
|---|---|---|---|---|
| dna_10m | 428 [420 .. 460] | 511 [500 .. 558] | 0.817 [0.803 .. 0.867] | 0 |
| rep_10m | 364 [355 .. 397] | 500 [485 .. 543] | 0.694 [0.684 .. 0.727] | 0 |


## dna_50m, MEM (information only, ESA alone)

cold 0.126 [0.125 .. 0.140] s, steady 0.082 [0.080 .. 0.089] s, 802 MEMs; ESA build 5.36 [5.13 .. 5.40] s, pack 0.647 [0.640 .. 0.647] s; peak RSS 0.87 GB (no Ukkonen tree was built at 50 Mbp: about 9 GB estimated from 1.9 GB at 10 Mbp on a 15 GB machine, and the 1b.4a attempt at this size failed).


## Where the MEM time goes (`nav --op split`, 3 fresh processes per input)

**dna_10m (817 MEMs)**: plain navigator, steady call 61.4 ms, cold 103 ms, 33 KiB allocated per call.

| member | calls per 100,000 query starts | replay (ns per call, ms, % of the plain call) | timed run (net of the timer floor, ms) |
|---|---|---|---|
| GetSuffixLink | 200,000 | 93 ns, 18.6 ms, 30 % | 37.8 |
| TryGetChild (after a link / otherwise) | 359,193 / 9,226 | 87 ns, 32.1 ms, 52 % | 40.5 |
| GetChildren | 1 | 89 ns, 0.0 ms, 0 % | 0.0 |
| CollectLeaves | 81,836 | 33 ns, 2.7 ms, 4 % | 0.9 |
| GetEdgeSymbol | 273,213 | 14 ns, 3.7 ms, 6 % | - |
| LengthOf | 633,828 | 9 ns, 5.4 ms, 9 % | - |
| IsNull / IsRoot | 745,171 / 200,000 | inlined compares | - |

Replay total 64.0 ms (104 % of the plain call: the members overlap in the processor, so the percentages are not additive); timed run: the timed members add up to 79.9 ms of a timed run of 117 ms (the timers serialise the pipeline: only the shares of the timed run, not its absolute times, mean anything).

**rep_10m (160,979 MEMs)**: plain navigator, steady call 147.5 ms, cold 183 ms, 8196 KiB allocated per call.

| member | calls per 100,000 query starts | replay (ns per call, ms, % of the plain call) | timed run (net of the timer floor, ms) |
|---|---|---|---|
| GetSuffixLink | 200,000 | 81 ns, 16.2 ms, 11 % | 51.8 |
| TryGetChild (after a link / otherwise) | 214,597 / 1,448,512 | 36 ns, 60.5 ms, 41 % | 70.3 |
| GetChildren | 1,392,772 | 34 ns, 47.1 ms, 32 % | 39.5 |
| CollectLeaves | 1,765,823 | 18 ns, 32.7 ms, 22 % | 20.2 |
| GetEdgeSymbol | 3,288,984 | 18 ns, 58.7 ms, 40 % | - |
| LengthOf | 1,963,249 | 10 ns, 19.1 ms, 13 % | - |
| IsNull / IsRoot | 641,521 / 200,000 | inlined compares | - |

Replay total 234.9 ms (159 % of the plain call: the members overlap in the processor, so the percentages are not additive); timed run: the timed members add up to 182.1 ms of a timed run of 490 ms (the timers serialise the pipeline: only the shares of the timed run, not its absolute times, mean anything).

Reading: on `dna_10m` (a text without repeats, so the 100,000 query starts produce 200,000 links, 368,000 child lookups of which 359,000 follow a link, and almost no subtree enumeration) the two heavy members, the link kernel (about 100 ns per link in place, 430-510 ns on a random node) and the child lookups of the rescans (about 90 ns each, two dependent reads per candidate child), account for about 80-90 % of the call; the rest is the cheap members (about 15 % by replay) and the generic algorithm, whose own cost and allocations (33 KiB per call) are below the resolution of the method. On `rep_10m` (3 % substitutions in a repeated block, 160,979 MEMs) the cost moves to the enumeration of the sibling subtrees that the MEM algorithm walks (1.39 million `GetChildren`, 1.45 million child lookups off the rescan path, 1.77 million `CollectLeaves` with 5.5 million positions, 3.3 million `GetEdgeSymbol`) and the link kernel falls to about 11 % (replay) of the call. The 8.4 MB allocated per call are the result list and the lists of the generic algorithm.

## Experiments that were run and not adopted

* **Lazy `Pos` of the link target** (the handle of the target carries `Pos = -1`, the navigator reads `SA[Lb]` when an edge symbol or a leaf position is needed): saves the one `SA[Lb]` read per link, which costs about 120 ns on a random node (prototype with and without it: 600 and 480 ns), but MEM on `dna_10m` was unchanged (4 paired runs: cold 0.100 s against 0.100 s, steady 0.057 s against 0.055 s), so it was dropped (in place the read apparently hits the cache; not examined further).
* **Software prefetch in the ISA scatter** (16 ranks ahead): pack build 0.098 s without, 0.100 s with: no gain, dropped.
* **One endpoint instead of two for the link** (`Psi(l)` only: once the value `t = depth - 1` is known, the maximal run of LCP entries `>= t` around any rank of the target is the target, so `Psi(r)` is redundant; mutant P10 survived for exactly this reason): two dependent random reads fewer, but MEM on `dna_10m` was 5-8 % SLOWER in 5 paired runs (cold 0.113-0.123 s against 0.109-0.117 s, steady 0.078-0.083 s against 0.073-0.076 s), the two ISA reads being independent and overlapping in the memory system. Dropped; the LCA form of the ADR stays.
* **`AggressiveOptimization` on the hot members** (adopted): cold MEM on `dna_10m` 0.145 s -> 0.098 s (6 processes each: 0.142-0.162 s against 0.093-0.107 s), steady unchanged (0.063 s against 0.065 s); on the navigator members alone (without the pyramid and the pack) the cold time is 0.116 s. Without it the first call runs the pointer code unoptimised until the tiering delay of the runtime has passed.
* **Pointer reader against the safe indexer**: per link, with equal work, the production kernel takes 0.82x (`dna_10m`) and 0.69x (`rep_10m`) of the prototype's time, so the escape-aware pointer reader pays for itself; without the completing `SA[Lb]` read the prototype takes 0.90-0.93x of the production time on `dna_10m` and 1.03-1.07x on `rep_10m` (two single-process runs each, not in the raw file), i.e. the kernel itself is bound by a handful of dependent cache misses per link and the reader matters at the 10-30 % level only.
