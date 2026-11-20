# minimap2_chains

minimap2 seeding and chaining.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `minimap2_chains` |
| **Method ID** | `Minimap2Chaining.MapChains` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The seeding and chaining stage of minimap2 2.26 (450 related pairs: the seed array of `--print-seeds` equal anchor by anchor, the chains of `--print-chains` equal, chain scores = the PAF s1 / cm tags; the chaining dynamic program mg_lchain_dp, its backtracking and compaction are ported line by line): the (k, w)-minimizers of the targets and of the query (mm_sketch, k <= 28, w < 256; symmetric k-mers skipped, A/C/G/T/U only) are indexed, every query minimizer that occurs in a target gives anchors on the forward or reverse strand (target / query position of the LAST base of the seed, sorted by minimap2's in-place radix sort), and the anchors are chained: the score of a step is min(min(dq, dr), span) minus floor(0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)); a predecessor must lie within maxGap on both sequences and within bandwidth of the diagonal, at most 5000 predecessors are scanned per anchor and the scan stops after 25 non-improving ones; chains are cut where the score drops by more than the bandwidth below its peak and kept when they reach minChainScore and minAnchors. Returns the chains (score, target, strand, anchors in increasing order) ordered by the target position of their first anchor. Frequent minimizers (A4-7): a minimizer that occurs more than mid_occ times in the targets (mid_occ = max(minMidOcc, the (1 - midOccFraction) quantile of the occurrence counts + 1), capped by maxMidOcc; midOccFraction >= 1 is the limit itself, as minimap2 -f) is dropped except for the least frequent ones of a streak of frequent minimizers (mm_seed_select: one per occDist query bases, never above 4095 occurrences); query minimizers that occur more than mid_occ times in the query and in more than queryOccFraction of its minimizers are dropped first; repeatLength = the query bases covered by dropped minimizers (the PAF rl tag). Also ported: homopolymerCompressed (-H), sdustThreshold (-T, the sdust low-complexity masking of the query minimizers), allVsAll (-X: with targetNames / queryName, hits on the diagonal of a same-named, same-length target and targets that sort before the query are skipped; all chains are reported, no long join), forwardOnly / reverseOnly. mode = chains (default, as before: every chain after the dynamic program, ordered by the first anchor) or paf: the whole mapping up to the end of the chaining stage (mm_map_frag without -c): the long-join re-chaining of -r b,B (bandwidthLong > bandwidth: mg_lchain_rmq, a range-minimum-query AVL tree, when the first chains do not cover the query), rmq chaining (--rmq), the primary / secondary marking (mm_set_parent: maskLevel, maskLength, hardMaskLevel; mm_select_sub: primaryRatio -p, bestN -N, strand retention), the divergence estimate dv (mm_est_err) and the mapping quality, i.e. every field of the PAF lines (tp P / S, cm, s1, s2, mapq, dv, rl) of minimap2 without -c (6 000+ oracle cases of the installed program: seeds in six groups of 330 cases, PAF lines in four groups of 330). Not part of this tool: the base-level alignment minimap2 -c adds (ksw2 gap filling and extension, registered A3-9c2), splice mode and the short-read / paired-end modes. Options are validated like the program (A4-4): k outside 1-28 is refused (minimap2: "-k and -w must be positive", assertion k <= 28 in mm_sketch); every other value minimap2 accepts is accepted and gives the program's chains - a negative bandwidth leaves no chain, a negative gap limit means the bandwidth, minAnchors 0 / -1 and a negative minChainScore are no limit (measured on 30 option values, see option_boundary_oracle.json). A query / target pair whose minimizers have more than 8 388 608 anchors in all (SequenceAligner.MinimapMaxAnchors, e.g. poly-A or a short tandem repeat against itself) is refused with an ArgumentException before the anchors are stored; minimap2 itself dies at the allocation or drops the frequent minimizers (not ported).

## Core Documentation Reference

- Source: [Minimap2Chaining.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Minimap2Chaining.cs), [Minimap2Chaining.Seeds.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Minimap2Chaining.Seeds.cs), [Minimap2Chaining.Hits.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Minimap2Chaining.Hits.cs), [Minimap2Chaining.Rmq.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Minimap2Chaining.Rmq.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [generate_minimap2_chain_oracle.py](../../../../scripts/alignment/generate_minimap2_chain_oracle.py), [generate_minimap2_seed_oracle.py](../../../../scripts/alignment/generate_minimap2_seed_oracle.py) (frequent minimizers, -f / -U, -H, -T, -X), [generate_minimap2_map_oracle.py](../../../../scripts/alignment/generate_minimap2_map_oracle.py) (PAF primary / secondary, long join, rmq)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `targets` | array<string> | Yes | Target (reference) sequences (A/C/G/T/U; other symbols break k-mers), at least one. |
| `query` | string | Yes | Query sequence. |
| `k` | integer | No | k-mer size (minimap2 -k, 1-28, default 15). Default: `15`. |
| `w` | integer | No | Minimizer window size (minimap2 -w, 1-255, default 10). Default: `10`. |
| `maxGap` | integer | No | Largest distance between two anchors of a chain on both sequences (minimap2 -g, default 5000). Default: `5000`. |
| `bandwidth` | integer | No | Chaining bandwidth (minimap2 -r, first value, default 500). Default: `500`. |
| `minChainScore` | integer | No | Smallest chain score (minimap2 -m, default 40). Default: `40`. |
| `minAnchors` | integer | No | Fewest anchors of a chain (minimap2 -n, default 3). Default: `3`. |
| `mode` | string | No | chains (default: every chain after the dynamic program, ordered by the first anchor) | paf (the mapping up to the end of the chaining stage: long-join re-chaining, primary / secondary marking, dv, mapq; the hits are returned in the PAF order). Default: `chains`. |
| `homopolymerCompressed` | boolean | No | Homopolymer-compressed k-mers (minimap2 -H). Default: `false`. |
| `sdustThreshold` | integer | No | sdust threshold masking low-complexity query minimizers (minimap2 -T, hidden option; 0 = off, default). Default: `0`. |
| `midOccFraction` | number | No | minimap2 -f: below 1 the fraction of repeated minimizers to filter (default 0.0002), otherwise the occurrence limit mid_occ itself. Default: `0.0002`. |
| `minMidOcc` | integer | No | minimap2 -U lower bound of mid_occ (default 10). Default: `10`. |
| `maxMidOcc` | integer | No | minimap2 -U upper bound of mid_occ (default 1000000). Default: `1000000`. |
| `queryOccFraction` | number | No | minimap2 --q-occ-frac (default 0.01). Default: `0.01`. |
| `occDist` | integer | No | minimap2 -e: occ_dist, one frequent minimizer per this many query bases is kept (default 500; 0 = none). Default: `500`. |
| `allVsAll` | boolean | No | All-vs-all mode (minimap2 -X): needs targetNames and queryName; all chains are reported without primary selection and without long join. Default: `false`. |
| `targetNames` | array<string> | No | Names of the targets (for allVsAll), one per target. |
| `queryName` | string | No | Name of the query (for allVsAll; in paf mode its hash orders chains of equal score). |
| `forwardOnly` | boolean | No | Only hits on the forward strand (minimap2 --for-only). Default: `false`. |
| `reverseOnly` | boolean | No | Only hits on the reverse strand (minimap2 --rev-only). Default: `false`. |
| `bandwidthLong` | integer | No | Long-join bandwidth (minimap2 -r second value, default 20000; the re-chaining runs only when it exceeds bandwidth, paf mode). Default: `20000`. |
| `maxChainSkip` | integer | No | Stop the search for a predecessor after this many non-improving anchors (minimap2 --max-chain-skip, default 25). Default: `25`. |
| `maxChainIterations` | integer | No | Largest number of predecessors examined per anchor (minimap2 --max-chain-iter, default 5000). Default: `5000`. |
| `chainGapScale` | number | No | Gap penalty scale (minimap2 --chain-gap-scale, default 0.8). Default: `0.8`. |
| `chainSkipScale` | number | No | Skip penalty scale (minimap2 --chain-skip-scale, default 0). Default: `0`. |
| `bestN` | integer | No | paf mode: retain at most this many secondary hits (minimap2 -N, default 5). Default: `5`. |
| `primaryRatio` | number | No | paf mode: minimum secondary-to-primary score ratio (minimap2 -p, default 0.8; 0 = no selection). Default: `0.8`. |
| `maskLevel` | number | No | paf mode: minimap2 --mask-level (default 0.5). Default: `0.5`. |
| `maskLength` | integer | No | paf mode: minimap2 --mask-len (default none). Default: `2147483647`. |
| `hardMaskLevel` | boolean | No | paf mode: minimap2 --hard-mask-level. Default: `false`. |
| `noLongJoin` | boolean | No | paf mode: no long-join re-chaining (minimap2 --no-long-join). Default: `false`. |
| `rmq` | boolean | No | paf mode: chain with the RMQ algorithm (minimap2 --rmq=yes). Default: `false`. |
| `outputSecondary` | boolean | No | paf mode: return the secondary hits too (default true; false = minimap2 --secondary=no). Default: `true`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `chains` | array | Chains `{ score, target, strand, anchors: [{ target, strand, targetEnd, queryEnd, span }] }` ordered by the target position of the first anchor. |
| `seedCount` | integer | Number of seed hits (anchors) the chains were built from. |
| `hits` | array | paf mode only: the hits in the order of the PAF lines of minimap2 without -c `{ id, parent, primary, target, strand, queryStart, queryEnd, targetStart, targetEnd, matchLength, blockLength, mapQ, score (s1), secondaryScore (s2), anchorCount (cm), divergence (dv), anchors }`; `chains` then lists the same hits as chains. |
| `repeatLength` | integer | Query bases covered by frequent minimizers that were dropped (the PAF rl tag). |
| `midOcc` | integer | The occurrence limit mid_occ that was applied. |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one target sequence is required |
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: A 100 base match, a 3 base substitution block and 79 more bases (k = 11, w = 5)

**User Prompt:**
> Chain the minimap2 seeds of the query against the reference with k 11 and w 5.

**Expected Tool Call:**
```json
{"tool": "minimap2_chains", "arguments": {"targets": ["TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCT"], "query": "CCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCTTTCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAA", "k": 11, "w": 5}}
```

**Response:**
```json
{"chains": [{"score": 172, "target": 0, "strand": "forward", "anchors": "53 anchors, the first at target 31 / query 11, the last at target 199 / query 181, span 11"}], "seedCount": 53}
```

`minimap2 -k 11 -w 5` prints the same 53 seeds and one chain of 53 anchors (PAF `cm:i:53 s1:i:172`); the anchors are abbreviated here.

## See Also

- [local_align_top_n](local_align_top_n.md) - non-intersecting local alignments
- [ksw2_extend](ksw2_extend.md) - the ksw2 extension minimap2 uses between and after the chain
- [nucmer_clusters](../analysis/nucmer_clusters.md) - the nucmer anchor clustering
