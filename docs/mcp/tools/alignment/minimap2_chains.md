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

The seeding and chaining stage of minimap2 2.26 (450 related pairs: the seed array of `--print-seeds` equal anchor by anchor, the chains of `--print-chains` equal, chain scores = the PAF s1 / cm tags; the chaining dynamic program mg_lchain_dp, its backtracking and compaction are ported line by line): the (k, w)-minimizers of the targets and of the query (mm_sketch, k <= 28, w < 256; symmetric k-mers skipped, A/C/G/T/U only) are indexed, every query minimizer that occurs in a target gives anchors on the forward or reverse strand (target / query position of the LAST base of the seed, sorted by minimap2's in-place radix sort), and the anchors are chained: the score of a step is min(min(dq, dr), span) minus floor(0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)); a predecessor must lie within maxGap on both sequences and within bandwidth of the diagonal, at most 5000 predecessors are scanned per anchor and the scan stops after 25 non-improving ones; chains are cut where the score drops by more than the bandwidth below its peak and kept when they reach minChainScore and minAnchors. Returns the chains (score, target, strand, anchors in increasing order) ordered by the target position of their first anchor. Not part of this tool: the base-level alignment minimap2 -c adds (ksw2 gap filling and extension), the secondary / primary selection, and the handling of frequent minimizers (a minimizer of the query that occurs more than mid_occ = max(10, 0.9998 quantile + 1) times in the targets is an error here, as is homopolymer compression).

## Core Documentation Reference

- Source: [Minimap2Chaining.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/Minimap2Chaining.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [generate_minimap2_chain_oracle.py](../../../../scripts/alignment/generate_minimap2_chain_oracle.py)

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

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `chains` | array | Chains `{ score, target, strand, anchors: [{ target, strand, targetEnd, queryEnd, span }] }` ordered by the target position of the first anchor. |
| `seedCount` | integer | Number of seed hits (anchors) the chains were built from. |

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
- [nucmer_clusters](nucmer_clusters.md) - the nucmer anchor clustering
