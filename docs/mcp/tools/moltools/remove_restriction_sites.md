# remove_restriction_sites

Synonymously eliminate restriction sites from a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `remove_restriction_sites` |
| **Method ID** | `CodonOptimizer.RemoveRestrictionSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Rewrites codons to remove the listed restriction recognition sequences while preserving the encoded protein (only synonymous swaps are used). Site strings may be DNA or RNA and may contain IUPAC ambiguity codes; **both strands** are cleared (a non-palindromic site such as BsaI `GGTCTC` is also removed where its reverse complement `GAGACC` occurs). Each removed occurrence changes exactly one codon — the synonymous substitution with the highest usage frequency in the supplied table. The output is RNA-alphabet, trimmed to whole codons. Sites for which no synonymous alternative eliminates the match are left in place.

**Site limits.** Every site is at most **64 characters** long (`RestrictionAnalyzer.MaxRecognitionLength`, the same bound as [find_pattern_occurrences](find_pattern_occurrences.md) and the custom-enzyme tools; the scan costs sequence length × site length) and is made of the 15 IUPAC nucleotide codes (`U` is read as `T`); empty entries are skipped. A longer site, or one with any other character, is refused **whatever the coding sequence is** (a sub-codon sequence included — the check no longer depends on the data).

**Cost.** The time is linear in the sequence length: the scan restarts only where a changed codon can have created or exposed an occurrence, and a candidate is judged on a mutable buffer instead of a rebuilt sequence. Release timings for 96 000 nt (32 000 codons) before → after the fix (review 2026-09, B08 R14): the site `N` 10.2 s → 22 ms, 16 000 tandem `GAATTC` 20.4 s → 26 ms, a 64-character `N` site 10.3 s → 0.2 s; outputs are byte-identical (464 stored cases and 9 315 further random cases). Through [optimize_codons](optimize_codons.md) with a GC target (`BalancedOptimization`, `MinimizeSecondary`) each candidate is additionally scored against the whole sequence, so that path grows with sequence length × the number of occurrences removed.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L1321](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L1321)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `restriction_sites` | string[] | Yes | Recognition sequences to remove (≥1 entry; DNA or RNA, IUPAC codes allowed, **at most 64 characters each**; empty entries are skipped). |
| `target_organism` | object | Yes | Preset id or inline custom table. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard, identical to omitting it; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Synonymous codons are taken from this code's families, so the protein is preserved under it. One of 1–6, 9–16, 21–33. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `optimizedSequence` | string | RNA-alphabet sequence with sites removed where possible. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | At least one restriction site is required |
| 1003 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| 1004 | Restriction site '…' contains a character that is not an IUPAC nucleotide code (checked whatever the sequence is) |
| 1005 | A restriction site is N characters long; the longest supported is 64. (checked whatever the sequence is) |

## Examples

### Example 1: Remove EcoRI

`GAATTC` (Glu-Phe) with site `GAATTC` → `GAAUUU` (UUC→UUU, the highest-frequency E. coli substitution; same protein, no `GAAUUC`).

### Example 2: No site present → sequence returned unchanged (as RNA).

## See Also

- [optimize_codons](optimize_codons.md), [find_restriction_sites](find_restriction_sites.md)
