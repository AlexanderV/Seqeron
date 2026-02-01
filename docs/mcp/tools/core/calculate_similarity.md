# calculate_similarity

Calculate similarity between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `calculate_similarity` |
| **Method ID** | `GenomicAnalyzer.CalculateSimilarity` |
| **Version** | 1.0.0 |

## Description

Calculates similarity between two DNA sequences using k-mer Jaccard index. Returns value between 0 (no similarity) and 1 (identical k-mer composition).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L238](../../../../src/Seqeron/Seqeron.Genomics/GenomicAnalyzer.cs#L238)

## Input/Output

**Input:** `{ sequence1: string, sequence2: string, kmerSize?: integer }`

**Output:** `{ similarity: number }` (0-1 scale)

## Example

**User Prompt:**
> How similar are ATGCATGC and ATGCATTT?

**Tool Call:**
```json
{ "tool": "calculate_similarity", "arguments": { "sequence1": "ATGCATGC", "sequence2": "ATGCATTT" } }
```

**Response:**
```json
{ "similarity": 0.75 }
```
