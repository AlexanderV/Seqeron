---
type: concept
title: "Semi-global alignment (fitting / query-in-reference)"
tags: [alignment, algorithm]
sources:
  - docs/Evidence/ALIGN-SEMI-001-Evidence.md
  - docs/algorithms/Alignment/Semi_Global_Alignment.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-semi-001-evidence
      evidence: "Test Unit ID: ALIGN-SEMI-001 ... Algorithm: Semi-Global Alignment (Fitting / Query-in-Reference)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:global-alignment-needleman-wunsch
      source: align-semi-001-evidence
      evidence: "Hybrid method ... combines features of both global and local alignment ... uses the NW recurrence with linear gap cost but free end gaps in the reference (first row = 0), traceback from max of last row"
      confidence: high
      status: current
---

# Semi-global alignment (fitting / query-in-reference)

A **hybrid ("glocal", ends-free) alignment** that finds the best partial alignment of two
sequences: the shorter **query** is aligned end-to-end while the longer **reference** may
begin and/or end in free (unpenalized) gaps. Seqeron implements the **fitting** variant
(query-in-reference, Rosalind SIMS); its validation record is
[[align-semi-001-evidence|ALIGN-SEMI-001]].

## Where it sits in the alignment family

Semi-global is a distinct member of the family alongside global (Needleman–Wunsch) and local
(Smith–Waterman), selected by an `AlignmentType.SemiGlobal` mode. The family differs only in
**matrix border initialization** and **traceback start**:

| Mode | First row F(0,j) | First column F(i,0) | Traceback starts at |
|------|-------------------|---------------------|----------------------|
| Global (NW) | d·j | d·i | bottom-right F(m,n) |
| Local (SW) | 0 | 0 | global max (with 0 floor) |
| Fitting (this page) | **0** (free start gaps in ref) | d·i (query fully aligned) | **max of last row** max_j F(m,j) |

## Model

- **Recurrence:** the plain NW recurrence F(i,j) = max( F(i−1,j−1)+S(aᵢ,bⱼ), F(i−1,j)+d,
  F(i,j−1)+d ) — **no zero floor** (unlike Smith–Waterman), so scores may be negative.
- **Free end gaps:** first row is zero (reference may skip a leading prefix at no cost); the
  query column keeps the linear gap penalty *d* so the query is consumed in full.
- **Traceback:** starts from the maximum cell in the **last row** (best reference endpoint),
  runs back to the top row (full query coverage); reference bases past the endpoint are
  appended as unpenalized trailing gaps.
- **Complexity:** O(m·n) time and space for query length m, reference length n.

## Variant chosen

The semi-global family has several members that differ in *which* ends are free:

| Variant | Free end gaps | Use case | Rosalind |
|---------|---------------|----------|----------|
| **Query-in-reference (fitting)** — *implemented* | start & end of reference | short-read / primer mapping | SIMS |
| Overlap | end of seq1, start of seq2 | assembly overlap | OAP |
| Full semi-global | all four ends | general substring finding | SMGB |

Selecting the fitting variant is a deliberate design choice picking one well-defined member,
corresponding to Rosalind SIMS. Its defining invariants: the query is fully represented
(`RemoveGaps(aligned_query) == query`) and the aligned reference is a substring of the
reference, with score = max_j F(m,j).

See [[algorithm-validation-evidence]] for the evidence-artifact pattern behind this unit and
[[global-alignment-needleman-wunsch]] for the global mode it extends.
