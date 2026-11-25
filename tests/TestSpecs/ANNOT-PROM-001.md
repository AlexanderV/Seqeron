# Test Specification: ANNOT-PROM-001

**Test Unit**: Promoter Detection
**Algorithm Group**: Annotation
**Canonical Method**: `GenomeAnnotator.FindPromoterMotifs(string dnaSequence[, int maxMismatches35, int maxMismatches10, int minSpacer, int maxSpacer])`
**Date**: 2026-11-25 (rewritten in review-2026-09 B11 F7; previous version 2026-03-05)
**Status**: Final

---

## Summary

Bacterial σ70 promoter boxes. The method delegates to the canonical `MotifFinder.FindSigma70Promoters`
(Harley & Reynolds 1987: −35 ≈ `TTGACA`, −10 ≈ `TATAAT`, spacer 15–21 bp, defaults ≤ 2 mismatches per box) and
reports each −35 / −10 box of every forward-strand pair once, ordered by position (−35 before −10 on ties), with
score = 1 − d/6 (d = Hamming distance to the consensus).

---

## Evidence Sources

| Source | Reference | Key Information |
|--------|-----------|-----------------|
| Harley & Reynolds (1987) | NAR 15(5):2343-2361, PMC340638 (abstract opened) | 263 promoters; −35 TTGACA / −10 TATAAT highly conserved; spacer 17 ± 1 bp in 92 % |
| Hawley & McClure (1983) | NAR 11(8):2237-2255, PMC325881 (abstract opened) | 112-promoter consensus compilation used as H&R's reference list |
| Wikipedia: Promoter (genetics) | opened 2026-11-25 | −10 TATAAT, −35 TTGACA, optimal spacing 17 bp; per-position % table has no citation (not used) |
| MotifFinder.FindSigma70Promoters | MOTIF-REGULATORY-001 (B05 F30) | Canonical pairing; spacer 15–21 default |
| lacUV5 control region | MotifFinder_Sigma70Promoters_Tests fixture | −35 TTTACA at 68, 18-bp spacer, −10 TATAAT at 92 |
| Independent Python oracle | docs/Evidence/ANNOT-PROM-001-Evidence.md | 400 random cases / 1,131 boxes, 0 mismatches |

---

## Must Tests

| ID | Test | Expected |
|----|------|----------|
| M01 | Consensus promoter `GG TTGACA C17 TATAAT GG` | exactly (2, −35, TTGACA, 1.0), (25, −10, TATAAT, 1.0) |
| M02/M03 | Isolated TTGACA / TATAAT / adjacent TTGACATATAAT | empty |
| M04 | Spacer 14 / 15 / 17 / 21 / 22 | empty / pair / pair / pair / empty |
| M05 | Overload spacer range (15–22 admits 22; 15–16 rejects 17) | as stated |
| M06 | 1 mismatch per box (TTTACA / TATGAT) | scores 5/6, 5/6 |
| M07 | 2-mismatch −35 (TGGCCA) default vs mm35 = 3 | default: 2 boxes; mm 3: adds GTGGCC (3/6) at 1 |
| M08 | Shared −35 box with −10 at spacer 15 and 21 | −35 reported once; default adds AATTAT (4/6) at 24 |
| M09 | lacUV5 published fixture | (1,0): (68 TTTACA 5/6), (92 TATAAT 1.0); defaults: + (54 TAGGCA 4/6), (77 TATGCT 4/6) |
| M10 | Random AT-rich 50-mer | 6 boxes exactly as the Python oracle |
| M11 | 50 random 120-mers | identical to boxes of `MotifFinder.FindSigma70Promoters` (delegation) |
| M12 | Promoter on the reverse strand | empty |
| M13 | Lowercase | same as uppercase |
| M14 | Empty / 4 nt / 11 nt | empty |
| M15 | null / `N` / mm 7 / minSpacer > maxSpacer | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException ×2 |

MCP (`find_promoter_motifs`): consensus promoter → both boxes; unpaired boxes → empty; equals `GenomeAnnotator`
on the oracle 50-mer (6 boxes); null/empty → ArgumentException.

---

## Deviations from Literature

None for the consensus definition. Forward strand only (tuple has no strand; use `FindSigma70Promoters(bothStrands:
true)`); promoter strength is out of scope (use `PredictSigma70Promoters`, Promoter Calculator v1.0).

---

## Test Files

- `Unit/Annotation/GenomeAnnotator_PromoterMotif_Tests.cs` (M01–M15)
- `Differential/AnnotationDifferentialTests.cs` (inline brute force, 5 cases)
- `Combinatorial/AnnotationCombinatorialTests.cs` (mismatches × spacer × box, 24 + 2)
- `Snapshots/AnnotationSnapshotTests.cs` (lacUV5)
- Heavy tier: `Properties/AnnotationProperties.cs`, `Metamorphic/AnnotationMetamorphicTests.cs`, `Fuzzing/AnnotationFuzzTests.cs`
- MCP: `Seqeron.Mcp.Annotation.Tests/FindPromoterMotifsTests.cs`
