---
type: concept
title: "Prokaryotic gene prediction & Shine-Dalgarno RBS detection"
tags: [annotation, algorithm]
mcp_tools:
  - predict_genes
sources:
  - docs/Validation/reports/ANNOT-GENE-001.md
source_commit: fc49476951d4c26bf663220b231007d8327e59cf
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-gene-001-report
      evidence: "Validation report ANNOT-GENE-001 validates GenomeAnnotator.PredictGenes + FindRibosomeBindingSitesBothStrands (Stage A PASS / Stage B PASS / CLEAN); the test-unit-registry tracks the ANNOT-GENE-001 unit."
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:open-reading-frame-detection
      source: annot-gene-001-report
      evidence: "ANNOT-GENE-001 validates the annotation-layer GenomeAnnotator ORF/gene finder (prokaryotic starts, both strands, CDS labelling), which the open-reading-frame-detection page records as deliberately NOT contract-equivalent to the Analysis-layer GenomicAnalyzer.FindOpenReadingFrames — the two are alternative entry points for ORF finding."
      confidence: high
      status: current
---

# Prokaryotic gene prediction & Shine-Dalgarno RBS detection

**Prokaryotic gene prediction** here is ORF-based gene finding in the annotation layer:
`GenomeAnnotator.PredictGenes(dna, minOrfLength, prefix)` scans all three frames on **both
strands** for start→stop spans, filters by a minimum length, labels each as a CDS, and emits a
gene record whose `protein_length` **excludes the stop codon** (`ProteinSequence.TrimEnd('*').Length`,
i.e. `(End−Start)/3 − 1`). Its companion is **Shine-Dalgarno ribosome-binding-site detection**,
which locates the SD motif upstream of a start codon at a valid **aligned spacing**. Validated
under test unit **ANNOT-GENE-001** (Stage A PASS / Stage B PASS / State CLEAN, independently
re-derived and mutation-checked); the validation record is [[annot-gene-001-report]],
[[test-unit-registry]] tracks the unit, and [[validation-ledger]] is the live status board.

## Shine-Dalgarno RBS model

The SD is a feature of the **mRNA**: it base-pairs with the anti-SD 3' tail of the 16S rRNA
(`5'-YACCUCCUUA-3'`, Y = pyrimidine), located ~8 bases upstream of the start codon. Seqeron
scans for the consensus `AGGAGG` and its biologically-justified substrings `GGAGG / AGGAG / GAGG /
AGGA` (all substrings of the consensus; T4 early genes are GAGG-dominated).

| Method | Behaviour |
|---|---|
| `FindRibosomeBindingSites(...)` | forward strand only (legacy; `.Where(o => !o.IsReverseComplement)`) |
| `FindRibosomeBindingSitesBothStrands(...)` | forward + reverse passes; reports `(position, sequence, score, strand)` |
| `ScanStrandForShineDalgarno(...)` | the single shared single-strand scan reused by both passes (no logic duplication) |

- **Aligned spacing** = SD 3' end → the A of the AUG: `distanceToStart = orf.Start − genomicPos −
  motif.Length`, accepted when in `[minDistance, maxDistance]` (`maxDistance` included, above it
  excluded). Score = `motif.Length / 6.0`.
- The default `[4,15]` window is a **declared permissive heuristic search window** — the method
  reports candidates, it does **not** score optimality. The literature optimum is a **5 nt** aligned
  spacing (Chen et al. 1994, *NAR* 22(23):4953-4957, PMID 7528374).
- Edge cases: empty/null input → empty; no ORF / no upstream window → no RBS (SD is defined
  relative to a start codon).

## Reverse-strand SD coordinate mapping (the correctness core)

Because a reverse-strand gene's mRNA is the reverse complement of the forward genomic strand, its
SD is found by scanning the **reverse complement**, not the forward strand. The reverse pass
un-maps reverse ORFs back into revComp space (`o with { Start = len−o.End, End = len−o.Start }`),
scans `revComp`, then maps a motif hit back to forward coordinates via:

```
forwardPosition = len − hit.position − motifLen        // GenomeAnnotator.cs:341
```

Worked oracle (ANNOT-GENE-001, independently re-derived): for
`S = C×10 + AGGAGG + C×8 + (ATG + A×297 + TAA)`, `|S| = 327`, genomic input `RC(S)` — AGGAGG at
revPos 10 maps to `forwardPos = 327 − 10 − 6 = 311`, where `RC(S)[311..317] = "CCTCCT"` (the anti-SD
complement on the forward strand) and the aligned spacing in revComp space is `24 − 10 − 6 = 8`. The
validator's mutation (`forwardPosition = hit.position`) failed R1/R3/R4, confirming the mapping is
locked.

## Scope and siblings

This is the **annotation-layer** prokaryotic gene/RBS finder. It is an [[open-reading-frame-detection|alternative]]
entry point to the Analysis-layer six-frame `GenomicAnalyzer.FindOpenReadingFrames` (ATG-only,
nucleotide `minLength`) — the two are deliberately **not contract-equivalent**; callers pick the
entry point. The static SD consensus string also appears as one row of the fixed
[[regulatory-element-detection]] catalog, but that scanner matches the bare `AGGAGG` anywhere and
does **not** model RBS-to-start aligned spacing or strand-aware coordinate mapping, which is what
this unit owns. Coding-potential scoring of a candidate ORF is a separate step
([[coding-potential-hexamer-score]]).
