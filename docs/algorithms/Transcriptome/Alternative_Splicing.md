# Alternative Splicing — Event Classification and Percent Spliced In (PSI)

| Field | Value |
|-------|-------|
| Algorithm Group | Transcriptome |
| Test Unit ID | TRANS-SPLICE-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Alternative splicing produces multiple mRNA isoforms from one gene by including or excluding specific RNA segments. This unit provides two operations: `CalculatePSI`, which quantifies the inclusion level (Percent Spliced In, Ψ) of a splicing event from inclusion and skipping read counts; and `DetectAlternativeSplicing`, which compares transcript isoforms of a gene and classifies their structural differences into the five canonical alternative-splicing classes. PSI is an exact arithmetic ratio; event classification is a deterministic, coordinate-driven rule set. [1][2][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A splicing event is defined relative to two isoforms of the same gene that differ in exon structure. The canonical taxonomy recognizes five event classes: skipped exon (SE, also "cassette exon"), intron retention (RI), alternative 5′ splice site (A5SS), alternative 3′ splice site (A3SS), and mutually exclusive exons (MXE). [1][4]

### 2.2 Core Model

**Percent Spliced In (unnormalized):** Ψ is the expression of inclusion isoforms as a fraction of total expression: Ψ = I / (I + S), where I is the inclusion read count and S the skipping (exclusion) read count. Equivalently μ̃ = γᵢ/(γᵢ+γₑ) for normalized junction expressions γᵢ, γₑ. [2]

**Percent Spliced In (length-normalized, rMATS):** ψ̂ = (I/lᵢ) / (I/lᵢ + S/lₛ), where lᵢ and lₛ are the effective lengths (number of unique isoform-specific read positions) of the inclusion and skipping isoforms. This corrects the bias toward the longer isoform when lᵢ ≠ lₛ. [3]

**Event classification (per isoform pair):** given the ordered exon sets of two isoforms, the exons unique to each isoform determine the class — an extra exon spanning an intron between two exons of the other isoform is RI; an extra cassette exon is SE; a unique exon sharing one boundary (same start/different end) is A3SS, (same end/different start) is A5SS; one unique non-overlapping exon in each isoform is MXE. [1]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ Ψ ≤ 1 when I,S ≥ 0 and I+S > 0 | Ψ is a part/whole ratio of non-negative counts [2] |
| INV-02 | I+S = 0 ⇒ Ψ = NaN (undefined) | 0/0 is undefined; rMATS/MISO add pseudo-counts to avoid it [2][3] |
| INV-03 | S=0,I>0 ⇒ Ψ=1 ; I=0,S>0 ⇒ Ψ=0 | direct from Ψ = I/(I+S) [2] |
| INV-04 | both lengths > 0 ⇒ Ψ = (I/lᵢ)/(I/lᵢ+S/lₛ) | rMATS length-normalization function [3] |
| INV-05 | a detected event references two isoforms of one gene differing in structure | an AS event is defined per isoform pair of a gene [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| inclusionReads | double | required | Reads supporting the inclusion isoform (I) | ≥ 0 |
| exclusionReads | double | required | Reads supporting the skipping isoform (S) | ≥ 0 |
| inclusionEffectiveLength | double | 0 | Effective length of inclusion isoform (lᵢ) | ≥ 0; >0 with lₛ>0 enables normalization |
| exclusionEffectiveLength | double | 0 | Effective length of skipping isoform (lₛ) | ≥ 0 |
| isoforms | IEnumerable&lt;TranscriptIsoform&gt; | required | Isoforms with gene id and ordered exon coordinates | exons ascending, one strand |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (CalculatePSI) | double | Ψ in [0,1]; NaN when no supporting reads |
| (DetectAlternativeSplicing) | IEnumerable&lt;SplicingEvent&gt; | one event per detected isoform-pair difference; `EventType` is one of the five class names |

### 3.3 Preconditions and Validation

`CalculatePSI` throws `ArgumentOutOfRangeException` if either read count is negative; returns NaN for the documented 0/0 undefined case. Length normalization is applied only when both effective lengths are strictly positive; otherwise the unnormalized ratio is used. `DetectAlternativeSplicing` treats null input as an empty sequence, skips genes with fewer than two isoforms, and emits no event for structurally identical isoform pairs. Exon coordinates are inclusive [Start, End], 5′→3′ ascending, on one strand.

## 4. Algorithm

### 4.1 High-Level Steps

1. **PSI:** validate read counts; if both effective lengths > 0 compute rate-normalized ratio (rMATS), else compute I/(I+S); return NaN if the denominator is 0.
2. **Detection:** group isoforms by gene; for genes with ≥2 isoforms, compare every isoform pair.
3. For each pair, sort exons, take exons unique to each isoform, and classify by the unique-exon pattern (MXE / A5SS / A3SS / RI / SE).

### 4.2 Decision Rules

- Two unique exons (one per isoform), non-overlapping → MutuallyExclusiveExons.
- Two unique exons sharing a start, differing end → AlternativeThreePrimeSS; sharing an end, differing start → AlternativeFivePrimeSS.
- One isoform has extra exon(s), the other none: extra exon bridging an intron between two of the other's exons → RetainedIntron; else SkippedExon.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculatePSI | O(1) | O(1) | arithmetic |
| DetectAlternativeSplicing | O(g·k²·e) | O(e) | g genes, k isoforms/gene, e exons/isoform |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [TranscriptomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs)

- `TranscriptomeAnalyzer.CalculatePSI(inclusionReads, exclusionReads, inclusionEffectiveLength, exclusionEffectiveLength)`: computes Ψ.
- `TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms)`: classifies AS events per isoform pair.

### 5.2 Current Behavior

Classification compares exon coordinate sets rather than searching sequence text, so the repository suffix tree is **not** applicable (no substring/occurrence search is involved). PSI defaults to the unnormalized read-count ratio (the definition shared by [1][2] and SUPPA2 [5]); supplying both effective lengths switches to the rMATS length-normalized form [3]. Detected events carry `InclusionLevel = NaN` because the inclusion level requires read counts, obtained separately via `CalculatePSI`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Ψ = I/(I+S) [2]; rMATS ψ̂ = (I/lᵢ)/(I/lᵢ+S/lₛ) [3].
- Five-class event taxonomy SE/RI/A5SS/A3SS/MXE [1][4].

**Intentionally simplified:**

- Event detection is deterministic structural comparison of annotated exon coordinates; it does not run the rMATS/MISO probabilistic isoform-abundance model. **Consequence:** no statistical confidence or differential ΔPSI testing is produced here (see `DetectDifferentialSplicing` for ΔPSI between conditions).

**Not implemented:**

- Read-level isoform assignment from alignments; **users should rely on:** an upstream aligner/quantifier (rMATS, SUPPA2) to produce the read counts fed to `CalculatePSI`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| I+S = 0 | Ψ = NaN | 0/0 undefined [2] |
| S = 0, I > 0 | Ψ = 1 | full inclusion [2] |
| negative reads | ArgumentOutOfRangeException | counts are non-negative |
| < 2 isoforms for a gene | no event | event needs two isoforms [1] |
| identical isoforms | no event | no structural difference |
| null isoforms | empty result | tolerant input contract |

### 6.2 Limitations

Classification assumes correctly ordered exon annotations on a single strand and one structural difference per isoform pair for unambiguous single-class assignment; complex multi-difference pairs fall back to SkippedExon. It does not model strand inference, intron-coordinate conventions across annotation formats, or read-level quantification uncertainty.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// PSI from read counts: 80 inclusion, 20 skipping → 0.80
double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20);

// rMATS length-normalized: I=80,S=20,l_I=200,l_S=100 → 0.6667
double psiN = TranscriptomeAnalyzer.CalculatePSI(80, 20, 200, 100);

// Classify a skipped-exon isoform pair
var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs) — covers INV-01..INV-05
- Evidence: [TRANS-SPLICE-001-Evidence.md](../../../docs/Evidence/TRANS-SPLICE-001-Evidence.md)
- Related algorithms: [Expression_Quantification](Expression_Quantification.md)

## 8. References

1. Wang ET, Sandberg R, Luo S, Khrebtukova I, Zhang L, Mayr C, Kingsmore SF, Schroth GP, Burge CB. 2008. Alternative isoform regulation in human tissue transcriptomes. Nature 456(7221):470–476. https://doi.org/10.1038/nature07509
2. Challenges in estimating percent inclusion of alternatively spliced junctions from RNA-seq data. 2012. BMC Bioinformatics 13(Suppl 6):S11. https://pmc.ncbi.nlm.nih.gov/articles/PMC3330053/
3. Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. 2014. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. PNAS 111(51):E5593–E5601. https://doi.org/10.1073/pnas.1419161111
4. rMATS project documentation. https://rmats.sourceforge.io/
5. Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E. 2018. SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. https://pubmed.ncbi.nlm.nih.gov/29571299/
