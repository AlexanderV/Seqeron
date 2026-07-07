# Fusion Breakpoint Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-FUSION-003 |
| Related Projects | Seqeron.Genomics.Oncology, Seqeron.Genomics.Core |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Fusion Breakpoint Analysis classifies a gene-fusion junction and predicts the chimeric protein product. Given the breakpoint site categories on each partner and the coding-frame quantities at the junction, `AnalyzeBreakpoint` reports whether the junction joins two coding regions and whether it preserves the reading frame; `PredictFusionProtein` assembles the chimeric coding sequence (5' partner CDS up to the breakpoint, joined to the 3' partner CDS from the breakpoint), translates it with the standard genetic code, and truncates the peptide at the first stop codon. The behavior is specification-driven (Arriba output schema) and reference-implementation-driven (AGFusion); it is deterministic. Partner CDS sequences and breakpoint offsets are caller-supplied (the library bundles no genome/transcript annotation), so the unit is a Framework algorithm.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A gene fusion joins parts of two genes at a breakpoint. The clinical consequence depends on (a) where each breakpoint falls in its partner — coding sequence (CDS), UTR, intron, exon, or intergenic — and (b) whether the junction keeps the 3' partner's codons in phase with the 5' partner's reading frame. Arriba reports the breakpoint location per partner in its `site1`/`site2` columns and the junction frame in its `reading_frame` column [1]. AGFusion predicts the cDNA/CDS/protein product for each isoform combination and labels the functional effect (in-frame / out-of-frame) [2].

### 2.2 Core Model

**Site categories.** Each breakpoint is one of `5'UTR`, `3'UTR`, `UTR`, `CDS`, `exon`, `intron`, `intergenic` [1]. A reading frame can only be joined when both breakpoints fall in CDS.

**Reading-frame rule.** Reading frames are read as triplets of nucleotides [3]. The 3' partner stays in phase iff the coding bases the 5' partner contributes up to the breakpoint, minus the 3' partner's coding-start phase, is a multiple of three:

```
in-frame  ⟺  (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0
```

This is the codon-phase form of AGFusion's check that the chimeric CDS keeps the junction at a codon boundary [2].

**Chimeric CDS and protein (AGFusion).** With the 5' CDS prefix to the breakpoint and the 3' CDS suffix from the breakpoint [2]:

```
cds_5prime = transcript1.coding_sequence[0 : junction5]
cds_3prime = transcript2.coding_sequence[junction3 : ]
chimeric    = cds_5prime + cds_3prime
protein     = translate(chimeric); protein = protein[0 : protein.find("*")]
```

For an out-of-frame junction AGFusion trims the chimeric CDS to a whole number of codons before translating (`chimeric[0 : 3*(len//3)]`), so the 3' partner is read in its shifted frame [2]. Translation uses the standard genetic code, NCBI Table 1 [3] (UAA/UAG/UGA = stop).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A reading-frame call (InFrame/OutOfFrame) is made only when both breakpoints are CDS; else NotPredicted | Arriba `reading_frame = .` when peptide is not predictable [1] |
| INV-02 | InFrame ⟺ `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0` | Triplet reading frame [3]; AGFusion frame rule [2] |
| INV-03 | The predicted peptide contains no internal stop; translation ends at the first stop codon | AGFusion `protein[0:find("*")]` [2] |
| INV-04 | ChimericCds = 5' CDS prefix (length `junction5`) ++ 3' CDS suffix (from `junction3`) | AGFusion concatenation [2] |
| INV-05 | An out-of-frame chimeric CDS is trimmed to whole codons before translation | AGFusion out-of-frame branch [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| fusion | FusionBreakpoint | required | Partners, per-partner site categories, and codon-frame quantities | FivePrimeCodingBases ≥ 0; ThreePrimeStartPhase ∈ {0,1,2} for a frame call |
| transcripts | (string FivePrimeCds, string ThreePrimeCds) | required | Full partner CDS sequences (DNA, A/C/G/T) | non-null; offsets in range of their CDS |

The breakpoint offsets are taken from `fusion`: the 5' prefix length is `FivePrimeCodingBases`, the 3' suffix start is `ThreePrimeStartPhase`.

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| BreakpointAnalysis.BreakpointInCoding | bool | True iff both breakpoints are CDS |
| BreakpointAnalysis.FrameStatus | BreakpointFrameStatus | InFrame / OutOfFrame / StopCodon / NotPredicted |
| FusionProteinPrediction.ChimericCds | string | 5' CDS prefix ++ 3' CDS suffix (uppercased DNA) |
| FusionProteinPrediction.Peptide | string | Translated peptide, truncated at first stop codon |
| FusionProteinPrediction.Effect | BreakpointFrameStatus | InFrame / OutOfFrame |
| FusionProteinPrediction.HasPrematureStop | bool | True iff a stop codon was reached before the ORF end |

### 3.3 Preconditions and Validation

Inputs are case-normalized to uppercase DNA. CDS strings must be non-null (`ArgumentNullException`). Offsets must be within their CDS and the 3' start phase within {0,1,2} for the protein prediction (`ArgumentOutOfRangeException`); `AnalyzeBreakpoint` validates the frame quantities only when both breakpoints are CDS (delegating to `IsInFrame`). Translation is 0-based, frame 0, codons read 5'→3' in triplets; a trailing partial codon is not translated.

## 4. Algorithm

### 4.1 High-Level Steps

1. `AnalyzeBreakpoint`: if both sites are CDS, call the reading frame via the codon-phase rule; otherwise report NotPredicted.
2. `PredictFusionProtein`: take the 5' CDS prefix `[0:junction5]` and the 3' CDS suffix `[junction3:]`; concatenate.
3. Determine the frame effect (in-frame / out-of-frame) by the codon-phase rule.
4. Translate the chimeric CDS codon-by-codon with the standard genetic code, stopping at the first stop codon; report the peptide and a premature-stop flag.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Site categories follow Arriba `site1`/`site2` values [1].
- The genetic code is the shared `Seqeron.Genomics.Core.GeneticCode.Standard` (NCBI Table 1); the codon table is not duplicated here.
- `CodonLength = 3` (triplet reading frame) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| AnalyzeBreakpoint | O(1) | O(1) | constant-time site/frame checks |
| PredictFusionProtein | O(n) | O(n) | n = chimeric CDS length; one pass over codons |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.AnalyzeBreakpoint(FusionBreakpoint)`: site categories + junction reading-frame consequence.
- `OncologyAnalyzer.PredictFusionProtein(FusionBreakpoint, (string, string))`: chimeric CDS, translation, first-stop truncation.
- Reuses `OncologyAnalyzer.IsInFrame(int, int)` (ONCO-FUSION-001) and `Seqeron.Genomics.Core.GeneticCode.Standard.Translate`.

### 5.2 Current Behavior

The repository has no genome/GTF/transcript database, so partner CDS sequences and breakpoint offsets are passed in directly rather than resolved from an annotation as AGFusion does. The concatenation, frame rule, translation, and first-stop truncation match AGFusion exactly; only the source of the inputs differs (an API-shape decision, no output change). No substring search is performed (the chimeric CDS is built by slice/concat), so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Site categories from the Arriba `site1`/`site2` vocabulary [1].
- Reading-frame call only for coding-to-coding junctions; otherwise NotPredicted (Arriba `reading_frame = .`) [1].
- Codon-phase in-frame rule `(b − p) mod 3 == 0` [2][3].
- Chimeric CDS = 5' prefix ++ 3' suffix; translate; truncate at first stop; out-of-frame trims to whole codons [2].

**Intentionally simplified:**

- Frame effect labels: AGFusion also emits "in-frame (with mutation)" for mid-codon junctions whose fractional codon parts complement; **consequence:** such junctions are reported here as InFrame (the chimeric CDS is still a multiple of 3 and translates in-frame), without the separate "with mutation" sub-label.

**Not implemented:**

- Transcript/isoform resolution from a genome annotation, Pfam domain architecture, and the assembled `fusion_transcript`/`peptide_sequence` glyphs (`|`, lowercase, `___`); **users should rely on:** AGFusion / Arriba for full isoform-aware visualization.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Breakpoint not in CDS (UTR/intron/intergenic) | FrameStatus = NotPredicted | Arriba `reading_frame = .` [1] |
| Premature stop in chimeric ORF | Peptide truncated, HasPrematureStop = true | AGFusion `protein[0:find("*")]` [2] |
| Out-of-frame junction | CDS trimmed to whole codons; 3' read shifted | AGFusion out-of-frame branch [2] |
| Empty 3' CDS suffix | Peptide from 5' prefix only | Concatenation with empty suffix [2] |
| ThreePrimeStartPhase ∉ {0,1,2} (frame) | OutOfFrame for protein; ArgumentOutOfRangeException via IsInFrame in AnalyzeBreakpoint | Phase is a codon offset (0,1,2) |

### 6.2 Limitations

No isoform/annotation resolution, no protein-domain annotation, no read-level transcript assembly, and no handling of non-template insertions at the junction. Frame is determined solely from the supplied coding-base count and phase; if those are wrong the call is wrong (the documented pitfall is using gene-level rather than CDS-level lengths).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var bp = new OncologyAnalyzer.FusionBreakpoint(
    "EML4", "ALK",
    OncologyAnalyzer.BreakpointSite.Cds, OncologyAnalyzer.BreakpointSite.Cds,
    FivePrimeCodingBases: 6, ThreePrimeStartPhase: 0);

var analysis = OncologyAnalyzer.AnalyzeBreakpoint(bp);          // FrameStatus = InFrame
var protein  = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));
// protein.ChimericCds = "ATGAAAGATGGT"; protein.Peptide = "MKDG"; Effect = InFrame
```

**Numerical walk-through:** chimeric `ATGAAAGATGGT` → codons ATG(M) AAA(K) GAT(D) GGT(G) → `MKDG`. With 3' CDS `GATTAAGGT`: chimeric `ATGAAAGATTAAGGT` → ATG(M) AAA(K) GAT(D) TAA(stop) → peptide `MKD`, HasPrematureStop = true.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-FUSION-003-Evidence.md](../../../docs/Evidence/ONCO-FUSION-003-Evidence.md)
- Related algorithms: [Fusion_Gene_Detection](Fusion_Gene_Detection.md), [Known_Fusion_Database_Lookup](Known_Fusion_Database_Lookup.md)

## 8. References

1. Uhrig S, Ellermann J, Walther T, et al. 2021. Accurate and efficient detection of gene fusions from RNA sequencing data. *Genome Research* 31(3):448–460. https://doi.org/10.1101/gr.257246.119 (output schema: https://github.com/suhrig/arriba/wiki/05-Output-files).
2. Murphy C, Elemento O. 2016. AGFusion: annotate and visualize gene fusions. *bioRxiv* 080903. https://doi.org/10.1101/080903 (source: https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py).
3. Badger JH, Olsen GJ. 1999. CRITICA: Coding Region Identification Tool Invoking Comparative Analysis. *Mol Biol Evol* 16(4):512–524. https://doi.org/10.1093/oxfordjournals.molbev.a026133 (via https://en.wikipedia.org/wiki/Reading_frame).
</content>
