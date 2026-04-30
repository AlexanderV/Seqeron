# Gene Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-GENE-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Gene prediction in this repository is a prokaryote-oriented, ORF-first heuristic. `GenomeAnnotator.PredictGenes(...)` converts qualifying ORFs into `CDS` annotations, while `FindRibosomeBindingSites(...)` separately scans upstream regions for Shine-Dalgarno-like motifs. The implementation is deterministic for its fixed codon and motif sets, but it does not implement a trained promoter or coding-potential model and it does not resolve overlaps between competing ORFs. It is therefore best understood as annotation scaffolding rather than a substitute for specialized gene finders such as GLIMMER or GeneMark.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The current document models prokaryotic genes as continuous coding regions without introns. In that model, an annotated coding sequence is associated with a start codon, a stop codon, and often upstream regulatory signals such as promoter elements and a ribosome-binding site. The current document identifies the canonical bacterial promoter elements as the `-35` box (`TTGACA`) and the `-10` box (`TATAAT`), and it describes the Shine-Dalgarno (SD) ribosome-binding sequence as an upstream translation-initiation signal that base-pairs with the `3'` end of `16S` rRNA.

The same source material describes the bacterial SD consensus as `AGGAGG` (or `AGGAGGU` in the E. coli form cited in the document), with shorter variants such as `GGAGG`, `AGGAG`, `GAGG`, and `AGGA`. The cited range for SD placement is `4` to `15` nucleotides upstream of the start codon, with Chen et al. (1994) reporting an optimal aligned spacing of `5` nucleotides from the motif `3'` end to the start codon.

### 2.2 Core Model

The gene-finding model documented here combines two classical signals:

- A coding candidate is an ORF bounded by an accepted start codon and an in-frame stop codon.
- A potential translation-initiation signal is a Shine-Dalgarno-like motif located a short distance upstream of the start codon.

For a start codon at index $s$ and a stop codon beginning at index $t$ in the same frame, the candidate coding interval spans $[s, t + 3)$. The upstream SD distance is measured from the motif `3'` end to the first nucleotide of the start codon, so only motifs with aligned spacing in the supported interval are considered valid initiation-site candidates.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The target gene structure is prokaryotic and can be approximated as a continuous ORF without introns | Eukaryotic genes with introns or alternative splicing are not represented correctly |
| ASM-02 | Translation initiation can be approximated by short Shine-Dalgarno-like sequence motifs upstream of the start codon | Genes that rely on different initiation signals may not receive supporting RBS hits |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every predicted gene is derived from an ORF that begins with `ATG`, `GTG`, or `TTG` and ends with `TAA`, `TAG`, or `TGA` | `PredictGenes(...)` delegates to the canonical ORF finder with `requireStartCodon: true` |
| INV-02 | Every predicted gene has valid genomic bounds and a frame attribute in `{1, 2, 3}` | `PredictGenes(...)` emits `GeneAnnotation` records from `OpenReadingFrame` values and preserves the positive frame index |
| INV-03 | Every reported RBS hit lies between `minDistance` and `maxDistance` from the associated start codon | `FindRibosomeBindingSites(...)` checks aligned spacing before emitting a hit |
| INV-04 | The current RBS scores are bounded between `4/6` and `1.0` for the supported motif library | The implementation normalizes motif length against the 6-base consensus length |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dnaSequence` | `string` | required | DNA sequence analyzed by both `PredictGenes(...)` and `FindRibosomeBindingSites(...)` | Case-insensitive; null or empty input yields no results in these two helpers |
| `minOrfLength` | `int` | `100` | Minimum ORF length in amino acids for `PredictGenes(...)` | Applied to ORFs before annotation records are created |
| `prefix` | `string` | `"gene"` | Prefix used when generating sequential gene IDs | The implementation appends `_{number:D4}` |
| `upstreamWindow` | `int` | `20` | Upstream search window for `FindRibosomeBindingSites(...)` | Used to bound the motif scan before each qualifying ORF start |
| `minDistance` | `int` | `4` | Minimum aligned distance from SD motif `3'` end to start codon | Must be less than or equal to `maxDistance` for meaningful scans |
| `maxDistance` | `int` | `15` | Maximum aligned distance from SD motif `3'` end to start codon | Values beyond this range are filtered out |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `PredictGenes(...)` | `IEnumerable<GeneAnnotation>` | Sequence of gene annotations derived from qualifying ORFs |
| `GeneId` | `string` | Sequential identifier such as `gene_0001` or `test_0001` |
| `Start` | `int` | 0-based inclusive start coordinate of the ORF-derived gene annotation |
| `End` | `int` | 0-based exclusive end coordinate of the ORF-derived gene annotation |
| `Strand` | `char` | `+` for forward-strand genes and `-` for reverse-strand genes |
| `Type` | `string` | Always `CDS` in the current implementation |
| `Product` | `string` | Always `hypothetical protein` in the current implementation |
| `Attributes` | `IReadOnlyDictionary<string, string>` | Includes `frame`, `protein_length`, and `translation` |
| `FindRibosomeBindingSites(...)` | `IEnumerable<(int position, string sequence, double score)>` | Sequence of upstream SD-like motif hits |
| `position` | `int` | 0-based genomic start of the detected motif |
| `sequence` | `string` | Matched SD-like motif string |
| `score` | `double` | Length-normalized motif score in the current implementation |

### 3.3 Preconditions and Validation

Both `PredictGenes(...)` and `FindRibosomeBindingSites(...)` return an empty sequence for null or empty input. Codon and motif matching are case-insensitive because the implementation uppercases scanned substrings internally. Neither helper performs full alphabet validation, so non-`A/C/G/T` characters are only relevant when they prevent exact start, stop, or motif matches. `PredictGenes(...)` returns ORF-derived coordinates in the repository's internal 0-based, end-exclusive convention. `FindRibosomeBindingSites(...)` first discovers ORFs and then scans only forward-strand ORFs for upstream motifs, so sequences without a qualifying ORF yield no RBS hits even if an SD-like motif is present.

## 4. Algorithm

### 4.1 High-Level Steps

1. Find ORFs across both strands using the repository's canonical ORF finder.
2. For `PredictGenes(...)`, keep ORFs whose translated length is at least `minOrfLength`, order them by genomic start coordinate, and emit `GeneAnnotation` records with sequential IDs and fixed metadata.
3. For `FindRibosomeBindingSites(...)`, find ORFs with a minimum length of `30` amino acids, inspect the upstream region before each forward-strand ORF, and scan that region for exact SD-like motifs.
4. Emit only those motif hits whose aligned distance from the motif `3'` end to the ORF start lies within `[minDistance, maxDistance]`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Category | Values / Rule |
|----------|---------------|
| Start codons | `ATG`, `GTG`, `TTG` |
| Stop codons | `TAA`, `TAG`, `TGA` |
| SD-like motif library | `AGGAGG`, `GGAGG`, `AGGAG`, `GAGG`, `AGGA` |
| RBS scoring | `score = motif.Length / 6.0` |
| SD distance measurement | Aligned spacing from motif `3'` end to the first nucleotide of the start codon |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictGenes(...)` | `O(n + m log m)` | `O(m)` | `n` = sequence length, `m` = ORFs materialized and then sorted by start coordinate |
| `FindRibosomeBindingSites(...)` | `O(n + m × upstreamWindow)` | `O(m)` | Includes ORF discovery plus motif scans over a fixed motif library |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.PredictGenes(string, int, string)`: Converts canonical ORFs into `GeneAnnotation` records.
- `GenomeAnnotator.FindRibosomeBindingSites(string, int, int, int)`: Scans upstream regions of qualifying ORFs for SD-like motifs.
- `GenomeAnnotator.FindOrfs(string, int, bool, bool)`: Underlying ORF finder used by both helpers.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `PredictGenes(...)` calls `FindOrfs(...)` with `searchBothStrands: true` and `requireStartCodon: true`, then sorts the resulting ORFs by `Start` before generating gene IDs.
- Every qualifying ORF is emitted as its own `CDS` annotation, including overlapping or nested candidates; there is no best-model selection or overlap suppression step.
- Every emitted annotation has `Type = "CDS"` and `Product = "hypothetical protein"`.
- The emitted attributes are `frame`, `protein_length`, and `translation`; `protein_length` trims the terminal `*` from the translated protein, while `translation` preserves the raw translated sequence.
- `FindRibosomeBindingSites(...)` internally uses `FindOrfs(dnaSequence, minLength: 30)` and then filters to forward-strand ORFs before scanning upstream windows.
- The current RBS score is based only on motif length, not on the literature-derived spacing optimum or a full affinity model.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ORF-based identification of continuous coding candidates.
- Recognition of the start codons `ATG`, `GTG`, and `TTG` and the stop codons `TAA`, `TAG`, and `TGA`.
- Shine-Dalgarno-like motif scanning in an upstream distance-constrained window.

**Intentionally simplified:**

- `PredictGenes(...)` uses ORF structure only and does not incorporate promoter boxes or RBS scores into the generated gene annotations; **consequence:** the gene list is not ranked or filtered by upstream regulatory evidence.
- `PredictGenes(...)` emits every qualifying ORF after start-coordinate sorting; **consequence:** overlapping or nested coding candidates remain separate predictions instead of being reconciled into one preferred model.
- `FindRibosomeBindingSites(...)` uses a fixed motif library and length-normalized scores; **consequence:** the reported scores do not model the literature's detailed initiation-strength differences.
- The helper scans only forward-strand ORFs for RBS hits; **consequence:** reverse-strand genes can be predicted without a corresponding RBS record from this method.

**Not implemented:**

- Promoter `-10` / `-35` scoring integrated into gene prediction; **users should rely on:** [Promoter_Detection.md](Promoter_Detection.md) for the repository's separate promoter-motif helper.
- Codon-bias, organism-specific training, or comparative gene-finding models; **users should rely on:** no current alternative in this repository.
- Eukaryotic intron-aware gene prediction; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | RBS detection is separate from gene prediction | Deviation | Users must combine `PredictGenes(...)` and `FindRibosomeBindingSites(...)` results themselves if they want joint interpretation | accepted | `PredictGenes(...)` does not call the RBS helper |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty sequence | Returns no genes and no RBS hits | Both helpers short-circuit through the underlying ORF scan |
| No valid ORF | Returns no genes and no RBS hits | RBS scanning is ORF-driven |
| Lowercase or mixed-case input | Handled identically to uppercase DNA | Codon and motif scans uppercase substrings internally |
| Reverse-strand ORF | Gene coordinates are remapped to the original forward coordinate system and `Strand = '-'` | `PredictGenes(...)` uses remapped `OpenReadingFrame` coordinates |
| ORFs below `minOrfLength` | Excluded from `PredictGenes(...)` | Length filtering happens before annotation generation |
| Full `AGGAGG` motif outside `[minDistance, maxDistance]` | Filtered out even if a shorter motif at the same locus still qualifies | Aligned spacing is measured after subtracting the specific motif length |

### 6.2 Limitations

The repository implements a simple ORF-based predictor rather than a trained gene-finding pipeline. It does not score promoters, does not model codon bias or organism-specific translation initiation, does not resolve overlapping candidates, and does not support eukaryotic introns or splice-aware gene models. The RBS helper is likewise a motif scanner, not a full translation-initiation model.

## 8. References

1. Wikipedia contributors. Gene prediction. https://en.wikipedia.org/wiki/Gene_prediction
2. Wikipedia contributors. Shine-Dalgarno sequence. https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence
3. Wikipedia contributors. Ribosome-binding site. https://en.wikipedia.org/wiki/Ribosome-binding_site
4. Shine J, Dalgarno L. Determinant of cistron specificity in bacterial ribosomes. Nature. 1975;254:34-38.
5. Chen H, Bjerknes M, Kumar R, Jay E. Determination of the optimal aligned spacing between the Shine-Dalgarno sequence and the translation initiation codon. Nucleic Acids Research. 1994;22(23):4953-4957.
6. Laursen BS, Sorensen HP, Mortensen KK, Sperling-Petersen HU. Initiation of protein synthesis in bacteria. Microbiology and Molecular Biology Reviews. 2005;69(1):101-123.
7. Stormo GD, Schneider TD, Gold L, Ehrenfeucht A. Characterization of translational initiation sites in E. coli. Nucleic Acids Research. 1982;10(9):2971-2996.
