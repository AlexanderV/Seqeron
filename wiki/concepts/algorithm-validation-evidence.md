---
type: concept
title: "Algorithm validation evidence artifacts"
tags: [validation, testing]
sources:
  - docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
  - docs/Evidence/ALIGN-MULTI-001-Evidence.md
  - docs/Evidence/ALIGN-SEMI-001-Evidence.md
  - docs/Evidence/ALIGN-STATS-001-Evidence.md
  - docs/Evidence/ANNOT-CODING-001-Evidence.md
  - docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
  - docs/Evidence/ANNOT-REPEAT-001-Evidence.md
  - docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
  - docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
  - docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
  - docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
  - docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
  - docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
  - docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
  - docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
  - docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
  - docs/Evidence/CHROM-ANEU-001-Evidence.md
  - docs/Evidence/CHROM-CENT-001-Evidence.md
  - docs/Evidence/CHROM-KARYO-001-Evidence.md
  - docs/Evidence/CHROM-SYNT-001-Evidence.md
  - docs/Evidence/CHROM-TELO-001-Evidence.md
  - docs/Evidence/CODON-CAI-001-Evidence.md
  - docs/Evidence/CODON-ENC-001-Evidence.md
  - docs/Evidence/CODON-OPT-001-Evidence.md
  - docs/Evidence/CODON-RARE-001-Evidence.md
  - docs/Evidence/CODON-RSCU-001-Evidence.md
  - docs/Evidence/CODON-STATS-001-Evidence.md
  - docs/Evidence/CODON-USAGE-001-Evidence.md
  - docs/Evidence/COMPGEN-ANI-001-Evidence.md
  - docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
  - docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
  - docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
  - docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
  - docs/Evidence/COMPGEN-RBH-001-Evidence.md
  - docs/Evidence/COMPGEN-REARR-001-Evidence.md
  - docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
  - docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md
  - docs/Evidence/DISORDER-LC-001-Evidence.md
  - docs/Evidence/DISORDER-MORF-001-Evidence.md
  - docs/Evidence/DISORDER-PRED-001-Evidence.md
source_commit: 05fff695e889b79023301d7319afbc8a24e0bec4
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-global-001-evidence
      evidence: "Evidence Artifact: ALIGN-GLOBAL-001 ... Online Sources ... Test Dataset ... Deviations and Assumptions ... References"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:definition-of-done
      source: align-global-001-evidence
      evidence: "Deviations and Assumptions: None ... the implementation follows the standard Needleman–Wunsch linear gap penalty model exactly as given in the Wikipedia pseudocode"
      confidence: medium
      status: current
---

# Algorithm validation evidence artifacts

Each algorithm [[test-unit-registry|test unit]] has a per-unit **Evidence artifact** under
`docs/Evidence/<UnitID>-Evidence.md`. These are the literature-traced source record behind
the [[definition-of-done]]'s "Evidence documented" criterion and the
[[validation-and-testing]] campaign: they pin exactly which external references and worked
examples the implementation and its tests are validated against.

## Templated structure

Every evidence file follows the same shape:

1. **Header** — Test Unit ID, algorithm name, date collected.
2. **Online sources** — Wikipedia / primary-literature URLs with access dates and the key
   extracted points (definitions, recurrences, complexity, worked examples).
3. **Test dataset** — the canonical worked example(s) with exact parameters and expected
   outputs, used as the oracle / differential test fixture.
4. **Deviations and assumptions** — where the implementation departs from (or exactly
   follows) the reference, plus API-contract behaviours (null/empty handling) that sit
   outside the algorithm spec.
5. **References** — primary literature and encyclopedic citations.

Because these files are near-templated across the ~213 documented units, the wiki keeps
**one** shared page for the pattern (this page) plus a concise per-file source summary
(e.g. [[align-global-001-evidence]], [[align-multi-001-evidence]],
[[align-semi-001-evidence]], [[align-stats-001-evidence]], [[annot-coding-001-evidence]],
[[annot-codonusage-001-evidence]], [[annot-repeat-001-evidence]],
[[assembly-consensus-001-evidence]], [[assembly-correct-001-evidence]],
[[assembly-cover-001-evidence]], [[assembly-dbg-001-evidence]],
[[assembly-merge-001-evidence]], [[assembly-olc-001-evidence]],
[[assembly-scaffold-001-evidence]], [[assembly-stats-001-evidence]],
[[assembly-trim-001-evidence]], [[chrom-aneu-001-evidence]],
[[chrom-cent-001-evidence]], [[chrom-karyo-001-evidence]], [[chrom-synt-001-evidence]],
[[chrom-telo-001-evidence]], [[codon-cai-001-evidence]],
[[codon-enc-001-evidence]], [[codon-opt-001-evidence]], [[codon-rare-001-evidence]],
[[codon-rscu-001-evidence]], [[codon-stats-001-evidence]],
[[codon-usage-001-evidence]], [[compgen-ani-001-evidence]],
[[compgen-cluster-001-evidence]], [[compgen-compare-001-evidence]],
[[compgen-dotplot-001-evidence]], [[compgen-ortho-001-evidence]],
[[compgen-rbh-001-evidence]], [[compgen-rearr-001-evidence]],
[[compgen-reversal-001-evidence]], [[compgen-synteny-001-evidence]],
[[disorder-lc-001-evidence]], [[disorder-morf-001-evidence]],
[[disorder-pred-001-evidence]]). An
individual algorithm gets its own concept page only when it is itself distinct and wiki-worthy
— for example [[global-alignment-needleman-wunsch]], [[multiple-sequence-alignment]],
[[semi-global-alignment-fitting]], [[alignment-statistics]],
[[coding-potential-hexamer-score]], [[relative-synonymous-codon-usage]],
[[codon-adaptation-index]] (the CAI index in the codon-usage family),
[[effective-number-of-codons]] (the reference-free ENC/Nc measure in the codon-usage family),
[[codon-optimization]] (the CDS-rewriting operation of the codon-usage family),
[[rare-codon-analysis]] (the thresholded-frequency + %MinMax/Sherlocc cluster-detection unit of the codon-usage family),
[[codon-usage-comparison]] (the raw codon-count table + TVD distribution-comparison unit of the codon-usage family),
[[repetitive-element-detection]] (the anchor for the repeats/tandem family),
[[consensus-sequence]] (the anchor for the assembly CONSENSUS family),
[[kmer-spectrum-error-correction]] (the anchor for the assembly CORRECT family),
[[coverage-depth-calculation]] (the anchor for the assembly COVER family),
[[de-bruijn-graph-assembly]] (the anchor for the assembly DBG family),
[[contig-merge-overlap-collapse]] (the anchor for the assembly MERGE family),
[[overlap-layout-consensus-assembly]] (the anchor for the assembly OLC family), or
[[scaffolding]] (the anchor for the assembly SCAFFOLD family),
[[assembly-statistics]] (the anchor for the assembly STATS family), or
[[quality-trimming-running-sum]] (the anchor for the assembly TRIM family), or
[[aneuploidy-detection]] (the anchor for the chromosome-analysis copy-number/ploidy family), or
[[centromere-analysis]] (the anchor for the chromosome centromere / alpha-satellite family), or
[[karyotype-analysis]] (the anchor for the chromosome karyotyping / ploidy-detection family), or
[[synteny-and-rearrangement-detection]] (the shared anchor for the chromosome + comparative-genomics synteny/rearrangement family), or
[[telomere-analysis]] (the anchor for the chromosome telomere family), or
[[average-nucleotide-identity]] (the anchor for the comparative-genomics ANI genome-similarity family), or
[[conserved-gene-clusters-common-intervals]] (the comparative-genomics common-interval / conserved-cluster unit), or
[[genome-comparison-core-dispensable]] (the comparative-genomics end-to-end genome-comparison pipeline — core/dispensable partition + syntenic fraction), or
[[dot-plot-word-match]] (the comparative-genomics word-match / k-tuple dot-matrix visual sequence-comparison unit), or
[[ortholog-detection-reciprocal-best-hits]] (the comparative-genomics RBH ortholog + within-genome in-paralog detection unit, shared anchor for COMPGEN-RBH-001), or
[[genome-rearrangement-breakpoint-distance]] (the comparative-genomics signed-permutation / breakpoint-distance rearrangement-detection unit, the alternative formulation to the block-signal [[synteny-and-rearrangement-detection]]), or
[[protein-low-complexity-seg]] (the anchor for the protein disorder / features family — SEG low-complexity region detection, the first ingested DISORDER-* unit), or
[[morf-prediction-dip-in-disorder]] (the MoRF dip-in-disorder prediction unit of the protein-disorder family, the second ingested DISORDER-* unit), or
[[intrinsic-disorder-prediction-top-idp]] (the TOP-IDP `PredictDisorder` sliding-window intrinsic-disorder anchor of the protein-disorder family, the third ingested DISORDER-* unit that MoRF + region detection sit on).
