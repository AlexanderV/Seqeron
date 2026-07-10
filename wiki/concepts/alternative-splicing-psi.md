---
type: concept
title: "Alternative splicing: event classification + Percent-Spliced-In (PSI / Ψ, ΔPSI)"
tags: [transcriptome, splicing, algorithm]
mcp_tools:
  - detect_alternative_splicing
sources:
  - docs/Evidence/TRANS-SPLICE-001-Evidence.md
  - docs/Validation/reports/TRANS-SPLICE-001.md
source_commit: 82e3e03992f6e370559efdde3124a4b870a57893
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-splice-001-evidence
      evidence: "Test Unit ID: TRANS-SPLICE-001, Algorithm: Alternative Splicing — Event Classification and Percent-Spliced-In (PSI); Methods TranscriptomeAnalyzer.CalculatePSI / DetectAlternativeSplicing"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:differential-expression
      source: trans-splice-001-evidence
      evidence: "Both TRANS-SPLICE-001 and TRANS-DIFF-001 are TranscriptomeAnalyzer units of the Transcriptome/RNA-seq family; PSI/ΔPSI answers 'which exon is spliced in' as the splicing-level counterpart to the anchor's gene-level differential expression."
      confidence: high
      status: current
---

# Alternative splicing: event classification + Percent-Spliced-In (PSI / Ψ)

**Quantify how a gene's exons are spliced together from RNA-seq reads.** Where
[[differential-expression]] asks whether a *gene's total level* changes, this unit asks *which
isoform / which exon* is used — it (1) **classifies** an alternative-splicing (AS) event between two
isoforms of one gene into a canonical class, and (2) estimates **Percent-Spliced-In (PSI, Ψ)** — the
fraction of transcripts that *include* an alternatively spliced segment — from inclusion vs exclusion
read support. A **Transcriptome / RNA-seq family** sibling of the anchor [[differential-expression]],
distinct in method (read-ratio splicing quantification, not a two-group mean test). Validated under
test unit **TRANS-SPLICE-001**; the record is [[trans-splice-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern. The
two-stage validation verdict is recorded in [[trans-splice-001-report]] — **Stage A
PASS-WITH-NOTES · Stage B FAIL → fixed, State ✅ CLEAN** (the swapped A5SS/A3SS event labels
in `ClassifyIsoformPair` were corrected in-session; full suite 6501/0).

Impl `TranscriptomeAnalyzer` (`Seqeron.Genomics.Annotation`, the Annotation server):
`CalculatePSI(...)` → the inclusion fraction, and `DetectAlternativeSplicing(...)` → the AS event
class for a pair of isoforms.

> **Not the genomic splice-site predictors.** This is **read-quantification of splicing usage** (PSI
> from RNA-seq reads). It is a different problem from [[splice-donor-site-prediction]] /
> [[splice-acceptor-site-prediction]] / [[gene-structure-prediction-intron-exon]], which score the
> `GT`/`AG` **sequence motif** at a candidate splice site. Same biology, orthogonal task: those
> predict *where a splice site is*; this quantifies *how often an exon is spliced in*.

## 1. Percent-Spliced-In (PSI / Ψ)

PSI is the **inclusion fraction** — inclusion read support over the total inclusion + exclusion
support (Wang et al. 2008; BMC Bioinformatics S11 / PMC3330053; SUPPA2):

```
Ψ = I / (I + S)            I = inclusion reads,  S = skipping (exclusion) reads
```

**Read classification** (PMC3330053, SUPPA): reads on the alternative exon *or* on its junctions with
the adjacent constitutive exons support **inclusion**; reads on the junction *between* the two
adjacent constitutive exons (which skips the alternative exon) support **exclusion**.

### Length-normalized PSI (rMATS, opt-in)

Reads accumulate in proportion to isoform length, so rMATS (Shen et al. 2014) divides each count by
the isoform's **effective length** (the number of unique isoform-specific read positions) before
forming the ratio:

```
ψ̂ = (I / lᵢ) / (I / lᵢ + S / lₛ)          lᵢ, lₛ = inclusion / skipping effective lengths
```

with the binomial model `I | ψ ∼ Binomial(n = I+S, p = lᵢψ / (lᵢψ + lₛ(1−ψ)))`. `CalculatePSI`
**defaults to the unnormalized ratio** `I/(I+S)` (the definition shared by Wang 2008, PMC3330053, and
SUPPA) and switches to the rMATS form **only when both effective lengths are supplied (> 0)** —
ASSUMPTION: an API-shape choice; both behaviors are source-backed, and they differ whenever lᵢ ≠ lₛ.

### ΔPSI — differential splicing

Comparing PSI between two conditions gives **ΔPSI** (delta-PSI = Ψ₁ − Ψ₂), the splicing-level analog
of the gene-level log2 fold change: it detects an exon whose *inclusion rate* shifts between groups
(SUPPA2 across conditions; rMATS replicate model). ΔPSI is the natural differential-splicing readout
layered on the per-condition PSI this unit computes.

## 2. Event classification — the five canonical AS classes

An AS event is defined **relative to two isoforms of the same gene** that differ in exon structure
(Wang 2008 — a single isoform defines no event). `DetectAlternativeSplicing` compares the two exon
coordinate lists and names the class (rMATS event codes in parentheses):

| Class | rMATS | How the two isoforms differ |
|-------|-------|-----------------------------|
| **SkippedExon** (cassette exon) | SE | a middle exon is present in one isoform, absent in the other |
| **RetainedIntron** | RI | an intron is retained as exon body in one isoform |
| **AlternativeFivePrimeSS** | A5SS | shared 5′ start, different 3′ end ⇒ alternative **donor** (5′ splice site = exon END) |
| **AlternativeThreePrimeSS** | A3SS | shared 3′ end, different 5′ start ⇒ alternative **acceptor** (3′ splice site = exon START) |
| **MutuallyExclusiveExons** | MXE | exactly one of two alternative middle exons is used |

The A5SS/A3SS distinction hinges on **which boundary moves** (rMATS-turbo coordinate columns, + strand):
alternative **donor** = downstream (3′/END) boundary of the upstream exon; alternative **acceptor** =
upstream (5′/START) boundary of the downstream exon. AS is pervasive — 92–94% of multi-exon human
genes (Wang 2008).

## Invariants and edge cases

- **INV:** `0 ≤ PSI ≤ 1` for any non-negative I, S (a part over the whole).
- **No exclusion reads (S = 0)** → PSI = **1** (fully included); **no inclusion reads (I = 0)** →
  PSI = **0** (fully excluded). Direct from Ψ = I/(I+S).
- **Zero total reads (0/0)** → PSI **undefined → NaN**. rMATS/the Gaussian model add a pseudo-count
  to avoid the division; here an event with no support yields no estimate.
- **Length normalization changes the estimate** whenever lᵢ ≠ lₛ (unnormalized vs rMATS form differ).
- **Fewer than two isoforms** for a gene → **no event**; **identical isoforms** → no structural
  difference → no event.

Worked oracles (from [[trans-splice-001-evidence]]):
- **Unnormalized PSI** I=80, S=20 → `80/(80+20)` = **0.80**.
- **rMATS length-normalized PSI** I=80, S=20, lᵢ=200, lₛ=100 → I/lᵢ=0.40, S/lₛ=0.20 →
  `0.40/(0.40+0.20)` = **0.6666…** (differs from 0.80 because lᵢ ≠ lₛ).
- **Classification** — coordinate pairs for a gene classify as SkippedExon / RetainedIntron /
  AlternativeFivePrimeSS / AlternativeThreePrimeSS / MutuallyExclusiveExons per the table above.

## Relationship to the rest of the RNA-seq family

- **[[differential-expression]]** (TRANS-DIFF-001, the family anchor) — the **gene-level** two-group
  test (log2FC + Welch-t + BH). PSI/ΔPSI is the **splicing-level** counterpart: a gene can be
  differentially *spliced* (an exon's inclusion shifts) with **no** change in its total expression,
  and vice versa — the two units answer complementary questions on the same RNA-seq data.
- **[[expression-quantification]]** (TRANS-EXPR-001) — supplies within-/cross-sample normalized
  expression; PSI is instead a self-normalizing **ratio of reads at one event**, so it does not need
  the TPM/FPKM library-size correction (the inclusion and exclusion reads share the same library).
- **Distinct from the genomic splice family** ([[splice-donor-site-prediction]] /
  [[splice-acceptor-site-prediction]] / [[gene-structure-prediction-intron-exon]]): those score the
  `GT`/`AG` splice-site **motif** from genomic sequence; this counts RNA-seq **reads** to quantify
  how often an already-annotated exon is spliced in. Cross-linked only to disambiguate.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the PSI definition and the
five-class event taxonomy. It is the **read-count ratio + coordinate-comparison** layer — **not** the
full rMATS/SUPPA2 statistical framework (replicate binomial/beta modeling, uncertainty-aware ΔPSI
significance, FDR across events). Two source-backed assumptions: length normalization is opt-in;
forward strand / ascending coordinates. All formulas match their primary sources (Wang 2008 taxonomy;
PMC3330053 / SUPPA Ψ = I/(I+S); Shen 2014 rMATS length-normalized Ψ) with **no source
contradictions** — the sources differ only in whether PSI is length-normalized, which the API exposes
as an option. **Not for clinical use.**
