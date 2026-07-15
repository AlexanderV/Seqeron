---
type: concept
title: "Fusion breakpoint reading-frame consequence + fusion protein prediction"
tags: [oncology, structural-variant, algorithm]
sources:
  - docs/Evidence/ONCO-FUSION-003-Evidence.md
  - docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md
source_commit: 5465dd6bff485c29023a0d3c5020b0dbb6a04c62
created: 2026-07-10
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:gene-fusion-detection-read-evidence
      source: onco-fusion-003-evidence
      evidence: "ONCO-FUSION-001 defers the premature-stop / transcript-reconstruction scope to ONCO-FUSION-003; both share the exon-phase frame rule (5' coding bases − 3' start phase) mod 3 == 0 (ONCO-FUSION-003-Evidence.md Assumption 0)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-fusion-003-evidence
      evidence: "Test Unit ID: ONCO-FUSION-003, Algorithm: Fusion Breakpoint Analysis (junction reading-frame consequence + fusion protein prediction)"
      confidence: high
      status: current
---

# Fusion breakpoint reading-frame consequence + fusion protein prediction

The **sixteenth ingested Oncology unit** (ONCO-FUSION-003) and the third member of the fusion
trio: the **protein-consequence** layer that both siblings explicitly deferred to. Given a
detected fusion's breakpoint it decides the **reading-frame status** of the chimeric transcript
and reconstructs the **fusion protein**. Validated under test unit **ONCO-FUSION-003**
([[onco-fusion-003-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

Where [[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001) decides *whether* a fusion is
detected — and stops at a binary exon-phase in/out-of-frame check, explicitly **not** scanning
for a premature stop — this unit computes the full **functional consequence**. The evidence model
is **Arriba** (Uhrig 2021, the `reading_frame` / `site` / `peptide_sequence` output columns) and
**AGFusion** (Murphy & Elemento 2016, the `model.py` protein-prediction reference); see
[[onco-fusion-003-evidence]] for the source-by-source trace.

## Four-state reading-frame status

`BreakpointFrameStatus` is a **four-state** call, richer than ONCO-FUSION-001's binary result:

- **`InFrame`** — the 3′ gene is read in its **native** reading frame.
- **`OutOfFrame`** — the 3′ gene is read frameshifted (its protein domains are not preserved).
- **`StopCodon`** — the junction or downstream sequence yields a premature stop terminating the
  chimeric ORF (Arriba's `stop-codon` value).
- **`NotPredicted`** — the peptide cannot be predicted (Arriba's `.`), e.g. a breakpoint outside a
  partner's coding region.

The in/out decision reuses the exon-phase rule from ONCO-FUSION-001 (reading frames are triplets):

```
in-frame  ⇔  (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0
```

with `fivePrimeCodingBases` = the coding bases the 5′ partner contributes up to the breakpoint and
`threePrimeStartPhase ∈ {0,1,2}` = the 3′ partner's coding-start phase at the breakpoint.

## Arriba two-way model, not AGFusion's three-way class

The unit deliberately follows **Arriba's two-way** in/out call — "whether the gene at the 3′ end
is fused in-frame or out-of-frame", i.e. whether the 3′ gene is in its **native** frame — rather
than **AGFusion's three-way** class (`in-frame` / `in-frame (with mutation)` / `out-of-frame`).
AGFusion's middle class labels a fusion whose two segment lengths' fractional codon parts
complement so the *contiguous* chimeric ORF length is a multiple of 3, **but the 3′ gene is still
read frameshifted**; the repo maps that case to **`OutOfFrame`**, because a frameshifted 3′ partner
does not preserve the functional (e.g. kinase-domain) protein. This is a documented model choice,
not a defect. The divergence is exactly the `ATGA|AAGGT` oracle below: `(4−0) mod 3 = 1` →
`OutOfFrame` here even though the 9-base chimeric CDS translates cleanly.

## Breakpoint-site classification gates the call

A frame status is only defined when **both** breakpoints sit in coding context. Each breakpoint is
placed as `5'UTR` / `3'UTR` / `UTR` / `CDS` / `exon` / `intron` / `intergenic` (Arriba `site1` /
`site2`, `+splice-site` at an exon boundary). A breakpoint in a UTR, intron, non-coding exon, or
intergenic region does **not** join two coding frames → status `NotPredicted`.

## Fusion protein prediction (`PredictFusionProtein`)

Following AGFusion's `model.py` exactly:

1. **5′ CDS prefix** — the 5′ partner contributes its coding sequence from the start up to the
   breakpoint offset (`cds[0:junction_5prime]`).
2. **3′ CDS suffix** — the 3′ partner contributes its coding sequence from the breakpoint offset to
   the end (`cds[junction_3prime:]`).
3. **Concatenate** the two segments into the chimeric CDS.
4. **Translate** with the standard genetic code (NCBI transl_table=1).
5. **Truncate at the first stop codon** (`protein[0:find("*")]`) — so even an in-frame fusion with
   an early stop yields a truncated peptide.

For an **out-of-frame** fusion the chimeric CDS is first **trimmed to a whole number of codons**
(`cds[0:3*int(len/3)]`) before translating, and the 3′ partner is read in its shifted frame.

Worked oracles (from [[onco-fusion-003-evidence]], standard genetic code):

- `ATGAAA | GATGGT` → `MKDG` (in-frame; 6%3=0, 6%3=0).
- `ATGAAA | GATTAAGGT` → `MKD` (in-frame, premature stop downstream of the junction → truncated).
- `ATGA | AAGGT`, 3′ suffix at native phase 0 → `MKG` translates cleanly but is **`OutOfFrame`**
  ((4−0) mod 3 = 1 — the Arriba-vs-AGFusion divergence point).
- `ATGA | AAGGT`, 3′ suffix at native phase 1 → **in-frame** ((4−1) mod 3 = 0).
- `ATGAA | GATGGT` → out-of-frame (5%3=2, 6%3=0), trim to 9 → `MKM`.

## Implementation surface

The algorithm spec (`Fusion_Breakpoint_Analysis.md`, ONCO-FUSION-003) pins the two entry
points on `OncologyAnalyzer` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`):

- **`AnalyzeBreakpoint(FusionBreakpoint)`** → `BreakpointAnalysis { BreakpointInCoding, FrameStatus }`.
  `BreakpointInCoding` is true iff **both** sites are `CDS`; the frame call is delegated to the shared
  `IsInFrame(int, int)` (ONCO-FUSION-001) and is evaluated **only** when both breakpoints are coding,
  else `NotPredicted`. `O(1)` — constant-time site/frame checks.
- **`PredictFusionProtein(FusionBreakpoint, (string FivePrimeCds, string ThreePrimeCds))`** →
  `FusionProteinPrediction { ChimericCds, Peptide, Effect, HasPrematureStop }`. The 5′ prefix length
  is `FivePrimeCodingBases` and the 3′ suffix start is `ThreePrimeStartPhase`, both taken from the
  `FusionBreakpoint`. `O(n)` in the chimeric CDS length — one pass over codons. Translation reuses
  `Seqeron.Genomics.Core.GeneticCode.Standard.Translate` (NCBI table 1); the codon table is **not**
  duplicated in the oncology layer. `HasPrematureStop` is true iff a stop is reached before the ORF end.

**Validation / normalization.** Inputs are case-normalized to uppercase DNA. CDS strings must be
non-null (`ArgumentNullException`); offsets must lie within their CDS and `ThreePrimeStartPhase ∈ {0,1,2}`
for a protein prediction (`ArgumentOutOfRangeException`). Translation is 0-based, frame 0, read
5′→3′ in triplets — a trailing partial codon is not translated. `Implementation Status: Framework`:
the repo bundles no genome/GTF, so partner CDS and offsets are caller-supplied and no substring search
(hence no suffix tree) is involved — the chimeric CDS is built purely by slice/concat.

## Scope and relation to the fusion trio

A [[research-grade-limitations|research-grade]] method, **not for clinical use**. It is the
**protein-consequence** third member of the fusion trio: downstream of the read-evidence caller
[[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001), whose binary in-frame check it elaborates
into a four-state consequence with actual transcript/peptide reconstruction, and orthogonal to the
HGNC naming unit [[gene-fusion-nomenclature-known-fusion-lookup]] (ONCO-FUSION-002). Because the
repository bundles no genome/GTF/transcript database, `PredictFusionProtein` takes the partner CDS
strings and breakpoint offsets **directly** from the caller (an API-shape assumption) — the
concatenation, frame rule, translation, and first-stop truncation are exactly AGFusion's; only the
input source differs. Orthogonal to the copy-number / clonal-structure and clinical-interpretation
ONCO units that would consume a called and characterized fusion. No source contradictions — the
two-way-vs-three-way framing is a documented, deliberate model selection.
