---
type: source
title: "Evidence: ONCO-FUSION-003 (fusion breakpoint reading-frame consequence + fusion protein prediction)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-FUSION-003-Evidence.md
sources:
  - docs/Evidence/ONCO-FUSION-003-Evidence.md
source_commit: 7b40b36e8b46c1f03926ea3cca120370653e3af0
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-FUSION-003

The validation-evidence artifact for test unit **ONCO-FUSION-003** — **Fusion Breakpoint
Analysis** (junction reading-frame consequence + fusion protein prediction). The **sixteenth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized
in [[fusion-breakpoint-frame-and-protein-prediction]]; it is the **protein-consequence** third
member of the fusion trio, downstream of the read-evidence caller
[[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001) and orthogonal to the naming unit
[[gene-fusion-nomenclature-known-fusion-lookup]] (ONCO-FUSION-002). Both siblings explicitly
deferred the premature-stop / transcript-reconstruction scope to this unit.
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (two reference implementations + the codon-triplet basis):**
  - **Arriba output-files spec** (Uhrig et al. 2021, *Genome Research*, rank 3) — the
    `reading_frame` column states whether the 3′ gene is fused **in-frame** or **out-of-frame**,
    with four values `in-frame` / `out-of-frame` / `stop-codon` / `.` (dot = peptide not
    predictable); `site1` / `site2` place each breakpoint in `5'UTR` / `3'UTR` / `UTR` / `CDS` /
    `exon` / `intron` / `intergenic` (+ `splice-site` at an exon boundary); `fusion_transcript`
    and `peptide_sequence` carry the assembled transcript / translated peptide with `|` marking
    the breakpoint; `type` = translocation / duplication / inversion / deletion / read-through.
  - **AGFusion `model.py`** (Murphy & Elemento 2016, rank 3) — the reference **protein
    prediction**: 5′ partner contributes its CDS **prefix** up to the breakpoint offset
    (`coding_sequence[0:junction_5prime]`), 3′ partner contributes its CDS **suffix** from the
    breakpoint (`coding_sequence[junction_3prime:]`), the two are **concatenated**, **translated**,
    and **truncated at the first stop codon** (`protein_seq[0:find("*")]`); an out-of-frame fusion
    first trims the chimeric CDS to whole codons (`seq[0:3*int(len/3)]`) before translating.
  - **Wikipedia "Reading frame"** (citing Badger & Olsen 1999, rank 4) — the codon-triplet /
    modulo-3 basis for the junction-phase test.

- **The two models, and the repo's choice:** AGFusion uses a **three-way** class — `in-frame`
  (both segments codon-aligned) / `in-frame (with mutation)` (segment fractional codon parts
  complement so the contiguous ORF length is a multiple of 3, but the 3′ gene is read
  frameshifted) / `out-of-frame`. The repo instead models **Arriba's two-way** call: `InFrame` /
  `OutOfFrame` / `StopCodon` / `NotPredicted`, where in-frame means the **3′ gene is read in its
  NATIVE reading frame**. So AGFusion's `in-frame (with mutation)` maps to the repo's
  **`OutOfFrame`** (the 3′ gene is frameshifted, so its protein domains are not preserved). The
  frame test is therefore the exon-phase-compatibility rule
  **`(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0`** — the same rule ONCO-FUSION-001
  uses, sourced there to Genomics England / Wikipedia.

- **Documented corner cases / failure modes:**
  - **Peptide not predictable** → `reading_frame = .`: a frame status is only defined when both
    breakpoints sit in coding context.
  - **Breakpoint not in CDS** (5′UTR / 3′UTR / intron / non-coding exon / intergenic): the junction
    does not join two coding frames, so no in/out call is made.
  - **Premature stop at junction** → `reading_frame = stop-codon`: the junction or downstream
    sequence yields a stop codon terminating the chimeric ORF; even a numerically **in-frame**
    fusion with an early stop yields a **truncated** peptide (translation always terminates at the
    first `*`).
  - **Out-of-frame frameshift**: the 3′ partner is read in a shifted frame; the CDS is trimmed to
    whole codons and translated, typically truncated early at the first downstream stop.

- **Datasets (deterministic, hand-derived from the AGFusion rule + standard genetic code
  transl_table=1):**
  - **Junction reading-frame truth table** over `b` (5′ coding bases) and `p` (3′ start phase):
    `(b−p) mod 3 == 0` → in-frame (9/0, 10/1, 11/2), else out-of-frame (10/0, 9/1, 9/2).
  - **Fusion protein prediction** worked vectors: `ATGAAA|GATGGT` → `MKDG` (in-frame, 6%3=0/6%3=0);
    `ATGAAA|GATTAAGGT` → `MKD` (in-frame, premature stop downstream of junction); `ATGA|AAGGT`
    (5′ 4 bases phase 1, 3′ suffix native phase 0) → `(4−0) mod 3 = 1` → **out-of-frame** under
    the Arriba 3′-native-frame model even though the 9-base chimeric CDS translates cleanly to
    `MKG` (this is exactly where the two models diverge); same bases with the 3′ suffix at native
    phase 1 → `(4−1) mod 3 = 0` → in-frame; `ATGAA|GATGGT` (5%3=2, 6%3=0) → out-of-frame, trim to
    9 → `MKM`.

- **Coverage recommendations (8 items):** MUST — in-frame call across all three residues of p;
  out-of-frame when `(b−p) mod 3 != 0`; breakpoint-site classification (CDS vs UTR/intron/intergenic)
  drives whether a call is made (non-CDS → not predicted); `PredictFusionProtein` concatenate +
  translate + truncate at first stop (`ATGAAA|GATGGT`→`MKDG`, `ATGAAA|GATTAAGGT`→`MKD`);
  out-of-frame trims to whole codons and reads the 3′ partner in a shifted frame. SHOULD —
  premature-`stop-codon` flag; null/empty/invalid-phase (outside {0,1,2}) raise the documented
  exceptions. COULD — stable/deterministic record shape.

## Deviations and assumptions

- **MODEL — Arriba two-way call, not AGFusion's three-way class.** `BreakpointFrameStatus` is
  `InFrame` / `OutOfFrame` / `StopCodon` / `NotPredicted`; AGFusion's `in-frame (with mutation)`
  (contiguous ORF multiple of 3 but 3′ gene frameshifted) maps to **`OutOfFrame`**. A deliberate
  model choice — a frameshifted 3′ partner does not preserve the functional (e.g.
  kinase-domain-preserving) protein, so it is not "in-frame" in the biological sense — not a defect.
- **ASSUMPTION — caller supplies CDS sequences and junction offsets.** The repository has no
  genome/GTF/transcript database, so `PredictFusionProtein` takes the partner CDS strings and the
  breakpoint CDS offsets directly (5′ coding-base count + 3′ coding-start offset) rather than
  resolving them from an Ensembl/RefSeq annotation as AGFusion does. The concatenation, frame rule,
  translation, and first-stop truncation are exactly AGFusion's; only the input source differs — no
  computed output changes for given CDS/offsets.

No source contradictions — Arriba and AGFusion agree on concatenate/translate/truncate; the
two-way-vs-three-way framing is a documented, deliberate model selection.
