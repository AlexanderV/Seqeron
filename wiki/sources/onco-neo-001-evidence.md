---
type: source
title: "Evidence: ONCO-NEO-001 (neoantigen candidate peptide window generation — somatic missense mutant/WT agretope windows)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-NEO-001-Evidence.md
sources:
  - docs/Evidence/ONCO-NEO-001-Evidence.md
source_commit: 643a974d3e132ca5be1beedd823ade8c4e535528
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-NEO-001

The validation-evidence artifact for test unit **ONCO-NEO-001** — **Neoantigen Candidate Peptide Window
Generation** for a somatic missense mutation (`GenerateNeoantigenPeptides`). The **twenty-fifth ingested
unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own concept,
[[neoantigen-peptide-generation]]; [[test-unit-registry]] tracks the unit. It is the **upstream partner**
of the affinity gate [[mhc-peptide-binding-prediction]] (ONCO-MHC-001).

## What this file records

Given a mutant protein and a 1-based missense position, the unit returns **every length-`k` window of the
mutant protein that spans the mutated residue**, each paired with its matched **wild-type (germline)**
peptide at the same coordinates and the mutation offset within the window — the mutant/WT **agretope**.
Class I neoantigen search uses lengths **8–11**.

- **Online sources (mutually consistent; no contradictions):**
  - **Hundal et al. (2020)** *Cancer Immunol Res* 8(3):409–420 — **pVACtools**: class I MHC peptides
    **8–11-mer** (class II 13–25); supported variant mechanisms (point mutation, in-frame / frameshift
    indel, fusion); ranking by **agretopicity index** (differential binding of mutant vs matched wild-type),
    filtering to mutant peptides that bind stronger than their WT counterpart (rank 1 / 3).
  - **pVACtools readthedocs** — `flanking_sequence_length` controls residues on either side of the mutation;
    the alteration is **centred in the returned protein sequence "if possible"**; class I benchmarking run
    over lengths **8–11** (rank 3).
  - **Li et al. (2020)** *BMC Med Genomics* 13:52 — **ProGeo-neo**: each coding missense mutation translated
    into a **21-mer** with **10 flanking residues on each side** of the substituted amino acid; the 8–11-mer
    epitope windows overlapping the mutation are extracted from it (rank 1).
  - **Jurtz et al. (2017)** *J Immunol* 199(9):3360–3368 + **NetMHCpan-4.1** DTU service — length options
    8–14-mer; "most presented MHC class I ligands are of length 9 amino acids", preference varies by allele
    → motivates enumerating **8–11** rather than 9 only (rank 1 / 3).
  - **Wells et al. (2020)** *Cell* 183(3):818–834 — **TESLA**: the **differential agretopic index (DAI)**
    (ratio of germline to matched-mutant peptide binding affinity) is a strong point-mutation immunogenicity
    signal; confirms the mutant peptide must be paired with its matched WT at the same coordinates — the
    agretope this unit produces (rank 1).

- **Window enumeration rule (recorded verbatim from the file):** for length `k`, valid 0-based start
  `s ∈ [max(0, p−1−k+1), min(p−1, L−k)]`, count `= last − first + 1`. Interior mutation → exactly `k`
  windows per length; WT window = mutant window with the mutant residue reverted at the mutation offset.

- **Documented corner cases / failure modes:** mutation near a terminus → fewer than `k` windows (truncated,
  built "if possible"); non-coding / synonymous variant → no candidate peptide; length outside 8–11 → not
  enumerated; protein shorter than a requested `k` → that length skipped, shorter windows still returned;
  null protein → `ArgumentNullException`; empty protein / mutant == wildtype / invalid length range /
  out-of-range position → exceptions.

- **Datasets (deterministic hand-derived oracles):**
  - **Interior Y5C** on WT `MKTAYIAKQRSTVWLNDEFGH` (L = 21) → mutant `MKTACIAKQRSTVWLNDEFGH`; each of
    k = 8/9/10/11 yields **5** windows (starts 1..5) → **20** candidate peptides total; e.g. k = 8 mutant
    windows `MKTACIAK, KTACIAKQ, TACIAKQR, ACIAKQRS, CIAKQRST`. Each WT peptide = mutant with `C`→`Y`.
  - **Terminal M1V** (position 1) → k = 9 yields exactly **1** window: mutant `VKTAYIAKQ`, WT `MKTAYIAKQ`,
    mutation offset 0.

- **Coverage recommendations:** MUST — interior missense produces all k∈[8,11] mutation-spanning windows
  (5 per length, 20 total) with correct mutant/WT pairing + offsets; mutant differs from WT only at the
  offset and the WT residue is the original; terminal mutation returns only the fitting (truncated) windows;
  single length k = 9 gives exactly `k` windows for an interior mutation. SHOULD — null / empty / mutant ==
  wildtype / invalid length range / out-of-range position raise exceptions. COULD — protein shorter than a
  requested `k` skips that length but still returns shorter windows.

## Deviations and assumptions

- **ASSUMPTION — single-residue substitution only.** Windowing is implemented for a somatic missense SNV
  (one amino-acid substitution). Frameshift / indel / fusion neopeptides need a different translation step
  (pVACtools supports them separately) and are out of scope — documented, not a correctness gap for the
  missense case.
- **ASSUMPTION — binding affinity out of scope.** IC50 / binding-rank scoring needs a trained MHC model
  (NetMHCpan weights) and is caller-supplied via [[mhc-peptide-binding-prediction]] (ONCO-MHC-001); not
  fabricated here. A source-backed scoping decision, not a correctness assumption.

No source contradictions — pVACtools, ProGeo-neo, NetMHCpan, and TESLA all agree on the
mutation-spanning-window enumeration and the matched-WT agretope pairing.
