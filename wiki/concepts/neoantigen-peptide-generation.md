---
type: concept
title: "Neoantigen candidate peptide window generation (somatic missense → mutant/WT agretope windows)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-NEO-001-Evidence.md
  - docs/algorithms/Oncology/Neoantigen_Peptide_Generation.md
source_commit: ce89ed9dbabb6aab5d19f8c05bd6b602f7a50b7e
created: 2026-07-10
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-neo-001-evidence
      evidence: "Test Unit ID: ONCO-NEO-001, Algorithm: Neoantigen Candidate Peptide Window Generation (somatic missense mutation)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:mhc-peptide-binding-prediction
      source: onco-neo-001-evidence
      evidence: "Binding affinity out of scope: IC50 / binding-rank scoring requires a trained MHC model (NetMHCpan weights) and is caller-supplied (ONCO-MHC-001); the mutant peptide is a candidate only if MHC scores it a sufficiently strong binder."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:hla-nomenclature-and-allele-specific-loh
      source: onco-neo-001-evidence
      evidence: "Wells 2020 (TESLA): the differential agretopic index pairs the mutant peptide with its matched germline counterpart at the same coordinates — the agretope this unit produces — which is presented on a tumour HLA allele; HLA LOH removes that presentation platform."
      confidence: medium
      status: current
---

# Neoantigen candidate peptide window generation

The **twenty-fifth ingested Oncology unit** (**ONCO-NEO-001**, `GenerateNeoantigenPeptides`) and the
**upstream partner** of the affinity gate [[mhc-peptide-binding-prediction]] (ONCO-MHC-001). Given a
somatic **missense SNV** (a single amino-acid substitution), it enumerates the tumour-specific candidate
peptides that must then be tested for MHC presentation. The literature-traced record is
[[onco-neo-001-evidence]]; [[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]]
describes the evidence-artifact pattern.

## What the unit produces

For a mutant protein and a mutation at 1-based position `p`, and for each requested peptide length `k`,
it returns **every length-`k` window of the mutant protein that spans the mutated residue**, each paired
with its matched **wild-type (germline)** peptide at the same coordinates plus the offset of the mutation
within the window. The mutant/WT pair at identical coordinates is the **agretope** — the object on which
downstream immunogenicity signals such as the **differential agretopic index (DAI)** are computed (Wells
et al. 2020, TESLA; the DAI is the ratio of the germline peptide's binding affinity to the matched mutant
peptide's).

Class I neoantigen search uses lengths **8–11** (pVACtools; NetMHCpan supports 8–14, most presented
ligands are 9-mers, so all of 8–11 are enumerated rather than 9 alone). A **21-mer with the substitution
centred** (10 flanking residues each side, ProGeo-neo's `±10` construction) contains exactly every 8–11-mer
window that overlaps the mutated residue — this unit produces those windows directly.

## Window enumeration rule

For peptide length `k`, protein length `L`, and 1-based mutation position `p`, the valid **0-based** start
offsets are:

```
s ∈ [ max(0, p−1−k+1),  min(p−1, L−k) ]        count = last − first + 1
```

Every such window includes the mutated residue. For a **sufficiently interior** mutation the count is
exactly `k` (there are `k` distinct offsets at which a length-`k` window can cover one fixed position). The
matched wild-type peptide equals the mutant peptide with the mutant residue reverted to the wild-type
residue at the window's mutation offset.

### Worked example (interior missense)

Wild-type `MKTAYIAKQRSTVWLNDEFGH` (L = 21), missense **Y5C** → mutant `MKTACIAKQRSTVWLNDEFGH`.

| k | first start (1-based) | last start (1-based) | window count | example mutant windows |
|---|---|---|---|---|
| 8 | 1 | 5 | 5 | `MKTACIAK` … `CIAKQRST` |
| 9 | 1 | 5 | 5 | `MKTACIAKQ` … `CIAKQRSTV` |
| 10 | 1 | 5 | 5 | `MKTACIAKQR` … `CIAKQRSTVW` |
| 11 | 1 | 5 | 5 | `MKTACIAKQRS` … `CIAKQRSTVWL` |

Total over k = 8..11 → **20** candidate peptides; each WT peptide equals its mutant with `C`→`Y` at the
mutation offset.

### Corner cases

- **Terminal mutation (truncated windows):** when the substitution is within fewer than `k−1` residues of
  either protein end, fewer than `k` windows of length `k` fit. Example **M1V** at k = 9 → exactly **one**
  window, mutant `VKTAYIAKQ` / WT `MKTAYIAKQ`, mutation offset 0. The set is whichever windows still fit
  inside the protein while spanning the mutation (ProGeo-neo / pVACtools centre the peptide "if possible").
- **No amino-acid change:** a non-coding or synonymous variant produces **no** candidate peptide — peptides
  come only from protein-altering variants.
- **Length outside 8–11:** outside the class I canonical search band; not enumerated. A protein shorter
  than some requested `k` simply skips that length and still returns the shorter windows.

## Scope and relation to the neoantigen pipeline

- **Missense-only.** This unit implements windowing for a single-residue substitution. Frameshift / indel /
  fusion neopeptides need a different translation step (pVACtools supports them separately) and are out of
  scope — documented, not a correctness gap for the missense case.
- **Binding affinity is out of scope.** IC50 / %Rank scoring needs a trained MHC model (NetMHCpan /
  MHCflurry weights) and is **caller-supplied** via [[mhc-peptide-binding-prediction]] (ONCO-MHC-001), the
  downstream **affinity gate**: a mutant window is an actual neoantigen candidate only once MHC scores it a
  sufficiently strong binder to one of the tumour's HLA alleles. ONCO-MHC-001's resolved class I length
  window (8–14, `MhcClassIMaxPeptideLength`) propagates to this generator's default length range.
- [[hla-nomenclature-and-allele-specific-loh]] (ONCO-HLA-001) is the presentation-platform sibling: HLA LOH
  **removes** an allele, so any neoantigen restricted to a lost allele can no longer be presented — immune
  escape — even though this unit would still generate the peptide.

## Implementation surface

The algorithm spec (`Neoantigen_Peptide_Generation.md`, unit **ONCO-NEO-001**, status **Framework**)
fixes the API contract in `OncologyAnalyzer.cs`:

- Entry point `OncologyAnalyzer.GenerateNeoantigenPeptides(string wildTypeProtein, char mutantResidue,
  int mutationPosition, int minLength = 8, int maxLength = 11)` → `IReadOnlyList<NeoantigenPeptide>`
  ordered by **length ascending, then start ascending** (INV-06).
- `NeoantigenPeptide` is a record struct with fields `Length`, `StartPosition` (1-based),
  `MutantPeptide`, `WildTypePeptide` (the agretope at identical coordinates), and `MutationOffset`
  (0-based offset of the substituted residue in the window). Constants `MhcClassIMinPeptideLength` /
  `MhcClassIMaxPeptideLength` = 8 / 11 supply the class I defaults.
- **Validation / exception contract:** null protein → `ArgumentNullException`; empty protein,
  `mutantResidue` equal to the wild-type residue (not a real substitution), `minLength < 1`, or
  `maxLength < minLength` → `ArgumentException`; `mutationPosition` outside `[1, L]` →
  `ArgumentOutOfRangeException`. A requested `k > L` is silently skipped; sequences are opaque
  one-letter strings (no alphabet validation, case preserved).
- Six invariants `INV-01`…`INV-06` (length in band, every peptide spans the mutation, equal-length
  mutant/WT pair differing only at the offset, mutant residue placed at the offset, interior-mutation
  count = k, deterministic ordering) are covered by
  `OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs`.
- **Deviation:** no suffix tree — "find windows spanning a position" is a bounded arithmetic range over
  one short protein, not a multi-query exact-match search, so the repository suffix tree does not apply
  (correctness unaffected). The class is named `OncologyAnalyzer` (project layout), superseding the
  checklist's `NeoantigenPredictor` placeholder.

## Sources and rigor

Source-traced against **pVACtools** (Hundal et al. 2020, class I 8–11-mer windowing + agretopicity),
**ProGeo-neo** (Li et al. 2020, the 21-mer ±10-flank construction), **NetMHCpan-4.0/4.1** (Jurtz et al.
2017 / DTU service, length band + 9-mer dominance), and **Wells et al. 2020** (TESLA, DAI agretope
pairing). The worked windows are derived by hand from the windowing definition, not from code output. A
[[scientific-rigor|research-grade]] correctness reference, **not for clinical or diagnostic use**. No
source contradictions — all four references agree on the mutation-spanning-window + matched-WT-agretope
framing.
