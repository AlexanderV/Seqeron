# Seqeron.Genomics — Validated Limitations & Operating Envelope

**Date:** 2026-06-16   **Library:** Seqeron.Genomics (mission-critical)
**Basis:** the independent Phase-1 + Phase-2 validation campaign (234 ☑ Registry units;
see [VALIDATION_LEDGER.md](VALIDATION_LEDGER.md), [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md),
and `docs/Validation/reports/`).

This document consolidates **documented, sourced limitations** — places where an algorithm is a
faithful but *simplified / subset* realisation of a fuller published method, or where it consumes a
caller-supplied input rather than computing it upstream. **None of these are defects:** every item
below was validated as correct *for its stated contract* and is `BY-DESIGN` in the findings register.
The point of this file is to make the library's honest operating envelope visible in one place so
downstream users are not surprised.

Genuine algorithm defects found during validation were **all fixed in-session** (13 in Phase 2);
they are recorded in the ledger, not here. Deferred *enhancements* that would lift some of these
limitations are tracked separately as the "Deferred BIG fixes" backlog in the ledger and the
`C. NOT-POSSIBLE (radical)` section of the findings register.

---

## 1. Algorithmic simplifications (a fuller method exists; the simpler one is correct for its scope)

| Unit(s) | Implemented | Fuller method not implemented | Note |
|---------|-------------|-------------------------------|------|
| ALIGN-MULTI-001 | Star (center-star) MSA **and** guide-tree progressive MSA | Iterative refinement (consistency-based / MUSCLE-style re-alignment) | **Both implemented (C4).** Star `MultipleAlign` is byte-for-byte unchanged; the added `MultipleAlignProgressive` is the Feng-Doolittle method (pairwise-NW identity distances → UPGMA guide tree → profile–profile NW with "once a gap, always a gap"). Progressive alignment remains single-pass and inherits the standard "once a gap, always a gap" limitation (early gap-placement errors are not later corrected). |
| RESTR-DIGEST-001 | Linear **and** circular-molecule digest | — | **Both topologies now implemented (C5).** `Digest(DnaSequence, params string[])` (linear) is byte-for-byte unchanged; the added `Digest(DnaSequence, MoleculeTopology, params string[])` overload with `MoleculeTopology { Linear, Circular }` handles closed circles — `k` sites → `k` fragments with a wrap-around origin-spanning fragment (length `(len − lastCut) + firstCut`), 0 sites → one uncut circular fragment, 1 site → one linearized fragment. Default stays `Linear`. No remaining gap. |
| RNA-STRUCT-001 | Greedy dot-bracket structure | Nussinov/Zuker DP traceback structure | The DP **energy/score values** themselves are correct; only the returned structure is the greedy approximation. |
| PHYLO-COMP-001 | Rooted-clade **and** unrooted-bipartition Robinson–Foulds | — | **Both metrics now implemented (C3).** `RobinsonFouldsDistance` (rooted clades) is unchanged; `CalculateUnrootedRobinsonFoulds` / `CalculateNormalizedUnrootedRobinsonFoulds` add the original Robinson & Foulds (1981) split-symmetric-difference metric, which is root-invariant. No remaining gap. |
| Phylogenetics (PHYLO-NEWICK/TREE/COMP) | Strictly **binary** trees (`PhyloNode` Left/Right) | True N-ary (multifurcating) trees | Newick now **throws** instead of silently truncating a multifurcation (fixed mid-campaign). |
| META-CLASS-001 | **Kraken taxonomy-tree + k-mer-LCA + RTL classifier** | — | **Now implemented (C1, approved breaking change); no remaining gap.** Replaced the former flat best-hit classifier with the faithful Kraken algorithm (Wood & Salzberg 2014): a `TaxonomyTree` with LCA, a `BuildKmerDatabase` that maps each shared canonical k-mer to the **LCA of its owning taxa**, and a `ClassifyReads` that assigns the leaf of the **maximum-scoring root-to-leaf path** (LCA-of-leaves on ties), no-hit → root/unclassified. Confidence is Kraken's C/Q. The unit is no longer `🔧 LIMITED`. |
| CRISPR-GUIDE/OFF-001 | **Doench 2014 "Rule Set 1"** on-target linear model **and** **MIT/Hsu 2013** single-hit + aggregate off-target specificity score (both real published models, C7) | Doench **"Rule Set 2" / Azimuth** (on-target) and **CFD** (off-target) | **Two published models now implemented (C7).** `CalculateOnTargetDoench2014` reproduces the Doench et al. 2014 (Nat Biotechnol 32:1262) logistic linear model over the 30-nt context — intercept + GC-count term + per-position single/di-nucleotide weights, coefficients transcribed verbatim from the reference `doenchScore.py`; reproduces that file's own worked examples (`TATAGCTGCGATCTGAGGTAGGGAGGGACC` → 71.309, `TCCGCACCTGTCACGGTCGGGGCTTGGCGC` → 1.898) to <1e-4. `CalculateMitHitScore` / `CalculateMitSpecificityScore` reproduce the Hsu et al. 2013 (Nat Biotechnol 31:827) single-hit (Π(1−W[i]) × distance term × 1/nmm²) and aggregate (100/(100+Σ)·100) formulas with the published 20-position W vector. The pre-existing honest-heuristic `EvaluateGuideRna` / `FindOffTargets` / `CalculateSpecificityScore` are unchanged and additive. **Residual:** (1) **Rule Set 2 / Azimuth is a gradient-boosted-tree model, NOT a coefficient table** — it cannot be faithfully reproduced from published numbers without the trained model file, so it is deliberately NOT implemented (no fabrication). (2) **CFD** (Doench 2016) is NOT implemented this session: its ~400-entry mismatch + 16-entry PAM penalty matrix is distributed only as a binary pickle (`mismatch_score.pkl` / `pam_scores.pkl`) and was not cleanly retrievable as authoritative verbatim text this session; deferred rather than approximated. |
| PRIMER-TM-001 | Honest heuristic / NN-thermodynamics | — | NN-thermodynamics (ΔH/ΔS/Tm) **is** validated (SEQ-THERMO-001); the primer-design convenience scoring is heuristic. |

## 2. "Threshold / aggregation / framework" layers — they classify or combine caller-supplied inputs, they do not predict upstream

These units are correct implementations of a published **rule or formula**, but they sit *downstream*
of a model/measurement the caller must provide. They are **decision-support computations, not
validated clinical-grade predictors.**

| Unit | Computes (validated) | Caller must supply / not modelled |
|------|----------------------|-----------------------------------|
| ONCO-SIG-002/003/004 | NNLS signature **refitting**, bootstrap CIs, aetiology mapping | **De-novo signature extraction (NMF)** is not implemented — reference signatures are an input. |
| ONCO-MHC-001 | Strong/weak-binder **classification** by IC50 / %Rank cutoffs | The **affinity / %Rank prediction** (NetMHCpan-style learned model) — supplied by caller. |
| ONCO-HRD-001 | HRD score = LOH + TAI + LST, ≥42 cutoff | The three component scores are **inputs** (per-segment derivation is ONCO-LOH/CNA). |
| ONCO-PURITY/PLOIDY/CCF/CLONAL | Purity/ploidy/CCF formulas + clonal rule | Allele-specific CN segments, multiplicity, VAF — supplied; ploidy-WGD uses supplied-segment length, not a chromosome-size table. |
| EPIGEN-AGE-001 | Horvath `anti.trafo` + linear predictor | The **353-CpG coefficient table** is a caller-supplied input (framework, no fabricated coefficients). |
| ONCO-CHIP-001 | CHIP filter (gene panel + ≥2% VAF + WBC subtraction) | Origin call uses a **gene+VAF heuristic** where matched-WBC data is absent (over-removes vs strict matched-WBC origin). |
| ONCO-MRD-001 | ≥2-of-N positivity, IMAF, Poisson LoD | IMAF is read-pooled, **without** INVAR-style background subtraction / tumour-AF weighting. |
| ONCO-SIG-003 | Bootstrap exposure CIs | Uses sigminer fixed-N **multinomial** resample (not Senkin's Poisson variant). |
| COMPGEN-ANI-001 | Goris ANI mean-identity formula | **Ungapped, single-direction** fragment placement (not gapped BLASTN, not reciprocal). |

## 3. Convention divergences (documented, internally consistent)

| Area | Convention used | Differs from |
|------|-----------------|--------------|
| Variant units (VARIANT-*, SV-*, ONCO somatic) | 0-based internal `Position` | VCF 1-based (1-based only on serialization, e.g. `ToVcfLines`). |
| GC / skew / composition outputs | Percentage (×100) in several units | Biopython fraction [0,1]. |
| Default thresholds | Configurable defaults that may differ from a named tool's published default (e.g. MSI 20% (MSIsensor2) vs MSIsensor ~3.5; ConsensusThreshold 0.5 vs Biopython 0.7; CNA inner cutoffs) | The named tool — but the value is sourced and the parameter is exposed. |
| DUST / compression complexity | `k=3` (DUST) and the LZ base-2 clamp are exactly sourced | `k≠3` is a documented, non-exact-asserted extrapolation. |
| Non-ACGT / ambiguous bases | Skipped or tracked as "Other" in several composition/thermo units | Tools that throw or partially count (e.g. Biopython S/W handling). |

---

## How to read this

- For **research / pipeline** use, these limitations are the normal scope boundaries of a from-first-
  principles library and are individually sourced.
- For **clinical or decision-grade** use, treat §2 especially as **heuristics / threshold layers that
  require an external validated predictor and clinical sign-off** — the library computes the rule, not
  the trained model behind it.
- Each row traces to a per-unit validation report (`git show <validate-commit>:docs/Validation/reports/<ID>.md`)
  and a `BY-DESIGN` / PASS-WITH-NOTES entry in the findings register.
