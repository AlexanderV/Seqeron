---
type: concept
title: "Protein domain prediction (PROSITE pattern + Plan7 profile-HMM) & signal-peptide prediction"
tags: [analysis, algorithm, protein, motif]
sources:
  - docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
source_commit: d7b350cc67d3684eb247a261ab387137e91e7fba
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-domain-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-DOMAIN-001 ... Algorithm: Domain Prediction & Signal Peptide Prediction (FindDomains / FindDomainsByHmm / PredictSignalPeptide)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:common-protein-motifs
      source: protmotif-domain-001-evidence
      evidence: "Both are ProteinMotif-family units on ProteinMotifFinder scanning an amino-acid sequence against PROSITE signatures; FindCommonMotifs matches a fixed short-motif dictionary while FindDomains matches PROSITE domain patterns (PS00028/PS00017/PS00678) and adds an opt-in Plan7 profile-HMM engine for profile-only families (SH3/PDZ/WD40) — distinct algorithms in one family"
      confidence: medium
      status: current
---

# Protein domain prediction & signal-peptide prediction

Test unit **PROTMOTIF-DOMAIN-001** validates two related protein-feature capabilities on
`ProteinMotifFinder`: **domain prediction** (`FindDomains`, plus an opt-in profile-HMM engine
`FindDomainsByHmm`) and **signal-peptide prediction** (`PredictSignalPeptide`). It is a unit of the
**ProteinMotif family**, a sibling of the fixed-dictionary [[common-protein-motifs]] and the windowed
[[coiled-coil-prediction]]. The validation record is [[protmotif-domain-001-evidence]] and
[[test-unit-registry]] tracks the unit; see [[algorithm-validation-evidence]] for the artifact pattern.

## `FindDomains` — deterministic PROSITE-pattern domain scan

A **protein domain** is a self-stabilizing, independently folding region (~50–250 aa; Wetlaufer 1973).
`FindDomains` detects domains that have an exact **PROSITE PATTERN** (a regex-like signature, not a
profile/HMM), so the match is deterministic and reproducible. Each hit carries the domain name,
PROSITE accession, Pfam cross-reference, and 0-based inclusive start/end.

| Domain | PROSITE | Pattern (verbatim) | Pfam |
|--------|---------|--------------------|------|
| Zinc finger C2H2 | PS00028 | `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` | PF00096 |
| P-loop / Walker A (ATP/GTP-binding motif A) | PS00017 | `[AG]-x(4)-G-K-[ST]` | PF00069 |
| WD40 repeat (Trp-Asp) | PS00678 | `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]` (14 elements over **15** residues) | PF00400 |

PROSITE → regex translation follows the ScanProsite rules: `x`→`.`, `[..]` kept, `{X}`→`[^X]`,
`e(n)`/`x(n)`→`{n}`, `x(m,n)`→`{m,n}`, `<`→`^`, `>`→`$`, separators dropped. Real positive control:
the PS00678 regex matches **GBB1_HUMAN (UniProt P62873)**, a 7-bladed WD40 β-propeller, at 0-based
starts 69/156/284 (each a 15-residue window; e.g. `LVSASQDGKLIIWDS`).

**Honest residual (profile-only families):** **SH3 (PS50002)** and **PDZ (PS50106)** are PROSITE
**PROFILEs** (weight matrices), not patterns — there is **no** deterministic signature to reproduce as
a regex. They are therefore **not** detected by `FindDomains`; the previously shipped ad-hoc SH3/PDZ
regexes (unsourced) were removed. These families are instead covered by the Plan7 engine below.

## `FindDomainsByHmm` — opt-in Plan7 profile-HMM engine (SH3, PDZ, WD40)

To cover the profile-only families, an **opt-in** Plan7 profile-HMM scorer (`Plan7ProfileHmm` +
`ProteinMotifFinder.FindDomainsByHmm` / `ScoreDomainHmm`) reproduces the **HMMER3** engine over three
bundled **CC0 (public-domain) Pfam** profiles: `PF00018.35` SH3 (LENG 48), `PF00595.30` PDZ (LENG 81),
`PF00400.39` WD40 (LENG 39). The exact-PROSITE-pattern `FindDomains` path and its defaults are unchanged.

The engine was built up across several sessions, each verified against **pyhmmer 0.12.1** ground truth
and an independent from-scratch Python re-derivation of the retrieved HMMER recurrences:

- **Viterbi / Forward log-odds** (Durbin et al. §5.4 recurrences; insert emissions hardwired to
  background → insert log-odds ≈ 0). Verified exact (1e-9) on a hand-built 2-symbol HMM
  (`ln(0.7/0.6)+ln0.9+ln(0.8/0.4)+ln0.8 = 0.5187937934…` nats).
- **hmmsearch-parity `pre_score`** — HMMER's **local-multihit** Forward bit score, scoring emissions
  against the Swiss-Prot **`bg->f`** background (not the profile `COMPO` line) → SH3 `pre = 68.7097`,
  PDZ `84.8629`, WD40 `213.4120` bits, all within ~1e-5 bits of hmmsearch.
- **null2 biased-composition correction** (`p7_GNull2_ByExpectation`) → per-seq `bias`; final
  `score = (fwd − (null + seqbias))/log2`. `omega = 1/256`.
- **E-value / P-value layer** — **Gumbel** survival for MSV/Viterbi, **exponential** tail for Forward
  (Eddy 2008; λ ≈ log2), `E = P·Z`, read from the profile `STATS LOCAL` lines. Verified to 1e-9 vs a
  hand-derived pin.
- **Multi-domain envelope decomposition** (`p7_domaindef` posterior heuristics rt1=0.25/rt2=0.10/
  rt3=0.20) → GBB1/PF00400 resolves to the **7** β-propeller-blade envelopes exactly (coords match
  hmmsearch; scores to ~1e-3 bits, HMMER uses float32 vs this engine's float64); SH3 → 1 envelope 3..50.
- **Stochastic-traceback clustering** (`region_trace_ensemble` + single-linkage `p7_spensemble_Cluster`,
  fixed-seed Easel LCG seed 42) for **closely-overlapping** tandem domains — reproduces hmmsearch
  envelope **coordinates** exactly on truncated tandem-SH3 constructs; opt-in via
  `FindDomains(seq, clusterOverlapping:true)`.

**Honest residuals:** exact `hmmsearch`-reported E-value pipeline parity is **not** reproduced — the
MSV/bias prefilters (which only gate which sequences reach Forward, not a hit's bit score) are not
reimplemented; **bit-for-bit** trace-ensemble identity would need HMMER's float32 arithmetic; and only
the **3 bundled CC0 profiles** are shipped (any other family is a caller-supplied `.hmm`). A
[[research-grade-limitations|research-grade]] implementation.

## `PredictSignalPeptide` — von Heijne tripartite model + −1,−3 rule

A **signal peptide** is an N-terminal 16–30 aa targeting sequence with a **tripartite** structure
(von Heijne 1985/1986):

- **n-region** (1–5 residues, positively charged K/R),
- **h-region** (7–15 hydrophobic residues, a single α-helix — "necessary and sufficient for membrane
  targeting"),
- **c-region** (3–7 polar residues) ending at a cleavage site obeying the **−1,−3 rule**: small neutral
  residues **{A, G, S}** at positions −1 and −3.

Scoring: `score = (nScore + 2·hScore + cScore)/4` (h-region double weight, per its dominant role).
Detection uses two **literature-derived constraints** (no arbitrary threshold): the n-region must carry
positive charge (`nScore > 0`) and the h-region must be predominantly hydrophobic (`hScore ≥ 0.5`).
`Probability = Score` directly (no scaling). Returns the three regions and the cleavage position.

## Invariants, oracles, corner cases

- Domain hits: 0-based inclusive `Start`/`End`; empty/null input → no domains; tandem-repeat domains
  (multiple zinc fingers, 7 WD40 blades) all detected; only a match fully contained in another is
  suppressed.
- Signal peptide: sequences shorter than ~15 aa cannot hold n+h+c → none; charged/cytoplasmic
  sequences → null (no signal); case-insensitive.
- Synthetic oracles: C2H2 `AAAACXXCXXXLXXXXXXXXHXXXHAAA` → domain at 4..24 (PF00096); P-loop
  `AAAAGXXXXGKSAAAA` → 4..11 (PF00069); signal peptide `MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF` → detected,
  cleavage ≈ 25.

## Design decisions (previously assumptions — all resolved)

1. **Signal-peptide weights 1:2:1** — h-region double weight, justified by von Heijne (1985) "hydrophobic
   core … necessary and sufficient".
2. **Evidence-based constraints** replace the old arbitrary 0.4 threshold (n-region charge + h-region
   hydrophobicity).
3. **`Probability = Score`** — the old `min(1, score×1.2)` scaling was eliminated.
4. **Small-residue set strictly {A, G, S}** — the canonical −1,−3 set (von Heijne 1983); T/C/V/L
   occasionally seen but excluded.
5. **PROSITE-derived regex patterns** for domain detection is a deliberate scope decision (PROSITE itself
   uses consensus patterns; Hulo 2006), not an assumption.
6. **Method naming `FindDomains`** (not `PredictDomains`) — pattern search, not probabilistic prediction.

## References

von Heijne 1983 (Eur J Biochem 133:17–21, −1,−3 rule) & 1985/1986 (J Mol Biol 184:99–105, tripartite
h-core); Walker 1982 (P-loop); Krishna 2003 & Pabo 2001 (C2H2 zinc finger); PROSITE PS00028/PS00017/
PS00678 + ScanProsite syntax; Pfam PF00096/PF00400/PF00018/PF00595/PF00069 (El-Gebali 2019, Mistry
2021, CC0 licence); HMMER User's Guide v3.4 (Eddy 2023), Durbin et al. 1998 §5.4, Eddy 2008/2011,
HMMER/Easel source (`generic_fwdback.c`, `generic_null2.c`, `p7_domaindef.c`, `p7_spensemble.c`,
`esl_gumbel.c`), pyhmmer 0.12.1. Full citations in [[protmotif-domain-001-evidence]] (do not duplicate).
