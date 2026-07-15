---
type: concept
title: "Protein domain prediction (PROSITE pattern + Plan7 profile-HMM) & signal-peptide prediction"
tags: [analysis, algorithm, protein, motif]
mcp_tools:
  - find_protein_domains
  - predict_signal_peptide
sources:
  - docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
  - docs/Evidence/PROTMOTIF-SP-001-Evidence.md
source_commit: 1513221ea012107594feec09dd5b4e850e3d5f37
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
      object: concept:test-unit-registry
      source: protmotif-sp-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-SP-001 ... Algorithm: Signal Peptide Cleavage-Site Prediction (von Heijne 1986 weight matrix) on ProteinMotifFinder.PredictSignalPeptide"
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

Two related protein-feature capabilities on `ProteinMotifFinder`: **domain prediction** (`FindDomains`,
plus an opt-in profile-HMM engine `FindDomainsByHmm`), validated by test unit **PROTMOTIF-DOMAIN-001**
([[protmotif-domain-001-evidence]]); and **signal-peptide prediction** (`PredictSignalPeptide`),
validated by test unit **PROTMOTIF-SP-001** ([[protmotif-sp-001-evidence]]) against the von Heijne 1986
weight matrix. These are units of the **ProteinMotif family**, siblings of the fixed-dictionary
[[common-protein-motifs]] and the windowed [[coiled-coil-prediction]]. [[test-unit-registry]] tracks the
units; see [[algorithm-validation-evidence]] for the artifact pattern.

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

## `PredictSignalPeptide` — von Heijne 1986 weight matrix (EMBOSS `sigcleave`)

The **current** `PredictSignalPeptide` implements the **von Heijne (1986) log-odds weight-matrix**
cleavage-site method — the algorithm of EMBOSS `sigcleave`. Test unit **PROTMOTIF-SP-001**
([[protmotif-sp-001-evidence]]) validates it. This **superseded** an earlier tripartite n/h/c + −1,−3
scoring (constants 0.95/0.825, an `NRegion`/`Probability` record) that was found **fabricated** and
removed; see the "Superseded" note below.

A **signal peptide** is an N-terminal 16–30 aa targeting sequence. Prediction is by a
position-specific **weight matrix** rather than the tripartite heuristic:

- At each candidate site the **score** is `Σ ln(count/expect)` over positions **−13..+2** (15 matrix
  columns; the `Expect` column is not scored). Natural log. Zero counts become `1.0e-10` at the
  conserved columns **−3 and −1** (a strong penalty) and `1.0` elsewhere before the log.
- The single prediction is the **argmax** of the score over all positions; cleavage lies between `−1`
  and `+1`, and `+1` (1-based `CleavagePosition`) is the first residue of the mature protein.
- **`IsLikelySignalPeptide ⇔ Score ≥ 3.5`** (`minWeight` default 3.5; ~95% sensitivity/specificity,
  cleavage-site correct in only 75–80% of cases — heuristic).
- Two 20×16 count matrices from von Heijne (1986): **eukaryotic** (161 aligned sequences, default) and
  **prokaryotic** (36, via `prokaryote:true`), each column summing to its sample size.

The record returns `CleavagePosition`, `Score`, `SignalSequence` (residues 1..`CleavagePosition−1`),
`WindowSequence` (the 15-residue window), and `IsLikelySignalPeptide`. Inputs shorter than one full
15-residue window return `null`; non-standard residues (X/B/Z/*) contribute 0; case-insensitive.

**Worked oracle:** **ACH2_DROME (UniProt P17644**, 576 aa) → maximum `Score` **13.739**,
`CleavagePosition` **42**, window `LLVLLLLCETVQA` (residues 29–41); reproduced exactly by an
independent Python re-derivation.

> **Superseded (2026-06-14):** DOMAIN-001 ([[protmotif-domain-001-evidence]]) described this same
> method as a **tripartite n/h/c + −1,−3-rule** model with `score = (nScore + 2·hScore + cScore)/4`
> and `Probability = Score`. That model's fields and constants were fabricated and have been replaced
> by the weight-matrix method above. Treat the tripartite description as **historical** for
> `PredictSignalPeptide`.

## Invariants, oracles, corner cases

- Domain hits: 0-based inclusive `Start`/`End`; empty/null input → no domains; tandem-repeat domains
  (multiple zinc fingers, 7 WD40 blades) all detected; only a match fully contained in another is
  suppressed.
- Signal peptide: `Score` = argmax over sites of `Σ ln(count/expect)` (−13..+2); `CleavagePosition`
  1-based mature start; `IsLikelySignalPeptide ⇔ Score ≥ 3.5`; inputs shorter than the 15-residue
  window return `null`; non-standard residues contribute 0; case-insensitive.
- Synthetic oracles: C2H2 `AAAACXXCXXXLXXXXXXXXHXXXHAAA` → domain at 4..24 (PF00096); P-loop
  `AAAAGXXXXGKSAAAA` → 4..11 (PF00069); signal-peptide ACH2_DROME (P17644) → `Score` 13.739,
  `CleavagePosition` 42, window `LLVLLLLCETVQA`.

## Design decisions (previously assumptions — all resolved)

1. **Signal-peptide method = von Heijne 1986 weight matrix** (EMBOSS `sigcleave`), replacing the
   earlier fabricated tripartite n/h/c heuristic (weights 1:2:1, `Probability = Score`, {A,G,S} −1,−3
   set) — see the Superseded note and [[protmotif-sp-001-evidence]].
2. **Natural-log log-odds + −3/−1 zero-count `1.0e-10` penalty** — the exact EMBOSS transform;
   base-10 or a missing penalty changes every score.
3. **`minWeight` default 3.5** governs `IsLikelySignalPeptide` only; the argmax site is always returned
   (no intrinsic cutoff).
4. **Eukaryotic matrix default; prokaryotic via `prokaryote:true`** — both von Heijne (1986).
5. **PROSITE-derived regex patterns** for domain detection is a deliberate scope decision (PROSITE itself
   uses consensus patterns; Hulo 2006), not an assumption.
6. **Method naming `FindDomains`** (not `PredictDomains`) — pattern search, not probabilistic prediction.

## References

von Heijne 1983 (Eur J Biochem 133:17–21, −1,−3 rule) & 1985/1986 (J Mol Biol 184:99–105, tripartite
h-core); Walker 1982 (P-loop); Krishna 2003 & Pabo 2001 (C2H2 zinc finger); PROSITE PS00028/PS00017/
PS00678 + ScanProsite syntax; Pfam PF00096/PF00400/PF00018/PF00595/PF00069 (El-Gebali 2019, Mistry
2021, CC0 licence); HMMER User's Guide v3.4 (Eddy 2023), Durbin et al. 1998 §5.4, Eddy 2008/2011,
HMMER/Easel source (`generic_fwdback.c`, `generic_null2.c`, `p7_domaindef.c`, `p7_spensemble.c`,
`esl_gumbel.c`), pyhmmer 0.12.1. Signal-peptide method: von Heijne 1986 (Nucleic Acids Res 14:4683–4690)
via EMBOSS 6.6.0 `sigcleave` (`sigcleave.c`, `data/Esig.euk`/`Esig.pro`; Rice 2000) + UniProt P17644.
Full citations in [[protmotif-domain-001-evidence]] and [[protmotif-sp-001-evidence]] (do not duplicate).
