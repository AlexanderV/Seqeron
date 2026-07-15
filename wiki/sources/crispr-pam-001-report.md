---
type: source
title: "Validation report: CRISPR-PAM-001 (PAM site detection — locating protospacer-adjacent motifs on both strands for gRNA design, CrisprDesigner.FindPamSites / GetSystem)"
tags: [validation, primer, governance]
doc_path: docs/Validation/reports/CRISPR-PAM-001.md
sources:
  - docs/Validation/reports/CRISPR-PAM-001.md
source_commit: 13507add3b861dc2ba03395b72f3913eda8ef054
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CRISPR-PAM-001

The two-stage **validation write-up** for test unit **CRISPR-PAM-001** (CRISPR **PAM site detection** —
locating protospacer-adjacent motifs such as SpCas9 `NGG` on **both** strands, the upstream step that
feeds candidate-guide extraction), area **MolTools**, validated 2026-06-24. This is the *report* artifact
that feeds one row of the [[validation-ledger]]; it records the validator's independent **verdict** on
both the algorithm description (Stage A) and the shipped code (Stage B), inside the wider
[[validation-and-testing]] campaign. PAM finding is the geometric front end of guide design: it resolves
each CRISPR system's motif/orientation and returns the PAM sites at which the composition heuristic then
extracts spacers — synthesized as **Layer 0** of the concept [[crispr-guide-rna-design]] (the
MolTools / `CrisprDesigner` anchor, sibling to the primer/probe units). It is the PAM-geometry sibling of
the on-target report [[crispr-guide-001-report]] and the off-target report [[crispr-off-001-report]], all
three on the same `CrisprDesigner` class. [[test-unit-registry]] defines the unit. Source under test:
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End state: ✅ CLEAN.** No defect; no production-code change
required. Every PAM string, IUPAC expansion, and 3′/5′ placement is correct against the cited literature;
the both-strand hand cross-check reproduces exactly; **58/58** PAM tests green (build 0 warnings /
0 errors). The single Stage-B "note" is a documented coordinate-system caveat on reverse-strand hits
(consistent with the field's own XML doc and the spec invariant) — **not** a bug.

## Canonical methods & source under test

In `CrisprDesigner.cs`:

- `FindPamSites(DnaSequence, CrisprSystemType)` (`:40–46`) — null-guarded
  (`ArgumentNullException.ThrowIfNull`), delegates to core.
- `FindPamSites(string, …)` overload (`:51–60`) — null/empty ⇒ empty; `ToUpperInvariant`
  (case-insensitive). Returns **identical** sites to the `DnaSequence` overload.
- `GetSystem(CrisprSystemType)` (`:18–28`) — resolves each of the 7 systems' Name / PAM / GuideLength /
  `PamAfterTarget` literals; unknown type → `ArgumentException`.
- `FindPamSitesCore` (`:62–135`) — the forward + reverse-complement two-pass scan.
- `MatchesPam` / `MatchesIupac` (`:137–154`) — delegate to `IupacHelper.MatchesIupac` (15-code table).

## Stage A — description (biology faithfulness)

Grounded against primary literature opened this session (Jinek 2012, Hsu 2013, Ran 2015 *Genome Biol*
16:215, Zetsche 2015, Liu 2019 CasX/Cas12e) plus Wikipedia *Protospacer adjacent motif*. All PAM
sequences, IUPAC expansions, and 3′/5′ placements confirmed — **findings: none, PASS.**

### Canonical PAM table (validated, 7 systems)

| System | PAM | IUPAC | Position vs protospacer | Guide len | Source |
|--------|-----|-------|-------------------------|-----------|--------|
| SpCas9 | `NGG` | N∈{A,C,G,T} | 3′ (after) | 20 | Jinek 2012; Wikipedia |
| SpCas9-NAG | `NAG` | secondary, lower activity | 3′ (after) | 20 | Hsu 2013 |
| SaCas9 | `NNGRRT` | R∈{A,G} | 3′ (after) | 21 (21–23) | Ran 2015 |
| Cas12a/Cpf1 | `TTTV` | V∈{A,C,G} | 5′ (before) | 23 | Zetsche 2015 |
| AsCas12a | `TTTV` | V∈{A,C,G} | 5′ (before) | 23 | Zetsche 2015 |
| LbCas12a | `TTTV` | V∈{A,C,G} | 5′ (before) | 24 (23–25) | Zetsche 2015 |
| CasX/Cas12e | `TTCN` | N∈{A,C,G,T} | 5′ (before) | 20 | Liu 2019 |

The `GetSystem` literals match this table exactly.

### IUPAC, strand, and coordinates

- **IUPAC:** N=any, R=A/G, V=A/C/G confirmed against the NC-IUB 1984 standard and against
  `IupacHelper.MatchesIupac` (N/R/V verified by reading the source).
- **Both strands scanned.** A guide can target either strand, so the PAM is searched on the **forward**
  strand *and* its **reverse complement**. A forward-strand `CCN` is the reverse-complement of a
  reverse-strand `NGG` — both strands must be scanned to find all sites.
- **Coordinates: 0-based** index into the supplied sequence; reported `Position` is the **0-based start of
  the PAM, always on the forward strand**. For reverse-strand hits `Position = len − i − pamLen` (forward
  coordinate) and `PamSequence = revcomp(...)`.

### Independent hand cross-check (both-strand)

`seq = "CCAACGTACGTACGTACGTACGTACGTACGT"` (len 31), SpCas9:
- **Forward** scan → no `NGG` (no `GG` dinucleotide forward) ⇒ **0** forward hits.
- **Reverse:** `revComp = "ACGTACGTACGTACGTACGTACGTACGTTGG"`; only `NGG` is `TGG` at revComp index 28.
  Back to forward: `Position = 31 − 28 − 3 = 0`; `PamSequence = revcomp("TGG") = "CCA"` (= the `CCN`
  reverse-strand PAM at forward position 0); target length 20.

Recomputed independently in Python this session — matches biology and the **M8** test (`Count==1`,
`Position==0`, `PamSequence=="CCA"`, target len 20). Further hand traces reproduced: forward `NGG` @20 with
20-nt target @0 (M6/M7); Cas12a `TTTV` accepts A/C/G rejects `TTTT` (M13/M14); SaCas9 `NNGRRT` accepts
R∈{A,G} rejects C at R (M15); overlap `AGGTGG` → two sites @20,@23 (S3); boundary exclusion both
directions (S1).

## Stage B — implementation

- `FindPamSitesCore` (`:62–135`) does the forward + reverse-complement scan with **per-enzyme placement**:
  `PamAfterTarget` ⇒ target is `[i−guide, i−1]`; else target is `[i+pamLen, i+pamLen+guide−1]`. A
  **boundary check** `targetStart ≥ 0 && targetEnd < len` drops sites whose spacer would run off either
  end. Reverse hits report the forward-coordinate `Position = len − i − pamLen` and the reverse-complemented
  `PamSequence`.
- The reverse-strand coordinate conversion (`len − i − pamLen`) and reverse-complemented `PamSequence` were
  reproduced exactly by the independent Python trace above. Formula realised correctly.
- **Variant/delegate consistency:** the `string` and `DnaSequence` overloads return identical sites
  (`FindPamSites_StringOverload_WorksIdentically`); AsCas12a (23 nt) and LbCas12a (24 nt) exercised
  end-to-end with correct target length/content.
- **Test-quality audit (HARD gate): PASS, no green-washing.** **58 tests** assert **exact** sourced values —
  exact PAM strings, exact positions (0/20/23), exact target content/length, exact counts (overlap=2,
  reverse=1), **positive and negative** IUPAC cases (`TTTV` accepts A/C/G rejects T; `NNGRRT` rejects C at
  R), and edge cases (empty, null, both-direction out-of-bounds). Deterministic, non-tautological.
  **58/58 pass.**

### The one documented note (not a defect)

For **reverse-strand** hits, `PamSite.TargetStart` is an index into the *reverse-complement* string (it is
used to slice `TargetSequence`), whereas `Position` is a *forward-strand* coordinate — the two fields are in
**different coordinate systems** for reverse hits. This is now explicitly documented in the
`PamSite.TargetStart` XML doc (`CrisprDesigner.cs :1035–1041`), and the spec's invariant 3 only requires
`Position` to be forward-strand (which it is). `TargetSequence` is sliced correctly from the reverse
complement, so the record is internally consistent. Flagged for downstream consumers that might assume
`TargetStart` is always forward-strand — this is the same coordinate caveat surfaced by the guide-design
Layer-1 note in [[crispr-guide-rna-design]] (`Position` copied from `pamSite.TargetStart` for reverse
designs).

## Findings & follow-ups

- **No code defect, no production-code change (State CLEAN).** Every PAM string, IUPAC expansion, and 3′/5′
  placement matches the cited literature; the both-strand hand cross-check (`NGG` forward / `CCN` reverse)
  reproduces exactly; 58/58 PAM tests green.
- One documented **coordinate-system note** for reverse-strand `PamSite.TargetStart` (consistent with the
  spec invariant and the field's own XML doc) — a downstream-consumer caveat, not a bug.
- Not-yet-ingested MolTools algorithm doc `docs/algorithms/MolTools/PAM_Site_Detection.md` (backlog slug
  `pam-site-detection`) describes the same `FindPamSites` surface at the algorithm level; left **pending** —
  this report validates the unit, not that algorithm doc. Research-grade, not for clinical use.
</content>
</invoke>
