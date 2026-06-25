# Validation Report: CRISPR-PAM-001 — PAM Site Detection

- **Validated:** 2026-06-24   **Area:** Molecular Tools (MolTools)
- **Canonical method(s):** `CrisprDesigner.FindPamSites(DnaSequence, CrisprSystemType)`, `FindPamSites(string, …)` overload, `CrisprDesigner.GetSystem(CrisprSystemType)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm

| Source (opened this session) | Confirms |
|--------|----------|
| Wikipedia, *Protospacer adjacent motif* | SpCas9 PAM = **5'-NGG-3'**, located **immediately 3' (downstream)** of the protospacer; the guide targets sequences "having a PAM at the 3'-end". "N" = any of A/T/G/C. PAM is read on the strand complementary to the guide RNA target (i.e. the non-target / protospacer strand). Cas12a/Cpf1 PAM is T-rich (5'-TTTN / 5'-TTTV), located **5' (upstream)** of the protospacer. |
| WebSearch — SaCas9 (Ran 2015; Genome Biology 16:215) | SaCas9 recognizes **5'-NNGRRT-3'** PAM (R = A or G); robust activity at spacer lengths 20–24, highest editing at 21–23 nt. |
| WebSearch — Cas12a (Zetsche 2015) | Cas12a/Cpf1 PAM = **TTTV** (V = A, C, or G), 5' of target. AsCas12a / LbCas12a orthologs. |
| WebSearch — CasX/Cas12e (Liu 2019) | DpbCasX cleaves dsDNA with **5'-TTCN** PAM; 20-nt guide. |
| (carried from prior validation, cb113ce) Jinek 2012, Hsu 2013 | SpCas9 NGG PAM + 20-nt guide; NAG as a secondary SpCas9 PAM with reduced activity. |

### Canonical PAM sequences and positions (validated)

| System | PAM | IUPAC expansion | Position vs protospacer | Guide len | Source |
|--------|-----|-----------------|-------------------------|-----------|--------|
| SpCas9 | NGG | N∈{A,C,G,T} | 3' (after) | 20 | Jinek 2012; Wikipedia PAM |
| SpCas9-NAG | NAG | secondary, lower activity | 3' (after) | 20 | Hsu 2013; Wikipedia Cas9 |
| SaCas9 | NNGRRT | R∈{A,G} | 3' (after) | 21 (21–23 range) | Ran 2015 |
| Cas12a/Cpf1 | TTTV | V∈{A,C,G} | 5' (before) | 23 | Zetsche 2015 |
| AsCas12a | TTTV | V∈{A,C,G} | 5' (before) | 23 | Zetsche 2015 |
| LbCas12a | TTTV | V∈{A,C,G} | 5' (before) | 24 (23–25 range) | Zetsche 2015 |
| CasX/Cas12e | TTCN | N∈{A,C,G,T} | 5' (before) | 20 | Liu 2019 |

### IUPAC, strand, coordinates

- IUPAC interpretation N=any, R=A/G, V=A/C/G confirmed against the NC-IUB 1984 standard, and matches `IupacHelper.MatchesIupac` (15-code table, N/R/V verified by reading the source).
- A guide can target either strand, so PAM is searched on **both** the forward strand and its reverse complement. A forward-strand **CCN** is the reverse-complement of a reverse-strand **NGG** — both strands must be scanned.
- Coordinate base: **0-based** index into the supplied sequence; reported `Position` is the 0-based start of the PAM, **always on the forward strand**.

### Independent hand cross-check (designed sequence with NGG on forward + CCN on reverse)

`seq = "CCAACGTACGTACGTACGTACGTACGTACGT"` (len 31), SpCas9.
- **Forward** scan: no `NGG` (no `GG` dinucleotide on the forward strand) → 0 forward hits.
- **Reverse**: `revComp = "ACGTACGTACGTACGTACGTACGTACGTTGG"`; the only `NGG` is `TGG` at revComp index 28. Converting back to forward: `Position = len − i − pamLen = 31 − 28 − 3 = 0`; `PamSequence = revcomp("TGG") = "CCA"`, which equals the forward bases at position 0 (the `CCN` reverse-strand PAM). Target length 20.

Recomputed independently in Python this session — matches both biology and the M8 test (`Count==1`, `Position==0`, `PamSequence=="CCA"`, target len 20). Additional hand traces reproduced: forward NGG @20 with 20-nt target @0 (M6/M7); Cas12a TTTV accepts A/C/G, rejects TTTT (M13/M14); SaCas9 NNGRRT accepts R∈{A,G}, rejects C at R (M15); overlap `AGGTGG` → two sites @20,@23 (S3); boundary exclusion both directions (S1).

### Findings / divergences (Stage A)

None. Every PAM string, IUPAC expansion, and 3'/5' placement is correct against the cited literature. **PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`
- `GetSystem` (L18–28): all 7 systems' Name/PAM/GuideLength/PamAfterTarget literals match the validated table exactly; unknown type → `ArgumentException`.
- `FindPamSites(DnaSequence)` (L40–46): null-guarded (`ArgumentNullException.ThrowIfNull`), delegates to core.
- `FindPamSites(string)` (L51–60): null/empty ⇒ empty; `ToUpperInvariant` (case-insensitive).
- `FindPamSitesCore` (L62–135): forward + reverse-complement scan; per-enzyme target placement (`pamAfterTarget` ⇒ target `[i−guide, i−1]`, else `[i+pamLen, i+pamLen+guide−1]`); boundary check `targetStart ≥ 0 && targetEnd < len`; reverse hits report `Position = len − i − pamLen` (forward coord) and `PamSequence = revcomp(...)`.
- `MatchesPam` / `MatchesIupac` (L137–154): delegates to `IupacHelper.MatchesIupac`; N/R/V verified correct.

### Formula realised correctly?

Yes. The forward/reverse two-pass scan with the per-enzyme upstream/downstream placement matches the validated description. The reverse-strand coordinate conversion (`len − i − pamLen`) and the reverse-complemented `PamSequence` were reproduced exactly by the independent Python trace above.

### Cross-verification table recomputed vs code

All hand traces (forward NGG @20; reverse TGG→fwd 0 / PAM "CCA"; Cas12a TTTV±; SaCas9 NNGRRT±; overlap @20,@23; both-direction boundary exclusion) reproduce the code's outputs and the passing assertions.

### Variant/delegate consistency

- String and `DnaSequence` overloads return identical sites (`FindPamSites_StringOverload_WorksIdentically`).
- AsCas12a (23 nt) and LbCas12a (24 nt) exercised end-to-end with correct target length/content.

### Test quality audit

58 tests assert **exact** sourced values — exact PAM strings, exact positions (0/20/23), exact target content/length, exact counts (overlap=2, reverse=1), positive **and** negative IUPAC cases (TTTV accepts A/C/G, rejects T; NNGRRT rejects C at R), and edge cases (empty, null, both-direction out-of-bounds). Deterministic, non-tautological. **58/58 pass** this session.

### Findings / defects (Stage B)

No defect against the validated description. One documented **note** (not a bug):

- For **reverse-strand** hits, `PamSite.TargetStart` is an index into the *reverse-complement* string (used to slice `TargetSequence`), whereas `Position` is a *forward-strand* coordinate — the two fields are in different coordinate systems for reverse hits. This is now **explicitly documented in the `PamSite.TargetStart` XML doc** (CrisprDesigner.cs L1035–1041) and the spec's invariant 3 only requires `Position` to be forward-strand (which it is). `TargetSequence` is sliced correctly from revComp, so the record is internally consistent. Flagged for downstream consumers that might assume `TargetStart` is always forward-strand.

## Verdict & follow-ups

- **Stage A: PASS** — every PAM sequence, IUPAC expansion, and 3'/5' placement confirmed against primary literature (Jinek 2012, Ran 2015, Zetsche 2015, Hsu 2013, Liu 2019) and Wikipedia.
- **Stage B: PASS-WITH-NOTES** — implementation faithfully realises the validated definitions; the both-strand hand cross-check (NGG forward / CCN reverse) reproduces exactly; 58/58 PAM tests green. One documented coordinate-system note for reverse-strand `TargetStart` (consistent with the spec invariant and the field's own XML doc — not a defect).
- **State: CLEAN** — no code change required; project builds and the unit's tests pass.
