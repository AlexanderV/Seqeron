# Validation Report: CRISPR-PAM-001 — PAM Site Detection

- **Validated:** 2026-06-12   **Area:** Molecular Tools (MolTools)
- **Canonical method(s):** `CrisprDesigner.FindPamSites(DnaSequence, CrisprSystemType)`, string overload, `CrisprDesigner.GetSystem(CrisprSystemType)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia, *Protospacer adjacent motif* | SpCas9 PAM = 5'-NGG-3' located **immediately 3' (downstream)** of the protospacer; "Cas9 will not bind/cleave unless the target is followed by the PAM". Cas12a/Cpf1 PAM is T-rich (5'-TTTN-3' / 5'-TTTV-3') located **5' (upstream)** of the protospacer. IUPAC: N=any, V=A/C/G (not T), Y=C/T. |
| Wikipedia, *Cas9* | Canonical PAM NGG; NAG and NGA tolerated as alternative motifs **with lower cleavage activity**; 20-nt spacer-complementary region of the gRNA; PAM is on the non-target strand 3 nt downstream of the cleavage site (i.e. 3' of protospacer). |
| Wikipedia, *Cas12a* | Cas12a recognizes T-rich PAM "5'-TTTV-3', where V is A, C, or G"; PAM is **5' upstream** of the target; crRNA shorter than Cas9 sgRNA; *Acidaminococcus* and *Lachnospiraceae* orthologs are active genome editors. |
| Jinek et al. 2012, Science 337:816 | SpCas9 NGG PAM, 20-nt guide segment of the sgRNA. |
| Ran et al. 2015, Nature 520:186 | SaCas9 recognizes **NNGRRT** PAM (R = A or G); ~21–23 nt guide optimal. |
| Zetsche et al. 2015, Cell 163:759 | Cas12a/Cpf1 T-rich (TTTV) PAM 5' of spacer; AsCas12a / LbCas12a orthologs; ~23 nt spacer. |
| Hsu et al. 2013, Nat Biotechnol 31:827 | NAG as secondary SpCas9 PAM with reduced activity. |
| Liu et al. 2019, Nature 566:218 | DpbCasX (Cas12e) has a 20-nt guide and recognizes a **5'-TTCN** PAM. |

### Canonical PAM sequences and positions (validated)

| System | PAM | IUPAC expansion | Position vs protospacer | Guide len | Source |
|--------|-----|-----------------|-------------------------|-----------|--------|
| SpCas9 | NGG | N∈{A,C,G,T} | 3' (after) | 20 | Jinek 2012; Wikipedia Cas9/PAM |
| SpCas9-NAG | NAG | secondary, lower activity | 3' (after) | 20 | Hsu 2013; Wikipedia Cas9 |
| SaCas9 | NNGRRT | R∈{A,G} | 3' (after) | 21 (21–23 range) | Ran 2015 |
| Cas12a/Cpf1 | TTTV | V∈{A,C,G} | 5' (before) | 23 | Zetsche 2015; Wikipedia Cas12a |
| AsCas12a | TTTV | V∈{A,C,G} | 5' (before) | 23 | Zetsche 2015 |
| LbCas12a | TTTV | V∈{A,C,G} | 5' (before) | 24 (23–25 range) | Zetsche 2015 |
| CasX/Cas12e | TTCN | N∈{A,C,G,T} | 5' (before) | 20 | Liu 2019 |

### IUPAC, strand, coordinates

- IUPAC interpretation N=any, R=A/G, V=A/C/G confirmed against NC-IUB 1984 standard and matches `IupacHelper.MatchesIupac`.
- PAM is biologically searched on **both strands** (a guide can target either strand); the implementation searches forward and reverse-complement strands.
- Coordinate base: **0-based** half-open indices into the supplied sequence (implementation/spec convention; the literature is biology, not coordinate-specific). Reported `Position` is the 0-based start of the PAM on the forward strand.

### Hand-computed worked examples

1. **Forward NGG (M6/M7).** `seq = "ACGTACGTACGTACGTACGTAGG"` (len 23), SpCas9. NGG matches at i=20 ("AGG"); `pamAfterTarget`⇒ target = seq[0..19] = "ACGTACGTACGTACGTACGT" (20 bp). Position=20, TargetStart=0. ✓ matches spec & test.
2. **Overlapping (S3).** `seq = "...ACGTAGGTGG"` — NGG at 20 ("AGG") and 23 ("TGG"); both reported, counts=2. ✓
3. **Reverse strand (M8).** `seq = "CCAACGTACGTACGTACGTACGTACGTACGT"` (len 31). No forward GG. revComp ends in "TGG" at revComp index 28 → NGG match. Forward position = 31−28−3 = **0**; PamSequence = revcomp("TGG") = **"CCA"** = forward bases at pos 0; target len 20. ✓ exactly matches the test's `Position==0`, `PamSequence=="CCA"`, `len==20`.
4. **Cas12a TTTV (M13/M14).** `TTTA…` matches (V=A); `TTTT…` does **not** (V excludes T). ✓
5. **SaCas9 NNGRRT (M15).** "AAGAAT","AGGAGT","AAGGGT" match (R∈{A,G}); "AAGCCT" rejected (C∉R). ✓
6. **Boundary (S1).** SpCas9 PAM at pos 0 ⇒ target would start at −20 ⇒ excluded. Cas12a TTTA with only 10 bp ⇒ target would run past end ⇒ excluded. ✓

### Findings / divergences (Stage A)

None biological. All PAM strings, IUPAC expansions, and 3'/5' placements are correct against the cited primary literature and Wikipedia. PASS.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`
- `GetSystem` (lines 18–28): all 7 systems' Name/PAM/GuideLength/PamAfterTarget literals match the validated table exactly; unknown type throws `ArgumentException`.
- `FindPamSites` (DnaSequence) lines 40–46: null-guarded, delegates to core.
- `FindPamSites` (string) lines 51–60: empty⇒empty; upper-cases (case-insensitive).
- `FindPamSitesCore` lines 62–135: forward + reverse-complement scan; per-enzyme upstream/downstream target placement; boundary check `targetStart>=0 && targetEnd<len`.
- `MatchesPam` / `MatchesIupac` lines 137–154: delegates to `IupacHelper.MatchesIupac` (15-code standard, correct N/R/V).

### Formula realised correctly?

Yes. Forward scan iterates `0..len-pamLen`, applies IUPAC PAM match, and for `pamAfterTarget` computes target = `[i-guide, i-1]`, else `[i+pamLen, i+pamLen+guide-1]`. Both-strand search done via `GetReverseComplementString`. Reverse PAM reported back in forward coordinates: `Position = len - i - pamLen` and `PamSequence = revcomp(revComp.Substring(i, pamLen))`. Hand traces (above) reproduce the exact code outputs.

### Cross-verification table recomputed vs code

All six worked examples above were recomputed by hand and match the code's reported values and the passing tests (Position, PamSequence, TargetSequence content/length, exclusion at boundaries, both-strand counts).

### Variant/delegate consistency

- String and DnaSequence overloads produce identical results (verified by `FindPamSites_StringOverload_WorksIdentically`).
- AsCas12a (23 bp) and LbCas12a (24 bp) exercised end-to-end through `FindPamSites` with correct target lengths/content.

### Test quality audit

58–59 tests assert **exact** sourced values: exact PAM strings, exact positions (20, 0, 23), exact target content, exact counts (overlap=2, reverse=1), positive+negative IUPAC cases (TTTV accepts A/C/G rejects T; NNGRRT rejects C at R), and edge cases (empty, null, both out-of-bounds directions). Assertions are deterministic and not tautological. PASS for test realism.

### Findings / defects (Stage B)

No defect against the spec. One documented **note** (not a bug):

- For **reverse-strand** hits, the `PamSite.TargetStart` field holds an index into the *reverse-complement* string (e.g. 8 in the M8 trace), whereas `Position` is a *forward-strand* coordinate. `TargetSequence` is correctly sliced from revComp using that index, so the record is internally consistent, but `TargetStart` and `Position` are in different coordinate systems for reverse hits. Spec invariant 3 and test S2 only require **`Position`** to be converted to forward coordinates (which it is, correctly); no test asserts reverse `TargetStart` as a forward coordinate, and `GetReverseComplement` consumers (guide design / off-target) only read `TargetSequence`. This is a pre-existing design characteristic, fully covered by the spec's narrowed invariant, so no code change is made. Flagged here for downstream consumers that might assume `TargetStart` is always forward-strand.

## Verdict & follow-ups

- **Stage A: PASS** — every PAM sequence, IUPAC expansion, and 3'/5' placement confirmed against primary literature (Jinek 2012, Ran 2015, Zetsche 2015, Hsu 2013, Liu 2019) and Wikipedia.
- **Stage B: PASS-WITH-NOTES** — implementation faithfully realises the validated definitions; all worked examples reproduce; 59 PAM tests + full suite (4461) green. One documented coordinate-system note for reverse-strand `TargetStart` (consistent with the spec's invariant, not a defect).
- **State: CLEAN** — no code change required; project builds and all tests pass.
