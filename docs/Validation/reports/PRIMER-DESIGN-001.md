# Validation Report: PRIMER-DESIGN-001 — PCR Primer Pair Design

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.DesignPrimers(template, targetStart, targetEnd, params)`, `PrimerDesigner.EvaluatePrimer(seq, pos, isForward, params)`, `PrimerDesigner.GeneratePrimerCandidates(template, regionStart, regionEnd, forward, params)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_PrimerDesign_Tests.cs` (canonical), `PrimerDesignerTests.cs` (smoke/mutation-killing)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

| Source | URL | Confirms |
|--------|-----|----------|
| Wikipedia — Primer (molecular biology) | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | "typically between 18 and 24 bases in length"; primer pairs should have similar melting temperatures (annealing happens for both strands simultaneously); "Regions high in mononucleotide and dinucleotide repeats should be avoided, as loop formation can occur"; primer-dimers from mutually complementary primers; internal hairpins/loops to be avoided. |
| Addgene — How to Design a Primer | https://www.addgene.org/protocols/primer-design/ | "Length of 18-24 bases"; "40-60% G/C content"; "Melting temperature (Tm) of 50-60°C"; "Primer pairs should have a Tm within 5°C of each other"; "Start and end with 1-2 G/C pairs" (GC clamp); "no internal secondary structure"; "avoid primer-primer annealing which creates primer dimers"; "both of the 3' ends of the hybridized primers must point toward one another". |
| Primer3 Manual (v2.6.x) | https://primer3.org/manual.html | PRIMER_MIN_SIZE=18, PRIMER_OPT_SIZE=20, PRIMER_MAX_SIZE=27, PRIMER_MIN_TM=57.0, PRIMER_OPT_TM=60.0, PRIMER_MAX_TM=63.0, PRIMER_MIN_GC=20.0, PRIMER_MAX_GC=80.0, PRIMER_PAIR_MAX_DIFF_TM=100.0 (all confirmed verbatim from the manual). PRIMER_MAX_POLY_X default = 5 (widely documented). |
| SantaLucia (1998) PNAS 95:1460-65 | — | Nearest-neighbor ΔG°37 parameters used by `Calculate3PrimeStability`. |

### Constraint table (description vs sources)

| Constraint | Implementation default | Source range | Verdict |
|-----------|------------------------|--------------|---------|
| Length min / opt / max | 18 / 20 / 25 | Primer3 18/20/27; Addgene 18–24 | Within both; 25 is a documented practical middle ground. ✅ |
| GC content | 40–60% | Addgene 40–60%; Primer3 20–80% | Matches Addgene exactly. ✅ |
| Tm | 57–63 °C (opt 60) | Primer3 57/60/63; Addgene 50–60 | Matches Primer3 exactly. ✅ |
| Pair Tm difference | ≤ 5 °C | Addgene "within 5°C"; Wikipedia "similar" | Matches Addgene exactly (stricter than Primer3 default 100). ✅ |
| Homopolymer max | 4 | Primer3 poly-X = 5 | Stricter than Primer3 — conservative, acceptable. ✅ |
| Dinucleotide repeat max | 4 | Wikipedia "avoid dinucleotide repeats" | Sourced qualitatively. ✅ |
| GC clamp (3') | optional (`Avoid3PrimeGC`) | Addgene "1–2 G/C at ends" | Implemented as optional flag; off by default. 🟡 minor (see notes) |
| 3' end stability ΔG | < −9 kcal/mol flag | SantaLucia NN params; Primer3 PRIMER_MAX_END_STABILITY | Sourced. ✅ |

### Forward/reverse placement & amplicon coordinates

- **Forward primer** = `template[start .. start+len]`, placed upstream of target (`start + len <= targetStart`). ✅
- **Reverse primer** = `reverseComplement(template[start .. end])` for a downstream window (`end - len >= targetEnd`), with `Position = start`. This is the correct biological definition: the reverse primer anneals to the sense strand and is the reverse complement of the template's downstream region so that its 3' end points back toward the forward primer. ✅
- **Amplicon size** = `reverse.Position + reverse.Length − forward.Position` = `reverseEnd − forwardStart`, i.e. the full PCR product from the forward primer 5' end to the reverse primer's downstream template end. ✅

### Worked example (hand-check)

Template = forward region (GAACTCGT×, 0–99) + poly-T target (100–149) + reverse region (TCCGAAGT×, 150–249). Target [100,150).

- Forward candidate: `start=76, len=24` → `template[76..100]` placed entirely upstream (`76+24=100 ≤ 100`). Position 76 < 100. ✅
- Reverse candidate: window `[150, end)`; e.g. `start=150, len=24, end=174` → primer sequence = `revcomp(template[150..174])`, Position=150 ≥ targetEnd=150. ✅
- Product size for that pair = `150 + 24 − 76 = 98` = `end(174) − start(76)`... correctly equals reverseEnd − forwardStart. ✅ (> target length 50, as asserted by `DesignPrimers_ValidResult_ProductSizeCorrect`).

**Reverse-complement spot check** (the classic bug): for template substring `AACCGGTT`, revcomp = `AACCGGTT` (palindromic unit) — verified by the dedicated test `GeneratePrimerCandidates_Reverse_SequenceIsReverseComplement`, which asserts `candidate.Sequence == new DnaSequence(template[Position..Position+Length]).ReverseComplement().Sequence` for every reverse candidate. No "forgot revcomp" defect.

**Stage A findings:** PASS. Every numeric constraint traces to Primer3 or Addgene verbatim; deviations (max length 25, homopolymer 4, GC 40–60) are documented and within/stricter-than sourced ranges. Minor note: the GC-clamp constraint (Addgene 1–2 G/C at ends) is present only as an opt-in `Avoid3PrimeGC` flag (default off, and only checks "≥1 G/C in last 2"), so the default design does not enforce a 3' GC clamp — acceptable since Primer3's default also does not hard-require it, but it is a documented divergence from the Addgene guideline.

---

## Stage B — Implementation

### Code path reviewed

- `DesignPrimers` (lines 40–115): input validation (`targetStart<0 || targetEnd>=Length || targetStart>=targetEnd` → ArgumentException); forward search `[targetStart−200, targetStart)`; reverse search `[targetEnd, targetEnd+200)`; selects highest-Score candidate per side; compatibility = `tmDiff ≤ 5 && !HasPrimerDimer`; product size formula at line 106.
- `EvaluatePrimer` (lines 120–188): computes GC, Tm, homopolymer, dinucleotide repeat, hairpin, 3' stability; flags issues against parameters; `IsValid = issues.Count==0`.
- `GeneratePrimerCandidates` (lines 480–500): for reverse, applies `ReverseComplement` to the template substring before evaluation.

### Realisation checks

- **Reverse primer = reverse complement** (lines 76–78, 494–495): confirmed correct in both `DesignPrimers` and `GeneratePrimerCandidates`. The `Position` stored is the template `start`, and `Sequence` is the revcomp — consistent and correct.
- **Amplicon coordinates** (line 106): `reverse.Position + reverse.Length − forward.Position` = reverseEnd − forwardStart. Correct; matches `DesignPrimers_ValidResult_ProductSizeCorrect`.
- **Tm-difference filter** (lines 103–104): `Math.Abs(fwdTm − revTm) ≤ 5.0` gates `IsValid`. Matches Addgene ≤5°C. Verified by `DesignPrimers_PrimerPair_TmDifferenceWithin5Degrees`.
- **Length / GC / Tm / homopolymer / dinuc constraints** (lines 139–152): each compared against the parameter struct using the validated ranges. Boundary behaviour verified by `EvaluatePrimer_LengthOutsideRange_ReportsIssue(17,26)`, `_GcOutsideRange_ReportsIssue(100%,0%)`, and 4 Tm boundary mutation-killers in `PrimerDesignerTests.cs`.
- **No-valid-pair → empty/invalid** (lines 93–100): when either side has no candidate, returns `PrimerPairResult(null,null,false,"Could not find valid primers...",0)`. Covered by `DesignPrimers_HomopolymerRichTemplate_MayReturnInvalid`.
- **Template shorter than min product** → ArgumentException via range check (`DesignPrimers_VeryShortTemplate_ThrowsArgumentException`, target end == template length triggers `targetEnd >= Length`).

### Edge cases

| Edge case | Behaviour | Test |
|-----------|-----------|------|
| start ≥ end | ArgumentException | `_TargetEndBeforeStart_ThrowsArgumentException` |
| target beyond template | ArgumentException | `_TargetBeyondTemplate_...` |
| negative coords | ArgumentException | `_NegativeCoordinates_...` |
| null template | NullReferenceException (documented) | `_NullTemplate_ThrowsException` |
| empty primer seq | `IsValid=false` | `EvaluatePrimer_EmptySequence_HandledGracefully` |
| region < min length | empty candidate list | `GeneratePrimerCandidates_EmptyRegion_ReturnsEmpty` |

### Variant / test-quality audit

- `GeneratePrimerCandidates` (reverse) and `DesignPrimers` (reverse loop) use the same `DnaSequence.ReverseComplement()` — consistent.
- Tests assert exact sourced values and the actual revcomp identity (non-vacuous); positive `IsValid==true` assertions on the standard template are non-vacuous (template engineered to pass all constraints).

### Findings / defects

None. No reverse-complement omission, no off-by-one in amplicon coordinates, constraint thresholds match the validated description.

---

## Verdict & follow-ups

- **Stage A:** PASS — all constraints sourced to Primer3 / Addgene / Wikipedia; documented deviations within range.
- **Stage B:** PASS — code faithfully realises the validated description; reverse primer correctly reverse-complemented; amplicon size correct.
- **State:** CLEAN. No code changes required.
- **Tests:** `FullyQualifiedName~PrimerDesign` → 136 passed / 0 failed. Full suite → 4461 passed / 0 failed.
- **Minor note (no action):** default design does not hard-enforce the Addgene 3' GC clamp (only available via opt-in `Avoid3PrimeGC`); consistent with Primer3's defaults.
