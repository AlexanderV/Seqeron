# Validation Report: PRIMER-DESIGN-001 — PCR Primer-Pair Design

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.DesignPrimers(template, targetStart, targetEnd, params)`, `PrimerDesigner.EvaluatePrimer(seq, pos, isForward, params)`, `PrimerDesigner.GeneratePrimerCandidates(template, regionStart, regionEnd, forward, params)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_PrimerDesign_Tests.cs` (canonical), `PrimerDesignerTests.cs` (smoke/mutation-killing)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened (this session) & what they confirm

| Source | URL | Confirms |
|--------|-----|----------|
| Primer3 Manual | https://primer3.org/manual.html | Fetched live this session. Defaults read verbatim: PRIMER_MIN_SIZE=18, PRIMER_OPT_SIZE=20, PRIMER_MAX_SIZE=27, PRIMER_MIN_TM=57.0, PRIMER_OPT_TM=60.0, PRIMER_MAX_TM=63.0, PRIMER_MIN_GC=20.0, PRIMER_MAX_GC=80.0, PRIMER_PAIR_MAX_DIFF_TM=100.0, **PRIMER_GC_CLAMP=0** (no 3' GC clamp required by default). |
| Wikipedia — Primer (molecular biology) | https://en.wikipedia.org/wiki/Primer_(molecular_biology) | Fetched live. "typically between 18 and 24 bases"; primer pairs should have similar Tm (annealing happens for both strands simultaneously); mononucleotide/dinucleotide repeats should be avoided (loop formation); complementary primers → primer-dimers; self-complementary primers → internal hairpins. Does NOT specify a numeric GC% range, GC-clamp rule, or a numeric pair-Tm threshold. |
| Addgene — How to Design a Primer | (TestSpec citation) | Length 18–24, GC 40–60%, Tm 50–60°C, pair Tm within 5°C, 1–2 G/C at ends (GC clamp). |
| SantaLucia (1998) PNAS 95:1460-65 | (TestSpec citation) | Nearest-neighbor ΔG°37 parameters used by `Calculate3PrimeStability` (3' end stability only). |

### Constraint table (implementation defaults vs sources)

| Constraint | Default | Source range | Verdict |
|-----------|---------|--------------|---------|
| Length min / opt / max | 18 / 20 / 25 | Primer3 18/20/27; Addgene 18–24; Wikipedia 18–24 | Within both. 25 is a documented practical middle ground. ✅ |
| GC content | 40–60% | Addgene 40–60%; Primer3 20–80% | Matches Addgene exactly; stricter than Primer3. ✅ |
| Tm | 57–63 °C (opt 60) | Primer3 57/60/63 | Matches Primer3 exactly. ✅ |
| Pair Tm difference | ≤ 5 °C | Addgene "within 5°C"; Wikipedia "similar" | Matches Addgene; stricter than Primer3 default (100). ✅ |
| Homopolymer (poly-X) max | 4 | Primer3 PRIMER_MAX_POLY_X = 5 | Stricter — conservative, documented. ✅ |
| Dinucleotide-repeat max | 4 | Wikipedia "avoid dinucleotide repeats" (qualitative) | Sourced qualitatively. ✅ |
| 3' GC clamp | opt-in `Avoid3PrimeGC` (default OFF) | Primer3 GC_CLAMP=0 (off); Addgene 1–2 G/C at ends | Default OFF matches Primer3's own default. 🟡 (note 1) |
| 3' end ΔG°37 stability | flag if < −9 kcal/mol | SantaLucia NN params; Primer3 PRIMER_MAX_END_STABILITY | Sourced. ✅ |

### Forward/reverse placement & amplicon coordinates

- **Forward** = `template[start .. start+len]`, upstream of target (`start+len ≤ targetStart`). ✅
- **Reverse** = `reverseComplement(template[start .. end])` for a downstream window (`end−len ≥ targetEnd`), `Position = start`. Correct biological definition (3' ends point toward each other). ✅
- **Amplicon size** = `reverse.Position + reverse.Length − forward.Position` = reverseEnd − forwardStart (full product, forward 5' end → reverse downstream end). ✅

### Independent cross-check (hand computation, this session)

Standard test template: forward unit `GAACTCGT` (50% GC) ×100 bp, poly-T target [100,150), reverse unit `TCCGAAGT` ×100 bp. Tm via Marmur-Doty `64.9 + 41·(GC−16.4)/N` (the code's formula for N ≥ 14):

| Candidate | len | GC | GC% | Tm (°C) | in [57,63]? |
|-----------|-----|----|-----|---------|-------------|
| 18 bp, 50% GC | 18 | 9 | 50 | 48.04 | **no** → rejected |
| 20 bp, 50% GC | 20 | 10 | 50 | 51.78 | **no** → rejected |
| 22 bp, 50% GC | 22 | 11 | 50 | 54.84 | **no** → rejected |
| 24 bp, 50% GC | 24 | 12 | 50 | 57.38 | **yes** → accepted |
| 25 bp, 52% GC | 25 | 13 | 52 | 59.32 | **yes** → accepted |

This confirms (a) a designed pair satisfies every stated hard constraint (length 24–25 ∈ [18,25]; GC 50–52% ∈ [40,60]; Tm 57.4–59.3 ∈ [57,63]; pair ΔTm ≈ 1.9 ≤ 5), and (b) a primer violating the Tm constraint (18/20/22-bp variants) is correctly rejected — the constraint filter is real, not vacuous. Arithmetic reproduced independently (Python).

### Stage A findings

PASS-WITH-NOTES. Every numeric constraint traces to Primer3/Addgene verbatim; deviations (max len 25, homopolymer 4, GC 40–60, pair Tm ≤5) are documented and within/stricter than sourced ranges.

- **Note 1 (GC clamp):** the default design does not enforce the Addgene 3' GC clamp; it is only available via opt-in `Avoid3PrimeGC` (and even then merely requires ≥1 G/C in the last 2 bases). This is acceptable because Primer3's own default `PRIMER_GC_CLAMP=0` also imposes no clamp — confirmed from the live manual this session. Documented divergence, not a defect.
- **Note 2 (TestSpec imprecision):** the TestSpec "Evidence Summary" lists SantaLucia (1998) as the basis for *Tm calculation*. The actual `CalculateMeltingTemperature` uses Wallace (N<14) / Marmur-Doty (N≥14); SantaLucia NN parameters are used only by `Calculate3PrimeStability`. This is a labelling imprecision in the spec; Tm-formula correctness is owned by PRIMER-TM-001 / PRIMER-STRUCT-001, not this unit. The hard-constraint enforcement validated here is unaffected.

---

## Stage B — Implementation

### Code path reviewed

- `DesignPrimers` (lines 40–115): validation (`targetStart<0 || targetEnd≥Length || targetStart≥targetEnd` → ArgumentException); forward search `[targetStart−200, targetStart)`; reverse search `[targetEnd, targetEnd+200)`; picks highest-`Score` candidate per side; compatibility = `tmDiff ≤ 5 && !HasPrimerDimer`; product size at line 106.
- `EvaluatePrimer` (lines 120–188): computes GC/Tm/homopolymer/dinuc/hairpin/3'-stability; appends an issue per violated constraint (lines 139–167); `IsValid = issues.Count==0`. Hard constraints enforced exactly with the validated thresholds.
- `GeneratePrimerCandidates` (lines 480–500): reverse branch applies `ReverseComplement` to the template substring before evaluation — consistent with the `DesignPrimers` reverse loop (lines 76–78).

### Realisation checks

- **Hard-constraint filter** (lines 139–152): length, GC%, Tm, homopolymer, dinuc each compared against the parameter struct with the validated ranges. ✅
- **Pair Tm filter** (lines 103–104): `Math.Abs(fwdTm−revTm) ≤ 5.0` gates `IsValid`. ✅
- **Reverse = reverse complement** (lines 76–78, 494–495): `Position` = template start, `Sequence` = revcomp. ✅
- **Amplicon size** (line 106): `reverse.Position + reverse.Length − forward.Position`. ✅
- **Ranking is a declared heuristic**: `CalculatePrimerScore` (lines 594–619) — 100 minus weighted deviations from optimal length/Tm/50%-GC/homopolymer plus a +5 GC-clamp bonus — selects the best candidate but does NOT relax any hard constraint (only `IsValid` candidates enter the pool, lines 61–62, 79–80). Heuristic ranking is acceptable per protocol; hard constraints verified enforced independently of score.
- **No-valid-pair path** (lines 93–100): returns `(null,null,false,"Could not find valid primers…",0)`. ✅

### Cross-verification vs code

The Stage-A hand table was reproduced against the actual code path: the standard template's accepted primers (24/25 bp) satisfy all constraints and the rejected shorter variants fail the Tm gate exactly as computed. `dotnet test --filter ~PrimerDesign` → **250 passed / 0 failed** this session.

### Edge cases (traced / tested)

| Edge case | Behaviour | Test |
|-----------|-----------|------|
| start ≥ end | ArgumentException | `_TargetEndBeforeStart_ThrowsArgumentException` |
| target beyond template | ArgumentException | `_TargetBeyondTemplate_…` |
| negative coords | ArgumentException | `_NegativeCoordinates_…` |
| null template | NullReferenceException (documented) | `_NullTemplate_ThrowsException` |
| empty primer seq | `IsValid=false` | `EvaluatePrimer_EmptySequence_HandledGracefully` |
| region < min length | empty candidate list | `GeneratePrimerCandidates_EmptyRegion_ReturnsEmpty` |
| length boundary 17/26 | issue reported | `EvaluatePrimer_LengthOutsideRange_ReportsIssue(17,26)` |
| GC boundary 100%/0% | issue reported | `EvaluatePrimer_GcOutsideRange_ReportsIssue` |
| Tm exactly at min/max | no Tm issue | mutation-killer boundary tests |

### Test-quality audit

- Positive `IsValid==true` assertions are non-vacuous (template engineered to pass all constraints; verified via the Tm hand table that this is genuinely achievable, not trivially true).
- `GeneratePrimerCandidates_Reverse_SequenceIsReverseComplement` asserts the actual revcomp identity per candidate.
- Boundary/mutation-killer tests assert exact thresholds (17/26 length, 0/100% GC, min/max Tm).

### Findings / defects

None. No reverse-complement omission, no amplicon off-by-one, all hard-constraint thresholds match the validated description.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — constraints sourced to Primer3/Addgene/Wikipedia; deviations documented & within range. Two minor notes: (1) default design omits the Addgene 3' GC clamp, consistent with Primer3 default `GC_CLAMP=0`; (2) TestSpec mislabels Tm basis as SantaLucia (actual: Wallace/Marmur-Doty; SantaLucia used only for 3'-end ΔG).
- **Stage B:** PASS — code faithfully realises the validated description; hard constraints enforced exactly; ranking is a declared heuristic that never relaxes a hard constraint.
- **State:** CLEAN. No code changes required.
- **Tests:** `--filter ~PrimerDesign` → 250 passed / 0 failed (no code touched; full suite not re-run since nothing changed).
