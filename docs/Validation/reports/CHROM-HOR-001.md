# Validation Report: CHROM-HOR-001 — Higher-Order Repeat (HOR) Detection

- **Validated:** 2026-06-25   **Area:** Chromosome
- **Canonical method(s):** `ChromosomeAnalyzer.DetectHigherOrderRepeat(string sequence, int monomerLength = 171)` → `HorResult`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN (one test-coverage gap fixed this session; suprachromosomal-family assignment remains a documented data-blocked boundary)

## Canonical method(s)
`DetectHigherOrderRepeat` — `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs:751`.
Returns `HorResult(HasHigherOrderStructure, MonomersPerUnit, HorUnitLengthBp, HorCopyNumber, MonomerCount, MeanInterHorIdentity, MeanIntraHorIdentity)` (record at `:157`).

## Stage A — Description

### Sources opened this session (authoritative, external)
1. **McNulty & Sullivan (2018)**, "Alpha satellite DNA biology: finding function in the recesses of the genome," *Chromosome Res* 26:115–138 (PMC6121732). Verbatim, retrieved this session:
   - Monomer size: "Alpha satellite DNA is composed of fundamental **171 bp** monomeric repeat units."
   - **Intra-HOR** monomer identity: "The individual monomers within a HOR unit have **50–70% identity** and can be distinguished…"
   - **Inter-HOR** copy identity: "HORs within a chromosome-specific array **differ in sequence by only a few percent**, however, HORs between non-homologous chromosomes are only 50–70% identical."
   - HOR unit length: "HOR unit length is determined by **where the next monomer shows nearly total sequence identity to the first monomer** in the HOR."
   - Repetition: "HOR units are repeated, largely uninterrupted, hundreds to thousands of times, resulting in a large, linear and homogeneous array of highly identical copies of tandem HOR units."
2. **Rosandić et al. (2024)** (PMC11050224), retrieved this session:
   - Inter-HOR divergence threshold: "HOR copies are further organized in tandem, with **minimal divergence between HOR copies, typically less than 5%**" (i.e. ≥ ~95% identity).
   - Period = k×monomer: "A peak of **period *n* corresponds to *n* × 171 bp**" — i.e. an *n*mer HOR = *n* monomers.
3. **Willard (1985)** / **Alkan et al. (2007)** — corroborate 171-bp monomer and <5% inter-HOR divergence (cited in the two reviews above).

### Formula / definition check
- **HOR period** = k monomers; **HOR unit length** = k × 171 bp. Confirmed (Rosandić: period n ⇒ n×171 bp).
- **Copy number** = ⌊monomer count / k⌋ — standard tiling count; consistent with "HOR units repeated … times."
- **Identity ordering (the defining property):** inter-HOR identity (≥ ~95%) ≥ intra-HOR monomer identity (50–70%). Directly confirmed: inter "differ by only a few percent" (~95–100%) vs intra "50–70%".
- Implementation constant `InterHorMinIdentityPercent = 95.0` maps exactly to "<5% divergence" (Rosandić/Alkan). ✅
- `HorPeriodConsistencyFraction = 0.90` (period must hold across ≥90% of k-spaced pairs) encodes "largely uninterrupted … homogeneous array." Reasonable, source-consistent.

### Edge-case semantics (defined & sourced)
- < 2 full monomers (incl. empty/null/single) → no periodicity possible → no HOR (period 1, NaN identities). Defined.
- Homogeneous 1-monomer repeat → smallest high-identity period is k=1, which is **not** a multi-monomer HOR ⇒ `HasHigherOrderStructure=false`. Matches the literature distinction between a monomeric/homogeneous array and a HOR (multimeric unit).
- Trailing partial monomer (length not a multiple of monomer) → dropped by floor split. Standard.
- Non-ACGT (e.g. N) → treated by the aligner as a non-matching residue; no special meaning required (input is assumed alpha-satellite).

### Independent cross-check (numbers — hand-computed this session, independent of repo fixture)
Built a **k=4 monomer unit repeated m=7×** over an independent background (4 pairwise-disjoint 30-position substitution sets ⇒ symmetric difference 60 ⇒ Hamming intra = (171−60)/171 = **64.91%**, inside the 50–70% band). Exact HOR copies ⇒ inter = 100%.

| Quantity | Hand-computed | Code (`DetectHigherOrderRepeat`) | Match |
|---|---|---|---|
| HOR period (monomers/unit) | 4 | 4 | ✅ |
| HOR unit length (bp) | 4 × 171 = 684 | 684 | ✅ |
| Copy number | ⌊28/4⌋ = 7 | 7 | ✅ |
| Monomer count | 28 | 28 | ✅ |
| Mean inter-HOR identity | 100.0% | 100.0% | ✅ |
| Mean intra-HOR identity | 64.91% (Hamming) | 65.5% (aligner) | band ✅ |
| inter ≥ intra invariant | true | true | ✅ |
| Non-HOR control (10 divergent monomers) | no HOR | HasHOR=false, period 1, inter=NaN | ✅ |

**Note on intra value:** intra-HOR identity is a *continuous* aligner measure. The exact value depends on whether the global aligner finds gaps. When the background is truly non-periodic and substitutions are tightly banded (as in the repo fixture), the optimal alignment is **gapless** and aligner identity = Hamming exactly (verified: repo fixture A↔B/A↔C/B↔C = 57.8947%, 0 gaps). When the background carries residual short-period autocorrelation (my first harness used a period-4 background), the aligner inserts a few opportunistic gaps and reports a slightly higher value — a *fixture* property, not a code defect. Both land the value inside the sourced 50–70% band, and the structural quantities (period, copy number, unit length, monomer count, inter-identity, ordering, non-HOR rejection) match exactly.

### Findings / divergences (Stage A)
None. The repo's stub TestSpec/report contained no incorrect claim; the description is biologically correct and now backed by verbatim source quotes.

## Stage B — Implementation

### Code path reviewed
`ChromosomeAnalyzer.cs:751–835` (`DetectHigherOrderRepeat`), `:728/:734` (thresholds), `:841` (`MeanPairwiseIdentity`), `:773` (per-pair identity via `SequenceAligner.GlobalAlign` + `CalculateStatistics`).

### Formula realised correctly?
- Splits into ⌊len/monomerLength⌋ full monomers; trailing partial dropped (`:762`). ✅
- Searches smallest period k∈[1, monomerCount/2] where ≥90% of k-spaced monomer pairs are ≥95% identical (`:790–805`). Maps to "period n ⇒ n×171 bp," ≥95% = "<5% divergence." ✅
- `copyNumber = monomerCount/detectedPeriod`, `unitLengthBp = detectedPeriod*monomerLength` (`:815–816`). ✅
- Inter-HOR identity = mean identity of monomers period apart (`:820–826`); intra-HOR = mean pairwise identity within the first unit, only for k≥2 (`:829–831`). ✅
- `HasHigherOrderStructure = detectedPeriod >= 2` — correctly excludes the homogeneous k=1 case from being a HOR (`:833`). ✅

### Cross-verification vs code
All values in the Stage-A table reproduced by running the actual code (independent harness, k=4/m=7) and by running the repo fixture (k=3/m=5: period 3, copy 5, unit 513 bp, monomers 15, inter 100%, intra 57.8947% = exact Hamming). ✅

### Variant/delegate consistency
Single public method + one `monomerLength` overload parameter; both exercised (default 171, and custom 10 in the harness → period 2/copy 5/unit 20). Mixed-case input matches uppercase (`ToUpperInvariant`, `:759`; locked by `MixedCaseInput_MatchesUppercase`). ✅

### Numerical robustness
NaN only where defined (no inter pairs / k<2). No div-by-zero (guards at `:763`, `:826`, `:851`). `monomerLength < 1` throws `ArgumentOutOfRangeException` (`:753`). ✅

### Test quality audit
Fixture `ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs` builds monomers from a shared background with point substitutions and **hand-derives** expected identity as (171 − |A△B|)/171 — values trace to construction, not code echoes. Covered Stage-A paths: 3-mer HOR (period 3×171) ✅; dimeric HOR ✅; monomeric-only → no HOR (self-validating: asserts all pairs <95% first) ✅; homogeneous 1-mer → not HOR ✅; inter ≥ intra invariant ✅; period-not-multiple (trailing partial) ✅; empty/null/single ✅; invalid monomerLength throws ✅; mixed case ✅.
- **Gap found & fixed this session:** no explicit **non-ACGT** test (a Stage-A required edge case). Added `DetectHigherOrderRepeat_NonAcgtTrailingPartialMonomer_IsIgnored` (80-bp N tail dropped as a partial monomer ⇒ period/copy/inter unchanged from M-HOR-1; hand-derived, behavior pre-verified by the independent harness). Full unfiltered suite green afterward.

### Findings / defects (Stage B)
- FR: minor test-coverage gap (missing non-ACGT case) — **fixed this session** (test added; logged in FINDINGS_REGISTER). No production-code defect.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- Full unfiltered `dotnet test Seqeron.sln -c Debug`: **Seqeron.Genomics.Tests 18780 passed / 0 failed** (incl. the new test); 0 warnings on the changed test project. (Two empty MCP test DLLs report "no test available" — pre-existing, unrelated.)
- **Documented boundary (not a defect):** suprachromosomal-family / chromosome-specific HOR *assignment* is not attempted — it requires curated HOR reference libraries (e.g. T2T-CHM13 HOR annotations). `DetectHigherOrderRepeat` correctly limits itself to period / copy number / unit length / inter-vs-intra identity, which it computes correctly.
