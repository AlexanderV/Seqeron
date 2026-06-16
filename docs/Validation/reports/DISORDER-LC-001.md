# Validation Report: DISORDER-LC-001 — Low-Complexity Region Detection in Protein Sequences (SEG)

- **Validated:** 2026-06-16   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.PredictLowComplexityRegions(string, int triggerWindow, double triggerThreshold, double extensionThreshold, int minLength)` (+ private `CalculateShannonEntropy`, `ClassifyLowComplexityType`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo's own Evidence)
- **GCG/Weizmann SEG help** — https://bip.weizmann.ac.il/education/materials/gcg/seg.html (WebFetch). Confirmed verbatim: `-WINdow=12`, `-LOWcut=2.2` (K1), `-HIGhcut=2.5` (K2); complexity in **bits/residue**; max amino-acid complexity = log₂(20) = **4.322** ("each character would add 4.322 bits … log(base 2) 20"); Stage 1 trigger = "complexity ≤ cutoff"; Stage 2 extension into overlapping segments ≤ cutoff.
- **`ncbi-seg` Ubuntu manpage** — https://manpages.ubuntu.com/manpages/bionic/man1/ncbi-seg.1.html (WebFetch). Confirmed W=12, K1=2.2 bits, K2=2.5 bits. Important clarification quoted verbatim: complexity "is defined by **equation (3) of Wootton & Federhen (1993)**"; trigger windows of length W and complexity ≤ K1 are "extended into a contig … by merging with extension windows … ≤ K(2)"; and a **third optimization stage** reduces each raw segment to "a single optimal low-complexity segment … the lowest value of the probability **P(0)** (equation (5))".
- **NCBI `blast_seg.c`** — https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_seg.c (WebFetch). Confirmed verbatim constants: `const int kSegWindow = 12; const double kSegLocut = 2.2; const double kSegHicut = 2.5;`. The file carries an `lnfact[]` factorial-log table → the reference tool uses the multinomial/factorial complexity (WF eq 3/5) for optimization, plus a Shannon-entropy (`s_Entropy`) trigger measure.
- **WebSearch** (Wootton & Federhen 1993, Comput. Chem. 17(2):149–163) — confirmed the paper defines "local compositional complexity" and that SEG-family tools "use a complexity measure based on Shannon entropy"; also the documented caveat that "Shannon entropy is not a good complexity measure for protein sequences" (the reason real SEG adds the P0 optimization).

### Formula check
Complexity per window = Shannon entropy of the residue composition, bits/residue:
H = −Σᵢ pᵢ·log₂(pᵢ), pᵢ = nᵢ/L. Independently hand-computed reference values (L=12):

| Composition | Distinct | H (bits) |
|---|---|---|
| 12×1 (`QQQ…`) | 1 | 0.000000 |
| 6+6 (`AAAAAALLLLLL`) | 2 | 1.000000 |
| 3+3+3+3 (`AAABBBCCCDDD`) | 4 | 2.000000 |
| 12 distinct (`ACDEFGHIKLMN`) | 12 | 3.584963 (= log₂12) |
| 20 equal | 20 | 4.321928 (= log₂20) |

All match the repo's described values and the GCG max-complexity statement.

### Edge-case semantics
Sequence shorter than W → no full window → empty (sourced). Window of W distinct residues → H=log₂W≈3.585 > K2 → never flagged (sourced). Homopolymer → H=0 → always triggers (sourced). All defined and source-backed.

### Independent cross-check
Hand-derived window entropies (above) reproduced from the sourced formula. Defaults W/K1/K2 verified verbatim in two independent reference sources (GCG help + NCBI `blast_seg.c`).

### Findings / divergences (→ PASS-WITH-NOTES)
1. **Complexity measure is the Shannon-entropy form, not WF equation (3) verbatim.** Real SEG's eq (3) is a multinomial complexity K₂=(1/L)·log_N(N!/Πnᵢ!); the canonical tool further reduces raw segments via the P0 significance test (eq 5). The repo uses the bits/residue Shannon entropy that NCBI `s_Entropy` and the GCG "bits/residue" docs describe, with the same 2.2/2.5 calibrated cutoffs, and **explicitly documents** the omission of WF eq-3 complexity + P0 optimization (algorithm doc §5.3 "Intentionally simplified / Not implemented"). This is a legitimate, disclosed simplification — boundaries can differ from NCBI SEG on mixed-complexity edges, but match on homopolymer/biased inputs. Note documented, not a Stage-A failure.
2. **`Type` label (`"X-rich"`/`"X/Y-rich"`)** is a repository presentation extension (>50% dominant rule), not part of SEG — already registered as Assumption #1. Does not affect coordinates.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`
- `CalculateShannonEntropy` (L280): counts A–Z, p=count/length, H−=p·log₂p — exact Shannon entropy normalized by full window length. Matches Stage A.
- `PredictLowComplexityRegions` (L497): null-guard → upper-case → short-circuit if len<W → Phase 1 trigger (≤K1) → collect spans → Phase 3 greedy extension (whole growing segment ≤K2) → merge adjacent (≤lastEnd+1) → minLength filter + classify. Realises the documented two-stage scan.

### Formula realised correctly?
Yes. I re-implemented the entire detector independently in Python (formula + spec only, not reading the C#) and reproduced **every** expected value in the test file:

| Case | Input | Expected (test) | Independent recompute |
|---|---|---|---|
| M1 | 26×Q | (0,25) | (0,25) ✓ |
| M2 | all-20 ×2 | empty | empty ✓ |
| M3 | AAABBBCCCDDD | (0,11) | (0,11) ✓ |
| M3′ | same, K1=0.5 | empty | empty ✓ |
| M5 | 12A+12L | (0,23) | (0,23) ✓ |
| M6 | 20Q+60spacer+20A | (0,34),(67,99) | (0,34),(67,99) ✓ |
| S1 | 12Q minLen=15 | empty | empty ✓ |
| S2 | 6Q | empty | empty ✓ |
| S3 | 12Q | (0,11) | (0,11) ✓ |

**M6 boundary entropies verified as genuine threshold crossings** (not arbitrary): H(0..34)=2.4939 ≤ 2.5 < H(0..35)=2.6077; H(67..99)=2.4250 ≤ 2.5 < H(66..99)=2.5452; spacer windows H=3.585 > K2 (gap preserved). These confirm the test's expected values trace to the sourced formula, not to a code echo.

### Variant/delegate consistency
Single public overload; defaults bound to `SegTriggerWindow=12 / SegTriggerComplexity=2.2 / SegExtensionComplexity=2.5` (verified verbatim against NCBI). No `*Fast` variant.

### Test quality audit (HARD gate)
Pre-existing fixture had 18 tests, all passing, with exact sourced boundaries/types and explicit entropy-evidence messages — **not** code echoes (each expected value independently reproduced above). Gate findings and fixes:
- **Coverage gaps closed this session** (logic not previously exercised):
  - `X/Y-rich` label branch (no residue >50%): added **M8** — dipeptide segment (0,23), A=L=50% → asserts `"A/L-rich"`.
  - `triggerWindow` custom parameter: added **M9** — W=4 on `AAABBBCCCDDD` (length-4 windows H∈{0.811,1.0} ≤ K1) → (0,11).
  - `extensionThreshold` custom parameter: added **M10** — K2=2.0 on the M6 input → (0,28),(71,99); boundaries verified as real crossings (H=1.877 ≤ 2.0 < 2.026).
- No green-washing: every added assertion uses **exact** sourced/hand-derived values (no Greater/Contains/range/widened tolerance); nothing skipped or weakened.
- Honest green: **full unfiltered suite 6612 passed / 0 failed** (1 pre-existing skipped benchmark, unrelated); changed test file builds warning-free.

Fixture now 21 tests; all MUST/SHOULD/COULD cases + INV-1/INV-2 properties + null/empty error paths + both label branches + all three tunable parameters covered.

### Findings / defects
None in the implementation. The only "deviations" are the disclosed Shannon-entropy-vs-WF-eq3 simplification and the cosmetic `Type` label — both documented in the algorithm doc and Assumption Register, neither affects the validated coordinate contract.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — biology/maths correct; defaults and entropy formula confirmed against two independent reference sources; one documented simplification (Shannon-entropy trigger in place of WF eq-3 complexity + P0 optimization), disclosed in the algorithm doc.
- **Stage B: PASS** — code faithfully realises the validated formula; all expected values independently reproduced from the spec; coverage gaps (X/Y-rich label, custom triggerWindow, custom extensionThreshold) closed with exact sourced values.
- **End-state: CLEAN / FIXED** — no implementation defect; test coverage gaps fully fixed; `dotnet build` 0 errors, full suite Failed: 0.
- **Test-quality gate: PASS.**
- No outstanding follow-ups for this unit. (A future enhancement — adding the WF eq-3 complexity + P0 optimization for NCBI-grade boundaries — is out of scope and already noted as a documented limitation, not a defect.)
