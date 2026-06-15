# Validation Report: PROTMOTIF-CC-001 — Coiled-Coil Prediction (heptad-repeat a/d hydrophobic-core detection)

- **Validated:** 2026-06-16   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.PredictCoiledCoils(string proteinSequence, int windowSize = 28, double threshold = 0.5)` → `IEnumerable<(int Start, int End, double Score)>`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (independently retrieved, not trusting the repo Evidence)

1. **Truebestein & Leonard (2016), "The Structure and Topology of α-Helical Coiled Coils", PMC7122542** — WebFetched.
   - Verbatim: *"the core-forming positions (a and d) are usually occupied by hydrophobic residues, whereas the remaining, solvent-exposed positions (b, c, e, f, and g) are dominated by hydrophilic residues."*
   - Residue preferences: dimers favour Ile/Val at **a**, Leu at **d**; tetramers reverse. Confirms I/L/V are the canonical a/d core residues.
2. **Wikipedia "Coiled coil" / "Heptad repeat"** (citing Mason & Arndt 2004; Chambers 1990) — WebSearch with verbatim extraction.
   - *"The positions in the heptad repeat are usually labeled abcdefg, where a and d are the hydrophobic positions, often being occupied by isoleucine, leucine, or valine."*
   - Also: *"mostly hydrophobic residues (Ala, Ile, Leu, Met, Val), or aromatic hydrophobic side chains (Phe, Trp and Tyr), are used at the a and d sites."* — i.e. {I,L,V} is the *named core* but not the *only* set.
3. **Lupas, Van Dyke, Stock (1991), Science 252:1162–1164** — WebSearch verbatim.
   - *"Because a window can be assigned seven different heptad repeat frames and a residue can occupy 28 positions in the gliding window, there are 196 preliminary scores for each residue, with the highest score being selected for each residue."*
   - *"gliding window of 28 residues"*, chosen because *"the shortest peptides still exhibiting a stable coiled-coil structure … are between four and five heptads long."*
   - COILS preliminary score = **geometric average** of per-position residue *relative frequencies* (a 21×20 PSSM), scaled against reference databases to a probability.

### Formula check

The documented model (doc §2.2, impl) is:
- heptad position `p(k,r) = (k − r) mod 7`; core positions `a = 0`, `d = 3` — matches abcdefg labelling.
- core residue set `{I, L, V}` — exactly the residues named verbatim in source (2).
- `occ(i,r) = (# a/d in window with residue ∈ {I,L,V}) / (# a/d in window)`; `score(i) = max over r∈{0..6} occ(i,r)` — the "7 frames, take the max" rule is verbatim Lupas (1991).
- window 28 = 4 heptads, min region 21 = 3 heptads — sourced (Lupas window; Mason & Arndt multi-heptad).

**Notes (why PASS-WITH-NOTES, not PASS):**
- **N-1 (simplification, documented):** This is NOT COILS. COILS scores a residue by the *geometric mean of position-specific residue frequencies* from a 21×20 PSSM and converts to a probability; this unit computes a plain a/d *occupancy fraction*. The doc (§2.5, §5.3) and Evidence state this openly and decline to fabricate the PSSM weights (which were not retrievable). This is a legitimate, source-grounded design choice, not a defect — but the algorithm name "coiled-coil prediction" should be read as the a/d-periodicity heuristic, not COILS.
- **N-2 (residue set, documented assumption):** the core set is restricted to {I,L,V}. Sources also name A, M, F, W, Y at a/d. The restriction is the verbatim-named set and avoids untraceable constants; documented as Assumption 1.

### Edge-case semantics check (all sourced)

| Case | Expected | Source |
|------|----------|--------|
| length < windowSize | empty (no full window) | Lupas window rule |
| no {I,L,V} residues | empty (occupancy 0 < 0.5) | core set {I,L,V} |
| region < 21 residues (3 heptads) | rejected | Mason & Arndt multi-heptad |
| off-frame coiled coil | found via 7-register max | Lupas seven frames |
| null / empty | empty | standard validation |

### Independent cross-check (numbers)

An independent Python reimplementation of the occupancy model (written this session, not derived from the C# code) reproduced every expected value — see Stage B table.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs:885–1034`:
- constants (885–910): `HeptadLength=7`, `HeptadPositionA=0`, `HeptadPositionD=3`, `CoiledCoilCoreResidues={I,L,V}`, window 28, threshold 0.5, `MinCoiledCoilRegion=21`.
- `PredictCoiledCoils` (924–977): null/empty/sub-window guard; uppercases; builds occupancy profile; single forward scan grouping contiguous ≥-threshold windows; emits regions (mid-loop drop branch + trailing branch).
- `BestHeptadOccupancy` (983–1012): max over 7 registers of (a/d ∈ {I,L,V})/(a/d count); skips register if no a/d in window.
- `BuildRegion` (1019–1027): maps window run `[first,last]` → residues `[first, last+W−1]`, rejects if span < 21.
- `Mod` (1030–1034): correct non-negative modulo.

Code realises the validated formula exactly (a=0, d=3, set {I,L,V}, 7 registers, max, peak-score per run, inclusive 0-based Start/End, min-region filter). No precision/overflow concerns (integer counts; one division by a positive count).

### Cross-verification table (independent Python model vs C# output — both agree)

| ID | Input | Independent expected | C# actual | Match |
|----|-------|----------------------|-----------|-------|
| M1 | `LAALAAA`×5 (35) | (0,34,1.0) | (0,34,1.0) | ✅ |
| M2 | `G`×40 | ∅ | ∅ | ✅ |
| M3 | `LAALAAA`×3 (21) | ∅ | ∅ | ✅ |
| M4 | `AA`+`LAALAAA`×5 (37) | (0,36,1.0) | (0,36,1.0) | ✅ |
| M5 | `LAALAAA`×4 (28) | (0,27,1.0) | (0,27,1.0) | ✅ |
| S1 | `LAAAAAA`×5 (35) | (0,34,0.5) | (0,34,0.5) | ✅ |
| S2 | S1, thr=0.5001 | ∅ | ∅ | ✅ |
| S5 | `LAALAAA`×5+`G`×35 (70) | (0,48,1.0) | (0,48,1.0) | ✅ |
| S3 | null | ∅ | ∅ | ✅ |
| S4 | "" | ∅ | ∅ | ✅ |
| C1 | `laalaaa`×5 | (0,34,1.0) | (0,34,1.0) | ✅ |
| **S6** (new) | `LAALAAA`×5+`G`×40+`LAALAAA`×5 (110) | (0,48,1.0),(58,109,1.0) | same | ✅ |
| **S7** (new) | `G`×7+`LAALAAA`+`G`×7 (21), W=7, thr=0.99 | ∅ (run [4..13] len 10 < 21) | ∅ | ✅ |
| **S8** (new) | `LAALAAA`×4 (28), W=14 | (0,27,1.0) | (0,27,1.0) | ✅ |

S5/S6 drop verified by the explicit occupancy profile: `1.0×8, 0.875×3, 0.75×4, 0.625×3, 0.5(idx 18–21), 0.375…→0`; last window ≥0.5 is index 21 → end 21+28−1 = 48.

### Variant/delegate consistency

Single public method; no `*Fast`/delegate variants. Private helpers `BestHeptadOccupancy`/`BuildRegion` are consistent with the canonical scan.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** all expected regions are exact closed-form values from the validated definition (`.EqualTo(...).Within(1e-10)` on Score, exact Start/End), independently reproduced in Python this session. A count-only or fixed-register or wrong-core-set implementation fails M1/M4/S1. Not tautologies.
- **No green-washing:** no weakened assertions, no widened tolerances, no skips, no expected-value adjustments to match output.
- **Coverage — gaps found and fixed (test-only, 0 code change):** the original 11 tests never exercised (a) the **multi-region** drop-and-continue branch (INV-03 non-overlap/increasing-Start), (b) a **non-default windowSize**, or (c) the **`BuildRegion` min-region rejection** path (unreachable at W=28>21). Added **S6** (two regions), **S7** (custom W=7 run below MinRegion → rejected), **S8** (custom W=14 still detects). Now every public-method parameter and every Stage-A branch (sub-window guard, no-core, off-frame max, threshold ≥/strict-below, single-window, multi-region, min-region reject, null/empty, lowercase) is covered.
- **Honest green:** full unfiltered `dotnet test` = **6582 passed, 0 failed** (1 pre-existing unrelated skip, `MFE_Benchmark_AllScenarios`); `dotnet build` 0 warnings / 0 errors. Fixture 11 → 14.

### Findings / defects

No algorithm or code defect. Two pre-existing Stage-A documentation notes (N-1 simplification of COILS, N-2 {I,L,V} restriction) are already disclosed in the doc/Evidence/TestSpec and require no change. The only action taken was closing three test-coverage gaps with sourced expected values.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — biology/maths correct; the model is an honestly-documented a/d-occupancy *simplification* of COILS (no PSSM) and restricts the core set to the verbatim-named {I,L,V}.
- **Stage B:** PASS-WITH-NOTES — code faithfully realises the validated description; test-coverage gaps (multi-region, custom window, min-region rejection) fixed in-session with externally-reproduced values.
- **End-state:** ✅ CLEAN — no defect; coverage gaps fixed; build + full suite green.
- **Test-quality gate:** PASS (exact sourced values, all branches covered, honest green 6582/0).
