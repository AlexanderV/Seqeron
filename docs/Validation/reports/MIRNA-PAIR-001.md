# Validation Report: MIRNA-PAIR-001 — MiRNA-Target Pairing Analysis

- **Validated:** 2026-06-15   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.AlignMiRnaToTarget(string,string)`, `GetReverseComplement(string)`, `CanPair(char,char)`, `IsWobblePair(char,char)` (+ private `CalculateDuplexEnergy`, `StackingEnergies`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one green-washed/contract-violating test + code gap found and fixed this session)

---

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **Wikipedia "Wobble base pair"** (cites Crick FHC 1966, J Mol Biol 19(2):548–555) and **ScienceDirect "Wobble base pair overview"** — confirmed: standard Watson-Crick pairs are **A·U and G·C**; the principal RNA wobble pair is **G·U**; G·U is explicitly *not* Watson-Crick and is "about as stable as an A-U base pair". → Confirms `CanPair` set {A-U, U-A, G-C, C-G, G-U, U-G} and `IsWobblePair` = {G-U, U-G} only.
2. **miRBase / RNAcentral search** — hsa-let-7a-5p (MIMAT0000062) mature sequence `UGAGGUAGUAGGUUGUAUAGUU`; seed (positions 2–8) `GAGGUAG`. Matches the spec's let-7a dataset. Watson-Crick reverse complement of `GAGGUAG` = `CUACCUC` (hand-verified: complement C-U-C-C-A-U-C, reversed → CUACCUC).
3. **Xia et al. (1998) Biochemistry 37:14719–14735** (Turner 2004 WC parameters), via the published Table-4 values — retrieved the full 10-stack Watson-Crick nearest-neighbor ΔG°37 set:
   AA/UU −0.93, AU/UA −1.10, UA/AU −1.33, CU/GA −2.08, CA/GU −2.11, GU/CA −2.24, GA/CU −2.35, CG/GC −2.36, GG/CC −3.26, GC/CG −3.42 (kcal/mol).
4. **NNDB Turner 2004 GU parameter page (search snippets)** — confirmed the unusual tandem-wobble values exist: 5′GU3′/3′UG5′ destabilising (~+1.29 in standard context), 5′UG3′/3′GU5′ (~+0.30/+0.47), with special stabilising context. (NNDB HTML pages 404 to the fetcher — same limitation as the authoring session; values cross-checked from search result text only, not re-opened in full.)

### Formula check

- **Pairing classification** (doc §2.2): Watson-Crick {A-U,U-A,G-C,C-G}; wobble {G-U,U-G}; else mismatch — matches sources 1.
- **Antiparallel orientation** (doc §2.2/§3.3): miRNA index `i` pairs target index `len−1−i` — matches Lewis (2005) reverse-complement seed convention. Verified the code uses exactly `target[target.Length-1-i]`.
- **Free energy** (doc §2.2): ΔG ≈ Σ Turner-2004 stacking over consecutive paired positions; magnitude declared "intentionally simplified" (no loop/bulge/init/terminal terms), only sign/ordering reliable. All 16 Watson-Crick `StackingEnergies` entries match the Xia-1998 values exactly (symmetric duplicates such as AA/UU = UU/AA are valid by strand-symmetry).

### Edge-case semantics

Empty/null → empty duplex / `""`; unequal lengths → overlap = min(len); DNA T → normalised to U; unknown base → `N`, never pairs. All sourced to the documented contract (§3.3, §6.1) and standard convention.

### Independent cross-check (numbers)

| Case | Source-derived expectation | Hand trace |
|------|----------------------------|-----------|
| RC(`GAGGUAG`) | `CUACCUC` (Watson-Crick RC, Lewis 2005) | `CUACCUC` ✓ |
| AlignMiRnaToTarget(AAAA,UUUU) | 4 WC, "\|\|\|\|" | ✓ |
| AlignMiRnaToTarget(GGGG,UUUU) | 4 wobble, "::::" (Crick) | ✓ |
| AlignMiRnaToTarget(AAAA,AAAA) | 4 mismatch (A pairs only U) | ✓ |
| AlignMiRnaToTarget(AGGU,AUCG) | 2 WC / 1 wobble / 1 mismatch, " \|:\|" | ✓ |
| GC-stem FreeEnergy sign | Σ(−3.42,−2.36)×alt = −20.76 ≤ 0 | ✓ (Python trace) |

### Findings / divergences (Stage A)

None substantive. The description is biologically and thermodynamically correct. Minor: the doc references the same Turner 2004 GU table; the GU wobble stacking numbers are context-dependent and not all independently re-opened this session (NNDB 404) — but the free-energy contract only claims sign/ordering, which is sound. **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`:
- `CanPair` L322, `IsWobblePair` L335, `GetReverseComplement` L294, `AlignMiRnaToTarget` L350, `CalculateDuplexEnergy` L415, `StackingEnergies` L478.

### Formula realised correctly?

Yes. Antiparallel pairing index, WC/wobble/mismatch classification, alignment symbols (`|`/`:`/space), ungapped overlap, and nearest-neighbor stacking sum all match the validated description. All 16 WC stacking constants match Xia 1998 exactly.

### Defect found & fixed this session

**FINDING (Stage B):** `CanPair`/`IsWobblePair` did **not** honour the documented "DNA T treated as U" contract (doc §3.1, §6.1; TestSpec M4 "A-T → true"). `CanPair('A','T')` returned **false** while `AlignMiRnaToTarget`/`GetReverseComplement` (which `.Replace('T','U')`) treat T as U — an inconsistency, and the MCP-exposed `CanPair` tool inherited the wrong behaviour.

The canonical test `CanPair_LowercaseAndDnaT_ReturnsTrue` was **green-washed**: its name + the spec promised DNA-T handling, but the body only asserted lowercase `a-u`/`g-c` and silently omitted the `A-T` assertion that would have failed. It would pass against the (wrong) implementation — exactly the gate-prohibited defect.

**Fix (code, matches the description — not the other way round):** added a private `NormalizeBase` (uppercase + T→U) used by `CanPair` and `IsWobblePair`, so they now honour the long-documented contract. Internal callers (`AnalyzeHairpin`, hairpin scan) pass pre-normalised sequences, so behaviour there is unchanged; the change only *adds* correct T-handling at the public/MCP boundary.

**Fix (tests):** `CanPair_LowercaseAndDnaT_ReturnsTrue` now asserts `CanPair('A','T')`, `('T','A')`, `('G','T')`, `('a','t')` == true; added `IsWobblePair_DnaT_TreatedAsU` (M4b: G-T/T-G wobble true, A-T wobble false). Also rewrote M15 (`AlignMiRnaToTarget_CountInvariant…`) which previously used a *perfect-complement* input while its comment claimed "mixed" — it now uses `AGGU`/`AUCG`, a genuinely mixed duplex, and asserts exact per-class counts (2 WC / 1 wobble / 1 mismatch) + alignment `" |:|"`, so the invariant test is no longer vacuous.

### Cross-verification table recomputed vs code

All rows in the Stage-A cross-check table reproduce against the actual code (canonical fixture green; GC-stem ΔG hand-trace −20.76 ≤ 0).

### Variant/delegate consistency

`AlignMiRnaToTarget` and `GetReverseComplement` normalise T→U; after the fix `CanPair`/`IsWobblePair` do too — now fully consistent. MCP `AnnotationTools.CanPair`/`IsWobblePair` delegate directly and inherit the corrected behaviour.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** exact values trace to Crick (wobble), Agarwal/Wikipedia (WC), Lewis (RC), Xia 1998 (stacking). ✓
- **No green-washing:** the one green-washed test (M4) was found and rewritten to the contract value; no assertions weakened; no tolerances widened; no skips. ✓
- **Coverage:** all four public methods + every Stage-A branch (WC/wobble/mismatch, empty miRNA/target, unequal length, DNA-T, RC involution, wobble⊆pair, FreeEnergy sign both directions). Free-energy *magnitude* is asserted by sign only — by design (NNDB pages 404; documented limitation), not a coverage gap in scope.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6543`; `dotnet build` 0 errors; changed files add 0 new warnings.

**Gate result: PASS.**

### Findings / defects

- FIXED: T-normalisation gap in `CanPair`/`IsWobblePair` + green-washed M4 test + vacuous M15 invariant test (all corrected and locked with sourced assertions).

---

## Verdict & follow-ups

- **Stage A:** PASS. **Stage B:** PASS-WITH-NOTES (defect found and completely fixed in-session).
- **End-state: CLEAN** — code now matches the documented contract, tests lock the sourced values, full suite green.
- Logged in FINDINGS_REGISTER. No outstanding work; free-energy magnitude remains intentionally simplified (sign-only validated) as documented.
