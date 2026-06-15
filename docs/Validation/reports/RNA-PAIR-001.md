# Validation Report: RNA-PAIR-001 — RNA Base Pairing

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CanPair(char,char)`, `RnaSecondaryStructure.GetBasePairType(char,char)`, `RnaSecondaryStructure.GetComplement(char)` (delegates to `SequenceExtensions.GetRnaComplementBase`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-coverage gap fixed in-session; no algorithm defect)

## Stage A — Description

### Sources opened this session (independent retrieval)

1. **Wikipedia — Wobble base pair** (https://en.wikipedia.org/wiki/Wobble_base_pair) — fetched 2026-06-16.
   - Verbatim: "The four main wobble base pairs are guanine–uracil (G–U), hypoxanthine–uracil (I–U), hypoxanthine–adenine (I–A), and hypoxanthine–cytosine (I–C)." Only G–U is over the standard {A,C,G,U} alphabet.
   - Wobble pairing table: **G pairs with C or U; U pairs with A or G.**
   - Verbatim: a wobble base pair "does not follow Watson–Crick base pair rules" — confirms G–U is a *distinct* type, never Watson-Crick.
   - Primary citation confirmed: Crick FH (Aug 1966), *J Mol Biol* **19**(2):548–555.

2. **Wikipedia — Base pair** (https://en.wikipedia.org/wiki/Base_pair) — fetched 2026-06-16.
   - Canonical Watson-Crick RNA pairs: A•U and G•C; in RNA uracil substitutes thymine.
   - A•U has 2 hydrogen bonds; G•C has 3.
   - Pairing is reciprocal/complementary (symmetric): A•U ≡ U•A, G•C ≡ C•G.

3. **Wikipedia — Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation) — fetched 2026-06-16.
   - Full IUPAC-IUB (1970) complement table extracted: A↔T, C↔G, G↔C, U→A; W↔W, S↔S, M↔K, K↔M, R↔Y, Y↔R, B↔V, V↔B, D↔H, H↔D, N↔N.

4. **ViennaRNA RNAlib** (independent reference implementation; WebSearch on RNAlib `md.pair` array / default pairing) — confirms the allowed-pair set used in production RNA folding: "RNA helices normally contain only 6 out of the 16 possible combinations: the Watson-Crick pairs **GC, CG, AU, UA**, and the somewhat weaker wobble pairs **GU and UG**." This is exactly the implementation's six seeded ordered pairs.

### Formula / model check

The model is a finite truth table, not a numeric formula. Independently reconstructed truth table:

| b1 | b2 | CanPair | Type        | Source |
|----|----|---------|-------------|--------|
| A  | U  | true    | WatsonCrick | Base pair |
| U  | A  | true    | WatsonCrick | Base pair (symmetry) |
| G  | C  | true    | WatsonCrick | Base pair |
| C  | G  | true    | WatsonCrick | Base pair (symmetry) |
| G  | U  | true    | Wobble      | Wobble / Crick 1966 / ViennaRNA |
| U  | G  | true    | Wobble      | Wobble / Crick 1966 / ViennaRNA |
| all other {A,C,G,U}² | | false | null | A pairs only with U; C only with G |

This matches the description in `docs/algorithms/RnaStructure/RNA_Base_Pairing.md` §2.2 and the TestSpec dataset exactly. The six-pair set is independently confirmed by ViennaRNA (a peer reviewed reference tool), not merely by the repo's own artifacts.

### Definitions & conventions

- Alphabet {A,C,G,U}, case-insensitive (non-correctness-affecting normalization). ✔ sourced/standard.
- T is a DNA base: not an RNA pairing input → `CanPair`/`GetBasePairType` false/null. Consistent with all pairing sources defining pairing over {A,C,G,U}. ✔
- RNA complement: A→U, U→A, G→C, C→G, T→A (T treated as U). Matches IUPAC table (A↔U in RNA variant) and Biopython `complement_rna("CGAUT")="GCUAA"`. ✔

### Edge-case semantics

- Symmetry (order independence): sourced from reciprocal pairing (Base pair). ✔
- G–U is Wobble, never WatsonCrick: sourced (Wobble article + Crick 1966). ✔
- Out-of-ASCII / non-base char → false/null, no exception: robustness, not in source (acceptable as defined contract). ✔
- Full IUPAC degenerate complement (W,S,M,K,B,D,H,V) is sourced from IUPAC-IUB (1970) — this session retrieved the complete table.

### Independent cross-check (numbers)

- ViennaRNA RNAlib default allowed-pair set = {GC, CG, AU, UA, GU, UG} — identical to `BuildPairLookup`.
- Biopython `complement_rna("CGAUT") → "GCUAA"` (C→G, G→C, A→U, U→A, T→A) — identical to `GetRnaComplementBase`.
- Full IUPAC complement self/cross map (W↔W, S↔S, M↔K, R↔Y, B↔V, D↔H, N↔N) — identical to `GetRnaComplementBase` arms.

### Findings / divergences

None. The biology and conventions are correct and externally confirmed by two independent reference implementations (ViennaRNA, Biopython) plus primary citation (Crick 1966). **Stage A: PASS.**

## Stage B — Implementation

### Code paths reviewed

- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
  - `BuildPairLookup()` (lines 427-434): seeds `A→U, U→A, G→C, C→G` = 1 (WatsonCrick), `G→U, U→G` = 2 (Wobble). All other cells 0. **Exactly the sourced six-pair set.**
  - `CanPair` (443-448): upper-cases both inputs, bounds-checks `(b1|b2)<128`, returns `PairLookup[..] != 0`. Symmetry holds by construction (table seeded both orders).
  - `GetBasePairType` (453-464): same lookup, maps 1→WatsonCrick, 2→Wobble, else→null. INV-03/04 hold by shared table.
  - `GetComplement` (469): thin delegate.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs`
  - `GetRnaComplementBase` (175-194): switch over uppercase+lowercase arms for A,U,G,C,T and all 11 IUPAC degenerate codes, `_ => nucleotide` pass-through. **Matches the externally-retrieved IUPAC table and Biopython exactly.**

### Formula realised correctly?

Yes. The lookup table contents equal the independently-reconstructed truth table cell-for-cell. The complement switch equals the independently-retrieved IUPAC/Biopython mapping arm-for-arm. T treated as U → A in complement (DNA T not an RNA pairing input). No approximation.

### Cross-verification table recomputed vs code

Traced by hand against the code and confirmed by the green test run (every row of the Stage-A truth table and complement table asserted with exact sourced values). All match.

### Variant/delegate consistency

`GetComplement` delegates to `GetRnaComplementBase`; verified the delegate target is the validated helper. `CanPair` and `GetBasePairType` share `PairLookup`, guaranteeing INV-03 by construction.

### Test quality audit (HARD gate)

Canonical test file: `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_CanPair_Tests.cs`.

- **Sourced, not code-echoes:** all assertions check exact externally-sourced values (true/false, WatsonCrick/Wobble/null, exact complement chars). A deliberately-wrong implementation (e.g. G–U classified WatsonCrick, or A→T complement) would FAIL these — not green-wash-prone.
- **No weakened assertions:** all use `Is.EqualTo`/`Is.True`/`Is.False`/`Is.Null` exact matches; no ranges/`Greater`/`Contains` where an exact value is known; no skips.
- **Coverage gap found and fixed:** the original S3 test exercised only N, R, Y of the IUPAC degenerate complement branches, leaving the *sourced* arms S, W, K, M, B, D, H, V (and the non-IUPAC pass-through, and lowercase degenerate) untested — i.e. real, source-defined logic with no assertion. This is a Stage-B test defect per the gate ("cover all the logic"). **No algorithm defect** (all arms verified correct vs the retrieved IUPAC table).
  - Strengthened S3 to assert all 11 degenerate complements: N→N, R→Y, Y→R, W→W, S→S, M→K, K→M, B→V, V→B, D→H, H→D (exact IUPAC-IUB / Wikipedia values).
  - Added `GetComplement_LowercaseIupacDegenerate_ReturnsExpected` (r→Y, n→N, k→M) — locks the lowercase arms.
  - Added `GetComplement_NonIupacChar_PassesThroughUnchanged` ('-'→'-', 'X'→'X') — locks the documented §6.1 pass-through contract.
- **Edge cases covered:** symmetry (M8 full 4×4 cross-product), CanPair⇔type consistency (M9), case-insensitivity (S1), DNA-T non-pairing (S2), out-of-range char (C1), wobble≠WatsonCrick (M5/INV-04).
- **Honest green:** full unfiltered suite `dotnet test` = **Failed: 0, Passed: 6586** (was 6584 before; +2 new test methods plus added assertions). `dotnet build` = 0 errors (4 pre-existing NUnit2007 warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, not touched).

### Findings / defects

No algorithm defect. One test-coverage gap (untested sourced IUPAC-complement branches) — fixed in-session with exact sourced values.

## Verdict & follow-ups

- **Stage A: PASS** — model independently confirmed against Crick 1966, Wikipedia Base pair / Wobble, ViennaRNA RNAlib (allowed-pair set), and Biopython / IUPAC-IUB 1970 (complement table).
- **Stage B: PASS-WITH-NOTES** — code realises the validated model exactly; test-coverage gap on IUPAC degenerate complement branches fixed (no algorithm change required).
- **End-state: CLEAN** — defect (test gap) completely fixed in-session; build + full suite green.
- **Test-quality gate: PASS** — exact sourced expectations, no green-washing, all logic now covered, honest green (6586/0).
