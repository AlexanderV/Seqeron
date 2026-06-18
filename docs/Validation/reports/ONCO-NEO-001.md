# Validation Report: ONCO-NEO-001 — Neoantigen Candidate Peptide Window Generation

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.GenerateNeoantigenPeptides(string wildTypeProtein, char mutantResidue, int mutationPosition, int minLength = 8, int maxLength = 11)` → `IReadOnlyList<NeoantigenPeptide>`; constants `MhcClassIMinPeptideLength = 8`, `MhcClassIMaxPeptideLength = 11`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Summary

ONCO-NEO-001 is the deterministic *windowing* step of MHC class I neoantigen prediction: given a single
somatic missense substitution (`P[p] → a`) in a protein `P` of length `L`, it enumerates every fixed-length
k-mer window (k = 8..11 by default) of the mutant protein `P'` that **spans** the mutated residue, each paired
with the wild-type k-mer at the same coordinates (the agretope). Binding affinity / IC50 (NetMHCpan) is
explicitly out of scope (caller-supplied, ONCO-MHC-001).

## Stage A — Description

### Sources opened this session (retrieved, not cited from memory)

1. **pVACtools (Hundal et al. 2020), Cancer Immunol Res 8(3):409–420** — WebSearch returned the verbatim
   indexed full-text: "predicts the strongest MHC binding peptides (8-11-mer for class I MHC and 13-25-mer for
   class II)"; supports "missense SNVs … along with frameshift variants and in-frame indels". Confirms the
   8–11 class I length band and the missense-substitution scope.
2. **ProGeo-neo (Li et al. 2020), BMC Med Genomics 13:52** — WebFetch of PMC7118832 returned verbatim:
   "Sequences corresponding to each of the coding missense mutations that would cause amino-acid substitutions
   were translated into a **21-mer amino acid fasta sequence, with 10 amino acids flanking the substituted
   amino acid on each side**"; "Class I MHC … binding peptide ligands of **8–11 amino acids**"; results are
   "peptides of 8–11 amino acids in length" and "**peptides without mutation sites were removed**". This is
   exactly the windowing rule: a 21-mer centred on the substitution (10 flanks each side) contains precisely
   every 8–11-mer that overlaps the central residue.
3. **NetMHCpan length basis (multiple, WebSearch)** — "Most MHCs prefer peptides of length 9"; NetMHCpan-4.0
   trained on 8–13-mer ligands; length preference varies by allele → motivates enumerating 8–11, not 9 only.
4. **Wells et al. 2020 (TESLA), Cell 183(3):818–834** — differential agretopic index = ratio of germline-peptide
   to matched mutated-peptide binding; confirms the mutant peptide must be paired with the wild-type peptide at
   the same coordinates (the agretope), which this unit produces.

### Formula check

Windowing rule (TestSpec §Neoantigen_Peptide_Generation.md §2.2, 1-based): a length-k window starting at `s`
spans the mutation iff `s ≤ p ≤ s+k−1`, i.e. `s ∈ [max(1, p−k+1), min(p, L−k+1)]`. The implementation's
0-based form (`firstStart = max(0, p0−k+1)`, `lastStart = min(p0, L−k)`, `p0 = p−1`) is the exact algebraic
translation. This is precisely the set of 8–11-mers contained in the ProGeo-neo 21-mer ±10-flank window.
**Confirmed against the ProGeo-neo verbatim definition.**

### Edge-case semantics (all sourced)

- **Terminal mutation** (within < k−1 of an end): fewer than k windows; only the windows that fit while spanning
  the mutation — ProGeo-neo "if possible" / pVACtools "centered … if possible". ✓
- **Non-substitution** (mutant == WT): not a missense variant → rejected. pVACtools generates peptides only from
  protein-altering variants. ✓
- **k outside 8–11**: outside the canonical class I band; the parameterised range respects it. ✓

### Independent cross-check (hand-computed from the definition, BEFORE reading code output)

Protein `MKTAYIAKQRSTVWLNDEFGH` (L=21): M1 K2 T3 A4 Y5 I6 A7 K8 Q9 R10 S11 T12 V13 W14 L15 N16 D17 E18 F19 G20 H21.

| Case | k | first/last start (1-based) | count | key window |
|------|---|----------------------------|-------|-----------|
| Y5C (p0=4) | 8 | 1 / 5 | 5 | s=1: mut `MKTACIAK`, WT `MKTAYIAK`, offset 4 |
| Y5C | 9 | 1 / 5 | 5 | last s=5: `CIAKQRSTV` |
| Y5C | 10 | 1 / 5 | 5 | — |
| Y5C | 11 | 1 / 5 | 5 | — |
| **Y5C total 8..11** | | | **20** | |
| V13A (p0=12, interior, ≥k−1 both ends) | 9 | 5 / 13 | **9** | exactly k windows |
| M1V (p0=0, N-term) | 9 | 1 / 1 | 1 | `VKTAYIAKQ` / `MKTAYIAKQ`, offset 0 |
| H21R (p0=20, C-term) | 8 | 14 / 14 | 1 | `WLNDEFGR` / `WLNDEFGH`, start 14, offset 7 |

Every value above was derived purely from the windowing definition and matches the Evidence dataset and the
test expectations. **No divergence.**

### Invariants

INV-1…INV-6 are all genuine consequences of the bound `s ∈ [max(1,p−k+1), min(p,L−k+1)]`: every peptide
is length k∈[min,max] and spans the mutation; mutant/WT have equal length and differ at exactly the offset
`p−s`; interior mutation gives exactly k windows (bound width `min(p,L−k+1) − max(1,p−k+1) + 1 = k` when
`p ≥ k` and `L−p ≥ k−1`). Verified algebraically.

### Stage A findings

No defect. The description, formula, edge-case semantics, and invariants are all confirmed verbatim against
retrieved primary sources (ProGeo-neo 21-mer ±10-flank; pVACtools 8–11mer; Wells agretope). **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:5123-5196` (`GenerateNeoantigenPeptides`),
record `NeoantigenPeptide` (5085), constants (5063/5070).

### Formula realised correctly?

Yes. Validation order: null → empty → `minLength < 1` → `maxLength < minLength` → position bounds → non-
substitution. Mutant protein materialised once by copying and replacing the residue at `mutationIndex`. The
loop ranges k = min..max, skips `k > proteinLength`, computes `firstStart = max(0, mutationIndex−k+1)`,
`lastStart = min(mutationIndex, proteinLength−k)`, emits `(k, start+1, mutantSubstr, wildTypeSubstr,
mutationIndex−start)`. This is the exact 0-based form of the sourced spanning bound. Output order is
k-ascending then start-ascending (INV-6) by loop structure. ✓

### Cross-verification table recomputed vs code

Ran the full suite (below). All 14 fixture cases pass; each independently-hand-computed value in the Stage A
table is reproduced by the code (counts 5/5/5/5 = 20 for Y5C; 9 interior for V13A; 1 each for M1V / H21R; exact
strings `MKTACIAK`/`MKTAYIAK`, `CIAKQRST`/`YIAKQRST`, `VKTAYIAKQ`/`MKTAYIAKQ`, `WLNDEFGR`/`WLNDEFGH`; offsets
4/0/0/7).

### Variant / delegate consistency

Single public method with two optional-parameter call shapes (default 8–11 range and explicit min/max). Both
shapes are exercised. `MhcClassIMinPeptideLength`/`MhcClassIMaxPeptideLength` are the defaults and match the
sourced 8/11 band.

### Test-quality audit (HARD gate)

- **Sourced, not code-echo:** every expected value is the hand-derived windowing value (exact strings, exact
  counts, exact offsets). M3 asserts mutant `'C'`/WT `'Y'` at the offset and `diffCount == 1`; M6/M7 assert exact
  terminal strings — none would pass against a deliberately-wrong implementation.
- **No green-washing:** counts and strings are exact equalities; no `Greater`/`AtLeast`/`Contains`/range where a
  value is known; no widened tolerance; no skipped/commented test.
- **Coverage:** all 7 MUST (M1, M2, M3, M3b, M4, M5, M6, M7), all 5 SHOULD (S1 null, S2 empty, S3 non-substitution,
  S4 position out of range ×2, S5 invalid length range ×2), both COULD (C1 protein shorter than all k → empty,
  C2 single-length subset). Both call shapes, all Stage-A branches (interior exact-k, N-terminal truncation,
  C-terminal truncation, k>L skip), and all documented error paths are exercised.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6661` (1 explicit benchmark not counted as a test);
  `dotnet build` 0 errors. The 4 NUnit2007 build warnings are pre-existing in an unrelated file
  (`ApproximateMatcher_EditDistance_Tests.cs`); no file touched this session emits a warning.

**Test-quality gate: PASS.**

### Stage B findings

No defect. Code faithfully realises the validated description; tests are exact, sourced, and fully cover the
Stage-A logic and edge cases. **Stage B: PASS.**

## Verdict & follow-ups

- **Stage A: PASS. Stage B: PASS. End-state: ✅ CLEAN.** No code, test, or spec change required.
- This session made **no** edits to code, tests, or specs; the unit was already correct and fully tested. Logged
  in the Findings Register under BY-DESIGN.
- Out-of-scope by design (documented, sourced): MHC binding/IC50 scoring (ONCO-MHC-001); indel/frameshift/fusion
  neopeptides; proteasomal cleavage / TAP; expression/VAF filtering.
