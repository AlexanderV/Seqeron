# Validation Report: VARIANT-INDEL-001 — Indel Detection

- **Validated:** 2026-06-15   **Area:** Variants
- **Canonical method(s):** `VariantCaller.FindInsertions(DnaSequence, DnaSequence)`, `VariantCaller.FindDeletions(DnaSequence, DnaSequence)`, `VariantCaller.FindIndels(DnaSequence, DnaSequence)` (delegate). Shared core: `CallVariants` → `SequenceAligner.GlobalAlign` → `CallVariantsFromAlignment`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **VCFv4.3 spec** (`https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex`, WebFetch 2026-06-15). Confirmed *verbatim*:
   - Single-base insertion: "A single base insertion of A after position 3 becomes REF=C, ALT=CA" ⇒ insertion ⇒ ALT longer than REF.
   - Single-base deletion: "A single base deletion of C at position 3 becomes REF=TC, ALT=T" ⇒ deletion ⇒ REF longer than ALT.
   - Padding-base rule for pure indels (left anchor base, reflected in POS; or base-after at contig position 1).
   - Microsatellite example `GTC → G,GTCT`: a 2-base deletion (TC) and a 1-base insertion (T), anchored at the preceding `G`.
   These directly support the variant-class definitions and INV-03/INV-04 directional-length invariants.
2. **minimal_representation `normalize.py`** (E. Minikel, raw GitHub, WebFetch 2026-06-15). Confirmed *verbatim* the two documented test tuples: `('7', 117199646, 'CTT', '-') → ('7', 117199644, 'ATCT', 'A')` and `('13', 32914438, 'T', '-') → ('13', 32914437, 'GT', 'G')`. Both corroborate deletion ⇒ len(REF) > len(ALT) and that pure (empty) indel alleles are padded with a left anchor — i.e. the directional-length invariant the repo asserts via the `"-"` sentinel form.
3. **Tan, Abecasis & Kang (2015)** — used (via the Evidence doc's verbatim quotes, re-confirmed against the DOI metadata) only to bound the *position* claim: an indel emitted from an arbitrary alignment column is not guaranteed to be in the left-aligned/parsimonious canonical form, so exact position is only deterministic where the optimal alignment is unique. This is the basis for ASM-02.

### Formula / model check
The model (pairwise global alignment; reference-gap column = insertion in query, query-gap column = deletion) is the standard VCF variant-class decomposition [VCFv4.3]. Insertion/deletion/SNP are distinct classes; an indel is length-changing while a SNP preserves length — matches the spec exactly.

### Edge-case semantics
- Identical sequences ⇒ no indels (INV-01): correct — no length-changing difference exists.
- Substitution-only (equal length) ⇒ no indels: correct per the distinct-class definition.
- Null inputs ⇒ `ArgumentNullException`: a reasonable, documented input-validation contract (not source-governed but standard).
- Multi-base indel ⇒ one per-base indel column per base (INV-05): consistent with the spec treating each extra/absent base as an event; this is a representation choice (per-base columns vs one merged allele), documented under ASM-01.

### Independent cross-check (numbers)
Re-derived a Needleman–Wunsch DP independently (match +1, mismatch −1, **linear** gap −1 — matching the repo's effective scoring; see Stage B note), enumerating **all** optimal alignments to test uniqueness:

| Case | ref / query | score | #optimal alignments | indels (type, refPos, REF, ALT) |
|------|-------------|------:|:---:|---|
| M1 | ATGCAT / ATGTCAT | 5 | 1 (unique) | INS @3 (`-`,`T`) |
| M2 | ATGTCAT / ATGCAT | 5 | 1 | DEL @3 (`T`,`-`) |
| M5 | ATGCAT / ATGTTCAT | 4 | 1 | INS @3, INS @3 |
| M6 | ATGTTCAT / ATGCAT | 4 | 1 | DEL @3, DEL @4 |
| C1 | ATGCAT / ATGTTTCAT | 3 | 1 | INS @3 ×3 |
| M3 | ATGCATGC / ATGTCATGG | 5 | 1 | INS @3 + SNP @7 |
| M4 | ATGTCATGC / ATGCATGG | 5 | 1 | DEL @3 + SNP @8 |
| S7 | ATGTCATGCAT / ATGCATGTCAT | 8 | 1 | DEL @3 + INS @8 |

Every position assertion in the spec/tests sits on a **provably unique** optimal alignment, so ASM-02's "exact position only on unique-alignment inputs" is honoured. Directional-length and per-base-column behaviour match the VCF spec and minimal_representation.

### Findings / divergences
None affecting correctness. ASM-01 (gap-sentinel in-memory form vs VCF padded form) and ASM-02 (no left-alignment/parsimony normalization) are explicitly documented; both are representation/position choices that do not change indel counts or types and are out of scope for the serialized-VCF rules (handled by `ToVcfLines`).

## Stage B — Implementation

### Code path reviewed
- `VariantCaller.cs:153-173` — `FindInsertions`/`FindDeletions`/`FindIndels` filter `CallVariants` by `VariantType`. Correct, trivial filters.
- `VariantCaller.cs:27-33,102-108` — `CallVariants` → `CallVariantsCore` → `SequenceAligner.GlobalAlign` + `CallVariantsFromAlignment`. Null checks at `:29-30`.
- `VariantCaller.cs:41-100` — `CallVariantsFromAlignment`: reference-gap ⇒ Insertion (`ReferenceAllele="-"`, ALT = query base), query-gap ⇒ Deletion (REF = ref base, `AlternateAllele="-"`), differing bases ⇒ SNP. `refPos` advances on reference-consuming columns (deletion/SNP/match), `queryPos` on query-consuming columns. Positions are 0-based (ASM-03).
- `SequenceAligner.cs:220-273` (`GlobalAlignCore`) + `:468-534` (`Traceback`): NW DP. **Note:** although `SimpleDna` declares `GapOpen=-2`, the DP uses only `GapExtend` (−1) for every gap step (`:255-257`), so the effective gap model is **linear −1** — exactly what the spec/doc state ("linear gap −1"). Verified the description matches the code, not the unused field.

### Formula realised correctly?
Yes. The independent NW (linear gap −1, repo traceback tie-order diag→up→left) reproduces the repo's exact aligned strings and variant decompositions for all 8 cases (table above). The classification logic matches the VCF variant-class definitions.

### Cross-verification table recomputed vs code
The full unfiltered suite (6482 tests) passes, including all 16 VARIANT-INDEL-001 tests; the independently computed counts/positions/alleles/types match the asserted values for M1–M8, S1–S7, C1.

### Variant/delegate consistency
`FindIndels` = union of Insertion+Deletion filters over the same `CallVariants` result as the individual finders — consistent by construction. Verified on S7 (one DEL + one INS).

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** counts, 0-based positions, exact alleles, and types are all traceable to the VCF spec single-base/microsatellite examples + minimal_representation directional length, and to the independently-enumerated unique optimal alignments. A deliberately-wrong implementation (e.g. swapping INS/DEL allele assignment, or merging multi-base indels) would fail M1/M2/M5/M6/C1.
- **No green-washing:** exact `Has.Count.EqualTo`, exact allele equality, exact `Position` equality where alignment is unique; no widened tolerances, no skips. The `All(...)` filter predicates in M3/M4/S7 are paired with `Is.Not.Empty`/`Any(...)`, so they are not vacuous; exact counts for those scenarios are independently locked by M1/M2/M5/M6.
- **Coverage:** all three public methods exercised; null reference + null query on both finders (S1–S4); identical-input (M7/M8); substitution-only (S5/S6); single, multi-base, and k-block indels; mixed indel+SNP filtering; insertion+deletion union. All Stage-A branches/edge cases covered.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6482`; `dotnet build` 0 errors. The 4 build warnings are pre-existing and in unrelated test files (none in `VariantCaller_FindIndels_Tests.cs`).

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS.** Description, formulas, invariants and edge cases match authoritative sources retrieved this session.
- **Stage B: PASS.** Code faithfully realises the validated model; all expected values independently confirmed; tests are real, exhaustive of the Stage-A surface, and honestly green.
- **End-state: ✅ CLEAN.** No defect found; no code or test changes required.
- Minor (non-defect) observation, already documented in the algorithm doc/spec: `SimpleDna.GapOpen` is unused by the global-alignment DP (effective gap model is linear −1). This is consistent with the description and does not affect indel detection.
