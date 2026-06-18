# Validation Report: VARIANT-CALL-001 — Variant Detection

- **Validated:** 2026-06-15   **Area:** Variants
- **Canonical method(s):** `VariantCaller.CallVariants(DnaSequence, DnaSequence)`,
  `VariantCaller.CallVariantsFromAlignment(string, string)`,
  `VariantCaller.ClassifyMutation(Variant)`, `VariantCaller.CalculateTiTvRatio(IEnumerable<Variant>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo artifacts)

1. **VCFv4.x specification (samtools/hts-specs).** The v4.3 PDF was fetched and saved this session;
   its FlateDecode text streams were decompressed locally and the field-definition + example text
   read directly. Cross-checked against the v4.1/v4.2 PDFs and a web search. Extracted verbatim:
   - **POS (1-based):** "POS | position: The reference position, with the 1st base having position 1."
   - **REF alphabet / case:** "REF | reference base(s): Each base must be one of A,C,G,T,N (case
     insensitive). Multiple bases are permitted. The value in the POS field refers to the position
     of the first base in the String."
   - **Indel padding base:** "...the REF and ALT Strings must include the base before the event
     (which must be reflected in the POS field), unless the event occurs at position 1 on the
     contig in which case it must include the base after the event; this padding base is not
     required..." — i.e. padding/1-based POS govern *serialized* VCF, not the in-memory model.
   - **Microsatellite example (§1.1 example record `20 1234567 microsat1 GTC G,GTCT ...`):**
     "a microsatellite with two alternative alleles, one a **deletion of 2 bases (TC)**, and the
     other an **insertion of one base (T)**." Confirms the M4/M5/M6 biology.
   - **Simple SNP example:** the record `20 14370 rs6054257 G A ...` is described as "a good simple
     SNP" — single-base REF/ALT, no padding. Confirms M2/M3 shape.
2. **Wikipedia — Transversion** (fetched). Verbatim: a transversion is "a point mutation in DNA in
   which a single (two ring) purine (A or G) is changed for a (one ring) pyrimidine (T or C), or
   vice versa." The four transversions listed: **A↔T, A↔C, G↔T, G↔C**. Transitions: **A↔G** and
   **C↔T**. Per base: 1 possible transition, 2 possible transversions.
3. **Wikipedia — Transition (genetics)** (fetched). Confirms **Ti/Tv = (#transitions)/(#transversions)**,
   and the documented transition bias ("transitions occur more often in genomes…", "~two out of
   three SNPs are transitions"), cited to Collins & Jukes (1994) / Ebersberger (2002).

### Formula / definition check

- **Variant classes** (SNP / insertion / deletion) — match VCF scope (Danecek 2011; VCF spec).
- **Transition iff {ref,alt}⊆{A,G} or ⊆{C,T}; transversion otherwise** — matches the verbatim
  purine/pyrimidine definitions above. INV-05 is exactly the sourced definition.
- **Ti/Tv numerator/denominator orientation** — matches Wikipedia (transitions over transversions).
- **Coordinate conventions** — the spec's 1-based POS + padding base apply to serialized VCF only;
  the in-memory `Variant` uses 0-based `Position` + `"-"` gap sentinel (ASM-01/ASM-02). This is an
  internal model choice not governed by a source; it is internally consistent and is the
  established sibling contract. Documented, not a defect.

### Edge-case semantics

- Identical sequences → 0 variants (variant = difference from reference): sourced.
- Empty alignment → empty; unequal aligned lengths → error: contract, sourced as structural.
- Case-insensitive classification: sourced (REF/ALT "case insensitive").
- Ti/Tv with zero transversions: mathematically undefined; repo maps to 0 (ASM-03). No source
  mandates a sentinel; documented as contract, not a correctness claim.

### Independent cross-check (hand computation, traced against the validated definitions)

- `GTC` vs `G` (M6 analogue): deletion of 2 bases TC — matches spec text exactly.
- `ATGC`/`ATTC`: 1 SNP G→T at 0-based pos 2 (= transversion, since G is purine, T pyrimidine).
- Ti/Tv on {A→G (Ti), A→C (Tv)} = 1/1 = 1.0; on {A→G, C→T (both Ti), A→C (Tv)} = 2/1 = 2.0;
  on {A→G, C→T} (0 Tv) = undefined → 0 by contract.

### Findings / divergences

None. Description is biologically and mathematically correct. The only divergences from VCF
(0-based Position, `"-"` sentinel, no left-align/parsimony normalization) are explicit, sourced
as out-of-scope for the in-memory detection contract, and recorded as ASM-01/02/03.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs`:
- `CallVariants` (L27–33): null-checks ref/query, delegates to `CallVariantsCore`.
- `CallVariantsCore` (L102–108): global-aligns then calls `CallVariantsFromAlignment`.
- `CallVariantsFromAlignment` (L41–100): column scan; ref-gap→Insertion, query-gap→Deletion,
  mismatch→SNP; correct 0-based ref/query position bookkeeping; empty→empty, unequal length→throw.
- `ClassifyMutation` (L184–198): non-SNP→Other; `ToUpperInvariant` both bases;
  `refPurine == altPurine ? Transition : Transversion` — exactly INV-05.
- `CalculateTiTvRatio` (L203–220): null-check; counts via `ClassifyMutation` over SNPs only;
  `transversions > 0 ? transitions/transversions : 0` — ASM-03 contract.

### Formula realised correctly?

Yes. Hand-traced each canonical case:
- M2/M3 `ATGC`/`ATTC` → 1 SNP, Position 2, REF G, ALT T. ✓
- M4 `AT-GC`/`ATTGC` → Insertion at refPos 2, queryPos 2, REF "-", ALT T. ✓
- M5 `ATTGC`/`AT-GC` → Deletion at refPos 2, REF T, ALT "-". ✓
- M6 `GTCAA`/`G--AA` → two Deletions at refPos 1 (T) and 2 (C). ✓ (matches spec "deletion of 2 bases TC")
- S5 `AAAA`/`TGTA` → SNPs at 0,1,2 with ALT T,G,T. ✓
- Classification + Ti/Tv values 1.0 / 2.0 / 0 reproduced by hand. ✓

### Cross-verification table recomputed vs code

| Input | Sourced expectation | Code result | Match |
|-------|--------------------|-------------|-------|
| `ATGC`/`ATTC` | 1 SNP pos2 G→T | same | ✓ |
| `GTCAA`/`G--AA` | 2 del (T@1, C@2) | same | ✓ |
| A→G | Transition | Transition | ✓ |
| C→T, G→A, T→C | Transition | Transition | ✓ |
| A→C/A→T/G→C/G→T | Transversion | Transversion | ✓ |
| a→g (lowercase) | Transition | Transition | ✓ |
| Deletion | Other | Other | ✓ |
| {A→G,A→C} | 1.0 | 1.0 | ✓ |
| {A→G,C→T,A→C} | 2.0 | 2.0 | ✓ |
| {A→G,C→T} | 0 (contract) | 0 | ✓ |

### Variant/delegate consistency

`CalculateTiTvRatio` reuses `ClassifyMutation`; both filter `Type == SNP` consistently. `CallVariants`
and `CallVariantsFromAlignment` share the same column-scan core. Out-of-scope siblings
(`FindSnps*`, `FindIndels`, `PredictEffect`, `ToVcfLines`, `CalculateStatistics`) are owned by other
units and were not re-validated here.

### Test quality audit (canonical file `VariantCaller_CallVariants_Tests.cs`)

- **Sourced, not code-echoes:** exact positions (0,1,2), exact alleles, exact Ti/Tv values (1.0/2.0/0)
  — each traceable to the VCF spec text / transversion definition / Ti/Tv definition retrieved this
  session, not to whatever the code happens to return. A deliberately-wrong implementation (e.g.
  swapping Ti/Tv orientation, or off-by-one position) would fail these.
- **No green-washing:** no `GreaterThan`/`Any`/range substituted where an exact value is known; tight
  `Within(1e-10)` tolerances on ratios; no skips/ignores.
- **Coverage:** all four canonical methods + all Stage-A branches: identical→empty, SNP via both
  entry points, insertion column, deletion column, two-base deletion, multiple SNPs, all four
  transition pairs, all four transversions, case-insensitive, non-SNP→Other, Ti/Tv 1:1 / 2:1 /
  zero-denominator, and error/edge cases (null ref, null query, null variants, unequal length,
  empty input) + a position-bounds/SNP-distinctness property test (C1).
- **Minor strengthening opportunity (not a defect):** M11 exercises case-insensitivity only on a
  lowercase transition (`a→g`). The `ToUpperInvariant` path is symmetric and fully exercised; a
  lowercase transversion case would add marginal coverage. Left as-is — does not affect the verdict.

### Findings / defects

None. `dotnet build` = 0 errors (the 4 NUnit2007 warnings are in an unrelated file,
`ApproximateMatcher_EditDistance_Tests.cs`, pre-existing and out of scope). Full unfiltered suite:
**Failed: 0, Passed: 6482, Skipped: 0**.

## Verdict & follow-ups

- **Stage A: PASS** — description matches primary sources (VCF spec, transition/transversion
  definitions, Ti/Tv definition), all retrieved this session.
- **Stage B: PASS** — code faithfully realises the validated description; cross-checks reproduce all
  sourced values; tests assert exact sourced values and cover every branch and edge case.
- **End-state: CLEAN** — no defect; no code or test change required.
- **Test-quality gate: PASS** — sourced (not echoed), no green-washing, full branch/edge coverage,
  honest green (full suite Failed: 0, build warning-free for changed files; none changed).
- No defects logged in FINDINGS_REGISTER.
