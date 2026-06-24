# Validation Report: PROTMOTIF-DOMAIN-001 — Protein Domain Identification

- **Validated:** 2026-06-24   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindDomains(string)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs:1320–1369`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_DomainPrediction_Tests.cs` (12 tests)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Scope reconciliation (important)

The TestSpec (last updated 2026-03-19) and Evidence doc described a **combined** "Domain Prediction
& Signal Peptide Prediction" unit, with `PredictSignalPeptide(sequence, maxLength)` returning a
tripartite n/h/c heuristic (`Score ∈ (0,1]`, `Probability = Score`, `NRegion/HRegion/CRegion`).

That signal-peptide model **no longer exists in the code.** Commit `8a6d3890`
("feat(ProteinMotif/Signal_Peptide_Prediction): von Heijne weight matrix") replaced it with the
von Heijne (1986) **position-specific weight-matrix / EMBOSS `sigcleave`** method — different
signature `PredictSignalPeptide(sequence, prokaryote, minWeight)` and different result record
(`CleavagePosition, Score [log-odds], SignalSequence, WindowSequence, IsLikelySignalPeptide`).
Signal-peptide validation was **split into a separate unit `PROTMOTIF-SP-001`**
(`tests/TestSpecs/PROTMOTIF-SP-001.md` + `ProteinMotifFinder_PredictSignalPeptide_Tests.cs`).

The **current DOMAIN-001 test file tests only `FindDomains`** and explicitly references SP-001 for
signal peptides. So the actual unit-under-test is correct; only the spec/evidence prose was stale.
**Fix applied:** retitled the TestSpec to domains-only and added a scope note marking the
signal-peptide material as superseded/moved to PROTMOTIF-SP-001 (doc-only; no source change). This
avoids "validating to a wrong spec".

---

## Stage A — Description (FindDomains)

### Sources opened & what they confirm
1. **PROSITE PS00028** (https://prosite.expasy.org/PS00028) — fetched. Pattern returned verbatim:
   `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`. Matches the code regex
   `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` exactly.
2. **PROSITE PS00017** (https://prosite.expasy.org/PS00017) — fetched. Pattern returned verbatim:
   `[AG]-x(4)-G-K-[ST]`. Matches the code regex `[AG].{4}GK[ST]` exactly.
3. **Pfam PF00018 (SH3) — InterPro** — confirms PF00018 is defined by a curated alignment + profile
   **HMM**, not a short consensus pattern. Same holds for PF00400 (WD40, β-propeller, ~40 aa/repeat)
   and PF00595 (PDZ, ~80–90 aa). The code's WD40/SH3/PDZ regexes are therefore deliberate
   **simplifications** of full Pfam HMMs, honestly declared in spec §1.4.1 ("Simplified pattern"),
   Evidence Design Decision 5, and inline test comments.

### Pattern check
- C2H2 zinc finger (PS00028) and Walker-A/P-loop (PS00017) are reproduced verbatim from PROSITE.
- WD40/SH3/PDZ carry the correct Pfam accessions (PF00400/PF00018/PF00595) but are short ad-hoc
  regexes, not the Pfam HMMs; they can produce false positives. This is transparently declared, not
  misadvertised (the API name `FindDomains` and XML doc say "signature patterns").

### Edge-case semantics
- Null/empty → empty enumerable (defined).
- No match (e.g. "AAAEEE") → empty (defined).
- Case insensitivity: `FindMotifByPattern` uses `ToUpperInvariant` + `RegexOptions.IgnoreCase`.
- Overlapping occurrences via lookahead wrapper `(?=(pattern))`; Start/End are 0-based inclusive of
  the captured substring; bounded 2 s regex timeout guards catastrophic backtracking.

### Independent cross-check (hand-computed, 0-based)
- Zinc finger `AAAACAACAAALEEEEEEEEHAAAHAAAA`: C@4, x(2)=5,6, C@7, x(3)=8,9,10, L@11, x(8)=12..19,
  H@20, x(3)=21,22,23, H@24 → **Start=4, End=24** ✓ (test M1).
- P-loop `AAAAGAAEAGKSAAAA`: [AG]@4, x(4)=5,6,7,8, G@9, K@10, [ST]@11 → **Start=4, End=11** ✓ (M2).

### Findings / divergences (NOTES)
- Honest-scope NOTE: only PS00028 and PS00017 are authoritative PROSITE patterns; WD40/SH3/PDZ are
  declared simplified regexes (not Pfam HMMs). → PASS-WITH-NOTES, not FAIL.

---

## Stage B — Implementation (FindDomains)

### Code path reviewed
- `FindDomains` (ProteinMotifFinder.cs:1320–1369) runs 5 regexes via `FindMotifByPattern`, wrapping
  each match into `ProteinDomain(Name, Accession, Start, End, Score, Description)`.
- `FindMotifByPattern` (:178–…) overlapping lookahead matching, IgnoreCase, timeout-guarded.

### Realised correctly?
- The two PROSITE patterns are literal copies of the validated strings; positions reported as
  0-based inclusive, Start ≤ End. Hand-computed M1/M2 positions match code and tests.

### Cross-verification table (recomputed vs code/tests)
| Case | Expected | Code/test result | Match |
|------|----------|------------------|-------|
| Zinc finger M1 | start 4, end 24, PF00096 | 4 / 24 / PF00096 | ✓ |
| P-loop M2 | start 4, end 11, PF00069 | 4 / 11 / PF00069 | ✓ |
| WD40 S1 | start 4, end 14, PF00400 | 4 / 14 / PF00400 | ✓ |
| SH3 S2 | start 4, end 16, PF00018 | 4 / 16 / PF00018 | ✓ |
| PDZ S3 | start 4, end 18, PF00595 | 4 / 18 / PF00595 | ✓ |
| No match S4 ("AAAEEE") | empty | empty | ✓ |
| Multi-domain C1 | ZF + Kinase both | both | ✓ |

### Variant/delegate consistency
- `FindDomains` reuses the same `FindMotifByPattern` engine as `FindCommonMotifs`/PROSITE search —
  consistent matching and Start/End/Score semantics across the class.

### Test quality audit
- 12 tests assert exact sourced values (positions, accessions, names), null/empty/no-match edges,
  case-insensitivity (incl. Score within 1e-10), metadata non-emptiness, Start ≤ End invariant, and
  multi-domain detection. No "no-throw"-only or tautological assertions. Signal-peptide tests
  correctly absent (covered by PROTMOTIF-SP-001).

### Findings / defects
- None for `FindDomains`. Build clean (0 warnings); 12 domain tests pass; full suite 18213/0.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — PS00028/PS00017 match PROSITE verbatim; WD40/SH3/PDZ are
  honestly-declared simplified regexes (Pfam uses HMMs). Spec/evidence prose for the (now removed)
  signal-peptide heuristic was stale → corrected to point at PROTMOTIF-SP-001.
- **Stage B:** PASS — `FindDomains` faithfully realises the validated patterns; all worked examples
  recomputed and matched.
- **End state:** CLEAN — no code defect; only a documentation-sync fix to the TestSpec.

**Code changed:** none. **Docs changed:** `tests/TestSpecs/PROTMOTIF-DOMAIN-001.md` (scope retitled
to domains-only; signal-peptide content marked superseded/moved to PROTMOTIF-SP-001).

**Sources:**
- PROSITE PS00028: https://prosite.expasy.org/PS00028
- PROSITE PS00017: https://prosite.expasy.org/PS00017
- Pfam/InterPro PF00018 (SH3, profile HMM): https://www.ebi.ac.uk/interpro/entry/pfam/PF00018
- Finn et al. (2014) Pfam: the protein families database, NAR 42:D222 — HMM-based families
- Krishna et al. (2003) NAR 31:532–550; Walker et al. (1982) EMBO J 1:945–951
