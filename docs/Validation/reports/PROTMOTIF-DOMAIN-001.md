# Validation Report: PROTMOTIF-DOMAIN-001 — Protein Domain Identification

## Update 2026-06-25 — Plan7 profile-HMM engine + bundled Pfam SH3/PDZ/WD40 (limitation fix); Status stays ☐

The data-blocked residual ("SH3/PDZ and any full Pfam HMM have no deterministic pattern") is
now **addressed for SH3/PDZ/WD40** with a faithful, opt-in profile-HMM scorer. The exact
PROSITE-pattern `FindDomains` path and its defaults are **unchanged**.

- **Profiles bundled (CC0 / public domain):** HMMER3/f ASCII profiles `PF00018.35` (SH3_1, LENG 48),
  `PF00595.30` (PDZ, LENG 81), `PF00400.39` (WD40, LENG 39), retrieved 2026-06-25 from the EMBL-EBI
  InterPro web API (`.../interpro/wwwapi/entry/pfam/<ACC>/?annotation=hmm`) and embedded verbatim
  under `Seqeron.Genomics.Analysis/Resources/`. Licence verbatim: "Pfam is freely available under
  the Creative Commons Zero ('CC0') licence" (InterPro/Pfam docs).
- **Engine:** new `Plan7ProfileHmm` (HMMER3/f parser + glocal Viterbi/Forward log-odds DP) per the
  HMMER User's Guide v3.4 file format ("negative natural-log probabilities; '*' = zero") and the
  Durbin et al. (1998) §5.4 / Eddy (2011) Plan7 recurrences. New opt-in methods
  `ProteinMotifFinder.FindDomainsByHmm(seq, minBitScore)` and `ScoreDomainHmm(seq, accession)`.
- **Exact DP verification:** a hand-built 2-match-state HMM scores "AC" via B→M1→M2→E to
  **0.5187937934151676 nats** = ln(0.7/0.6)+ln 0.9+ln(0.8/0.4)+ln 0.8; engine matches to 1e-9 (test H1).
- **Real detection (ranking):** SRC_HUMAN SH3 core ≈ +60 bits vs PF00018; DLG4_HUMAN PDZ1 (res 61–151)
  ≈ +83 bits vs PF00595; GBB1_HUMAN ≈ +36 bits vs PF00400; a low-complexity negative is strongly
  negative for all three and reported by none. Cross-domain specificity confirmed (SH3 core vs
  PF00400 ≈ −25 bits). 17 new tests (H1–H12), full suite green.
- **Honest residual:** exact `hmmsearch` bit-score / E-value parity (MSV/bias filters, null2,
  Gumbel/exponential calibration) and the full Pfam library beyond these 3 domains are out of scope.
- Status stays **☐** in the root registry (independent re-validation; this is not a ☑ self-claim).

---

- **Validated:** 2026-06-24   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindDomains(string)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs` (Domain Finding region)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_DomainPrediction_Tests.cs` (14 tests)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Update 2026-06-24 — exact PROSITE patterns (limitation fix); Status reset to ☐

The prior PASS-WITH-NOTES rested on WD40/SH3/PDZ being honestly-declared *simplified* ad-hoc
regexes. This limitation fix removes that residual where a deterministic PROSITE pattern exists:

- **WD40 → EXACT PROSITE PATTERN PS00678** (WD_REPEATS_1), retrieved verbatim 2026-06-24 from
  https://prosite.expasy.org/PS00678 and https://prosite.expasy.org/PDOC00574:
  `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]`.
  Translated by the official ScanProsite syntax rules
  (https://prosite.expasy.org/scanprosite/scanprosite_doc.html) to
  `[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]`
  (14 elements, fixed 15 residues). The prior ad-hoc `[LIVMFYWC]-x(5,12)-[WF]-D` regex was replaced.
  Validated against the real WD40 β-propeller GBB1_HUMAN (UniProt P62873): matches at 0-based
  starts 69, 156, 284.
- **SH3 (PROSITE PS50002) and PDZ (PROSITE PS50106) are PROFILES (weight matrices), not patterns**
  (confirmed via PDOC50002 / PDOC50106 + EBI InterPro). No deterministic pattern exists, so the
  previously shipped *unsourced* ad-hoc SH3/PDZ regexes were **removed** rather than fabricating a
  signature. `FindDomains` no longer reports SH3/PDZ — an **honest residual**. Pfam HMM profiles
  (PF00018, PF00595, the full PF00400 family) are trained models and are not bundled.
- **Unchanged exact patterns:** zinc finger C2H2 **PS00028** and Walker-A / P-loop **PS00017**.
- **Tests:** now 14 (was 12) — added M7 (PS00678 real-segment match, Start=4/End=18/15 aa),
  M8 (near-miss: invariant Trp removed → no match), M9 (GBB1_HUMAN hits {69,156,284}),
  M10 (verbatim PS00678 string through the PROSITE→regex translator vs hand-traced 15-mer),
  S5 (Src SH3 core is NOT reported as SH3/PDZ). Removed the old SH3/PDZ positive tests.
- Full suite green; Status reset ☑→☐ in the root registry for independent re-validation.

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
