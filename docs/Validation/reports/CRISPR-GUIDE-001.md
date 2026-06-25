# Validation Report: CRISPR-GUIDE-001 — CRISPR guide RNA on-target efficacy scoring

- **Validated:** 2026-06-24 (independent re-confirmation of the Doench-2014 "Rule Set 1" fix, commit 129c2ca, plus the Doench-2016 Rule Set 2 / Azimuth addition, commit 57730b9)   **Area:** MolTools
- **Canonical method(s):**
  - `CrisprDesigner.CalculateOnTargetDoench2014(string context30Mer)` — Rule Set 1 (linear logistic model). `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs:521`; constants/table at lines 473–502.
  - `CrisprDesigner.CalculateOnTargetRuleSet2(...)` — Rule Set 2 / Azimuth (GBRT). Engine `AzimuthRuleSet2.cs`; embedded `Resources/azimuth_rs2_{nopos,full}.bin`; extractor `scripts/azimuth/extract_azimuth_model.py`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN — no defect found; no code or test change required.

This session independently re-confirms the prior detailed validation. The Doench 2014 model was
re-grounded from scratch against the freshly re-downloaded CRISPOR reference; the Rule Set 2 /
Azimuth pipeline was re-confirmed at the level of provenance, oracle externality, and green tests.

## Stage A — Description

### Sources opened this session (independent of the repo)
- **Doench et al. 2014**, "Rational design of highly active sgRNAs…", *Nat Biotechnol* 32:1262–1267, PMID 25184501 — Rule Set 1 is a logistic-regression linear model over a fixed 30-nt context.
- **Reference `doenchScore.py`** (CRISPOR, Haeussler et al. 2016, *Genome Biol* 17:148). **Re-downloaded the raw 54-line file this session** to `/tmp/doenchScore.py` from
  `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/doenchScore.py` and read the
  constants and `params` list **directly** — did not trust the repo's copied table.
- **Microsoft Research Azimuth** (BSD-3-Clause) — provenance of the Rule Set 2 trained pickles, as recorded in the repo's extractor header.

### Rule Set 1 — model layout (confirmed against the reference, read directly)
- **30-mer window:** `[4 nt 5′] + [20 nt protospacer] + [3 nt PAM] + [3 nt 3′] = 30`. GC term over `seq[4:24]`; features indexed by raw 0-based offset (`seq[pos:pos+len]`). Confirmed.
- **Intercept `0.59763615`** (reference line 30) — matches repo `DoenchIntercept` exactly.
- **`gcLow = -0.2026259`, `gcHigh = -0.1665878`** (reference lines 31–32) — match repo `DoenchGcLow`/`DoenchGcHigh` exactly.
- **GC term:** `score += abs(10 - gcCount) * gcWeight`, `gcWeight = gcLow if gcCount<=10 else gcHigh` (reference lines 38–43). Boundary `<= 10 → gcLow`. Confirmed.
- **Output:** `1/(1+exp(-score))` → (0,1); repo scales ×100. Confirmed.

### Coefficient table — full independent diff (this session)
Parsed both tables programmatically: the reference `params` via `exec` of the downloaded file, and
the repo's `DoenchParams` via regex over `CrisprDesigner.cs`. Result:

```
ref count: 70   repo count: 70
ordered exact match: True
in ref not repo: ∅
in repo not ref: ∅
```

Byte-for-byte identical, including the reference's intentional quirks where `(24,'AG'/'CG'/'TG')`
dinucleotides reuse the `(24,'A'/'C'/'T')` weights and `(26,'GT')` = `(27,'T')` = 0.11787758.

### Independent cross-checks (worked examples reproduced from scratch)
Re-implemented the reference formula in Python from the constants and ran it using the **repo's
parsed table** (so the check ties the repo table to the published expected outputs):

| 30-mer | Computed ×1 | Reference published | Δ |
|---|---|---|---|
| `TATAGCTGCGATCTGAGGTAGGGAGGGACC` | 0.7130893469547 | 0.713089368437 | 2.1e-08 |
| `TCCGCACCTGTCACGGTCGGGGCTTGGCGC` | 0.0189838431933 | 0.0189838463593 | 3.2e-09 |
| `AAAAAAAAAAAAAAAAAAAAAAAAAGGAAA` (M-003) ×100 | 4.4338168085440 | test oracle 4.4338168085 | <1e-9 |

The ~1e-8 deltas are float-print precision in the reference's quoted literals; the model reproduces exactly.

### Edge-case semantics
- **Wrong length / null / empty / non-ACGT:** repo throws (`ArgumentException` / `ArgumentNullException`) — stricter, defensible typed contract over the reference's guard-free function. Sourced as the model's documented 30-mer requirement.
- **Lowercase:** upper-cased before scoring → identical value.
- **NGG PAM guard (offsets 25–26 == GG):** repo enforces SpCas9 specificity; the bare reference does not. A documented input guard, not a scoring-model change. Does not alter scores for valid inputs.

### Rule Set 2 / Azimuth — description (re-confirmed, not re-derived this session)
Rule Set 2 is a trained scikit-learn GBRT, not a coefficient table; it cannot be reproduced from
published numbers, only from the shipped pickles. The repo reproduces it sklearn-free. This session
re-confirmed: (a) both `azimuth_rs2_{nopos,full}.bin` and `scripts/azimuth/extract_azimuth_model.py`
are present; (b) the oracle CSVs (`scripts/azimuth/oracle/{nopos,full}_oracle.csv`, 947 rows) carry
both `ref_score` (the verified reference prediction) **and** `upstream` (the Microsoft Azimuth fixture)
with an `agrees` flag — confirming the test oracles are externally derived, not read off the C# code;
(c) the extractor header documents the full provenance and the known upstream-fixture drift (~38% rows),
consistent with the prior multi-session report. The detailed sklearn-free recovery, 1e-13 featurizer
match, real-Biopython Tm cross-check, CPython-2.7 column order, and bit-identical scikit-learn traversal
were established in the prior report and are accepted as still valid (model files and tests unchanged).

**Stage A verdict: PASS** — Rule Set 1 layout, intercept, GC term form/boundary, output transform, and
the full 70-entry coefficient table all match the primary reference exactly; three worked examples
reproduced independently. Rule Set 2 provenance and oracle externality re-confirmed.

## Stage B — Implementation

- **Code path (Rule Set 1):** `CrisprDesigner.cs:521–568`. `score = intercept`; GC over `seq[4..24)` with `Math.Abs(10-gc)*gcWeight`, `gc<=10→gcLow`; feature loop adds `weight` when `string.CompareOrdinal(seq,pos,modelSeq,0,len)==0` (≡ `seq[pos:pos+len]==modelSeq`); `1/(1+e^-score)*100`. Line-for-line faithful to the reference.
- **Code path (Rule Set 2):** `AzimuthRuleSet2.cs` reader + featurizer + GBRT traversal; thin wrappers in `CrisprDesigner.cs`. Additive; no existing signature/test changed.
- **Cross-verification vs the actual C#:** test classes `CrisprDesigner_Doench2014_Tests` (9 tests), `CrisprDesigner_RuleSet2_Tests` (13 tests), and `CrisprDesigner_GuideRNA_Tests` (heuristic evaluator) run **Passed: 54, Failed: 0**. M-001/M-002 lock the reference worked examples ×100 (71.3089368437 / 1.89838463593, tol 1e-4); M-003 locks the independently-recomputed 4.4338168085 — all three reproduced by hand this session.
- **Variant/delegate consistency:** Doench 2014 is a single static method; Rule Set 2 nopos/full wrappers delegate to one engine. Heuristic `EvaluateGuideRna`/`DesignGuideRnas` (TestSpec scoring model) unchanged and green.
- **Test-quality audit (HARD gate):** PASS, no green-washing. M-001/M-002/M-003 are exact externally-sourced cross-checks; per the prior report and the model's sensitivity, a wrong intercept (0.59763615→0.69763615) shifts the score ~2 points ≫ tol and fails them. Edge tests cover wrong length (29/31/short), null/empty, non-ACGT (`N`), lowercase equivalence, range/ordering, and non-NGG PAM rejection. Rule Set 2 tests assert against the verified reference and upstream-agreeing subset (counts locked), never against C# output.
- **Numerical robustness:** sigmoid bounded; GBRT is finite-tree sum; no overflow/precision concern over the coefficient/threshold ranges.

**Stage B verdict: PASS** — the code faithfully realises both validated reference models; tests assert exact externally-sourced values.

## Verdict & follow-ups
- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.** Re-downloaded the CRISPOR reference and re-confirmed the intercept (0.59763615), gcLow (−0.2026259), gcHigh (−0.1665878), the 30-mer layout (4+20+3+3), the GC-term boundary, and the full 70-entry table as byte-identical; reproduced three worked-example scores to ≤2e-8. Rule Set 2 / Azimuth provenance and oracle externality re-confirmed.
- **No code or test change required this session.**
- **Tests run:** `CrisprDesigner_Doench2014_Tests | CrisprDesigner_RuleSet2_Tests | CrisprDesigner_GuideRNA_Tests` → 54 passed, 0 failed (build: 0 warnings, 0 errors). Full unfiltered suite not re-run because no code changed.
