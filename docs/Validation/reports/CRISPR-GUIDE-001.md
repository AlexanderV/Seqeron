# Validation Report: CRISPR-GUIDE-001 — On-target guide efficacy (Doench 2014 "Rule Set 1")

- **Validated:** 2026-06-17 (Phase-3 independent re-validation of the Doench-2014 enhancement, commit 129c2ca)   **Area:** MolTools
- **Canonical method(s):** `CrisprDesigner.CalculateOnTargetDoench2014(string context30Mer)` (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs:521`); coefficient table `DoenchParams` + `DoenchIntercept`/`DoenchGcLow`/`DoenchGcHigh` (lines 473–502).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN — no defect found; no code/test change required.

This re-validation is scoped to the **Doench 2014 "Rule Set 1" on-target model** added in commit
129c2ca. The pre-existing honest-heuristic `CalculateGuideScore`/`DesignGuides` paths and the
MIT/Hsu off-target additions (CRISPR-OFF-001) were not the subject of this session.

## Stage A — Description

### Sources opened this session (independent of the repo)
- **Doench, Hartenian, Graham, et al.**, "Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation." *Nat Biotechnol* 32:1262–1267 (2014), PMID 25184501 — the model is a logistic-regression linear model over a fixed 30-nt context.
- **Reference implementation `doenchScore.py`** (CRISPOR, Haeussler et al. 2016, *Genome Biol* 17:148; github.com/maximilianh/crisporWebsite). **Re-downloaded the raw file this session** to `/tmp/doenchScore.py` (54 lines) and read it directly — did NOT rely on the repo's copied constants.

### Model layout (confirmed against the reference)
- **30-mer window:** `[4 nt 5′ context] + [20 nt protospacer] + [3 nt PAM] + [3 nt 3′ context] = 30`. The reference computes the GC term over `seq[4:24]` (= the 20-nt protospacer) and indexes all features by raw 0-based offset into the 30-mer (`subSeq = seq[pos:pos+len(modelSeq)]`). Confirmed.
- **Intercept:** `0.59763615` (reference line 30). ✔
- **GC-count term:** `score += abs(10 - gcCount) * gcWeight` with `gcWeight = gcLow (-0.2026259)` if `gcCount <= 10` else `gcHigh (-0.1665878)` (reference lines 31–32, 38–43). The form is `|GC − 10|` (absolute value) and the boundary is `<= 10 → gcLow`. Confirmed: the repo uses `Math.Abs(10 - gcCount)` and `gcCount <= 10 ? DoenchGcLow : DoenchGcHigh` (CrisprDesigner.cs:556–557) — identical.
- **Output:** `1.0/(1.0+exp(-score))` → probability in (0,1). The repo multiplies by 100 to expose a 0–100 scale. Confirmed and documented.

### Coefficient table — full independent diff
The repo's 70-entry `DoenchParams` table was parsed from the C# source and compared **tuple-by-tuple, in order** against the reference `params` list read from `/tmp/doenchScore.py`:

```
ref count: 70   repo count: 70
EXACT MATCH (ordered): True
In ref not repo: ∅
In repo not ref: ∅
```

There is **no transcription error**: every `(pos, seq, weight)` is byte-for-byte identical, including the reference's own intentional quirk where the dinucleotide tuples `(24,'AG'/'CG'/'TG')` reuse the single-nucleotide weights `(24,'A'/'C'/'T')` and `(26,'GT',0.11787758)` equals `(27,'T',0.11787758)`. Spot-checked weights explicitly against the reference file: `(22,'T',-0.8770074)`, `(23,'C',-0.8762358)`, `(11,'GG',-1.5169074)`, `(20,'GG',-0.7822076)`, `(21,'TC',-1.029693)`, `(28,'GG',-0.69774)`, `(29,'G',0.38634258)`, `(6,'C',-0.7411813)`, `(17,'G',-0.6780964)` — all match.

### Independent cross-check (worked examples reproduced from scratch)
Implemented the reference formula in Python from the coefficients (no trust of the repo's test literals) and ran the reference's own two self-test 30-mers:

| 30-mer | My computation ×1 | Reference published value | |
|---|---|---|---|
| `TATAGCTGCGATCTGAGGTAGGGAGGGACC` | 0.713089346955 | 0.713089368437 | Δ 2.2e-08 |
| `TCCGCACCTGTCACGGTCGGGGCTTGGCGC` | 0.018983843193 | 0.0189838463593 | Δ 3.2e-09 |

The tiny deltas are float-print precision in the reference's quoted literal; the model reproduces exactly.

### Edge-case semantics
- **Wrong length:** the reference assumes a 30-mer (no guard). The repo adds an explicit `length != 30 → ArgumentException` guard — a stricter, defensible contract for a typed API. Sourced as the model's documented input requirement.
- **Lowercase:** the repo upper-cases before scoring (`ToUpperInvariant`), so lowercase is accepted and gives the identical value. Reasonable extension.
- **Non-ACGT:** repo throws `ArgumentException`. The reference would silently produce a degenerate score; throwing is the safer typed contract.
- **0–100 vs 0–1:** the repo deliberately scales reference×100; documented in code and tests.

### Stage A divergences (all documented, none a defect)
- The repo adds an **NGG PAM check** (offsets 25–26 must be `GG`) that the reference does NOT enforce. This is a defensible SpCas9-specificity guard and is documented; it does not alter scoring for valid inputs. (Note: it restricts inputs the bare reference would still score; it is a guard, not a model change.)

**Stage A verdict: PASS** — layout, intercept, GC term form/boundary, output transform, and the full coefficient table all match the primary reference exactly; worked examples reproduced independently.

## Stage B — Implementation

- **Code path reviewed:** `CrisprDesigner.CalculateOnTargetDoench2014` (CrisprDesigner.cs:521–568) and `DoenchParams`/constants (473–502).
- **Formula realised correctly:** `score = intercept`; GC over `seq[4..24)` with `Math.Abs(10-gc)*gcWeight`, `gc<=10→gcLow`; feature loop adds `weight` when `string.CompareOrdinal(seq, pos, modelSeq, 0, len)==0` (= `seq[pos:pos+len]==modelSeq`); `1/(1+e^-score)*100`. Matches the reference line-for-line.
- **Cross-verification vs the actual C#:** the 11 Doench tests pass (`--filter ~Doench2014` → Passed 11, Failed 0). M-001/M-002 lock the reference worked examples ×100 (71.3089368437 / 1.89838463593, tol 1e-4). M-003 all-A-with-PAM expected 4.4338168085 — independently reproduced from the coefficients in Python: `4.4338168085440035`. ✔
- **Test-quality audit (HARD gate):** PASS, no green-washing.
  - M-001/M-002 are exact sourced cross-checks, not tautologies. Mutation test: perturbing the intercept by +0.001 shifts the score by ~0.020 (≫ tol 1e-4) — a wrong intercept/weight **would** fail these tests. (Implementer's own commit notes a 0.59763615→0.69763615 flip kills both; reconfirmed the sensitivity here.)
  - Edge tests genuinely cover wrong length (29/31/short), null/empty (`ArgumentNullException`), non-ACGT (`N`), lowercase equivalence, range [0,100] + ordering, and non-NGG PAM rejection.
- **Numerical robustness:** sigmoid is bounded, no overflow/precision concern over the coefficient range.

**Stage B verdict: PASS** — the code faithfully realises the validated reference model; tests assert exact externally-sourced values and a wrong constant fails.

## Verdict & follow-ups
- **End-state: ✅ CLEAN.** No transcription error in the 70-coefficient table, intercept, or GC terms; both reference worked examples reproduced independently; tests are real (exact sourced values, mutation-sensitive). No code or test change required this session.
- Residual (faithfully scoped, not a defect): Doench Rule Set 2 / Azimuth (GBT model, no coefficient table) and CFD (binary pickle) remain intentionally unimplemented — documented in FINDINGS_REGISTER C7 and LIMITATIONS.
- **Full unfiltered suite:** 6772 passed, 0 failed; build 0 errors (4 pre-existing NUnit2007 warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).
