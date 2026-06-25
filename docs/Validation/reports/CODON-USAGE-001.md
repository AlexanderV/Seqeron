# Validation Report: CODON-USAGE-001 — Codon Usage Analysis

- **Validated:** 2026-06-24   **Area:** Codon
- **Canonical method(s):** `CodonOptimizer.CalculateCodonUsage(string)` (per-codon counts),
  `CodonOptimizer.CompareCodonUsage(string, string)` (Total-Variation-Distance similarity)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN

---

## Scope note (important)

The session prompt frames CODON-USAGE-001 around the full Kazusa/EMBOSS-cusp output set:
per-codon **count**, **frequency per 1000**, **fraction within the synonymous (amino-acid)
family**, and RSCU. The **authoritative unit definition**, however — per
`tests/TestSpecs/CODON-USAGE-001.md`, `docs/Evidence/CODON-USAGE-001-Evidence.md`,
`ALGORITHMS_CHECKLIST_V2.md` (line 1291), and the named test file
`CodonOptimizer_CodonUsage_Tests.cs` — is the narrower pair:

| Method | Class | Role | Output |
|--------|-------|------|--------|
| `CalculateCodonUsage(seq)` | `CodonOptimizer` | Canonical | `Dictionary<string,int>` of in-frame codon **counts** |
| `CompareCodonUsage(seq1, seq2)` | `CodonOptimizer` | Comparison | TVD similarity ∈ [0,1] |

This unit does **NOT** compute frequency-per-1000, per-family fraction, or RSCU. Those richer
Kazusa/cusp outputs live in **other** units:
- **RSCU** → CODON-RSCU-001 / `CodonUsageAnalyzer.CalculateRscu` (TestSpec D1 records this
  separation explicitly: "Tests focus on `CodonOptimizer` methods, not `CodonUsageAnalyzer`").
- **Per-codon frequency / fraction tables** → CODON-STATS-001 / SEQ-CODON-FREQ-001.

So the prompt's "frequency per 1000 + fraction within family + RSCU" emphasis is the general
Kazusa table format, not this unit's deliverable. Verdict for Stage A is PASS-WITH-NOTES purely
because of this scope mismatch; the unit's own description is mathematically correct.

---

## Stage A — Description

### Sources opened & what they confirm

1. **EMBOSS `cusp` manual** (fetched, bioinformatics.nl). Confirms the canonical codon-usage
   table columns: **Number** = "observed number of codons in the input sequences" (raw count);
   **Frequency** = "expected number of codons … per 1000 bases"; **Fraction** = "the proportion
   of usage of the codon among its redundant set" (i.e. within the synonymous family). The manual
   does not specify trailing-incomplete-codon or frame handling.
2. **Kazusa Codon Usage Database** format (search + readme). Each row is
   `[triplet] [aa] [fraction] [frequency:per-thousand] ([number])`, e.g. E. coli
   `UUU F 0.57 19.7 (101)` — fraction = 0.57 of Phe codons, 19.7 per 1000, 101 total. Confirms the
   four output semantics and the per-thousand normalisation = count / total · 1000.
3. **Wikipedia — "Codon usage bias"** (per Evidence doc): degeneracy (64 codons → 20 aa + 3 stop),
   bias = differing synonymous-codon frequencies; lists CAI/Fop/RCA/Nc/RSCU as bias measures.

### Formula check

- **Counts (this unit's `CalculateCodonUsage`).** Count(c) = number of in-frame triplets equal to
  codon c. Matches cusp "Number" / Kazusa "(number)". The richer frequency-per-1000
  (= count/total·1000) and per-family fraction columns are real and correctly defined per the
  sources, but are **not** part of this unit (see Scope note).
- **Comparison (`CompareCodonUsage`).** `Similarity = 1 − (Σ_c |f₁(c) − f₂(c)|) / 2`, with
  f_i(c) = count_i(c)/total_i and the sum over the union of observed codons. This is exactly
  Total Variation Distance similarity (TVD = ½·L¹ between two discrete distributions), giving a
  value in [0,1]. The listed invariants (identity = 1, symmetry, range [0,1], disjoint → 0) are
  genuine mathematical facts of TVD.

### Edge-case semantics

- Empty / null sequence → empty dictionary (convention "no data"). Sourced as standard.
- Incomplete trailing 1–2 nt → ignored (Kazusa / cusp convention; cusp manual silent but this is
  the universal practice).
- Reading frame: **frame 0 only** — the unit reads triplets from index 0 with no frame selection
  or ORF detection. This is the standard assumption (input is a coding sequence already in frame),
  matching cusp/Kazusa which operate on supplied CDS.
- DNA T treated as RNA U (biological equivalence).
- TVD: empty / one-empty → 0 similarity (no comparable data) — a defensible convention.

### Independent cross-check (hand computation)

- **Counts** on ORF `ATGGCTGCTTAA` (DNA → `AUGGCUGCUUAA`, frame 0): codons AUG, GCU, GCU, UAA ⇒
  {AUG:1, GCU:2, UAA:1}, total 4. (If one applied the Kazusa columns this unit doesn't emit:
  per-1000 → AUG 250, GCU 500, UAA 250, sum 1000; GCU fraction within Ala family = 1.0.)
- **TVD M9** `AUGAUGCCCUUU` vs `AUGUUUUUUCCC`: f₁ = {AUG 0.5, CCC 0.25, UUU 0.25};
  f₂ = {AUG 0.25, UUU 0.5, CCC 0.25}. Σ|Δ| = 0.25 + 0 + 0.25 = 0.5 ⇒ sim = 1 − 0.25 = **0.75**. ✓
- **TVD M7** `CUGCUGCUGCUA` vs `CUACUACUACUG`: f₁={CUG 0.75, CUA 0.25}, f₂={CUA 0.75, CUG 0.25};
  Σ|Δ| = 0.5 + 0.5 = 1 ⇒ sim = **0.5**. ✓
- **TVD S6** `AUGGCUAUG` vs `AUGUUUAUG`: f₁={AUG 2/3, GCU 1/3}, f₂={AUG 2/3, UUU 1/3};
  Σ|Δ| = 0 + 1/3 + 1/3 = 2/3 ⇒ sim = 1 − 1/3 = **2/3**. ✓

### Findings / divergences

PASS-WITH-NOTES: scope mismatch only (prompt's frequency-per-1000 / fraction / RSCU framing
belongs to adjacent units). The unit's counts + TVD definitions are correct against EMBOSS cusp,
Kazusa, and TVD theory. No formula is wrong.

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`
- `CalculateCodonUsage` (634–652): `string.IsNullOrEmpty` → empty; `ToUpperInvariant()` then
  `Replace('T','U')`; `SplitIntoCodons` → in-frame triplets; counts into a dictionary.
- `CompareCodonUsage` (657–681): per-codon frequency = count / per-sequence total; returns
  `1 − Σ|f₁−f₂|/2`; returns 0 when union empty or either total = 0.
- `SplitIntoCodons` (687–695): loop `for (i = 0; i + 2 < len; i += 3)` ≡ `i ≤ len−3`, i.e. only
  complete in-frame codons from frame 0; trailing 1–2 nt dropped. Verified for len 3 (last codon at
  0), len 5 (one codon), len 9 (three codons), len 192 (the all-64 test, last codon at 189).

### Formula realised correctly?

Yes. `CalculateCodonUsage` realises Count(c) exactly (raw integer counts, no normalisation —
correctly matching the unit's stated output, not the full Kazusa table). `CompareCodonUsage`
realises TVD similarity exactly: each sequence normalised by its **own** total, ½·L¹, `1 − TVD`.

### Cross-verification table recomputed vs code (all confirmed by passing tests)

| Case | Input | Expected | Spec | Result |
|------|-------|----------|------|--------|
| Counts | `AUGGCUGCU` | {AUG:1, GCU:2} | M1 | ✓ |
| Incomplete | `AUGGC` | {AUG:1} | M3 | ✓ |
| All 64 | each codon ×1 | 64 keys, each 1 | M5 | ✓ |
| DNA→RNA | `ATGGCTTAA` | {AUG:1,GCU:1,UAA:1}, no T | S1 | ✓ |
| Sum invariant | 18 nt | total = 6 | S3 | ✓ |
| TVD 0.5 | CUG/CUA mirror | 0.5 | M7 | ✓ |
| TVD 0.75 | M9 3-codon | 0.75 | M9 | ✓ |
| TVD 0.0/0.25/0.75 | M10 cases | 0.0/0.25/0.75 | M10 | ✓ |
| TVD 0.0 disjoint | UUU vs GGG | 0.0 | S5 | ✓ |
| TVD 2/3 | shared AUG | 0.666… | S6 | ✓ |
| Empty / one-empty | "" / "",seq | 0.0 | M8,S4 | ✓ |

### Variant/delegate consistency

`CompareCodonUsage` delegates to `CalculateCodonUsage` for both inputs — frequencies derived from
the same counting path. Consistent.

### Test quality audit

`CodonOptimizer_CodonUsage_Tests.cs` (22 instances) asserts **exact** values
(`Is.EqualTo(...).Within(1e-10)` for TVD; exact ints for counts), covering all M1–M10 / S1–S6 /
edge (null, too-short, single-codon, mixed-case, T→U, sum invariant, symmetry, disjoint→0,
empty/one-empty). Deterministic; no tautological "no-throw" assertions. Real tests.

### Findings / defects

None. Code faithfully realises the validated formulas. Frame handling (frame 0, trailing-codon
drop) matches the Kazusa/cusp convention.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — counts + TVD definitions are correct against EMBOSS cusp, Kazusa,
  and TVD theory. The only note is the scope mismatch: the prompt's frequency-per-1000 / per-family
  fraction / RSCU framing is the general Kazusa table format and belongs to adjacent units
  (CODON-RSCU-001, CODON-STATS-001, SEQ-CODON-FREQ-001), not to this unit.
- **Stage B:** PASS — `CalculateCodonUsage` (counts, frame 0, T→U, case-insensitive,
  trailing-codon drop) and `CompareCodonUsage` (TVD similarity) both match authoritative sources
  and the hand cross-checks; tests assert the sourced exact values.
- **State:** CLEAN — no defect; no code change required. Build green; CODON-USAGE-001 filter:
  22/22 pass; full `Seqeron.Genomics.Tests` suite: **18208 passed, 0 failed**.

### Sources
- EMBOSS `cusp` manual: https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
- Kazusa Codon Usage Database (format + readme): https://www.kazusa.or.jp/codon/ ,
  https://www.kazusa.or.jp/codon/readme_codon.html
- Wikipedia, "Codon usage bias": https://en.wikipedia.org/wiki/Codon_usage_bias
- Total Variation Distance (½·L¹ between discrete distributions) — standard probability theory.
