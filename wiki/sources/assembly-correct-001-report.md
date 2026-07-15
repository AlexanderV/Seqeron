---
type: source
title: "Validation report: ASSEMBLY-CORRECT-001 (k-mer spectrum two-sided read error correction)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-CORRECT-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-CORRECT-001.md
source_commit: 02551f587247480dd3dff1cbb59fba61ed5ffae2
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-CORRECT-001

The two-stage **validation write-up** for test unit **ASSEMBLY-CORRECT-001** (k-mer spectrum,
two-sided read error correction — the error-correction preprocessing step ahead of assembly),
validated 2026-06-15 in a fresh context. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm description and
the shipped code. The algorithm itself is summarized in the concept
[[kmer-spectrum-error-correction]]; the two-stage methodology is the [[validation-protocol]], and the
overall strategy is [[validation-and-testing]]. Distinct from [[assembly-correct-001-evidence]] (the
pre-implementation evidence artifact) — this is the independent re-validation verdict.

Canonical method under test:
`SequenceAssembler.ErrorCorrectReads(reads, kmerSize, minKmerFrequency)` (private helpers
`BuildKmerSpectrum`, `CorrectRead`, `IsPositionTrusted`, `AllCoveringKmersTrusted`,
`CoveringKmerStarts`, `KmerAt`).

## Verdict

**Stage A: 🟡 PASS-WITH-NOTES · Stage B: ✅ PASS · State: ✅ CLEAN · Test-quality gate: PASS.**
One code-echoing test (M4) was found and rewritten (no algorithm defect). Full unfiltered suite
**6535 passed, 0 failed**; `dotnet build` 0 errors (the 4 pre-existing NUnit2007 warnings in
`ApproximateMatcher_EditDistance_Tests.cs` surface only on `--no-incremental` and are out of scope).

## Stage A — description (algorithm faithfulness)

Theory checked against primary sources opened this session, quoted verbatim, not the repo's own
assertions:

- **Musket** (Liu, Schmidt & Maskell 2013, *Bioinformatics* 29(3):308) — the two-sided rule:
  trusted iff multiplicity **exceeds a coverage cut-off** (auto-chosen from "the smallest density
  around the valley"), a base is trusted "if a base is covered by any trusted k-mer," two-sided
  correction "aims to find a unique alternative base that makes all k-mers that cover position *i*
  trusted," evaluates **leftmost + rightmost** covering k-mers "significantly improving speed,"
  keeps the base **unchanged on ambiguity**, and "conservatively assumes … at most one substitution
  error in any k-mer."
- **Quake** (Kelley, Schatz & Salzberg 2010, *Genome Biology* 11:R116) — corroborates high- vs
  low-coverage = trusted/untrusted, localizes the error region by intersection-then-union of a
  read's untrusted k-mers, and searches "the maximum likelihood set of corrections that makes all
  k-mers overlapping the region trusted" (single-base nucleotide-edit model).
- **Spectral Alignment Problem / EULER** (Pevzner et al. 2001) and **Song & Florea 2018**
  (PMC6311904) — same idea under *solid* vs *weak* k-mers; SAP = "find a sequence with minimal
  corrections such that each k-mer on the corrected sequence is trusted"; spectrum
  `S = {s ∈ Kk | abundance(s) > m}`.

**Formula check.** The doc's core model (trusted iff `mult(x) ≥ t`; trusted position iff some
covering k-mer trusted; two-sided correction = unique alternative base making covering k-mers
trusted; ambiguity ⇒ unchanged; substitution-only) matches the cited sources exactly; INV-1..INV-5
are genuine properties.

**Independent cross-check (numbers).** Hand + Python reimplementation of the sourced rule. Worked
example `k=3`, cut-off `2`, reads `ACGTACGT`×3 + `ACGTTCGT`×1 → spectrum `ACG=7, CGT=8, GTA=3,
TAC=3` trusted vs `GTT=1, TTC=1, TCG=1` untrusted; in the error read only index 4 is untrusted;
candidates there: A→{GTA,TAC,ACG} all trusted (unique), C→GTC absent, G→GTG absent ⇒ unique
correction to **A** ⇒ `ACGTACGT`. The reimplementation reproduced M1/M2/M3/S1/S2/C1/C2 exactly.

**Notes (why Stage A is PASS-WITH-NOTES — documented, not defects).**

- **N1 — "all-covering" vs "leftmost+rightmost".** Musket, *for speed*, evaluates only the leftmost
  and rightmost covering k-mers; the repo evaluates **all** covering k-mers. A stricter,
  equally-valid realisation of the same "make all covering k-mers trusted" goal; under the
  ≤1-error-per-k-mer assumption the two agree, and checking interior windows is conservative.
- **N2 — `≥` vs `>` threshold.** SAP/EULER write `abundance > m`, Musket says "exceeds a cut-off";
  the repo uses `mult ≥ minKmerFrequency`. A parameter-naming convention choice, documented in the
  contract ("k-mers with multiplicity ≥ this value are trusted").
- **N3 — out-of-scope, honestly declared** (§5.3): auto-cut-off, q-mer/quality weighting,
  multi-stage iteration, read trimming/discarding.

## Stage B — implementation (code review)

Code path `SequenceAssembler.cs:1104-1268`. `ErrorCorrectReads` validates (`ArgumentNullException`
on null reads, `ArgumentOutOfRangeException` on `kmerSize<1`), builds the spectrum once
(`BuildKmerSpectrum`, case-insensitive, null elements skipped, reads shorter than k contribute
nothing), then corrects each read left-to-right. `CorrectRead` returns the upper-cased read
unchanged when `len < k`; otherwise, per untrusted position, it counts candidate bases (fixed order
A,C,G,T) making **all** covering k-mers trusted and applies the substitution only when exactly one
qualifies. Spectrum is fixed during the pass; an applied correction is visible to later positions in
the same read.

**Formula realised correctly?** Yes — trusted-base rule (`IsPositionTrusted`), two-sided
all-covering-trusted test (`AllCoveringKmersTrusted`), unique-alternative rule (`candidateCount ==
1`), ambiguity (>1 ⇒ restore), substitution-only (StringBuilder, length unchanged), determinism
(fixed candidate order, single spectrum). Cross-verification table (reference reimplementation, all
match):

| Case | Input (k, cut-off=2) | Expected (sourced/hand) | Code |
|------|----------------------|-------------------------|------|
| M1 | ACGTACGT×3 + ACGTTCGT, k=3 | all four → ACGTACGT | ✅ |
| M2 | ACGTACGT×3, k=3 | unchanged | ✅ |
| M3 | A,A,C,C,G,T, k=1 | T unchanged (A,C both valid → ambiguous) | ✅ |
| M4 | AAAAAAAA×3 + AACCAAAA, k=4 | AACCAAAA unchanged (pos 2,3 untrusted, 0 valid candidates) | ✅ |
| C1 | "acg", k=5 | "ACG" (no covering k-mers) | ✅ |
| C2 | [], k=3 | [] | ✅ |

**Variant/delegate consistency.** Single public method; no `*Fast`/instance variant. Defaults
(`kmerSize=15`, `minKmerFrequency=2`) are non-behavioral for the tested contract (every behavioral
test passes both explicitly).

## Findings

- **A22 (test-quality, fixed this session — no algorithm defect).** The original M4 used
  `[ACGTACGT, ACGTACGT, TTTTTTTT]` claiming to test the "no correcting set" branch, but `TTT` has
  multiplicity 6 ≥ 2, so **every** position of `TTTTTTTT` is covered by a trusted k-mer and is
  skipped by the trusted-base rule (same path as M2/INV-3) — the candidate-search /
  no-valid-correction branch (INV-4) was never reached, a code-echoing blind spot that would still
  pass against a broken no-correction branch. **Fix:** rewrote M4 to a genuine no-valid-correction
  case — genome `AAAAAAAA`×3 (only `AAAA` trusted, mult=16) plus a two-adjacent-error read
  `AACCAAAA` (k=4); positions 2 and 3 are each covered solely by untrusted 4-mers and the candidate
  search returns **zero** valid bases, so the read stays `AACCAAAA`. Verified by hand and by the
  reference reimplementation (`pos2={}`, `pos3={}`); the test now exercises the INV-4 no-correction
  path and fails against a broken branch. No assertion weakened, no tolerance widened, no skip, no
  expected value tuned to code output.
- **No algorithm defect.** Code faithfully realises the validated rule.
- **Test coverage** spans M1 (unique correction), M2 (trusted-base rule), M3 (ambiguity), M4 (no
  valid correction), M5 (count/length INV-1/2), M6/M7 (null / k<1 exceptions), S1 (error-free
  passthrough), S2 (case-insensitive), S3 (determinism), C1 (read < k), C2 (empty list), P1
  (fixed-seed property) — all public-method branches and Stage-A edge/error cases.
- **No follow-ups.**
