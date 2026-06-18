# Validation Report: ASSEMBLY-CORRECT-001 — K-mer spectrum (two-sided) read error correction

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.ErrorCorrectReads(reads, kmerSize, minKmerFrequency)`
  (private helpers `BuildKmerSpectrum`, `CorrectRead`, `IsPositionTrusted`, `AllCoveringKmersTrusted`,
  `CoveringKmerStarts`, `KmerAt`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (after fixing one code-echoing test)
- **End-state:** ✅ CLEAN
- **Test-quality gate:** PASS (one defect found and fixed)

## Stage A — Description

### Sources opened this session (with extracted text)

1. **Musket** — Liu, Schmidt & Maskell (2013), *Bioinformatics* 29(3):308-315
   (https://academic.oup.com/bioinformatics/article/29/3/308/257257, fetched 2026-06-15). Confirmed verbatim:
   - Trusted vs untrusted: *"the k-mers whose multiplicity exceeds a coverage cut-off are deemed to be trusted and otherwise, untrusted."*
   - Cut-off estimation: *"the multiplicity corresponding to the smallest density around the valley is an appropriate estimation of the cut-off."*
   - Trusted-base rule: *"If a base is covered by any trusted k-mer, the base is deemed to be trusted and untrusted, otherwise."*
   - Two-sided goal: *"our two-sided correction aims to find a unique alternative base that makes all k-mers that cover position i trusted."*
   - Speed shortcut: *"Musket chooses to only evaluate both the leftmost and the rightmost k-mers that cover position i on the read, significantly improving speed."*
   - Ambiguity: *"if more than one alternative is found to make both the leftmost and the rightmost k-mers trusted, the base will keep unchanged as a result of ambiguity."*
   - At-most-one-error: *"our two-sided correction conservatively assumes that there is at most one substitution error in any k-mer of a read."*

2. **Quake** — Kelley, Schatz & Salzberg (2010), *Genome Biology* 11:R116
   (https://pmc.ncbi.nlm.nih.gov/articles/PMC3156955/, fetched 2026-06-15). Confirmed: high-coverage = trusted,
   low-coverage = untrusted; error region localized by intersection then union of a read's untrusted k-mers;
   *"search for the maximum likelihood set of corrections that makes all k-mers overlapping the region trusted"*;
   single-base nucleotide-edit model.

3. **Spectral Alignment Problem / EULER (Pevzner et al. 2001)** and **Song & Florea (2018), PMC6311904**
   (web search 2026-06-15): *"a solid k-mer does not likely contain errors, while a weak k-mer most likely contains errors"*;
   SAP = *"take a set of trusted kmers and attempt to find a sequence with minimal corrections such that each k-mer on the corrected sequence is trusted"*;
   spectrum `S = {s ∈ Kk | abundance(s) > m}`.

### Formula / definitions check

The doc's core model (trusted iff `mult(x) ≥ t`; trusted position iff some covering k-mer trusted;
two-sided correction = unique alternative base making covering k-mers trusted; ambiguity ⇒ unchanged;
substitution-only) matches the cited sources exactly. INV-1…INV-5 are genuine properties of this rule.

### Notes / divergences (minor, documented — hence PASS-WITH-NOTES)

- **N1 — "all covering" vs "leftmost+rightmost".** Musket, *for speed*, evaluates only the leftmost and rightmost
  covering k-mers; the repo evaluates **all** covering k-mers. This is a stricter, equally-valid realisation of the
  same "make all covering k-mers trusted" goal, and the doc/Evidence describe it accurately as "all k-mers covering
  position i". For k-mers and reads where the at-most-one-error assumption holds, the two agree on the leftmost and
  rightmost windows; checking the interior windows too is conservative, not wrong. Documented, not a defect.
- **N2 — `≥` vs `>` threshold convention.** SAP/EULER write `abundance > m` and Musket says "exceeds a cut-off";
  the repo uses `mult ≥ minKmerFrequency`. This is purely a parameter-convention choice (off-by-one in how the cut-off
  is named) and is documented in the contract ("k-mers with multiplicity ≥ this value are trusted"). Consistent and
  internally sound.
- **N3 — Auto-cut-off and q-mer/quality weighting, multi-stage iteration, read trimming/discarding** are out of scope
  and honestly declared in §5.3 "Intentionally simplified / Not implemented".

### Independent cross-check (numbers)

Hand + Python reimplementation of the sourced rule. Worked example k=3, cut-off=2, reads `ACGTACGT`×3 + `ACGTTCGT`×1:
spectrum `ACG=7, CGT=8, GTA=3, TAC=3, GTT=1, TTC=1, TCG=1` (recomputed, matches doc); trusted = {ACG,CGT,GTA,TAC};
in the error read only index 4 is untrusted; candidates at index 4: A→{GTA,TAC,ACG} all trusted (unique), C→GTC absent,
G→GTG absent ⇒ unique correction to **A** ⇒ `ACGTACGT`. Reimplementation reproduced M1/M2/M3/S1/S2/C1/C2 outputs exactly.

## Stage B — Implementation

### Code path reviewed

`SequenceAssembler.cs:1104-1268`. `ErrorCorrectReads` validates (`ArgumentNullException` on null reads,
`ArgumentOutOfRangeException` on `kmerSize<1`), builds the spectrum once (`BuildKmerSpectrum`, case-insensitive,
null elements skipped, reads shorter than k contribute nothing), then corrects each read left-to-right. `CorrectRead`
returns the upper-cased read unchanged when `len < k`; otherwise, per untrusted position, it counts candidate bases
(fixed order A,C,G,T) making **all** covering k-mers trusted and applies the substitution only when exactly one
qualifies. Spectrum is fixed during the pass; an applied correction is visible to later positions in the same read.

### Formula realised correctly?

Yes. Trusted-base rule (`IsPositionTrusted`), two-sided all-covering-trusted test (`AllCoveringKmersTrusted`),
unique-alternative rule (`candidateCount == 1`), ambiguity (>1 ⇒ restore), substitution-only (StringBuilder, length
unchanged), determinism (fixed candidate order, single spectrum). Matches the Stage-A description.

### Cross-verification table recomputed vs code (reference reimplementation, all match)

| Case | Input (k, cut-off=2) | Expected (sourced/hand) | Code |
|------|----------------------|-------------------------|------|
| M1 | ACGTACGT×3 + ACGTTCGT, k=3 | all four → ACGTACGT | ✅ |
| M2 | ACGTACGT×3, k=3 | unchanged | ✅ |
| M3 | A,A,C,C,G,T, k=1 | T unchanged (A,C both valid → ambiguous) | ✅ |
| **M4** | AAAAAAAA×3 + AACCAAAA, k=4 | AACCAAAA unchanged (pos 2,3 untrusted, 0 valid candidates) | ✅ |
| C1 | "acg", k=5 | "ACG" (no covering k-mers) | ✅ |
| C2 | [], k=3 | [] | ✅ |

### Variant/delegate consistency

Single public method; no `*Fast`/instance variant. Defaults (`kmerSize=15`, `minKmerFrequency=2`) are non-behavioral
for the tested contract (every behavioral test passes both explicitly).

### Test quality audit (the gate)

- **Defect found & fixed (A22).** Original M4 used `[ACGTACGT, ACGTACGT, TTTTTTTT]` and claimed to test the
  "no correcting set" branch. But `TTT` has multiplicity 6 ≥ 2, so **every** position of `TTTTTTTT` is covered by a
  trusted k-mer and is skipped by the **trusted-base rule** (the same path as M2/INV-3). The candidate-search /
  no-valid-correction branch (INV-4) was never reached, so M4 would still pass against an implementation whose
  no-correction branch was broken — a code-echoing blind spot. Confirmed by spectrum recomputation
  (`TTT=6` trusted) and by tracing every position of `TTTTTTTT` as trusted.
- **Fix:** rewrote M4 to a genuine no-valid-correction case — genome `AAAAAAAA`×3 (only `AAAA` trusted, mult=16) plus
  a two-adjacent-error read `AACCAAAA` (k=4). Positions 2 and 3 are each covered solely by untrusted 4-mers and the
  candidate search returns **zero** valid bases (the remaining error keeps a covering window weak), so the read stays
  `AACCAAAA`. Verified by hand and by the reference reimplementation (`pos2={}`, `pos3={}`). The test now exercises the
  INV-4 no-correction path and fails against a broken branch.
- **No green-washing:** no assertion weakened, no tolerance widened, no skip; no expected value tuned to code output —
  M4's expected value is the externally-derived no-correction outcome.
- **Coverage:** M1 (unique correction), M2 (trusted-base rule), M3 (ambiguity), M4 (no valid correction), M5
  (count/length INV-1/2), M6/M7 (null / k<1 exceptions), S1 (error-free passthrough), S2 (case-insensitive), S3
  (determinism), C1 (read < k), C2 (empty list), P1 (fixed-seed property INV-2/INV-3). All public-method branches and
  Stage-A edge/error cases are exercised.

### Findings / defects

- **A22** (test-quality, fixed this session) — see above. No algorithm defect.

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** (N1 all-covering vs leftmost+rightmost; N2 `≥` vs `>` convention; N3 documented scope
  reductions — all accurate and disclosed).
- **Stage B: PASS** — code faithfully realises the validated rule; the one code-echoing test (M4) was rewritten to a
  genuine, source-traced no-correction case.
- **End-state: ✅ CLEAN.** Build 0 errors (changed test file warning-free); full unfiltered suite **6535 passed,
  0 failed**.
- Pre-existing, unrelated: 4 NUnit2007 warnings in `ApproximateMatcher_EditDistance_Tests.cs` (surface only on
  `--no-incremental`); not in scope for this unit and not introduced here.
