# Validation Report: COMPGEN-ANI-001 — Average Nucleotide Identity (ANI)

- **Validated:** 2026-06-16   **Area:** Comparative Genomics
- **Canonical method(s):** `ComparativeGenomics.CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction)` (+ private `BestUngappedFragmentMatch`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Goris et al. 2007** (IJSEM 57:81-91, DOI:10.1099/ijs.0.64483-0) — method text retrieved via
   WebSearch (reproduced verbatim in the secondary literature / pyani docs). Confirms, verbatim:
   - Query genome "cut into consecutive 1020 nt fragments" to mirror ~1 kb DDH fragmentation.
   - Fragments searched against the reference whole genome with **BLASTN**; best match saved.
   - BLAST settings: `X=150`, `q=-1` (mismatch penalty), `F=F`.
   - **ANI = "mean identity of all BLASTN matches that showed more than 30 % overall sequence
     identity (recalculated to an identity along the entire sequence) over an alignable region of
     at least 70 % of their length."**
   - "ANI values of approximately 95 % correspond to the 70 % DNA-DNA hybridization standard."
2. **Konstantinidis & Tiedje 2005** (PNAS 102:2567-2572) — WebSearch confirms "ANI values of
   ≈94 % corresponded to the traditional 70 % DNA–DNA reassociation standard," and that newer
   analyses favour the ≈95–96 % species boundary.
3. **pyani ANIb source** (`pyani/anib.html`, WebFetch of the module) — confirms the exact identity
   and coverage convention used by the canonical reference implementation:
   - `ani_pid = ani_alnids / qlen` → **identity recalculated over the whole query fragment length**,
     not over the aligned region.
   - `ani_coverage = ani_alnlen / qlen` → coverage = aligned length / fragment length.
   - filter: `data[(ani_coverage > 0.7) & (ani_pid > 0.3)]` → **>70 % coverage, >30 % identity**.

### Formula check
The repo's documented model (`Average_Nucleotide_Identity.md` §2.2) is:
`id(f_i) = matching bases / L`, `cov(f_i) = alignable length / L`, qualify iff `id > 0.30 ∧ cov ≥ 0.70`,
`ANI = mean of qualifying id(f_i)`. This matches Goris 2007 and the pyani convention exactly
(identity recalculated over the full fragment length; strict `>30 %`; `≥/>70 %` coverage). Defaults
`L=1020`, `0.30`, `0.70` all trace to Goris 2007. ✔

### Edge-case semantics
- Trailing partial fragment dropped (consecutive non-overlapping fragmentation). ✔ (Goris "consecutive 1020 nt fragments")
- Non-conserved fragments (below identity or coverage cut-off) discarded, not counted as 0. ✔ (Goris; pyani filter)
- Empty / single-fragment / no-qualifying-fragment → 0 (mean over empty set reported as 0). ✔ (definition; convention)
- Direction-dependence (query is the fragmented genome) acknowledged. ✔

### Independent cross-check (numbers, this session)
Re-derived every deterministic test value from the Goris mean-identity formula in **Python**
(`/tmp/ani_check.py`, independent of the C# code), reproducing the matching-bases / fragment-length
identity and the >30 %/≥70 % filter:

| Case | Inputs (fragLen 4) | Re-derived ANI | TestSpec |
|------|--------------------|----------------|----------|
| M1 identical | R,R | 1.0 | 1.0 ✔ |
| M2 1 sub (`…TTTA`) | last frag 3/4 | 0.9375 | 0.9375 ✔ |
| M3 half (`…AATT`) | last frag 2/4 | 0.875 | 0.875 ✔ |
| M4 id cut-off | `AAAACGTC` vs `AAAAAAAA`, frag2=0/4 excluded | 1.0 | 1.0 ✔ |
| M5 align cut-off | `AAAA` vs `AA` (ref<frag) | 0.0 | 0.0 ✔ |
| M6 trailing | `AAAACCCCXX` (2 full frags) | 1.0 | 1.0 ✔ |
| S2 short query | `AAA` < fragLen | 0.0 | 0.0 ✔ |
| C1 default | `AAAAACGT` vs `AAAAAAAA`, frag2=0.25 excluded | 1.0 | 1.0 ✔ |
| C1 minId 0.20 | frag2=0.25 now kept | 0.625 | 0.625 ✔ |

All match. Also confirmed `AATT` best placement against R is 2/4 at several offsets (no better window exists).

### Findings / divergences (Stage A)
The unit implements **single-direction, ungapped** ANI rather than the full reciprocal gapped-BLASTN
ANIb. This is **explicitly and honestly documented** as a simplification (algorithm doc §5.3/§5.4,
Evidence ASSUMPTION 1): for substitution-only divergence the ungapped "matching bases / fragment
length" identity equals the Goris recalculated-over-fragment identity; indels and reciprocal
averaging are out of scope. The description does not overclaim taxonomy-grade ANI. The abstract
maths (the mean-identity formula and cut-offs) is correct and sourced. → **Stage A PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:1069-1135`.

- Validation: null/empty → 0 (`:1076`); `fragmentLength ≤ 0` → `ArgumentOutOfRangeException` (`:1078`). ✔
- Fragmentation: `for (start=0; start+fragmentLength <= length; start += fragmentLength)` — consecutive,
  non-overlapping, trailing partial dropped (`:1086`). ✔ Matches INV-04.
- Per fragment: `BestUngappedFragmentMatch` slides the fragment over every full-length offset of the
  reference, counts matching bases at the best offset, identity = bestMatches / fragLen (`:1108-1135`). ✔
- Coverage: returns `1.0` whenever ref ≥ fragLen, else `(0,0)` (ref shorter → excluded). This is the
  correct consequence of the ungapped full-length design. ✔
- Filter: `identity > minIdentity && alignableFraction >= minAlignableFraction` (`:1092`) — strict `>`
  on identity (matches Goris "more than 30 %"), `>=` on coverage; mean of kept identities, else 0 (`:1096`). ✔

### Formula realised correctly?
Yes. The code computes the validated mean-identity formula with identity recalculated over the full
fragment length and the >30 %/≥70 % cut-offs. The earlier LCS-length proxy defect (noted in §5.2 of
the doc) has been replaced by a true matching-base count; M2/M3/M4 lock this in.

### Cross-verification vs code
The full unfiltered test suite (which runs the C# `CalculateANI` on exactly the inputs above)
passes; every value equals the independently Python-re-derived value in the Stage-A table. ✔

### Variant/delegate consistency
Single public method, single signature with defaults. No `*Fast`/delegate variants. The MCP wrapper
(`AnalysisTools.cs`) forwards to the same method. No divergence.

### Test quality audit (HARD gate)
File: `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs` (11 tests).

- **Sourced, not code-echoed:** every deterministic expected value is hand-derived from the Goris
  formula and independently reproduced in Python this session. M4 explicitly notes a wrong
  "average-all-fragments" implementation would return 0.5 and thus fail — a real discriminating
  oracle, not a tautology. M2/M3 separate true identity from an LCS proxy. ✔
- **No green-washing:** exact equality (`Within(1e-10)` on exact fractions) for every known value;
  no weakened `Greater/AtLeast/Contains`; no widened tolerances; no skipped/ignored tests. The only
  `Is.InRange` (S1) is the range *invariant* INV-01 over arbitrary pairs — appropriate, not a
  weakened exact assertion. ✔
- **Coverage of logic/branches:** identity recalc (M2/M3), identity cut-off strict-`>` (M4/C1),
  coverage cut-off (M5), non-overlapping fragmentation + trailing drop (M6), null/empty (M7),
  invalid fragmentLength (M8), range invariant incl. default-1020-length run (S1), query<fragLen
  (S2), custom `minIdentity` (C1). All Stage-A branches/edge cases covered. ✔
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6605`; build `0 Error(s)`. ✔

**Gate result: PASS.**

### Findings / notes (Stage B)
1. **`minAlignableFraction` parameter is structurally inert (not a defect).** Because the matcher is
   ungapped full-length, coverage is always exactly `1.0` (ref ≥ fragLen) or `0.0` (ref < fragLen).
   Any `minAlignableFraction ∈ (0,1]` therefore yields identical behaviour, so no test varies it.
   This is an inherent consequence of the documented ungapped simplification, not a coding error; it
   cannot be meaningfully exercised without introducing gapped/partial alignment. Noted as a
   limitation, consistent with the algorithm doc.
2. **Simplification vs full ANIb (documented).** Single-direction, ungapped, no IUPAC/reverse-
   complement, no reciprocal averaging — accurately disclosed in the doc (§5.3/§6.2) and Evidence.
   Not taxonomy-grade; the well-defined averaging formula is what is under test and it is correct.

## Verdict & follow-ups
- **Stage A: PASS** — description matches Goris 2007 / pyani verbatim; thresholds and species
  boundary independently confirmed; simplification honestly scoped.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated (simplified) formula; tests
  are exact, sourced, discriminating, and cover all branches; full suite green. Notes: the
  `minAlignableFraction` parameter is inert under the ungapped design, and the method is a documented
  non-taxonomy-grade simplification of full ANIb.
- **End-state: CLEAN** — no defect found; no code or test changes required. Working tree only gains
  this report + ledger/register rows.
