# Validation Report: MOTIF-REGULATORY-001 — Regulatory Element Scan

- **Validated:** 2026-06-16   **Area:** Matching (Motif Discovery)
- **Canonical method(s):** `MotifFinder.FindRegulatoryElements(DnaSequence)`; `MotifFinder.KnownMotifs.*` constants; underlying `MotifFinder.FindDegenerateMotif(DnaSequence, string)` IUPAC scan
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Summary

`FindRegulatoryElements` scans a DNA sequence against a fixed library of 12 published regulatory
consensus motifs (3 eukaryotic core-promoter elements, 2 prokaryotic σ70 hexamers, Kozak context,
Shine-Dalgarno, poly(A) signal, and 4 transcription-factor binding sites). Each entry is matched as
an exact (or, for the E-box, IUPAC-degenerate) consensus string and every 0-based occurrence is
reported. All 12 consensus strings were independently confirmed against external authoritative
sources retrieved this session (see table below). No defect was found; two SHOULD-level tests were
strengthened from `Contains`/partial assertions to exact sourced values. Full suite green.

## Stage A — Description

### Sources opened this session (independent retrieval, not the repo's own citations)

| Element | Repo pattern | Source-confirmed value | Source (retrieved 2026-06-16) |
|---|---|---|---|
| TATA Box | `TATAAA` | consensus `5'-TATA(A/T)A(A/T)-3'`, canonical core `TATAAA`, ~25–35 bp upstream of TSS | en.wikipedia.org/wiki/TATA_box |
| -10 Box (Pribnow) | `TATAAT` | "a sequence of *TATAAT* of six nucleotides" | en.wikipedia.org/wiki/Pribnow_box |
| -35 Box | `TTGACA` | "-35 sequence (5′-TTGACA-3′)" σ70 promoter | WebSearch: ASM/NAR σ70 promoter element refs |
| CAAT Box | `CCAAT` | "CCAAT pentanucleotide present in ∼30% of promoters" (Bucher 1990); extended `RRCCAATSR` | WebSearch: NAR 26(5):1135 survey of 178 CCAAT boxes / Bucher 1990 |
| GC Box | `GGGCGG` | "consensus sequence GGGCGG … binding site for zinc finger proteins" (Sp1) | en.wikipedia.org/wiki/GC_box |
| Kozak | `GCCGCCACCATGG` | `gccRccAUGG` (R=purine, A more frequent); -3 purine, +4 G ⇒ most-preferred DNA string `GCCGCCACCATGG` | en.wikipedia.org/wiki/Kozak_consensus_sequence |
| Shine-Dalgarno | `AGGAGG` | "six-base consensus sequence is AGGAGG"; complementary to 3' 16S rRNA | en.wikipedia.org/wiki/Shine–Dalgarno_sequence |
| Poly(A) Signal | `AATAAA` | hexanucleotide `A(A/U)UAAA`, canonical `AAUAAA` (DNA `AATAAA`); ~10–30 nt upstream of cleavage | WebSearch: Proudfoot & Brownlee 1976 / Genes Dev review |
| E-box | `CANNTG` | "DNA sequence CANNTG (where N can be any nucleotide), palindromic canonical CACGTG" | en.wikipedia.org/wiki/E-box |
| AP-1 (TRE) | `TGACTCA` | "consensus AP-1 binding site is the palindrome TGA(C/G)TCA"; canonical TRE 7-bp `TGACTCA` | WebSearch: Lee/Mitchell/Tjian 1987; PMC2693225 |
| NF-κB | `GGGACTTTCC` | Sen & Baltimore 1986 bound `GGGACTTTCC` (Ig κ enhancer); consensus `GGGRNWYYCC` | WebSearch: Cell "30 years of NF-κB"; PMC NF-κB κ enhancer |
| CREB (CRE) | `TGACGTCA` | "8-base palindrome 5'-TGACGTCA-3'" in somatostatin CRE (Montminy et al. 1986) | WebSearch: Montminy 1986; wikidoc CRE; PubMed 7912577 |

### Formula / definition check

The "core model" (§2.2 of the algorithm doc) is a per-position IUPAC window scan: report every start
`0 ≤ i ≤ n−m` where `S[i+j]` ∈ IUPAC-set(`P[j]`) for all `j`. This is the standard exact/degenerate
consensus-matching definition; coordinates are 0-based inclusive start; only the given strand is
scanned (documented limitation). All conventions are explicit and standard.

### Edge-case semantics

- Null sequence → `ArgumentNullException` (contract) — sourced as a contract, standard.
- Empty sequence → empty result (no offset satisfies `0 ≤ i ≤ n−m`) — correct window definition.
- Multiple/overlapping occurrences → all reported (exhaustive scan) — correct.
- E-box `CANNTG` IUPAC-degenerate (N = any) → `CACGTG`, `CAGCTG`, … all match — correct per source.

### Independent cross-check (hand-computed positions)

- `GGGTATAAAGGG` → `TATAAA` starts at index 3 (G0 G1 G2 T3…). ✓
- `TTGCCGCCACCATGGAA` → `GCCGCCACCATGG` (len 13) at index 2, ends index 14. ✓
- `AATAAACGAATAAA` → `AATAAA` at indices 0 and 8 (A8 A9 T10 A11 A12 A13). ✓
- `TATAAAAGGAGG` → `TATAAA` at 0; `AGGAGG` at 6 (A6 G7 G8 A9 G10 G11). ✓

### Findings / divergences (Stage A)

Two documented representative-site simplifications (both source-backed, recorded as assumptions):
1. **NF-κB** scanned as the strong reference site `GGGACTTTCC` rather than expanding the degenerate
   consensus `GGGRNWYYCC`. The exact string is the Sen & Baltimore (1986) Ig κ enhancer site — not
   fabricated. Consequence: weaker/variant κB sites are not reported. Documented as a limitation.
2. **Kozak** scanned as the single most-preferred-base string `GCCGCCACCATGG` rather than expanding
   the -3 purine / +4 G degeneracy. The -3 = A is the more frequent purine per Kozak (1987), so the
   exact string is the canonical optimal context. Consequence: suboptimal Kozak contexts not reported.
   Documented as a limitation.

These are honest, explicitly-scoped simplifications, not biological errors. Note on AP-1: the general
consensus is degenerate `TGA(C/G)TCA`; the repo scans the single canonical TRE `TGACTCA`. The prior
repo value `TGAGTCA` (the other allowed member) was changed to the canonical `TGACTCA` — defensible
since `TGACTCA` is the archetypal/most-studied TRE. (A fuller implementation could use `TGASTCA`; the
exact-string choice is a documented simplification, consistent with NF-κB/Kozak.)

**Stage A verdict: PASS.** Every consensus string matches an independently-retrieved authoritative
source; definitions, conventions, edge cases, and invariants are correct.

## Stage B — Implementation

### Code path reviewed

- `MotifFinder.cs:685–717` — `FindRegulatoryElements`: null-checks, builds the 12-entry library from
  `KnownMotifs.*`, scans each via `FindDegenerateMotif`, yields `RegulatoryElement(Name, Position,
  Sequence, Pattern, Description)`.
- `MotifFinder.cs:639–678` — `KnownMotifs` constants.
- `MotifFinder.cs:90–119` — `FindDegenerateMotif(DnaSequence, string)`: per-position IUPAC window scan
  using `IupacCodes[motifChar].Contains(seqChar)`.
- `MotifFinder.cs:46–63` — `IupacCodes` table (N = "ACGT", R = "AG", etc.). All entries correct.

### Formula realised correctly?

Yes. The scan loops `for i in [0, len−m]` and tests each position against the IUPAC set — exactly the
Stage-A definition. `MatchedSequence = seq.Substring(i, m)` and `Pattern` is the (upper-cased) library
string, so INV-01 (length equality) and INV-02 (IUPAC match) hold by construction; the linear window
gives INV-03/INV-04 (exhaustiveness). Constants match the validated values verbatim (INV-04 / no
fabrication). Null → `ArgumentNullException` (`ArgumentNullException.ThrowIfNull`, line 687/92).

### Cross-verification table (sourced expectation vs actual code)

All 12 probes + edge cases run via the canonical test class and pass: TATA@3, -10@2, -35@2, CAAT@2,
GC@2, Kozak@2, SD@2, polyA@2, E-box `CACGTG`@2, AP-1 `TGACTCA`@2, NF-κB@2, CREB@2; AP-1 `TGAGTCA`→0
hits; null→throws; empty→empty; repeated polyA→{0,8}; multi-element→TATA@0+SD@6; no-consensus→empty.
Every expected value traces to the sourced consensus, not to code output.

### Variant/delegate consistency

`FindRegulatoryElements` uses the `DnaSequence` overload of `FindDegenerateMotif` (lines 90–119); the
`CancellationToken` overload (`FindDegenerateMotifCore`, 148–201) uses an equivalent `switch`-based
IUPAC test producing identical match semantics. Both are consistent for the codes used here.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** every M1–M16 assertion uses the externally-confirmed consensus string,
  exact 0-based position, and (where applicable) exact matched substring. M11 is a genuine negative
  regression (`TGAGTCA` ⇒ 0 AP-1 hits) that would fail a wrong implementation. M16 locks all 12
  constants to their sourced values.
- **No green-washing found in the as-shipped MUST tests** (exact `Is.EqualTo`, no widened tolerances,
  none skipped/ignored).
- **Strengthened two SHOULD tests this session** (test-quality improvements, not code fixes):
  - **S2** (`MultipleDistinctElements`): was `Does.Contain(name)` only → now asserts exactly one
    TATA box at position 0 and exactly one Shine-Dalgarno at position 6 (hand-sourced positions).
  - **S3** (`NoConsensusPresent`): asserted only that GC Box and TATA Box subsets were empty → now
    asserts the **entire** scan result `Is.Empty`, so any spurious match across all 12 elements is
    caught.
- **Coverage:** all 16 MUST, 3 SHOULD, 1 COULD cases present; the single public method, the
  constants, null/empty/boundary, the degenerate E-box path, and the AP-1 negative case are all
  exercised. C1's `Is.GreaterThan(0)` is a sanity precondition guarding a stronger `All(...)`
  invariant assertion — acceptable.
- **Honest green:** FULL unfiltered suite = **6606 passed, 0 failed, 1 skipped** (unrelated
  `MFE_Benchmark_AllScenarios`); `dotnet build` 0 errors, no new warnings.

**Test-quality gate: PASS** (after strengthening S2 and S3).

### Findings / defects (Stage B)

No code defect. The implementation faithfully realises the validated description. Two pre-existing
SHOULD tests were under-asserting; both strengthened to exact sourced values this session.

## Verdict & follow-ups

- **Stage A: PASS** — all 12 consensus strings independently confirmed against external sources.
- **Stage B: PASS** — code realises the formula; constants and edge cases correct; tests now lock
  exact sourced values throughout.
- **End-state: CLEAN.** No defect; the two test-quality gaps found were completely fixed and the full
  suite is green.
- **Non-blocking notes (documented limitations, not defects):** single-strand scan; NF-κB and Kozak
  (and effectively AP-1) use representative exact strings rather than the full degenerate consensus —
  all explicitly documented in the algorithm doc §5.3/§6.2 and the spec Assumption Register.
