# Validation Report: RESTR-FILTER-001 — Restriction Enzyme Filtering

- **Validated:** 2026-06-15   **Area:** MolTools
- **Canonical method(s):** `RestrictionAnalyzer.GetEnzymesByCutLength(int,int)`, `GetEnzymesByCutLength(int)`, `GetBluntCutters()`, `GetStickyCutters()`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one Stage-B test-quality defect fixed in-session; no code defect)

## Stage A — Description

The unit selects subsets of a curated Type II enzyme library by (a) recognition-sequence length
(exact or inclusive range) and (b) end type (blunt vs sticky), the latter derived from cut-position
equality (`cf == cr` ⇒ center cut ⇒ blunt; otherwise staggered ⇒ overhang).

### Sources opened this session (independent of the repo artifacts)
- **Wikipedia "Sticky and blunt ends"** (definitions) — blunt = both strands terminate in a base pair; sticky = a 5'/3' overhang. Every end is one or the other ⇒ total partition.
- **Wikipedia "Restriction enzyme"** — Type II undivided sites "4–8 nucleotides"; center cut → blunt, staggered → sticky.
- **NEB R0103 HincII** (WebSearch) — `GTY^RAC`, blunt.
- **NEB R0167 RsaI** (WebSearch) — `GT^AC`, blunt.
- **NEB R0176 DpnI** + Wikipedia DpnI + addgene/UTech (WebSearch) — `Gm6A^TC`, **blunt** (cut at the center of the palindrome; offset 2/2). The first WebFetch model rendering ("5' overhang") was a misread and was corrected against three further sources.
- **ScaI / StuI / SwaI** (WebSearch: ScaI `AGT^ACT` blunt; StuI `AGG^CCT` blunt; SwaI `ATTT^AAAT` blunt — "cleaves at its center to generate blunt-ended fragments", Stachowiak et al.).
- **NdeI** (WebSearch) — `CA^TATG`, 5' overhang (sticky). **PacI** — `TTAAT^TAA`, 3' overhang (sticky).
- **SfiI** (Springer/PMC548270 + NEB R0123, WebSearch) — interrupted palindrome `5'-GGCCNNNN^NGGCC-3'`, 3 nt 3' overhang, "eight nucleotide specificity"; the recognition *string* is 13 nt.

### Formula / definition check
- Blunt criterion `cf == cr` matches the cited "center cut ⇒ blunt" rule (Wikipedia, Restriction enzyme).
- Recognition length = string length; the 4–8 nt range applies to *undivided* sites. SfiI (13-nt interrupted palindrome) legitimately lies outside [4,8]. The 8-nt-specificity-vs-13-nt-string distinction is a documented modeling choice and is internally consistent.
- Range inclusivity is an API-shape assumption (documented); recognition lengths are all source-backed.

### Independent cross-check (numbers)
Library enumerated and recomputed locally (39 enzymes):
- **Length histogram:** {4: 9, 6: 24, 8: 5, 13: 1}. Only SfiI (13) is outside [4,8] ⇒ M5 `count == total − 1` (38) is structurally correct.
- **Blunt set (cf==cr), 10 enzymes:** AluI, DpnI, EcoRV, HaeIII, HincII, RsaI, ScaI, SmaI, StuI, SwaI — **every one independently confirmed blunt** against REBASE/NEB/Wikipedia this session.
- **Sticky set, 29 enzymes** — spot-checked NdeI (5'), PacI (3'), KpnI/PstI (3'), EcoRI (5'), all confirmed.

Stage A verdict: **PASS** — biology and conventions are correct and sourced.

## Stage B — Implementation

- Code path: `RestrictionAnalyzer.cs:84-115` (four filters), `RestrictionEnzyme.IsBluntEnd` (`:480`), `RecognitionLength` (`:494`).
- The filters are exactly `Where(len==L)`, `Where(min<=len<=max)`, `Where(IsBluntEnd)`, `Where(!IsBluntEnd)` — faithful to the validated description. Total methods, never throw/return-null; inverted/non-positive ranges yield empty (INV-05, C1) verified.
- Library cut-position entries each cross-checked: all 10 blunt entries and the spot-checked sticky entries match REBASE/NEB. **No code defect found.**

### Cross-verification (recomputed vs code)
| Check | External value | Code |
|-------|---------------|------|
| `[6,6]` count | 24 6-cutters | 24 ✓ |
| `[4,4]` count | 9 4-cutters | 9 ✓ |
| `[4,8]` count | 39 − 1 (SfiI) = 38 | 38 ✓ |
| `[9,10]` | ∅ | ∅ ✓ |
| Blunt set | exactly {AluI,DpnI,EcoRV,HaeIII,HincII,RsaI,ScaI,SmaI,StuI,SwaI} | matches ✓ |
| Blunt ∩ Sticky / ∪ | ∅ / library | ✓ |

### Test quality audit — one defect fixed
- **Defect (Stage-B test weakness):** M1's `GetBluntCutters_KnownEnzymes_*` checked only 4 of 10 blunt members present + 2 sticky absent. A wrong cut-position for any of the 6 unchecked blunt enzymes (RsaI, ScaI, StuI, SwaI, HincII, DpnI) would misclassify it as sticky and **the test would still pass** — i.e. it would survive a deliberately-wrong implementation.
- **Fix:** added `GetBluntCutters_ReturnsExactlyTheDocumentedBluntSet` asserting `Is.EquivalentTo` the full externally-sourced 10-enzyme blunt set. With M3 (partition) this also pins the 29-enzyme sticky set by composition. All 10 names traced to REBASE/NEB/Wikipedia cut sites retrieved this session — sourced values, not code echoes.
- Other tests (M3 partition disjoint+union+count, M4/M7 `All(len==N)`+membership, M5 `count==total−1`+SfiI absent, S1/S2/C1 boundaries) are real, exact, deterministic. No green-washing, no weakened assertions.

### Test-quality gate: **PASS** (after fix)
Filter fixture 11 → 12 tests. Full unfiltered suite: **6569 passed, 0 failed, 0 skipped-of-interest**; `dotnet build` 0 errors (4 pre-existing warnings in unrelated files, untouched).

## Verdict & follow-ups
- Stage A PASS; Stage B PASS-WITH-NOTES (test-quality defect fixed in-session, no code change).
- **End-state: CLEAN.** Algorithm fully functional; tests now lock the complete sourced blunt set.
- Logged in FINDINGS_REGISTER as a test-quality fix.
