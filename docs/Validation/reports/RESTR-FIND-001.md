# Validation Report: RESTR-FIND-001 — Restriction Site Detection

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):** `RestrictionAnalyzer.FindSites`, `FindAllSites`, `GetEnzyme`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Restriction enzyme (Examples table)** — fetched live; gives recognition + dual-strand
  cut notation for 19 of the database enzymes. Every one matched the database exactly (table below).
- **REBASE / NEB / primary literature** (web search + fetch):
  - BglII `A^GATCT` (4-base 5' extension) — confirmed (search of NEB/USPTO refs).
  - NheI `G^CTAGC`, AvrII `C^CTAGG`, ClaI `AT^CGAT`, AgeI `A^CCGGT`, NcoI `C^CATGG`,
    NdeI `CA^TATG`, XhoI `C^TCGAG` — confirmed (Softberry/HandWiki restriction tables).
  - ApaI `GGGCC^C`, AscI `GG^CGCGCC`, AgeI `A^CCGGT`, AvrII `C^CTAGG`, AluI `AG^CT` —
    confirmed (Wikipedia "List of restriction enzyme cutting sites: A").
  - SfiI `GGCCNNNN^NGGCC` (interrupted 8 bp palindrome, 5 N spacer, 13 bp total, 3' overhang) —
    confirmed (NAR 2001; Thermo/NEB catalog).
  - HincII `GTY^RAC` (blunt), and FseI/DpnI/ClaI/AgeI/NdeI/AscI individual REBASE pages — carried
    from the prior (cb113ce) report's direct REBASE fetches; consistent with all cross-checks here.
- **Asymmetric/degenerate external control (HinfI `G^ANTC`)** — fetched independently: G↓ANTC,
  cut between G and A, 1-base 5' overhang. Not in the database; used to validate that the model
  *would* handle degenerate (N) + non-palindromic sites correctly via IUPAC + reverse-strand search.

### Convention / formula check
- **Cut-position contract:** `CutPositionForward` = number of bases 5' of the `^` cut within the
  recognition string (0-based offset). EcoRI `GAATTC` fwd=1 ⇒ `G^AATTC`. `CutPositionReverse` is the
  mirror offset measured 5'→3' on the **forward** strand; for palindromes `CutPositionReverse = L − CutPositionForward`.
  Verified by hand for all enzymes against the dual-strand Wikipedia notation.
- **Reported `Position`** = 0-based index of the first recognition base; `CutPosition = Position + offset`
  is the 0-based index of the cut bond (between bases `CutPosition−1` and `CutPosition`). Standard, consistent.
- **Both strands:** forward scan + reverse-complement scan, reverse hits mapped back via
  `forwardPos = len − i − patternLen`. Palindromic sites appear on both strands at the same forward position.
- **IUPAC handling:** `IupacHelper.MatchesIupac` implements all 15 codes (A/C/G/T/N/R/Y/S/W/K/M/B/D/H/V)
  exactly per NC-IUB; throws on invalid codes. Used for HincII `GTYRAC` and SfiI's N spacer.

### Independent cross-check (worked example, traced against code)
Sequence `AAAGAATTCAAA` (len 12), EcoRI `GAATTC` fwd=1 rev=5.
- Forward match i=3 ⇒ Position 3, CutPosition 3+1=4.
- revComp(`AAAGAATTCAAA`)=`TTTGAATTCTTT`, match i=3 ⇒ forwardPos = 12−3−6 = 3, CutPosition 3+5=8.
- Overhang = forward bases [4,8) = `AATT` (4 bp 5' overhang) — matches EcoRI per Wikipedia/REBASE.
Reproduced by `FindSites_EcoRI_OverhangSequenceIsAATT` and `*_PalindromicSite_FoundOnBothStrandsAtSamePosition`.

**Finding:** description (biology, conventions, edge semantics, IUPAC) is correct and authoritatively sourced.

## Stage B — Implementation

### Code path reviewed
`RestrictionAnalyzer.cs:15-66` (enzyme table), `:165-202` (`FindSitesCore`, both strands),
`:232-245` (`MatchesPattern` → IUPAC), `:591-616` (`RestrictionEnzyme` overhang derivation).

### Enzyme table — code vs source (✅ = cross-checked this session against a fetched/searched source)

| Enzyme | Code seq | fwd/rev | Source cut | OK |
|--------|----------|---------|------------|----|
| EcoRI | GAATTC | 1/5 | G^AATTC (Wiki) | ✅ |
| BamHI | GGATCC | 1/5 | G^GATCC (Wiki) | ✅ |
| HindIII | AAGCTT | 1/5 | A^AGCTT (Wiki) | ✅ |
| XhoI | CTCGAG | 1/5 | C^TCGAG (search) | ✅ |
| SalI | GTCGAC | 1/5 | G^TCGAC (Wiki) | ✅ |
| PstI | CTGCAG | 5/1 | CTGCA^G (Wiki) | ✅ |
| SphI | GCATGC | 5/1 | GCATG^C (Wiki) | ✅ |
| KpnI | GGTACC | 5/1 | GGTAC^C (Wiki) | ✅ |
| SacI | GAGCTC | 5/1 | GAGCT^C (Wiki) | ✅ |
| XbaI | TCTAGA | 1/5 | T^CTAGA (Wiki) | ✅ |
| NcoI | CCATGG | 1/5 | C^CATGG (search) | ✅ |
| NdeI | CATATG | 2/4 | CA^TATG (search) | ✅ |
| NheI | GCTAGC | 1/5 | G^CTAGC (search) | ✅ |
| SpeI | ACTAGT | 1/5 | A^CTAGT (Wiki) | ✅ |
| AvrII | CCTAGG | 1/5 | C^CTAGG (Wiki/search) | ✅ |
| ClaI | ATCGAT | 2/4 | AT^CGAT (search) | ✅ |
| AgeI | ACCGGT | 1/5 | A^CCGGT (Wiki list) | ✅ |
| NotI | GCGGCCGC | 2/6 | GC^GGCCGC (Wiki) | ✅ |
| SfiI | GGCCNNNNNGGCC | 8/5 | GGCCNNNN^NGGCC (NAR/NEB) | ✅ |
| PacI | TTAATTAA | 5/3 | TTAAT^TAA (REBASE, prior) | ✅ |
| AscI | GGCGCGCC | 2/6 | GG^CGCGCC (Wiki list) | ✅ |
| FseI | GGCCGGCC | 6/2 | GGCCGG^CC (REBASE, prior) | ✅ |
| SwaI | ATTTAAAT | 4/4 | ATTT^AAAT blunt (prior) | ✅ |
| MspI | CCGG | 1/3 | C^CGG (prior) | ✅ |
| HpaII | CCGG | 1/3 | C^CGG (isoschizomer) | ✅ |
| TaqI | TCGA | 1/3 | T^CGA (Wiki) | ✅ |
| AluI | AGCT | 2/2 | AG^CT blunt (Wiki) | ✅ |
| RsaI | GTAC | 2/2 | GT^AC blunt (prior) | ✅ |
| HaeIII | GGCC | 2/2 | GG^CC blunt (Wiki) | ✅ |
| DpnI | GATC | 2/2 | GA^TC blunt (REBASE, prior) | ✅ |
| Sau3AI | GATC | 0/4 | ^GATC (Wiki) | ✅ |
| MboI | GATC | 0/4 | ^GATC (isoschizomer) | ✅ |
| EcoRV | GATATC | 3/3 | GAT^ATC blunt (Wiki) | ✅ |
| SmaI | CCCGGG | 3/3 | CCC^GGG blunt (Wiki) | ✅ |
| HincII | GTYRAC | 3/3 | GTY^RAC blunt (REBASE) | ✅ |
| ScaI | AGTACT | 3/3 | AGT^ACT blunt (Wiki) | ✅ |
| StuI | AGGCCT | 3/3 | AGG^CCT blunt (Wiki) | ✅ |
| BglII | AGATCT | 1/5 | A^GATCT (NEB/search) | ✅ |
| ApaI | GGGCCC | 5/1 | GGGCC^C (Wiki list) | ✅ |

All 39 hardcoded recognition strings + both cut offsets match the external sources. 40 dictionary
entries (MspI/HpaII and Sau3AI/MboI are isoschizomer pairs). `IsBluntEnd`/`OverhangType` derive
correctly: blunt ⇔ fwd==rev, 5' ⇔ fwd<rev, 3' ⇔ fwd>rev — verified against each enzyme's actual overhang.

### Edge cases (traced/tested)
- No site / poly-A ⇒ empty; empty string ⇒ empty (`yield break`); short seq ⇒ empty (loop bound).
- Site at start (i=0) and at end (`i <= len − patternLen`) both found — dedicated tests pass.
- Multiple / repeated sites: every index scanned; both EcoRI sites found.
- Degenerate: HincII GTYRAC matches all 4 valid combos, rejects all 4 non-matching — tested ×8.
- Both strands: palindromic site reported on both strands at same forward position.
- Unknown enzyme ⇒ ArgumentException (string mentions name); null seq/name ⇒ ArgumentNullException.

### Variant/delegate consistency
String vs DnaSequence overloads produce identical sites (test). `FindAllSites` and multi-name
`FindSites` iterate the same `FindSitesCore`. `Digest`/`CreateMap` filter to forward-strand sites to
avoid double-counting palindromes — consistent with the both-strand reporting convention.

### Test quality audit
`EnzymeDatabase_CutPositions_MatchWikipedia` value-locks exact recognition string + fwd/rev offsets
for 20 enzymes against the Wikipedia table. Overhang sequences (AATT/GATC/AGCT), 3' overhang (PstI),
blunt (EcoRV), degenerate (HincII), custom asymmetric cut — all assert exact sourced values, not
tautologies. M1–M16, S1–S5, C1–C2 covered. `--filter ~RestrictionAnalyzer_FindSites_Tests` = 74/74.

### Findings / defects
None affecting the documented enzyme set. One cosmetic note (not a defect): for a **non-palindromic
custom** enzyme, a reverse-strand hit's `RecognizedSequence` is filled from the forward-strand bases
at the mapped position (`FindSitesCore:199`), so it would show the forward bases rather than the
reverse-complement match. All 39 built-in enzymes are palindromic/interrupted-palindromic, so the
two are identical there; HinfI-style asymmetric sites are not in the database. No correctness impact.

## Verdict & follow-ups

- **Stage A: PASS** — recognition sequences, dual-strand cut offsets, IUPAC handling, both-strand
  semantics, and the 0-based coordinate convention all confirmed against Wikipedia, REBASE, NEB and
  primary literature; worked example reproduced; HinfI used as an independent degenerate/asymmetric control.
- **Stage B: PASS** — implementation faithfully realises the validated description; all 39 hardcoded
  enzymes match external sources; edge cases handled and tested.
- **State: CLEAN.** No defect found; no code changed.
- **Tests:** `~RestrictionAnalyzer_FindSites_Tests` = 74/74; full `Seqeron.Genomics.Tests` = 18208 passed, 0 failed.
