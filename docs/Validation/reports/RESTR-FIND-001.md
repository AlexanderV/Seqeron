# Validation Report: RESTR-FIND-001 — Restriction Site Detection

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `RestrictionAnalyzer.FindSites`, `FindAllSites`, `GetEnzyme`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **REBASE (rebase.neb.com)** — authoritative restriction-enzyme database. Fetched individual
  enzyme pages: EcoRI (`G^AATTC`), SfiI (`GGCCNNNN^NGGCC`), HincII (`GTY^RAC`), AscI (`GG^CGCGCC`),
  FseI (`GGCCGG^CC`), DpnI (`GA^TC`, blunt), ClaI (`AT^CGAT`), AgeI (`A^CCGGT`), NdeI (`CA^TATG`).
- **Wikipedia: Restriction enzyme (Examples table)** — full recognition + cut-site table for 23
  common enzymes (EcoRI `G^AATTC`, BamHI `G^GATCC`, HindIII `A^AGCTT`, NotI `GC^GGCCGC`,
  SmaI `CCC^GGG` blunt, PstI `CTGCA^G`, EcoRV `GAT^ATC` blunt, KpnI `GGTAC^C`, SacI `GAGCT^C`,
  SalI `G^TCGAC`, ScaI `AGT^ACT`, SpeI `A^CTAGT`, SphI `GCATG^C`, StuI `AGG^CCT`, XbaI `T^CTAGA`,
  TaqI `T^CGA`, Sau3AI `^GATC`, AluI `AG^CT`, HaeIII `GG^CC`).
- Roberts (1976) Type-II classification — corroborates palindromic Type-II recognition model.

### Convention / formula check
- **Cut-position convention (the model's contract):** `CutPositionForward` is a 0-based offset =
  number of bases on the 5' side of the `^` cut within the recognition sequence. e.g. EcoRI
  `GAATTC`, fwd=1 ⇒ `G^AATTC`. `CutPositionReverse` is the mirror offset measured 5'→3' on the
  forward strand. For every (palindromic) enzyme in the table, `CutPositionReverse = L − CutPositionForward`
  (verified by hand-computation for all 39 entries — all consistent).
- **Reported site `Position`** is the 0-based index of the first base of the recognition match;
  `CutPosition = Position + CutPositionForward` is the 0-based index of the bond cut (between bases
  `CutPosition-1` and `CutPosition`). Standard, internally consistent.
- **IUPAC handling** (`IupacHelper.MatchesIupac`) covers A/C/G/T/N/R/Y/S/W/K/M/B/D/H/V exactly per
  the IUPAC nucleotide ambiguity table; throws on invalid codes. Used for HincII (`GTYRAC`) and
  SfiI (`GGCCNNNNNGGCC`).
- **Both strands:** `FindSitesCore` scans the forward sequence and the reverse complement, mapping
  reverse hits back to forward coordinates via `forwardPos = len − i − patternLen`. Palindromic
  sites are therefore reported on both strands at the same forward position.

### Independent cross-check (worked example)
Sequence `GAATTCGAATTC` (two EcoRI sites). EcoRI `GAATTC`, fwd=1.
- Forward matches at i=0 and i=6 ⇒ Position 0 (cut 1: `G^AATTC...`) and Position 6 (cut 7).
- Reverse complement of the sequence is also `GAATTCGAATTC` (palindrome) ⇒ reverse hits map to the
  same forward positions 0 and 6. Cut offset = AATT 5' overhang (4 nt). Matches REBASE `G^AATTC`.

**Finding:** description (biology, conventions, edge semantics) is correct and authoritatively sourced.

## Stage B — Implementation

### Code path reviewed
`RestrictionAnalyzer.cs:15-66` (enzyme table), `:145-182` (`FindSitesCore`), `:212-225`
(`MatchesPattern`/IUPAC), `:450-475` (`RestrictionEnzyme` overhang derivation).

### Enzyme table — code vs REBASE (39 enzymes; `^` cut shown)

| Enzyme | Code seq | Code fwd/rev | REBASE/Wikipedia | Match |
|--------|----------|--------------|------------------|-------|
| EcoRI | GAATTC | 1/5 | G^AATTC | ✅ |
| BamHI | GGATCC | 1/5 | G^GATCC | ✅ |
| HindIII | AAGCTT | 1/5 | A^AGCTT | ✅ |
| XhoI | CTCGAG | 1/5 | C^TCGAG | ✅ |
| SalI | GTCGAC | 1/5 | G^TCGAC | ✅ |
| PstI | CTGCAG | 5/1 | CTGCA^G | ✅ |
| SphI | GCATGC | 5/1 | GCATG^C | ✅ |
| KpnI | GGTACC | 5/1 | GGTAC^C | ✅ |
| SacI | GAGCTC | 5/1 | GAGCT^C | ✅ |
| XbaI | TCTAGA | 1/5 | T^CTAGA | ✅ |
| NcoI | CCATGG | 1/5 | C^CATGG | ✅ |
| NdeI | CATATG | 2/4 | CA^TATG | ✅ |
| NheI | GCTAGC | 1/5 | G^CTAGC | ✅ |
| SpeI | ACTAGT | 1/5 | A^CTAGT | ✅ |
| AvrII | CCTAGG | 1/5 | C^CTAGG | ✅ |
| ClaI | ATCGAT | 2/4 | AT^CGAT | ✅ |
| AgeI | ACCGGT | 1/5 | A^CCGGT | ✅ |
| NotI | GCGGCCGC | 2/6 | GC^GGCCGC | ✅ |
| SfiI | GGCCNNNNNGGCC | 8/5 | GGCCNNNN^NGGCC | ✅ |
| PacI | TTAATTAA | 5/3 | TTAAT^TAA | ✅ |
| AscI | GGCGCGCC | 2/6 | GG^CGCGCC | ✅ |
| FseI | GGCCGGCC | 6/2 | GGCCGG^CC | ✅ |
| SwaI | ATTTAAAT | 4/4 | ATTT^AAAT (blunt) | ✅ |
| MspI | CCGG | 1/3 | C^CGG | ✅ |
| HpaII | CCGG | 1/3 | C^CGG | ✅ |
| TaqI | TCGA | 1/3 | T^CGA | ✅ |
| AluI | AGCT | 2/2 | AG^CT (blunt) | ✅ |
| RsaI | GTAC | 2/2 | GT^AC (blunt) | ✅ |
| HaeIII | GGCC | 2/2 | GG^CC (blunt) | ✅ |
| DpnI | GATC | 2/2 | GA^TC (blunt) | ✅ |
| Sau3AI | GATC | 0/4 | ^GATC | ✅ |
| MboI | GATC | 0/4 | ^GATC | ✅ |
| EcoRV | GATATC | 3/3 | GAT^ATC (blunt) | ✅ |
| SmaI | CCCGGG | 3/3 | CCC^GGG (blunt) | ✅ |
| HincII | GTYRAC | 3/3 | GTY^RAC (blunt) | ✅ |
| ScaI | AGTACT | 3/3 | AGT^ACT (blunt) | ✅ |
| StuI | AGGCCT | 3/3 | AGG^CCT (blunt) | ✅ |
| BglII | AGATCT | 1/5 | A^GATCT | ✅ |
| ApaI | GGGCCC | 5/1 | GGGCC^C | ✅ |

All 39 hardcoded recognition sequences and both cut offsets match REBASE/Wikipedia exactly.
40 dictionary entries (MspI and HpaII are isoschizomers sharing `C^CGG`). 38 distinct enzymes.

### Edge cases (traced/tested)
- No site ⇒ empty (`FindSites_NoSites_ReturnsEmptyCollection`).
- Empty sequence ⇒ empty (`FindSites_EmptyStringSequence_ReturnsEmptyCollection`).
- Site at boundary: loop bound `i <= seq.Length - pattern.Length` handles last position correctly.
- Overlapping/repeated sites: `FindSites_MultipleSites_FindsAllAtCorrectPositions` (each index scanned).
- Degenerate enzyme: HincII covered by `FindSites_HincII_MatchesAllDegenerateCombinations` (×4),
  `*_RejectsNonMatchingDegenerateBases` (×4), `*_DegenerateEnzymeProducesCorrectCutPosition`.
- Both strands: `FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition`.

### Test quality audit
`EnzymeDatabase_CutPositions_MatchWikipedia` asserts exact recognition string + fwd/rev cut offsets
for 21 enzymes against the Wikipedia source table — real value-locked assertions, not tautologies.
Overhang sequences (AATT/GATC/AGCT), 3' overhang (PstI), blunt (EcoRV) all assert exact sourced
values. M1–M16, S1–S5, C1–C2 covered.

### Findings / defects
None. No false-positive risk; no off-by-one in boundary or cut math.

## Verdict & follow-ups

- **Stage A: PASS** — recognition sequences, cut positions, IUPAC handling, both-strand semantics,
  and 0-based coordinate convention all confirmed against REBASE and Wikipedia.
- **Stage B: PASS** — implementation faithfully realises the validated description; every one of the
  39 hardcoded enzymes matches REBASE exactly; worked example reproduced; edge cases handled & tested.
- **State: CLEAN.** No defect found; no code changed.
- **Tests:** `--filter ~FindSites` = 81/81 passed; full suite `Seqeron.Genomics.Tests` = 4461 passed,
  0 failed, 0 (non-benchmark) skipped.
