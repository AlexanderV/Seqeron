# B08 — restriction slice of the round-4 worklist (H6, H7, H13, H24, H25, H27, H28, H30–H36, H45, H52, H54-restriction)

Companion to `B08.md`. The round-4 worklist in `B08.md` is implemented by two sessions in parallel:
the orchestrator takes the CRISPR / Azimuth items, this session takes the restriction-enzyme items.
Entries here are numbered **F100+**; the orchestrator uses F91–F99. `B08.md` itself is touched only
to tick this slice's own `- ☐ **Hnn` lines.

**Owned files:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer*.cs`, their
tests / TestSpecs / Evidence, `docs/algorithms/MolTools/Restriction_*.md`, and the restriction
wrappers in `src/Seqeron/Mcp/Seqeron.Mcp.MolTools/MolToolsTools.cs`.
**Frozen, not touched:** `MolTools/Resources/rebase_enzymes.tsv`,
`scripts/rebase/generate_rebase_enzymes.py` (licence question pending).

## Oracle environment (this session's container)

| Oracle | Version | How it was obtained |
|---|---|---|
| Biopython `Bio.Restriction` | 1.88 | `pip install biopython==1.88` |
| EMBOSS `restrict` | 6.6.0.0 | `apt-get install -y emboss` |
| REBASE `allenz.104` | — | only via the frozen `rebase_enzymes.tsv` already in the repo; `rebase.neb.com` is not used (TLS verification fails) |

EMBOSS ships no `embossre.enz` (REBASE is licensed), so `restrict` was driven with an enzyme file
generated from Biopython 1.88's `Restriction_Dictionary`
(`c1 = fst5`, `c2 = size + fst3`, `c3 = scd5`, `c4 = size + scd3`, `blunt = 1` iff `c1 == c2`,
`ncuts` 2 or 4; the `blunt` column is an **integer**, not `Y`/`N`, or `restrict` silently loads
nothing). The file must sit at `$EMBOSS_DATA/REBASE/embossre.enz` beside an (empty) `embossre.equ`.
**Calibration** on `GGAATTCCGGTCTCAAAATGGCCATTTTTACGTACGTAAA`: EMBOSS's `5prime` column equals
Biopython's 0-based top-strand bond for every forward hit —
EcoRI 2/2, BsaI 15/15, HauII 35/35, SnaBI 31,35 / 31,35, BsiWI 31/31, AciI none/none.
The only extra EMBOSS row is HauII on the minus strand, which is F101 (H7) itself.

## Confirmed fixes

### F100 — RESTR-FIND-001 (H6): the circular search skipped every molecule shorter than the recognition sequence

- **Defect.** `RestrictionAnalyzer.cs` applied the linear `seq.Length < pattern.Length` guard
  *before* the topology dispatch, so `FindSitesCircular`'s wrap-around extension was never built.
  On a circle a window longer than `n` still exists — it reads past the origin and back in.
- **Reproduced (Biopython 1.88, `Enzyme.search(Seq(m), linear=False)`).** Seqeron returned
  **0 sites for all 10** probe rows; Biopython finds a site in 9 of them:

  | molecule | enzyme | site (len) | Biopython (1-based) | reduced, 0-based | Seqeron before |
  |---|---|---|---|---|---|
  | `ACGT` | BsiWI, Pfl23II, PspLI, SplI | `CGTACG` (6) | `[3]` | 2 | **0 sites** |
  | `ACATGA` | MslI, RseI, SmiMI | `CAYNNNNRTG` (10) | `[1]` | 0 | **0 sites** |
  | `ACATGA` | DrdV | `CATGNAC` (7) | `[13]` | 0 | **0 sites** |
  | `TATTGCAAAA` | AlfI | `GCANNNNNNTGC` (12) | `[5, 19]` | 4, 8 | **0 sites** |

  (Biopython's circular `_drop` normalises only a prefix/suffix of its own sorted result list, so it
  can return 13 for a 6-mer and 19 for a 10-mer; both oracles are compared reduced modulo `n`.)
  An existing test, `FindSites_Circular_MoleculeShorterThanSite_ReportsNothing`, **pinned the
  defect**, arguing in its XML doc that such a reading "wraps a short molecule more than once".
  That premise is false: `ATCGAT` on the 4-bp circle `GATC` starts at index 1 and needs two bases
  past the origin — exactly Biopython's own one-pass `data + data[1:size]` extension.
- **Fix.** The length guard is scoped to `MoleculeTopology.Linear`; on a circle only an empty
  molecule has no window. `FindSitesCircular` builds its extension with a new `WrapExtend` helper
  that repeats the molecule to `n + len - 1` bases instead of slicing `seq[..(len - 1)]` (which
  threw for `n < len - 1`). For every `n >= len - 1` `WrapExtend` returns the identical string the
  old expression did, so the change is provably a no-op at every length that already worked.
- **Post-fix oracle agreement.** All 9 rows above now match Biopython exactly. Differential over
  **all 754 REBASE enzymes × 84 molecules = 63 336 rows** (every length 1–13 in five compositions,
  plus 4 random molecules each at 20/25/40/97/200, plus the two H6 molecules), comparing the
  circular `Digest` top-strand cut set against `Enzyme.search(..., linear=False)` reduced mod `n`:
  **0 divergences at every length that can hold the site** (0 missing, 0 mismatched).
- **Declared divergence (superset, 44 of the 63 336 rows, all at `n < len`).** Biopython's
  `data[1:size]` slice truncates to the molecule, so it can only reach start indices
  `0 .. n + min(n, len-1) - len` and never tests the rest. Those windows are physically present on
  the circle. Example locked by a test: `SnaBI` (`TACGTA`) on circular `ACGT` is the window at start
  3; `Bio.Restriction.SnaBI.search(Seq("ACGT"), linear=False) == []`, Seqeron reports it with cut 2.
  EMBOSS `restrict -plasmid` reports nothing at all for `n < len` (verified on all three H6
  molecules), so it does not arbitrate; it does agree with Biopython on an ordinary wrap
  (`ACGAAAAAAACGT` + BsiWI → EMBOSS `5prime = 11`, Biopython `[12]` 1-based).
- **Tests.** `RestrictionAnalyzer_FindSites_Tests.cs`:
  `FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite` (replaces the defect-pinning
  test, locks Biopython's `ClaI.search(Seq("GATC"), linear=False) == [4]`),
  `Digest_Circular_MoleculeShorterThanSite_MatchesBiopython` (9 cases, the table above),
  `FindSites_Circular_ShortMolecule_ReportsWindowsBiopythonTruncatesAway` (the declared superset),
  `FindSites_Circular_EmptyMolecule_ReportsNothing` (the guard that remains).
- **Commit.** see `git log --grep "B08 H6"`.

### F101 — RESTR-DIGEST-001 (H7): a window reading as the site on both strands lost one of its two cuts

- **Defect.** `CollectForwardCuts` (the cut collector behind `Digest`, `GetDigestSummary` and
  `PartialDigest`) and the mirror filter in `CreateMapCore` dropped a reverse-strand record whenever
  the same window also read as the site on the forward strand. That is Biopython's rule, and it is
  exact only while the two readings cut the same bonds. `FindSites` reported both readings all along.
- **The judgement call: EMBOSS is followed, not Biopython.** Re-derived from both tools' source in
  this container:
  - *Biopython 1.88* `Restriction.py`: `Palindromic._search` builds its result from `_modify` alone,
    and `NonPalindromic._search` loops `for start, group in iterator:` and takes `modif(start)`
    whenever `group(s)` — the **forward** named group of the `(?=fwd)|(?=rev)` compsite — matched,
    falling back to `_rev_modify` only otherwise. So one reading per *window*, forward preferred,
    whatever the cut.
  - *EMBOSS 6.6.0* `nucleus/embpat.c`, `embPatRestrictScan` (source obtained with
    `apt-get source emboss`, 6.6.0+dfsg-12ubuntu2): the function runs **two independent scan loops**
    — `forward = ajTrue` over `substr`/`binstr`, then `forward = ajFalse` over `revstr`/`binrev` —
    each pushing hits, and afterwards keeps one hit per **distinct `cut1`**
    (`if(m->cut1 != v) { push; v = m->cut1; }`), preferring the forward hit when two coincide.
    So EMBOSS dedupes by *cut*, not by window.
  - *Chemistry.* A palindromic recognition window is bound by the enzyme's recognition domain in
    either orientation; for a distal cutter the cleavage then falls on opposite sides of the window,
    so both cuts are real. Biopython's one-reading-per-window rule is an optimisation that happens to
    be lossless only for a symmetric cutter — which is the case for 562 of its 564 palindromic
    enzymes, and for all 39 of this library's built-ins. **EMBOSS's rule is the one that matches the
    chemistry and it is adopted**: a reverse record is dropped only when its cut set duplicates the
    forward record's (`ReverseReadingCutsAlike`).
- **Reproduced (EMBOSS 6.6.0 `restrict`, datafile as described above).** Over 170 rows
  (34 molecules × the 5 enzymes, linear): **44 divergent, every one a cut Seqeron was missing, never
  an extra.** Per enzyme: MspJI 19, SgeI 19, LpnPI 2, SgrTI 2, HauII 2.
  Worked single-window case — molecule `A×16 + CAAG + A×20` (40 nt), one `CNNG` window at 0-based 16,
  positioned so both readings' cuts are internal:

  | | cuts (0-based top-strand bond) | fragments |
  |---|---|---|
  | Seqeron before | `29` | 2: 29, 11 |
  | Biopython 1.88 `SgeI.search(…, linear=True)` | `[30]` 1-based = `29` | — |
  | EMBOSS `restrict` | `-` row `5prime 3, 3prime 7`; `+` row `5prime 29, 3prime 33` → `3, 29` | 3: 3, 26, 11 |
  | **Seqeron after** | **`3, 29`** | **3: 3, 26, 11** |

  `CreateMap` on the same molecule: before `TotalSites = 1`, `UniqueCutters = [SgeI]`; after
  `TotalSites = 2`, `UniqueCutters = []` — an enzyme that binds once but cuts twice is not a single
  cutter, which is what the MCP description's stated use ("plan cloning around single-cutter
  enzymes") needs.
- **Post-fix oracle agreement.**
  - vs **EMBOSS**, the 5 enzymes, 118 molecules = **590 rows → 0 divergences** (was 44/170).
  - vs **Biopython**, all 754 enzymes × 118 molecules = **88 972 rows → 76 divergences**, involving
    **exactly the 5 enzymes** (MspJI 31, SgeI 31, LpnPI 5, SgrTI 5, HauII 4) and in **all 76 rows
    Seqeron is a strict superset** — it never misses a cut Biopython reports. The other **749
    enzymes keep the F59 agreement with Biopython exactly**.
  - circular: 63 336 rows; the rows that changed relative to the pre-H7 run are exactly the same 5
    enzymes, and there is still no row where Seqeron misses a Biopython cut.
- **Tests.** New region *Both-strand readings of one window* in `RestrictionAnalyzer_Digest_Tests.cs`:
  `Digest_WindowReadingOnBothStrands_CutsAtBothReadings` (SgeI + MspJI, the table above),
  `CreateMap_WindowReadingOnBothStrands_CountsBothCleavages`,
  `Digest_PalindromicSymmetricCutter_StillCollapsesTheReverseReading` (EcoRI unchanged).
  `Digest_ByName_EveryRebaseEnzyme_CutsAtBiopythonSearchPositions` now takes its expectation for
  those 5 enzymes from a **new EMBOSS fixture**,
  `tests/.../TestData/Restriction/emboss_restrict_bothstrand.tsv`, generated by the **new**
  `scripts/rebase/generate_emboss_restrict_oracle.py` over the same molecules as
  `biopython_rebase_search.tsv`; the test asserts the EMBOSS rows are a superset of the Biopython
  rows and that exactly 5 enzymes take the EMBOSS oracle. **The frozen
  `scripts/rebase/generate_rebase_enzymes.py` and `Resources/rebase_enzymes.tsv` were not touched.**

### F102 — RESTR-DIGEST-001 (H45): the region comment's count of affected enzymes was wrong

- **Defect.** `RestrictionAnalyzer.cs` said the rule "differs for the 4 palindromic distal cutters
  (HauII, SgeI, …)". That is the stated justification for F101's rule, and "4" matches neither
  reading of the set.
- **Computed over the 754 REBASE enzymes with a defined cut** (both in a scratch harness from
  Biopython's `Restriction_Dictionary` and, independently, in the new test from the library's own
  table — the two agree):
  - **564** have a palindromic recognition *string*;
  - **9** of those have `fst5 + fst3 != 0`, i.e. an asymmetric first cut pair:
    AlfI, BdaI, BplI, FalI, HauII, NmeDI, RdeGBIII, SgeI, UcoMSI;
  - **2** of those 9 have a cut *set* that differs under the reverse reading: **HauII, SgeI**. The
    other 7 are two-sided (Type IIB) and their two cut pairs swap into each other, so the set
    survives;
  - across the whole set **166** enzymes have a differing cut set, and **exactly 5** of those admit a
    window matching the pattern on both strands: **HauII** `TGGCCA`, **SgeI** `CNNG`, **MspJI**
    `CNNR`, **LpnPI** `CCDG`, **SgrTI** `CCDS`.
- **Fix.** The region comment now carries all of these numbers, states both oracles' rules with their
  source locations, says which is followed and why, and records the differential result. The same
  stale "4 palindromic distal cutters" wording in
  `RestrictionAnalyzer_RebaseSet_Tests.cs` was corrected with it.
- **Tests.** `ReverseReadingScope_OverTheWholeRebaseSet_MatchesTheDocumentedCounts` (754 / 564 / the
  9 names / the 2 names / 166) and `ReverseReadingScope_ShortSites_AffectExactlyFiveEnzymes`
  (enumerates every concrete window of every site of at most 8 positions and returns exactly
  HauII, LpnPI, MspJI, SgeI, SgrTI) — computed, not asserted from the comment.

## Test results

| Scope | Result |
|---|---|
| `dotnet build Seqeron.sln -c Debug` | 0 warnings, 0 errors |
| targeted tier after F100 (`~Restriction`, fast filter) | 621 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F100 | 471 passed / 0 failed |
| targeted tier after F101/F102 (`~Restriction\|~Digest\|~Enzyme\|~MolTools`, fast filter) | 3131 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F101/F102 | 471 passed / 0 failed |

## Leftovers

(pending — the slice is not finished)
