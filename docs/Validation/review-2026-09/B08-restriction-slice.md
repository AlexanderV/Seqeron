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

### F103 — RESTR-FIND-001 (H13): `GetEnzyme` could not resolve `CviKI-1`, the enzyme's real REBASE name

- **Defect.** `GetEnzyme` did two exact dictionary lookups against a table keyed by Biopython's
  class-name spelling. Biopython replaces a REBASE hyphen with an underscore because a Python class
  name cannot contain one, so the shipped enzyme is keyed `CviKI_1` while REBASE and NEB (R0710)
  spell it **`CviKI-1`**. No trimming either, although `MethylaseSets.FromName` trims.
- **Reproduced.** `GetEnzyme("CviKI-1")` → `null`; `GetEnzyme("cviki-1")` → `null`;
  `GetEnzyme(" EcoRI")` → `null`; `GetEnzyme("EcoRI ")` → `null`. Every consuming member then failed:
  `FindSites(seq, "CviKI-1")`, `GetOverhangSequence("CviKI-1")` and `AreCompatible("CviKI-1", "AluI")`
  all threw `ArgumentException: Unknown enzyme: CviKI-1`. Over MCP that is
  `find_restriction_sites` erroring and `get_enzyme` answering "not a REBASE enzyme with a defined
  cut" for a real NEB catalogue enzyme.
- **Fix, in the consuming code only — the frozen `rebase_enzymes.tsv` is not touched.**
  `GetEnzyme` trims the name and, when the trimmed name contains a hyphen and is otherwise unknown,
  retries with `-` replaced by `_`. The table's own spelling is still what the resolved enzyme
  reports back. Exactly one shipped enzyme is affected: of the 5 underscore-keyed names in
  Biopython's dictionary (`CviKI_1`, `Hpy99XIV_mut1`, `HpyAXVI_mut1`, `HpyAXVI_mut2`,
  `HpyUM032XIII_mut1`) only `CviKI_1` has a defined cut and so appears in the frozen table.
- **Post-fix oracle agreement (Biopython 1.88).** `CviKI-1`, `CviKI_1`, `cviki-1` and `CVIKI_1` all
  resolve to the entry whose recognition sequence is `RGCY`; `GetOverhangSequence` → `""`, matching
  `Bio.Restriction.CviKI_1.ovhgseq == ''` (`is_blunt() == True`); the digest of `AGCTAGCT` cuts at
  the 0-based bonds **2 and 6**, matching `CviKI_1.search(Seq("AGCTAGCT"), linear=True) == [3, 7]`;
  `AreCompatible("CviKI-1", "AluI")` → `true` (both blunt). A genuinely unknown name, and a hyphen
  in a position the table does not use (`Eco-RI`), still answer `null`.
- **Tests.** `GetEnzyme_RebaseHyphenSpelling_ResolvesToTheSameEnzyme` (4 spellings × 5 assertions,
  Biopython numbers above), `GetEnzyme_SurroundingWhitespace_IsIgnored`,
  `GetEnzyme_StillNullForAnUnresolvableName` — in `RestrictionAnalyzer_FindSites_Tests.cs`.

### F104 — RESTR-FIND-001/RESTR-DIGEST-001 (H24): one enzyme-name argument shape, four meanings

- **Defect, reproduced.** For the same `params string[] enzymeNames` shape:

  | call | before | after |
  |---|---|---|
  | `FindSites(seq, Linear, [])` | **0 sites** | `ArgumentException` |
  | `FindSites(seq, Linear, methylases, [])` | `ArgumentException` | unchanged |
  | `Digest(seq, [])` / `PartialDigest` / `GetDigestSummary` | `ArgumentException` | unchanged |
  | `CreateMap(seq)` (no names) | every common enzyme | unchanged, now the single **documented** exception |
  | `FindSites(seq, "NoSuchEnzyme1", "EcoRI")` | **built fine, threw on `.ToList()`** | `ArgumentException` at the call |
  | `Digest(seq, "NoSuchEnzyme1")` | **built fine, threw on `.ToList()`** | `ArgumentException` at the call |

  The lazy rejection contradicted the class's own name-resolution contract and the `AreCompatible`
  remark, which both state that a misspelled name is rejected eagerly — and the raw-string
  `FindSites(string, string, MoleculeTopology)` already was eager ("an unknown enzyme name is still
  rejected, eagerly, before the sequence is looked at").
- **Fix.** `FindSitesCore(DnaSequence, string[], MoleculeTopology)` is no longer an iterator: it
  resolves the names through a new `RequireEnzymeNames` helper (null array → `ArgumentNullException`,
  empty array → `ArgumentException`, each name through `RequireEnzyme`) and then streams. The same
  rule now holds for the methylation overload of `FindSites`, `Digest`, `PartialDigest` and
  `GetDigestSummary`, whose `ResolveEnzymes` is materialised rather than lazily projected.
  `CreateMap`'s name overload keeps "empty = every common enzyme" because that is its parameter's
  documented meaning, and `RequireEnzymeNames`' doc names it as the one deliberate exception.
- **Documentation.** `RequireEnzyme` now states the whole contract in one place — a null **or empty**
  name is `ArgumentNullException`, an unresolvable name is `ArgumentException`, both eager — and the
  7 `<exception cref="ArgumentNullException">` tags on name-taking public members that said only
  "is null" now say "is null or empty", which closes the undocumented `GetOverhangSequence("")` and
  `AreCompatible("", …)` throws. The exception *type* for an empty name is left as
  `ArgumentNullException`: it is pinned by an existing test across the whole restriction surface, and
  changing it is H19's subject, not this item's.
- **Tests.** `EmptyEnzymeNameArray_IsAnErrorEverywhereExceptCreateMap` (7 members throw; `CreateMap`
  still considers all 39 common enzymes) and `UnknownEnzymeNameInAnArray_IsRejectedAtTheCall`
  (6 members throw without enumerating); the existing
  `FindSites_EmptyOrNullEnzymeName_ThrowsArgumentNullException` gained the eagerness assertion.

### F105 — RESTR-FIND-001 (H25): `FindPatternOccurrences` threw where its own `<param>` doc promised an empty result

- **Defect, reproduced.** `<param name="sequence">` promised "a null or empty sequence yields no
  occurrences"; the body began with `ArgumentNullException.ThrowIfNull(sequence)`, and the method's
  own `<exception cref="ArgumentNullException">` tag documented that throw — so two doc statements
  contradicted each other and one contradicted the code. Measured:
  `FindPatternOccurrences(null, "GAATTC")` → `ArgumentNullException`, while the sibling
  `FindSites((string)null, "EcoRI")` → `0` and `FindSites("", "EcoRI")` → `0`.
- **Fix.** The code follows the `<param>` doc and the sibling: `sequence` is `string?` and a null or
  empty sequence yields no occurrences (`string.IsNullOrEmpty`, exactly the sibling's expression).
  The stale `<exception>` tag is gone. The **pattern** and `startIndex` are still validated eagerly,
  before the sequence is looked at — `FindPatternOccurrences(null, "")`,
  `FindPatternOccurrences(null, "GAZTC")` and `FindPatternOccurrences(null, "GAATTC", -1)` all still
  throw. The nullable annotation is source- and binary-compatible.
  This also matters for F106 (H27): the MCP wrappers route an IUPAC-bearing sequence through this
  method.
- **Tests.** `FindPatternOccurrences_NullOrEmptySequence_YieldsNoOccurrences` (new);
  `FindPatternOccurrences_InvalidArguments_Throw` updated — its null-sequence row now asserts the
  empty result and its doc records the supersession.

### F106 — RESTR-FILTER-001 (H30): the compatibility "shared end" is the queried enzyme's own end

- **Defect, reproduced.** `AreCompatible`'s `out compatibleEnd` is set to the *first* enzyme's
  overhang window, so the value is order-dependent, while `FindCompatibleEnzymes`' `<returns>` called
  it "the shared end" and the MCP `enzymes_compatible` description said it "echoes the shared end".
  Measured: `GetOverhangSequence("BsmI") == "CN"`, `GetOverhangSequence("BsrI") == "GN"`;
  `FindCompatibleEnzymes("BsmI", Rebase)` reports BsrI with `CN` and
  `FindCompatibleEnzymes("BsrI", Rebase)` reports BsmI with `GN`. **BsrI never leaves `CN`.**
- **Scope, computed over all 754 REBASE enzymes.** Of the sticky compatible pairs, **18 541** anneal
  because the two windows intersect directly, and exactly **12** anneal only through the reverse
  complement — all of them `GN` against `CN`: Bse1I/BsmI, Bse1I/LmnI, Bse1I/Mva1269I, Bse1I/PctI,
  BseNI/BsmI, BseNI/LmnI, BseNI/Mva1269I, BseNI/PctI, BsmI/BsrI, BsrI/LmnI, BsrI/Mva1269I, BsrI/PctI.
  For those 12 **no single shared string exists**, so no value of one string can be truthful.
- **Fix.** The value is left as it is — the queried/first enzyme's own end, which is the useful one —
  and every place that described it is corrected to say so, naming the 12 pairs. Additively, the MCP
  output records gain `end2`, the *partner's* own window (`"blunt"` for a blunt pair), so a caller
  sees both ends of a pair in one call instead of having to call `get_enzyme` again:
  `EnzymeCompatibilityResult(bool Compatible, string? CompatibleEnd, string? End2)` and
  `EnzymeCompatibilityPair(string Enzyme1, string Enzyme2, string CompatibleEnd, string? End2)`, both
  with an optional parameter so existing construction sites keep compiling.
- **Post-fix.** `enzymes_compatible("BsmI", "BsrI")` → `compatibleEnd "CN"`, `end2 "GN"`; the reversed
  call → `"GN"` / `"CN"`; `("EcoRV", "SmaI")` → `"blunt"` / `"blunt"`; `("EcoRI", "BamHI")` →
  `compatible false` with both ends `null`.
- **Tests.** `FindCompatibleEnzymes_ReverseComplementOnlyPair_ReportsTheQueriedEnzymesOwnEnd`
  (`RestrictionAnalyzer_FindSites_Tests.cs`) and `EnzymesCompatible_ReverseComplementOnlyPair_ReportsBothEnds`
  (`EnzymeSetToolsTests.cs`). Two existing MCP tests that built their expectation with the 3-argument
  constructor were updated to carry `end2`.
- **Docs.** `FindCompatibleEnzymes(string)`'s `<returns>` and summary; MCP descriptions of
  `enzymes_compatible` and `compatible_enzymes`; `docs/mcp/tools/moltools/enzymes_compatible.{md,mcp.json}`
  (output schema + a worked BsmI/BsrI example) and `compatible_enzymes.{md,mcp.json}`.

### F107 — RESTR-FIND-001 (H31): `FindSites` promised a global ascending order it never had

- **Defect, reproduced.** The `<returns>` doc said "forward-strand hits first, then reverse-strand
  hits, each in ascending `Position` order". The reverse-strand scan runs 5'→3' along the reverse
  complement and `Position = n - i - len` *falls* as the scan index rises. Measured, EcoRI on
  `GAATTCAAGAATTCAAGAATTC`: `p0F p8F p16F p16R p8R p0R`. On a circle it is not monotone at all —
  Sau3AI on circular `GATCGATCGATC` gives positions `0, 4, 8, 8, 4, 0`.
- **Second defect in the same block.** It claimed `FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI")`
  "reports the single site at position 4". One *window* survives the linear cut-feasibility filter,
  but it is reported as **two records** (forward cut 4, reverse cut 8) by the
  two-records-per-palindromic-site rule the same doc states three paragraphs earlier.
- **Fix.** The `<returns>` doc now states the real order (forward ascending, then reverse descending
  on a linear molecule, non-monotone on a circle), says to sort by `Position` for a global order, and
  notes that `CrisprDesigner.FindOffTargets` behaves the same way; the Sau3AI example counts records
  and names the top-strand cut that Biopython reports (`Sau3AI.search(Seq("GATCGATCGATC")) == [5]`).
- **Tests.** `FindSites_ReverseStrandHits_ComeOutInDescendingPositionOrder` (both topologies, the
  exact sequences above) and `FindSites_Sau3AI_OneSurvivingWindow_IsTwoRecords`.

### F108 — RESTR-DIGEST-001 (H32): `CreateMap`'s two overloads documented `TotalSites` contradictorily

- **Defect, reproduced.** One overload said `TotalSites`/`UniqueCutters` "count forward-strand
  records only"; the other said "physical sites"; the code implements neither phrasing exactly and
  the MCP description took the first one, so its stated use — "plan cloning around single-cutter
  enzymes" — read wrong for every Type IIS enzyme, i.e. every Golden Gate enzyme. Measured on
  `A×10 + GGTCTC + A×10 + GAGACC + A×10`: one forward-strand record and one reverse-strand record, so
  the code gives `TotalSites = 2` and `UniqueCutters = []` where "forward records only" implies 1 and
  a unique cutter.
- **Fix.** One statement, in all four places (both XML docs, the MCP tool description, and
  `docs/mcp/tools/moltools/restriction_map.{md,mcp.json}`): the two members count **cleavages** — the
  rule F101 settled — so a palindromic site with symmetric cuts counts once, a non-palindromic
  enzyme's reverse-strand sites count on their own, a two-sided (Type IIB) site's two cut records
  count once, and a window admitting two different cleavages counts twice. `UniqueCutters` therefore
  means "cleaves exactly once", which is what the cloning use needs.
- **Checked unchanged by F101:** the measured 754-enzyme example in `restriction_map.md`
  (`AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA` → 14 enzymes / 740 non-cutters / `totalSites 41` /
  `uniqueCutters ["SetI"]`) still holds, because on that molecule the five H7 enzymes' reverse
  readings cut outside the molecule and are dropped by the linear feasibility filter.
- **Tests.** `CreateMap_TypeIISEnzyme_CountsTheReverseStrandSiteAndIsNotAUniqueCutter` and the MCP
  `RestrictionMap_TypeIISEnzyme_CountsBothStrandsCleavages`.

### F109 — RESTR-FILTER-001 (H33): the cited Biopython method was misdescribed

- **Defect, reproduced.** `FindCompatibleEnzymes(string)`'s remark claimed the reference "filters its
  candidates with `x.is_5overhang()`, so it answers the empty list for every blunt and every
  3'-overhang enzyme". Biopython 1.88 has **three** `compatible_end` implementations, read from the
  installed source: `Blunt.compatible_end` returns every `x.is_blunt()`, `Ov5.compatible_end` every
  `x.is_5overhang() and x % cls`, `Ov3.compatible_end` every `x.is_3overhang() and x % cls`.
  Measured: `EcoRV.compatible_end()` → **137** (every blunt enzyme), `PstI.compatible_end()` → **13**,
  `SacI.compatible_end()` → **9**, `EcoRI.compatible_end()` → **51**; each includes the queried
  enzyme itself.
- **Second defect.** The summary said "the enzyme itself included" without qualification; that is
  false when the query lies outside the chosen set. Measured:
  `FindCompatibleEnzymes("BsaI", EnzymeSet.Common)` does not list BsaI (not one of the 39), while the
  REBASE query does.
- **Also recorded:** Biopython's `compatible_end` ignores its own `batch` argument and always
  iterates `AllEnzymes`; this library honours its `EnzymeSet`.
- **Fix.** The remark now states all three reference implementations with their predicates and the
  four measured counts, and the summary scopes self-inclusion to membership of the chosen set.
- **Tests.** `FindCompatibleEnzymes_BluntAndThreePrimeEnzymes_AreNotEmpty` (137 blunt REBASE cutters,
  cross-checking Biopython's own `EcoRV.compatible_end()` count, plus a non-empty 3'-overhang result
  containing NsiI, which is in Biopython's `PstI.compatible_end()` list) and
  `FindCompatibleEnzymes_QueryOutsideTheChosenSet_DoesNotListItself`.

### F110 — docs (H34, H35, H36): algorithm-doc invariants brought to the shipped code

- **H34, reproduced.** `Restriction_Site_Detection.md` INV-01 asserted
  `0 <= Position <= Length - recognitionLength` unconditionally, while §2.3 and INV-07 do scope their
  claims by topology. Measured: EcoRI on circular `TTCAAAAAAAGAA` (13 nt, 6-nt site) reports
  `Position = 10`, where `Length - 6 = 7`; Biopython agrees
  (`EcoRI.search(..., linear=False) == [12]`). INV-01 is now scoped to Linear and a new **INV-01c**
  gives the circular bound `0 <= Position <= Length - 1`.
  **Also fixed in the same table, which the audit did not list but F100 invalidated:** INV-08 ("No
  site is reported when `Length < RecognitionLength`") is now scoped to Linear, with a new
  **INV-08c** for the circular case that F100 created.
- **H35, reproduced.** Three places said a circle with `k` cuts gives `k` fragments. Measured on a
  circular 20-mer with EcoRI: k=0 → **1** fragment (20), k=1 → 1 (20), k=2 → 2 (10, 10) — only k=0
  diverges, and the behaviour is deliberate and already pinned by
  `RestrictionAnalyzer_Digest_Tests.cs`, with the carve-out stated correctly on `DigestCircular`'s own
  XML remark. So four statements of one rule disagreed. `Restriction_Digest_Simulation.md` INV-01 is
  now scoped to Linear, INV-05 carries the `k = 0` carve-out with the measured numbers, and §2.1's
  bare `k + 1` is qualified.
- **H36, reproduced.** `Restriction_Enzyme_Filtering.md` INV-07 claimed reflexivity "for every
  enzyme whose overhang is determined (all 39 built-ins)", while the same document's §2.2 names
  `SfiI` → `NNN` explicitly. Measured: of the 39 built-ins exactly **one**, `SfiI`
  (`GGCCNNNNNGGCC`), has a non-ACGT overhang window. Reflexivity survives; the parenthetical did not,
  and is replaced by the fact.
- **Also brought up to date (F101's consequence, not a separate audit item):**
  `Restriction_Digest_Simulation.md` §2.5 step 3 and §4, which stated the Biopython
  one-reading-per-window rule, and `docs/Evidence/RESTR-DIGEST-001-Evidence.md` §"Fix", which stated
  the same; plus `docs/Evidence/RESTR-FIND-001-Evidence.md` §3, whose "Deliberate divergences"
  paragraph asserted that no site is reported when `Length < RecognitionLength` on either topology and
  argued, wrongly, that such a reading wraps the molecule more than once. Each now states the shipped
  rule with the measured differential.

### F111 — RESTR-FILTER-001 (H52): the three isoschizomer queries' null guards were untested

- **Defect.** `FindIsoschizomers` / `FindEquischizomers` / `FindNeoschizomers` each document
  `ArgumentNullException`, and the sibling guards of `GetEnzyme`, `GetOverhangSequence` and
  `AreCompatible` are tested in the same file — but the one guard test there passes `"Nope"`, which
  exercises the `ArgumentException` arm. Neither arm was covered for these three.
- **Fix.** Two tests in `RestrictionAnalyzer_Filter_Tests.cs`:
  `IsoschizomerQueries_NullOrEmptyName_ThrowArgumentNullException` covers all **six** overloads
  (three methods × with and without `EnzymeSet`) for both `null` and `""`, the contract F104 states
  on `RequireEnzyme`; `IsoschizomerQueries_UnknownName_ThrowArgumentExceptionEagerly` covers the
  other arm for the same six, asserts it is eager (the throw happens at the call, not on first
  enumeration), and checks a resolvable name still works — including the hyphen spelling F103 added.

### F112 — RESTR-FIND-001/RESTR-FILTER-001 (H54, restriction rows): tests that could not fail

The `CrisprDesigner_GuideRNA_Tests.cs` rows of H54 belong to the orchestrator and are untouched here.

| Row | Was | Now |
|---|---|---|
| `RestrictionAnalyzer_FindSites_Tests.cs:115` — `Enzymes_ContainsAtLeast30CommonEnzymes` asserted `>= 30` where 39 is the documented fact, so nine enzymes could have vanished silently | `Is.GreaterThanOrEqualTo(30)` + 9 `ContainsKey` spot checks | renamed `Enzymes_ContainsExactlyThe39CommonEnzymes`; `Has.Count.EqualTo(39)` **and** the exact 39-name membership list, so neither a loss nor an unannounced addition passes |
| `:141-151` — `GetEnzymesByCutLength_ReturnsEnzymesOfCorrectLength` asserted `Is.Not.Empty` plus a soundness predicate against an invariant that says "exactly" | 3 lengths, soundness only | renamed `…_ReturnsExactlyTheEnzymesOfThatLength`; the exact set per length for 4 (9 enzymes), 6 (24), 8 (5), 13 (SfiI) and the two empty lengths 5 and 7, plus a new `GetEnzymesByCutLength_PerLengthSets_PartitionTheTable` asserting the length histogram `{4:9, 6:24, 8:5, 13:1}` and that the per-length sets sum to 39 |
| `:158-186` — `GetBluntCutters_ReturnsOnlyBluntEndEnzymes` / `GetStickyCutters_ReturnsOnlyStickyEndEnzymes` passed for an implementation returning only `AluI` | `Is.Not.Empty` + `All(...)` | counts pinned (10 blunt, 29 sticky), plus the blunt cut-offset identity and an explicit partition assertion against `Enzymes.Keys` |
| `RestrictionAnalyzer_Filter_Tests.cs:143-158` — `GetEnzymesByCutLength_Range6To6_ReturnsOnly6Cutters` | soundness + a 5-name spot check | the exact 24-name set, plus `ByLength(6) == ByLength(6,6)` (INV-04) |
| `RestrictionAnalyzer_Filter_Tests.cs:201-215` — `GetEnzymesByCutLength_Range4To4_ReturnsOnly4Cutters` | soundness + a 4-name spot check | the exact 9-name set, plus `ByLength(4) == ByLength(4,4)` |

Every number above is computed from the shipped table, not transcribed: 39 built-ins, lengths
`{4: 9, 6: 24, 8: 5, 13: 1}`, 10 blunt and 29 sticky. The blunt count cross-checks against
Biopython 1.88, whose `EcoRV.compatible_end()` returns the 137 blunt REBASE cutters (F109), of which
these 10 are the common-set share.

### F113 — RESTR-FIND-001/RESTR-DIGEST-001 (H27): the documented IUPAC-tolerant scan was unreachable over MCP

- **Defect, reproduced.** All seven sequence-construction sites in the restriction MCP wrappers built
  `new DnaSequence(sequence)`, and `DnaSequence` accepts only A/C/G/T
  (`Core/DnaSequence.cs:19-31` → `ValidateSequence`). So scanning, digesting or mapping an
  N-containing contig or a soft-masked assembly was an **error over MCP, not a result** — while
  `RestrictionAnalyzer` documents the raw-string `FindSites` overloads as handling exactly that
  (`FindSites("GGCCNNNNNGGCC", "SfiI")` → the site at 0) and `FindPatternOccurrences` as "exposed for
  callers that hold … a sequence that `DnaSequence` would reject". Probe:
  `new DnaSequence("nnnGAATTCnnnGGATCCNNN")` throws `ArgumentException`, and every one of
  `find_restriction_sites`, `find_all_restriction_sites`, `restriction_digest`, `digest_summary`,
  `partial_digest`, `restriction_map` failed on it.
- **Fix.** `DnaSequence` is B01's file and is not touched. Instead:
  1. The **five private cores** that still took a `DnaSequence` — `CollectForwardCuts`, `DigestCore`,
     `DigestCircular`, `CreateMapCore`, `PartialDigestCore` — now take the `string` they already
     reduced it to (each used only `.Sequence` and `.Length`). Every `DnaSequence` entry point passes
     `sequence.Sequence`, so **no behaviour changes**: see the regression evidence below.
  2. A new partial file `RestrictionAnalyzer.RawSequence.cs` adds the raw-string public entry points,
     each forwarding to the same core — `FindAllSites`, `FindSites` (names and enzyme records),
     `Digest` (both), `GetDigestSummary` (both), `PartialDigest` (both), `CreateMap` (both). A null or
     empty sequence yields the natural empty result, as on the pre-existing
     `FindSites(string, string)`; the enzyme arguments are still validated eagerly (F104); the
     ASCII-only upper-case fold of F49 is preserved.
  3. The six MCP wrappers route through them, and each `sequence` parameter description now states
     that IUPAC codes are accepted and what the matching rule is.
- **Post-fix, over MCP, on `nnnGAATTCnnnGGATCCNNN`** (EcoRI at 0-based 3, BamHI at 12):
  `find_restriction_sites` → 4 records (2 windows × 2 strands, forward positions 3 and 12);
  `find_all_restriction_sites` → 10; `restriction_digest` → fragments 4, 17 (G^AATTC, 0-based bond 4);
  `digest_summary` → 2 fragments; `partial_digest` → 6 ((k+1)(k+2)/2 for k=2);
  `restriction_map` → `totalSites` 2. The masked runs create **no** extra site: a sequence ambiguity
  code blocks a concrete site code rather than widening it, so `AluI` still finds nothing and
  `find_pattern_occurrences("GAANTC", "GAATTC")` is empty — Biopython's rule, which
  `find_pattern_occurrences("GGCCNNNNNGGCC", "GGCCNNNNNGGCC") == [0]` shows from the other side.
- **Regression evidence that the core refactor is a no-op.** The three oracle differentials were
  re-run against the refactored build and are **byte-identical** to the pre-refactor runs:
  linear 88 972 rows vs Biopython (76 divergent, the same 5 H7 enzymes), circular 63 336 rows vs
  Biopython (116 extra, 0 missing), and 590 rows vs EMBOSS on the 5 H7 enzymes (0 divergent).
  A test also pins every raw-string overload against its `DnaSequence` sibling on a pure-ACGT
  molecule, both topologies, with and without methylation.
- **Tests.** `RestrictionAnalyzer_FindSites_Tests.cs`:
  `RawStringOverloads_OnPureAcgt_EqualTheDnaSequenceOverloads` (2 topologies × 11 members),
  `RawStringOverloads_AcceptAnIupacSequenceDnaSequenceRejects` (including the `DnaSequence` throw that
  is the defect), `RawStringOverloads_NullOrEmptySequence_YieldTheEmptyResult`.
  New `Seqeron.Mcp.MolTools.Tests/RestrictionIupacAndCustomEnzymeTests.cs`:
  `DnaSequence_StillRejectsTheMaskedMolecule`, `EveryRestrictionTool_AcceptsAnIupacSequence`,
  `IupacCodesInTheSequence_DoNotCreateExtraSites`.
- **Note.** `find_degenerate_motif`, which H27 mentions as having the same guard, is on the Analysis
  server (`MotifFinder`, batch B05) and is not this slice's file; it is recorded under
  *Cross-batch requests* below.

### F114 — RESTR-FIND-001/RESTR-FILTER-001 (H28): three C# capabilities had no MCP expression

| Capability | Before | After |
|---|---|---|
| ranged `GetEnzymesByCutLength(min, max, set)` | only the exact-length overload was exposed, so "all 4–6 bp cutters" took three calls | `enzymes_by_cut_length` gained an optional `max_length`; measured `(4)` → 9, `(6)` → 24, `(5)` → 0, `(4, max 6)` → **33** in one call, `(4, "rebase", 6)` → **597**, and an inverted range → empty (the doc's INV-05) |
| every `IEnumerable<RestrictionEnzyme>` overload | no MCP expression, so a caller-supplied enzyme, or one of the 334 REBASE enzymes whose cut REBASE leaves undefined (`GetEnzyme("Dpi3084I")` → `null`), could not be used at all | two new tools, `find_sites_custom_enzyme` and `restriction_digest_custom_enzyme`: recognition sequence + the two cut offsets, optional Type IIB second pair, optional methylation, IUPAC accepted |
| `FindPatternOccurrences` | no MCP expression (H27 lists it here) | new tool `find_pattern_occurrences` |

- **Post-fix numbers.** `find_sites_custom_enzyme("AAAGAATTCAAA", "GAATTC", 1, 5)` returns records
  **identical** to `find_restriction_sites(..., ["EcoRI"])` — `p3F cut 4`, `p3R cut 8`; the Type IIS
  geometry `("GGTCTC", 7, 11)` gives `p10F cut 17`; BaeI's two-sided geometry
  `("ACNNNNGTAYC", -10, -15, second (23, 18))` yields two forward records, one flagged as the second
  cut. `restriction_digest_custom_enzyme` on the EcoRI-by-hand definition gives 4, 8 linear and a
  single 12-nt fragment circular (one cut on a circle linearises it). Half a second cut pair, an
  inconsistent stagger and a non-IUPAC recognition sequence are all errors.
- **Bookkeeping.** `scripts/skills/gen-catalog.py`'s `EXPECTED_COUNTS["moltools"]` 56 → 59 and
  `EXPECTED_TOTAL` 483 → 486 (the script calls this table its own source of truth; without the bump
  `--check` warns). Catalog and the ten generated files regenerated with the script;
  `gen-catalog.py --check` and `sync-github-skills.py --check` both clean. Three new doc pairs under
  `docs/mcp/tools/moltools/`, whose `description` and `inputSchema` were generated from the live MCP
  schema rather than typed, and the seven existing restriction `.mcp.json` files re-synced the same
  way (`partial_digest`'s doc-pair contract test enforces this). `ToolRegistrationContractTests`
  gained the three names. The `bio-moldesign` skill's `reference/tool-map.md` gained the three rows
  and the IUPAC note, and the `.github/skills` mirror was resynced; the "483 tools" figure was
  updated in the 7 skill/strategy documents that state it.
- **Tests.** `EnzymesByCutLength_RangedForm_IsOneCall`,
  `FindSitesCustomEnzyme_EcoRIByHand_EqualsTheNamedTool`,
  `RestrictionDigestCustomEnzyme_MatchesTheNamedToolAndValidates`,
  `IupacCodesInTheSequence_DoNotCreateExtraSites` (the `find_pattern_occurrences` contract).

## Cross-batch requests

### B05 — MCP `find_degenerate_motif` rejects non-ACGT (verified; must be registered by the orchestrator)

**Verified against the current code, per the README's "Verify first" rule.** H27 names this tool as
having "the same guard, so there is no workaround".

- `AnalysisTools.cs:1313` calls `RequireDna(sequence, nameof(sequence))`, and `RequireDna`
  (`AnalysisTools.cs:15-22`) does `DnaSequence.TryCreate` and throws
  `ArgumentException("Invalid DNA sequence")` when it fails. `DnaSequence` accepts only A/C/G/T
  (`Core/DnaSequence.cs:19-31`).
- **Minimal repro:** `find_degenerate_motif("nnnGAATTCnnnGGATCCNNN", "GAATTC")` → `ArgumentException`,
  i.e. a degenerate-motif scan over an N-containing contig or a soft-masked assembly is an error over
  MCP, exactly as the restriction family was before F113.
- **Not already resolved there:** B05's report records no such item.
- **The fix is one line, but it is a behaviour change.** `MotifFinder` already has the raw-string
  overload `FindDegenerateMotif(string sequence, string motif, CancellationToken)`
  (`Analysis/MotifFinder.cs:118`), so the wrapper only has to call it instead of `RequireDna`. But the
  tool would go from *erroring* to *returning results*, and that overload's own doc states a
  different ambiguity rule from the restriction family's — "any non-ACGT sequence character
  (including IUPAC ambiguity symbols, `U` and gaps) matches no motif code" — so whether a sequence
  `N` should be matched by a motif `N` (as Biopython's restriction regex does) is a decision for
  B05's owner. By the README that makes it **large → register**, not "small → do it yourself".
- **This session cannot register it:** its brief forbids editing `QUEUE.md` (the parallel session is
  writing there). **Orchestrator: please add this row** to
  `QUEUE.md` → *Requests to finished batches*:

  | # | Target | File:line / test | Evidence | From | Status |
  |---|---|---|---|---|---|
  | R11 | B05 | `src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs:1313` (`find_degenerate_motif` → `RequireDna`) | `find_degenerate_motif("nnnGAATTCnnnGGATCCNNN", "GAATTC")` → `ArgumentException "Invalid DNA sequence"`; `DnaSequence` is ACGT-only (`Core/DnaSequence.cs:19-31`). `MotifFinder.FindDegenerateMotif(string, string, CancellationToken)` (`Analysis/MotifFinder.cs:118`) already exists, so the wrapper only needs to route through it — but the tool changes from erroring to answering, and that overload treats a sequence `N` as matching nothing, where B08's restriction scan lets a pattern `N` match it (Biopython's rule). B05 decides which semantics the tool should advertise | B08 restriction slice (H27) | open |

## Notes for the orchestrator (not defects of this slice)

- **One MCP test is red on the branch and it is a CRISPR doc, not this slice's.**
  `B07ToolDocsContractTests.DocPair_SchemaAndExamples_MatchTheTool("calculate_on_target_rule_set2")`
  fails because the tool's `[Description]` and its `.mcp.json` `description` disagree about the PAM
  clause: the attribute says "3 nt PAM (must be NGG) + 3 nt downstream, A/C/G/T only; the NGG
  requirement is …" and the json says "3 nt PAM + 3 nt downstream (the PAM is not audited, as in
  CRISPOR; …)". **Pre-existing on `origin/claude/stoic-maxwell-0olr2z` at `3bcd8e07`, before this
  slice's commit:** both strings are byte-identical there, and
  `git diff origin/claude/stoic-maxwell-0olr2z -- MolToolsTools.cs` from this slice's tree contains
  **0** lines touching that text. It looks like the H22 attribute change and the json's new wording
  landed out of step. This session's brief forbids touching any CRISPR doc, so it is reported rather
  than fixed — one field in
  `docs/mcp/tools/moltools/calculate_on_target_rule_set2.mcp.json`, or the attribute, needs to match
  the other. Every restriction row of that test, and all four `AuditedSourceStems` of
  `DocRef_LineAnchor_PointsAtMethodDeclaration`, are green.
- **`sync_moltools_docref.py` still reports 11 drifted anchors, all in B02's files**
  (`CodonOptimizer.cs` ×10, `CodonUsageAnalyzer.cs` ×1) and so outside the H29 guard's
  `AuditedSourceStems` = PrimerDesigner / ProbeDesigner / CrisprDesigner / RestrictionAnalyzer. They
  predate this slice (none of those files is touched here) and no test covers them; widening the
  guard is H29's call, not this slice's.

- `B08.md:162` (an earlier F-entry's test list) names
  `FindSites_Circular_MoleculeShorterThanSite_ReportsNothing`, which F100 replaced with
  `FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite`. Per the coordination rules this
  session does not edit `B08.md` beyond ticking its own items, so the rename is recorded here.
  Likewise H54 renamed `Enzymes_ContainsAtLeast30CommonEnzymes` →
  `Enzymes_ContainsExactlyThe39CommonEnzymes` and
  `GetEnzymesByCutLength_ReturnsEnzymesOfCorrectLength` →
  `GetEnzymesByCutLength_ReturnsExactlyTheEnzymesOfThatLength`; no `docs/Validation/reports/**`,
  `docs/Evidence/**` or `docs/algorithms/**` file referenced the old names (checked by grep).
- `H29` (not this slice's item) notes that `partial_digest`'s `docRef` cannot be matched by the
  guard's `(?<file>\w+\.cs)` regex because `\w` rejects the dot in
  `RestrictionAnalyzer.PartialDigest.cs`. All 11 restriction anchors are themselves correct; the
  regex widening belongs to H29.

## Heavy-tier inventory (listed, not run — the orchestrator runs the heavy tier once for B08)

Every `Fuzzing/`, `Properties/` and `Metamorphic/` test that exercises a method this slice changed,
with whether its expectation could be affected. **Conclusion: none needs updating.** The reasons are
structural, so they hold without running the tier:

- the five enzymes F101 changes (`HauII`, `SgeI`, `MspJI`, `LpnPI`, `SgrTI`) are **REBASE-only** —
  none is one of the 39 built-ins — and every heavy-tier test below draws its enzymes from
  `RestrictionAnalyzer.Enzymes` or names a built-in literally;
- F100 changes only circular molecules **shorter than the recognition sequence**, and every
  heavy-tier molecule is far longer;
- F113's core refactor is proven byte-identical on three differentials (see F113), and its new public
  members are additive;
- F103/F104/F105 widen or harden argument handling; the only tightening is "an empty enzyme-name
  array now throws", which the orchestrator's own `fdb10a0a` already wrote the fuzz test for.

| Heavy-tier test | Touches | Affected? |
|---|---|---|
| `Fuzzing/MolToolsFuzzTests.FindSites_EmptyEnzymeList_ThrowsArgumentException_NotAnEmptyScan` | F104 | **No — it is the new contract.** Asserts `ThrowExactly<ArgumentException>` for `Array.Empty<string>()` and `ArgumentNullException` for `(string[])null!`; `RequireEnzymeNames` throws exactly those. Written by the orchestrator's `fdb10a0a` for H24 |
| `…FindSites_UnknownEnzymeName_TypedOverload_ThrowsArgumentException` | F104 | No. `FindSites(sequence, "NotARealEnzyme")` binds the single-name `FindSites(DnaSequence, string)` overload, which was already eager; F104 only made the `params` form match it |
| `…FindSites_UnknownEnzymeName_RawOverload_ThrowsArgumentExceptionOnEnumeration` | F104 | No. Its lambda calls `.ToList()`, so it passes whether the throw is eager or lazy. Its XML comment still says "only on ENUMERATION"; for the raw-string overload that is still literally true (the name is resolved before the sequence is read, but the enumerable is returned lazily) |
| `…FindSites_NullOrEmptyEnzymeName_TypedOverload_ThrowsArgumentNullException` | F104 | No. `RequireEnzyme` keeps `ArgumentNullException` for null **and** empty; F104 documented that rather than changing it |
| `…FindSites_NullDnaSequence_ThrowsArgumentNullException` | F113 | No. The new raw-string sibling made a bare `null!` ambiguous only where the test wrote one; the two call sites in `RestrictionAnalyzer_MethylationEntryPoints_Tests` were given explicit `(DnaSequence)null!` casts and this fuzz file compiles unchanged |
| `…FindSites_NonDnaSequence_TypedOverload_RejectedAtConstruction` | F113 | No — and it is the defect's own witness: `DnaSequence` still rejects non-ACGT, which is why the raw-string entry points exist |
| `…FindSites_NonDnaSequence_RawOverload_NoSitesNoCrash` | F113 | No. The raw-string rule is unchanged; the new overloads apply the same `ToUpperAscii` fold and the same matcher |
| `…FindSites_DegenerateIupacEnzyme_MatchesViaAmbiguityCodes`, `…FindSites_JunkInsideRecognitionWindow_SuppressesThatSite`, `…FindSites_KnownEcoRISite_DetectedAtCorrectPosition`, `…FindSites_EnzymeNameLookup_IsCaseInsensitive`, `…FindSites_EmptySequence_IsEmptyAndDoesNotThrow` | F100, F103, F113 | No. Linear, built-in enzymes, molecules longer than the site; F103 only *widens* name resolution (a trimmed or hyphenated spelling now also resolves), it rejects nothing it used to accept |
| `…Digest_KnownEcoRISites_YieldCorrectFragmentCountAndConserveMass`, `…Digest_LinearNoCutSites_YieldsSingleWholeFragment`, `…Digest_CircularWithTwoSites_YieldsTwoFragments_MassConserved`, `…Digest_CircularZeroCutSites_YieldsSingleCircularFragment_NoWrapAroundCrash` | F100, F101, F110 | No. EcoRI only; the circular rows are the `k=2 → 2` and `k=0 → 1` behaviour F110 *documented* rather than changed |
| `…Digest_EmptySequence_YieldsSingleEmptyFragment` | F113 | No — and it pins the behaviour the new raw-string overloads reproduce for a null or empty string (one zero-length fragment) |
| `…Digest_OneHundredPlusEnzymes_CompletesAndConservesMass` | F101 | No. Its catalog is `RestrictionAnalyzer.Enzymes.Values` (the 39 built-ins, all palindromic with symmetric cuts), so the F101 rule collapses every reverse record exactly as before |
| `Properties/RestrictionProperties` — all 18 tests (`FindSites_*`, `Digest_*`, `DigestSummary_*`, `FilteredSites_*`) | F100, F101, F104, F113 | No. All linear, all built-in enzymes (`EcoRI`/`HindIII`/the table), all molecules ≥ 16 nt. `FindSites_Positions_WithinBounds` asserts `Position < Length`, which is the **linear** bound F110 scoped INV-01 to; no property test uses `MoleculeTopology.Circular` |
| `Fuzzing/RestrictionFilterFuzzTests` — all 8 tests | F112 (tests only), F110 (doc only) | No. They assert against values computed from the shipped table (`GetBluntCutters`, `GetEnzymesByCutLength`), not against the doc prose F110 corrected, and this slice changed neither method |
| `Properties/CodonProperties.RemoveRestrictionSites_PreservesProtein…`, `Metamorphic/CodonMetamorphicTests` (CAI/usage row), `Metamorphic/B02FinisherMetamorphicTests` (3 rows) | F105 | No. They go through `CodonOptimizer.RemoveRestrictionSites`, whose only contact with this slice is `FindPatternOccurrences(ToMatchableDna(sequence), pattern, from)` — F105 only added tolerance for a **null** sequence, and that caller never passes null |
| `Metamorphic/RepeatsMetamorphicTests.FindPalindromes_KnownRestrictionSite…` | — | No. It names a restriction site as test data for `RepeatFinder`; it does not call `RestrictionAnalyzer` |
| `Metamorphic/MolToolsMetamorphicTests` — all 40 tests | — | No. PAM/guide/Tm/primer/probe only; none calls `RestrictionAnalyzer` |

**New behaviour that deserves heavy-tier cover and has none yet** (a suggestion for the single heavy
run, not a leftover of this slice — every item is already locked by unit tests, and the brief assigns
the heavy tier to the orchestrator): a property that the raw-string overloads agree with the
`DnaSequence` overloads on any generated pure-ACGT molecule and both topologies (the unit test
`RawStringOverloads_OnPureAcgt_EqualTheDnaSequenceOverloads` does this for one fixture), and a
property that a circular molecule of any length ≥ 1 yields exactly one site record per start index
that matches, which is F100's invariant.

## Test results

| Scope | Result |
|---|---|
| `dotnet build Seqeron.sln -c Debug` | 0 warnings, 0 errors |
| targeted tier after F100 (`~Restriction`, fast filter) | 621 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F100 | 471 passed / 0 failed |
| targeted tier after F101/F102 (`~Restriction\|~Digest\|~Enzyme\|~MolTools`, fast filter) | 3131 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F101/F102 | 471 passed / 0 failed |
| targeted tier after F103–F105 (same filter) | 3163 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F103–F105 | 472 passed / 0 failed |
| targeted tier after F106–F110 (`~Restriction\|~CodonOptimizer`, fast filter) | 910 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F106–F110 | 477 passed / 0 failed |
| `python3 scripts/skills/gen-catalog.py --check` after F106–F110 | fresh: 483 tools, 10 generated files up to date |
| targeted tier after F111/F112 (`~Restriction\|~CodonOptimizer`, fast filter) | 916 passed / 0 failed |
| targeted tier after F113/F114 (same filter) | 920 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F113/F114 | 483 passed / 0 failed |
| `gen-catalog.py --check` / `sync-github-skills.py --check` after F113/F114 | fresh: 486 tools, 10 generated files up to date / 21 skills byte-identical |

The counts rise between rows because the orchestrator is adding CRISPR tests to the same branch in
parallel; every row is a run of the branch as it stood at that commit, with 0 failures.

## Leftovers

**none.**

All seventeen items of this slice are implemented, each with the oracle number reproduced before the
fix and confirmed after it, a test that locks it, and a commit:

| Item | Entry | Oracle | Reproduced → after |
|---|---|---|---|
| H6 | F100 | Biopython 1.88 | 0 sites for all 10 probe rows → all 9 Biopython rows match; 63 336-row differential, 0 divergences at any length that can hold the site |
| H7 | F101 | EMBOSS 6.6.0 (adopted) vs Biopython 1.88 | 44 of 170 rows missing a cut → 0; 88 972-row Biopython regression, 76 rows, exactly the 5 named enzymes, all supersets |
| H45 | F102 | the shipped table | comment said "4"; computed 564 / 9 / 2 / 166 / 5 |
| H13 | F103 | Biopython 1.88 | `GetEnzyme("CviKI-1")` → null → resolves; cuts 2 and 6 match `CviKI_1.search` |
| H24 | F104 | the class's own contract | 4 meanings for one shape, lazy where the doc says eager → one meaning, eager |
| H25 | F105 | the method's own `<param>` doc + its sibling | `ArgumentNullException` → no occurrences |
| H30 | F106 | the shipped table | `CN`/`GN` order dependence; 12 revcomp-only pairs enumerated; `end2` added |
| H31 | F107 | measured | `p0F p8F p16F p16R p8R p0R`; Sau3AI is 2 records, not 1 |
| H32 | F108 | Biopython 1.88 (14 580 rows) | `TotalSites` 2 / no unique cutter where "forward records only" implies 1 |
| H33 | F109 | Biopython 1.88 source | 3 `compatible_end` implementations; 137 / 13 / 9 / 51 |
| H34 | F110 | Biopython 1.88 | circular `Position = 10` on a 13-mer with a 6-nt site |
| H35 | F110 | measured | circular k=0 → 1, k=1 → 1, k=2 → 2 |
| H36 | F110 | the shipped table | exactly one built-in (SfiI) has a non-ACGT window |
| H52 | F111 | the class's own contract | both arms untested → all six overloads covered |
| H54 (restriction rows) | F112 | the shipped table | `>= 30` → exactly 39; per-length sets; 10 blunt / 29 sticky |
| H27 | F113 | the library's own documented rule | `new DnaSequence("nnnGAATTC…")` throws and all 6 tools failed → all 6 answer; refactor byte-identical on 3 differentials |
| H28 | F114 | the C# surface | 3 capabilities with no MCP expression → `max_length` plus 3 new tools (486) |

Nothing is BLOCKED and nothing is deferred. The two items that are **not** this slice's to finish are
recorded above rather than here, because their target files belong to other sessions:

1. **B05 — MCP `find_degenerate_motif` rejects non-ACGT** (named by H27). Verified against the
   current code with a minimal repro; it is a behaviour change on another batch's public MCP tool, so
   by the README it is *register, not do*, and this session's brief forbids editing `QUEUE.md`. The
   exact register row is under *Cross-batch requests* above, ready to paste.
2. **The `calculate_on_target_rule_set2` doc pair is red on the branch tip** from before this slice
   (attribute and `.mcp.json` disagree on the PAM clause; verified byte-identical at `3bcd8e07` and
   untouched by this slice's diff). It is a CRISPR doc this session's brief forbids editing. Details
   under *Notes for the orchestrator*.

Two other observations, recorded because the completeness auditor greps for the words and they are
**not** leftovers of this slice:

- `RestrictionAnalyzer.PartialDigest.cs:17` and `:97` say the fragment **yield** is "out of scope".
  Pre-existing, not on the round-4 worklist, and a real modelling boundary rather than a
  simplification: the relative abundance of partial products depends on enzyme kinetics and digestion
  time, not on the sequence, so there is nothing sequence-derived to compute. The methods compute
  which fragments *can* arise, which is what they claim.
- `sync_moltools_docref.py` still reports 11 drifted anchors, all in B02's `CodonOptimizer.cs` /
  `CodonUsageAnalyzer.cs` and all outside the H29 guard's scope. Pre-existing; see
  *Notes for the orchestrator*.
