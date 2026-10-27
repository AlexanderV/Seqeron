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

## Test results

| Scope | Result |
|---|---|
| `dotnet build Seqeron.sln -c Debug` | 0 warnings, 0 errors |
| targeted tier after F100 (`~Restriction`, fast filter) | 621 passed / 0 failed |
| `Seqeron.Mcp.MolTools.Tests` after F100 | 471 passed / 0 failed |

## Leftovers

(pending — the slice is not finished)
