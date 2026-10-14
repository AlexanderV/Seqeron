# Evidence Artifact: RESTR-DIGEST-001

**Test Unit ID:** RESTR-DIGEST-001
**Algorithm:** Restriction Digest Simulation (`Digest` both overloads, `DigestCore`, `DigestCircular`,
`GetDigestSummary` both overloads, `DigestFragment`, MCP `restriction_digest` / `digest_summary`)
**Date Collected:** 2026-10-14
**Reviewer note:** every number in this file was computed in the 2026-10-14 review session by running
the code and the reference implementation. Nothing is carried over from an earlier pass. Out of
scope and therefore not re-measured here: site finding (RESTR-FIND-001) and enzyme
filtering/compatibility (`FindCompatibleEnzymes` / `AreCompatible` / the enzyme-level overhang,
RESTR-FILTER-001 — re-validated separately on the same day; `GetOverhang` was folded into the shared
`OverhangWindow` geometry by that unit, so `CutOverhangWindow` is now the only molecule-side reader
of it).

---

## Reference Oracle

### Biopython 1.88 — `Bio.Restriction`

**Package:** `biopython` 1.88, installed locally
(`/usr/local/lib/python3.13/dist-packages/Bio/Restriction/Restriction.py`), run offline.
**Accessed / run:** 2026-10-14.

Reference code read this session, with line numbers in that file:

| Member | Lines | What it establishes |
|---|---|---|
| `RestrictionType.search` | 514-539 | `linear=False` wraps the molecule; returns 1-based **top-strand** cut positions |
| `AbstractCut._drop` | 542-568 | Linear: drop a cut unless `1 < cut_on_watson <= length` **and** `1 < cut_on_crick <= length`. Circular: fold out-of-range cuts back by ±`length` **in place, without re-sorting** |
| `Defined._modify` | 789-814 | Top-strand cut `= location + fst5` |
| `catalyse` / `catalyze` | 1131-1178 | Linear: fragments are `d[1:r[0]]`, `d[r[x]:r[x+1]]`, `d[r[-1]:]`. Circular: `d[r[-1]:] + d[1:r[0]]` **first**, then the rest. `if not r: return (d[1:],)` — an enzyme that does not cut returns the whole molecule as one fragment, circular included |
| `ovhgseq`, `is_5overhang`, `is_3overhang`, `is_blunt` | 480-540 | Per-enzyme overhang sequence and 5'/3'/blunt classification |

Derived conventions used throughout:

- Biopython's 1-based top-strand cut `p` == Seqeron forward-strand `CutPosition` `p-1`.
- `catalyse` partitions the **top strand** at `search()`'s top-strand cuts. Seqeron therefore keeps
  only `IsForwardStrand` site records as fragment boundaries; the reverse-strand record of the same
  palindromic site carries the bottom-strand cut of that same site and is reported as the fragment
  **end overhang** instead of as a second boundary.
- Fragment order on a circle: Biopython emits the origin-spanning fragment **first**, Seqeron emits
  it **last** (ascending `StartPosition`). The multisets are identical; the orders differ by a
  rotation of one. All comparisons below rotate Biopython's circular output by one before comparing
  element-wise, and this convention is now stated in the public XML docs, the algorithm doc and both
  MCP tool docs.

### Oracle caveat found this session — Biopython's own circular `catalyse` is unsafe

`_drop` (L542-568) folds a circular out-of-range cut back by ±`length` **in place and does not
re-sort**, while `catalyse` (L1146-1176) indexes `r[0]` / `r[-1]` as if the list were ascending.
Measured:

```
AluI.search(Seq("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"), linear=False)    == [11, 2]     # unsorted
AluI.catalyse(Seq("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"), linear=False)
    == ("CTGTTATAGCTAAGAGCGAACCGGAGTATAGCTGTTATAG", "")                           # 40 bases, 31 bp molecule
```

**Frequency measured:** of 4000 pseudo-random circular sequences with a planted origin-spanning site
(seed 20260914), **55 (1.4%)** return an unsorted `search(..., linear=False)` and in **exactly the
same 55** cases `catalyse(..., linear=False)` violates the fragment-sum invariant.

Consequence for this validation: for circular cases the oracle is
`sorted(search(..., linear=False))` expanded into fragments by Biopython's own documented rule, and
`catalyse` is used directly only as a cross-check of that expansion (agreeing on every case outside
the 55). Seqeron sorts its cut positions (`SortedSet<int>`), so its fragments partition the molecule
in all 1452 comparison cases. A regression test locks the `AluI` case above.

---

## 1. Enzyme-table confound check — 39/39 agreement

Before comparing digests, the 39 built-in enzyme definitions were re-derived from
`enzyme.elucidate()` (padding-corrected via `enzyme.site`) and compared against the literals in
`RestrictionAnalyzer.cs`.

**Result: 39 entries, 39 matched on `(recognition sequence, CutPositionForward, CutPositionReverse)`,
0 differences.** The enzyme table is therefore not a confound in any number below.

## 2. Overhang chemistry — 39/39 agreement with `ovhgseq`

Seqeron's rule: at a cut with top-strand cut `watson`, the bottom-strand cut is
`crick = watson + (CutPositionReverse − CutPositionForward)`; the single-stranded extension spans
the top-strand window `[min(watson, crick), max(watson, crick))`, of length
`|CutPositionReverse − CutPositionForward|`; `r > f` ⇒ 5' overhang, `r < f` ⇒ 3', `r == f` ⇒ blunt.

Taking that window over the **recognition sequence** and comparing with Biopython:
**39 enzymes, 39 agreed on both the overhang sequence (`ovhgseq`) and the 5'/3'/blunt class
(`is_5overhang()` / `is_3overhang()` / `is_blunt()`), 0 differences.** Spot values:

| Enzyme | fwd / rev | Seqeron window | `ovhgseq` | Type (Seqeron / Biopython) |
|---|---|---|---|---|
| EcoRI | 1 / 5 | `AATT` | `AATT` | 5' / 5' |
| BamHI | 1 / 5 | `GATC` | `GATC` | 5' / 5' |
| Sau3AI | 0 / 4 | `GATC` | `GATC` | 5' / 5' |
| HindIII | 1 / 5 | `AGCT` | `AGCT` | 5' / 5' |
| NdeI | 2 / 4 | `TA` | `TA` | 5' / 5' |
| PstI | 5 / 1 | `TGCA` | `TGCA` | 3' / 3' |
| KpnI | 5 / 1 | `GTAC` | `GTAC` | 3' / 3' |
| SacI | 5 / 1 | `AGCT` | `AGCT` | 3' / 3' |
| PacI | 5 / 3 | `AT` | `AT` | 3' / 3' |
| SfiI | 8 / 5 | `NNN` | `NNN` | 3' / 3' |
| EcoRV | 3 / 3 | `` (empty) | `` | blunt / blunt |
| AluI | 2 / 2 | `` | `` | blunt / blunt |

Seqeron reports this window over the **molecule**, not over the pattern, so for a degenerate enzyme
it gives the real bases where Biopython's `ovhgseq` gives `NNN` (SfiI) — strictly more informative,
same rule.

---

## 3. Digest comparison vs Biopython — before & after

**Harness** (scratch, not committed): a .NET console program referencing
`Seqeron.Genomics.MolTools`, reading `sequence \t topology \t enzymes` lines and emitting the
fragment list as JSON; a Python driver building the Biopython expectation per case. "Before" is
commit `55579583` (the state this unit inherited), built from a `git worktree` of that commit;
"after" is the state committed by this unit. **The same 1452 cases, fixed seed 20260914, were run
against both.**

Case mix: 32 designed cases (16 sequences × both topologies: origin-spanning sites, uncut molecules,
edge cutters Sau3AI/MboI, 3'-overhang cutters, adjacent/overlapping 4-cutter sites, double digests),
700 uniform-random ACGT sequences 6–150 bp (single and 2–3 enzyme), 700 with an enzyme site
deliberately split across the origin (single and multi-enzyme), 120 over the `GATC` alphabet with
dense 4-cutter sets. Comparison is element-wise on the **fragment sequence list in order**.

| Group | Cases | Mismatches BEFORE | Mismatches AFTER |
|---|---|---|---|
| linear, single-enzyme | 363 | 0 | **0** |
| linear, multi-enzyme | 349 | 0 | **0** |
| circular, single-enzyme | 405 | **164** | **0** |
| circular, multi-enzyme | 335 | **149** | **0** |
| **Total** | **1452** | **313** | **0** |

Every one of the 313 pre-fix mismatches is circular; the linear path was already exact, single- and
multi-enzyme. Representative pre-fix failures (got → expected):

| Case | Before | Biopython (= after) |
|---|---|---|
| `TTCAAAAAAAGAA` circ, EcoRI | `["TTCAAAAAAAGAA"]` (uncut) | `["AATTCAAAAAAAG"]` |
| `GATCGATCGATC` circ, Sau3AI | `["GATCGATCGATC"]` (uncut) | `["GATC","GATC","GATC"]` |
| `GATCGATCGATC` circ, MboI+DpnI | `["TC","GA","TCGA","TCGA"]` | `["GA","TC","GA","TC","GA","TC"]` |
| `TTCGGATCCAAAAAAAGAA` circ, EcoRI+BamHI | `["GATCCAAAAAAAGAATTCG"]` | `["GATCCAAAAAAAG","AATTCG"]` |
| `ATTCAAAAAAAAGA` circ, EcoRI | `["ATTCAAAAAAAAGA"]` (uncut) | `["AATTCAAAAAAAAG"]` |

### 3b. Fragment-end overhang comparison (post-fix; the fields did not exist before)

A second suite (seed 777001, 1300 cases: 800 wrap-planted circular/linear pairs over all 39 enzymes,
300 uniform-random 10–150 bp with 1–3 enzymes, 200 dense `GATC`-alphabet cases with 1–3 enzymes
including the 3'-overhang cutters) compares, per fragment and in order, the **fragment sequence plus
all four overhang fields** (`LeftOverhangType`, `LeftOverhang`, `RightOverhangType`,
`RightOverhang`) against an independent Python derivation from Biopython's cut positions and enzyme
geometry.

**1300 cases, 0 mismatches, 0 fragment-sum violations.**

---

## 4. Hand-traced worked values (locked by tests)

| Call | Result | Reference |
|---|---|---|
| `Digest("AAAGAATTCAAA", "EcoRI")` | `AAAG` (0,4) + `AATTCAAA` (4,8); cut ends 5' `AATT` | `EcoRI.catalyse(Seq("AAAGAATTCAAA")) == ("AAAG","AATTCAAA")`; `EcoRI.ovhgseq == "AATT"` |
| `Digest("CTGCAGAAACTGCAG", "PstI")` | `CTGCA` + `GAAACTGCA` + `G`; cut ends 3' `TGCA` | `PstI.catalyse(...) == ("CTGCA","GAAACTGCA","G")`; `PstI.ovhgseq == "TGCA"` |
| `Digest("TTCAAAAAAAGAA", Circular, "EcoRI")` | 1 fragment `AATTCAAAAAAAG`, start 11, len 13, EcoRI/EcoRI, both ends 5' `AATT` | `search(...,linear=False)==[12]`; `catalyse(...,linear=False)==("AATTCAAAAAAAG",)` |
| `Digest("TTCAAAAAAAGAA", Linear, "EcoRI")` | 1 uncut fragment, both enzymes null | `EcoRI.search(Seq("TTCAAAAAAAGAA")) == []` |
| `Digest("GATCGATCGATC", Linear, "Sau3AI")` | `GATC` + `GATCGATC` | `Sau3AI.catalyse(Seq("GATCGATCGATC")) == ("GATC","GATCGATC")` (`search == [5]`) |
| `Digest("GATCGATCGATC", Circular, "Sau3AI")` | `GATC`,`GATC`,`GATC` at 0/4/8 | `search(...,linear=False)==[1,5,9]`; `catalyse == ("GATC","GATC","GATC")` |
| `Digest("GATCGATCGATC", Circular, "MboI","DpnI")` | 6 × 2 bp: `GA`,`TC`,`GA`,`TC`,`GA`,`TC` at 0..10 | batch `search(...,linear=False) == {MboI:[1,5,9], DpnI:[3,7,11]}` |
| `Digest("TTCGGATCCAAAAAAAGAA", Circular, "EcoRI","BamHI")` | `GATCCAAAAAAAG` (4,13) + `AATTCG` (17,6) | batch `search(...,linear=False) == {EcoRI:[18], BamHI:[5]}` |
| `Digest("GCTGTTATAGCTAAGAGCGAACCGGAGTATA", Circular, "AluI")` | `CTGTTATAG` (1,9) + `CTAAGAGCGAACCGGAGTATAG` (10,22), sum 31 | `sorted(search(...,linear=False)) == [2,11]`; Biopython's own `catalyse` is wrong here (see oracle caveat) |
| `Digest("CCGGAAACCGG", Linear, "MspI","HpaII")` | `C`,`CGGAAAC`,`CGG` at 0/1/8 — the shared cut is one boundary, not two; 5' `CG` ends | batch `search == {MspI:[2,9], HpaII:[2,9]}`; `MspI.ovhgseq == "CG"` |
| `Digest("GATCAAAGATC", Circular, "MboI","Sau3AI","DpnI")` | `GA`,`TCAAA`,`GA`,`TC` at 0/2/7/9; `GATC` 5' at the `^GATC` cuts, blunt at the `GA^TC` cuts | batch `search(...,linear=False) == {MboI:[1,8], Sau3AI:[1,8], DpnI:[3,10]}` |
| `GetDigestSummary("GAATTCAAAGAATTCAAA", Linear, "EcoRI")` | 3 fragments, sizes `[9,8,1]` | cuts {1,10} from `search == [2,11]` |
| `GetDigestSummary("GAATTCAAAGAATTCAAA", Circular, "EcoRI")` | 2 fragments, sizes `[9,9]` | same cuts, circular rule k → k |
| `GetDigestSummary("TTCGGATCCAAAAAAAGAA", Circular, "EcoRI","BamHI")` | 2 fragments, sizes `[13,6]`, sum 19 | batch search above |

Same-coordinate enzymes (`MspI`/`HpaII`, and `MboI`/`Sau3AI` vs `DpnI`) are checked explicitly in
the last two rows: the duplicate coordinate is merged into one boundary (`SortedSet<int>`) and the
fragments are identical to the reference; only the reported representative enzyme name — the last
enzyme seen for that coordinate — is a choice the reference does not make at all.

The circular fragment-count rule itself (k cuts → k fragments) is sourced to Univ. of Illinois MolBio
"Restriction Digestion/Gel Electrophoresis" and Addgene Plasmids 101 (a single cutter linearizes a
plasmid), as already recorded for the pre-existing circular tests; the fragment-sum rule is Addgene's
Restriction Digest Protocol ("the sum of the individual fragments should equal the size of the
original fragment").

---

## 5. Zero-length-fragment / `FragmentNumber` reachability

`DigestCore` skipped a zero-length piece while numbering fragments by **boundary index**, so a skip
would have broken both "k cuts → k+1 fragments" and contiguous numbering. Checked this session: the
linear forward-strand cut positions `DigestCore` consumes all come from `FindSitesLinear`, which
admits a site only when `CutsAreInternal(watson, crick, n)` holds, i.e. the cut lies in
`[1, Length−1]`; the positions are also distinct (`SortedSet<int>`). Every interval
`[boundary_i, boundary_{i+1})` is therefore non-empty and a zero-length fragment is **unreachable
through the public API**. The fix keeps the guard (defensive) but numbers fragments with an emission
counter, so the numbering is contiguous by construction; a test asserts count = k+1 (linear) / k
(circular), contiguous `1..n` numbering, all lengths > 0, `Length == Sequence.Length`, ascending
`StartPosition` and the fragment-sum invariant over 7 dense cut patterns × both topologies.

---

## 6. Coverage of the public surface

| Entry point | Topology reachable before | Topology reachable after |
|---|---|---|
| `Digest(DnaSequence, params string[])` | Linear only (by design) | unchanged (Linear) |
| `Digest(DnaSequence, MoleculeTopology, params string[])` | both, but Circular used the linear site finder | both, correct |
| `GetDigestSummary(DnaSequence, params string[])` | Linear only | unchanged (Linear) |
| `GetDigestSummary(DnaSequence, MoleculeTopology, params string[])` | **did not exist** | added |
| MCP `restriction_digest` | Linear only | `topology` parameter, default `Linear` |
| MCP `digest_summary` | Linear only | `topology` parameter, default `Linear` |

`DigestSummary` gained a `Topology` field (trailing, defaulted to `Linear`, so source-compatible) and
`DigestFragment` gained `LeftOverhangType` / `LeftOverhang` / `RightOverhangType` / `RightOverhang`
(trailing, defaulted to `null`). MCP tool counts are unchanged (moltools 53, total 480:
`scripts/skills/gen-catalog.py --check` passes) because only parameters and descriptions changed.

---

## 7. Sources opened this session

1. `/usr/local/lib/python3.13/dist-packages/Bio/Restriction/Restriction.py` (Biopython 1.88) —
   L514-568 (`search`, `_drop`), L789-825 (`_modify`, `_rev_modify`), L1131-1178 (`catalyse`),
   plus `ovhgseq` / `is_5overhang` / `is_3overhang` / `is_blunt` via the live objects.
2. Biopython's REBASE snapshot `Bio/Restriction/Restriction_Dictionary.py` via the live enzyme
   objects (`site`, `elucidate()`, `ovhgseq`) — used for the 39/39 table and overhang checks.
3. `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs` (whole file).
4. `src/Seqeron/Mcp/Seqeron.Mcp.MolTools/Tools/MolToolsTools.cs` (restriction tool region),
   `src/Seqeron/Mcp/Seqeron.Mcp.MolTools/Models/Results.cs`.
5. `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs`,
   `tests/Seqeron/Seqeron.Mcp.MolTools.Tests/RestrictionDigestTests.cs`,
   `tests/Seqeron/Seqeron.Mcp.MolTools.Tests/DigestSummaryTests.cs`.
6. `tests/TestSpecs/RESTR-DIGEST-001.md`, `docs/algorithms/MolTools/Restriction_Digest_Simulation.md`,
   `docs/mcp/tools/moltools/restriction_digest.{md,mcp.json}`,
   `docs/mcp/tools/moltools/digest_summary.{md,mcp.json}`.
7. `docs/Validation/VALIDATION_PROTOCOL.md`.

**Not reachable from this environment**, verified by me this session: the agent proxy answers
**403 to CONNECT** for `en.wikipedia.org:443`, `rebase.neb.com:443` and `www.addgene.org:443`
(`curl` exit 56, "CONNECT tunnel failed, response 403" for all three), matching what
`RESTR-FIND-001-Evidence.md` already recorded. The Wikipedia, REBASE and Addgene citations inherited
from the pre-existing TestSpec could therefore not be re-fetched; they are carried forward as-is and
are **not** the basis of any number in this file, every one of which comes from Biopython 1.88
(including its own REBASE 404 snapshot) or from running the code. The sourced claims that rest on
them — the circular fragment-count rule and the fragment-sum rule — are independently confirmed here
by Biopython's `catalyse` behaviour, which produces exactly k fragments for k circular cuts and
exactly k+1 for k linear cuts in all 1452 comparison cases.

**Update (WP-SRC, 2026-10-14).** Those denials were re-read from the proxy's status endpoint, not
re-attempted, and still stand. The inherited Wikipedia/REBASE/Addgene citations are now
**explicitly annotated as unverifiable in this environment** in the TestSpecs that carry them, so a
reader cannot mistake them for checked; and REBASE's organism and methylation columns turn out to be
obtainable from an allowed host after all (REBASE `allenz.104`, vendored in `crisprVerse/crisprBase`).
See [`B08-Source-Provenance.md`](B08-Source-Provenance.md).
