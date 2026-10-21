# Evidence Artifact: RESTR-DIGEST-001

**Test Unit ID:** RESTR-DIGEST-001
**Algorithm:** Restriction Digest Simulation (`Digest` both overloads, `DigestCore`, `DigestCircular`,
`GetDigestSummary` both overloads, `DigestFragment`, MCP `restriction_digest` / `digest_summary`;
`PartialDigest` + MCP `partial_digest` added 2026-10-15, §8)
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
  **end overhang** instead of as a second boundary. *(Refined by §9, F55: a reverse-strand record is
  skipped only when its window also matches on the forward strand; the reverse-strand sites of a
  non-palindromic enzyme are boundaries at their own top-strand cut.)*
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
fragments are identical to the reference. Which enzyme names are reported for such a boundary is a
choice the reference does not make at all; since review 2026-09 B08 F48 the fragment lists ALL of
them (`LeftEnzymes`/`RightEnzymes`, ordinal order of the table names) and the single-name fields are
the first of that list — previously the last enzyme seen, i.e. argument-order dependent. Built-in
collisions confirmed with Biopython 1.88 `search` (same 1-based position): MspI/HpaII `[5]` on
`AAACCGGAAA`; MboI/Sau3AI/BamHI `[5]` on `AAAGGATCCAAA`; ClaI/TaqI `[6]` on `AAAATCGATAAA`;
StuI/HaeIII `[7]` on `AAAAGGCCTAAA`; ScaI/RsaI `[7]` on `AAAAGTACTAAA`; FseI `[11]` / HaeIII `[7, 11]` on
`AAAAGGCCGGCCAAAA` (BglII/MboI·Sau3AI by the same nesting). All share `ovhgseq`/`ovhg` except
FseI (`ovhg = 4`, `CCGG`, 3') vs HaeIII (`ovhg = 0`, blunt); for that pair the reported end is the
representative's (FseI's) alone — the combined digest's physical end is not modelled (it also
depends on HaeIII's second cut, which falls on FseI's bottom-strand cut).

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

---

## 8. Partial digestion (`PartialDigest`, B08 D8 / F53, 2026-10-15) — CONSTRUCTED ORACLE

**No reference tool implements partial digestion** (Biopython 1.88 `Bio.Restriction` has only
complete `catalyse`; no other tool was available). The checks below are therefore a *constructed
oracle*: an explicit, independent derivation from the complete digest, plus the reference numbers
that do exist (Biopython cut positions and per-subset digests, and a textbook sample).

**Definition used (opened 2026-10-15).** Rosalind problem PDPL "Creating a Restriction Map"
(<https://rosalind.info/problems/pdpl/>, figure "Taken from Jones & Pevzner, *An Introduction to
Bioinformatics Algorithms*"), quoted: *"For a set X containing numbers, the difference multiset of X
is the multiset ΔX defined as the collection of all positive differences between elements of X. …
If X contains n elements, then ΔX will contain one element for each pair of elements from X, so that
ΔX contains (n choose 2) elements"*, and *"we can easily obtain a multiset L corresponding to the
distances between restriction sites on a chromosome"*; sample dataset `2 2 3 3 4 5 6 7 8 10`, sample
output `0 2 4 7 10`. X includes both molecule ends (0 and 10), so for k cut positions n = k+2 and the
partial digest has C(k+2, 2) = (k+1)(k+2)/2 fragments. English Wikipedia has no article on partial
digestion (checked: "Restriction digest", "Restriction map", "Restriction enzyme", "Restriction
fragment" contain no definition; "Partial digest problem", "Partial digestion", "Turnpike problem"
do not exist; full-text search for "partial digest" restriction returns only a citation title), so
Rosalind/Jones & Pevzner is the cited definition.

**Circular extension (this library's convention, consistent with F23–F27).** Every arc between two
cuts, cyclically (k start cuts × 0…k−1 missed cuts; k−1 = the one-cut linearised molecule of the
complete circular digest), plus the intact circle when all k sites are missed, reported exactly as
`Digest` reports an uncut circle — k²+1 fragments. Ends, overhangs and `LeftEnzymes`/`RightEnzymes`
(F25, F48) are taken from the complete digest's cut records.

**Reference numbers (Biopython 1.88, run 2026-10-15 in a fresh venv):**

| Input | Biopython | Seqeron `PartialDigest` |
|---|---|---|
| `RestrictionBatch([AluI, Sau3AI, TaqI]).search(Seq("AGCTGATCGA"))` | AluI [3], Sau3AI [5], TaqI [8] (1-based) → 0-based cuts {2,4,7}, i.e. Rosalind's X = {0,2,4,7,10} | lengths `{2,2,3,3,4,5,6,7,8,10}` = Rosalind's sample ΔX |
| linear, union over all 8 enzyme subsets of Biopython's fragment rule (single-enzyme subsets checked with `catalyse` itself) | `AG, CT, CGA, GAT, AGCT, CTGAT, GATCGA, AGCTGAT, CTGATCGA, AGCTGATCGA` | identical (10) |
| circular (`linear=False`), same union, empty subset = the uncut circle | `CT, GAT, CGAAG, CTGAT, CGAAGCT, GATCGAAG, CTGATCGAAG, GATCGAAGCT, CGAAGCTGAT, AGCTGATCGA` | identical (10 = 3²+1); single-enzyme `catalyse(linear=False)` gives the three linearised molecules |

Each enzyme cuts once in this molecule, so enzyme subsets *are* site subsets — this is the one place
a reference implementation can produce partial-digest fragments directly. Locked by
`PartialDigest_Linear_LengthsAreTheRosalindDifferenceMultiset_AndMatchBiopythonSubsetUnion` and
`PartialDigest_Circular_MatchesBiopythonSubsetUnion_AndEndsWithTheIntactCircle`.

**Constructed oracle on random sequences** (`PartialDigest_EqualsUnionOfCompleteDigestsOverAllCutSubsets`,
seed 531, 250 cases per topology, 4–89 bp, 1–6 enzymes from a 13-enzyme panel incl. FseI/HaeIII and
isoschizomer pairs, k ≤ 10, max k reached ≥ 6 asserted): for every one of the 2^k subsets of the
complete digest's cuts the subset digest is built from the complete digest's pieces — consecutive
pieces between retained cuts concatenated, left end from the first piece, right end from the last;
the empty circular subset is `Digest`'s zero-cut circle. The union (keyed on every member except
`FragmentNumber`, plus the missed-cut count = pieces − 1) must equal `PartialDigest`'s output as a set,
with no duplicate, and the count must be (k+1)(k+2)/2 / k²+1. **Result: 500/500 cases agree, 0
mismatches.** The oracle never uses PartialDigest's interval/wrap arithmetic — only the complete
digest and concatenation.

**Property** (`PartialDigest_ZeroMissedCutSubset_IsExactlyTheCompleteDigest`, seed 20261015, 300 cases
per topology): the `MissedCuts == 0` fragments equal `Digest(...)` member for member, numbering
included; `maxMissedCuts = 0` returns exactly `Digest(...)`; `maxMissedCuts = m` (1–3) is the unlimited
result's prefix with `MissedCuts ≤ m`. 0 failures.

**Defect caught by the tests before commit:** the first circular fragment-count formula
`k·min(m+1, k)` overflowed for the unbounded case (`m = long.MaxValue`, `m+1` wraps negative) so the
100,000-fragment guard never fired on a circle; replaced by `k·(min(m, k−1)+1)`, and the guard test
asserts both topologies (125,250 linear / 250,001 circular refused, 1497 / 1500 with `maxMissedCuts = 2`).

**Out of scope:** the yield (relative abundance / band intensity) of each partial product — it
depends on enzyme units, time and per-site cleavage rates, and no source opened gives a model; only
which fragments *can* arise is computed.

**Sources opened for §8:** `https://rosalind.info/problems/pdpl/` (2026-10-15); en.wikipedia.org pages
listed above (2026-10-15, reachable this session); Biopython 1.88 (`pip install biopython==1.88` into a
scratchpad venv) `RestrictionBatch.search`, `catalyse`.

---

## 9. Distal (single-sided Type IIS) cutters and non-palindromic reverse-strand sites (B08 WP11a / D6 part 1, F55, 2026-10-15)

**Gaps found.** (1) `Digest`/`PartialDigest` accepted enzyme **names** only, so no enzyme outside
the 39-entry table — no Type IIS cutter at all — could be digested. (2) The cut collector kept
`IsForwardStrand` records only (the convention stated in "Derived conventions" above). That is
right for a palindromic site, whose reverse record is the same site, but a **non-palindromic**
enzyme's reverse-strand occurrence is a different site; Biopython `NonPalindromic._search` adds it
via `_rev_modify` (`location − fst3`). Every reverse-strand site of a non-palindromic enzyme was
therefore never cut. (3) The end-overhang window was positioned from `CutPosition`, which for a
reverse record is the bottom-strand cut.

**Fix.** New overloads `Digest(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>)` and
`PartialDigest(DnaSequence, MoleculeTopology, int?, IEnumerable<RestrictionEnzyme>)` (eager
validation: null / empty / null element / non-IUPAC pattern). The collector keys each site by its
top-strand cut (forward: `CutPosition`; reverse: `CutPosition + CutPositionForward −
CutPositionReverse`, mod `n` on a circle) and skips a reverse record only when its window also reads
as the site on the forward strand — Biopython's "one cut per window, forward branch first"
(`Palindromic._search` searches one strand; `NonPalindromic._search`'s regex alternation takes the
forward group). The overhang window is read relative to the top-strand cut, which is valid for both
orientations because `crick − watson = CutPositionReverse − CutPositionForward` either way. The
name-based entry points are unchanged in behaviour (all 39 built-ins are palindromic, so every
reverse record is skipped exactly as before); `Digest_EnzymeOverload_BuiltInEnzymeRecord_EqualsNameOverload`
locks that.

**Differential** (harness described in RESTR-FIND-001 Evidence §2.6; every OneCut enzyme of the
Biopython 1.88 snapshot as a custom record, both topologies; "before" = the new overloads over the
old collector):

| check (N = 60 run) | compared | before (mismatching cases) | after |
|---|---|---|---|
| fragment (start, length), distal | 53,224 fragments / 18,360 cases | **8,424** | 0 |
| end overhang type + bases, distal | 106,448 ends | **2,350** | 0 |
| partial-digest multiset, distal | 331,080 fragments | **2,850** | 0 |
| fragment (start, length), in-site non-palindromic control | 29,508 / 13,824 cases | **301** | 0 |
| end overhangs / partial, in-site control | 59,016 / 22,874 | **67** / **100** | 0 / 0 |

Plus 5 cases of the documented TscAI/TspRI divergence (RESTR-FIND-001 Evidence §2.6 / §3), excluded
from both columns. The overhang oracle uses `ovhg = fst5 − (size + fst3)`, which equals Biopython's
stored `ovhg` for 727/729 OneCut enzymes (not for TscAI/TspRI). Circular fragments are compared
against the sorted, reduced `search` list, not Biopython's circular `catalyse` (oracle caveat above;
574 circular `search` results in the N = 60 run contained a raw position outside `[1, n]`, all on
circles shorter than the enzyme's cut offset).

**Locked values** (`RestrictionAnalyzer_DistalCutter_Tests`, all Biopython 1.88 output): BsaI linear
`ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT` → `search` `[10, 29]`, `catalyse` 9/19/13, insert
`AATGCCCCCCCCCCCCCCG` with 5' ends `AATG`/`AATG`; same circular → (9, 19), (28, 22); BsaI site across
the origin of `TCAAATGCCCCCCCCATGGTC` → `[4]`, one 21-nt fragment; FokI circular reverse site at 0 →
`[30]`; BsmBI 7/15/11, BbsI 9/11/13, SapI 10/7/13 (3-nt `GCA`/`TTT`), MlyI 15/4/13 blunt, AcuI 22/10
(3' `AA`), FokI 16/8; partial digest of the BsaI molecule = intervals over {0, 9, 28, 41} (linear) and
{19, 22, 41, 41, 41} (circular).

**Not in this unit (→ WP11b — all closed there, see §10):** two-sided Type IIB cutters (BaeI/BcgI; 25 in the snapshot);
enzyme-level overhang / compatibility of distal ends (`GetOverhangSequence` → `null`,
`AreCompatible` → `false`); custom enzymes in `GetDigestSummary` / `CreateMap`; MCP exposure.

## 10. Two-sided (Type IIB) cutters and custom enzymes in summary / map (B08 WP11b / D6 part 2, F56, F58, 2026-10-15)

- `RestrictionEnzyme.SecondCut` (F56): each site of a two-cut enzyme contributes two top-strand
  boundaries, so the site is excised as its own fragment; every end's overhang is read with that cut's
  own pair. Differential vs Biopython 1.88 over all 25 `TwoCuts` enzymes (RESTR-FIND-001 Evidence §2.7):
  fragments 1 988 → **0** mismatching cases, partial-digest multisets 680 → **0**, 17 718 end overhangs
  0 mismatches. Biopython's `catalyse` on palindromic two-cutters with ≥ 2 sites mis-slices its unsorted
  cut list (empty fragment, lengths not summing to the molecule) — 111 cases, not reproduced.
- Locked values: BaeI `A×20 ACAAAAGTACC T×20` → 10/33/8 with 3' ends `AAAAA` / `TTTTT`; its reverse
  complement → 13/33/5; site at 0 → only the downstream cut (`[24]`, 23/8); circular origin-spanning
  site → cuts 17/25; UcoMSI (palindromic) → 13/18/25; partial-digest lengths {8, 10, 33, 41, 43, 51}.
- F58: `GetDigestSummary(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>)` and
  `CreateMap(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>)` (and multi-enzyme
  `FindSites`) take custom records; `CreateMap` counts physical sites (forward windows + reverse windows
  that do not also read forward), identical to the old forward-record count on the palindromic
  built-ins (locked by `EnzymeCollectionOverloads_BuiltInRecords_EqualTheNameOverloads`).


---

## 11. Methylation-blocked digest (B08 round 2 E2 / F68, 2026-10-21)

`Digest`, `GetDigestSummary` and `CreateMap` accept an optional methylase set (name form `(sequence, topology,
methylases, params enzymeNames)`, record form `(sequence, topology, enzymes, methylases)`); the mask is computed once per call
and filters `FindSitesCore`, so boundaries, overhangs and fragment order are those of the unblocked sites. The oracle for the
site set is EMBOSS `restrict -methylation` (RESTR-FIND-001 Evidence §7); EMBOSS has no digest, so fragment values are derived:
on `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` (51 bp) EMBOSS keeps EcoRI 0 and 18 and blocks MboI 10 and 36 with Dam, so
`Digest(EcoRI, MboI)` is 1/18/32 (cuts 1 and 19; unblocked: 5 fragments, cuts 1, 10, 19, 36) and the circular digest is 18/33. `PartialDigest` is
carried since B08 round 3 G3 (§12). Locked in `RestrictionAnalyzer_Methylation_Tests` (`Digest_BlockedSitesDoNotCut`, `Digest_Circular_AndEnzymeListOverload`,
`GetDigestSummary_AndCreateMap_CarryTheSet`) and the MCP `restriction_digest` / `digest_summary` / `restriction_map` tests.

## 12. Methylation-blocked partial digest (B08 round 3 G3 / F83, 2026-10-21)

`PartialDigest(sequence, topology, maxMissedCuts, methylases, params enzymeNames)` and the record form `(sequence, topology,
maxMissedCuts, enzymes, methylases)` thread the same mask `DigestWithMask` uses into `CollectForwardCuts` (`PartialDigestCore` gained an
optional `bool[]? blocked`); MCP `partial_digest` gained `methylation` (`none|dam|dcm|dam+dcm|emboss`, `MethylaseSets.FromName`). F68 listed
`PartialDigest` as "not carried" with no proof; there is none: the partial digest is a function of the cut list alone (§8), so the
methylation-aware cut list is the whole change — a blocked site is neither a cut nor a missed cut, and the `maxMissedCuts` / 100 000 guard
counts surviving cuts only.

**Oracle (independent of the library).** The surviving cut set is built from two executed tools, not from `FindSites`: Biopython 1.88
`Bio.Restriction.<Enzyme>.search(Seq, linear=...)` gives every cut of every palindromic enzyme, and a cut is kept when EMBOSS 6.6
`restrict` (table generated from Biopython, one sequence per call, with `-methylation` and the chosen `-mfile`) still reports that
enzyme's site at `cut - fst5`. The partial-digest fragments are then the F53 constructed set over those cuts: linear = every `[X_i, X_j)`
over `X = {0, cuts, L}` with `j - i - 1` missed cuts; circular = for each cut, the arcs to the next 1..k cuts (`k-1 .. 0` missed) plus the
intact circle. Surviving cuts (0-based bases 5' of the cut): L2 (51 bp, `EcoRI, MboI, SmaI, TaqI`) none `{1,10,19,30,36,39,45}` (36 fragments), Dam /
Dam+Dcm / emboss `{1,19,30,39,45}` (21); L1 (55 bp, `MboI, Sau3AI, BamHI, BclI`) none `{5,15,24}`, Dam `{}` (the intact molecule, `MissedCuts 0`);
L3 (56 bp, `KpnI, Acc65I, EcoRI`) none `{8,12,22,36,40,48,52}`, emboss `{22,48,52}`; a seeded 500-bp random molecule (8 enzymes) none 16 cuts
(153 linear / 257 circular fragments), emboss 10 cuts (66 / 101), Dam only 12 (91 / 145), Dcm only 15 (136 / 226). Linear and circular cut lists
agree (no methylase site spans the origin there). Library = oracle member for member (start, length, missed cuts) in both topologies, by enzyme
name and by enzyme record; the `MissedCuts = 0` prefix equals `Digest(sequence, topology, methylases, names)` (invariant 8c holds under
methylation). Default off (`null`, empty, `none`) is byte-identical for the six `{topology} x {maxMissedCuts}` cases checked.
Locked in `RestrictionAnalyzer_MethylationEntryPoints_Tests` and the MCP `RestrictionMethylationToolsTests` (`PartialDigest_*`).
