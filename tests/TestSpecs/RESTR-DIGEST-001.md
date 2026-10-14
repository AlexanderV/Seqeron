# Test Specification: RESTR-DIGEST-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | RESTR-DIGEST-001 |
| **Area** | MolTools |
| **Title** | Restriction Digest Simulation |
| **Canonical Class** | `RestrictionAnalyzer` |
| **Canonical Methods** | `Digest` (both overloads), `GetDigestSummary` (both overloads), `CreateMap`, `AreCompatible`, `FindCompatibleEnzymes` |
| **Complexity** | O(n + k log k) where k = number of cut sites |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-10-14 |
| **Re-validated** | 2026-10-14 (RESTR-DIGEST-001 unit, batch B08) — see `docs/Evidence/RESTR-DIGEST-001-Evidence.md` |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Restriction digest | Academic | Fragment generation, gel electrophoresis verification, fragment sum = original length |
| Wikipedia: Restriction enzyme | Academic | Type II enzymes, palindromic recognition, overhang types, enzyme Examples table with cut positions |
| Wikipedia: Restriction map | Academic | Site positions, unique cutters identification |
| Wikipedia: Sticky and blunt ends | Academic | Overhang compatibility rules: same type + complementary sequence required; blunt ends always compatible |
| Addgene: Restriction Digest Protocol | Protocol | Fragment sum = original length, validation by gel electrophoresis |
| Roberts RJ (1976) | Research | Restriction endonuclease classification |
| REBASE | Database | Enzyme data, overhang sequences |
| Biopython 1.88 `Bio.Restriction` (`Restriction.py`) | Reference implementation | `search`/`_drop` L514-568, `catalyse` L1131-1178, `ovhgseq`/`is_5overhang`/`is_3overhang`/`is_blunt`: top-strand fragment boundaries, circular (`linear=False`) search, uncut-circle behaviour, per-enzyme overhang sequence and class |
| Univ. of Illinois MolBio, "Restriction Digestion/Gel Electrophoresis" | Course material | Circular molecule with k cuts → k fragments |
| Addgene, Plasmids 101 | Protocol | A single cutter linearizes a plasmid |

---

## Invariants

### Digest Invariants

1. **Fragment Count (linear)**: With k unique top-strand cut positions, produces k+1 fragments (Source: Wikipedia Restriction digest; Biopython `catalyse` emits `d[1:r[0]]`, the k−1 inner slices and `d[r[-1]:]`)
1b. **Fragment Count (circular)**: With k unique top-strand cut positions, produces exactly k fragments; k=0 leaves the circle whole and k=1 linearizes it (Source: Univ. of Illinois MolBio; Addgene Plasmids 101; Biopython `catalyse(..., linear=False)`)
2. **Fragment Sum**: Sum of all fragment lengths == original sequence length (Source: Addgene Protocol — "The sum of the individual fragments should equal the size of the original fragment")
3. **Fragment Order**: Fragments are numbered sequentially from 5' to 3' (Source: Standard molecular biology convention)
4. **Fragment Positions**: Fragment start positions are sorted in ascending order (Source: Standard molecular biology convention)
5. **Boundary Fragments**: First fragment has LeftEnzyme=null, last has RightEnzyme=null (Source: DNA ends at sequence termini have no enzyme cut)
6. **Non-Empty Fragments**: All fragments have Length > 0 (Source: cut positions are distinct and, on a linear molecule, confined to [1, L−1] by the cut-feasibility filter, so every interval is non-empty)
6b. **Contiguous Numbering**: `FragmentNumber` runs 1..n over the returned fragments without a gap (Source: API contract; emission counter, not a boundary index)
7. **No Cuts Returns Whole**: No cut sites returns original sequence as single fragment, in **both** topologies (Source: Standard — no double-strand break = no fragmentation; Biopython `catalyse`: `if not r: return (d[1:],)`, circular included)
7b. **Top-Strand Boundaries**: Fragment boundaries are the top-strand (forward-strand) cuts only; the reverse-strand record of a palindromic site carries the bottom-strand cut of the *same* site and must not be used as a second boundary (Source: Biopython `catalyse` slices at `search()`, which returns top-strand cuts)
7c. **Fragment Order**: Ascending `StartPosition`, so on a circle the origin-spanning fragment is **last**. Biopython emits it first; the two orders differ by a rotation of one (Source: Biopython `catalyse` L1166)
7d. **End Chemistry**: Each end created by a cut reports `OverhangType` (Blunt / FivePrime / ThreePrime) and the extension's bases, read 5'→3' on the strand carrying it and taken from the molecule; a linear molecule's original termini and an uncut circle report `null`. The window is `[min(watson, crick), max(watson, crick))` with `crick = watson + (CutPositionReverse − CutPositionForward)` (Source: Biopython `ovhgseq` / `is_5overhang` / `is_3overhang` / `is_blunt`; 39/39 agreement measured 2026-10-14)
7e. **Ligation Symmetry**: At one cut, the upstream fragment's `RightOverhang` is the reverse complement of the downstream fragment's `LeftOverhang` — identical for a palindromic enzyme (Source: Wikipedia Sticky and blunt ends — overhangs must be complementary to ligate)

### Summary Invariants

8. **Fragment Sizes Sorted**: FragmentSizes in GetDigestSummary are sorted descending (Source: Convention for gel electrophoresis comparison)
9. **Size Bounds**: SmallestFragment ≤ AverageFragmentSize ≤ LargestFragment (Source: Mathematical)
10. **Enzyme List**: EnzymesUsed contains all input enzyme names (Source: API contract)
10b. **Topology Echo**: `DigestSummary.Topology` echoes the topology digested; the topology-free overload is exactly `Linear` (Source: API contract)
10c. **Fragment Sum (both topologies)**: `FragmentSizes` sum to the molecule length for linear and circular alike (Source: Addgene Restriction Digest Protocol)

### Map Invariants

11. **Non-Cutters Identified**: Enzymes with zero sites listed in NonCutters (Source: Wikipedia Restriction map)
12. **Unique Cutters Identified**: Enzymes with exactly one forward-strand site listed in UniqueCutters (Source: Wikipedia Restriction map — forward-strand only to avoid double-counting palindromic sites)
13. **Site Count Consistency**: TotalSites counts forward-strand sites only (Source: Palindromic sites appear on both strands at same position; counting both would double-count)

### Compatibility Invariants

14. **Blunt Compatibility**: All blunt-end enzymes are compatible with each other (Source: Wikipedia Sticky and blunt ends — "blunt ends are always compatible with each other")
15. **Overhang Compatibility**: Enzymes with identical overhang **type** (5' or 3') AND identical overhang **sequence** are compatible (Source: Wikipedia Sticky and blunt ends — "overhangs have to be complementary in order for the ligase to work")
15b. **Cross-Type Incompatibility**: A 5' overhang and a 3' overhang are NOT compatible, even if the overhang sequence string is the same palindrome (Source: Wikipedia — strand geometry prevents base-pairing between same-strand extensions)
16. **Symmetry**: AreCompatible(A, B) == AreCompatible(B, A) (Source: Mathematical)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Single cut produces two fragments | Invariant #1 | Wikipedia |
| M2 | Fragment lengths sum to original sequence length | Invariant #2 | Addgene |
| M3 | No cut sites returns original sequence as single fragment | Invariant #7 | Wikipedia |
| M4 | Multiple cuts produce correct number of fragments | Invariant #1 | Wikipedia |
| M5 | First fragment has LeftEnzyme=null | Invariant #5 | Molecular biology |
| M6 | Last fragment has RightEnzyme=null | Invariant #5 | Molecular biology |
| M7 | Fragment numbers are sequential (1, 2, 3...) | Invariant #3 | Convention |
| M8 | Fragment start positions increase monotonically | Invariant #4 | Convention |
| M9 | Multiple enzymes produce combined cut sites | API contract | Implementation |
| M10 | No enzymes provided throws ArgumentException | API contract | Implementation |
| M11 | Null sequence throws ArgumentNullException | API contract | Implementation |
| M12 | GetDigestSummary fragment sizes sorted descending | Invariant #8 | Convention |
| M13 | GetDigestSummary LargestFragment ≥ SmallestFragment | Invariant #9 | Mathematical |
| M14 | GetDigestSummary EnzymesUsed contains input enzymes | Invariant #10 | API contract |
| M15 | CreateMap identifies non-cutters | Invariant #11 | Wikipedia |
| M16 | CreateMap identifies unique cutters (forward-strand) | Invariant #12 | Wikipedia |
| M17 | CreateMap TotalSites counts forward-strand only | Invariant #13 | Palindromic convention |
| M18 | CreateMap null sequence throws ArgumentNullException | API contract | Implementation |
| M19 | Blunt enzymes compatible with each other | Invariant #14 | Wikipedia |
| M20 | Same-type same-overhang enzymes compatible (BamHI/BglII) | Invariant #15 | Wikipedia |
| M21 | Different-overhang enzymes not compatible | Invariant #15 | Wikipedia |
| M22 | Unknown enzyme in AreCompatible returns false | API contract | Implementation |
| M23 | FindCompatibleEnzymes includes known pairs | Invariant #15 | Wikipedia |
| M24 | Cross-type overhangs NOT compatible (HindIII/SacI) | Invariant #15b | Wikipedia |
| M25 | Cross-type overhangs NOT compatible (SphI/NcoI) | Invariant #15b | Wikipedia |
| M26 | GATC overhang family all compatible (MboI/Sau3AI/BamHI/BglII) | Invariant #15 | Wikipedia |
| M27 | CTAG overhang family all compatible (NheI/XbaI/SpeI/AvrII) | Invariant #15 | Wikipedia |
| M28 | TCGA overhang family compatible (SalI/XhoI) | Invariant #15 | Wikipedia |
| M29 | Enzyme database cut positions match Wikipedia Examples table | Data | Wikipedia |
| M30 | UniqueCutters correctly identifies palindromic single-site enzymes | Invariant #12 | Wikipedia |
| M31 | UniqueCutters excludes multi-site enzymes | Invariant #12 | Wikipedia |
| M32 | FindCompatibleEnzymes excludes cross-type pairs | Invariant #15b | Wikipedia |
| M33 | Circular digest cuts an **origin-spanning** site (`TTCAAAAAAAGAA` + EcoRI → 1 fragment `AATTCAAAAAAAG` at 11; linear is uncut) | Invariant #1b | Biopython `catalyse(..., linear=False)` |
| M34 | Circular digest applies no linear cut-feasibility filter, so an edge cutter cuts every site (`GATCGATCGATC` + Sau3AI → `GATC`×3 at 0/4/8; linear → `GATC` + `GATCGATC`) | Invariant #1b | Biopython `search`/`_drop` |
| M35 | Circular **double** digest with an origin-spanning site (`TTCGGATCCAAAAAAAGAA` + EcoRI/BamHI → `GATCCAAAAAAAG`(4) + `AATTCG`(17)) | Invariant #1b | Biopython `RestrictionBatch.search(..., linear=False)` |
| M36 | Circular multi-enzyme edge cutters (`GATCGATCGATC` + MboI/DpnI → 6 × 2 bp) | Invariant #1b | Biopython batch search |
| M37 | Circular fragments partition the molecule even where Biopython's own `catalyse` does not (unsorted-`search` case `GCTGTTATAGCTAAGAGCGAACCGGAGTATA` + AluI → 2 fragments summing to 31) | Invariant #2 | Biopython `_drop` L542-568 defect, 55/4000 measured |
| M38 | `FragmentNumber` contiguous 1..n, no zero-length fragment, count k+1 / k, ascending starts, sum = length, over 7 dense cut patterns × both topologies | Invariants #6, #6b, #1, #1b, #2 | API contract + cut-feasibility filter |
| M39 | Fragment cut ends report sourced overhang type + sequence (EcoRI 5' `AATT`, PstI 3' `TGCA`, EcoRV blunt empty) and linear termini report `null` | Invariant #7d | Biopython `ovhgseq` |
| M40 | On a circle every end is a cut (both flanks carry enzyme + overhang); an uncut circle reports `null` ends; ligation symmetry holds | Invariants #7d, #7e | Biopython + Wikipedia |
| M41 | The overhang window wraps past the origin (`TTCAAAAAAAGAA` + EcoRI circular → `AATT` from seq[11],seq[12],seq[0],seq[1]) | Invariant #7d | Hand-traced + `ovhgseq` |
| M42 | `GetDigestSummary` topology overload agrees with `Digest`, echoes `Topology`, and the topology-free overload is exactly Linear | Invariants #10b, #1, #1b | API contract |
| M43 | `GetDigestSummary(Circular)` sees an origin-spanning site (sizes `[13, 6]` vs linear `[15, 4]`) | Invariant #1b | Biopython batch search |
| M44 | MCP `restriction_digest` exposes `topology` (default Linear) and reaches the circular origin-spanning case | API contract | Biopython |
| M45 | MCP `restriction_digest` returns the fragment-end overhang fields | Invariant #7d | Biopython `ovhgseq` |
| M46 | MCP `digest_summary` exposes `topology` (default Linear) and echoes it | Invariants #10b, #1b | API contract |
| M47 | Enzymes cutting at the same coordinate produce ONE boundary (`MspI`+`HpaII` on `CCGGAAACCGG` → 3 fragments; `MboI`+`Sau3AI`+`DpnI` on circular `GATCAAAGATC` → 4) with one representative enzyme name and that enzyme's overhang | Invariants #1, #1b, #7d | Biopython `RestrictionBatch.search` |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Fragment sequence content is correct substring | Verify fragment extraction | Convention |
| S2 | Very short sequence (shorter than recognition) returns whole | Edge case | Convention |
| S3 | Adjacent cut sites handled correctly | Edge case | Convention |
| S4 | AreCompatible is symmetric (A,B) == (B,A) | Invariant #16 | Mathematical |
| S5 | CreateMap with no enzymes specified searches all | API contract | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | Performance with many cut sites | Complexity verification | Implementation |
| C2 | Fragment properties include correct enzyme names | Detailed verification | Implementation |

---

## Deviations and Assumptions

**None.** All behavior is derived from external sources:
- Fragment generation, counting, and sizing: Wikipedia (Restriction digest), Addgene Protocol
- Enzyme cut positions and recognition sequences: Wikipedia (Restriction enzyme) Examples table, REBASE
- Overhang compatibility rules: Wikipedia (Sticky and blunt ends)
- Unique cutter identification: Wikipedia (Restriction map)
- Forward-strand-only counting for palindromic sites: standard molecular biology convention for Type II palindromic enzymes

---

## Audit Results

### Test Coverage (RestrictionAnalyzer_Digest_Tests.cs)

| Test | Status | Coverage |
|------|--------|----------|
| Digest_SingleCut_ReturnsTwoFragmentsWithCorrectSum | ✓ | M1, M2 |
| Digest_NoCuts_ReturnsWholeSequenceAsSingleFragment | ✓ | M3 |
| Digest_MultipleCuts_ReturnsCorrectFragmentCount | ✓ | M4 |
| Digest_FragmentsHaveCorrectProperties | ✓ | M5, M6, M7, M8 |
| Digest_FragmentSequenceContent_MatchesExpectedSubstring | ✓ | S1 |
| Digest_MultipleEnzymes_CutsWithBoth | ✓ | M9, C2 |
| Digest_NoEnzymes_ThrowsArgumentException | ✓ | M10 |
| Digest_NullSequence_ThrowsArgumentNullException | ✓ | M11 |
| Digest_SequenceShorterThanRecognition_ReturnsWholeSequence | ✓ | S2 |
| Digest_AdjacentCutSites_ProducesSmallFragment | ✓ | S3 |
| GetDigestSummary_ReturnsCorrectSummaryWithInvariants | ✓ | M12, M13, M14 |
| CreateMap_ReturnsCorrectMapWithAllFields | ✓ | M13, M17, M31 |
| CreateMap_IdentifiesUniqueCutters | ✓ | M16, M30 |
| CreateMap_IdentifiesNonCutters | ✓ | M15 |
| CreateMap_TotalSites_CountsForwardStrandOnly | ✓ | M17 |
| CreateMap_NullSequence_ThrowsArgumentNullException | ✓ | M18 |
| CreateMap_NoEnzymesSpecified_SearchesAll | ✓ | S5 |
| AreCompatible_BluntEnzymes_AreCompatible | ✓ | M19 |
| AreCompatible_SameOverhang_AreCompatible_BamHI_BglII | ✓ | M20 |
| AreCompatible_SameOverhang_AreCompatible_SalI_XhoI | ✓ | M28 |
| AreCompatible_GATC_OverhangFamily_AllCompatible (3 cases) | ✓ | M26 |
| AreCompatible_CTAG_OverhangFamily_AllCompatible (3 cases) | ✓ | M27 |
| AreCompatible_CrossTypeOverhang_NotCompatible_HindIII_SacI | ✓ | M24 |
| AreCompatible_CrossTypeOverhang_NotCompatible_SphI_NcoI | ✓ | M25 |
| AreCompatible_DifferentOverhangs_NotCompatible | ✓ | M21 |
| AreCompatible_UnknownEnzyme_ReturnsFalse | ✓ | M22 |
| AreCompatible_IsSymmetric (5 cases) | ✓ | S4 |
| EnzymeDatabase_MatchesWikipediaData (20 cases) | ✓ | M29 |
| FindCompatibleEnzymes_FindsKnownPairs | ✓ | M23, M32 |
| FindCompatibleEnzymes_AllReturnedPairsAreActuallyCompatible | ✓ | M23 |
| Digest_Circular_OriginSpanningSite_IsCut_MatchesBiopython | ✓ | M33 |
| Digest_Circular_EdgeCutter_CutsEverySite_MatchesBiopython | ✓ | M34 |
| Digest_Circular_DoubleDigest_OriginSpanningSite_MatchesBiopython | ✓ | M35 |
| Digest_Circular_MultiEnzymeEdgeCutters_MatchesBiopython | ✓ | M36 |
| Digest_Circular_BiopythonUnsortedCutCase_FragmentsStillPartitionMolecule | ✓ | M37 |
| Digest_FragmentNumbering_IsContiguous_AndNoZeroLengthFragments | ✓ | M38 |
| Digest_FragmentEnds_ReportSourcedOverhangTypeAndSequence | ✓ | M39 |
| Digest_Circular_EveryEndIsACut_UncutCircleHasNoEnds | ✓ | M40 |
| Digest_Circular_OverhangWindow_WrapsPastTheOrigin | ✓ | M41 |
| GetDigestSummary_TopologyOverload_MatchesDigest_AndDefaultsToLinear | ✓ | M42 |
| GetDigestSummary_Circular_SeesOriginSpanningSite | ✓ | M43 |
| RestrictionDigest_Topology_CircularFindsOriginSpanningSite (MCP) | ✓ | M44 |
| RestrictionDigest_ReportsFragmentEndOverhangs (MCP) | ✓ | M45 |
| DigestSummary_Topology_CircularSeesOriginSpanningSite (MCP) | ✓ | M46 |
| Digest_EnzymesSharingACutCoordinate_ProduceOneBoundary | ✓ | M47 |

### Coverage Classification Summary

All Must (M1–M47), Should (S1–S5), and Could (C2) test cases are covered.
No missing, weak, or duplicate tests remain. M33–M47 were added by the 2026-10-14 re-validation; all
their expected values come from Biopython 1.88 runs recorded in
`docs/Evidence/RESTR-DIGEST-001-Evidence.md`.

---

## Open Questions

None for digestion. Note for `RESTR-FILTER-001`: `GetOverhang(RestrictionEnzyme)` (used only by
`AreCompatible`/`FindCompatibleEnzymes`, out of scope here) derives the overhang from the recognition
**pattern** with `Substring(CutPositionForward, |cutDiff|)` for 5' and
`Substring(CutPositionReverse, |cutDiff|)` for 3' — for the 39 built-in enzymes this coincides with
the window `[min(fwd, rev), max(fwd, rev))` validated here, but it is a second, independent
derivation of the same quantity and is a duplication candidate for that unit.
