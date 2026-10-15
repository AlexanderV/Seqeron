# Test Specification: RESTR-FILTER-001

**Test Unit ID:** RESTR-FILTER-001
**Area:** MolTools
**Algorithm:** Restriction Enzyme Filtering (recognition-site length, blunt/sticky cutters, enzyme-level overhang, end compatibility, isoschizomers)
**Status:** ☑ Re-validated (review-2026-09 batch B08, unit RESTR-FILTER-001)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sticky and blunt ends (Wikipedia) | 4 | https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | 2026-06-13 |
| 2 | Restriction enzyme (Wikipedia) | 4 | https://en.wikipedia.org/wiki/Restriction_enzyme | 2026-06-13 |
| 3 | List of restriction enzyme cutting sites (Wikipedia) | 4 | https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | 2026-06-13 |
| 4 | KpnI (NEB R0142 / REBASE) | 3 | https://www.neb.com/en/products/r0142-kpni | 2026-06-13 |
| 5 | EcoRI-HF (NEB R3101 / REBASE) | 3 | https://www.neb.com/en/products/r3101-ecori-hf | 2026-06-13 |
| 6 | SfiI interrupted palindrome (PMC548270 / REBASE) | 1 | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC548270/ | 2026-06-13 |
| 7 | Biopython 1.88 `Bio.Restriction` (reference implementation, installed locally) — per enzyme `.site`, `.size`, `.fst5`, `.fst3`, `.ovhg`, `.ovhgseq`, `.elucidate()`, `.is_blunt()`, `.is_5overhang()`, `.is_3overhang()`, `.isoschizomers(batch)`, `.equischizomers(batch)`, `.neoschizomers(batch)`, and the `%` compatible-ends operator (`RestrictionType.__mod__` → `Blunt/Ov5/Ov3._mod1` → `Defined/Ambiguous._mod2`) | 2 | `python3 -c "import Bio.Restriction"` (local package; `rebase.neb.com`, `en.wikipedia.org` and `addgene.org` are BLOCKED from this container — the agent proxy answers 403 to CONNECT) | 2026-10-14 |

> **Provenance note (WP-SRC, 2026-10-14).** Sources **1-5** are on hosts denied by the
> organization's egress policy (`en.wikipedia.org`, `www.neb.com` — agent proxy answers 403 to
> CONNECT), and were accessed 2026-06-13, before this campaign. **No unit of this campaign has
> opened them**, so they are **unverifiable in this environment** and must not be read as checked.
> They are retained, not deleted, because the substance of each is independently re-sourced:
> the blunt / 5' / 3' definitions behind sources 1-2 are *executed* against source 7
> (`is_blunt()` / `is_5overhang()` / `is_3overhang()`, 39/39), and the cutting-site data behind
> sources 3-5 is re-sourced from the pinned REBASE **404 (2024)** snapshot and REBASE
> **`allenz.104`**. The NEB product-specific content of sources 4-5 (HF variant, buffers) is **not**
> re-sourced, and nothing in this spec depends on it. See
> [`docs/Evidence/B08-Source-Provenance.md`](../../docs/Evidence/B08-Source-Provenance.md) §1, §7c.

### 1.2 Key Evidence Points

1. A blunt end means "both strands terminate in a base pair"; a sticky/cohesive end is "a stretch of unpaired nucleotides" (a 5' or 3' overhang). Every end is one or the other. — Source 1.
2. Type II enzymes "can either cleave at the center of both strands to yield a blunt end, or at a staggered position leaving overhangs called sticky ends." — Source 2.
3. Type II recognition sites are "4–8 nucleotides in length." — Source 2.
4. SmaI (CCCGGG) is blunt; EcoRI (GAATTC) is a 5' overhang (sticky); KpnI (GGTACC) and PstI (CTGCAG) are 3' overhangs (sticky). — Sources 2, 4, 5.
5. Length categories: EcoRI/BamHI/HindIII/PstI = 6 bp; AluI/HaeIII/TaqI = 4 bp; NotI = 8 bp. — Source 3.
6. SfiI is an interrupted palindrome `5'-GGCCNNNN^NGGCC-3'` (13-nt recognition string), a divided site outside the 4–8 nt undivided range. — Source 6.
7. **Enzyme-level overhang (`ovhgseq`).** The overhang an enzyme is *defined* to leave is the
   top-strand bases between its two strand cuts, written 5'→3': `EcoRI` `"AATT"` (`ovhg -4`),
   `PstI` `"TGCA"` (`ovhg +4`), `NdeI` `"TA"` (`-2`), `TaqI` `"CG"` (`-2`), `PacI` `"AT"` (`+2`),
   `""` for every blunt cutter (`ovhg 0`). Where the site is degenerate *inside* that window the
   reference keeps the IUPAC code: `SfiI` → `"NNN"` (`ovhg +3`). `HincII` (`GTYRAC`) is blunt, so
   its degeneracy never reaches an overhang. — Source 7.
8. **Ligation rule.** Two ends anneal only when their single-stranded extensions have the same
   length, the same polarity (both 5' or both 3') and are reverse complements of each other; blunt
   ends ligate to blunt ends. One cut offers BOTH `W` (the downstream fragment's end, read on the
   top strand) and `revcomp(W)` (the upstream fragment's end, read on the bottom strand), so two
   enzymes are compatible iff they share a polarity and `W1 ∈ { W2, revcomp(W2) }`. Derived from
   the `ovhg`/`ovhgseq` primitives of Source 7 and cross-checked against that source's own `%`
   operator (see §1.5). — Sources 1, 7.
9. **Isoschizomers.** Same recognition site; *equischizomer* = same site **and** same cut
   positions; *neoschizomer* = same site, different cut position. — Source 7
   (`equischizomers`/`neoschizomers` docstrings).

### 1.5 Reference audit performed this session (numbers)

| What was compared | Comparisons | Disagreements |
|---|---|---|
| All 39 built-in entries: `.site` string, `.size`, both cut offsets (`fst5` and `size + fst3`), `is_blunt`, 5'/3' polarity, `ovhgseq` | 39 enzymes × 7 fields | **0** |
| `AreCompatible(name, name)` vs Biopython `a % b or b % a`, all ordered pairs of the 39 built-ins | 1 521 | **0** |
| Same rule vs Biopython `%` over the 389 REBASE enzymes whose `ovhgseq` is blunt or fully A/C/G/T, all unordered pairs | 75 466 | **0** |
| `FindIsoschizomers`/`FindEquischizomers`/`FindNeoschizomers` vs `isoschizomers/equischizomers/neoschizomers(batch = the 39)` | 39 enzymes × 3 queries | **0** |
| `FindCompatibleEnzymes(name)` vs the symmetrised `%` set, per enzyme | 39 | **0** |

**Where Biopython is NOT usable as an oracle** (so the rule is derived from the primitives instead
of copied): `%` is **asymmetric** for 2 386 of 283 881 REBASE pairs, because `_mod1` dispatches on
the left operand and only `Defined._mod2` falls back to the other side; and `Ambiguous._mod2` builds
a broken regular expression when an ambiguity code repeats in the overhang (it splits its own
character class on that code), so `StyI % PscI` (`"CWWG"` vs `"CATG"`) is `False` although the two
overhangs clearly intersect. 12 REBASE pairs are compatible **only** through the
reverse-complement branch that plain `ovhgseq` equality cannot see (`BsmI` `"CN"` with `BsrI`
`"GN"`, and the `CN`/`GN` family generally).

### 1.3 Documented Corner Cases

- A length range excluding all of 4–8 nt returns no enzymes; a range covering 4–8 returns the whole library (Source 2, 3).
- Blunt and sticky sets are complementary and disjoint over the library — the partition is total (Source 1).

### 1.4 Known Failure Modes / Pitfalls

1. Treating a palindromic-but-staggered cutter (KpnI, PstI) as blunt because the overhang string is symmetric — geometry, not the string, defines blunt vs sticky. — Source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetEnzymesByCutLength(int minLength, int maxLength)` | RestrictionAnalyzer | **Canonical** | New range overload; inclusive bounds |
| `GetBluntCutters()` | RestrictionAnalyzer | **Canonical** | Enzymes with blunt ends |
| `GetStickyCutters()` | RestrictionAnalyzer | **Canonical** | Enzymes with 5'/3' overhangs |
| `GetEnzymesByCutLength(int length)` | RestrictionAnalyzer | **Delegate** | Existing single-length overload; smoke + consistency with range overload |
| `GetOverhangSequence(string)` / `GetOverhangSequence(RestrictionEnzyme)` | RestrictionAnalyzer | **Canonical** (new, F28) | REBASE/Biopython `ovhgseq`; `null` from the record overload when the enzyme cuts outside its own site |
| `AreCompatible(string, string)` | RestrictionAnalyzer | **Canonical** | Ligation rule; **unknown name now throws** (F31) |
| `AreCompatible(RestrictionEnzyme, RestrictionEnzyme)` | RestrictionAnalyzer | **Canonical** (new, F30) | Reaches the degenerate and non-palindromic cases the built-in table cannot express |
| `TryAreCompatible(string, string, out bool)` | RestrictionAnalyzer | **Delegate** (new, F31) | Lenient probe kept for callers that cannot guarantee the names exist |
| `FindCompatibleEnzymes()` / `FindCompatibleEnzymes(string)` | RestrictionAnalyzer | **Canonical** | Whole-table pairs / per-enzyme query incl. itself (F32) |
| `FindIsoschizomers` / `FindEquischizomers` / `FindNeoschizomers` | RestrictionAnalyzer | **Canonical** (new, F33) | Biopython `isoschizomers`/`equischizomers`/`neoschizomers` restricted to the built-in table |
| `GetEnzyme(string)` | RestrictionAnalyzer | **Delegate** | `null` for an unknown name; `ArgumentNullException(nameof(name))` for `null` (F31) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Blunt-cutter set and sticky-cutter set are disjoint and their union equals the full enzyme library (total partition). | Yes | Source 1 (every end is blunt or an overhang) |
| INV-2 | `GetBluntCutters()` returns exactly enzymes whose forward and reverse cut positions are equal (center cut → blunt). | Yes | Source 2 (center cut → blunt) |
| INV-3 | `GetEnzymesByCutLength(min,max)` returns exactly enzymes with `min ≤ RecognitionLength ≤ max` (inclusive). | Yes | Source 2/3 (recognition lengths) + ASSUMPTION (bound inclusivity) |
| INV-4 | `GetEnzymesByCutLength(L)` equals `GetEnzymesByCutLength(L,L)` for any L (overload consistency). | Yes | Derivation |
| INV-5 | When `min > max`, the range filter returns an empty set. | Yes | Empty-interval logic |
| INV-6 | `GetOverhangSequence(e)` has length `|cf − cr|` and is empty exactly when `e.IsBluntEnd`. | Yes | Source 7 (`ovhg`/`ovhgseq`) |
| INV-7 | `AreCompatible` is symmetric, and reflexive for every enzyme whose overhang is determined (all 39 built-ins; a distal cutter has no enzyme-level overhang and is compatible with nothing, itself included). | Yes | Derivation (the rule is symmetric by construction) + Source 7 |
| INV-8 | Same overhang letters with opposite polarity are NEVER compatible. | Yes | Sources 1, 7 (`Ov5._mod1`/`Ov3._mod1` require the same class) |
| INV-9 | `FindIsoschizomers = FindEquischizomers ⊎ FindNeoschizomers`, and the queried enzyme is in none of the three. | Yes | Source 7 (`isoschizomers` = equi ∪ neo, `del r[i]`) |
| INV-10 | Every equischizomer leaves the same end, hence is compatible. | Yes | Derivation from INV-6 + the definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | BluntCutters_IncludesKnownBlunt_ExcludesSticky | Blunt set contains SmaI, EcoRV, AluI, HaeIII; excludes EcoRI, KpnI | SmaI/EcoRV/AluI/HaeIII ∈ result; EcoRI/KpnI ∉ result | Source 1, 2 (SmaI blunt) |
| M2 | StickyCutters_IncludesKnownSticky_ExcludesBlunt | Sticky set contains EcoRI, KpnI, PstI; excludes SmaI, EcoRV | EcoRI/KpnI/PstI ∈ result; SmaI/EcoRV ∉ result | Source 2 (EcoRI 5', KpnI/PstI 3'), Source 4 |
| M3 | BluntAndSticky_PartitionLibrary | Blunt ∪ Sticky = all enzymes; Blunt ∩ Sticky = ∅; counts sum to total | Disjoint, union == full library | Source 1 (total partition) |
| M4 | EnzymesByCutLength_Range6_ReturnsAll6Cutters | Range [6,6] returns all and only 6-bp recognition enzymes (EcoRI, BamHI, PstI ∈; AluI, NotI ∉) | All results length 6; EcoRI/BamHI/PstI present; AluI/NotI absent | Source 3 (6-bp set) |
| M5 | EnzymesByCutLength_Range4to8_ReturnsAllUndividedSites_ExcludesSfiI | Range [4,8] returns the full library except the 13-nt interrupted palindrome SfiI | count == total − 1; SfiI absent; all results 4–8 nt | Source 2 (4–8 nt undivided); Source 6 (SfiI interrupted palindrome) |
| M6 | EnzymesByCutLength_Range9to10_ReturnsEmpty | Range above all recognition lengths returns nothing | empty | Source 2 (max 8 nt) |
| M7 | EnzymesByCutLength_Range4_ReturnsOnly4Cutters | Range [4,4] returns only 4-bp enzymes (AluI, HaeIII, TaqI ∈; EcoRI ∉) | All results length 4 | Source 3 (4-bp set) |
| M8 | GetOverhangSequence_MatchesBiopythonOvhgseq | 34 parametric rows: the enzyme-level overhang equals `ovhgseq` | EcoRI `AATT`, PstI `TGCA`, NdeI `TA`, PacI `AT`, SfiI `NNN`, HincII `""`, … | Source 7 |
| M9 | GetOverhangSequence_WholeTable_LengthEqualsCutStaggerAndEmptyOnlyWhenBlunt | All 39: `|overhang| == |cf − cr|`, empty iff blunt | holds for 39/39 | INV-6, Source 7 |
| M10 | EnzymeTable_EveryEntry_MatchesTheReferenceRow | All 39 rows locked: site, size, both cut offsets, polarity, blunt flag, `ovhgseq` | 39/39 match | Source 7 |
| M11 | AreCompatible_SameFivePrimeOverhang_IsCompatible | 11 rows: `GATC` (BamHI/BglII/MboI/Sau3AI), `TCGA` (SalI/XhoI), `CTAG` (NheI/XbaI/SpeI/AvrII), `CG` (ClaI/TaqI/MspI/HpaII) | all `true`, shared end asserted | Sources 1, 2, 7 |
| M12 | AreCompatible_SameLettersOppositePolarity_IsNotCompatible | All four such pairs in the table: HindIII/SacI `AGCT`, NcoI/SphI `CATG`, NotI/ApaI `GGCC`, AgeI/FseI `CCGG` | all `false`, both directions | INV-8, Sources 1, 7 |
| M13 | AreCompatible_BluntAndSticky_FollowTheEndChemistry | blunt+blunt `true`; blunt+5' and blunt+3' `false` | as stated | Source 1, Biopython `Blunt._mod1` |
| M14 | AreCompatible_DifferentOverhangLength_IsNotCompatible | 6 rows incl. SfiI(3 nt)/PstI(4 nt) | all `false` (Biopython `SfiI % PstI` is `False`) | Source 7 |
| M15 | AreCompatible_DegenerateOverhang_MatchesCodeWiseNotAsAString | SfiI `NNN` vs itself and vs a 3-nt 3' `GTC` cutter | both `true` (`SfiI % SfiI` is `True`) | Source 7 (`Ambiguous._mod2`) |
| M16 | AreCompatible_ReverseComplementaryOverhangs_AreCompatible | 5' `CA` with 5' `TG`; 3' `CN` with 3' `GN` (the REBASE BsmI/BsrI shape); 5' `CA` not with 5' `GA` nor with 3' `TG` | as stated | Evidence point 8 (derivation) + the 12 REBASE pairs |
| M17 | FindCompatibleEnzymes_WholeTable_MatchesTheReferenceGrouping | 64 pairs: 45 blunt, 6 `GATC`, 6 `CTAG`, 6 `CG`, 1 `TCGA`, no other end | exact counts | Source 7 (1521/1521) |
| M18 | FindCompatibleEnzymes_PerEnzyme_MatchesTheReferenceSet | 8 rows incl. `SfiI → {SfiI}`, `EcoRV →` the 10 blunt cutters | exact sets | Source 7 |
| M19 | FindCompatibleEnzymes_PerEnzyme_AgreesWithAreCompatible_AndIsSymmetric | 39 per-enzyme sets + all 1521 ordered pairs for symmetry | consistent | INV-7 |
| M20 | FindIsoschizomers_MatchesBiopythonWithinTheBuiltInTable | MspI↔HpaII equi; MboI↔Sau3AI equi; DpnI neo of both; EcoRI/HincII/SfiI none | exact sets | Source 7 (`batch=` the 39) |
| M21 | FindIsoschizomers_WholeTable_PartitionAndEndConsistency | iso = equi ⊎ neo over all 39; equischizomers share the end and are compatible | holds | INV-9, INV-10 |
| M22 | AreCompatible_UnknownEnzyme_Throws | `AreCompatible("EcoRI","UnknownEnzyme")` → `ArgumentException(enzyme2Name)`; `TryAreCompatible` reports `false` instead | throws / `false` | API contract (F31) |
| M23 | AreCompatible_EcoRIAndMfeI_ShareTheFivePrimeAattEnd | The textbook compatible pair. MfeI (`CAATTG`, `C^AATT_G`) is not a built-in, so it is synthesised from its REBASE record as `(1, 5)` — Biopython `charac` `(1, -1, …)` equals EcoRI's. 5' overhang `AATT` for both; EcoRI↔MfeI and MfeI↔MfeI `true` both directions; MfeI↔BamHI `false` both directions | as stated | Source 7 (`MfeI.ovhgseq == "AATT"`, `EcoRI % MfeI`, `MfeI % EcoRI` `True`; `MfeI % BamHI`, `BamHI % MfeI` `False`) — B08 WP4/D5, F42 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | EnzymesByCutLength_MinGreaterThanMax_ReturnsEmpty | min=8, max=4 (inverted) → empty | empty | INV-5 |
| S2 | EnzymesByCutLength_RangeMatchesSingleLengthOverload | Range [L,L] equals single-length overload for L=4,6,8 | set equality | INV-4 |
| S3 | GetOverhangSequence_UnknownEnzyme_Throws | unknown name → `ArgumentException("Unknown enzyme: Nope")` | throws | API contract |
| S4 | GetOverhangSequence_DistalCutter_IsThePatternWindowPaddedWithN; `GetOverhangSequence_DistalEnzyme_MatchesBiopythonOvhgseq` (6), `GetOverhangSequence_TwoCutEnzyme_IsTheFirstPairsWindow`, `AreCompatible_DistalAndTwoCutEnds_MatchBiopythonMod` (13 × 2 orders) | `BsaI`-shaped custom enzyme (`GGTCTC`, cuts 7/11) → `"NNNN"`, `AreCompatible(it, it)` → `true`; BsmI → `CN`, TspRI → `NNCASTGNN`, BaeI → `NNNNN`; BsaI % EcoRI true, BsaI % PstI / SapI / EcoRV false | as stated | Biopython 1.88 `ovhgseq` + `%` (BsaI is `Ambiguous`, not `NotDefined` — superseded F28 by F57) |
| S5 | AreCompatible_NullEnzymeDefinition_Throws | `null` enzyme / `null` name → `ArgumentNullException` | throws | API contract |
| S6 | FindCompatibleEnzymes_PerEnzyme_UnknownName_Throws, FindIsoschizomers_UnknownName_Throws | unknown name → `ArgumentException` | throws | API contract |
| S7 | AreCompatible_UniqueOverhang_OnlyCompatibleWithItself | `FindCompatibleEnzymes("EcoRI")` → `{EcoRI}` | exact set | Source 7. Built-ins only; the off-table partner MfeI is M23 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | EnzymesByCutLength_ZeroOrNegativeBounds_ReturnsEmpty | min=-1, max=0 → empty | empty | No site has length ≤ 0 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_FindSites_Tests.cs` (RESTR-FIND-001) — site finding only.
- `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Digest_Tests.cs` (RESTR-DIGEST-001) — digest/map/compatibility.
- No existing tests cover `GetEnzymesByCutLength`, `GetBluntCutters`, or `GetStickyCutters`. No prior range overload exists.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Filter_Tests.cs` — all filtering tests for RESTR-FILTER-001.
- **Remove:** none (no pre-existing filter tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| RestrictionAnalyzer_Filter_Tests.cs | Canonical for RESTR-FILTER-001 | 33 test methods (12 original + 20 added in review-2026-09 for M8–M22, S3–S7, + 1 for M23 by B08 WP4/D5) |
| RestrictionAnalyzer_Digest_Tests.cs | Holds the pre-existing `AreCompatible`/`FindCompatibleEnzymes` cases from RESTR-DIGEST-001's spec; `AreCompatible_UnknownEnzyme_ReturnsFalse` was renamed to `…_Throws` for F31 | — |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact membership assertions |
| M2 | ✅ Covered | Exact membership assertions |
| M3 | ✅ Covered | Disjoint + union == full library |
| M4 | ✅ Covered | Exact length + membership |
| M5 | ✅ Covered | count == library size − 1; SfiI (interrupted palindrome) excluded |
| M6 | ✅ Covered | empty |
| M7 | ✅ Covered | exact length + membership |
| S1 | ✅ Covered | empty on inverted range |
| S2 | ✅ Covered | set equality with single-length overload |
| C1 | ✅ Covered | empty on non-positive bounds |

**Total in-scope cases:** 10 | **✅:** 10

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Range bounds of `GetEnzymesByCutLength(min,max)` are inclusive on both ends (API-shape convention; recognition-length values themselves are source-backed). | INV-3, M4–M7, S1, S2, C1 |

**Not assumptions, but derivations with the reasoning written out** (see Evidence point 8): the
reverse-complement branch of the compatibility rule and the per-position IUPAC intersection for a
degenerate overhang. Both follow from `ovhg`/`ovhgseq` and the antiparallel base-pairing geometry;
neither could be checked against a published compatible-ends table because `rebase.neb.com`,
`en.wikipedia.org` and `addgene.org` are blocked from this container (403 from the agent proxy on
CONNECT, confirmed three times across this batch's units).

---

## 7. Open Questions / Decisions

1. The single assumption is API-shape only (range bound inclusivity); all biological/recognition-length values are source-backed.
2. **Decided (F31):** `AreCompatible(string, string)` throws on an unknown enzyme name instead of
   returning `false`. The whole class already resolves names that way (`FindSites`, `Digest`,
   `CreateMap`, and now `GetOverhangSequence`/`FindCompatibleEnzymes(name)`/the isoschizomer
   queries), and "not compatible" is the answer that silently loses a real cloning option. The only
   in-repo caller was the MCP `enzymes_compatible` tool, which now surfaces the error — consistent
   with `find_restriction_sites`. `TryAreCompatible` keeps the lenient path available explicitly.
3. **Decided:** `GetEnzymesByCutLength` keeps measuring the full recognition string, so `SfiI` is a
   13 and not an 8 — this is Biopython `.size` and the previous behaviour; the XML doc and the MCP
   tool description now say so, because "8-cutter" in catalogue prose means SfiI's eight *specified*
   bases.
