# Evidence Artifact: RESTR-FILTER-001

**Test Unit ID:** RESTR-FILTER-001
**Algorithm:** Restriction Enzyme Filtering (recognition-site length, blunt/sticky cutters, enzyme-level overhang, end compatibility, isoschizomers)
**Date Collected:** 2026-06-13; **re-validated and extended 2026-10-14** (review-2026-09 batch B08, unit RESTR-FILTER-001 — see the addendum at the end)

---

## Online Sources

### Sticky and blunt ends (Wikipedia)

**URL:** https://en.wikipedia.org/wiki/Sticky_and_blunt_ends
**Accessed:** 2026-06-13 (fetched via WebFetch of the article URL)
**Authority rank:** 4 (Wikipedia citing primary literature; used for definitions, not numeric parameters)

**Key Extracted Points:**

1. **Blunt end definition:** Verbatim — "In a blunt-ended molecule, both strands terminate in a base pair." (No unpaired overhang nucleotides.)
2. **Blunt-end compatibility:** Verbatim — "blunt ends are always compatible with each other."
3. **Overhang / sticky end definition:** Verbatim — "An overhang is a stretch of unpaired nucleotides in the end of a DNA molecule." and "Longer overhangs are called cohesive ends or sticky ends."
4. **Overhang strand:** Verbatim — "These unpaired nucleotides can be in either strand, creating either 3' or 5' overhangs." A sticky end is therefore a non-blunt end (a 5' or a 3' overhang).

### Restriction enzyme (Wikipedia)

**URL:** https://en.wikipedia.org/wiki/Restriction_enzyme
**Accessed:** 2026-06-13 (fetched via WebFetch of the article URL)
**Authority rank:** 4 (Wikipedia citing primary literature; recognition-length range cross-checked below)

**Key Extracted Points:**

1. **Type II recognition-site length:** Verbatim — Type II enzymes have "recognition sites that are usually undivided and palindromic and 4–8 nucleotides in length." This bounds the meaningful range for length-based filtering.
2. **Blunt vs sticky from cut geometry:** Type II enzymes "can either cleave at the center of both strands to yield a blunt end, or at a staggered position leaving overhangs called sticky ends." A symmetric (center) cut → blunt; a staggered cut → sticky.
3. **EcoRI example (5' overhang / sticky):** Recognizes GAATTC and cuts `5'---G  AATTC---3'`, producing a 5' overhang (sticky end).
4. **SmaI example (blunt):** Recognizes CCCGGG and produces blunt ends `5'---CCC  GGG---3'` (cut at the center of the palindrome).
5. **KpnI example (3' overhang / sticky):** Recognizes GGTACC and cuts `5'---GGTAC  C---3'`, leaving a 3' overhang (sticky end).
6. **PstI example (3' overhang / sticky):** Recognizes CTGCAG, produces a 3' overhang `5'---CTGCA  G---3'`.

### List of restriction enzyme cutting sites (Wikipedia)

**URL:** https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites
**Accessed:** 2026-06-13 (fetched via WebFetch of the article URL)
**Authority rank:** 4

**Key Extracted Points:**

1. **Recognition-length categories:** EcoRI, BamHI, HindIII, XhoI, SacI, PstI, NdeI recognize 6-bp sequences; TaqI, AluI, HaeIII recognize 4-bp sequences; NotI recognizes an 8-bp sequence. This corroborates the 4/6/8-bp lengths used by the length-filter test cases.
2. **End type:** Enzymes "produce either blunt or overhung ends" (blunt ends marked with asterisks in the source's resource lists).

### SfiI — interrupted palindrome (PMC / REBASE, search result text)

**URL:** retrieved via WebSearch query `SfiI GGCCNNNNNGGCC interrupted palindrome recognition site 13 bp REBASE NEB` (results: ncbi.nlm.nih.gov/pmc/articles/PMC548270, medical-dictionary "interrupted palindrome")
**Accessed:** 2026-06-13
**Authority rank:** 1/5 (peer-reviewed PMC article) / 5 (REBASE)

**Key Extracted Points:**

1. **SfiI recognition site:** "SfiI … recognizes the interrupted palindromic sequence 5'GGCCNNNN^NGGCC3'" — a 13-nt string with a 5-base degenerate (NNNNN) spacer. This is a **divided / interrupted** palindrome, not the undivided 4–8 nt class of Source 2; its recognition-string length is 13 and therefore lies outside the inclusive [4,8] length filter.

### KpnI — NEB / REBASE (search result text)

**URL:** retrieved via WebSearch query `KpnI GGTAC^C 3' overhang sticky end recognition site NEB REBASE` (top results: neb.com/en/products/r0142-kpni, rebase.neb.com)
**Accessed:** 2026-06-13
**Authority rank:** 3/5 (NEB product page / REBASE database)

**Key Extracted Points:**

1. **KpnI cut/overhang:** "KpnI recognizes the palindromic … sequence 5′-GGTACC-3′ and generates a 3′, 4-base overhang … notated as G_GTAC^C." Confirms KpnI is a sticky (3' overhang) cutter, not blunt.

### EcoRI — NEB / REBASE (search result text)

**URL:** retrieved via WebSearch query `REBASE EcoRI GAATTC cut site G^AATTC blunt or 5' overhang NEB` (top result: neb.com/en/products/r3101-ecori-hf, rebase.neb.com)
**Accessed:** 2026-06-13
**Authority rank:** 3/5 (NEB / REBASE)

**Key Extracted Points:**

1. **EcoRI cut/overhang:** EcoRI "recognizes the palindromic sequence GAATTC and cuts between the G and the A on both … strands, leaving … a sticky end" — a 4-base 5' AATT overhang. Confirms EcoRI is a sticky (5' overhang) cutter.

---

## Documented Corner Cases and Failure Modes

### From Restriction enzyme (Wikipedia)

1. **Recognition length outside 4–8 nt:** A length-filter range that excludes all of 4–8 (e.g. min=9) must return no enzymes; a range covering 4–8 returns the whole library. Lengths below 4 or above 8 are not meaningful Type II recognition sites.

### From Sticky and blunt ends (Wikipedia)

1. **Blunt/sticky partition is total:** Every enzyme end is either blunt (both strands terminate in a base pair) or a sticky end (a 5' or 3' overhang). There is no third category, so the blunt-cutter set and sticky-cutter set are complementary and disjoint over the library.

---

## Test Datasets

### Dataset: Library enzymes with known end type and recognition length

**Source:** Restriction enzyme (Wikipedia) §"Types"; KpnI NEB/REBASE; EcoRI NEB/REBASE; List of restriction enzyme cutting sites (Wikipedia).

| Enzyme | Recognition | Length (bp) | Cut geometry | End type |
|--------|-------------|-------------|--------------|----------|
| SmaI | CCCGGG | 6 | center (CCC^GGG) | Blunt |
| EcoRV | GATATC | 6 | center (GAT^ATC) | Blunt |
| AluI | AGCT | 4 | center (AG^CT) | Blunt |
| HaeIII | GGCC | 4 | center (GG^CC) | Blunt |
| EcoRI | GAATTC | 6 | staggered, 5' overhang (G^AATTC) | Sticky (5') |
| KpnI | GGTACC | 6 | staggered, 3' overhang (GGTAC^C) | Sticky (3') |
| PstI | CTGCAG | 6 | staggered, 3' overhang (CTGCA^G) | Sticky (3') |
| NotI | GCGGCCGC | 8 | staggered | Sticky |
| TaqI | TCGA | 4 | staggered | Sticky |
| SfiI | GGCCNNNN^NGGCC | 13 (interrupted) | staggered, 3' overhang | Sticky (divided palindrome; outside 4–8 nt) |

---

## Assumptions

1. **ASSUMPTION: Range filter is inclusive on both bounds** — The checklist names `GetEnzymesByCutLength(minLength, maxLength)`. No authoritative source defines the inclusivity of a software filter's range bounds. The implemented contract treats both bounds as inclusive (`minLength ≤ len ≤ maxLength`), which is the conventional meaning of a min/max range. This is an API-shape decision, not a biological-correctness parameter: the recognition-length values themselves are source-backed; only the boundary convention of the helper is conventional. Documented in the algorithm doc and TestSpec.

---

## Recommendations for Test Coverage

1. **MUST Test:** `GetBluntCutters()` returns only enzymes whose two strands cut symmetrically (blunt: SmaI, EcoRV, AluI, HaeIII present) and excludes sticky cutters (EcoRI, KpnI). — Evidence: Sticky and blunt ends (Wikipedia), Restriction enzyme (Wikipedia) SmaI example.
2. **MUST Test:** `GetStickyCutters()` returns only overhang-producing enzymes (EcoRI, KpnI, PstI present) and excludes blunt cutters (SmaI). — Evidence: Restriction enzyme (Wikipedia) EcoRI/KpnI/PstI examples; KpnI NEB/REBASE.
3. **MUST Test:** Blunt and sticky sets are disjoint and together equal the full library (total partition). — Evidence: Sticky and blunt ends (Wikipedia) — every end is blunt or an overhang.
4. **MUST Test:** `GetEnzymesByCutLength(min,max)` returns exactly the enzymes with recognition length in `[min,max]` inclusive (e.g. 6→6 yields all 6-cutters; 9→10 yields none). The `[4,8]` range yields the full library **except** the interrupted-palindrome SfiI (recognition string `GGCCNNNNNGGCC`, length 13), which lies outside the 4–8 nt undivided-site range and is correctly excluded. — Evidence: Restriction enzyme (Wikipedia) 4–8 nt range (undivided sites); List of restriction enzyme cutting sites (Wikipedia) length categories; SfiI PMC/REBASE (interrupted palindrome).
5. **SHOULD Test:** Range with `min > max` returns empty (no enzyme length can satisfy an empty interval). — Rationale: boundary/error behavior of a numeric range.
6. **SHOULD Test:** Existing single-length overload `GetEnzymesByCutLength(length)` agrees with the range overload at `min == max == length`. — Rationale: consistency between the two overloads.
7. **COULD Test:** Negative or zero bounds return empty (no recognition site has length ≤ 0). — Rationale: defensive boundary.

---

## References

1. Wikipedia. 2026. Sticky and blunt ends. https://en.wikipedia.org/wiki/Sticky_and_blunt_ends
2. Wikipedia. 2026. Restriction enzyme. https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia. 2026. List of restriction enzyme cutting sites. https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites
4. New England Biolabs / REBASE. KpnI (R0142). https://www.neb.com/en/products/r0142-kpni ; http://rebase.neb.com
5. New England Biolabs / REBASE. EcoRI-HF (R3101). https://www.neb.com/en/products/r3101-ecori-hf ; http://rebase.neb.com
6. SfiI interrupted palindrome `5'-GGCCNNNN^NGGCC-3'` (homology model article). PMC. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC548270/

---

---

# Addendum — review-2026-09 / B08, unit RESTR-FILTER-001 (2026-10-14)

Everything in this addendum was computed in that session; no number is carried over.

## A1. Blocked hosts — what could NOT be opened, and what was used instead

Three hosts that the 2026-06 evidence above cites are unreachable from this container. Verified
this session: `curl` exits with status `000`, and `$HTTPS_PROXY/__agentproxy/status` lists the
denials verbatim —

```
2026-10-14T19:18:52.434Z connect_rejected rebase.neb.com:443   gateway answered 403 to CONNECT (policy denial or upstream failure)
2026-10-14T19:18:52.700Z connect_rejected en.wikipedia.org:443 gateway answered 403 to CONNECT (policy denial or upstream failure)
2026-10-14T19:18:53.016Z connect_rejected www.addgene.org:443  gateway answered 403 to CONNECT (policy denial or upstream failure)
```

So **no published compatible-ends table (REBASE, NEB, Addgene) and no Wikipedia page was opened in
this session.** The quotations in the 2026-06 sections above are left exactly as they were
recorded then and are **not** re-asserted here.

**Update (WP-SRC, 2026-10-14).** The denials above still stand and were not re-attempted. What has
changed is that REBASE itself is partly obtainable *without* `rebase.neb.com`: besides the EMBOSS
**404 (2024)** snapshot in Biopython, REBASE **`allenz.104`** (version 104, Copyright Dr. Richard J.
Roberts 2021, released Mar 31 2021) is vendored in the Bioconductor package `crisprVerse/crisprBase`
and served from the allowed host `raw.githubusercontent.com` — which supplies the organism and
methylation columns the Biopython snapshot omits. That does **not** yield a compatible-ends table
(REBASE publishes none), so B1 stays blocked. Datasets, sha256 hashes, mirror cross-check numbers
and a per-citation reachability verdict for this file's sources 1-5 are in
[`B08-Source-Provenance.md`](B08-Source-Provenance.md). What this session *could* read is the locally
installed reference implementation:

**Oracle: Biopython 1.88 `Bio.Restriction`** (`python3 -c "import Bio; print(Bio.__version__)"` →
`1.88`). Primitives used, per enzyme: `.site`, `.size`, `.fst5`, `.fst3`, `.ovhg`, `.ovhgseq`,
`.elucidate()`, `.is_blunt()`, `.is_5overhang()`, `.is_3overhang()`, `.isoschizomers(batch)`,
`.equischizomers(batch)`, `.neoschizomers(batch)`, and the `%` compatible-ends operator
(`Restriction.py`: `RestrictionType.__mod__` → `Blunt._mod1` / `Ov5._mod1` / `Ov3._mod1` →
`Defined._mod2` / `Ambiguous._mod2` / `NotDefined._mod2`), all read by `inspect.getsource`.

## A2. Enzyme-table audit — all 39 built-ins, seven fields each

Method: the `_enzymes` initialiser in `RestrictionAnalyzer.cs` was parsed mechanically and each row
compared to the same enzyme in `Bio.Restriction`. Biopython's bottom-strand cut offset is
`size + fst3`; its top-strand offset is `fst5`.

| Field compared | Seqeron value | Biopython value | Entries | Wrong |
|---|---|---|---|---|
| Recognition string | `RecognitionSequence` | `.site` | 39 | **0** |
| Site size class | `RecognitionLength` | `.size` | 39 | **0** |
| Top-strand cut | `CutPositionForward` | `.fst5` | 39 | **0** |
| Bottom-strand cut | `CutPositionReverse` | `.size + .fst3` | 39 | **0** |
| Blunt vs sticky | `IsBluntEnd` | `.is_blunt()` | 39 | **0** |
| 5' vs 3' polarity | `OverhangType` | `.is_5overhang()` / `.is_3overhang()` | 39 | **0** |
| Overhang sequence | `GetOverhangSequence` | `.ovhgseq` | 39 | **0** |

Class counts, both sides identical: **10 blunt**, **21 with a 5' overhang**, **8 with a 3'
overhang**. Distinct ends in the table: 5' `AATT, AGCT, CATG, CCGG, CG, CGCG, CTAG, GATC, GGCC, TA,
TCGA`; 3' `AGCT, AT, CATG, CCGG, GGCC, GTAC, NNN, TGCA`; blunt `""`.

Representative rows (Biopython `.elucidate()` in brackets):

| Enzyme | Site | cf / cr | `ovhg` | `ovhgseq` | `elucidate()` |
|---|---|---|---|---|---|
| EcoRI | GAATTC | 1 / 5 | −4 | `AATT` | `G^AATT_C` |
| BamHI | GGATCC | 1 / 5 | −4 | `GATC` | `G^GATC_C` |
| Sau3AI | GATC | 0 / 4 | −4 | `GATC` | `N^GATC_N` |
| NdeI | CATATG | 2 / 4 | −2 | `TA` | `CA^TA_TG` |
| TaqI | TCGA | 1 / 3 | −2 | `CG` | `T^CG_A` |
| PstI | CTGCAG | 5 / 1 | +4 | `TGCA` | `C_TGCA^G` |
| KpnI | GGTACC | 5 / 1 | +4 | `GTAC` | `G_GTAC^C` |
| PacI | TTAATTAA | 5 / 3 | +2 | `AT` | `TTA_AT^TAA` |
| SfiI | GGCCNNNNNGGCC | 8 / 5 | +3 | `NNN` | `GGCCN_NNN^NGGCC` |
| HincII | GTYRAC | 3 / 3 | 0 | `""` | `GTY^_RAC` |
| DpnI | GATC | 2 / 2 | 0 | `""` | `GA^_TC` |
| SwaI | ATTTAAAT | 4 / 4 | 0 | `""` | `ATTT^_AAAT` |

**On `SfiI` specifically:** `"NNN"` is Biopython's own `ovhgseq`, not an artefact of slicing the
pattern. The enzyme-level overhang of a divided palindrome genuinely *is* three unspecified bases;
the per-site bases are obtainable only from a molecule, which is what the digest's
`CutOverhangWindow` reports (`DigestFragment.LeftOverhang` / `RightOverhang`). The pattern-derived
value was therefore **not** a defect — the inherited note that it was "the literal `NNN` rather than
the real per-site overhang" does not reproduce as a wrong value.

**On `GetEnzymesByCutLength`:** the length measured is the full recognition string, i.e. `.size`.
`SfiI` is therefore a 13, not an 8, which is Biopython's classification too (`SfiI.size == 13`). The
XML doc and the MCP tool description now state this, because catalogue prose calls SfiI an
"8-base-specificity" enzyme.

## A3. The ligation rule, derived from `ovhg` / `ovhgseq`

With `W` = the top-strand bases between the two strand cuts (`ovhgseq`):

1. A cut leaves two ends. The downstream fragment's end carries the extension `W` read 5'→3' on the
   strand that protrudes; the upstream fragment's end carries `revcomp(W)` on the opposite strand.
   They are complementary and antiparallel — which is exactly the geometry T4 ligase needs — so they
   re-ligate. **An enzyme therefore offers both `W` and `revcomp(W)`.**
2. Two ends anneal iff the extensions have equal length, equal polarity and one is the reverse
   complement of the other. Opposite polarity cannot anneal: both extensions would lie on the same
   strand of the joint.
3. Hence **E1 and E2 are compatible iff they share a polarity and `W1 ∈ { W2, revcomp(W2) }`**;
   blunt ligates to blunt.
4. An ambiguity code denotes a set of bases, independently per position, so `expand(W)` is a
   Cartesian product and `expand(revcomp(W)) = revcomp(expand(W))`. The condition in 3 becomes
   `expand(W1) ∩ expand(W2) ≠ ∅` **or** `expand(W1) ∩ expand(revcomp(W2)) ≠ ∅`, which is exactly
   per-position code intersection. A degenerate verdict is therefore "**possibly** compatible —
   some pair of real sites ligates".

Implemented as `OverhangsAnneal` + `CodesIntersect`, the latter composed from the canonical
`IupacHelper.MatchesIupac` over `A/C/G/T` (no second ambiguity table in `RestrictionAnalyzer`).

## A4. Cross-check of that rule against Biopython `%`

| Comparison | Pairs | Disagreements |
|---|---|---|
| All ordered pairs of the 39 built-ins, `AreCompatible(a,b)` vs `a % b or b % a` | 1 521 | **0** |
| All unordered pairs of the 389 REBASE enzymes whose `ovhgseq` is blunt or fully A/C/G/T, vs `a % b` | 75 466 | **0** |
| Same 75 466 pairs with the reverse-complement branch disabled | 75 466 | **0** (no defined REBASE pair is an rc-partner without also being equal) |
| `FindCompatibleEnzymes(name)` vs the symmetrised `%` set | 39 enzymes | **0** |
| `FindIsoschizomers` / `FindEquischizomers` / `FindNeoschizomers` vs `*(batch = the 39)` | 117 queries | **0** |

**Why Biopython's `%` is used as a cross-check and not copied.** Two defects were found in it this
session, both by reading `Restriction.py` and reproducing them:

- **Asymmetry.** `RestrictionType.__mod__` delegates to `cls._mod1(other)`, i.e. it dispatches on
  the *left* operand, and only `Defined._mod2` contains a `other._mod2(cls)` fallback —
  `Ambiguous._mod2` has none. So an `Ambiguous` × `Ambiguous` pair is judged from one side only:
  `PciSI % Psp5II` is `True` while `Psp5II % PciSI` is `False` (`"NNN"` vs `"GWC"`).
  **2 386 of 283 881 REBASE pairs disagree with their own transpose.**
- **Broken regex for a repeated ambiguity code.** `Ambiguous._mod2` expands the first occurrence of
  a code into a character class and then, on the next iteration of the same code, calls
  `se.split(base)` on a string that now contains that code *inside the class it just inserted*
  (`matching['W'] == 'ABDHKMNRTWVY'`). `"CWWG"` becomes
  `C[ABDHKMNRT[ABDHKMNRTWVY]VY][ABDHKMNRT[ABDHKMNRTWVY]VY]G`, which matches nothing — so
  `StyI % PscI` and `PscI % StyI` are both `False` although `"CWWG"` and `"CATG"` plainly intersect.
  (`"NNN"` survives because `matching['N']` is substituted via the separate `"."` branch.)

Both make `%` unsound exactly where ambiguity is involved, so the rule above is derived instead.
The implemented rule is symmetric by construction and agrees with `%` on every pair where `%` is
sound (the 75 466 defined pairs and all 1 521 built-in pairs).

## A5. Confirmed compatible / incompatible pairs

Compatible (`AreCompatible` `true`; Biopython `%` `True`):

| Pair | Shared end | Polarity |
|---|---|---|
| BamHI / BglII / MboI / Sau3AI (all 6 pairs) | `GATC` | 5' |
| SalI / XhoI | `TCGA` | 5' |
| NheI / XbaI / SpeI / AvrII (all 6 pairs) | `CTAG` | 5' |
| ClaI / TaqI / MspI / HpaII (all 6 pairs) | `CG` | 5' |
| any two of the 10 blunt cutters (all 45 pairs) | `blunt` | — |
| SfiI with itself | `NNN` | 3' |

Incompatible, and **why**:

| Pair | Ends | Reason |
|---|---|---|
| HindIII / SacI | 5' `AGCT` vs 3' `AGCT` | **opposite polarity**, same letters |
| NcoI / SphI | 5' `CATG` vs 3' `CATG` | opposite polarity |
| NotI / ApaI | 5' `GGCC` vs 3' `GGCC` | opposite polarity |
| AgeI / FseI | 5' `CCGG` vs 3' `CCGG` | opposite polarity |
| SfiI / PstI, SfiI / ApaI, SfiI / KpnI | 3' `NNN` (3 nt) vs 3' 4-nt ends | **length** (Biopython agrees: `SfiI % PstI` is `False`) |
| EcoRI / anything but itself | 5' `AATT` | no other built-in leaves `AATT` (REBASE `MfeI`, `CAATTG`, would — it is not a built-in) |
| EcoRV / EcoRI, SmaI / PstI | blunt vs sticky | a blunt end has no extension to anneal |

All four same-letters-opposite-polarity pairs that exist in the table are asserted in
`AreCompatible_SameLettersOppositePolarity_IsNotCompatible`, in both argument orders.

**The reverse-complement branch** is not reachable from the built-in table: every built-in overhang
window is its own reverse complement, because a palindromic site cut symmetrically always yields a
palindromic window (verified: 0 of the 39 windows is non-palindromic). It is real in REBASE,
though — **17 REBASE enzymes have a non-palindromic `ovhgseq`** (`BbvCI` `TCA` 5'; `BsiI`/`BauI`/
`BssSI`/`Bst2BI` `ACGA` 5'; `PspFI`/`BseYI` `CCAG` 5'; `SimI` `GTC` 5'; `GsaI` `CCAG` 3'; `BscXI`
`CAGG` 3'; `BsmI`/`Mva1269I`/`PctI`/`LmnI` `CN` 3'; `BsrI`/`BseNI`/`Bse1I` `GN` 3') and **12 REBASE
pairs are compatible only through it**: `{BsmI, Mva1269I, LmnI, PctI}` (`CN`) × `{BsrI, BseNI,
Bse1I}` (`GN`). Worked instantiation: `BsmI` can leave a 3' `GG` (its `revcomp(CN)` end with N = C)
and `BsrI` a 3' `CC` (its `revcomp(GN)` end with N = G); `revcomp(GG) == CC`, so they ligate. It is
reachable in this library through the `AreCompatible(RestrictionEnzyme, RestrictionEnzyme)` overload
and is tested with that shape.

## A6. Isoschizomers (Biopython definitions, `batch` restricted to the 39 built-ins)

| Enzyme | isoschizomers | equischizomers (same cut) | neoschizomers (different cut) |
|---|---|---|---|
| MspI (CCGG 1/3) | HpaII | HpaII | — |
| HpaII (CCGG 1/3) | MspI | MspI | — |
| MboI (GATC 0/4) | DpnI, Sau3AI | Sau3AI | DpnI |
| Sau3AI (GATC 0/4) | DpnI, MboI | MboI | DpnI |
| DpnI (GATC 2/2) | MboI, Sau3AI | — | MboI, Sau3AI |
| every other built-in | — | — | — |

Identical on both sides for all 39 enzymes and all three queries (117/117).

## A7. Change log of this unit (see `docs/Validation/review-2026-09/B08.md` F28–F34)

- **F28** one canonical overhang-window geometry (`OverhangWindow`), public
  `GetOverhangSequence`; the duplicate pattern-slicing derivation inside `GetOverhang` is gone and
  `CutOverhangWindow` now reads the same geometry off the molecule.
- **F29/F30** the derived ligation rule, with the reverse-complement branch and per-position IUPAC
  intersection; public `AreCompatible(RestrictionEnzyme, RestrictionEnzyme)`.
- **F31** `AreCompatible(string, string)` and `GetEnzyme(null)` exception contracts.
- **F32** `FindCompatibleEnzymes(string)` per-enzyme query.
- **F33** `FindIsoschizomers` / `FindEquischizomers` / `FindNeoschizomers`.
- **F34** the whole-table audit locked in a test.

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-10-14**: Addendum A1–A7 — Biopython 1.88 oracle, 39-enzyme seven-field audit (0 wrong
  entries), the derived ligation rule and its 1 521 + 75 466-pair cross-check, two defects found in
  Biopython's own `%` operator, the confirmed compatible/incompatible pair table, and the
  isoschizomer sets. Blocked hosts recorded instead of cited.
