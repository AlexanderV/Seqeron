# Evidence Artifact: RESTR-FILTER-001

**Test Unit ID:** RESTR-FILTER-001
**Algorithm:** Restriction Enzyme Filtering (by recognition-site length, blunt cutters, sticky cutters)
**Date Collected:** 2026-06-13

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

## Change History

- **2026-06-13**: Initial documentation.
