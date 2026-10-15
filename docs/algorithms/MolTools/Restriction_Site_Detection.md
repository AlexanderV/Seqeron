# Restriction Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | N/A |
| Related Projects | N/A |
| Implementation Status | Reference-verified |
| Test Unit ID (actual) | RESTR-FIND-001 |
| Last Reviewed | 2026-10-14 |
| Reference Oracle | Biopython 1.88 `Bio.Restriction` (enzyme table = REBASE EMBOSS files version 404, 2024) |

## 1. Overview

Restriction site detection identifies DNA sequences recognized by restriction endonucleases and computes the corresponding cut positions. In this repository, the implementation scans both strands for recognition-sequence matches, supports IUPAC ambiguity codes in enzyme patterns, and returns site records that include the enzyme, strand, cut position, and matched sequence. The current implementation is sequence-based and centered on a built-in enzyme catalog.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Restriction enzymes cleave DNA at specific recognition sequences, and Type II enzymes are the standard molecular-biology workhorses because they cut at or near their recognition sites. Recognition sequences are commonly 4-8 bases long, and many are palindromic so the same site appears on both strands. The original document also classifies cleavage products into 5' overhangs, 3' overhangs, and blunt ends. Sources: Wikipedia (Restriction enzyme, Restriction site, EcoRI), Roberts (1976), IUPAC-IUB (1970).

### 2.2 Core Model

For each enzyme, the algorithm slides a window of the recognition-sequence length across the input sequence and checks character-by-character matches under IUPAC ambiguity rules. All reported coordinates are 0-based and in **forward-strand** space, for both strands; `CutPosition` is the number of forward-strand bases 5' of the cut, i.e. the cut lies between bases `CutPosition - 1` and `CutPosition`, **on the strand named by `IsForwardStrand`**. Biopython `Bio.Restriction.<enzyme>.search(seq)` reports the same top-strand cut 1-based, so

$$
biopython = CutPosition_{forward} + 1 .
$$

A forward-strand match at position `i` has cuts

$$
watson = i + CutPositionForward, \qquad crick = i + CutPositionReverse ,
$$

and the forward-strand record reports `watson`.

Reverse-strand matches are found by scanning the reverse complement and converting the recognition-site start back to forward coordinates,

$$
forwardPos = sequenceLength - i - patternLength .
$$

For such a site the enzyme is bound the other way round, so with `L` the recognition length its cuts are

$$
watson = forwardPos + L - CutPositionReverse, \qquad crick = forwardPos + L - CutPositionForward ,
$$

and the reverse-strand record reports `crick` — the cut on the strand it read — together with the matched window **as read on that strand** (the reverse complement of the forward-strand window). For a palindromic enzyme `CutPositionForward + CutPositionReverse = L`, so `crick = forwardPos + CutPositionReverse` and the two formulations coincide; every built-in enzyme is palindromic, so this distinction is visible only for a custom non-palindromic enzyme. (Reference: Biopython `Defined._modify` / `Defined._rev_modify`, `Palindromic._search`, `NonPalindromic._search`.)

### 2.2a Distal (single-sided Type IIS) cutters

Nothing in §2.2 assumes the cut lies inside the window: the offsets are signed and may be negative (cut upstream of the site) or exceed `L` (cut downstream). The record maps onto Biopython's `fst5`/`fst3` as

$$
CutPositionForward = fst5, \qquad CutPositionReverse = size + fst3, \qquad ovhg = CutPositionForward - CutPositionReverse ,
$$

and REBASE `SITE(n/m)` onto `(L + n, L + m)`: `BsaI` GGTCTC(1/5) = `("GGTCTC", 7, 11)`, `BsmBI` = `("CGTCTC", 7, 11)`, `BbsI` GAAGAC(2/6) = `("GAAGAC", 8, 12)`, `SapI` GCTCTTC(1/4) = `("GCTCTTC", 8, 11)`, `FokI` GGATG(9/13) = `("GGATG", 14, 18)`, `MlyI` GAGTC(5/5) = `("GAGTC", 10, 10)` (blunt). §2.2's formulas are then exactly Biopython's `location + fst5` (forward) and `location − fst3` (reverse), so for a non-palindromic site the reverse-strand occurrence cuts on the **other** side of the window (BsaI on `GAGACC` cuts 5 nt 5' of it on the top strand). Validated by a differential over all 153 single-sided distal enzymes of the Biopython 1.88 snapshot × random sequences × both topologies: 0 mismatches after the fix of review B08 F54 (the circular reverse-strand cut was reduced with C#'s signed `%` and could be negative; 2,575 affected cases before). Enzymes that cut on **both** sides (Type IIB `BaeI`, `BcgI`; Biopython `TwoCuts`) carry their second pair in the additive `RestrictionEnzyme.SecondCut = (scd5, size + scd3)` (BaeI = `("ACNNNNGTAYC", -10, -15) { SecondCut = (23, 18) }`); every site is then reported once per cut (`RestrictionSite.IsSecondCut`), with Biopython's `_modify` (`location + scd5`) / `_rev_modify` (`location − scd3`), each cut filtered on its own on a linear molecule. Both pairs must leave the same stagger (Biopython stores one `ovhg`; true for all 25). Differential over all 25 two-cut enzymes × both topologies: 0 mismatches (review B08 F56).

One documented divergence: Biopython 1.88 stores `TscAI`/`TspRI` (CASTG, `fst5 = 7`, `fst3 = −7`) with `ovhg = 10`, contradicting their own `fst5`/`fst3` and 9-nt `ovhgseq` `NNCASTGNN`; its `_drop` therefore rejects a linear site whose bottom-strand cut is the (existing) bond after base 0 (`TspRI.search(Seq("AAACAGTGAAAAAA")) == []`). This library uses the offsets (`ovhg = 9`) and keeps the site.

### 2.3 Cut Feasibility and Topology

A cut can only fall on a phosphodiester bond that exists. On a **linear** molecule of length `n` a site is therefore reported only when

$$
1 \le watson \le n - 1 \quad\text{and}\quad 1 \le crick \le n - 1 ,
$$

which is Biopython's `AbstractCut._drop` filter (`1 < cut_on_watson <= length` together with the same test on `cut_on_crick = cut_on_watson - ovhg`, `ovhg = CutPositionForward - CutPositionReverse`) expressed in 0-based coordinates. For an in-site cutter it excludes only an enzyme that cuts at the very edge of its own recognition window — in the built-in table `Sau3AI` and `MboI` (`^GATC` / `GATC^`) — sitting at the extreme 5' or 3' end; a distal cutter loses every site whose cut would fall beyond either molecule end (`FokI` needs 13 bases after `GGATG`), as in the reference.

On a **circular** molecule the recognition window may span the origin, so the search runs over `seq + seq[0 .. L-2]` (Biopython's `FormattedSeq.finditer` appends `data[1:size]`), cut positions are reduced into `[0, n)` with a non-negative modulo (a cut after the last base is the bond at the origin; a distal offset may be negative or exceed `n`), and nothing is dropped. Biopython's `_drop` folds a circular cut back by one length only and does not re-sort, so on a circle shorter than the cut offset it can report positions outside `[1, n]`; the reduced positions agree.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= Position <= sequence.Length - recognitionLength` for every reported site | Sites are yielded only from valid window starts |
| INV-02 | `RecognizedSequence.Length == Enzyme.RecognitionSequence.Length` | The source slices exactly the pattern length |
| INV-03 | Enzyme-name lookup is case-insensitive | The built-in enzyme dictionary uses `StringComparer.OrdinalIgnoreCase` |
| INV-04 | Empty raw-string input yields no sites | The raw-string overload short-circuits on null or empty input |
| INV-05 | `RecognizedSequence.Length == Enzyme.RecognitionLength` and the window matches the pattern under IUPAC rules on the strand named by `IsForwardStrand` | The source slices exactly the pattern length from the strand it matched |
| INV-06 | Linear molecules: `1 <= CutPosition <= Length - 1` for every reported site, on both strands | The cut-feasibility filter (§2.3) |
| INV-07 | Circular molecules: `0 <= CutPosition <= Length - 1` | Cut positions are reduced modulo the length |
| INV-08 | No site is reported when `Length < RecognitionLength` | A recognition window longer than the molecule cannot occur |
| INV-09 | `CreateMap`: `SitesByEnzyme.Keys` and `NonCutters` partition the considered enzyme set, both keyed by canonical table name | Both are derived from one case-insensitive, de-duplicated resolved set |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to scan | Null `DnaSequence` input throws; empty raw-string input yields no results |
| `enzymeName` | `string` | required | Name of the restriction enzyme to search | Unknown names throw `ArgumentException` |
| `enzymeNames` | `string[]` | required for multi-enzyme overload | Set of enzyme names to scan | Multi-enzyme search yields the union of per-enzyme results |
| `enzyme` | `RestrictionEnzyme` | required for custom-enzyme overload | Enzyme definition to apply directly | Null input throws `ArgumentNullException`; an empty or non-IUPAC `RecognitionSequence` throws `ArgumentException` eagerly |
| `topology` | `MoleculeTopology` | `Linear` | Whether the molecule is a linear fragment or a closed circle | See §2.3 |
| `recognitionPattern` | `string` | required for `FindPatternOccurrences` | Raw IUPAC pattern | Non-empty, IUPAC codes only (case-insensitive) |
| `startIndex` | `int` | `0` | First index `FindPatternOccurrences` considers | Must be non-negative |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Start position of the recognition sequence |
| `Enzyme` | `RestrictionEnzyme` | Enzyme definition that matched |
| `IsForwardStrand` | `bool` | `true` for forward-strand matches and `false` for reverse-strand matches |
| `CutPosition` | `int` | Computed cut position |
| `RecognizedSequence` | `string` | Actual sequence segment that matched the enzyme pattern |

### 3.3 Preconditions and Validation

The `DnaSequence` overloads require non-null sequence input. The raw-string overload returns no results for null or empty strings but still validates the enzyme name **eagerly**, at the call. Unknown enzyme names raise `ArgumentException`, while enzyme lookup itself is case-insensitive. A custom enzyme's recognition sequence is validated eagerly too: empty, or containing a character outside the 15 IUPAC nucleotide codes, raises `ArgumentException` (previously such a pattern either threw `ArgumentOutOfRangeException` from inside the IUPAC matcher or, when empty, reported one "site" per position). Pattern matching uses IUPAC ambiguity codes through the canonical `IupacHelper.MatchesIupac(...)`; ambiguity is interpreted in the **pattern** only, so a sequence position holding anything other than A/C/G/T blocks a match rather than widening it.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the enzyme definition from the built-in catalog or use the supplied custom enzyme.
2. Scan the forward strand for IUPAC-aware recognition-sequence matches.
3. Yield a forward-strand `RestrictionSite` for each match.
4. Reverse-complement the sequence and repeat the recognition search.
5. Convert reverse-strand positions back to forward coordinates and yield reverse-strand `RestrictionSite` records.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Recognition-sequence categories called out in the original document:

| Category | Examples | Recognition Length |
|----------|----------|-------------------|
| 4-cutters | AluI, MspI, TaqI | 4 bp |
| 6-cutters | EcoRI, BamHI, HindIII | 6 bp |
| 8-cutters | NotI, PacI, AscI | 8 bp |

Overhang types:

| Type | Description | Example |
|------|-------------|---------|
| 5' overhang | `CutPositionForward < CutPositionReverse` | EcoRI |
| 3' overhang | `CutPositionForward > CutPositionReverse` | PstI |
| Blunt end | `CutPositionForward == CutPositionReverse` | EcoRV |

Cognate methylation site (B08 D7 / F52). Each built-in carries REBASE `allenz.104` field `<6>`
verbatim as `RebaseMethylationField` (e.g. EcoRI `"3(6)"`, HpaII `"2(5)"`, BamHI `"5(4)"`, NcoI
`"?(4)"`, DpnI `""`) and its parse `CognateMethylationSites`, per the file's format header:
`X(Y)` where X is the 1-based base within the recognition sequence (negative = complementary strand,
numbered from that strand's 5' base; `?` = not determined) and Y is `(6)` N6-methyladenine, `(5)`
5-methylcytosine, `(4)` N4-methylcytosine; several entries are comma-separated. It is the site the
enzyme's **cognate** methyltransferase modifies to protect its own DNA — not the enzyme's sensitivity
to Dam/Dcm/CpG methylation, which is not modelled and is not read by site finding. The 39 values were
transcribed programmatically and are locked per enzyme (`docs/Evidence/RESTR-FIND-001-Evidence.md` §5).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Single-enzyme `FindSites` | `O(n × m)` | `O(1)` streaming | `n` is sequence length and `m` is recognition-sequence length |
| Multi-enzyme `FindSites` / `FindAllSites` | `O(n × k × m)` | `O(1)` streaming | `k` is the number of enzymes scanned |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs), [IupacHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/IupacHelper.cs)

- `RestrictionAnalyzer.FindSites(DnaSequence, string)`: Searches a validated DNA sequence for one enzyme.
- `RestrictionAnalyzer.FindSites(DnaSequence, params string[])`: Aggregates per-enzyme searches.
- `RestrictionAnalyzer.FindSites(string, string)`: Raw-string overload that uppercases input and yields no results for empty input.
- `RestrictionAnalyzer.FindAllSites(DnaSequence)`: Scans the 39 common enzymes (the default set).
- `RestrictionAnalyzer.GetEnzyme(string)`: Case-insensitive enzyme lookup against the **full REBASE set** (754 enzymes, §5.2a); every name-taking member resolves through it.
- `RestrictionAnalyzer.RebaseEnzymes` / `GetEnzymes(EnzymeSet)` and the `EnzymeSet` overloads of `GetEnzymesByCutLength`, `GetBluntCutters`, `GetStickyCutters`, `FindAllSites`, `FindCompatibleEnzymes(string, …)`, `FindIsoschizomers`/`FindEquischizomers`/`FindNeoschizomers`: choose the enumerated set explicitly (`EnzymeSet.Common`, the default, or `EnzymeSet.Rebase`). `CreateMap`/`GetDigestSummary`/`FindSites`/`Digest` over the full set: pass `RebaseEnzymes.Values` to their `IEnumerable<RestrictionEnzyme>` overloads.
- `RestrictionAnalyzer.FindSites(DnaSequence, string, MoleculeTopology)` / `FindSites(DnaSequence, RestrictionEnzyme, MoleculeTopology)` / `FindSites(string, string, MoleculeTopology)` / `FindSites(DnaSequence, MoleculeTopology, params string[])` / `FindAllSites(DnaSequence, MoleculeTopology)` / `CreateMap(DnaSequence, MoleculeTopology, params string[])`: the same searches on a linear fragment or a closed circle (§2.3).
- `RestrictionAnalyzer.FindPatternOccurrences(string sequence, string recognitionPattern, int startIndex = 0)`: the raw IUPAC pattern scan underlying the above — 0-based ascending start indices, overlapping occurrences included, one strand only, no enzyme/cut arithmetic and no feasibility filter. Exposed for callers that hold a pattern rather than an enzyme and a sequence `DnaSequence` would reject (e.g. one still carrying IUPAC ambiguity codes).

### 5.2 Current Behavior

The implementation ships with an in-memory catalog of 39 entries (37 distinct enzymes; `MspI`/`HpaII` share `CCGG` and `DpnI`/`Sau3AI`/`MboI` share `GATC`), including 4-cutters, 6-cutters, 8-cutters, blunt cutters, and enzymes with degenerate recognition motifs such as `SfiI` and `HincII`. Every entry's recognition sequence, both cut offsets and blunt/sticky classification were audited against REBASE EMBOSS version 404 (2024) via Biopython 1.88 on 2026-10-14: **39/39 matched on all four fields, 0 wrong entries** (test `EnzymeDatabase_AllEntries_MatchRebase404`). Both strands are reported, so palindromic sites appear twice at the same genomic position with different `IsForwardStrand` values; each record carries the cut made on its own strand and the window as read on that strand (§2.2). Overlapping sites are all reported, matching the reference (whose site regexes are lookaheads).

### 5.2a The vendored REBASE set (B08 D4 part 1, F59)

`RestrictionAnalyzer.RebaseEnzymes` holds every enzyme of the Biopython 1.88 REBASE 404 snapshot that has a
defined cut: **754 of 1088** (729 Biopython `OneCut` — 576 in-site + 153 distal — and 25 `TwoCuts`). It is
the embedded resource `Resources/rebase_enzymes.tsv`, generated by
`scripts/rebase/generate_rebase_enzymes.py` from two sha256-pinned inputs: Biopython's
`Restriction_Dictionary.py` (site, `CutPositionForward = fst5`, `CutPositionReverse = size + fst3`,
`SecondCut = (scd5, size + scd3)`) and REBASE `allenz.104` (organism field `<3>` verbatim with strain,
methylation field `<6>`), joined by name; 7 names absent from `allenz.104` have both empty. The 39 common
enzymes are represented by their curated instances, so nothing computed with them changes.

**Excluded — 334 enzymes with no defined cut** (Biopython `NoCut`/`NotDefined`, `fst5 is None`, e.g.
`Dpi3084I` CGRAG): a record must carry integer cut offsets, an invented offset would make every fragment a
false statement, and such an enzyme cannot be digested. `GetEnzyme` returns `null` for them.

**Which set a member uses.** Name lookups always resolve against the full set (`Digest(seq, "BsaI")` works).
Members that *enumerate* enzymes without being given a list — `FindAllSites`, `CreateMap` with no names,
`GetEnzymesByCutLength`, `GetBluntCutters`, `GetStickyCutters`, `FindCompatibleEnzymes`, the isoschizomer
queries — keep the 39 common enzymes by default (their output would otherwise grow ~20-fold), and take an
explicit `EnzymeSet` where the full set is wanted.

**Readings vs. windows (F60).** `FindSites` reports every strand reading of a window; Biopython `search`
reports one cut set per window (forward reading first). The two agree once a reverse record whose window
also reads as the site on the forward strand is set aside — the rule `Digest`/`CreateMap`/`PartialDigest`
apply. It is invisible for the common enzymes and matters for the 4 palindromic distal cutters (HauII,
SgeI, …) and for degenerate sites whose two readings can share a window (MspJI, LpnPI, SgrTI).

**Differential** (`RestrictionAnalyzer_RebaseSet_Tests`): every one of the 754 enzymes, resolved **by name**,
on a molecule with its site planted in both orientations, linear and circular: `FindSites` top-strand cuts
(with the rule above) = Biopython 1.88 `search`, and `Digest` fragment lengths = the gaps between the
sorted distinct `search` cuts — **0 mismatches** over 1 508 + 1 508 comparisons.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Recognition-sequence scanning for a catalog of Type II-style enzymes, audited entry by entry against REBASE 404.
- IUPAC ambiguity handling in enzyme patterns, through the canonical `IupacHelper`.
- Strand-aware reporting of recognition and cut positions, including the correct reverse-strand cut for a non-palindromic custom enzyme.
- Overlapping-site reporting.
- Linear-molecule cut feasibility (Biopython `AbstractCut._drop`).
- Circular-molecule (origin-spanning) search with modulo cut positions (Biopython `search(..., linear=False)`).
- The raw IUPAC pattern scan as a public API (`FindPatternOccurrences`).

**Intentionally simplified:**

- ~~The implementation relies on a built-in enzyme catalog; **consequence:** only enzymes present in that catalog are available through name-based lookup.~~ *(Superseded 2026-10-15, B08 D4 part 1 F59: name lookup covers all 754 REBASE 404 enzymes with a defined cut, §5.2a.)* The 334 REBASE enzymes without a defined cut position are not available — they cannot be digested. Any other enzyme can still be supplied through the `RestrictionEnzyme` overloads.
- The built-in catalog contains only enzymes whose cuts fall inside their site; **consequence:** distal cutters must be supplied as custom records. *(Corrected 2026-10-15, review B08 WP11a F54/F55: the earlier claim that Type IIS cutters "cannot be expressed" was wrong — the record's offsets are signed and unrestricted, so a single-sided distal cutter is `CutPositionForward = fst5`, `CutPositionReverse = size + fst3`, e.g. `BsaI` = `("GGTCTC", 7, 11)`, `FokI` = `("GGATG", 14, 18)`; see §2.2a.)* Two-sided Type IIB cutters (`BaeI`, `BcgI`) are expressible through `RestrictionEnzyme.SecondCut` (F56).
- Site finding is purely sequence-based; **consequence:** methylation sensitivity and experimental context are not modeled (`DpnI` requires Dam methylation; `HpaII`/`MspI` differ only in methylation sensitivity).

**Not implemented:**

- Methylation sensitivity (whether Dam/Dcm/CpG methylation overlapping a site blocks the enzyme). The
  record's `CognateMethylationSites` (REBASE `<6>`) is descriptive data only and must not be read as a
  blocking rule; no source reachable in this campaign tabulates sensitivity for the built-ins.
- Contextual enzyme activity effects such as methylation-state filtering; **users should rely on:** external enzymology references or downstream validation when those matter.

**Documented divergence from the oracle:**

- When the molecule is **shorter** than the recognition site and circular, Biopython still reports a site, because it pads a circular sequence with `data[1:size]` unconditionally and therefore wraps a short molecule more than once (`ClaI.search(Seq("GATC"), linear=False) == [4]` for the 6-bp site `ATCGAT` on a 4-bp circle). This library reports nothing: a 6-bp recognition site cannot be present on a 4-bp circle.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty raw-string input | Returns no sites | Explicit short-circuit in source |
| Unknown enzyme name | Throws `ArgumentException` | Name lookup must succeed before scanning |
| Palindromic recognition site | May appear on both strands at the same position | Both strands are searched independently |
| Degenerate recognition sequence | Matched with IUPAC ambiguity rules | `IupacHelper.MatchesIupac(...)` handles the pattern |
| IUPAC code in the *sequence* | Matched by a site `N` only; blocks every other site code | Biopython 1.88's rule (site `N` → regex `.`, other codes → their A/C/G/T class): `SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]`, `EcoRI.search(Seq("AGAANTCA")) == []`. Composed from `IupacHelper.IsNucleotideCode` / `MatchesIupac`. `DnaSequence` rejects such input, so this is reachable only via the raw-string and `FindPatternOccurrences` entry points (B08 F40) |
| Non-IUPAC character in the *sequence* | Matches nothing, not even a site `N` | Deliberate narrowing: Biopython lets `.` match any letter (`SfiI.search(Seq("GGCCXXXXXGGCC")) == [9]`) |
| `Sau3AI`/`MboI` site at the extreme 5' or 3' end of a linear molecule | Not reported | One of its two cuts would fall outside the molecule (§2.3) |
| Recognition site spanning the origin of a circular molecule | Reported, with cut positions modulo the length | `topology = Circular` (§2.3) |
| Molecule shorter than the recognition site | Returns no sites | The window cannot occur |
| Custom enzyme with an empty or non-IUPAC recognition sequence | Throws `ArgumentException` at the call | Eager validation |
| Overlapping occurrences of a self-overlapping site (e.g. `PacI` TTAATTAA in TTAATTAATTAA) | All reported | The scan advances one base at a time, like the reference's lookahead regex |

### 6.2 Limitations

The current implementation supports in-site, single-sided distal (Type IIS) and two-sided (Type IIB, `SecondCut`) cutters, and is limited to sequence-only matching. The 1088-enzyme Biopython 1.88 snapshot has 153 single-sided distal OneCut enzymes and 25 two-sided ones (Biopython `TwoCuts`), all expressible, 0 differential mismatches (F54, F56); the 334 enzymes without a defined cut remain out of scope. It does not model methylation dependence, assay conditions, or experimental digestion efficiency, and the two records reported per palindromic site must be interpreted downstream when unique cut counts are required (`RestrictionMap.TotalSites` and `UniqueCutters` count physical sites: forward-strand windows plus reverse-strand windows that do not also read forward — for the palindromic built-ins exactly the forward records; F58).

## 8. References

1. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
2. Wikipedia: Restriction site - https://en.wikipedia.org/wiki/Restriction_site
3. Wikipedia: EcoRI - https://en.wikipedia.org/wiki/EcoRI
4. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry.
5. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/ — used here as the **REBASE EMBOSS files version 404 (2024)** snapshot shipped in Biopython 1.88 `Bio/Restriction/Restriction_Dictionary.py` (that generated file states its own REBASE version); rebase.neb.com is not reachable from the validation environment.
6. IUPAC-IUB (1970) - Nucleic acid notation.
7. Biopython 1.88, `Bio.Restriction` (`Bio/Restriction/Restriction.py`) — the reference oracle for site positions, cut positions, cut feasibility on linear molecules and circular search. Relevant members: `FormattedSeq.finditer`, `AbstractCut._drop`, `Defined._modify`, `Defined._rev_modify`, `Palindromic._search`, `NonPalindromic._search`.
8. REBASE `allenz.104` (REBASE version 104, R. J. Roberts, 2021), field `<6>` METHYLATION SITE and its format header — as vendored in `crisprVerse/crisprBase` `inst/enzymes/allenz.txt` / `notation.txt` (sha256 pinned in `docs/Evidence/B08-Source-Provenance.md` §3a).
9. `scripts/rebase/generate_rebase_enzymes.py` — generator of the vendored table (`Resources/rebase_enzymes.tsv`) and the Biopython `search` oracle (`TestData/Restriction/biopython_rebase_search.tsv`); inputs pinned by sha256.
