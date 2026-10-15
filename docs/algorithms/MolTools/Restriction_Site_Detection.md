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

### 2.3 Cut Feasibility and Topology

A cut can only fall on a phosphodiester bond that exists. On a **linear** molecule of length `n` a site is therefore reported only when

$$
1 \le watson \le n - 1 \quad\text{and}\quad 1 \le crick \le n - 1 ,
$$

which is Biopython's `AbstractCut._drop` filter (`1 < cut_on_watson <= length` together with the same test on `cut_on_crick = cut_on_watson - ovhg`, `ovhg = CutPositionForward - CutPositionReverse`) expressed in 0-based coordinates. It excludes only an enzyme that cuts at the very edge of its own recognition window — in the built-in table `Sau3AI` and `MboI` (`^GATC` / `GATC^`) — sitting at the extreme 5' or 3' end.

On a **circular** molecule the recognition window may span the origin, so the search runs over `seq + seq[0 .. L-2]` (Biopython's `FormattedSeq.finditer` appends `data[1:size]`), cut positions are reduced modulo `n` (a cut after the last base is the bond at the origin), and nothing is dropped.

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
- `RestrictionAnalyzer.FindAllSites(DnaSequence)`: Scans the full built-in enzyme set.
- `RestrictionAnalyzer.GetEnzyme(string)`: Case-insensitive enzyme lookup.
- `RestrictionAnalyzer.FindSites(DnaSequence, string, MoleculeTopology)` / `FindSites(DnaSequence, RestrictionEnzyme, MoleculeTopology)` / `FindSites(string, string, MoleculeTopology)` / `FindSites(DnaSequence, MoleculeTopology, params string[])` / `FindAllSites(DnaSequence, MoleculeTopology)` / `CreateMap(DnaSequence, MoleculeTopology, params string[])`: the same searches on a linear fragment or a closed circle (§2.3).
- `RestrictionAnalyzer.FindPatternOccurrences(string sequence, string recognitionPattern, int startIndex = 0)`: the raw IUPAC pattern scan underlying the above — 0-based ascending start indices, overlapping occurrences included, one strand only, no enzyme/cut arithmetic and no feasibility filter. Exposed for callers that hold a pattern rather than an enzyme and a sequence `DnaSequence` would reject (e.g. one still carrying IUPAC ambiguity codes).

### 5.2 Current Behavior

The implementation ships with an in-memory catalog of 39 entries (37 distinct enzymes; `MspI`/`HpaII` share `CCGG` and `DpnI`/`Sau3AI`/`MboI` share `GATC`), including 4-cutters, 6-cutters, 8-cutters, blunt cutters, and enzymes with degenerate recognition motifs such as `SfiI` and `HincII`. Every entry's recognition sequence, both cut offsets and blunt/sticky classification were audited against REBASE EMBOSS version 404 (2024) via Biopython 1.88 on 2026-10-14: **39/39 matched on all four fields, 0 wrong entries** (test `EnzymeDatabase_AllEntries_MatchRebase404`). Both strands are reported, so palindromic sites appear twice at the same genomic position with different `IsForwardStrand` values; each record carries the cut made on its own strand and the window as read on that strand (§2.2). Overlapping sites are all reported, matching the reference (whose site regexes are lookaheads).

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

- The implementation relies on a built-in enzyme catalog; **consequence:** only enzymes present in that catalog are available through name-based lookup. Any other enzyme can still be supplied through the `RestrictionEnzyme` overloads.
- A `RestrictionEnzyme` expresses one cut per strand **inside** the recognition window, so Type IIS / Type IIG enzymes that cut at a distance (e.g. `BsaI` GGTCTC(1/5), `FokI`) cannot be expressed; **consequence:** the catalog contains only enzymes whose cuts fall inside their site. Adding the reference's `fst5`/`fst3` sign convention (which allows cuts outside the site) would be an additive record change.
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

The current implementation is limited to enzymes whose cuts fall inside their recognition window (so Type IIS / Type IIG distal cutters cannot be expressed as a `RestrictionEnzyme`) and to sequence-only matching. It does not model methylation dependence, assay conditions, or experimental digestion efficiency, and the two records reported per palindromic site must be interpreted downstream when unique cut counts are required (`RestrictionMap.TotalSites` and `UniqueCutters` already count forward-strand records only).

## 8. References

1. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
2. Wikipedia: Restriction site - https://en.wikipedia.org/wiki/Restriction_site
3. Wikipedia: EcoRI - https://en.wikipedia.org/wiki/EcoRI
4. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry.
5. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/ — used here as the **REBASE EMBOSS files version 404 (2024)** snapshot shipped in Biopython 1.88 `Bio/Restriction/Restriction_Dictionary.py` (that generated file states its own REBASE version); rebase.neb.com is not reachable from the validation environment.
6. IUPAC-IUB (1970) - Nucleic acid notation.
7. Biopython 1.88, `Bio.Restriction` (`Bio/Restriction/Restriction.py`) — the reference oracle for site positions, cut positions, cut feasibility on linear molecules and circular search. Relevant members: `FormattedSeq.finditer`, `AbstractCut._drop`, `Defined._modify`, `Defined._rev_modify`, `Palindromic._search`, `NonPalindromic._search`.
8. REBASE `allenz.104` (REBASE version 104, R. J. Roberts, 2021), field `<6>` METHYLATION SITE and its format header — as vendored in `crisprVerse/crisprBase` `inst/enzymes/allenz.txt` / `notation.txt` (sha256 pinned in `docs/Evidence/B08-Source-Provenance.md` §3a).
