# Restriction Enzyme Filtering

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | RESTR-FILTER-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-10-21 (review-2026-09 batch B08, unit RESTR-FILTER-001) |

## 1. Overview

Restriction enzyme filtering selects subsets of a restriction-enzyme library by the properties that matter for cloning design: the length of the recognition sequence, the type of DNA end the enzyme produces (blunt, 5' overhang, 3' overhang), the **overhang sequence** it leaves, which other enzymes leave a **ligatable (compatible)** end, and which other enzymes recognise the **same site** (isoschizomers). These are pure, specification-driven set operations over a curated enzyme table — no sequence input is involved. They let a user, for example, list all 6-cutters, pick blunt-end producers for a particular ligation strategy, or ask what an SfiI end can be ligated to.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Type II restriction endonucleases recognize short, usually palindromic DNA sequences and cleave at or near them. Their undivided recognition sites are "4–8 nucleotides in length" [2] (a few enzymes, e.g. SfiI, use a longer interrupted palindrome with a degenerate spacer [6]). Cleavage either happens "at the center of both strands to yield a blunt end, or at a staggered position leaving overhangs called sticky ends" [2]. A blunt-ended molecule is one where "both strands terminate in a base pair" [1]; a sticky (cohesive) end carries "a stretch of unpaired nucleotides" — a 5' or 3' overhang [1].

### 2.2 Core Model

For an enzyme `e` with recognition sequence `r(e)` and per-strand cut positions `cf(e)` (forward) and `cr(e)` (reverse, measured from the 5' end of the recognition site):

- **Recognition length:** `len(e) = |r(e)|`.
- **End type:** `e` is **blunt** iff `cf(e) = cr(e)` (the two strands are cut at the same offset, i.e. a center/symmetric cut [2]); otherwise `e` is **sticky** (a staggered cut leaving a 5' or 3' overhang [1][2]).
- **Length filter (range):** `ByLength(min, max) = { e : min ≤ len(e) ≤ max }` (both bounds inclusive).
- **Blunt set:** `Blunt = { e : cf(e) = cr(e) }`; **Sticky set:** `Sticky = { e : cf(e) ≠ cr(e) }`.
- **Polarity:** `pol(e) = Blunt` if `cf = cr`, `5'` if `cf < cr`, `3'` if `cf > cr` — the sign of
  REBASE/Biopython `ovhg` [7].
- **Overhang window / sequence:** the bases lying between the two strand cuts,
  `win(e) = [min(cf, cr), max(cf, cr))`; the enzyme-level overhang sequence is that window read off
  the recognition pattern, `W(e) = r(e)[min(cf,cr) .. max(cf,cr))`, written 5'→3' on the top strand.
  `|W(e)| = |ovhg(e)|`, and `W(e) = ""` exactly when `e` is blunt. This is precisely Biopython's
  `ovhgseq` [7] — including the degenerate case: where the site carries an IUPAC code inside the
  window, `W` carries that code (`SfiI` → `NNN`), because the enzyme-level overhang of a divided
  palindrome genuinely is unspecified bases. The *per-site* bases are a property of a molecule, not
  of the enzyme, and are reported by the digest instead (`DigestFragment.LeftOverhang` /
  `RightOverhang`).
- **Ligation compatibility.** A cut leaves two ends: the downstream fragment's extension reads `W`
  on the protruding strand, the upstream fragment's reads `revcomp(W)` on the other strand. They are
  complementary and antiparallel, so they re-ligate — an enzyme therefore offers **both**. Two ends
  anneal iff their extensions have equal length, **equal polarity** and one is the reverse
  complement of the other; opposite polarity can never anneal, because both extensions would lie on
  the same strand of the joint. Hence

  `Compatible(e1, e2) ⟺ (pol(e1) = pol(e2) = Blunt) ∨ (pol(e1) = pol(e2) ≠ Blunt ∧ W(e1) ≈ W(e2)) ∨ (pol(e1) = pol(e2) ≠ Blunt ∧ W(e1) ≈ revcomp(W(e2)))`

  where `a ≈ b` means equal length and, at every position, the two IUPAC codes admit a common base
  (an ambiguity code denotes a set, independently per position, so `expand(W)` is a Cartesian
  product and the set condition `expand(W1) ∩ expand(W2) ≠ ∅` reduces to per-position intersection).
  A verdict involving a degenerate overhang is therefore "**possibly** compatible — some pair of
  real sites ligates".
- **Isoschizomers.** `Iso(e) = { f ≠ e : r(f) = r(e) }`, split into
  `Equi(e) = { f ∈ Iso(e) : cf(f) = cf(e) ∧ cr(f) = cr(e) }` (same site, same cut — interchangeable
  in a digest) and `Neo(e) = Iso(e) \ Equi(e)` (same site, different cut — different ends). These are
  Biopython's `equischizomers` / `neoschizomers` / `isoschizomers` [7].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Blunt ∩ Sticky = ∅` and `Blunt ∪ Sticky = Library` (total partition). | Every DNA end is either blunt or an overhang [1]; the predicate `cf = cr` is a Boolean split of the library. |
| INV-02 | `GetBluntCutters()` returns exactly enzymes with `cf = cr`. | Center cut → blunt [2]; implemented via `IsBluntEnd`. |
| INV-03 | `ByLength(min,max)` returns exactly enzymes with `min ≤ len ≤ max`. | Direct inclusive-range predicate; recognition lengths are 4–8 nt [2][3]. |
| INV-04 | `ByLength(L) = ByLength(L,L)` for any L. | Single-length filter is the degenerate equal-bounds range. |
| INV-05 | `min > max ⇒ ByLength(min,max) = ∅`. | No integer length satisfies an empty interval. |
| INV-06 | `|W(e)| = |cf(e) − cr(e)|`, and `W(e) = ""` iff `e ∈ Blunt`. | `W` is the window between the two cuts; a center cut leaves none [7]. |
| INV-07 | `Compatible` is symmetric, and reflexive for every enzyme whose overhang is determined (all 39 built-ins). | `≈` is symmetric by construction, and `W ≈ W` trivially holds; an enzyme that cuts outside its own site has the `N`-padded window `W` (F57), so it is only "possibly compatible" under §2.2. |
| INV-08 | `pol(e1) ≠ pol(e2) ⇒ ¬Compatible(e1, e2)`, same overhang letters notwithstanding. | Two extensions of opposite polarity sit on the same strand of the joint and cannot base-pair [1][7]. |
| INV-09 | `Iso(e) = Equi(e) ⊎ Neo(e)` and `e ∉ Iso(e)`. | By construction; matches the reference's own `isoschizomers = equi ∪ neo`, `del r[i]` [7]. |
| INV-10 | `f ∈ Equi(e) ⇒ W(f) = W(e) ⇒ Compatible(e, f)`. | Same site and same cut offsets give the same window. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `length` | `int` | required | Exact recognition-sequence length to match (bp) | any int; meaningful 4–8 |
| `minLength` | `int` | required | Inclusive lower bound on recognition length (bp) | any int |
| `maxLength` | `int` | required | Inclusive upper bound on recognition length (bp) | any int |

`GetBluntCutters()`, `GetStickyCutters()` and `FindCompatibleEnzymes()` take no parameters.

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `enzymeName`, `enzyme1Name`, `enzyme2Name` | `string` | required | Built-in enzyme name, resolved case-insensitively | must be in the table, else `ArgumentException`; `null`/empty → `ArgumentNullException` |
| `enzyme`, `enzyme1`, `enzyme2` | `RestrictionEnzyme` | required | Enzyme definition; may be caller-supplied, i.e. not a built-in | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IEnumerable<RestrictionEnzyme>` | The matching enzymes from the built-in library (deferred LINQ enumeration; order follows dictionary insertion). |

### 3.3 Preconditions and Validation

The four **set filters** (`GetEnzymesByCutLength` ×2, `GetBluntCutters`, `GetStickyCutters`) are
total: they never throw and never return null. They enumerate the static built-in enzyme dictionary
and apply a predicate. A range with `minLength > maxLength`, or non-positive bounds, yields an empty
sequence rather than an error. There is no sequence input, so no alphabet/normalization concerns
apply.

The **name-taking** members (`GetEnzyme`, `GetOverhangSequence(string)`,
`AreCompatible(string, string)`, `FindCompatibleEnzymes(string)`, `FindIsoschizomers` /
`FindEquischizomers` / `FindNeoschizomers`) share one contract with the rest of the class
(`FindSites`, `Digest`, `CreateMap`): a `null` or empty name is an `ArgumentNullException` and a name
outside the table is an `ArgumentException("Unknown enzyme: …")`. The one exception is `GetEnzyme`
itself, which is the lookup primitive and returns `null` for an unknown name (but still throws
`ArgumentNullException(nameof(name))` for `null`).

**Behaviour change, 2026-10-14 (B08/F31):** `AreCompatible(string, string)` previously returned
`false` for an unknown name. A misspelling was therefore indistinguishable from a genuinely
incompatible pair — and "not compatible" is the answer that silently loses a real cloning option. It
now throws. `TryAreCompatible(string, string, out bool)` was added for callers that cannot guarantee
the names exist: it returns `false` when either name is unknown and leaves the verdict `false`.

`AreCompatible(RestrictionEnzyme, RestrictionEnzyme)` and
`GetOverhangSequence(RestrictionEnzyme)` accept caller-supplied enzymes. If an enzyme cuts **outside**
its own recognition sequence (a Type IIS distal cut, or a two-sided Type IIB cutter), its overhang
window is read off the recognition pattern **padded with `N`** — the bases outside the site belong to
the flanking molecule: BsaI GGTCTC(1/5) → `NNNN`, BsmI → `CN`, TspRI → `NNCASTGNN`, BaeI → `NNNNN`.
That is Biopython 1.88's own `ovhgseq` for all 178 distal / two-cut enzymes of its snapshot, which it
classes `Ambiguous` (not `NotDefined`), and the ligation rule above then applies unchanged: an `N`
intersects every base, so BsaI is "possibly compatible" with every 4-nt 5'-overhang end
(`BsaI % EcoRI` is True) and with nothing 3' or blunt. Differential over all 236 740 ordered pairs that
involve one of the 178 (× all 754 expressible enzymes) against the symmetrised `%`: 0 disagreements
apart from the 24 reverse-complement-only pairs of F30 (review 2026-09 B08 F57; before: 28 072
disagreements, every distal end judged incompatible and `GetOverhangSequence` `null`).

## 4. Algorithm

### 4.1 High-Level Steps

1. Enumerate the static enzyme library values.
2. Apply the predicate: recognition-length equality, inclusive-range membership, or `IsBluntEnd` / `!IsBluntEnd`.
3. Yield the matching enzymes (lazy LINQ `Where`).

For the overhang and compatibility members:

4. `OverhangWindow(e)` returns the window as `(start, length) = (min(cf, cr), |cf − cr|)` — the
   **single** definition of that geometry in the class. `GetOverhangSequence` reads it off the
   recognition pattern; the digest's `CutOverhangWindow` reads the same window off the molecule,
   translated so the enzyme's forward cut lands on the site's reported cut position. There is no
   second derivation.
5. `AreCompatible` applies the rule of §2.2: blunt×blunt → true; different polarity → false;
   otherwise `OverhangsAnneal(W1, W2)`, which tests equal length and then per-position IUPAC
   intersection against `W2` and against `revcomp(W2)` (via the canonical
   `DnaSequence.GetReverseComplementString` and `IupacHelper.MatchesIupac` — no second ambiguity
   table lives in this class).
6. `FindCompatibleEnzymes()` enumerates unordered pairs of distinct enzymes;
   `FindCompatibleEnzymes(name)` enumerates the whole table against one enzyme and **includes that
   enzyme itself**, so an enzyme with a unique overhang is not invisible.
7. The three isoschizomer queries share one private same-recognition-sequence scan and then split on
   cut-position equality.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The library is a fixed `Dictionary<string, RestrictionEnzyme>` of Type II enzymes; each `RestrictionEnzyme` stores its recognition sequence and per-strand cut positions, from which `RecognitionLength` and `IsBluntEnd` (`cf == cr`) are derived. The blunt/sticky predicate uses only the cut-position equality, matching the center-vs-staggered criterion [2]. End-type classification of representative enzymes (SmaI/EcoRV/AluI/HaeIII = blunt; EcoRI = 5' overhang; KpnI/PstI = 3' overhang) is corroborated against [2], [4], [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Any filter | O(e) | O(1) extra | `e` = library size; single linear pass, deferred enumeration. No search over sequences, so the repository suffix tree is not applicable. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs)

- `RestrictionAnalyzer.GetEnzyme(string name)`: lookup primitive; `null` for an unknown name.
- `RestrictionAnalyzer.Enzymes`: the whole table as a read-only dictionary.
- `RestrictionAnalyzer.GetEnzymesByCutLength(int length)`: enzymes whose recognition length equals `length`.
- `RestrictionAnalyzer.GetEnzymesByCutLength(int minLength, int maxLength)`: enzymes whose recognition length is in the inclusive range.
- `RestrictionAnalyzer.GetBluntCutters()`: enzymes with `IsBluntEnd == true`.
- `RestrictionAnalyzer.GetStickyCutters()`: enzymes with `IsBluntEnd == false`.
- `RestrictionAnalyzer.GetOverhangSequence(string)` / `(RestrictionEnzyme)`: the enzyme-level
  overhang, REBASE/Biopython `ovhgseq`.
- `RestrictionAnalyzer.AreCompatible(string, string)` / `(RestrictionEnzyme, RestrictionEnzyme)` and
  `TryAreCompatible(string, string, out bool)`: the ligation rule.
- `RestrictionAnalyzer.FindCompatibleEnzymes()` / `(string)`: all compatible pairs / one enzyme's
  compatible set (itself included).
- `RestrictionAnalyzer.FindIsoschizomers(string)` / `FindEquischizomers(string)` /
  `FindNeoschizomers(string)`: same-site enzymes, split by cut position.

**MCP surface (Seqeron.Mcp.MolTools):** `get_enzyme` (now also reports `overhangSequence`,
`equischizomers`, `neoschizomers`), `enzymes_by_cut_length`, `blunt_cutters`, `sticky_cutters`,
`compatible_enzymes` (now takes an optional `enzyme_name` for the per-enzyme query),
`enzymes_compatible` (now also reports `compatibleEnd`, and an unknown name is an error). No tool
was added or removed by this unit. Since review B08 WP12b (F61) the listing tools take an optional
`enzyme_set` (`"common"` default = the 39, byte-identical to before; `"rebase"` = all 754 through
the `EnzymeSet` overloads): `enzymes_by_cut_length`, `blunt_cutters`, `sticky_cutters`,
`get_enzyme`'s schizomer lists and `compatible_enzymes` (per-enzyme form and, since F72, the paged
all-pairs form with `max_pairs` / `offset`); an unknown value is an `ArgumentException` on
`enzyme_set`. Names are always resolved against all 754.

### 5.2 Current Behavior

Filters operate over the enzyme library only (never over a user sequence): the parameterless forms read the 39-enzyme common set, and every filter and listing member has an `EnzymeSet` overload (`GetEnzymesByCutLength` ×2, `GetBluntCutters`, `GetStickyCutters`, `FindIsoschizomers` / `FindEquischizomers` / `FindNeoschizomers`, `FindCompatibleEnzymes(string, EnzymeSet)` and `FindCompatibleEnzymes(EnzymeSet)`) that switches them to `EnzymeSet.Rebase` — all 754 REBASE 404 enzymes with a defined cut (B08 F61, F72). `IsBluntEnd` is a record-derived property: `CutPositionForward == CutPositionReverse`. The methods return deferred `IEnumerable<>`; callers materialize with `ToList()` when a snapshot is needed. The unit is not a search/matching operation, so the repository suffix tree was evaluated and found inapplicable (no text to search).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Blunt vs sticky split by center (`cf = cr`) vs staggered cut [2]; blunt = "both strands terminate in a base pair", sticky = an overhang [1].
- Inclusive recognition-length range filter over Type II 4–8 nt sites [2][3].
- Total partition of the library into blunt and sticky (INV-01) [1].

- The enzyme-level overhang sequence as REBASE/Biopython `ovhgseq`, degenerate codes included [7];
  audited base for base against that oracle for all 39 built-ins (39/39).
- The full ligation rule of §2.2, including the polarity restriction (INV-08), the
  reverse-complement branch and per-position IUPAC intersection; cross-checked against Biopython's
  `%` operator on 1 521 built-in pairs and 75 466 REBASE pairs with a fully-defined overhang
  (0 disagreements each).
- Isoschizomer / equischizomer / neoschizomer reporting [7]; 117/117 queries identical to the
  reference restricted to the same 39 enzymes.

**Intentionally simplified:**

- (none)

**Not implemented:**

- Nothing in this unit's scope. Filtering *directly* by overhang direction or overhang sequence has
  no dedicated method, but it is one `Where` over the public `RestrictionEnzyme.OverhangType` /
  `GetOverhangSequence`, and the compatibility members answer the question these filters would be
  used for.
- ~~Out of this unit's scope by construction: the table is a curated 39-enzyme subset, not REBASE; and a `RestrictionEnzyme` carries one cut offset per strand *relative to the site*, so an enzyme that cuts distally (Type IIS `BsaI` GGTCTC(1/5), `FokI`) has no enzyme-level overhang — `GetOverhangSequence` returns `null` and `AreCompatible` returns `false` rather than guessing. Supporting such enzymes means adopting REBASE's signed `fst5`/`fst3` convention on the record.~~ *(Withdrawn, review B08 F54-F59/F61: the record's offsets are signed and unrestricted, so a distal cutter is expressible — `BsaI` = `("GGTCTC", 7, 11)`, two-sided Type IIB cutters through `RestrictionEnzyme.SecondCut`; `EnzymeSet.Rebase` reaches all 754 REBASE enzymes with a defined cut, not only the 39; and since F57 `GetOverhangSequence` of a distal cutter is the `N`-padded window (`BsaI` → `NNNN`) with the "possibly compatible" rule of §2.2, as in Biopython. Same strike as `Restriction_Site_Detection.md` §5.3.)*  Still out of scope: the 334 REBASE enzymes without a defined cut position, which cannot be digested and are not in either set.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Range bounds inclusive | Assumption | Determines whether boundary lengths are returned | accepted | API-shape convention; recognition-length values are source-backed [2][3]. |
| 2 | Recognition length is the full site string, so `SfiI` is 13 and not 8 | Convention | `GetEnzymesByCutLength(8)` does not return `SfiI` | accepted | Matches Biopython `.size == 13` [7]; documented on the method and in the MCP tool description, because catalogue prose calls SfiI an "8-base-specificity" enzyme. |
| 3 | A degenerate overhang yields a "possibly compatible" verdict | Derivation | `SfiI` is reported compatible with any 3-nt 3' overhang | accepted | Follows from the IUPAC set semantics (§2.2); it is also the reference's intent (`Ambiguous._mod2` expands the code to a regex) [7]. |
| 4 | Biopython's `%` operator is a cross-check, not the specification | Deviation | Two `%` results are deliberately not reproduced | accepted | `%` dispatches on the left operand and only `Defined._mod2` falls back, so it is **asymmetric for 2 386 of 283 881 REBASE pairs**; and `Ambiguous._mod2` builds a broken regex when an ambiguity code repeats in the overhang (`StyI % PscI`, `"CWWG"` vs `"CATG"`, is `False` in both directions). The rule of §2.2 is symmetric by construction and agrees with `%` wherever `%` is sound. Both defects were reproduced this session; see the Evidence addendum §A4. |
| 5 | `AreCompatible(string, string)` throws for an unknown name | Behaviour change (2026-10-14) | A caller relying on the old `false` now sees an exception | accepted | Consistency with every other name-taking member; `TryAreCompatible` preserves the lenient path. Only in-repo caller was the MCP `enzymes_compatible` tool, updated with it. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `ByLength(min,max)` with `min > max` | empty | Empty interval (INV-05) |
| `ByLength(4,8)` | full library except SfiI | Undivided sites are 4–8 nt [2]; SfiI is a 13-nt interrupted palindrome `GGCCNNNN^NGGCC` [6], excluded |
| `ByLength(9,10)` | empty | Above the 8-nt maximum for undivided sites [2] |
| Non-positive bounds | empty | No recognition site has length ≤ 0 |
| `GetOverhangSequence` of a blunt cutter | `""` | `ovhg == 0`; there is no window (INV-06) |
| `GetOverhangSequence("SfiI")` | `"NNN"` | The enzyme-level overhang of a divided palindrome is unspecified bases; Biopython reports the same [7] |
| `GetOverhangSequence` of a distal (Type IIS-style) cutter | the `N`-padded window (`BsaI` GGTCTC(1/5) → `"NNNN"`); never `null` (F57) | The bases outside the site belong to the flanking molecule; Biopython's `ovhgseq` for all 178 distal / two-cut enzymes [7] |
| `AreCompatible` with a distal cutter | `true` against every end whose polarity and length fit the `N` window (BsaI vs EcoRI: "possibly compatible"); `false` against 3' and blunt ends | `N` intersects every base (§2.2, F57); `BsaI % EcoRI` is True in Biopython [7] |
| 5' and 3' overhangs with the same letters (HindIII/SacI, NcoI/SphI, NotI/ApaI, AgeI/FseI) | `false` | INV-08 — opposite polarity cannot anneal |
| `AreCompatible` with an unknown name | `ArgumentException` | F31; use `TryAreCompatible` to probe |
| `FindCompatibleEnzymes()` and an enzyme with a unique overhang | that enzyme appears in **no** pair | Only pairs of *distinct* enzymes are enumerated; 15 of the 39 are in this position — use `FindCompatibleEnzymes(name)` |

### 6.2 Limitations

The default (parameterless) library is the curated 39-enzyme common set, so its filters and compatible-pair lists reflect only those enzymes; the `EnzymeSet.Rebase` overloads cover all 754 REBASE enzymes with a defined cut. On the common set `FindCompatibleEnzymes("EcoRI")` returns EcoRI alone because REBASE's `MfeI` (`CAATTG`, 5' `AATT`) is not one of the 39; ~~"`EcoRI` is reported compatible with itself alone"~~ is **no longer true of the REBASE set** — `FindCompatibleEnzymes("EcoRI", EnzymeSet.Rebase)` (MCP `compatible_enzymes("EcoRI", "rebase")`) returns `MfeI` among the 51 enzymes of Biopython's `EcoRI.compatible_end(...)` (MCP test `Rebase_CompatibleEnzymes_SingleEnzyme_MatchesBiopythonCompatibleEnd`), and the all-pairs listing `FindCompatibleEnzymes(EnzymeSet.Rebase)` contains `("EcoRI", "MfeI", "AATT")` (`FindCompatibleEnzymesSet_Rebase_MatchesBiopythonModOperatorOnTheDefinedOverhangEnzymes`, F42/F63/F72); `AreCompatible(EcoRI, MfeI)` is asserted in `AreCompatible_EcoRIAndMfeI_ShareTheFivePrimeAattEnd`. Blunt/sticky classification and the overhang sequence depend on the stored cut positions being correct; they do not re-derive cleavage from sequence. The stored values were audited entry by entry against Biopython 1.88 on 2026-10-14 (39 enzymes × 6 fields, **0 wrong entries**). An enzyme that cuts outside its own site has an `N`-padded enzyme-level overhang (§3.3, F57) rather than none.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// All blunt-end producers in the library
var blunt = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToList();

// All 6-cutters
var sixCutters = RestrictionAnalyzer.GetEnzymesByCutLength(6, 6).ToList();

// Every undivided 4-8 nt site (excludes the 13-nt interrupted palindrome SfiI)
var undivided = RestrictionAnalyzer.GetEnzymesByCutLength(4, 8).ToList();

// The enzyme-level overhang (REBASE / Biopython ovhgseq)
RestrictionAnalyzer.GetOverhangSequence("EcoRI");   // "AATT"  (5', ovhg -4)
RestrictionAnalyzer.GetOverhangSequence("PstI");    // "TGCA"  (3', ovhg +4)
RestrictionAnalyzer.GetOverhangSequence("SfiI");    // "NNN"   (3', ovhg +3)
RestrictionAnalyzer.GetOverhangSequence("EcoRV");   // ""      (blunt)

// End compatibility
RestrictionAnalyzer.AreCompatible("BamHI", "BglII");   // true  — both 5' GATC
RestrictionAnalyzer.AreCompatible("SalI", "XhoI");     // true  — both 5' TCGA
RestrictionAnalyzer.AreCompatible("HindIII", "SacI");  // FALSE — 5' AGCT vs 3' AGCT
RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI");    // true  — blunt + blunt
RestrictionAnalyzer.AreCompatible("EcoRI", "Nope");    // throws ArgumentException
RestrictionAnalyzer.TryAreCompatible("EcoRI", "Nope", out bool c); // false, c == false

// What can an SfiI end be ligated to? (it is in no pair of FindCompatibleEnzymes())
RestrictionAnalyzer.FindCompatibleEnzymes("SfiI").Select(t => t.Enzyme);   // { "SfiI" }
RestrictionAnalyzer.FindCompatibleEnzymes("BamHI").Select(t => t.Enzyme);  // BamHI, BglII, MboI, Sau3AI

// Same site, interchangeable or not
RestrictionAnalyzer.FindEquischizomers("MspI").Select(e => e.Name);  // { "HpaII" }  (CCGG 1/3 both)
RestrictionAnalyzer.FindNeoschizomers("MboI").Select(e => e.Name);   // { "DpnI" }   (GATC, 0/4 vs 2/2)
```

### 7.2 Reference cross-check (2026-10-14)

| Comparison against Biopython 1.88 `Bio.Restriction` | Count | Disagreements |
|---|---|---|
| 39 built-ins × {site, `size`, `fst5`, `size+fst3`, `is_blunt`, polarity, `ovhgseq`} | 39 × 7 | **0** |
| `AreCompatible` vs `a % b or b % a`, all ordered built-in pairs | 1 521 | **0** |
| Same rule vs `a % b` over the 389 REBASE enzymes with a fully-defined overhang | 75 466 | **0** |
| Isoschizomer queries vs `*(batch = the 39)` | 117 | **0** |

Class counts, identical on both sides: 10 blunt, 21 with a 5' overhang, 8 with a 3' overhang.
`FindCompatibleEnzymes()` yields **64** pairs: 45 blunt×blunt, 6 `GATC`, 6 `CTAG`, 6 `CG`, 1 `TCGA`.

**`FindCompatibleEnzymes(EnzymeSet)` (2026-10-21, B08 F72).** The all-pairs listing takes an `EnzymeSet` like its
siblings; `Common` reproduces `FindCompatibleEnzymes()` exactly (same 64 pairs, same order) and `Rebase` yields
**27 869** of the 283 881 unordered pairs of the 754 enzymes, lazily, in table order. Over the 389 enzymes with a
blunt or fully A/C/G/T overhang the listing is 11 200 pairs (9 316 blunt) and equals Biopython 1.88 `a % b`
exactly (sha256 `a38477bb…d030` over the sorted `name1,name2,overhang` lines); `(EcoRI, MfeI, AATT)` is among
them. The other 16 669 pairs involve an ambiguous (distal Type IIS) overhang and follow the derived rule of §2.2.
The MCP `compatible_enzymes` tool pages this listing (`max_pairs`, `offset`; `totalPairs`, `truncated`).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RestrictionAnalyzer_Filter_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Filter_Tests.cs) — covers `INV-01`–`INV-10`
- Tests (pre-existing compatibility cases, from RESTR-DIGEST-001's spec): [RestrictionAnalyzer_Digest_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Digest_Tests.cs)
- Evidence: [RESTR-FILTER-001-Evidence.md](../../../docs/Evidence/RESTR-FILTER-001-Evidence.md)
- Related algorithms: [Restriction_Digest_Simulation](../MolTools/Restriction_Digest_Simulation.md)

## 8. References

1. Wikipedia. 2026. Sticky and blunt ends. https://en.wikipedia.org/wiki/Sticky_and_blunt_ends
2. Wikipedia. 2026. Restriction enzyme. https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia. 2026. List of restriction enzyme cutting sites. https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites
4. New England Biolabs / REBASE. KpnI (R0142). https://www.neb.com/en/products/r0142-kpni
5. New England Biolabs / REBASE. EcoRI-HF (R3101). https://www.neb.com/en/products/r3101-ecori-hf
6. SfiI interrupted palindrome `5'-GGCCNNNN^NGGCC-3'` (homology model article). PMC. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC548270/
7. Biopython 1.88, `Bio.Restriction` (`Restriction.py`) — reference implementation, read locally:
   per-enzyme `.site`, `.size`, `.fst5`, `.fst3`, `.ovhg`, `.ovhgseq`, `.elucidate()`,
   `.is_blunt()`, `.is_5overhang()`, `.is_3overhang()`, `.isoschizomers/equischizomers/neoschizomers(batch)`,
   and the `%` compatible-ends operator (`RestrictionType.__mod__` → `Blunt/Ov5/Ov3._mod1` →
   `Defined/Ambiguous/NotDefined._mod2`). **References 1–6 were NOT re-opened on 2026-10-14**:
   `rebase.neb.com`, `en.wikipedia.org` and `www.addgene.org` are blocked from this container (the
   agent proxy answers 403 to CONNECT). See `docs/Evidence/RESTR-FILTER-001-Evidence.md` §A1.
