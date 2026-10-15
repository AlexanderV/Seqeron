# Restriction Digest Simulation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | N/A |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-10-15 |

## 1. Overview

Restriction digest simulation models DNA fragmentation after cleavage by one or more restriction enzymes. In this repository, digestion is performed by locating cut sites with the topology-aware site finder, taking top-strand (forward-strand) cut positions to avoid double-counting palindromic sites, and emitting fragments bounded by those cuts together with the blunt/5'/3' chemistry of each cut end. Both linear and circular (plasmid) molecules are supported on every public entry point and over MCP. The same source surface also provides digest summaries, restriction maps, and overhang-compatibility checks.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A restriction digest cleaves DNA at enzyme recognition sequences and produces fragments whose sizes can be checked experimentally by gel electrophoresis. The original document emphasizes the fragment-sum principle: the sum of fragment lengths must equal the original sequence length. It also states that two enzymes produce compatible ligatable ends when both are blunt cutters or when they generate identical sticky-end overhangs. Sources: Wikipedia (Restriction digest, Restriction enzyme, Restriction map), Addgene Restriction Digest Protocol, REBASE.

### 2.2 Core Model

Given a sequence of length `L` and forward-strand cut positions `c1 < c2 < ... < ck`, the digest fragments are the half-open intervals:

$$
[0, c_1), [c_1, c_2), \ldots, [c_k, L)
$$

This yields `k + 1` fragments when `k` cuts are present. `GetDigestSummary(...)` sorts fragment sizes in descending order, and `CreateMap(...)` groups sites by enzyme name while treating unique cutters as enzymes with exactly one forward-strand site.

For a **circular** molecule of length `L` the site search itself wraps past the origin (so an origin-spanning recognition site is found and cut), and the fragments are the cyclic intervals between consecutive cuts:

$$
[c_1, c_2), [c_2, c_3), \ldots, [c_k, L) \cup [0, c_1)
$$

which yields exactly `k` fragments for `k` cuts; `k = 0` leaves the circle whole and `k = 1` linearizes it into one full-length fragment. The origin-spanning fragment is reported last, so `StartPosition` increases monotonically over the returned fragments in both topologies.

### 2.3 Fragment End Chemistry

A cut by an enzyme with forward cut offset `f` and reverse cut offset `r` places the top-strand cut at `watson` and the bottom-strand cut at `crick = watson + (r - f)`. The single-stranded extension spans the top-strand window `[min(watson, crick), max(watson, crick))`, of length `|r - f|`:

| Condition | End |
|-----------|-----|
| `r > f` | 5' overhang |
| `r < f` | 3' overhang |
| `r == f` | blunt |

Each fragment end created by a cut reports this extension, written 5'→3' on the strand that carries it: at one cut the upstream fragment's right overhang and the downstream fragment's left overhang are reverse complements of each other, which is why they re-ligate (and are the same string for a palindromic enzyme). Cross-checked against Biopython 1.88 `Bio.Restriction`: over the recognition sequence this window equals `enzyme.ovhgseq` and the type equals `is_5overhang()`/`is_3overhang()`/`is_blunt()` for all 39 built-in enzymes (39/39, 2026-10-14) — EcoRI `AATT` 5', PstI `TGCA` 3', EcoRV `""` blunt. Taken over the molecule rather than the pattern, a degenerate enzyme reports its real bases (SfiI yields the three actual bases, not Biopython's `NNN`).

### 2.3a Partial Digestion (`PartialDigest`, B08 F53)

A partial digest leaves any subset of the cut sites uncut. It needs no data beyond the complete digest's cut list. For a linear molecule of length `L` with cut positions `c_1 < … < c_k`, the fragments are all intervals `[a, b)`, `a < b`, over `X = {0, c_1, …, c_k, L}`. Their lengths form the difference multiset ΔX of the partial digest problem, with `C(k+2, 2) = (k+1)(k+2)/2` elements (Rosalind PDPL "Creating a Restriction Map", after Jones & Pevzner: "the collection of all positive differences between elements of X"; sample `X = {0,2,4,7,10}` → `ΔX = {2,2,3,3,4,5,6,7,8,10}`).

A circular molecule with `k ≥ 1` cuts gives every cut-to-cut arc, wrapping past the origin, with `0…k−1` sites missed. An arc with `k−1` missed sites is the full-length molecule linearised at its start cut. The intact circle (all `k` missed) is reported as `Digest` reports an uncut circle, so the total is `k² + 1`. Each fragment carries `MissedCuts`, the number of cut positions strictly inside it. Its ends carry the complete digest's enzyme lists and overhangs.

Fragments are ordered by `MissedCuts` and then by `StartPosition`, and `FragmentNumber` is the position in that order. The `MissedCuts = 0` prefix is therefore exactly `Digest`'s output. The optional `maxMissedCuts` keeps the prefix with `MissedCuts ≤ m`. A request for more than `MaxPartialDigestFragments` (100,000) fragments throws `ArgumentException`, and the caller can lower `maxMissedCuts` to stay under it. The model is combinatorial only: the **yield** of each partial product is out of scope. No reference tool implements partial digestion, so validation uses a constructed oracle: the union of complete digests over all `2^k` cut subsets (Evidence §8).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `k` forward-strand cut positions produce `k + 1` fragments | `Digest(...)` inserts boundaries `[0, cuts..., sequence.Length]` |
| INV-02 | The sum of fragment lengths equals the original sequence length | Adjacent boundaries partition the sequence |
| INV-03 | `FragmentSizes` in `DigestSummary` are sorted descending | `GetDigestSummary(...)` orders sizes descending before constructing the record |
| INV-04 | `AreCompatible(A, B) == AreCompatible(B, A)` | Compatibility is determined from bluntness or identical overhang type and sequence |
| INV-05 | `k` top-strand cut positions on a **circle** produce exactly `k` fragments | `DigestCircular(...)` emits one fragment per cut, the last one joining `[lastCut, L)` with `[0, firstCut)` |
| INV-06 | Fragment lengths sum to the molecule length in **both** topologies | Linear and cyclic intervals both partition the molecule |
| INV-07 | `FragmentNumber` runs 1..n contiguously in the returned order | Emission counter, not a boundary index |
| INV-08 | At one cut, the upstream fragment's `RightOverhang` is the reverse complement of the downstream fragment's `LeftOverhang` | Both are read 5'→3' on opposite strands of the same top-strand window |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to digest or map | Null input throws `ArgumentNullException` |
| `enzymeNames` | `string[]` | required for digestion | Enzymes used for site discovery and fragment generation | `Digest(...)` throws `ArgumentException` when none are supplied |
| `topology` | `MoleculeTopology` | `Linear` | Whether the molecule is a linear fragment or a closed circle (plasmid) | Additive overload on `Digest(...)`, `GetDigestSummary(...)` and `CreateMap(...)`; the topology-free overloads are exactly `Linear` |
| `enzyme1Name`, `enzyme2Name` | `string` | required | Enzyme names for compatibility analysis | Unknown names produce `false` in `AreCompatible(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DigestFragment` | record | Top-strand sequence, start position, length, the left/right enzyme name per fragment boundary (the ordinally first one when several enzymes cut there) plus the complete ordinal-ordered lists `LeftEnzymes`/`RightEnzymes`, fragment number, and the overhang type + overhang sequence of each end created by a cut (`null` for an original terminus or an uncut circle) |
| `DigestSummary` | record | Total fragment count, sorted sizes, largest/smallest fragment, average size, enzyme list, and the topology digested |
| `RestrictionMap` | record | Sequence length, sites, grouped positions, forward-strand total site count, unique cutters, and non-cutters |
| `compatible` | `bool` | Whether two enzymes generate compatible ends |

### 3.3 Preconditions and Validation

`Digest(...)` requires a non-null sequence and at least one enzyme name. When no cuts are found, it returns a single fragment equal to the full sequence. `CreateMap(...)` accepts zero enzyme names and then scans the full built-in enzyme catalog. `AreCompatible(...)` returns `false` if either enzyme name is not present in the catalog.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input sequence and enzyme list.
2. Find restriction sites for each requested enzyme **with the requested topology**, so a circular molecule's origin-spanning sites are found and no linear-only cut-feasibility filter is applied to a closed circle.
3. Keep only top-strand (forward-strand) cut positions to avoid double-counting palindromic sites; the reverse-strand record of the same site carries the bottom-strand cut, which is reported as the end overhang rather than as a second boundary.
4. Sort and deduplicate cut positions.
5. Emit fragments between adjacent boundaries — linear: `[0, c1), ..., [ck, L)`; circular: `[c1, c2), ..., [ck, L) ∪ [0, c1)` — each with the overhang type and overhang sequence of its cut ends.
6. Build summaries or maps by sorting fragment sizes and grouping site positions by enzyme.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Compatibility rules preserved from the original document and confirmed in source:

| Rule | Result |
|------|--------|
| Both enzymes produce blunt ends | Compatible |
| Overhang types differ | Not compatible |
| Overhang types match and overhang sequences match | Compatible |
| Overhang sequences differ | Not compatible |

Examples called out in the original document:

| Pair | Outcome |
|------|---------|
| BamHI and BglII | Compatible (`GATC` overhang) |
| EcoRV and SmaI | Compatible (both blunt) |
| EcoRI and PstI | Not compatible |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Digest` | `O(n + k log k)` | `O(k)` | Site discovery plus sorted cut handling |
| `GetDigestSummary` | `O(n + k log k)` | `O(k)` | Materializes digest fragments and sorts sizes |
| `CreateMap` | `O(n + k log k)` | `O(k)` | Groups positions by enzyme and identifies unique/non-cutters |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs)

- `RestrictionAnalyzer.Digest(DnaSequence, params string[])`: Generates digest fragments of a linear molecule.
- `RestrictionAnalyzer.Digest(DnaSequence, MoleculeTopology, params string[])`: Generates digest fragments of a linear or circular molecule.
- `RestrictionAnalyzer.GetDigestSummary(DnaSequence, params string[])`: Produces sorted fragment statistics for a linear molecule.
- `RestrictionAnalyzer.GetDigestSummary(DnaSequence, MoleculeTopology, params string[])`: Produces sorted fragment statistics for a linear or circular molecule.
- MCP: `restriction_digest(sequence, enzyme_names, topology)` and `digest_summary(sequence, enzyme_names, topology)` delegate to the topology overloads.
- `RestrictionAnalyzer.PartialDigest(DnaSequence, [MoleculeTopology], [int? maxMissedCuts], params string[])` ([RestrictionAnalyzer.PartialDigest.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.PartialDigest.cs), a partial-class file so that existing line anchors into `RestrictionAnalyzer.cs` stay valid): every partial-digest fragment as `PartialDigestFragment(DigestFragment Fragment, int MissedCuts)`; MCP `partial_digest(sequence, enzyme_names, topology, max_missed_cuts)` delegates to it.
- `RestrictionAnalyzer.CreateMap(DnaSequence, params string[])`: Builds a restriction map from requested or all enzymes.
- `RestrictionAnalyzer.AreCompatible(string, string)`: Tests end compatibility.
- `RestrictionAnalyzer.FindCompatibleEnzymes()`: Enumerates compatible enzyme pairs from the built-in catalog.

### 5.2 Current Behavior

The current implementation records forward-strand cut positions only when digesting so that palindromic sites are not double-counted. If no cut positions are found, the digest yields a single fragment equal to the original sequence. When multiple enzymes cut at the same coordinate, the cut position is one boundary and the digest lists every such enzyme once in `LeftEnzymes`/`RightEnzymes` (ordinal order of the table names, so the result is independent of the argument order); `LeftEnzyme`/`RightEnzyme` are the first of that list (review 2026-09 B08 F48). `GetDigestSummary(...)` sorts fragment sizes descending before constructing the record. `CreateMap(...)` counts `TotalSites` on the forward strand only and defines `UniqueCutters` from distinct forward-strand positions, while `NonCutters` are requested enzymes that do not appear in the grouped map.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Fragment generation from ordered restriction cut positions.
- Restriction-map construction by grouping sites per enzyme.
- Compatibility rules for blunt ends and identical sticky-end overhangs.
- Circular (plasmid) digestion, including origin-spanning recognition sites, the uncut circle and the single-cutter linearization.
- Per-fragment-end blunt/5'/3' overhang chemistry taken from the molecule.

**Reference-conformant by design (not a simplification):**

- Digest boundaries are top-strand cut positions only. This is Biopython's own convention: `enzyme.catalyse(seq)` slices the molecule at `enzyme.search(seq)`, the list of top-strand cuts. The reverse-strand record of a palindromic site carries the bottom-strand cut of the *same* site, so using it as a second boundary would split one cut into two. The bottom-strand cut is not discarded: it defines the end overhang reported on both adjoining fragments.

**Not implemented:**

- Gel electrophoresis simulation and other downstream laboratory effects; **users should rely on:** external analysis or visualization tools for those steps.
- Partial (incomplete) digestion and methylation-dependent blocking; **users should rely on:** laboratory controls — the enzyme record carries no cut-efficiency or methylation-sensitivity data to model them from.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No cut sites found | Returns one fragment equal to the original sequence | Explicit special case in `Digest(...)` |
| First fragment | `LeftEnzyme = null` | No enzyme cuts before position 0 |
| Last fragment | `RightEnzyme = null` | No enzyme cuts after the final boundary |
| Zero-length fragments | Cannot arise | The linear cut-feasibility filter confines every cut to `[1, L-1]`, so every interval is non-empty; the emission guard is defensive and `FragmentNumber` counts emitted fragments so the numbering stays contiguous regardless |
| Circular, 0 cut sites | One full-length uncut circular fragment, no ends (overhang fields `null`) | Matches Biopython `catalyse(..., linear=False)`, which returns the whole sequence for a non-cutter |
| Circular, 1 cut site | One full-length linearized fragment starting at the cut | Addgene Plasmids 101: a single cutter linearizes the plasmid |
| Circular, origin-spanning site | Found and cut; the fragment that crosses the origin is reported last | Biopython `search(..., linear=False)` searches `data + data[1:size]` |

### 6.2 Limitations

The digest model is a sequence-partitioning simulation. It does not simulate gel electrophoresis or methylation-dependent blocking; incomplete (partial) digestion is enumerated combinatorially by `PartialDigest` (§2.3a) but the yield of each partial product is not modelled, and it relies on the enzyme definitions available in the built-in restriction catalog (all 39 of which match Biopython 1.88's recognition sequence and both cut offsets exactly, verified 2026-10-14). Circular-DNA behaviour **is** modelled — see §2.2 and §6.1. When different enzymes cut at the same coordinate, all of them are listed, but the overhang reported for that boundary is the first-listed enzyme's own end. Among the built-ins this matters only for FseI/HaeIII (FseI 3' `CCGG`, HaeIII blunt, both cutting after base 6 of `GGCCGGCC`): the physical end of such a combined digest is not modelled.

**Reference-implementation caveat.** Biopython 1.88's `catalyse(..., linear=False)` is not safe to use as a fragment oracle on a circle: `AbstractCut._drop` (Restriction.py L542-568) folds an out-of-range circular cut back by ±length *in place without re-sorting*, so `search(..., linear=False)` can return an unsorted list (e.g. `AluI.search(Seq("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"), linear=False) == [11, 2]`), and `catalyse` (L1146-1176) indexes `r[0]`/`r[-1]` assuming ascending order. The result is an empty fragment plus an over-long one whose lengths sum to more than the molecule (40 vs 31 bases in that example). Measured on 4000 circular sequences with a planted origin-spanning site, 55 (1.4%) trigger it. Seqeron sorts its cut positions, so its fragment lengths always sum to the molecule length; the Seqeron↔Biopython comparisons in `docs/Evidence/RESTR-DIGEST-001-Evidence.md` therefore use `sorted(search(...))` as the oracle for circular cases.

## 8. References

1. Wikipedia: Restriction digest - https://en.wikipedia.org/wiki/Restriction_digest
2. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia: Restriction map - https://en.wikipedia.org/wiki/Restriction_map
4. Addgene: Restriction Digest Protocol - https://www.addgene.org/protocols/restriction-digest/
5. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry.
6. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/
7. Biopython 1.88, `Bio.Restriction` (`Restriction.py`): `AbstractCut._drop` L542-568, `catalyse`/`catalyze` L1131-1178, `ovhgseq`/`is_5overhang`/`is_3overhang`/`is_blunt`. Reference implementation for cut coordinates, circular search and overhang chemistry.
8. Univ. of Illinois MolBio, "Restriction Digestion/Gel Electrophoresis" - fragment count on a circle.
9. Addgene, Plasmids 101: a single cutter linearizes a plasmid - https://blog.addgene.org/plasmids-101-restriction-cloning
10. Rosalind, PDPL "Creating a Restriction Map" (after Jones NC & Pevzner PA, *An Introduction to Bioinformatics Algorithms*, MIT Press 2004) - https://rosalind.info/problems/pdpl/ — partial digest = difference multiset ΔX.
