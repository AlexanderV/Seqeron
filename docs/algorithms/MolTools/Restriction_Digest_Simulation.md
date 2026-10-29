# Restriction Digest Simulation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | N/A |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-10-15 |

## 1. Overview

Restriction digest simulation models DNA fragmentation after cleavage by one or more restriction enzymes. In this repository, digestion is performed by locating cut sites with the topology-aware site finder, taking each site's top-strand cut position (one per recognition window, so a palindromic site is not double-counted; the reverse-strand sites of a non-palindromic enzyme cut on their own), and emitting fragments bounded by those cuts together with the blunt/5'/3' chemistry of each cut end. Both linear and circular (plasmid) molecules are supported on every public entry point and over MCP. The same source surface also provides digest summaries, restriction maps, and overhang-compatibility checks.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A restriction digest cleaves DNA at enzyme recognition sequences and produces fragments whose sizes can be checked experimentally by gel electrophoresis. The original document emphasizes the fragment-sum principle: the sum of fragment lengths must equal the original sequence length. It also states that two enzymes produce compatible ligatable ends when both are blunt cutters or when they generate identical sticky-end overhangs. Sources: Wikipedia (Restriction digest, Restriction enzyme, Restriction map), Addgene Restriction Digest Protocol, REBASE.

### 2.2 Core Model

Given a sequence of length `L` and forward-strand cut positions `c1 < c2 < ... < ck`, the digest fragments are the half-open intervals:

$$
[0, c_1), [c_1, c_2), \ldots, [c_k, L)
$$

This yields `k + 1` fragments when `k` cuts are present **on a linear molecule**; see §2.2 for the circular count, which is `k` for `k >= 1` and `1` for `k = 0`. `GetDigestSummary(...)` sorts fragment sizes in descending order, and `CreateMap(...)` groups sites by enzyme name while counting a cleavage per cut, so `UniqueCutters` are the enzymes that cut the molecule exactly once (B08 H32/H35).

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
| INV-01 | **Linear** molecules: `k` top-strand cut positions produce `k + 1` fragments | `Digest(...)` inserts boundaries `[0, cuts..., sequence.Length]` |
| INV-02 | The sum of fragment lengths equals the original sequence length | Adjacent boundaries partition the sequence |
| INV-03 | `FragmentSizes` in `DigestSummary` are sorted descending | `GetDigestSummary(...)` orders sizes descending before constructing the record |
| INV-04 | `AreCompatible(A, B) == AreCompatible(B, A)` | Compatibility is determined from bluntness or identical overhang type and sequence |
| INV-05 | **Circular** molecules: `k >= 1` top-strand cut positions produce exactly `k` fragments; `k = 0` produces **1** (the intact circle) | `DigestCircular(...)` emits one fragment per cut, the last one joining `[lastCut, L)` with `[0, firstCut)`; with no cut there is no boundary to open the circle at, so the whole molecule is returned as one fragment. Measured on a circular 20-mer with EcoRI: k=0 → 1 fragment of 20, k=1 → 1 fragment of 20, k=2 → 2 fragments of 10 (B08 H35; the carve-out is also stated on `DigestCircular`'s own XML remark and pinned by `RestrictionAnalyzer_Digest_Tests`) |
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
| `enzyme1Name`, `enzyme2Name` | `string` | required | Enzyme names for compatibility analysis | Unknown or empty names throw `ArgumentException` in `AreCompatible(...)` (`TryAreCompatible(...)` answers `false`); a null name throws `ArgumentNullException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DigestFragment` | record | Top-strand sequence, start position, length, the left/right enzyme name per fragment boundary (the ordinally first one when several enzymes cut there) plus the complete ordinal-ordered lists `LeftEnzymes`/`RightEnzymes`, fragment number, and the overhang type + overhang sequence of each end created by a cut (`null` for an original terminus or an uncut circle) |
| `DigestSummary` | record | Total fragment count, sorted sizes, largest/smallest fragment, average size, enzyme list, and the topology digested |
| `RestrictionMap` | record | Sequence length, sites, grouped positions, total cleavage count (`TotalSites`, F108: a palindromic site counts once, a two-sided Type IIB site counts once, a Type IIS enzyme with sites on both strands counts each), unique cutters (exactly one cleavage), and non-cutters |
| `compatible` | `bool` | Whether two enzymes generate compatible ends |

### 3.3 Preconditions and Validation

`Digest(...)` requires a non-null sequence and at least one enzyme name. When no cuts are found, it returns a single fragment equal to the full sequence. `CreateMap(...)` accepts zero enzyme names and then scans the full built-in enzyme catalog. `AreCompatible(...)` throws `ArgumentException` if either enzyme name is not a REBASE enzyme with a defined cut (since B08 F31; `TryAreCompatible(...)` is the lenient form that returns `false`).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input sequence and enzyme list.
2. Find restriction sites for each requested enzyme **with the requested topology**, so a circular molecule's origin-spanning sites are found and no linear-only cut-feasibility filter is applied to a closed circle.
3. Key every site by its top-strand cut, de-duplicating by **cut** rather than by window. A reverse-strand record is skipped only when reading its window backwards cleaves the same bonds as reading it forwards — which is the case for every palindromic enzyme with symmetric cuts, so a palindromic site is not double-counted and its bottom-strand cut is reported as the end overhang rather than as a second boundary. Any other reverse-strand record contributes its top-strand cut `Position + L − CutPositionReverse` (Biopython `location − fst3`): the reverse-strand sites of a non-palindromic enzyme, and also the second cleavage of a window that reads as the site on both strands and cuts differently each way. That last case is EMBOSS 6.6.0's rule, not Biopython's — `embPatRestrictScan` scans `substr` and then `revstr` and keeps one hit per distinct `cut1`, while Biopython reports one reading per window — and it affects exactly 5 of the 754 REBASE enzymes: HauII, SgeI, MspJI, LpnPI, SgrTI (review 2026-09 B08 F55, revised by H7/F101).
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
- `RestrictionAnalyzer.PartialDigest(DnaSequence, MoleculeTopology, int?, IReadOnlyList<Methylase>?, params string[])` and `PartialDigest(DnaSequence, MoleculeTopology, int?, IEnumerable<RestrictionEnzyme>, IReadOnlyList<Methylase>?)` (B08 round 3 G3 / F83; [RestrictionAnalyzer.Methylation.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.Methylation.cs)): the same partial digest over the methylation-aware cut list (§2.3a, `Restriction_Site_Detection.md` §5.2b).
- `RestrictionAnalyzer.Digest(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>)` and `RestrictionAnalyzer.PartialDigest(DnaSequence, MoleculeTopology, int?, IEnumerable<RestrictionEnzyme>)` (review 2026-09 B08 F55): the same digests with caller-supplied enzyme records, e.g. single-sided Type IIS cutters (`BsaI` = `("GGTCTC", 7, 11)`); the enzymes are validated eagerly. `Digest` with records is exposed over MCP as `restriction_digest_custom_enzyme` (B08 H28); `PartialDigest` with records is not.
- `RestrictionAnalyzer.CreateMap(DnaSequence, params string[])`: Builds a restriction map from requested or all enzymes.
- `RestrictionAnalyzer.AreCompatible(string, string)`: Tests end compatibility.
- `RestrictionAnalyzer.FindCompatibleEnzymes()`: Enumerates compatible enzyme pairs from the built-in catalog.

### 5.2 Current Behavior

The current implementation records one top-strand cut per distinct cleavage when digesting, so palindromic sites with symmetric cuts are not double-counted while the reverse-strand sites of a non-palindromic (e.g. Type IIS) enzyme are cut, and a window that admits two different cleavages (HauII, SgeI, MspJI, LpnPI, SgrTI) contributes both. If no cut positions are found, the digest yields a single fragment equal to the original sequence. When multiple enzymes cut at the same coordinate, the cut position is one boundary and the digest lists every such enzyme once in `LeftEnzymes`/`RightEnzymes` (ordinal order of the table names, so the result is independent of the argument order); `LeftEnzyme`/`RightEnzyme` are the first of that list (review 2026-09 B08 F48). `GetDigestSummary(...)` sorts fragment sizes descending before constructing the record. `CreateMap(...)` counts `TotalSites` as the number of distinct cleavages — every forward-strand record, plus every reverse-strand record whose backwards reading cuts other bonds — and defines `UniqueCutters` as the enzymes contributing exactly one, so a Type IIS enzyme's reverse-strand site counts and an enzyme that binds one window but cuts it twice is not a unique cutter; `NonCutters` are requested enzymes that do not appear in the grouped map (B08 H32).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Fragment generation from ordered restriction cut positions.
- Restriction-map construction by grouping sites per enzyme.
- Compatibility rules for blunt ends and identical sticky-end overhangs.
- Circular (plasmid) digestion, including origin-spanning recognition sites, the uncut circle and the single-cutter linearization.
- Per-fragment-end blunt/5'/3' overhang chemistry taken from the molecule.

**Reference-conformant by design (not a simplification):**

- Digest boundaries are top-strand cut positions only. This is Biopython's own convention: `enzyme.catalyse(seq)` slices the molecule at `enzyme.search(seq)`, the list of top-strand cuts. The reverse-strand record of a palindromic site carries the bottom-strand cut of the *same* site, so using it as a second boundary would split one cut into two. The bottom-strand cut is not discarded: it defines the end overhang reported on both adjoining fragments.
- Distal (single-sided Type IIS) cutters (review 2026-09 B08 F55): the end window `[min(watson, crick), max(watson, crick))` with `crick = watson + (CutPositionReverse − CutPositionForward)` holds for a site on either strand and lies in the flanking molecule, so a `BsaI` end reports its four real bases (Biopython's enzyme-level `ovhgseq` is only `NNNN`). Linear sites whose cut would fall off the molecule are dropped (Biopython `_drop`). Differential vs Biopython 1.88 over all 153 single-sided distal enzymes (+ 576 in-site as control) × random sequences × both topologies — fragment start/length, both end overhangs and the partial-digest fragment multiset: **0 mismatches** after the fix; before it, 8,424 distal and 301 in-site non-palindromic cases had wrong fragments because reverse-strand sites were never cut (Evidence RESTR-DIGEST-001 §9). Two-sided Type IIB cutters (`BaeI`) are not expressible.
- Two-sided (Type IIB) cutters (review 2026-09 B08 F56): an enzyme with `SecondCut` contributes two top-strand cuts per site, so the recognition site is excised as its own fragment (BaeI on `A×20 ACAAAAGTACC T×20`: 10/33/8, Biopython `catalyse`); each end's window uses that cut's own pair. Boundaries are the sorted distinct cuts — Biopython's `Palindromic._search` does not sort a two-cut list, so its `catalyse` mis-slices with two or more sites (UcoMSI on two sites: an empty fragment, lengths not summing to the molecule); not reproduced. `GetDigestSummary` and `CreateMap` take `IEnumerable<RestrictionEnzyme>` as well (F58).
- Caller-supplied enzymes (B08 A6-3, A6-5, A7-2, A7-3): the two strand cuts of a custom `RestrictionEnzyme` may be at most `RestrictionAnalyzer.MaxCutStagger` = 64 bases apart (the stagger is the overhang length, materialised for every fragment end; REBASE's largest of the 754 enzymes is 9), its recognition sequence at most `MaxRecognitionLength` = 64 bases (REBASE: 15) and its name at most `MaxEnzymeNameLength` = 64 characters (REBASE: 9) - checked at the call by every member that takes an enzyme or enzyme list; the second cut pair's stagger is compared in `long` with the first's. Cut offsets themselves may be any `int` (`int.MinValue` included): on a circular molecule each cut is reduced modulo the length in `long` arithmetic, for both strands' records, and so is the test that decides whether the reverse reading of a site cuts the same bonds as the forward one (`ReverseReadingCutsAlike`, which on a circle compares the cuts modulo the length: offsets `(1, 5)` and `(13, 17)` cut a 12-mer at the same two bonds and count as one cleavage in the map, as they already did in the digest; its `int` arithmetic wrapped for an offset within the site length of `int.MinValue` and dropped a real cleavage from the circular digest, map and partial digest; `FindSites` was right). Exact arithmetic for `GAATTC` at 3 of a 12-mer, 2^31 = 8 (mod 12): cut offsets `(int.MinValue+4, int.MinValue+2)` give top-strand cuts 11 and 3, i.e. the fragments `3+8 GAATTCAA` and `11+4 AAAA`.

**Not implemented:**

- Gel electrophoresis simulation and other downstream laboratory effects; **users should rely on:** external analysis or visualization tools for those steps.
- ~~Partial (incomplete) digestion and methylation-dependent blocking~~ *(superseded: partial digestion — F53 `PartialDigest`; methylation-dependent blocking — 2026-10-21, B08 round 2 E2 / F68, opt-in `Digest(sequence, topology, methylases, enzymeNames)` / `Digest(sequence, topology, enzymes, methylases)`, `GetDigestSummary` and `CreateMap` likewise, EMBOSS `restrict -methylation` semantics: a recognition window containing a base methylated by the chosen set — `none`/`dam`/`dcm`/`dam+dcm`/`emboss` — is not a cut site; see `Restriction_Site_Detection.md` §5.2b. `PartialDigest` takes it since B08 round 3 G3 / F83: `PartialDigest(sequence, topology, maxMissedCuts, methylases, enzymeNames)` and the record form `(sequence, topology, maxMissedCuts, enzymes, methylases)`, MCP `partial_digest` `methylation`; a blocked site is neither a cut nor a missed cut.)* **Still not modelled:** a per-enzyme sensitivity table (which enzyme Dam/Dcm/CpG actually blocks — the mode is a window rule, so BamHI is blocked by Dam as in EMBOSS), cut efficiency, and kinetics; **users should rely on:** laboratory controls for those.

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

The digest model is a sequence-partitioning simulation. It does not simulate gel electrophoresis, cut efficiency or kinetics. Methylation-dependent blocking is **opt-in** (default off): `Digest`, `PartialDigest`, `GetDigestSummary` and `CreateMap` (and the MCP `restriction_digest`, `partial_digest`, `digest_summary`, `restriction_map` via `methylation` = `none`\|`dam`\|`dcm`\|`dam+dcm`\|`emboss`) take a set of methylase sites and drop every recognition window that contains a base the set methylates, EMBOSS `restrict -methylation` semantics (`Restriction_Site_Detection.md` §5.2b; B08 F68, F83). It is a **window rule, not a per-enzyme sensitivity table**: EMBOSS blocks BamHI by Dam although BamHI is insensitive at the bench, and `DpnI` (which needs methylation) is not modelled as such. Incomplete (partial) digestion is enumerated combinatorially by `PartialDigest` (§2.3a) but the yield of each partial product is not modelled, and it relies on the enzyme definitions available in the built-in restriction catalog or supplied by the caller as `RestrictionEnzyme` records (F55; two-sided Type IIB cutters such as `BaeI` via `RestrictionEnzyme.SecondCut`, F56) (all 39 built-ins match Biopython 1.88's recognition sequence and both cut offsets exactly, verified 2026-10-14). Circular-DNA behaviour **is** modelled — see §2.2 and §6.1. When different enzymes cut at the same coordinate, all of them are listed, but the overhang reported for that boundary is the first-listed enzyme's own end. Among the built-ins this matters only for FseI/HaeIII (FseI 3' `CCGG`, HaeIII blunt, both cutting after base 6 of `GGCCGGCC`): the physical end of such a combined digest is not modelled.

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
