# Evidence Artifact: ASSEMBLY-SCAFFOLD-001

**Test Unit ID:** ASSEMBLY-SCAFFOLD-001
**Algorithm:** Scaffolding (joining ordered contigs into scaffolds with N-gaps)
**Date Collected:** 2026-06-13

---

## Online Sources

### Jackman et al., ABySS 2.0 — "Scaffolding a genome sequence assembly using ABySS"

**URL:** http://sjackman.ca/abyss-scaffold-paper/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the paper page)
**Authority rank:** 1 (peer-reviewed; ABySS 2.0 published in Genome Research 2017, 27:768–777)

**Key Extracted Points:**

1. **Scaffold construction:** Verbatim — "The sequences of the vertices in a path are concatenated, interspersed with gaps represented by a run of the character N, whose length corresponds to the estimate of the distance between those two contigs." This defines the canonical scaffold = concatenation of ordered contigs separated by runs of `N`.
2. **Gap length = distance estimate:** The number of `N` characters between two contigs equals the estimated distance between them.
3. **Distance estimation:** Verbatim — "A maximum likelihood estimator is used to estimate the distance between the two contigs from the alignments of the paired reads to the contigs." (The estimator is upstream; this unit consumes a supplied gap estimate.)
4. **Negative estimate = overlap:** Verbatim — "It is possible that the distance estimate is negative, indicating that the two contigs should in fact overlap. If such an overlap is indeed found in the contig overlap graph, the two contigs are merged." (When overlap cannot be resolved, no positive run of `N` corresponds to a negative estimate.)

### NCBI AGP Specification v2.1 — gap representation in scaffolds

**URL:** https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
**Accessed:** 2026-06-13 (retrieved via WebFetch)
**Authority rank:** 2 (official NCBI/INSDC file-format specification)

**Key Extracted Points:**

1. **Gap component type:** Gaps are represented with component type `N` (gap with specified size) or `U` (gap of unknown size).
2. **Gap length must be positive:** Verbatim — "Gap lengths must be positive. Negative gaps and gap lines with zero length are not valid."
3. **Unknown / negative gap default:** Verbatim — "For negative gaps, or gaps of unknown size, use U as the component_type and 100 as the gap size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of unknown size." This fixes the placeholder gap length at **100** when the estimate is non-positive.

### Sahlin et al. (2012) — "Improved gap size estimation for scaffolding algorithms"

**URL:** https://academic.oup.com/bioinformatics/article/28/17/2215/246308
**Accessed:** 2026-06-13 (retrieved via WebFetch)
**Authority rank:** 1 (peer-reviewed, Bioinformatics 28(17):2215–2222)

**Key Extracted Points:**

1. **Gap definition:** Verbatim — "The (unknown) distance between c1 and c2 is given by d." The gap is the unknown distance `d` between two contigs in a scaffold.
2. **Negative gap = overlap, and is common:** Verbatim — "The negative gap case frequently occurs since a de Brujin-based assembler splits its contigs at a given node in the de Bruijn graph that leaves an overlap (negative gap) of one k-mer length." Confirms a non-positive gap is a real, expected input class.

### Pop, Kosack & Salzberg (2004) — "Hierarchical Scaffolding With Bambus" (via Wikipedia primary citation)

**URL:** https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics)
**Accessed:** 2026-06-13 (retrieved via WebFetch; used for the primary citation, not as authority itself)
**Authority rank:** 4 (Wikipedia) → primary: Genome Research 14(1):149–159

**Key Extracted Points:**

1. **Scaffolding definition:** "Link together a non-contiguous series of genomic sequences into a scaffold, consisting of sequences separated by gaps of known length." Scaffolding orders and orients contigs.
2. **Greedy linkage:** Bambus "joins together contigs with the most links first" — the canonical greedy scaffolder follows the strongest links; here links are supplied pre-ordered.

---

## Documented Corner Cases and Failure Modes

### From NCBI AGP Specification

1. **Zero-length gap:** Not valid as a gap line; a scaffold gap must have a positive length. A zero estimate therefore cannot map to "0 N"; it falls into the unknown-gap default (100).
2. **Negative gap:** Not valid as a positive run; mapped to the unknown-size default (100 N) unless overlap is resolved.

### From Jackman et al. (2017)

1. **Negative distance estimate:** Indicates contigs should overlap; if overlap is not found/resolved, a positive `N` run does not represent it.

---

## Test Datasets

### Dataset: Constructed scaffold from ordered contigs with a positive gap

**Source:** Jackman et al. (2017) scaffold construction rule (concatenate path contigs interspersed with a run of `N` of length = distance estimate).

| Parameter | Value |
|-----------|-------|
| contigs | `["ACGT", "TTGG", "CCAA"]` |
| links | `[(0,1,3), (1,2,2)]` |
| gapCharacter | `N` (default) |
| Expected scaffold | `ACGT` + `NNN` + `TTGG` + `NN` + `CCAA` = `ACGTNNNTTGGNNCCAA` |
| Expected length | 4 + 3 + 4 + 2 + 4 = 17 |
| Expected scaffold count | 1 |

### Dataset: Non-positive (negative) gap → unknown-gap default

**Source:** NCBI AGP Specification v2.1 (100 = GenBank/EMBL/DDBJ unknown-gap length); Jackman et al. (2017) (negative = overlap, here unresolved).

| Parameter | Value |
|-----------|-------|
| contigs | `["AAAA", "TTTT"]` |
| links | `[(0,1,-5)]` |
| Expected gap run | 100 × `N` |
| Expected scaffold | `AAAA` + (`N`×100) + `TTTT` |
| Expected length | 4 + 100 + 4 = 108 |

---

## Assumptions

1. **ASSUMPTION: Unresolved-overlap placeholder uses the AGP unknown-gap length (100).** Jackman et al. (2017) merge contigs when a negative estimate's overlap is *found*; this unit does not perform overlap resolution, so a non-positive estimate is emitted as a gap of unknown size. The chosen length (100 `N`) is the GenBank/EMBL/DDBJ standard for unknown-size gaps per the NCBI AGP Specification, so the constant itself is source-backed; the decision to fall back to it (rather than resolve the overlap) is the scoping assumption. This does not invent a numeric value.

---

## Recommendations for Test Coverage

1. **MUST Test:** Ordered contigs with positive gaps concatenate into a single scaffold with exact `N` runs (`ACGTNNNTTGGNNCCAA`). — Evidence: Jackman et al. (2017) scaffold construction.
2. **MUST Test:** A positive gap of size *g* emits exactly *g* gap characters; scaffold length = Σ|contig| + Σ gap. — Evidence: Jackman et al. (2017) ("length corresponds to the estimate").
3. **MUST Test:** A non-positive (zero / negative) estimate emits exactly 100 gap characters. — Evidence: NCBI AGP Specification v2.1.
4. **MUST Test:** A custom gap character is used verbatim instead of `N`. — Evidence: gap is "a run of the character N" → parameterized fill character.
5. **MUST Test:** Each contig appears in at most one scaffold; a link to an already-placed contig is skipped. — Evidence: scaffold = a path of distinct contigs (Bambus/ABySS path model).
6. **MUST Test:** Contigs with no links each become their own single-contig scaffold. — Evidence: scaffold is a path; an unlinked contig is a length-1 path.
7. **SHOULD Test:** Out-of-range / self link indices are ignored. — Rationale: malformed link robustness; indices reference contig positions.
8. **SHOULD Test:** Null `contigs` / null `links` throw `ArgumentNullException`. — Rationale: sibling input-validation convention (`MergeContigs`).
9. **COULD Test:** Empty contig list returns an empty result. — Rationale: trivial identity (no contigs → no scaffolds).

---

## References

1. Jackman SD, Vandervalk BP, Mohamadi H, Chu J, Yeo S, Hammond SA, Jahesh G, Khan H, Coombe L, Warren RL, Birol I (2017). ABySS 2.0: resource-efficient assembly of large genomes using a Bloom filter. *Genome Research* 27:768–777. https://genome.cshlp.org/content/27/5/768 (scaffold-construction text retrieved from http://sjackman.ca/abyss-scaffold-paper/)
2. NCBI (National Center for Biotechnology Information). AGP Specification v2.1 — Accessioned Golden Path file format. https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
3. Sahlin K, Street N, Lundeberg J, Arvestad L (2012). Improved gap size estimation for scaffolding algorithms. *Bioinformatics* 28(17):2215–2222. https://academic.oup.com/bioinformatics/article/28/17/2215/246308
4. Pop M, Kosack DS, Salzberg SL (2004). Hierarchical Scaffolding With Bambus. *Genome Research* 14(1):149–159. (Primary citation located via https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics))

---

## Change History

- **2026-06-13**: Initial documentation.
