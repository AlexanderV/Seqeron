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
**Accessed:** 2026-06-13 (retrieved via WebFetch; used for the primary citation, not as authority itself) · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (Wikipedia) → primary: Genome Research 14(1):149–159

**Key Extracted Points:**

1. **Scaffolding definition:** "Link together a non-contiguous series of genomic sequences into a scaffold, consisting of sequences separated by gaps of known length." Scaffolding orders and orients contigs.
2. **Greedy linkage:** Bambus "joins together contigs with the most links first" — the canonical greedy scaffolder follows the strongest links; here links are supplied pre-ordered.

---

### B10 review (2026-11-25) — sources opened and reference numbers

**Huson, Reinert & Myers (2002), "The greedy path-merging algorithm for contig scaffolding", J. ACM 49(5):603–615.**
PDF opened in full: https://publications.mpi-cbg.de/Huson_2002_6349.pdf (rank 1).
- §2.1: each contig gives two nodes (its ends) joined by a contig edge; a mate-edge joins ends of two
  contigs and "induces a relative orientation and approximate distance"; "l(e) may be negative, which
  indicates that the two contigs should overlap by approximately −l(e) base pairs". Weight = evidence (§2, §2.2).
- Algorithm 5.1: "Select all contig edges. for each mate-edge e in descending order of weight: … if P1 ≠ P2
  and we can merge P1 and P2 …"; invariant "every node is adjacent to at most two selected edges";
  merging that interleaves paths ("zipper", Algorithm 5.2) and the happy/unhappy test need σ(e).

**ABySS `MergePaths/MergeContigs.cpp`** (raw.githubusercontent.com/bcgsc/abyss/master, read): `mergeContigs`
asserts d < 0, overlap = −d, takes the last −d bases of the merged sequence and the first −d bases of the
next contig, `createConsensus` (equal ignoring case; N yields to the other base; IUPAC subset → OR; else
fail); on success collapses the overlap; else a gapped `alignOverlap`, else appends `'n'` + contig.

**SSPACE-Basic** (raw.githubusercontent.com/nsoranzo/sspace_basic/master, `SSPACE_Basic.pl`
`mergeContigs`, `bin/PairingAndScaffolding.pl` `computeLayout`, read): `reverseComplement($seq) if
($orient eq "r")`; negative estimate → overlap search (≥ −n 15) and merge, otherwise a single `n`;
positive → `"N" x distance`; ratio rule `-a 0.7` between the two best links, `-k 5` minimum links;
"Reverse all contigs if the whole scaffold is a reverse complement".

**Python oracle** (independent; Biopython 1.8x `Seq.reverse_complement`). Rule: links sorted by
descending weight (stable); link (c1,o1)→(c2,o2) joins end exit(c1) = 3′ if `+` else 5′ and entry(c2)
= 5′ if `+` else 3′; selected iff both ends free and c1, c2 on different paths; each path read from a
terminal, keep the reading with more forward contigs (tie: lowest-index contig forward); gaps: d > 0 →
d N; d < 0 and the −d bases (since the last gap) agree (case-insensitive, N wildcard) → collapse;
otherwise 100 N (AGP). Scaffolds ordered by lowest contig index.

| Case | contigs | links (c1,c2,gap[,rev1,rev2,weight]) | Oracle = C# |
|---|---|---|---|
| link into lower index | AA, TT | (1,0,2) | `TTNNAA` (old code: `AA`,`TT` — link lost) |
| predecessor of contig 0 | AA, CC, GG | (2,0,1),(0,1,1) | `GGNAANCC` (old: `AANCC`,`GG`) |
| second successor of a joined end | AA, CC, GG | (1,0,1),(1,2,1) | `CCNAA`, `GG` (old: `AA`,`CCNGG`) |
| heavier wins | AAAA, CCCC, GGGG | (0,1,2,w1),(0,2,3,w5) | `AAAANNNGGGG`, `CCCC` |
| cycle | AA, CC, GG | (2,0,1),(0,1,1),(1,2,1) | `GGNAANCC` |
| orientation 0+ 1- | AACG, TTGC | (0,1,3,F,T) | `AACGNNNGCAA` |
| orientation 0- 1+ | AACG, TTGC | (0,1,3,T,F) | `GCAANNNAACG` |
| orientation 0- 1- | AACG, TTGC | (0,1,3,T,T) | `TTGCNNNAACG` |
| oriented chain | ACCA, GGTT, CATG | (0,1,2,F,T),(1,2,4,T,F) | `ACCANNAACCNNNNCATG` |
| inconsistent orientation | ACCA, GGTT, CATG | (0,1,2,F,T),(1,2,4,F,F) | `ACCANNAACC`, `CATG` |
| overlap agrees | ACGTAC, TACGGA | (0,1,−3) | `ACGTACGGA` |
| overlap with N | ACGTNC, TACGGA | (0,1,−3) | `ACGTACGGA` |
| overlap disagrees | ACGTAC, TTCGGA | (0,1,−3) | `ACGTAC`+100N+`TTCGGA` |
| overlap on RC | GGATCC, TTTGG | (0,1,−2,F,T) | `GGATCCAAA` |
| overlap then gap | AAACG, CGTTT, GGGG | (0,1,−2),(1,2,3) | `AAACGTTTNNNGGGG` |
| singleton ordering | AAA, CCC, GGG | (2,1,4) | `AAA`, `GGGNNNNCCC` |
| MCP example | AAAA, CCCC, GGTT | (0,1,2),(0,2,3,F,T,w4) | `AAAANNNAACC`, `CCCC` |

Differential run: 5000 random cases (3000 unconstrained: 1–7 contigs over ACGT/N incl. empty, up to 3n
links with out-of-range/self indices, gaps in {0, −6..−1, 1..9}, random orientations, weights 1–4;
2000 genome-derived: 2–6 overlapping pieces, random RC and shuffling, gaps −o / mismatched −o / 0..5)
— C# `Scaffold` = oracle in **5000/5000** (600 cases fully collapse into one gap-free scaffold).
Oracle (verbatim):

<details><summary>scaffold_oracle.py</summary>

```python
"""Independent oracle for ASSEMBLY-SCAFFOLD-001.
Greedy (Huson, Reinert & Myers 2002, Alg. 5.1, end-to-end merges; Bambus 'most links first'):
mate-edges in descending weight (stable on input order); an edge joins two contig END nodes;
accepted iff both end nodes still have no selected mate-edge and they lie in different selected
paths. Paths are reverse-complemented with Biopython. Negative distance d: ABySS MergeContigs
(the -d bases of the tail must agree, N matching anything) -> collapse; otherwise AGP 'U' 100.
"""
from Bio.Seq import Seq

UNKNOWN = 100

def rc(s): return str(Seq(s).reverse_complement())

def scaffold(contigs, links, gap='N'):
    n = len(contigs)
    # link = (c1, c2, gap, rev1, rev2, weight)
    L = [l if len(l) == 6 else (l[0], l[1], l[2], False, False, 1) for l in links]
    order = sorted(range(len(L)), key=lambda i: -L[i][5])
    partner = {}  # (contig,'L'/'R') -> ((contig,end), gap)
    comp = list(range(n))
    def find(x):
        while comp[x] != x: x = comp[x]
        return x
    for i in order:
        c1, c2, g, r1, r2, w = L[i]
        if not (0 <= c1 < n and 0 <= c2 < n) or c1 == c2: continue
        a = (c1, 'L' if r1 else 'R'); b = (c2, 'R' if r2 else 'L')
        if a in partner or b in partner: continue
        if find(c1) == find(c2): continue
        partner[a] = (b, g); partner[b] = (a, g); comp[find(c1)] = find(c2)
    groups = {}
    for c in range(n): groups.setdefault(find(c), []).append(c)
    out = []
    for members in groups.values():
        # terminal: contig with a free end; walk from it with free end on the left
        def walk(start):
            path = []; c = start
            o = '+' if (c, 'L') not in partner else '-'
            g = None
            while True:
                path.append((c, o, g))
                ex = (c, 'R' if o == '+' else 'L')
                if ex not in partner: break
                (nc, ne), g = partner[ex]
                c, o = nc, ('+' if ne == 'L' else '-')
            return path
        terms = [c for c in members if (c, 'L') not in partner or (c, 'R') not in partner]
        cands = [walk(t) for t in terms]
        lo = min(members)
        def key(p):
            fwd = sum(1 for (_, o, _) in p if o == '+')
            lo_fwd = any(c == lo and o == '+' for (c, o, _) in p)
            return (-fwd, 0 if lo_fwd else 1)
        best = min(cands, key=key)
        s = ''; since = 0
        for (c, o, g) in best:
            seq = contigs[c] if o == '+' else rc(contigs[c])
            if g is None:
                s = seq; since = len(seq); continue
            if g > 0:
                s += gap * g + seq; since = len(seq); continue
            ov = -g if g < 0 else 0
            if g < 0 and ov <= since and ov <= len(seq) and all(
                    x.upper() == y.upper() or x.upper() == 'N' or y.upper() == 'N'
                    for x, y in zip(s[len(s)-ov:], seq[:ov])):
                tail = s[len(s)-ov:]
                merged = ''.join(y if x.upper() == 'N' else x for x, y in zip(tail, seq[:ov]))
                s = s[:len(s)-ov] + merged + seq[ov:]; since = since - ov + len(seq)
            else:
                s += gap * UNKNOWN + seq; since = len(seq)
        out.append((lo, s))
    return [s for _, s in sorted(out)]
```
</details>

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

- **2026-11-25** (B10 review): greedy path merging (Huson et al. 2002), orientation, weights, overlap
  collapse (ABySS MergeContigs); sources, oracle and reference table added above. The earlier
  "unresolved-overlap placeholder" assumption is narrowed to overlaps that are not confirmed.
- **2026-06-13**: Initial documentation.
