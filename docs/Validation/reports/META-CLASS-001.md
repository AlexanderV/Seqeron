# Validation Report: META-CLASS-001 — Taxonomic Classification (k-mer / Kraken-style)

- **Validated:** 2026-06-12   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.ClassifyReads(reads, kmerDatabase, k)`, `MetagenomicsAnalyzer.BuildKmerDatabase(referenceGenomes, k)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** 🔧 LIMITED (implements flat **best-hit** count rule; advertises Kraken/LCA but does not realise the Kraken weighted root-to-leaf-path algorithm nor true LCA assignment — documented honestly below)

---

## Stage A — Description

### Sources opened & what they confirm

1. **Wood & Salzberg (2014), Kraken — Genome Biology 15:R46** (PMC4053813, full text fetched).
   The canonical algorithm, quoted verbatim:
   - Database maps each k-mer to **"the LCA of all organisms whose genomes contain that k-mer."**
   - **"These LCA taxa and their ancestors in the taxonomy tree form … the classification tree, a pruned subtree."**
   - **"Each node in the classification tree is weighted with the number of k-mers in K(S) that mapped to the taxon associated with that node."**
   - **"Each root-to-leaf (RTL) path … is scored by calculating the sum of all node weights along the path."**
   - **"The maximum scoring RTL path … is the classification path, and S is assigned the label corresponding to its leaf."**
   - Unmatched: **"Sequences for which none of the k-mers in K(S) are found in any genome are left unclassified."**
   - Ties: **"If there are multiple maximally scoring paths, the LCA of all those paths' leaves is selected."**

2. **Kraken 2 Manual (GitHub wiki, fetched).** Confidence score **C/Q**:
   - C = **"the number of k-mers mapped to LCA values in the clade rooted at the label"** (i.e. winning label **and its descendant clade**).
   - Q = **"the number of k-mers in the sequence that lack an ambiguous nucleotide"** (non-ACGT k-mers excluded — they are not queried).
   - Optional confidence threshold pushes the assignment up the tree toward the root if not met; else unclassified.

3. **Wikipedia — Metagenomics (fetched).** Confirms similarity-based binning (BLAST/MEGAN) and that metagenomes are compared by k-mer profiles; LCA is the standard assignment idea (MEGAN/Kraken). Generic, supportive context only.

### The exact, sourced classification rule

Kraken is **NOT** simple best-hit and **NOT** plain LCA-of-all-hit-taxa. It is a **weighted root-to-leaf-path** rule over a per-read classification tree built from per-k-mer LCAs, with LCA used only as a **tie-breaker** among maximal paths. Because every node's weight is the count of k-mers at that taxon and a path score **sums ancestor weights**, ancestors accumulate the weight of all descendants — so a deep leaf wins only if its whole lineage's k-mer support beats competitors.

### Hand-computed worked examples (against the authoritative rule)

- k-mers map to {SpeciesA⊂GenusG (2 k-mers), SpeciesB⊂GenusG (1 k-mer)} → classification tree: root–…–GenusG(0)–SpeciesA(2) and –SpeciesB(1). Path scores (sum along path): GenusG itself = 3 (its own 0 + nothing; but ancestor accumulation gives GenusG-level score 3 via children), SpeciesA path = 2, SpeciesB path = 1. Max-scoring **leaf** is SpeciesA (2 > 1). If instead SpeciesA=1, SpeciesB=1 → two maximal leaf paths tie → **LCA = GenusG**. → assignment = **GenusG**.
- All k-mers map to SpeciesA only → SpeciesA path dominates → **SpeciesA**.
- Conflicting kingdoms (e.g. Bacteria 1, Archaea 1) → tied maximal paths → **LCA = root** (unclassified-equivalent / cellular root).

### Findings / divergences (Stage A)

The **TestSpec** and **Evidence** docs describe the algorithm's step 6 as *"Classify to taxon with most k-mer hits"* (best-hit) and label this *"Kraken (LCA)"* / *"LCA-like behavior"* (Evidence line 91; TestSpec S2). This conflates three distinct rules. Per the cited primary source they are not equivalent: best-hit ≠ LCA ≠ weighted-RTL. The C/Q formula text is faithfully quoted, but the **C definition** in the docs ("k-mers supporting winning taxon") drops Kraken's *"in the clade rooted at the label"* — i.e. omits descendant-clade accumulation. Hence Stage A is PASS-WITH-NOTES: the *quoted* formulas are correct, but the doc's prose equates best-hit with LCA/Kraken, which the primary source contradicts.

---

## Stage B — Implementation

### Code path reviewed

`MetagenomicsAnalyzer.cs:110-228`.

- **K-mer extraction / canonicalisation** (`ClassifyReads` 127-145, `GetCanonicalKmer` 210-214): extracts every length-k substring, `ToUpperInvariant`, **skips k-mers with any non-ACGT** (matches Kraken ambiguous filtering), canonicalises via `min(kmer, revComp)` ordinal — **correct**.
- **Q / TotalKmers** (137): incremented once per non-ambiguous k-mer → equals Kraken's Q. **Correct.**
- **DB build** (`BuildKmerDatabase` 180-208): canonical k-mer → **single TaxonId string** (first writer wins; no LCA merge of references). One taxon per k-mer, **not** the LCA of all genomes containing it.
- **Classification rule** (154-158): `taxonCounts.OrderByDescending(kv => kv.Value).First()` → **flat best-hit**: the single taxon string with the highest raw k-mer count. No classification tree, no ancestor/RTL summation, no LCA tie-break.
- **C / Confidence** (157-158): `C = bestTaxon.Value` (winning string's own count), `confidence = C / totalKmers`. Correct C/Q arithmetic for a **flat** DB; omits clade-rooted accumulation.
- **Unclassified** (117-122, 147-152): empty/`<k` → Unclassified, TotalKmers=0; no hits → Unclassified, TotalKmers>0. **Matches Kraken.**
- **Taxonomy parse** (216-228): `;` or `|` split into 7 fixed ranks. Implementation-defined but reasonable.

### Worked example recomputed vs the code

| Scenario | Kraken (sourced) | Seqeron code | Match? |
|---|---|---|---|
| SpeciesA=2, SpeciesB=1 (same GenusG) | SpeciesA (RTL leaf) | SpeciesA (best count) | ✓ coincidental |
| SpeciesA=1, SpeciesB=1 (same GenusG) | **GenusG** (LCA of tied leaves) | **arbitrary leaf** (unstable `OrderByDescending`, dictionary order) | ✗ |
| Bacteria=1, Archaea=1 (conflicting kingdoms) | **root / unclassified-equivalent** | **arbitrary kingdom** picked | ✗ |
| All k-mers → SpeciesA | SpeciesA | SpeciesA | ✓ |
| No k-mer matches | Unclassified | Unclassified | ✓ |

The code's behaviour matches Kraken **only** when one taxon has a strict plurality. On ties / cross-clade conflict it returns an arbitrary best-hit rather than the LCA. It never climbs the taxonomy tree.

### Cross-verification of the test-asserted numbers

All 27 tests' exact numbers were recomputed and are **arithmetically self-consistent** with the implemented best-hit + C/Q rule:
- M10: 21bp,k=14 → Q=8, C=1, conf=0.125 ✓
- M12: Q=2, C=1, conf=0.5 ✓
- S4: N at pos 9, k=10 → Q=1 (10 of 11 k-mers span N), C=1, conf=1.0 ✓
- S5: 4 k-mers, Taxon1=3 / Taxon2=1 → winner Taxon1, **C=3 (winning taxon only), Q=4, conf=0.75** ✓

S5 is the telling test: it locks C = winning-taxon-only count, which is correct for this **flat single-level** DB but is **not** Kraken's clade-rooted C in a hierarchical DB.

### Variant / delegate consistency

`ClassifyReads` and `BuildKmerDatabase` use the same `GetCanonicalKmer` and the same ACGT filter on both build and query sides — consistent. No `*Fast` variant.

### Test-quality audit

Tests assert **exact** sourced-style values (not just "no throw"), are deterministic, and cover the Stage-A edge cases (empty, `<k`, no-match, ambiguous-N, canonical RC lookup, mixed case, multi-taxon). **Gap:** no test exercises a **hierarchical** DB (k-mer → LCA at an internal rank) or a **tie** producing an LCA — precisely the cases where best-hit diverges from Kraken/LCA. The suite therefore cannot detect the divergence, which is why it is green.

### Findings / defects (Stage B)

**Divergence (documented, not a within-session-fixable defect):** `ClassifyReads` implements **flat best-hit** (`OrderByDescending(count).First()`), while the algorithm it advertises (Kraken, Wood & Salzberg 2014) uses a **weighted root-to-leaf-path** assignment with **LCA tie-breaking**, and its DB maps each k-mer to the **LCA of all genomes** containing it. The Seqeron DB maps each k-mer to a single (first-seen) taxon string with no LCA merge and no parent links. Confidence C omits Kraken's *clade-rooted* accumulation. Consequences: on ambiguous/tied/cross-clade reads it returns an arbitrary leaf rather than the LCA, and it cannot "pull up" a read to a higher rank.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — primary formulas (C/Q, ambiguous exclusion, unclassified) are correctly cited; the docs' prose wrongly equates *best-hit* with *LCA/Kraken weighted-RTL* and drops "clade rooted at the label" from C.
- **Stage B:** PASS-WITH-NOTES — the code **faithfully realises its own (simplified) description** and all 27 tests pass with exact values; but that description is a **best-hit** approximation of Kraken, not the weighted-RTL/LCA algorithm it names.
- **End-state: 🔧 LIMITED.**

**Why LIMITED (root cause + what's missing):** A complete fix is out of reach in this session because it requires a **data-model redesign**, not a local edit:
1. `BuildKmerDatabase` would need to map each canonical k-mer to the **LCA node of all reference genomes containing it** (requires a taxonomy DAG with parent pointers and an LCA routine over reference TaxonIds), replacing the current `Dictionary<string,string>` first-writer-wins leaf mapping.
2. `ClassifyReads` would need to build a per-read **classification tree**, weight nodes by k-mer counts, **sum weights along each RTL path**, pick the max-scoring leaf, and break ties by **LCA of tied leaves** — replacing the one-line `OrderByDescending(count).First()`.
3. Confidence C would need to become the **clade-rooted** count (winning label plus descendants), per Kraken 2.

Such a change alters the public method contracts and the meaning of every one of the 27 locked tests (S2/S5 explicitly assert best-hit semantics), so it constitutes a larger redesign rather than a defect fix, and must not be done by silently rewriting tests to a different rule. No code was changed; build + all 27 META-CLASS-001 tests remain green.

**Recommendation (future work):** Either (a) implement true Kraken weighted-RTL + LCA DB with a taxonomy tree input and re-derive the test expectations from the Kraken algorithm, or (b) keep the simplified best-hit method but **correct the TestSpec/Evidence/XML-doc wording** to state plainly "flat best-hit (highest-count taxon) over a k-mer→leaf-taxon database; C/Q confidence; this is a simplification of Kraken's weighted root-to-leaf-path/LCA algorithm" so the advertising matches the behaviour.

### Sources
- [Wood & Salzberg (2014), Kraken — PMC4053813](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053813/)
- [Kraken: Genome Biology 15:R46](https://link.springer.com/article/10.1186/gb-2014-15-3-r46)
- [Kraken 2 Manual (confidence C/Q)](https://github.com/DerrickWood/kraken2/wiki/Manual)
- [Wikipedia — Metagenomics](https://en.wikipedia.org/wiki/Metagenomics)
