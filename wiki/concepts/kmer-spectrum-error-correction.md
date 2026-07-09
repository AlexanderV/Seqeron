---
type: concept
title: "K-mer spectrum read error correction (two-sided, Musket/Quake)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
source_commit: 4167da78fdda4821db88c34df210958d7da43cde
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-correct-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-CORRECT-001 ... K-mer spectrum (two-sided) read error correction"
      confidence: high
      status: current
---

# K-mer spectrum read error correction

Correcting substitution errors in sequencing reads using the **k-mer spectrum**: real genomic
k-mers recur at high coverage, while k-mers carrying a sequencing error are rare. This is the
error-correction preprocessing step ahead of assembly, and the anchor for the assembly **CORRECT**
family. Validated under test unit **ASSEMBLY-CORRECT-001**; the validation record is
[[assembly-correct-001-evidence]], and [[test-unit-registry]] tracks the unit. See
[[algorithm-validation-evidence]] for the artifact pattern.

## Trusted vs untrusted (the spectrum model)

Traced verbatim to Musket (Liu et al. 2013) and Quake (Kelley et al. 2010):

- A **k-mer is trusted** iff its multiplicity in the read set exceeds a coverage **cut-off**
  (`minKmerFrequency`), otherwise **untrusted**. In the reference tools the cut-off is chosen
  automatically from the valley between the error peak and the genomic peak of the k-mer coverage
  histogram; here it is an explicit parameter.
- A **base is trusted** iff it is covered by *any* trusted k-mer. Bases inside a trusted k-mer are
  presumed correct and are never modified.
- Song & Florea (2018) name the same split *solid* (frequent, error-free) vs *weak* k-mers.

## Two-sided correction rule (per position)

1. **Locate the error** to positions covered *only* by untrusted k-mers (Quake: the intersection,
   then union, of the read's untrusted k-mers).
2. For a candidate position `i`, find an alternative base that makes **all** k-mers covering `i`
   trusted — evaluating **both** the leftmost and the rightmost covering k-mers (this two-sidedness
   is what disambiguates the fix).
3. **At most one substitution per k-mer** is assumed (conservative); the search is over the four
   nucleotides (substitution / single-base-edit model — no indels).
4. **Ambiguity → unchanged.** If more than one alternative base makes both the leftmost and
   rightmost covering k-mers trusted, the base is left **unchanged** (no guess).
5. **No correcting set → unchanged.** If no base makes the covering k-mers trusted, the base is
   left unchanged (Quake may additionally trim/discard — out of scope for this unit).

Because the model only substitutes bases, **read count and per-read length are preserved**; correct
reads pass through untouched. Reads shorter than `k` contribute no k-mers and are returned as-is.

## Limitations (source-documented)

- The frequency cut-off cannot perfectly separate erroneous from correct k-mers — a solid k-mer may
  contain an error, a weak k-mer may be error-free (a general limit of all k-spectrum methods).
- Multiple errors close together within a single k-mer violate the ≤1-error assumption and may be
  uncorrectable.

## Repository defaults (non-behavioral)

The library exposes both parameters with fixed defaults `kmerSize=15`, `minKmerFrequency=2` rather
than auto-selecting them from the coverage histogram. This is an **assumption record**, not a
deviation from the algorithm: every behavioral test passes `k` and the cut-off explicitly, so
correctness for a given `(k, cut-off)` is fully source-defined. Contract behaviour: null reads ⇒
`ArgumentNullException`; `kmerSize < 1` ⇒ `ArgumentOutOfRangeException`; input is upper-cased
(case-insensitive) per the analyzer convention.
