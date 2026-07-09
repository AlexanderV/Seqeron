---
type: source
title: "Evidence: EPIGEN-BISULF-001 (Bisulfite sequencing methylation calling)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-BISULF-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-BISULF-001-Evidence.md
source_commit: 3fdfb015426652fb931b21c72bd92c5b0a214a7e
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-BISULF-001

The validation-evidence artifact for test unit **EPIGEN-BISULF-001** — **bisulfite sequencing
analysis** (in-silico conversion, per-CpG methylation calling, weighted profile aggregation). This is
the **second ingested unit of the Epigenetics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm is synthesized in its own
concept, [[bisulfite-methylation-calling]]; [[test-unit-registry]] tracks the unit. Sibling of
[[epigenetic-age-horvath-clock]] — this unit *produces* the per-CpG β-values that an age clock consumes.

## What this file records

- **Online sources (all four mutually consistent, no contradictions):**
  - **Frommer et al. (1992)** "A genomic sequencing protocol that yields a positive display of
    5-methylcytosine residues", *PNAS* 89(5):1827–1831 (authority rank 1, primary paper) — bisulfite
    "converts cytosine to uracil, but 5-methylcytosine remains nonreactive"; uracil reads/amplifies as
    **thymine**, 5mC reads as **cytosine** (a remaining C marks methylation); conversion is
    **strand-specific** (one strand at a time, each a separate molecule).
  - **Krueger & Andrews (2011)** "Bismark", *Bioinformatics* 27(11):1571–1572 (rank 1 / 3) — the call
    rule: a **read C at a reference C = methylated** (protected), a **read T = unmethylated**
    (converted); Bismark discriminates CpG / CHG / CHH context.
  - **Bismark User Guide v0.15.0** (Babraham, rank 3) — the percentage formula
    "**methylation % = 100 · count_methylated / (count_methylated + count_unmethylated)**" (= fraction
    `meth/(meth+unmeth)`); call symbols `z/Z` CpG, `x/X` CHG, `h/H` CHH (lowercase unmeth / uppercase
    meth); percentage computed **individually per context**.
  - **Schultz et al. (2012)** 'Leveling' the playing field, *Trends Genet.* 28:583–585 (rank 1) —
    **weighted methylation level** = total methylated reads / total (methylated + unmethylated) reads
    = `Σ(levelᵢ·coverageᵢ)/Σ(coverageᵢ)` over a context's sites.

- **Documented corner cases / failure modes:** non-C bases (A/G/T) unaffected — leave unchanged;
  single-strand conversion only (no reverse-complement merge); **zero-coverage CpG excluded** (denominator
  0 undefined); read bases outside the reference / past its end ignored; a CpG needs both C and the
  following G (last reference base cannot start a CpG); a reference-C read base that is neither C nor T
  (A/G mismatch) is not a valid bisulfite call and is excluded.

- **Datasets (documented oracles):**
  - *Conversion (Frommer chemistry)* — input `ACGTCGAA`, methylated `{1}` → **`ACGTTGAA`** (C@1 protected
    stays C; C@4 unmethylated → T; non-C unchanged).
  - *Calling (Bismark rule)* — reference `ACGTACGT` (CpG at index 1 and 5); at index 1 one `C` + one `T`
    read → level **1/2 = 0.5**, coverage 2; at index 5 one `T` read → level **0/1 = 0.0**, coverage 1.
  - *Weighted profile (Schultz)* — site A (level 1.0, coverage 10) + site B (level 0.0, coverage 30) →
    weighted CpG methylation `(1.0·10 + 0.0·30)/(10+30) = 10/40 = ` **0.25** ≠ unweighted mean **0.5**.

- **Test-coverage recommendations:** MUST — unmethylated C→T / methylated C protected / non-C unchanged;
  level = meth/(meth+unmeth) per CpG with coverage = valid calls and zero-coverage sites excluded;
  per-context weighted level = Σ(level·coverage)/Σcoverage giving 0.25 ≠ 0.5. SHOULD — empty/null →
  empty; lowercase handling; T at a CpG counts as unmethylated. COULD — reads past the reference end
  ignored; a reference A/G at the read position not counted.

## Deviations and assumptions

- **ASSUMPTION (coordinate base):** positions are **0-based** offsets into the supplied sequence
  (matches the sibling `FindCpGSites`/`FindMethylationSites`); Frommer/Bismark are silent on the internal
  index base — an API convention, not correctness-affecting for the values.
- **ASSUMPTION (single-strand conversion):** `SimulateBisulfiteConversion` converts **only** the supplied
  strand and does not synthesise/merge the complement (Frommer's strand-specific protocol — each strand
  is a separate molecule; a two-strand merge is a different operation, out of scope).
- **API-contract deviation (from the algorithm doc):** the registry signature
  `CalculateMethylationFromBisulfite(bsSeq, refSeq)` is realised as `(referenceSequence, bisulfiteReads)`
  because per-site coverage needs per-read multiplicity a single converted string cannot carry
  (accepted).

No source contradictions.
