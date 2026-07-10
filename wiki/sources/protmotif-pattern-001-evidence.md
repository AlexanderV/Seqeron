---
type: source
title: "Evidence: PROTMOTIF-PATTERN-001 (Protein pattern-matching methods — PROSITE→regex + FindMotifByProsite + FindDomains)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md
source_commit: 1527877d257a6d630aba1236cf1b6d4c6e184832
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-PATTERN-001

The validation-evidence artifact for test unit **PROTMOTIF-PATTERN-001** — **Protein Pattern
Matching Methods** (`FindMotifByPattern`, `FindMotifByProsite`, `ConvertPrositeToRegex`,
`FindDomains`). It is a **second Evidence artifact over the same PROSITE→regex engine** as
[[protmotif-find-001-evidence]] (PROTMOTIF-FIND-001): same three primitives, same
information-content scoring, no contradictions. Its distinct contributions are (1) the
end-to-end `FindMotifByProsite` (PROSITE string → regex → matches) call path, (2) a sharper
statement of PROSITE PA-line grammar corner cases, and (3) exact numeric IC oracles. The
engine model, PROSITE→regex table, overlap mechanism and scoring live in
[[protein-motif-pattern-search]]; `FindDomains` is covered by
[[protein-domain-and-signal-peptide-prediction]]. One instance of the templated
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]].

## What this file records

- **Online sources:**
  - **PROSITE pattern syntax** (ScanProsite documentation / PROSITE User Manual §IV.E, rank 2):
    the element-by-element PA-line grammar and its regex mapping — `x`→`.`, `[ALT]`→`[ALT]`,
    `{AM}`→`[^AM]`, `A(3)`/`x(3)`→`A{3}`/`.{3}`, `x(2,4)`→`.{2,4}`, `<`→`^`, `>`→`$`, `-`
    separator dropped, trailing period terminates the pattern.
  - **PROSITE worked-example entries (rank 5):** PS00001 ASN_GLYCOSYLATION `N-{P}-[ST]-{P}`
    →`N[^P][ST][^P]`; PS00005 PKC_PHOSPHO_SITE `[ST]-x-[RK]`→`[ST].[RK]`; PS00016 RGD `R-G-D`
    →`RGD`; PS00017 ATP_GTP_A / P-loop `[AG]-x(4)-G-K-[ST]`→`[AG].{4}GK[ST]`; PS00029
    LEUCINE_ZIPPER `L-x(6)-L-x(6)-L-x(6)-L`→`L.{6}L.{6}L.{6}L`.
  - **De Castro et al. 2006 (ScanProsite, rank 1):** PROSITE signatures are "patterns (regular
    expressions) … or generalized profiles (weight matrices)" — confirms patterns are realized
    as regexes.
  - **Schneider & Stephens 1990 (sequence logos, rank 1):** per-position information content
    `Rseq = log₂N − Σ pₙ log₂pₙ`, in bits; protein maximum `log₂20 ≈ 4.32` bits/site; a
    position allowing `k` of 20 residues uniformly gives `IC = log₂(20/k)`.

- **PROSITE PA-line corner cases (sharper than FIND-001):**
  - **Ranges only on `x`:** `x(2,4)` is valid but a range on a residue letter such as
    `A(2,4)` is **not** a valid PROSITE element (fixed counts like `A(3)` remain valid).
  - **Period terminator:** characters after the `.` that ends a PA line are not part of the
    pattern.
  - **Unsupported `*` Kleene star:** the `<{C}*>` metacharacter belongs to the ScanProsite
    *query* extension, not the standard PA-line grammar — a PA-line→regex converter must
    **reject it (FormatException), not silently treat `*` as a residue**.

- **Overlapping enumeration (ASSUMPTION):** `FindMotifByPattern` wraps the regex in a
  zero-width lookahead `(?=(...))` so every start position, including overlaps, is enumerated.
  Marked ASSUMPTION — a documented repository contract, not PROSITE-mandated; it changes only
  whether overlapping starts are all listed, not the set of start positions.

- **E-value (ASSUMPTION):** `CalculateEValue` uses `E = (N−L+1)·2^(−IC)`, the expected number
  of random matches under a uniform i.i.d. amino-acid background. Consistent with the IC
  definition but not the exact ScanProsite E-value (which uses Swiss-Prot residue
  frequencies); tested for its defining formula, not against a ScanProsite number.

## Test datasets (exact oracles)

| Dataset | PROSITE | Regex | Sequence | Expected |
|---------|---------|-------|----------|----------|
| PS00016 RGD (literal) | `R-G-D` | `RGD` | `AAARGDAAA` | start 3 (0-based), substring `RGD`; IC = 3·log₂20 ≈ 12.965784284662087 bits |
| PS00001 N-glyc (exclusion+class) | `N-{P}-[ST]-{P}` | `N[^P][ST][^P]` | `AANASAAANGTAAAA` | starts 2, 8; substrings `NASA`, `NGTA` |
| PS00017 P-loop (fixed range) | `[AG]-x(4)-G-K-[ST]` | `[AG].{4}GK[ST]` | — | conversion only |
| IC scoring | — | — | — | fixed letter k=1 → log₂20 ≈ 4.321928094887363; class `[ST]` k=2 → log₂10 ≈ 3.321928094887362; wildcard `.` k=20 → 0 |

## Recommended coverage

MUST: `ConvertPrositeToRegex` produces the exact regex for each worked example
(PS00001/05/16/17/29) and each syntax atom (`x`, `x(n)`, `x(n,m)`, `[..]`, `{..}`, `A(n)`,
`<`, `>`, trailing `.`); `FindMotifByPattern` returns exact 0-based start/end/substring for a
literal and a class pattern and `Score` equals `Σ log₂(20/kᵢ)`; `FindMotifByProsite`
end-to-end on PS00001 and PS00016 yields the dataset positions; `FindDomains` detects the
P-loop (PS00017-equivalent) at the correct position. SHOULD: overlapping enumeration via
lookahead; unsupported `*` rejected with FormatException. COULD: null/empty inputs return
empty (no exception) for every method.

## Relationship to PROTMOTIF-FIND-001

No contradictions. PROTMOTIF-FIND-001 already validated `FindMotifByPattern` /
`ConvertPrositeToRegex` / `CalculateMotifScore` / `CalculateEValue` and corrected two impl
patterns (PS00007, PS00018). PATTERN-001 revalidates the same engine over a partly overlapping
example set and adds the `FindMotifByProsite` and `FindDomains` surface plus the `*`-rejection
and range-only-on-`x` grammar contracts. Both cite Schneider & Stephens 1990 for IC and De
Castro 2006 for the pattern-is-regex claim; both keep 0-based `MotifMatch.Start`/`End` vs
ScanProsite's 1-based coordinates.
