---
type: source
title: "Evidence: SEQ-COMPOSITION-001 (sequence composition — nucleotide + amino-acid counts/fractions, GC content)"
tags: [validation, sequence-statistics, composition]
doc_path: docs/Evidence/SEQ-COMPOSITION-001-Evidence.md
sources:
  - docs/Evidence/SEQ-COMPOSITION-001-Evidence.md
source_commit: 2fa9affeb77d7240ffffd91ffd809647c4297484
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-COMPOSITION-001

The validation-evidence artifact for test unit **SEQ-COMPOSITION-001** — **sequence
composition**: exact **nucleotide** counts/fractions (A/T/G/C/U/N/Other partition of
Length), **GC content** `(G+C)/(A+T+G+C+U)`, and **amino-acid** residue counts. It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; [[test-unit-registry]] tracks the unit.

The formula/definition family (base counts, fractions, GC content, residue counts) is
synthesized on the concept [[base-composition]]; this page records only what the artifact
adds. The two **skew** members mentioned in the same doc (GC skew / AT skew) live on the
sibling concept [[nucleotide-composition-skew]].

**Duplicate/consolidated registry entry.** The Change History records that
SEQ-COMPOSITION-001 is a **duplicate Registry entry for the two composition methods already
delivered under SEQ-STATS-001** — consolidated rather than re-implemented (see TestSpec §7).
So this unit is a documentation/registry artifact over an existing implementation, not a new
algorithm. (SEQ-STATS-001 is not yet ingested — see follow-up.)

## What this file records

- **Online sources:**
  - **Biopython `Bio.SeqUtils`** (`gc_fraction`, `GC_skew`; rank 3, reference impl) —
    `gc_fraction` core `gc = sum(seq.count(x) for x in "CGScgs")`, returns `gc/length`, a
    float in **[0, 1]**; **empty sequence ⇒ 0** (length-zero denominator handled);
    **case-insensitive** (`"CGScgs"` includes lowercase). Default `ambiguous="remove"` counts
    only GCS and includes only ACTGSWU in the length.
  - **Wikipedia "GC skew"** (rank 4, tracing to Lobry 1996) — `GC skew = (G−C)/(G+C)`,
    `AT skew = (A−T)/(A+T)`, cumulative-skew peaks = replication terminus/origin. (Skew
    detail lives on [[nucleotide-composition-skew]].)
  - **IUPAC codes** (rank 2, nomenclature standard) — canonical nucleotides **A, C, G, T, U**;
    degenerate R/Y/S/W/K/M/B/D/H/V and **N = any base**; the **20** standard amino-acid
    single-letter codes (A C D E F G H I K L M N P Q R S T V W Y).
- **Datasets (hand-derived worked examples, arithmetic — no library run needed):**
  - Nucleotide: `ATGC` → A/T/G/C/U = 1/1/1/1/0, GC content 2/4 = **0.5**, GC skew 0, AT skew 0;
    `GGGC` → 0/0/3/1/0, GC content **1.0**, GC skew 2/4 = 0.5, AT skew **0** (`a+t=0`);
    `AAUUGGCC` → 2/0/2/2/2, GC content 4/8 = **0.5**, GC skew 0, AT skew 2/2 = **1.0**.
  - Amino-acid: `MKVLWA` → residues M,K,V,L,W,A, Length 6, each count = 1.
- **Corner cases / failure modes:** empty/null ⇒ **all-zero composition**; **no G or C**
  (zero denominator) ⇒ GC skew **0.0** (and AT skew 0 when `a+t=0`); **mixed case** result is
  case-insensitive; **non-canonical letters** (`N`, degenerate codes, `X`) are tracked
  separately from the four/five canonical bases (routed to `CountN`/`CountOther`), not into
  A/T/G/C/U.

## Deviations and assumptions

**One documented assumption — degenerate IUPAC codes are NOT counted toward composition
totals.** Biopython's `gc_fraction` counts `S` toward GC and `W` toward the length
denominator; the repository counts **only A/T/G/C/U** toward GC/AT totals and routes other
letters to `CountN`/`CountOther`. For sequences over the standard {A,T,G,C,U} alphabet (this
unit's scope) the two agree **exactly**; the difference manifests only on degenerate symbols.
Documented as an intentional simplification, not an invented constant. **No source
contradictions** — Biopython, Wikipedia, and IUPAC agree on the formulas and alphabet.

Recommended coverage (from the artifact): MUST — exact A/T/G/C/U/N/Other counts + Length
partition; GC content `(G+C)/(A+T+G+C+U)`; GC skew `(G−C)/(G+C)` incl. a negative case; AT
skew `(A−T)/(A+T)`; empty/null ⇒ all-zero; amino-acid exact residue counts + length. SHOULD —
case-insensitivity; zero-denominator skews ⇒ 0.
