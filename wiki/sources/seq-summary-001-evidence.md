---
type: source
title: "Evidence: SEQ-SUMMARY-001 (SummarizeNucleotideSequence — per-sequence summary record aggregating composition + GC + Shannon entropy + linguistic complexity + Tm)"
tags: [validation, sequence-statistics, composition]
doc_path: docs/Evidence/SEQ-SUMMARY-001-Evidence.md
sources:
  - docs/Evidence/SEQ-SUMMARY-001-Evidence.md
source_commit: 37c18482dfb0ee53a0be1fc073c2b4c6694012cc
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-SUMMARY-001

The validation-evidence artifact for test unit **SEQ-SUMMARY-001** — the **per-sequence summary
record** `SummarizeNucleotideSequence` that bundles, into one result object, the per-sequence
statistics computed by already-validated sibling methods:

| Summary field | Value | Canonical per-metric method (already validated) |
|---------------|-------|--------------------------------------------------|
| `Length` | sequence length | — |
| `Composition` | A/T/G/C/U/N counts | `CalculateNucleotideComposition` → [[base-composition]] |
| `GcContent` | `(G+C)/(A+T+G+C+U)` ∈ [0,1] | `CalculateGcContent` → [[base-composition]] |
| `Entropy` | Shannon `H = −Σ pᵢ log₂ pᵢ` (bits) | [[seq-entropy-profile-001-evidence]] / [[windowed-sequence-complexity-profile]] |
| `Complexity` | linguistic complexity (vocabulary-usage) | `CalculateLinguisticComplexity` → [[windowed-sequence-complexity-profile]] |
| `MeltingTemperature` | Wallace (len<14) / GC-Marmur-Doty (len≥14) | `CalculateMeltingTemperature` → [[primer-dimer-thermodynamics-tm]] (legacy Wallace/Marmur default) |

It is one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; [[test-unit-registry]] tracks the unit.

## A pure aggregation — the field-consistency contract

**SEQ-SUMMARY-001 introduces no new computation.** Its correctness is defined **field-by-field**:
every summary field MUST equal the value its canonical per-metric method returns on the *identical*
input. Any divergence (different rounding, alphabet handling, or formula selection) is a defect of
the summary, **independent** of whether each underlying metric is itself correct — those are the
scope of the sibling units (composition → [[base-composition]] / SEQ-COMPOSITION-001; entropy →
SEQ-ENTROPY-PROFILE-001; linguistic complexity → SEQ-COMPLEX-WINDOW-001; Tm → the SEQ-TM-001
`CalculateMeltingTemperature`, not yet separately ingested). So this page adds no new concept: every
aggregated method is already synthesized on an existing concept page.

**Scope note / cross-reference correction.** [[base-composition]] and [[seq-stats-001-evidence]]
previously described `SummarizeNucleotideSequence` as a *"thin aggregation wrapper re-exposing the
same GC content / counts."* This artifact shows it is broader than composition alone: the record
also carries the **Shannon entropy**, **linguistic complexity**, and **melting temperature** of the
sequence. Not a contradiction — the earlier pages named only the composition facet they cared about
— but the summary is a full per-sequence statistics record, cross-linking the composition, entropy,
complexity, and Tm families in one object.

## What this file records

- **Online sources** (each establishes the canonical definition of one aggregated field, so the
  per-field expectations are evidence-derived):
  - **Biopython `Bio.SeqUtils.gc_fraction`** (rank 3, reference impl) — `gc = Σ count(x) for x in
    "CGScgs"`, returns `gc/length` a float; `length = gc + Σ ATWUatwu`; **empty ⇒ 0**;
    case-insensitive. (GC-content detail on [[base-composition]].)
  - **Biopython `Bio.SeqUtils.MeltingTemp`** (`Tm_Wallace`, `Tm_GC`; rank 3) — Wallace rule
    `Tm = 4·(G+C) + 2·(A+T)` (rule of thumb for ~14–20 nt oligos; docstring example
    `Tm_Wallace('ACGTTGCAATGCCGTA') = 48.0`); GC general form `Tm = A + B·%GC − C/N + salt − D·%mismatch`,
    Marmur & Doty 1962 valueset `Tm = 69.3 + 0.41·%GC − 650/N`.
  - **Wikipedia "Entropy (information theory)"** (rank 4, citing Shannon 1948) — `H(X) = −Σ p(x)
    log p(x)`; **base 2 ⇒ bits**; **maximum `log₂ n`** at the uniform distribution over `n` outcomes.
  - **Wikipedia "Linguistic sequence complexity"** (rank 4, citing Trifonov 1990) — vocabulary usage
    `U` = actual / maximal possible vocabulary size; complexity combines as a **product**
    `C = U₁U₂…Uw`; `0 < C < 1` for DNA fragments.
- **Datasets (hand-derived from the cited formulas — no library run needed):**

  | Sequence | Len | Composition A,T,G,C,U,N | GcContent | Entropy (bits) | MeltingTemperature (°C) | Complexity |
  |----------|-----|--------------------------|-----------|----------------|--------------------------|------------|
  | `ATGCATGC` | 8 | 2,2,2,2,0,0 | **0.5** (4/8) | **2.0** (= log₂4, four symbols @ ¼) | **24.0** (len<14 → Wallace `2·4 + 4·4`) | equals `CalculateLinguisticComplexity` (**0.83968253968…**) |
  | `ATGCATGCATGCATGC` | 16 | GC count 8 | — | — | **43.375** (len≥14 → GC formula `64.9 + 41·(8−16.4)/16`) | — |

  The 16-mer verifies the **length-dispatch to the Marmur-Doty GC branch at length ≥ 14** (repo
  variant `64.9 + 41·(GC − 16.4)/N`, validated in SEQ-TM-001).
- **Corner cases / failure modes:** **empty/null ⇒ the degenerate zero summary** (Length 0,
  GcContent 0, Entropy 0, Complexity 0, Tm 0 — inheriting Biopython's guarded zero-length GC);
  **case-insensitive** (lowercase input gives an identical summary, since each per-metric method
  uppercases internally); Composition is the only nested structure — its A/T/G/C/U/N keys/counts
  must match `CalculateNucleotideComposition` (U and N counted for RNA / ambiguous inputs).

## Deviations and assumptions

**One documented ASSUMPTION — the Tm formula-selection threshold (`useWallaceRule: length < 14`).**
The summary passes `sequence.Length < 14` to `CalculateMeltingTemperature`; the 14 nt boundary is
the sibling SEQ-TM-001 convention (`ThermoConstants.WallaceMaxLength = 14`). Biopython documents
Wallace as a rule of thumb for ~14–20 nt without fixing an exact switch point. This is
**non-correctness-affecting for the summary**: its contract is "MeltingTemperature equals
`CalculateMeltingTemperature` with this flag," so the threshold belongs to the already-validated
SEQ-TM-001 unit, and the summary is tested only for **equality with that canonical method** on the
same input. **No source contradictions** — Biopython, Wikipedia (Shannon), and Trifonov agree on the
per-field formulas; the aggregation is a mechanical bundling of them.

Recommended coverage (from the artifact): MUST — each field equals its canonical per-metric method
on the same input; the exact `ATGCATGC` oracle (GcContent 0.5, Entropy 2.0, Tm 24.0); the GC/Marmur
branch for length ≥ 14 (43.375 for the 16-mer). SHOULD — Composition dict keys/counts match
`CalculateNucleotideComposition`; empty/null ⇒ zero summary. COULD — case-insensitivity.
