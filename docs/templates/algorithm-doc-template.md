# <Algorithm Name>

<!--
Canonical template for repository algorithm documentation.

Rules:
- Keep the section order below.
- Remove optional subsections that do not apply; do not leave empty headings.
- Use `## References` only; do not create a separate `Sources` section.
- Keep formal theory/specification separate from repository implementation details.
- Do not repeat the same fact in Overview, Scientific / Formal Basis, and Implementation Notes.
- Back every formula, threshold, invariant, biological claim, or file-format rule with an inline citation.

Link policy:
- Markdown resolves links relative to the current file, not the repo root.
- All in-repo links MUST be file-relative paths (e.g., `../../../src/Foo/Bar.cs`),
  NOT repo-rooted paths (`src/Foo/Bar.cs`). For docs at `docs/algorithms/<group>/<file>.md`
  the prefix to repo root is `../../../`.

Identifier policy:
- Modeling assumptions in 2.3 use IDs `ASM-NN` (zero-padded two-digit, e.g. ASM-01..ASM-99).
- Invariants in 2.4 use IDs `INV-NN` (same format).
- Deviation/Assumption rows in 5.4 may reference the above IDs in their `Notes` column.
- In combined-algorithm files (see below) prefix IDs with the algorithm token, e.g.
  `INV-UPGMA-01`, `ASM-NJ-02`, so identifiers remain unique across the file.

Implementation Status vocabulary (use exactly one in the header table):
- **Production**     — complete; matches the cited theory/spec for in-scope inputs; no known correctness gaps; tests cover the contract.
- **Simplified**     — implements the core mechanism but omits parameters or refinements that the cited source defines; gaps are listed under 5.3 "Intentionally simplified".
- **Reference**      — educational / illustrative implementation; correctness is the priority over performance, scale, or full parameter coverage.
- **Framework**      — provides API and extension points; users must supply domain-specific data (gene sets, signature matrices, scoring tables) for production use.
- **Experimental**   — under active development; contract or output may change without backwards compatibility.

Multi-algorithm files (UPGMA + NJ, Bray-Curtis + Jaccard, ESTIMATE + CIBERSORT, etc.):

The hierarchy is fixed — do not improvise:

- Sections **2** and **4** remain single `## ` (H2) umbrellas with their original
  numbered title. Per-algorithm splits live as H3 children, one per peer method —
  the same pattern extends to suites of three or more peer algorithms (e.g.,
  Tree_Comparison-style triples):
    `### 2.A <Algorithm 1 name>`
    `### 2.B <Algorithm 2 name>`
    `### 2.C <Algorithm 3 name>`  (and so on for D, E, … as needed)
  Inside each H3, use H4 **named, unnumbered** subsections — the position inside the
  block makes the role unambiguous:
    `#### Domain Context`
    `#### Core Model`
    `#### Modeling Assumptions`
    `#### Properties and Invariants`
    `#### Comparison with Related Methods` (optional, may also appear once at end of section 2)
  Section 4 follows the same pattern: `#### High-Level Steps`, `#### Decision Rules / Reference Tables`, `#### Complexity`.

- Section **5.3** remains a single `### ` (H3) umbrella titled `### 5.3 Conformance to Theory / Spec`.
  Per-algorithm splits live as H4 children:
    `#### 5.3.A <Algorithm 1 name>`
    `#### 5.3.B <Algorithm 2 name>`
  Each H4 contains the three required labeled bold blocks (Implemented / Intentionally simplified / Not implemented).

- The letter→algorithm mapping is declared once, as a single quoted line directly under
  the `## 2. Scientific / Formal Basis` heading, e.g.:
    `> A = UPGMA, B = Neighbor-Joining`
  or for 3+ peer suites:
    `> A = Robinson-Foulds, B = Quartet distance, C = Path-difference`
  The same letters then apply across 2.A/2.B/2.C, 4.A/4.B/4.C, 5.3.A/5.3.B/5.3.C.

- ID prefixes in combined files use the **algorithm token**, not the letter, so IDs stay
  readable from tests/evidence: `INV-UPGMA-01`, `ASM-NJ-02`.

- Sections **1, 3, 5.1, 5.2, 5.4, 6, 7, 8** are shared. Split a row inside one of these
  shared tables only when it is genuinely algorithm-specific; in that case prefix the
  row's first cell with the algorithm token in brackets, e.g.
  `[UPGMA] negative branch lengths impossible`.
-->

| Field | Value |
|-------|-------|
| Algorithm Group | <Area / subdomain> |
| Test Unit ID | <ID / N/A> |
| Related Projects | <Project names / N/A> |
| Implementation Status | <Production / Simplified / Reference / Framework / Experimental> |
| Last Reviewed | <YYYY-MM-DD> |

## 1. Overview

<3-6 sentences describing what the algorithm does, what problem it solves, when it should be used, and whether it is exact, heuristic, probabilistic, or specification-driven.>

## 2. Scientific / Formal Basis

<!-- Use biology, mathematics, or a formal specification depending on the domain. -->
<!-- Do not describe repository-specific behavior here. -->

### 2.1 Domain Context

<Only the background required to understand the algorithm. For file parsers and protocol handlers, summarize the governing specification instead of biological theory.>

### 2.2 Core Model

<Formal definition, recurrence, formula, scoring rule, state machine, parsing model, coordinate system, or decision criteria. Cite the originating paper or specification for every non-trivial formula.>

### 2.3 Modeling Assumptions (Optional)

<!--
Use this subsection for statistical, probabilistic, or biological models whose validity
depends on stated preconditions about the data or system being modeled. Examples:
random mating and no selection (Hardy-Weinberg), molecular clock (UPGMA), additive
distances (Neighbor-Joining), nearest-neighbor independence (RNA folding), linear
mixture (CIBERSORT), fixed temperature/buffer (thermodynamic models). Omit for purely
combinatorial or specification-driven algorithms.

IDs: `ASM-NN` (zero-padded). In combined-algorithm files prefix with the algorithm
token (e.g., `ASM-UPGMA-01`).
-->

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | <e.g., random mating> | <e.g., excess homozygotes; HWE rejected> |

### 2.4 Properties and Invariants

<Determinism, optimality guarantees, correctness conditions, monotonicity, value bounds, coordinate conventions, or other facts that must always hold for any valid input. State invariants in a form that can be checked by tests.>

<!--
IDs: `INV-NN` (zero-padded two-digit). In combined-algorithm files prefix with the
algorithm token (e.g., `INV-UPGMA-01`). Invariant IDs are referenceable from tests,
evidence docs, and the Notes column of 5.4.
-->

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | <e.g., 0 ≤ score ≤ 1> | <derivation or citation> |

### 2.5 Comparison with Related Methods (Optional)

<!-- Include only when contrast with a classical alternative is essential to understanding (e.g., Smith-Waterman vs Needleman-Wunsch, UPGMA vs Neighbor-Joining, Bray-Curtis vs Jaccard, NNLS vs ν-SVR). Compare on properties relevant to user choice, not exhaustive trivia. -->

| Aspect | <This algorithm> | <Alternative> |
|--------|------------------|---------------|
| <e.g., tree type> | <rooted> | <unrooted> |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| <param> | <type> | <default / required> | <meaning> | <range, units, indexing, alphabet, normalization> |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| <field> | <type> | <meaning and interpretation, including units and coordinate base where relevant> |

### 3.3 Preconditions and Validation

<How null, empty, malformed, or out-of-range input is handled. Specify indexing conventions (0-based vs 1-based), coordinate system (inclusive/exclusive), accepted alphabet (DNA / RNA / protein / IUPAC degenerate), case sensitivity, normalization (T↔U, lowercase→uppercase), and exception types raised for each failure mode.>

## 4. Algorithm

### 4.1 High-Level Steps

1. <Step 1>
2. <Step 2>
3. <Step 3>

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

<Use when operational details are needed beyond the formal model:
- thresholds, scoring penalties, dynamic programming states, windowing rules, parser states;
- **reference tables / lookup constants** that materially define behavior — IUPAC codes, codon tables, BLOSUM/PAM matrices, Turner 2004 thermodynamic parameters, organism-specific codon usage, restriction-enzyme libraries, RBS spacing constraints;
- core data structures (suffix arrays, trees, hash maps) whose layout affects complexity or output ordering.
Cite the origin of every numeric parameter table.>

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| <main operation> | <Big-O> | <Big-O> | <assumptions, dominant factors, output-size terms, alphabet size> |

## 5. Implementation Notes

<!-- This section documents the actual repository code. Repository-specific behavior only — do not restate Section 2. -->

### 5.1 Location and Entry Points

<!--
Use a path relative to *this* document, not the repo root.
For docs at `docs/algorithms/<group>/<file>.md` the prefix is `../../../`.
Example: `[SequenceAligner.cs](../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)`.
-->

**Implementation location:** [<FileName.cs>](../../../src/<path-from-repo-root>/<FileName.cs>)

- `<Type.Method(...)>`: <what it does>
- `<Type.Method(...)>`: <what it does>

### 5.2 Current Behavior

<Document only repository-specific behavior that is not already captured by the public contract: fallback rules, supported modes, storage/layout choices, caching, internal optimizations (SIMD, log-sum-exp, hybrid storage), and any observable semantics that differ from the abstract algorithm description.>

### 5.3 Conformance to Theory / Spec

<!--
This subsection is structurally enforced: it MUST contain the three labeled blocks
below. Each block is a bullet list; use "(none)" rather than deleting the block.
A free-text paragraph here is non-conforming.
-->

**Implemented (verbatim from the cited theory/spec):**

- <feature or formula realized exactly as in Section 2>

**Intentionally simplified:**

- <feature: <approximation>; **consequence:** <what the user observes that differs from the cited source>>

**Not implemented:**

- <feature: <out-of-scope item>; **users should rely on:** <pointer to the better-conforming class, external tool, or "no current alternative">>

### 5.4 Deviations and Assumptions (Optional)

<!-- Include this subsection only if there are real deviations from the cited theory/spec or implementation-level assumptions that affect correctness, interpretation, interoperability, or expected output. Omit it entirely when not needed. -->

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | <name> | Deviation / Assumption | <why it matters to users> | <open / accepted / fixed / blocked> | <resolution, mitigation, or pointer to the better-conforming class> |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| <edge case> | <result> | <specification clause, implementation contract, or biological reasoning> |

### 6.2 Limitations

<Document unsupported scenarios, ambiguity sources, heuristic weaknesses, performance cliffs, biological caveats (e.g., does not handle pseudoknots, no codon bias modelling, prokaryote-only), interpretation limits, and conditions under which output is undefined or unreliable.>

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

<!-- Either form is acceptable; both is rare and only when each illustrates something the other cannot. -->

**API usage example:**

```csharp
// Minimal example that demonstrates the canonical API and the expected shape of the result.
```

**Numerical / biological walk-through (optional):**

<Step-by-step trace on a small input that shows intermediate quantities: PWM scores per position, pairwise distances, frequency tables, DP cell values, etc. Useful when the formula in 2.2 is not visually obvious from a code snippet.>

### 7.2 Applications and Use Cases (Optional)

<!-- Concrete real-world applications when they are not already obvious from Section 1 (Overview) or Section 2.1 (Domain Context). Examples: forensic STR profiling, primer design, motif discovery, microsatellite-based diagnostics, codon optimisation for heterologous expression. Do not duplicate the Overview. -->

- **<Application>:** <what the algorithm enables in this context, with a citation when the link is non-obvious>

### 7.3 Related Tests, Evidence, or Documents

<!--
Cross-link tests, evidence/QA documents, and related algorithm docs. All paths
must be relative to *this* file (Markdown does not resolve repo-rooted paths).
For docs at `docs/algorithms/<group>/<file>.md` the prefix to repo root is `../../../`.
Tests/evidence may also reference invariant or assumption IDs from sections 2.3 / 2.4.
-->

- Tests: [<TestFile.cs>](../../../tests/<path>/<TestFile.cs>) — covers `INV-01`, `INV-02`
- Evidence: [<TEST-UNIT-ID>-Evidence.md](../../../docs/Evidence/<TEST-UNIT-ID>-Evidence.md)
- Related algorithms: [<Related_Algorithm>](../<group>/<Related_Algorithm>.md)

### 7.4 Change History (Optional)

| Date | Version | Changes |
|------|---------|---------|
| <YYYY-MM-DD> | <version> | <summary> |

## 8. References

<!--
Every entry must be a full bibliographic citation followed by a DOI or stable URL.
Format: <Author or specification owner>. <Year>. <Title>. <Journal / standard / publisher / site>. <DOI or URL>.
Prefer primary literature, official specifications, or canonical project documentation.
For Wikipedia entries, use the article title and stable URL; do not rely on Wikipedia as the only source for a numeric parameter or formula.
-->

1. <Author>. <Year>. <Title>. <Journal / standard / site>. <DOI or URL>
2. <Author>. <Year>. <Title>. <Journal / standard / site>. <DOI or URL>
