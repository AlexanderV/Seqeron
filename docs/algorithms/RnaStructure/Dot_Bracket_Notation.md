# Dot-Bracket (Extended WUSS) Notation

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-DOTBRACKET-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Dot-bracket notation is a one-line, specification-driven encoding of RNA secondary structure: matching brackets denote base pairs and dots denote unpaired residues [1][2]. This unit parses a notation string into the list of paired index positions (`ParseDotBracket`) and tests whether a string is well-formed (`ValidateDotBracket`). The implementation follows the extended ViennaRNA / WUSS convention, where four bracket families (`()`, `[]`, `{}`, `<>`) and uppercase/lowercase letter pairs are independent pairing systems, so crossing helices (pseudoknots) can be represented [1][2][3]. Both operations are exact, deterministic, single-pass string algorithms.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An RNA secondary structure is a set of base pairs over sequence positions `0..n-1`. Dot-bracket notation serializes such a structure as a string the same length as the sequence: each position is either paired (a bracket/letter) or unpaired (a dot or other single-stranded symbol) [1][2]. In standard notation the bracket structure is properly nested; the extended notation adds further bracket families and letter pairs whose helices may cross, which is how pseudoknots are written [1][2][3].

### 2.2 Core Model

- A base pair `(i, j)` with `i < j` is written as an opening symbol at position `i` and the matching closing symbol of the same family at position `j` [1][2].
- Bracket families and their open/close symbols: `(`/`)`, `[`/`]`, `{`/`}`, `<`/`>` [2][3]. Letter pairs use an uppercase letter as the 5' (opening) partner and the matching lowercase letter as the 3' (closing) partner, e.g. `AAAA....aaaa` [1][2].
- Within one family, closers match openers in last-opened-first-closed (nested) order; the choice of family symbol carries no meaning beyond "left and right partners must match up" [3][4].
- Each family / letter case-pair is an **independent** pairing system: families need not nest with one another, which permits crossing helices (pseudoknots), e.g. `<<<<[[[[....>>>>]]]]` [1][2].
- A dot `.` and other WUSS single-stranded symbols (`-`, `,`, `:`) denote unpaired residues [1][5].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every returned pair has `Position1 < Position2` | an opener is pushed before its closer is read; the closer's index is the larger [1][2] |
| INV-02 | Both indices of a returned pair belong to the same family / letter case-pair | each family has its own stack; a closer only pops its own family [3][4] |
| INV-03 | For a string accepted by `ValidateDotBracket`, `ParseDotBracket` returns exactly one pair per opening symbol | balanced ⇒ every opener is matched [1] |
| INV-04 | `ValidateDotBracket(s)` is true iff every closer matches an earlier unmatched opener of the same family and no opener is left unclosed | per-family stack balance [3][4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dotBracket` | `string` | required | dot-bracket / extended WUSS notation string | 0-based positions; null/empty allowed (treated as empty structure) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ParseDotBracket` → | `IEnumerable<(int Position1, int Position2)>` | base pairs as 0-based (5', 3') index tuples |
| `ValidateDotBracket` → | `bool` | true iff the string is well-formed (per-family balanced) |

### 3.3 Preconditions and Validation

Positions are 0-based indices into the notation string. Null or empty input is treated as a valid, pair-free structure: `ValidateDotBracket(null) == ValidateDotBracket("") == true` and `ParseDotBracket(null)`/`("")` yield no pairs (see §5.4). Letter recognition uses `char.IsLetter`/case via the invariant culture. Any character that is neither a recognized bracket nor a letter (dots, `-`, `,`, `:`, and others) is treated as unpaired and skipped [5]. `ParseDotBracket` is best-effort on malformed input: an unmatched closing symbol is dropped without throwing — callers should test `ValidateDotBracket` first.

## 4. Algorithm

### 4.1 High-Level Steps

1. Maintain one index stack per opening symbol (the four opening brackets and each uppercase letter encountered).
2. Scan left to right. On an opening bracket or uppercase letter, push its index onto that symbol's stack.
3. On a closing bracket, pop the matching opening-bracket stack and emit `(opener, current)`. On a lowercase letter, pop the matching uppercase-letter stack and emit the pair.
4. Any other character (dot or single-stranded WUSS symbol) is unpaired and skipped.
5. `ValidateDotBracket` runs the same scan: a closer with an empty/absent same-family stack ⇒ false; after the scan, any non-empty stack ⇒ false; otherwise true.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Bracket family map (opening ↔ closing), from ViennaRNA / WUSS [2][3]:

| Open | Close |
|------|-------|
| `(` | `)` |
| `[` | `]` |
| `{` | `}` |
| `<` | `>` |

Letter pairs: uppercase `X` opens, matching lowercase `x` closes [1][2]. Data structure: a `Dictionary<char, Stack<int>>` keyed by opening symbol, giving each independent pairing system its own LIFO stack.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ParseDotBracket` | O(n) | O(n) | single pass; stacks hold at most n open positions |
| `ValidateDotBracket` | O(n) | O(n) | single pass; per-family stacks |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.ParseDotBracket(string)`: returns the 0-based base pairs encoded by the notation.
- `RnaSecondaryStructure.ValidateDotBracket(string)`: returns whether the notation is well-formed (per-family balanced).

### 5.2 Current Behavior

Each opening bracket family and each uppercase letter gets its own stack; closers pop only their own family's stack. This is what distinguishes the implementation from a single-counter / single-stack approach and is required to (a) correctly parse crossing families such as `([)]` and (b) reject mismatched families such as `(]` during validation. Letters are matched with uppercase as the opener. Non-bracket, non-letter characters are treated as unpaired and ignored.

**Search reuse (suffix tree):** N/A. This is a single linear scan with a stack, not a substring-search / occurrence-enumeration task, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `()`/`.` base-pair / unpaired encoding and balanced-bracket well-formedness [1][2].
- Four independent bracket families `()`, `[]`, `{}`, `<>` for nested and crossing helices [2][3].
- Uppercase/lowercase letter pairs as independent pseudoknot pairing systems, uppercase = 5' opener [1][2].
- "Partners must match up" rule (closer matches opener of same family) in both parse and validation [3][4].
- Non-bracket WUSS symbols (`-`, `,`, `:`) and dots treated as single-stranded [5].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Semantic WUSS distinctions among unpaired symbols (e.g. hairpin-loop `_` vs bulge `-` vs external `:`): all unpaired symbols are treated identically as "unpaired"; **users should rely on:** the raw symbol if loop-type classification is needed — only pairing is decoded here.
- Conversion to/from a full thermodynamic structure object: see [RNA_Base_Pairing.md](RNA_Base_Pairing.md) and [Minimum_Free_Energy.md](Minimum_Free_Energy.md).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Best-effort parse of malformed input | Assumption | a stray closer is silently dropped rather than throwing | accepted | callers should `ValidateDotBracket` first |
| 2 | null/empty = valid empty structure | Assumption | no exception on null/empty | accepted | empty string is unambiguously balanced |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `""` / null | valid; no pairs | empty balanced structure |
| `".....".` (all unpaired) | valid; no pairs | dots are unpaired [1] |
| `([)]` | parse {(0,2),(1,3)}; valid | crossing families matched independently [1][3] |
| `(]` | invalid | mismatched families; partners must match up [3][4] |
| `())` | parse yields {(0,1)} only | best-effort; stray closer dropped (§5.4) |

### 6.2 Limitations

Does not classify loop types or convert to energy models; only decodes pairing. Letter recognition relies on `char.IsLetter`; non-ASCII letters would be treated as pairing symbols, which is outside the WUSS A–Z convention.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var pairs = RnaSecondaryStructure.ParseDotBracket("((((....))))").ToList();
// pairs == { (0,11), (1,10), (2,9), (3,8) }

bool ok = RnaSecondaryStructure.ValidateDotBracket("([)]"); // true (crossing families)
bool bad = RnaSecondaryStructure.ValidateDotBracket("(]");  // false (mismatched families)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_ParseDotBracket_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_ParseDotBracket_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [RNA-DOTBRACKET-001-Evidence.md](../../../docs/Evidence/RNA-DOTBRACKET-001-Evidence.md)
- Related algorithms: [Pseudoknot_Detection.md](Pseudoknot_Detection.md), [RNA_Base_Pairing.md](RNA_Base_Pairing.md)

## 8. References

1. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26. RNA Structure Notations: https://viennarna.readthedocs.io/en/latest/io/rna_structures.html
2. ViennaRNA Package. Dot-Bracket Notation of Secondary Structures. https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/dotbracket.html
3. ViennaRNA Package. Washington University Secondary Structure (WUSS) notation. https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/wuss.html
4. Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. https://doi.org/10.1093/bioinformatics/btt509
5. Rfam Documentation. Glossary (WUSS format). https://docs.rfam.org/en/latest/glossary.html
