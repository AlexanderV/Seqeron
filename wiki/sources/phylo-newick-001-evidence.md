---
type: source
title: "Evidence: PHYLO-NEWICK-001 (Newick I/O — tree serialize/parse)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-NEWICK-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-NEWICK-001-Evidence.md
source_commit: 24402115d77d909285d42a30e194bcf6ecfcfa59
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-NEWICK-001

The validation-evidence artifact for test unit **PHYLO-NEWICK-001** — **Newick I/O**
(`PhylogeneticAnalyzer.ToNewick()` / `ParseNewick()`): the **serialization layer** of the
phylogenetics family. It converts the rooted binary `PhyloNode` tree — the one built from an
[[evolutionary-distance-matrix]] (PHYLO-DIST-001) by UPGMA/NJ, resampled by
[[phylogenetic-bootstrap-support]] (PHYLO-BOOT-001), and read/queried by
[[tree-comparison-metrics]] (PHYLO-COMP-001) — to and from the **Newick text format** (nested
parentheses, leaf names, `:`-prefixed branch lengths, internal-node labels, `;` terminator). It is
one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. This is a **format serializer, not a new algorithm** — it warrants no dedicated concept
page; the tree semantics it round-trips live in the PHYLO concept pages above. See
[[test-unit-registry]] for how units are tracked. Research-grade
([[scientific-rigor|research-grade]]), not for clinical use.

## What this file records

- **Authoritative sources:** **Wikipedia — Newick format** (grammar, examples, escaping/label
  restrictions), **PHYLIP — The Newick tree format** (Felsenstein, one of the format authors:
  origin, underscore→blank convention, non-uniqueness), and **Gary Olsen — "Interpretation of
  Newick's 8:45 Tree Format" (1990)** (formal BNF grammar, quoting rules, root branch length/label).
- **Grammar (two consistent forms):** Wikipedia `Tree → Subtree ";"`, `Internal → "(" BranchSet ")"
  Name`, `Branch → Subtree Length`, `Length → empty | ":" number`; Olsen additionally allows the
  **root** to carry `[root_label]` and `[:branch_length]` (the TreeAlign `:0.0` convention). Both
  admit **multifurcation** (`BranchSet → Branch | Branch "," BranchSet`), unnamed nodes (`(,,(,));`),
  and single-taxon trees (`A;`).
- **Label rules:** unquoted labels may not contain blanks, `()`, `[]`, `'`, `:`, `;`, or `,`;
  underscore→blank is the PHYLIP convention.
- **Invariants (N1–N9):** string MUST end with `;` (N1); parsed leaf count = leaf names in string
  (N2); round-trip `ToNewick → ParseNewick` preserves **topology** (N3) and **leaf names** (N4);
  branch lengths are reals after `:` (N5); empty string → throw (N6); null node → empty string (N7,
  defensive); internal names emitted after `)` (N8); numbers use **`.` decimal separator via
  InvariantCulture**, locale-independent (N9). Parsing handles balanced parens (recursive descent),
  scientific notation (`e`/`E`/`+`/`-`), and the Olsen root branch length.

## Deviations and assumptions

No deviations — implementation is compliant on every checked spec row (semicolon termination,
internal-node names, invalid-label suppression per Olsen rules, InvariantCulture number format, root
branch-length parsing, balanced-paren recursive descent, scientific notation). **Documented
scope limitations** (each with a spec reference, all out-of-scope choices rather than bugs):
**binary trees only** (UPGMA/NJ produce bifurcating trees, though the grammar permits N children);
**no quoted `'…'` labels**; **no `[]` comments**; **no underscore→blank** rewrite; float precision
±0.00005 (F4, adequate for UPGMA/NJ; the spec imposes no precision limit). UPGMA/NJ labels
containing Newick metacharacters are omitted on output per the Olsen label rules. No source
contradictions.
