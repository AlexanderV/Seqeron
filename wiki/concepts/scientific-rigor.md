---
type: concept
title: "Scientific rigor by construction"
tags: [rigor]
sources:
  - README.md
  - docs/Validation/LIMITATIONS.md
source_commit: fc1c1454fa0b1bb8600eddb59af98471c2ff2ddc
created: 2026-07-09
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: readme
      evidence: "A runtime LimitationPolicy guards each algorithm's validated scope, and results are tool-computed with provenance"
      confidence: medium
      status: current
---

# Scientific rigor by construction

The principle that Seqeron results are trustworthy because the system is built to be honest at runtime — the assistant never fabricates a number. Every value is computed by a real algorithm and carries provenance, and each algorithm is guarded against use outside its validated scope.

## Mechanisms

- **`LimitationPolicy`** — a runtime guard that constrains each algorithm to its validated operating envelope, refusing (or flagging) use outside it. Its three modes (`Strict`/`Moderate`/`Permissive`) and honest scope catalog (`docs/Validation/LIMITATIONS.md`) are written up in [[operating-envelope-and-limitation-policy]].
- **Tool-only computation** — results come from tool/algorithm calls, not model guesses; the `bio-rigor` skill enforces this discipline (tool-only, 0-based coordinates, provenance on every result).
- **Provenance** — outputs are reproducible and cited back to the tool chain that produced them.

## Relationship to other pieces

Rigor is the *runtime* half of Seqeron's correctness story; [[validation-and-testing]] is the *development-time* half (how the algorithms were checked before release). It is also what makes the [[three-front-doors]] equivalence meaningful — the same guarded algorithm answers on every path — and it is delivered to plain-language users through the [[skill-layer]]. The honesty stance is bounded by real caveats: see [[research-grade-limitations]].
