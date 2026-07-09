---
type: concept
title: "Snapshot / approval testing (Verify) — golden-master regression capture"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/05_SNAPSHOT_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: snapshot-testing-checklist
      evidence: "Priority P1, framework Verify + VerifyNUnit. Serialises the full output to a committed .verified.txt golden master; any output change requires explicit developer approval — 'ловить ненавмисні регресії у складних структурованих виходах' (catches unintended regressions in complex structured outputs)."
      confidence: high
      status: current
---

# Snapshot / approval testing (Verify)

Snapshot (approval) testing serialises an algorithm's **entire output** and compares it against a
committed `.verified.txt` **golden-master** file. The first run records the golden file; every
later run diffs against it and **fails on any change** — and accepting a new snapshot is an
explicit, reviewed developer action. It is the right tool for **complex structured outputs**
where a regression is easy to introduce and hard to spot with hand-written assertions:
alignments, phylogenetic trees, annotation tables, parse results. Seqeron uses **Verify** with
`VerifyNUnit`. This is a **P1** member of the [[validation-and-testing]] program; the checklist
record is [[snapshot-testing-checklist]].

## The workflow

1. First run creates the `.verified.txt` file (the golden master).
2. Subsequent runs compare and **fail on any diff**.
3. On an intended behaviour change: review the diff and **accept** the new snapshot; on an
   unintended one: it's a regression to fix.

## Coverage — the one genuinely incomplete methodology

Snapshot testing is the outlier of the ten: **37 / 255 complete, 221 not started** at the
2026-03-19 checklist, with ~20 snapshot files under `Snapshots/` covering Alignment, Annotation,
Chromosome, Codon, CRISPR, Disorder, Epigenetics, FileIO, Metagenomics, MiRNA, MolTools,
PatternMatching, Phylogenetic, Population, PrimerProbe, ProteinMotif, Repeat, Restriction, RNA,
and Splicing — five new files still needed (Composition, Kmer, Translation, Oncology). This is
the concrete counter-example to any claim that the testing program is uniformly complete: unlike
[[property-based-testing]], [[metamorphic-testing]], and [[fuzzing]] (all at 258/258), snapshot
coverage is a real, tracked gap — one more entry in the [[research-grade-limitations|deep but
uneven]] caveat.

## Relation to characterization testing

Snapshot and [[characterization-testing]] overlap heavily — both freeze current output as a
golden master. The distinction the checklists draw: snapshots are a **standing regression guard**
run on every build, whereas characterization tests are **on-demand, generated specifically before
a risky refactor** and discarded after. In practice the `Snapshots/` files already do much of the
characterization job.
