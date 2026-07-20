---
type: gotcha
title: "predict_replication_origin assumes one circular bacterial chromosome — min-skew = origin is a heuristic"
tags: [sequence-composition, gotcha]
mcp_tools:
  - predict_replication_origin
  - cumulative_gc_skew
sources:
  - docs/algorithms/Sequence_Composition/Replication_Origin_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_replication_origin assumes one circular bacterial chromosome — min-skew = origin is a heuristic

**The trap.** The origin/terminus predictor maps the **cumulative GC-skew minimum directly to the
replication origin** and the maximum to the terminus. That mapping assumes a **single circular
chromosome with symmetric bidirectional replication**, and it does **not** correct for genome
rearrangements or horizontal gene transfer. It is the one *Simplified* step here.

**Why it bites.** Run it on a **linear** chromosome, a plasmid, a **draft assembly of many
contigs**, a **eukaryotic** chromosome (many origins), or a genome with large inversions / HGT, and
the min/max-skew coordinate is **not a biological origin** — you still get back a confident position,
so the wrongness is silent. The underlying GC-skew math ([[gc-skew]]) is fully sourced and correct;
only the origin/terminus *interpretation* is the heuristic.

**What to rely on instead.** Use `gc_skew` / `cumulative_gc_skew` as a descriptive strand-asymmetry
signal on any sequence. Treat `predict_replication_origin` as valid **only for a single circular
bacterial chromosome**, and confirm the coordinate with independent evidence (DnaA boxes,
Ori-Finder, experimental data) before relying on it. Full model:
[[replication-origin-cumulative-skew]]; the skew primitive is [[gc-skew]].
