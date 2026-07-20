---
type: gotcha
title: "count_kmers is strand-specific and count_kmers_both_strands is additive — neither is Jellyfish -C canonical collapsing"
tags: [kmer, gotcha]
mcp_tools:
  - count_kmers
  - count_kmers_both_strands
sources:
  - docs/algorithms/K-mer/K-mer_Counting.md
  - docs/algorithms/K-mer/Both_Strand_Kmer_Counting.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# count_kmers is strand-specific and count_kmers_both_strands is additive — neither is Jellyfish -C canonical collapsing

**The trap.** Seqeron's k-mer counters do **not** collapse a k-mer with its reverse complement into
a single canonical key (Jellyfish `-C` / Mash behaviour). `count_kmers` counts the **forward strand
only** — `w` and `RC(w)` are separate keys. `count_kmers_both_strands` uses the **additive** kPAL
"balance" convention (`count(w) = forward[w] + forward[RC(w)]`), which is **still not** the canonical
`occurrences of min(w, RC(w))`. Canonical-collapsing is explicitly **not implemented**.

**Why it bites.** If you compare k-mer spectra, distances, or profiles against a Jellyfish `-C` /
Mash pipeline — the usual reference for strand-agnostic counting — the numbers will **not match**,
and a k-mer table you assume is canonical will double-represent palindrome-asymmetric k-mers. On
double-stranded DNA where strand is arbitrary, `count_kmers` (forward-only) can also under-count the
biological occurrence of a motif that sits on the reverse strand.

**What to rely on instead.** For a strand-agnostic view pick the convention deliberately: use
`count_kmers_both_strands` for the additive kPAL view, or canonicalize as a **post-processing** step
(collapse `w`/`RC(w)` to `min(w, RC(w))`) if you need Jellyfish-`-C` keys. The raw-string / span
surfaces also do **not** enforce a DNA alphabet — validate input first. Full models:
[[k-mer-counting]], [[both-strand-kmer-counting]]; the LCA classifier that *does* canonicalize is
[[taxonomic-classification]].
