# bio-assembly — fuller recipes & parameter guidance

Progressive-disclosure detail for `SKILL.md`. Rigor rules (parse-with-a-tool, envelope, provenance,
cross-check, 0-based coords, alpha caveat) come from **[`bio-rigor`](../../bio-rigor/SKILL.md)** — not
repeated here.

> **Server note.** The assembly **engine** tools (`assemble_de_bruijn`, `assemble_olc`,
> `assembly_stats`, `calculate_coverage`, `find_overlap`, `find_all_overlaps`, `merge_contigs`,
> `scaffold_contigs`, `error_correct_reads`, `quality_trim_reads`, `compute_consensus`) live on the
> **Alignment** server as `SequenceAssembler.*`. This skill drives them; `bio-alignment` owns their
> reference docs. The **k-mer** (`KmerAnalyzer.*`) and **repeat/complexity**
> (`RepeatFinder.*`, `SequenceComplexity.*`, `GenomicAnalyzer.*`) tools are the Analysis + Core subset
> this skill owns.

## Assembly parameters (de Bruijn / OLC)

Both assemblers accept the same four knobs, but each path uses a different operative one:

| Param | Default | de Bruijn | OLC |
|---|---|---|---|
| `kmerSize` | `31` | **operative** — graph is built on (k−1)-mers; lower for short/low-coverage reads, raise to resolve short repeats | ignored (shared param) |
| `minOverlap` | `20` | shared param | **operative** — minimum suffix-prefix overlap (bp) to accept an edge |
| `minIdentity` | `0.9` | shared param | **operative** — minimum overlap identity ratio (0–1) |
| `minContigLength` | `100` | drop contigs shorter than this (bp) | same |

Guidance:
- **Many short reads / dense coverage / small `kmerSize`** → de Bruijn. Pick `kmerSize` from the
  `kmer_spectrum`: high enough that genuine k-mers separate from the low-count error peak, low enough
  that the graph stays connected. Odd `k` avoids palindromic k-mer collapse.
- **Few reads with long clear overlaps** → OLC; tune `minOverlap`/`minIdentity` to the read length and
  error rate. Inspect the graph first with `find_all_overlaps`.
- **Always report the knob set in provenance** — N50/contig set are only comparable under the same
  parameters.
- Reads are DNA-oriented (A/C/G/T). Validate the alphabet with a `bio-qc` / parser step first
  (`bio-rigor` rule); assemblers do not clean input for you — use `error_correct_reads` /
  `quality_trim_reads` upstream.

## N50 and assembly stats — pick the right call

- `assembly_stats(contigs, totalReads)` recomputes N50 / longestContig / totalLength from a contig
  set you already have (echoes contigs). Use it to **independently verify** the `n50` an assembler
  returned, or to score externally-produced contigs.
- N50 = the contig length *L* such that contigs ≥ *L* cover ≥ 50 % of `totalLength` (sorted
  descending, cumulative crosses half). `totalReads` only fills the read-accounting fields; it does
  not affect N50. See [`docs/algorithms/Assembly/Assembly_Statistics.md`](../../../../docs/algorithms/Assembly/Assembly_Statistics.md).
- **Cross-check:** engine `n50` (from `assemble_*`) must equal `assembly_stats` `n50` on the same
  contigs; `totalLength` must equal Σ|contig|. A mismatch means the contig set changed between steps.

## Coverage over a reference

`calculate_coverage(reference, reads, minOverlap=20)` maps each read to its best **ungapped** position
requiring ≥ `minOverlap` matching bases, and returns `coverage` — an int array of length
`|reference|` where `coverage[i]` = reads spanning **0-based** position *i*. Derive:
- **mean / median depth** = summary over the array;
- **assembly gaps** = maximal runs where `coverage[i] == 0`;
- **breadth** = fraction of positions with `coverage[i] > 0`.
Reads that don't reach `minOverlap` matches anywhere are dropped (not placed), so low `minOverlap`
increases spurious placements — report the value used. Details:
[`docs/algorithms/Assembly/Coverage_Calculation.md`](../../../../docs/algorithms/Assembly/Coverage_Calculation.md).

## k-mer QC for assembly

- **`kmer_spectrum(sequence, k)`** → frequency-of-frequencies map (occurrence count → # distinct
  k-mers). The **low-count peak** (usually count 1–2) is sequencing error; the **main peak** sits near
  the true coverage depth; a **second peak at ~2× / 3×** flags repeats or ploidy. Use it to choose
  `kmerSize` and an error cutoff before assembling. Empty when `k` > |sequence|.
- **`count_kmers` / `count_kmers_both_strands`** → exact per-k-mer counts (single vs canonical
  both-strand). Both-strand is the correct total for a double-stranded genome.
- **`most_frequent_kmers(sequence, k[, top])`** and **`kmers_with_min_count(sequence, k, minCount)`** →
  over-represented k-mers = repeats, adapters, contamination. Screen these before trusting a contig.
- **`unique_kmers(sequence, k)`** → k-mers occurring exactly once (candidate unique anchors / markers).
- **`find_clumps(sequence, k, windowSize, minOccurrences)`** → k-mers that occur ≥ `minOccurrences`
  times within some window of size `windowSize` (`windowSize` ≥ `k`); locates ori boxes, satellite
  arrays, and other local repeat structure. Set unordered.
- **`kmer_frequencies`, `kmer_positions`, `generate_all_kmers`, `analyze_kmers`, `kmer_distance`** —
  supporting k-mer profiling (positions, full frequency table, alphabet enumeration, composite summary,
  k-mer-profile distance between two sequences). Background:
  [`docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md`](../../../../docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md),
  [`docs/algorithms/K-mer/K-mer_Statistics.md`](../../../../docs/algorithms/K-mer/K-mer_Statistics.md).

## Repeat & low-complexity screen (assembly-relevant)

Repeats fragment assemblies and inflate coverage; screen inputs and contigs:
- **`find_tandem_repeats` / `find_microsatellites`** — tandem arrays / SSRs that collapse in graphs.
- **`find_direct_repeats` / `find_inverted_repeats` / `find_palindromes`** — dispersed / inverted
  structures that create spurious branches.
- **`tandem_repeat_summary`** — compact per-sequence tandem-repeat report.
- **`find_low_complexity_regions` + `dust_score`** — mask before assembly with `mask_low_complexity`
  to stop low-complexity reads from tangling the graph.
- **Independent cross-check** of the single biggest repeat: `find_longest_repeat` (DNA, with position)
  or `suffix_tree_lrs` (raw text, suffix-tree path) on **Core** — a second code path corroborating the
  RepeatFinder result. Background:
  [`docs/algorithms/Repeat_Analysis/Repeat_Detection.md`](../../../../docs/algorithms/Repeat_Analysis/Repeat_Detection.md),
  [`Tandem_Repeat_Detection.md`](../../../../docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md).

## Cross-checking playbook (satisfy `bio-rigor` cross-check rule)

| Result | Independent corroboration |
|---|---|
| assembler `n50` / `totalLength` | recompute with `assembly_stats`; `totalLength` == Σ|contig| |
| a contig's k-mer basis | `kmer_spectrum` — no anomalous high-count k-mer for a unique contig |
| coverage breadth / gaps | `coverage` zero-runs vs the contig boundaries |
| largest repeat (RepeatFinder) | `find_longest_repeat` / `suffix_tree_lrs` (Core, suffix-tree) |
| over-represented k-mer set | re-count with `count_kmers_both_strands`; check strand symmetry |

## Coordinates & units (report explicitly)

- `coverage` is indexed by **0-based** reference position; repeat/clump positions are **0-based** as
  documented per tool. State the base in output.
- Lengths, N50, `totalLength`, `longestContig` are in **bp**. `kmer_spectrum` keys are occurrence
  counts (integers), values are counts of distinct k-mers — not lengths.
- `minIdentity` is a fraction in [0,1]; do not confuse with a percentage.

## Envelope note

None of the 9 `LimitationPolicy`-guarded units cover assembly, k-mer, or repeat/complexity tools, so
no `MinimumMode` gate applies to a normal assemble / coverage / N50 / k-mer / repeat call. If a task
chains into a guarded unit (e.g. downstream chromosome binning `META-BIN-001`, or an ONCO/mirna path
in another skill), `bio-rigor`'s envelope rule takes over — check
[`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) and **STOP** on a
`SeqeronLimitationException` rather than forcing output.
