# bio-alignment — fuller recipes & parameter guidance

Progressive-disclosure detail for `SKILL.md`. Rigor rules (parse-with-a-tool, envelope, provenance,
cross-check, 0-based coords, alpha caveat) come from **`bio-rigor`** — not repeated here.

## Scoring parameters (shared by all pairwise + MSA tools)

`global_align`, `local_align`, `semi_global_align`, `multiple_align` take the same four knobs:

| Param | Default | Meaning |
|---|---|---|
| `match` | `+1` | reward for an identical column |
| `mismatch` | `-1` | penalty for a substitution column |
| `gapOpen` | `-2` | penalty applied when a gap starts |
| `gapExtend` | `-1` | penalty **per gap position** (linear gap model) |

Guidance:
- **DNA, closely related seqs:** defaults are fine. Raise `|mismatch|` (e.g. `-2`) to punish
  substitutions harder and favor gaps; lower it to tolerate SNP-dense inputs.
- **Favor few long gaps over many short gaps:** make `gapOpen` more negative relative to
  `gapExtend` (affine-style intent). The model here is linear per position; treat `gapOpen` as the
  one-time cost and `gapExtend` as the per-base cost.
- **Always report the scoring set in provenance** — identity/score are only comparable under the
  same scheme.
- These tools are **DNA-oriented** (A/C/G/T, uppercased internally). Validate the alphabet first
  (`bio-rigor` rule 5). For protein/BLOSUM-style scoring, confirm support via `seqeron-discovery`
  before assuming it — do not silently reuse the DNA match/mismatch scheme on protein.

## Identity denominators — pick deliberately

Three different "identity" answers exist; they are **not interchangeable**:

| Tool | Denominator | Use when |
|---|---|---|
| `alignment_statistics` (`SequenceAligner.CalculateStatistics`) | **full alignment length, incl. gap columns** (EMBOSS needle convention) | reporting identity/similarity/gap % of an alignment you already produced |
| `sequence_identity` (`SequenceAssembler.CalculateIdentity`) | **equal-length, gapless**, position-by-position; returns 0 if lengths differ | two already-equal-length seqs (e.g. two MSA rows, fixed-length barcodes) |
| `calculate_similarity` (`GenomicAnalyzer.CalculateSimilarity`) | **k-mer Jaccard**, not positional at all | fast alignment-free similarity screen |

So the same pair can read 83 % (gapped-length identity) vs a different % (gapless) vs a Jaccard
value — always state which you used.

`alignment_statistics.similarity` = (matches + similar substitutions)/length; for DNA it typically
equals identity unless a similarity/substitution grouping is in play — see
[`docs/algorithms/Alignment/Alignment_Statistics.md`](../../../../docs/algorithms/Alignment/Alignment_Statistics.md).

## Choosing an alignment mode

- **Global (`global_align`)** — two sequences you believe are homologous over their whole length
  (e.g. two alleles of the same gene). Penalizes end gaps; every base is aligned.
- **Local (`local_align`)** — locate a shared motif/domain inside otherwise dissimilar sequences.
  Output is only the high-scoring core with its 0-based `start/end` in each input; flanks dropped.
- **Semi-global (`semi_global_align`)** — free end gaps in seq2. Query (seq1, shorter) is fit into
  reference (seq2, longer); the reference is reproduced in full and `startPosition2..endPosition2`
  tells you where the query landed (0-based). Use for overlap detection, containment, or locating a
  primer/insert inside an amplicon.

## MSA + consensus + conservation

`multiple_align(sequences=[…])` → `alignedSequences` (equal-length rows), `consensus`
(majority vote), `totalScore` (sum-of-pairs). It is **anchor-based progressive** (center picked by
4-mer cosine similarity, suffix-tree anchoring) — good for a family of related DNA seqs; it is not a
guaranteed-optimal MSA. Details: [`docs/algorithms/Alignment/Multiple_Sequence_Alignment.md`](../../../../docs/algorithms/Alignment/Multiple_Sequence_Alignment.md).

Downstream:
- **Independent consensus check:** run `compute_consensus(alignedReads=alignedSequences)` — a
  different code path (Biopython dumb_consensus, threshold 0.5, ties → `N`). It should reproduce the
  `multiple_align` consensus; a mismatch flags a threshold/tie effect worth reporting.
- **Per-column conservation:** columns where the consensus is a definite base (not `N`) and all rows
  agree are conserved; `N` marks tied/ambiguous columns. `compute_consensus` requires **equal-length
  pre-aligned** rows (error 1003 otherwise) — feed it MSA rows, not raw reads.

## Alignment-free / quick comparisons

- **`edit_distance`** — minimum insert/delete/substitute edits; sequences may differ in length.
  Good for "how many changes separate these?". Case-insensitive.
- **`hamming_distance`** — substitutions only; **requires equal length** (error 1004 otherwise).
  Use for fixed-length barcodes / already-aligned columns.
- **`calculate_similarity`** — k-mer Jaccard in [0,1], optional `kmerSize`. Fast screen before
  committing to a full alignment; robust to small rearrangements but ignores position.
- **`find_longest_common_region`** (DNA, returns positions) or **`suffix_tree_lcs`** (raw text) —
  the longest exact shared block; a cheap, independent corroboration of a local alignment's core.

## Approximate pattern search (motif / primer occurrence)

- **`find_with_mismatches`** — Hamming, fixed-length window, up to `maxMismatches`; reports position,
  matched substring, distance, and 0-based mismatch positions.
- **`find_with_edits`** — Levenshtein, variable-length windows (`pattern.Length ± maxEdits`); use when
  indels are possible, not just substitutions.
- **`find_best_match`** — the single best (min-Hamming) window; leftmost on ties, exact match
  short-circuits.
- **`frequent_kmers_with_mismatches`** — most-frequent k-mers allowing d mismatches (neighborhood
  tally); motif discovery.

## Cross-checking playbook (satisfy `bio-rigor` rule 4)

| Result | Independent corroboration |
|---|---|
| local-alignment core region | `find_longest_common_region` positions should overlap it |
| `alignment_statistics` identity | equals `matches / alignmentLength` from the same call |
| MSA consensus | `compute_consensus` on the same rows should match |
| edit distance ≈ 0 | `hamming_distance` (if equal length) or `sequence_identity` ≈ 1 |
| approximate hit set | re-search reverse-complement / re-count with `count_approximate_occurrences` |

## Coordinates & units (report explicitly)

- All `startPosition*` / `endPosition*` and `find_*` positions are **0-based**, end-inclusive as
  documented per tool. State the base in output.
- `alignment_statistics` returns **percentages** (identity/similarity/gapPercent); `sequence_identity`
  and `calculate_similarity` return **fractions in [0,1]** — do not mix the scales.
- Lengths are in **bp**. Scores are unitless and only comparable under the same scoring scheme.

## Envelope note

The Alignment/Core tools used here are not among the 9 `LimitationPolicy`-guarded units, so no
`MinimumMode` gate applies to a normal pairwise/MSA/distance call. If a task chains into a guarded
unit (e.g. downstream popgen/phylo/ONCO), `bio-rigor`'s envelope rule takes over — check
[`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) and stop on a
`SeqeronLimitationException` rather than forcing output.
