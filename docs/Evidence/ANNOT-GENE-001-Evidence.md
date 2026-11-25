# Evidence Artifact: ANNOT-GENE-001

**Test Unit ID:** ANNOT-GENE-001
**Algorithm:** Gene prediction / ribosome binding sites (`GenomeAnnotator.PredictGenes`, `FindRibosomeBindingSites*`, `TrainRbsModel`, `ClassifyRbsStart`)
**Date Collected:** 2026-11-25 (review-2026-09 B11, F4–F6 and audit item A2 / F44)

---

## Online sources (what was actually opened)

### Hyatt D et al. (2010) Prodigal — BMC Bioinformatics 11:119
Algorithm description (nodes = candidate start/stop codons; log-odds scores for coding, RBS and start type; single-genome
self-training; dynamic programming over the node graph). Primary reference of the implemented method.

### Prodigal source V2.6.3 (specification of formulas and constants only)
`raw.githubusercontent.com/hyattpd/Prodigal/GoogleImport/{node.c,node.h,sequence.c,dprog.c,dprog.h,gene.c,training.c,training.h,main.c}`
fetched with a neutral User-Agent. Prodigal is GPL-3: the C# engine (`ProdigalStartModel.cs`) is a clean-room
re-implementation of the published algorithm using the source only as the specification of facts — bin definitions
(`shine_dalgarno_exact/_mm`), scoring formulas (`score_nodes`, `raw_coding_score`, `rbs_score`, `train_starts_sd`,
`calc_dicodon_gene`, `intergenic_mod`, `dprog`/`score_connection`, `tweak_final_starts`) and numeric constants
(st_wt 4.35, ±4/±5 clamps, 10 training iterations, sthresh 35, windows 6–20 nt / 15 nt spacer, MIN_GENE 90 …). No code text
and no trained file is bundled; the trained model is always estimated from the caller's own genome.

### Oracle genomes (raw.githubusercontent.com/valburkin/GenLabs.Genome, GenBank records, no usage restriction per NCBI policy)
* `escherichia_coli.txt` — E. coli K-12 MG1655, GenBank U00096.2, 4,639,675 bp
* `salmonella_enterica.txt` — S. enterica Typhi CT18, GenBank AL513382.1, 4,809,037 bp 

(E. coli file sha256 `b1d61ce0fac63311a301966a65d052c8061b6747afc537f879192027f14308f1`; S. enterica `8545234d1b00441965b1f546d6d8d1e14c480b51375081a25bdf9733da1f088a`.) NCBI/Ensembl were not needed.

### Oracle binary
`apt-get install prodigal` → Prodigal V2.6.3 (Ubuntu 1:2.6.3-6). Runs: `prodigal -i g.fa -t g.trn` (single-genome training,
binary training file), `prodigal -i g.fa -t g.trn -f gff -o g.gff -s g.starts` (genes + per-node start-score file).
The `.trn` file is the raw `struct _training` (doubles); `gc` @0, `st_wt` @16, `bias[3]` @24, `type_wt[3]` @48,
`uses_sd` @72, `rbs_wt[28]` @80 were read with Python `struct`.

---

## Stage A — definition

Previous behaviour (B11 F4–F6, kept as the default): RBS score = motif length / 6; start = most upstream ATG/GTG/TTG per
stop. The blocker recorded in earlier rounds ("needs genome-trained weights") is not proof of impossibility: Prodigal trains
the weights self-supervised on the input genome (single-genome mode, ≥ 20 kb), so the whole method can be implemented.

## Implemented (opt-in; defaults unchanged)

| API | Behaviour |
|---|---|
| `GenomeAnnotator.TrainRbsModel(genome)` | Prodigal single-genome training: nodes → GC-frame bias DP → dicodon model → 10 iterations of `train_starts_sd` → `RbsModel` (28 RBS-bin + 3 start-type log-odds, upstream composition, coding model, SD-usage decision). |
| `TrainRbsModel(genome, genes)` | Same estimator with a caller-supplied gene set as the "real" starts (no coding re-choice). |
| `ClassifyRbsStart(seq, start, strand, model?)` | Prodigal `rbs_score` bins (best exact + one-mismatch motif in the −20…−6 window → bin). |
| `FindRibosomeBindingSites*(…, RbsScoreMode.LogOdds, model?)` | Hit score = weight of the hit's motif/spacer bin; default `MotifLength` (= len/6) unchanged. |
| `PredictGenes(…, StartSelection.ProdigalScore, model?)` | Per stop, argmax of Prodigal start score (coding + RBS + type + upstream). |
| `PredictGenes(…, StartSelection.ProdigalDynamicProgramming, model?)` | Complete Prodigal gene model (DP + `eliminate_bad_genes` + `tweak_final_starts`). |

Not ported: Prodigal's non-SD motif trainer (`train_starts_nonsd`, used when `uses_sd` = 0 — e.g. Bacteroides, Mycoplasma-like
genomes) → `NotSupportedException`; metagenomic mode (`-p meta`); multi-record concatenation; masking (`-m`); `-c` closed ends;
translation tables other than 11. These are scope, not impossibility (Prodigal formulas are available).

## Stage B — numerical cross-check against the Prodigal V2.6.3 binary

**1. Training output (`-t`), full genomes.** All 28 `rbs_wt` and 3 `type_wt`: max |diff| ≤ 5e-5 (the 4-decimal print of the comparison
script) on both genomes; `uses_sd` = 1 on both. E. coli: type_wt = (0.9069, −1.5143, −2.9429); rbs_wt[27] (AGGAGG, spacer 5–10) = 3.269,
rbs_wt[0] (no motif) = −1.9086; S. enterica: type_wt = (0.8284, −1.3498, −2.6635), rbs_wt[27] = 2.9485, rbs_wt[0] = −1.8318.
On the 150 kb E. coli slice used by the unit tests all weights agree to < 1e-9 (locked in the test).

**2. Node scores (`-s` start file), E. coli, all 192,610 start nodes.** coding score (CodPot), start score (StrtSc), RBS score (RBSScr),
upstream score (UpsScr) and type score (TypeScr): 0 nodes differ by more than the file's 2-decimal rounding (max |diff| = 0.005 for
every column); 0 nodes missing.

**3. Gene predictions (`-f gff`), DP mode.** E. coli: 4314/4314 genes identical (begin, end, strand); S. enterica: 4676/4676 identical;
150 kb slice: 129/129 identical.

**4. Start-site agreement of the ORF-based modes** (genes sharing the Prodigal stop; ORFs ≥ 30 aa; model = genome-trained):

| Genome | stop-shared genes | `LongestOrf` identical start | `ProdigalScore` identical start |
|---|---:|---:|---:|
| E. coli K-12 | 4,310 | 2,771 (64.3 %) | 4,246 (98.5 %) |
| S. enterica CT18 | 4,673 | 2,992 (64.0 %) | 4,580 (98.0 %) |
| E. coli 1–150 kb | 128 | 81 (63.3 %) | 126 (98.4 %) |

`ProdigalScore` omits the DP's neighbour (intergenic/overlap) terms and the final adjustment, which explains the residual 1.5–2 %.

**5. Start-type frequency sanity check** (Prodigal-gene-model starts, DP mode). E. coli: ATG 90.7 % (3911), GTG 7.6 % (327), TTG 1.7 % (75), 1 partial gene;
S. enterica: ATG 89.9 %, GTG 7.7 %, TTG 2.4 %. These equal the Prodigal binary's own calls (same gene set). Annotated E. coli K-12 CDS use
≈ 83 % ATG / 14 % GTG / 3 % TTG (Blattner et al. 1997, Science 277:1453 — recalled, not opened); Prodigal systematically over-calls ATG
(its start-type prior learned from the training set is ATG-heavy), so the departure is a property of Prodigal, not of the port.

**6. Performance.** Training 3.2–3.5 s and DP 3.7–3.8 s for a 4.6–4.8 Mb genome (Prodigal binary: 2.6 s total).

## RBS bin formulas (cross-checked by an independent Python port of `shine_dalgarno_exact/_mm` for the unit-test cases)

Window of 6 bases compared to AGGAGG: A at positions 0,3 scores 2; G at 1,2,4,5 scores 3; motif score = −2 + Σ. Exact motifs of length 3–6 with no
mismatch; spacer class from `rdis = start − (motif end)`: `<5` → 2 (motif < 5) / 1 (≥ 5); `11–12` → 1 / 2; `≥ 13` → 3; else 0; `rdis > 15` or score < 6 rejected.
Bins (score, class): (6,2)→1 (6,3)→2 (8/9,3)→3 (6,1)→6 (11/12/14,3)→10 (8/9,2)→11 (8/9,1)→12 (6,0)→13 (8,0)→15 (9,0)→16 (11,2)→20 (11,1)→21 (11,0)→22
(12,2)→20 (12,1)→23 (12,0)→24 (14,2)→25 (14,1)→26 (14,0)→27. One-mismatch 5/6-mers: bins 2–5, 7–9, 14, 17–19. Weight = ln(freq in genes / freq among all candidate starts), clamped ±4.

## Files

Implementation `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/ProdigalStartModel.cs` (+ API in `GenomeAnnotator.cs`); tests
`tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/GenomeAnnotator_ProdigalStart_Tests.cs`; oracle data `TestData/Prodigal/*` (150 kb slice of U00096.2,
Prodigal training weights and gene set); MCP `find_ribosome_binding_sites.scoreMode`, `predict_genes.startSelection`.
