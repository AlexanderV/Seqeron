# Alignment significance statistics (Karlin-Altschul, BLAST+ tables, Spouge E-value, Z-score)

> Audit R14-7 (B09 F300). Unit ALIGN-STATS-001. Code: `AlignmentSignificance` (`AlignmentSignificance.cs`, `.Gapped.cs`, generated `.Tables.cs`) in `Seqeron.Genomics.Alignment`.

## What it computes

For a local alignment score S of a scoring system with integer scores:

- **Bit score** S' = (lambda S - ln K) / ln 2 (NCBI `bitscore`).
- **E-value**: the expected number of chance alignments of score >= S in the search space; **P-value** = 1 - exp(-E) (computed without cancellation: 1.24e-33 stays 1.24e-33).
- **Karlin-Altschul parameters** lambda, K, H: lambda is the positive root of sum_s p_s exp(lambda s) = 1 over the score distribution p_s = sum of p_i q_j over the residue pairs with score s (the safeguarded Newton iteration of NCBI `NlmKarlinLambdaNR`), H = lambda sum s p_s exp(lambda s) (`BlastKarlinLtoH`), K by the series of Karlin & Altschul 1990 (`BlastKarlinLHtoK`; the closed forms for steps of one). Background: Robinson & Robinson 1991 for protein (the NCBI list), uniform A / C / G / T for nucleotides, or any distribution (`FromScoreDistribution`). This is the "ideal" block of BLAST; BLAST+ composition-based statistics (blastp default) change lambda per alignment and are not modelled (the values equal `-comp_based_stats 0`).
- **Gapped parameters**: lambda, K, H, alpha, beta and the Gumbel block (C, alpha_v, sigma) cannot be derived and are tabulated by NCBI: `TryGappedProtein` (BLOSUM45 / 50 / 62 / 90, PAM30 / 70 / 250 and the ungapped rows; the table of BLOSUM80 belongs to the NCBI matrix whose A/A is 5, not to the EMBOSS / Biopython BLOSUM80 with A/A = 7, and `ForMatrix` refuses a matrix whose ungapped lambda differs from the table's) and `TryGappedNucleotide` (blastn reward / penalty 1/-5, 1/-4, 2/-7, 1/-3, 2/-5, 1/-2, 2/-3, 3/-4, 4/-5, 1/-1, 3/-2, 5/-4 with their gap costs; a pair is divided by its gcd; odd scores of the even-only tables are rounded down to an even number for the E-value only). The tables are generated from NCBI `blast_stat.c` by `scripts/alignment/generate_blast_stat_tables.py`.
- **E-value formulas**: for the protein matrices BLAST+ prints Spouge's finite-size-corrected value (Spouge 2002; `BLAST_SpougeStoE`): E = area K exp(-lambda S) with area = p1 p2 + c(y) P(m) P(n), where p = (L - l(y)) P + sqrt(v(y)) phi, l(y) = a y + b, the variances and covariance from alpha_v and sigma (`SpougeSignificance`); for blastn the classic length-adjusted formula E = (m - ell)(n - N ell) K exp(-lambda S) with the length adjustment ell of `BLAST_ComputeLengthAdjustment` (`KarlinAltschulSignificance`, `LengthAdjustment`). `ForMatrix` / `ForNucleotide` choose the formula BLAST+ uses; gap costs are in the BLAST convention (a gap of length k costs existence + k x extension = the library's `GapOpen` / `GapExtend` negated).
- **Z-score** (`ShuffleZScore`, MCP `local_align_zscore`): the score against `shuffles` seeded shuffles of the second sequence (EMBOSS `prss` / FASTA idea): z = (S - mean) / sd.

## Verification

BLAST+ 2.12.0 (`blastp -comp_based_stats 0 -seg no -subject`, `blastn -task blastn -dust no -subject`, `-outfmt 5`; `scripts/alignment/generate_significance_oracle.py`): the bit score (to 0.0006) and the E-value (relative 2e-5) of **2 283 HSPs** equal over 18 protein matrix / gap-cost combinations (1 285 HSPs) and 15 blastn reward / penalty / gap combinations (998 HSPs); the ungapped lambda / K / H computed from the score distribution equal the first row of each NCBI table (BLOSUM62 0.3176 / 0.134 / 0.4012; BLOSUM45, 50, 90, PAM30, 70, 250); the worked example of the NCBI comment (scores -2 / 0 / 3 with 0.7 / 0.1 / 0.2: lambda 0.330, K 0.154); the +-1 walk (lambda = ln 3, K = 1/3).

## MCP

`local_align`, `local_align_affine` and `water_align` take `significance` (+ `dbLength`, `dbSequences`) and return `significance` (bitScore, eValue, pValue, method, lambda, k, h, lengthAdjustment, searchSpace); gap costs outside the NCBI tables are refused with the table's reason. `local_align_zscore`.

## Sources actually opened

NCBI `blast_stat.c` (raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public, User-Agent `Seqeron-review/1.0`): the tables, `NlmKarlinLambdaNR`, `BlastKarlinLtoH`, `BlastKarlinLHtoK`, `BLAST_ComputeLengthAdjustment`, `BLAST_SpougeStoE`, `Blast_GumbelBlkLoadFromTables`, `s_GetNuclValuesArray`; the installed BLAST+ 2.12.0 as the oracle. Karlin S, Altschul SF (1990) PNAS 87:2264; Altschul et al. (2001) NAR 29:351; Spouge JL (2002) (finite-size correction), as cited in the NCBI source.
