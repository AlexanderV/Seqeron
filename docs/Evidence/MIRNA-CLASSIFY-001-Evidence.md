# Evidence Artifact: MIRNA-CLASSIFY-001 (MiPred, WP-R2-1)

**Test Unit ID:** MIRNA-CLASSIFY-001
**Algorithm:** Pre-miRNA classification — published MiPred random forest (opt-in `ClassifyPreMiRnaMiPred`)
**Date Collected:** 2026-11-25

---

## Primary sources opened

### MiPred package (original model and feature script)

**URLs (public mirror of the web-server package):**
`https://raw.githubusercontent.com/PeterScott/mirna-prediction-nodejs/master/progs/MiPred/model.RData` (3,043,828 bytes, SHA-256 `51dda943d61599de19663eafc063c3ebe07f6c3e039c9cca9638d5485d53e350`),
`.../microRNAcheck_parallel.pl` (17,996 bytes, SHA-256 `5b8fd4af9cfc6ee334479607a414cde01d637e12a3cf361a161c8532d1c0bb75`),
`.../Note.txt` ("assumes the RandomForest (R package) and RNAFold (ViennaRNA package) are installed").
**Accessed:** 2026-11-25 (neutral User-Agent; the GitHub API directory listing was refused by the environment, file names as supplied).
**Handling:** untrusted download in its own directory; the Perl was read as text only and never executed; `model.RData` was never loaded by R, it was parsed by a pure-Python R-serialization reader.

**Key extracted points (script):**

1. Clean-up: characters other than `atcgATCGUu` removed, upper case, T→U.
2. Fold with RNAfold; MFE read from the printed value.
3. Hairpin pre-judgement: `MFE <= -15 && pairs >= 18 && no multiloop && 51 <= length <= 137`; pairs = number of `(` unless the structure matches `\(\.+\)[\.,\(,\)]*\(\.+\)`; multiloop = a `(` after the first `)`.
4. Features (34): for positions 2..n-1 the centre nucleotide (A,U,G,C) × local structure (`(((`,`((.`,`(.(`,`(..`,`.((`,`.(.`,`..(`,`...`) counted, divided by the sequence length (closing brackets never counted; the `tr/)/(/` line is commented out); then MFE; then the p-value.
5. p-value: 999 randfold dinucleotide shuffles (triplet-swap, `rand`), folded with RNAfold; `(cpt_ega + cpt_inf)/(N + 1)` with `cpt_ega` starting at 1.
6. Prediction: `randomForest` `predict(model, test, type = "prob")[, 2] > 0.5` ⇒ real; confidence = 100·p (real) or 100 − 100·p.

**Key extracted points (model.RData):** randomForest classification forest, 1000 trees, mtry 3, classes N (168 pseudo) / P (163 human pre-miRNAs), 34 numeric variables V1..V34, nrnodes 103, 56,324 nodes (28,662 terminal), OOB error 0.0332 (N 0.0238, P 0.0429), formula `V35 ~ V1 + … + V34`, stored `train_data` (331 × 35).

### Jiang P, Wu H, Wang W, Ma W, Sun X, Lu Z. 2007. MiPred: classification of real and pseudo microRNA precursors using random forest prediction model with combined features. Nucleic Acids Res 35(Web Server issue):W339–W344

**URL:** doi:10.1093/nar/gkm368 (web search snippets only; full text not fetched)
Confirms: random forest on combined features (local contiguous structure-sequence composition, MFE, randomization-test p-value from 1000 dinucleotide shuffles), trained on the Triplet-SVM data set; MFE and p-value the most important features (consistent with the stored `importance`: V34 40.8, V33 23.1, V9 10.7, V32 10.4).
Note: Ng & Mishra 2007 (Bioinformatics 23:1321) is a different tool (miR-abela); the task text attributing MiPred to it is a mis-attribution.

### R randomForest 4.7-1.1 (local, trusted oracle)

`predict(rf, newdata, type = "vote", norm.votes = FALSE)` on forests trained locally (set.seed(20261125), ntree 60) — descent rule left iff `x[var] <= split point`; `getTree(rf, k, labelVar = FALSE)` export.

### ViennaRNA 2.8 (Python `RNA.fold`), numpy

Structures/MFE of the oracle hairpins and the independent re-derivation of the features/p-value (port of the Perl, numpy RNG).

## Test datasets

1. Oracle forest: 60 trees, 560 rows, 4 exact 30:30 ties, 160 rows with features equal to split points (MiPredOracleData.cs).
2. Reference hairpins: hsa-mir-21 (MI0000077), hsa-let-7a-1 (MI0000060), hsa-mir-16-1 (MI0000070), hsa-mir-155 (MI0000681), hsa-mir-30a (MI0000088) precursors; three shuffles of let-7a-1; three pri-miRNAs with flanks; 60 stored training rows (MiPredReferenceData.cs).

## Assumptions and documented differences

- Fold: canonical Turner-2004 MFE engine instead of RNAfold 1.x (Turner 1999); p-value feature uses `DinucleotideShuffle` (uShuffle-exact) with a fixed seed instead of randfold's triplet-swap shuffle with `rand` — statistical equivalence tested, draws not reproducible.
- Licence: neither the package nor the mirror states one; only numeric model parameters (split variables/points/classes) are bundled; flagged for the owner.
