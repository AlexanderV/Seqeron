# Evidence Artifact: PROTMOTIF-SP-001

**Test Unit ID:** PROTMOTIF-SP-001
**Algorithm:** Signal Peptide Cleavage-Site Prediction (von Heijne 1986 weight matrix)
**Date Collected:** 2026-06-14

---

## Online Sources

### EMBOSS 6.6.0 `sigcleave` application documentation

**URL:** https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html
**Accessed:** 2026-06-14 (fetched with `curl`; HTML saved and read line-by-line)
**Authority rank:** 3 (reference implementation in an established library)

**Key Extracted Points:**

1. **Method:** "sigcleave predicts the site of cleavage between a signal sequence and the mature exported protein using the method of von Heijne." It is "the EMBOSS implementation of the weight matrix method (von Heijne 1986)."
2. **Eukaryotic weight matrix (verbatim):** The page embeds the default data file "Amino acid counts for 161 Eukaryotic Signal Peptides, from von Heijne (1986), Nucl. Acids. Res. 14:4683-4690", header `# R -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 +1 +2 Expect`, `Sample: 161 aligned sequences`. Full 20×16 table copied into the implementation (see Dataset below).
3. **Cleavage site convention:** "The cleavage site is between +1 and -1."
4. **Acceptance threshold:** "-minweight should be at least 3.5. At this level, the method should correctly identify 95% of signal peptides, and reject 95% of non-signal peptides. The cleavage site should be correctly predicted in 75-80% of cases."
5. **Prokaryote option:** "-prokaryote … changes the default scoring data file name." (eukaryotic is the default).
6. **Worked example (ACH2_DROME):** Output block: `Maximum score 13.739`; hit `(1) Score 13.739 length 13 at residues 29->41`, `Sequence: LLVLLLLCETVQA`, `mature_peptide: NPDAKRLYDD…` — i.e. cleavage between residue 41 and 42, mature protein starts at residue 42.

### EMBOSS 6.6.0 source: `emboss/sigcleave.c`

**URL:** https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/sigcleave.c
**Accessed:** 2026-06-14 (fetched with `curl`; full file read)
**Authority rank:** 3 (reference-implementation source code)

**Key Extracted Points:**

1. **Scoring loop (`main`):** For each candidate index `i`, `weight = 0`; sum `weight += ajFloat2dGet(matrix, aa(j), ic)` over `j = i+pval … i+nval-1` with `ic` starting at `13+pval`, where `pval = -13` and `nval = 2`. The window therefore covers positions −13..+1 of the sequence into matrix columns 0..14 (the +2 / `Expect` column is not used in scoring).
2. **Best site selection:** `if(weight>maxweight){ maxweight=weight; maxsite=i; }` — the single returned prediction is the argmax of the weight over all positions.
3. **Reporting threshold:** `if(weight>minweight)` collects reportable hits; default `minweight = 3.5`.
4. **Cleavage position:** the end of the signal peptide is `maxsite-1` and the mature peptide begins at `maxsite` (`ajStrAssignSubC(&stmp, …, maxsite+pval, maxsite-1)` for the signal window).
5. **Matrix transform (`sigcleave_readSig`):** `expected = mat[i][d2-1]` (the `Expect` column). For each position column `j` of a residue row: if the count is zero it is replaced by `1.0e-10` when `j==10 || j==12` (columns −3 and −1), otherwise by `1.0`; then `mat[i][j] = log(mat[i][j] / expected)` — natural logarithm.

### EMBOSS 6.6.0 data file: `emboss/data/Esig.euk` and `emboss/data/Esig.pro`

**URL:** https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/data/Esig.euk (and `…/Esig.pro`)
**Accessed:** 2026-06-14 (fetched with `curl`; full files read)
**Authority rank:** 3 (reference-implementation data)

**Key Extracted Points:**

1. **Esig.euk:** `Sample: 161 aligned sequences`, "from von Heijne (1986)". The 20×16 count matrix (positions −13..+2 plus `Expect`) was copied verbatim into `EukaryoticCounts`/`EukaryoticExpect`.
2. **Esig.pro:** `Sample: 36 aligned sequences`, "from von Heijne (1986)". The 20×16 count matrix copied verbatim into `ProkaryoticCounts`/`ProkaryoticExpect`. All `Expect` values are > 0 (smallest is W = 0.4).
3. **Column sums:** Each position column sums to the sample size (eukaryotic = 161 column sum verified for several columns; prokaryotic column 0 and 14 sum to 36), matching the EMBOSS sanity check `fabs(sample - rt)`.

### UniProt entry P17644 (ACH2_DROME)

**URL:** https://rest.uniprot.org/uniprotkb/P17644.fasta
**Accessed:** 2026-06-14 (fetched; full FASTA read)
**Authority rank:** 5 (curated sequence database)

**Key Extracted Points:**

1. **Sequence:** The 576-aa sequence used by the EMBOSS worked example. Residues 29–41 are `LLVLLLLCETVQA` and residue 42 is `N` (mature-protein start), matching the EMBOSS output exactly.

### von Heijne G. (1986) — primary publication (metadata only)

**URL:** https://doi.org/10.1093/nar/14.11.4683
**Accessed:** 2026-06-14 (DOI verified via search; full text retrievable only as scanned PDF — the numeric matrix is obtained from the EMBOSS data files above, which reproduce it in text and cite this paper)
**Authority rank:** 1 (primary peer-reviewed paper)

**Key Extracted Points:**

1. **Origin of matrix:** The eukaryotic (161 sequences) and prokaryotic (36 sequences) count matrices used here are explicitly attributed to this paper by the EMBOSS data files.

---

## Documented Corner Cases and Failure Modes

### From EMBOSS `sigcleave` documentation

1. **All sites reported, no intrinsic cutoff:** "sigcleave may predict any number of cleavage sites in a protein sequence but not all of these will be biologically relevant." The `-minweight` threshold (default 3.5) only governs which sites are reported; the maximum-scoring site is always identified.
2. **Accuracy ceiling:** the cleavage site is correctly predicted in only 75–80% of cases (heuristic, not exact).

### From `sigcleave.c` source

1. **Short sequence:** a sequence shorter than the 15-column window still scores (off-window columns are skipped), but a minimum length is required for a meaningful score; this implementation requires at least 15 residues (one full window) and otherwise returns `null`.
2. **Non-standard residues:** characters outside the 20 standard residues are not present in the matrix; in this implementation they contribute 0 to the score (the matrix lookup is skipped), consistent with EMBOSS treating only mapped residues.

---

## Test Datasets

### Dataset: Eukaryotic count matrix (von Heijne 1986, 161 sequences)

**Source:** EMBOSS 6.6.0 `data/Esig.euk`; von Heijne (1986) Nucleic Acids Res. 14:4683-4690.

The full 20×16 integer count matrix (rows A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y; columns −13..+2 then `Expect`) is reproduced verbatim in `ProteinMotifFinder.EukaryoticCounts` / `EukaryoticExpect`. Selected values:

| Residue | −3 | −1 | +1 | Expect |
|---------|----|----|----|--------|
| A | 47 | 80 | 18 | 14.5 |
| L | 8  | 1  | 8  | 12.1 |
| G | 5  | 39 | 10 | 12.1 |

### Dataset: ACH2_DROME worked example (EMBOSS reference output)

**Source:** EMBOSS `sigcleave` documentation; UniProt P17644.

| Parameter | Value |
|-----------|-------|
| Input | UniProt P17644 (576 aa) |
| Maximum score | 13.739 |
| Best signal-peptide window (positions −13..−1) | `LLVLLLLCETVQA` (residues 29–41) |
| Cleavage | between residue 41 and 42 |
| Mature-protein start (CleavagePosition) | 42 |
| Mature peptide prefix | `NPDAKRLYDD…` |

Independent re-derivation of the score (Python, natural-log transform with the EMBOSS pseudocount rule) reproduced **13.739** at mature-start index 42, confirming the formula and matrix.

---

## Assumptions

1. **ASSUMPTION: Minimum input length = one full window (15 aa).** EMBOSS scores any length by skipping off-window columns, but a score over a truncated window is not meaningful for cleavage-site prediction. This implementation returns `null` below 15 residues. This affects only very short inputs (the in-scope eukaryotic/prokaryotic signal peptides are ≥ 15 aa) and does not change the score of any in-window candidate.

---

## Recommendations for Test Coverage

1. **MUST Test:** ACH2_DROME (P17644) returns `CleavagePosition == 42` and `Score == 13.739` (±1e-3) — Evidence: EMBOSS `sigcleave` worked example + independent re-derivation.
2. **MUST Test:** The score equals the sum of `ln(count/expect)` log-odds weights over positions −13..+2 for a hand-computed small window — Evidence: `sigcleave.c` scoring loop + matrix transform.
3. **MUST Test:** Best site is the argmax of the weight (a second known hit, e.g. residues 26→38 score 12.135, is strictly lower) — Evidence: EMBOSS output ordering.
4. **MUST Test:** `IsLikelySignalPeptide` is `true` iff `Score ≥ 3.5` — Evidence: EMBOSS `-minweight` default.
5. **MUST Test:** Null/empty/short (< 15 aa) inputs return `null` — Evidence: window requirement.
6. **SHOULD Test:** Case-insensitivity (upper vs lower produce identical results) — Rationale: implementation upper-cases input.
7. **SHOULD Test:** Prokaryotic matrix selection scores against `Esig.pro` and can differ from the eukaryotic score — Rationale: EMBOSS `-prokaryote` option.
8. **COULD Test:** Non-standard residues (X) contribute 0 and do not throw — Rationale: documented residue handling.

---

## References

1. von Heijne, G. (1986). A new method for predicting signal sequence cleavage sites. Nucleic Acids Research, 14(11):4683–4690. https://doi.org/10.1093/nar/14.11.4683
2. Rice, P., Longden, I., Bleasby, A. (2000). EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics, 16(6):276–277. `sigcleave` application + `data/Esig.euk`, `data/Esig.pro`. https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html
3. UniProt Consortium. ACH2_DROME (P17644). https://rest.uniprot.org/uniprotkb/P17644.fasta (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation. Resolves the prior ⛔ block: the von Heijne (1986) weight matrix and the deterministic single-site selection model were retrieved in text form from the EMBOSS reference implementation (data files + `sigcleave.c`), and the worked example was independently reproduced.
