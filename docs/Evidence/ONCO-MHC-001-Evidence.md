# Evidence Artifact: ONCO-MHC-001

**Test Unit ID:** ONCO-MHC-001
**Algorithm:** MHC-Peptide Binding (length filtering + affinity/%rank thresholds; matrix-based prediction)
**Date Collected:** 2026-06-14 (classification); 2026-06-25 (matrix-based prediction extension)

---

## Online Sources

### NetMHCpan-4.1 / NetMHCIIpan-4.0 (Reynisson et al. 2020)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7319546/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Nucleic Acids Research*)
**How retrieved:** WebSearch query `NetMHCpan-4.1 Reynisson 2020 strong binder weak binder %rank threshold 0.5 2 IC50` → opened the PMC full text with WebFetch.

**Key Extracted Points:**

1. **Class I %Rank thresholds (verbatim):** "by default, %Rank < 0.5% and %Rank < 2% thresholds are considered for detecting SBs and WBs for class I" (SB = strong binder, WB = weak binder).
2. **Class II %Rank thresholds (verbatim):** "and %Rank < 2% and %Rank < 10%, for SBs and WBs for class II".
3. **Class I peptide length range (verbatim):** "for class I, the length range goes from 8 to 14 amino acids, default is 8–11".
4. **Citation extracted from the page:** Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. *Nucleic Acids Res*. 2020;48(W1):W449–W454. doi:10.1093/nar/gkaa379.

---

### Sette et al. (1994) — class I affinity vs immunogenicity

**URL:** https://pubmed.ncbi.nlm.nih.gov/7527444/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Journal of Immunology*)
**How retrieved:** WebSearch query `Sette 1994 Journal of Immunology relationship affinity MHC class I immunogenicity 500 nM CTL response threshold` → opened the PubMed abstract page with WebFetch.

**Key Extracted Points:**

1. **Affinity threshold (verbatim from abstract):** "an affinity threshold of approximately 500 nM (preferably 50 nM or less) apparently determines the capacity" (of a peptide epitope to elicit a CTL response). This is the primary biological basis for the 50 nM (strong) and 500 nM (binder) IC50 cutoffs.
2. **Citation:** Sette A, Vitiello A, Reherman B, Fowler P, Nayersina R, Kast WM, Melief CJ, Oseroff C, Yuan L, Ruppert J, Sidney J, del Guercio MF, Southwood S, Kubo RT, Chesnut RW, Grey HM, Chisari FV. The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *J Immunol*. 1994 Dec 15;153(12):5586–92. PMID: 7527444.

---

### IEDB — selecting thresholds for MHC binding predictions

**URL:** https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated bioinformatics resource; IEDB)
**How retrieved:** WebSearch query `IEDB MHC class I binding prediction strong binder IC50 < 50 nM intermediate < 500 nM threshold` returned this IEDB help article; the article body is served behind a 403 to WebFetch, so the verbatim statement below is the text returned in the IEDB search-result snippet (the direct page fetch returned HTTP 403 — recorded here for auditability).

**Key Extracted Points:**

1. **IC50 affinity tiers (verbatim from snippet):** "Peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity."
2. **Binder cutoff (verbatim from snippet):** "an absolute binding affinity (IC50) threshold of 500 nM identifies strong binders".

---

### Roomp, Antes & Lengauer (2010) — corroborating the 500 nM binder cutoff

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2836306/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *BMC Bioinformatics*)
**How retrieved:** WebSearch (IEDB threshold query, above) listed this PMC article; opened it with WebFetch (after following the ncbi→pmc redirect) to corroborate the IC50 = 500 nM binder demarcation independently of IEDB.

**Key Extracted Points:**

1. **500 nM binder demarcation (verbatim):** "Any peptides annotated in IEDB as binders with IC50 values greater than 500 nM, and peptides annotated as non-binders with IC50 values less than 500 nM were discarded." — confirms 500 nM as the binder/non-binder demarcation in a peer-reviewed source.
2. **Citation:** Roomp K, Antes I, Lengauer T. Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics*. 2010 Feb 17;11:90. doi:10.1186/1471-2105-11-90.

---

### IEDB — MHC class II binding tool description (peptide length range)

**URL:** https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated bioinformatics resource; IEDB)
**How retrieved:** WebSearch query `IEDB MHC class II peptide length 13-25 binding prediction core 9-mer length range` returned this IEDB article; page body is 403 to WebFetch, so the verbatim statement below is from the search-result snippet (direct fetch returned HTTP 403 — recorded for auditability).

**Key Extracted Points:**

1. **Class II length range (verbatim from snippet):** "Peptides binding to MHC class II molecules can vary considerably, and typically range between 13 and 25 amino acids long."
2. **Binding core (verbatim from snippet):** "The peptide binding core ... is usually nine amino acids long (9-mer)."

---

### BIMAS — HLA peptide motif search scoring documentation (matrix-based prediction)

**URL:** http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html (server retired)
**Accessed:** 2026-06-25
**Authority rank:** 5 (curated NIH bioinformatics tool documentation, restating Parker 1994)
**How retrieved:** server `www-bimas.cit.nih.gov` no longer resolves (curl exit 6 / ECONNREFUSED). Retrieved via the Internet Archive with `curl http://web.archive.org/web/20041016022153id_/http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html` on 2026-06-25, then stripped HTML.

**Key Extracted Points:**

1. **Score meaning (HLA-A2, verbatim):** "In the case of HLA-A2, this score corresponds to the estimated **half-time of dissociation** of complexes containing the peptide at 37 oC at pH 6.5."
2. **Scoring rule (verbatim):** "The initial (running) score is set to 1.0. … The running score is then multiplied by the coefficient for that amino acid type, at that position … Using 9-mers, nine multiplications are performed … The resulting running score is **multiplied by a final constant** to yield an estimate of the half time of disassociation."
3. **Matrix shape (verbatim):** "the 181 values (**180 coefficient values plus one final constant**)" — a 20×9 coefficient matrix per allele plus a final constant.
4. **Unlisted / ambiguous residue (verbatim):** ambiguous chars "given a coefficient of **1.00** … the effect … is to leave the score unchanged" → neutral multiplicative coefficient = 1.0.
5. **Independent-binding model (verbatim):** "each amino acid in the peptide contributes independently to binding … Highly favorable amino acids have coefficients substantially greater than 1, and unfavorable amino acids have positive coefficients that are less than one."
6. **A2 table is published (verbatim):** "the coefficient tables used at this Web site have not been published elsewhere (**except for HLA-A2**)" — the HLA-A2 coefficients are those of Parker 1994.
7. **Coefficient values not embeddable:** the per-allele values are served only by a dynamic CGI (`/cgi-bin/molbio/hla_coefficient_viewing_page`); the Internet Archive captured the *input form* (the molecule menu lists `A_0201_standard`, `A1_standard`, …) but **not** the generated value table. The A_0201_standard numeric coefficients could not be retrieved.

---

### Parker, Bednarek & Coligan (1994) — primary paper (HLA-A2 coefficient table)

**URL:** https://pubmed.ncbi.nlm.nih.gov/8254189/
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed primary paper, *Journal of Immunology*)
**How retrieved:** WebSearch `Parker Bednarek Coligan 1994 "Scheme for ranking potential HLA-A2 binding peptides"` → opened the PubMed abstract with WebFetch. Full text is paywalled (journals.aai.org → academic.oup.com; no free PMC/PDF), so the 180-value table itself could not be retrieved.

**Key Extracted Points:**

1. **Matrix (verbatim abstract):** "the binding data from a set of 154 peptides were combined together to generate a table containing **180 coefficients (20 amino acids x 9 positions)**."
2. **Scoring rule (verbatim abstract):** the "'theoretical' binding stability (**calculated by multiplying together the corresponding coefficients**) matched the experimental binding stability to within a factor of 5." Confirms the product-of-coefficients rule; accuracy ≈ factor of 5.
3. **Half-life basis (verbatim):** binding assayed by "the rate of dissociation of beta 2m"; non/weak-binders had "a half-life of beta 2m dissociation of less than 5 min at 37 degrees C."

---

### SMM score ↔ IC50 linearisation (Peters & Sette 2005 / IEDB log50k)

**URL:** https://dmnfarrell.github.io/bioinformatics/mhclearning ; primary: https://doi.org/10.1186/1471-2105-6-132 (Peters & Sette 2005)
**Accessed:** 2026-06-25
**Authority rank:** 1 (primary: Peters & Sette 2005) / 3 (reference-tutorial restating the IEDB transform)
**How retrieved:** WebSearch `IEDB MHC binding "50000" log50k ic50 transform` → opened the dmnfarrell tutorial with WebFetch for the verbatim formula; Peters & Sette 2005 confirmed via WebSearch + Springer (full text login-gated, PMC reCAPTCHA-gated); GitHub `ykimbiology/smm` LICENSE inspected with `gh api`.

**Key Extracted Points:**

1. **Linearisation (verbatim):** "`log50k = 1-log(ic50) / log(50000)`". Inverting: log(IC50) = (1 − score)·log(50000) ⇒ **`IC50 = 50000^(1 − score)`**.
2. **Additive SMM matrix:** SMM prediction is a sum of position-specific contributions plus an intercept; fitted on log-transformed IC50 (Peters & Sette 2005).
3. **Derived anchor points (exact):** score 0 → 50000 nM; score 1 → 1 nM; score 0.5 → √50000 = 223.6067977499790 nM.
4. **Not redistributable:** GitHub `ykimbiology/smm` ships a **0-byte LICENSE** (all-rights-reserved); IEDB SMM matrices are non-commercial. The matrix is caller-supplied (like CIBERSORT LM22 in ONCO-IMMUNE-001).

---

### MHCflurry 2.0 — Class I pan-allele binding-affinity neural network (O'Donnell et al. 2020; Apache-2.0 source + weights)

**URL:** https://github.com/openvax/mhcflurry (source, Apache-2.0) ; primary paper https://doi.org/10.1016/j.cels.2020.06.010 ; licence https://github.com/openvax/mhcflurry/blob/master/LICENSE
**Accessed:** 2026-06-25
**Authority rank:** 1 (primary paper) / 3 (reference implementation: the MHCflurry source actually installed — `mhcflurry` 2.1.5 — and the `models_class1_pan` release 20200610 read in this session)

**How retrieved:** `pip3 install --user mhcflurry` (version 2.1.5) and `mhcflurry-downloads fetch models_class1_pan` to obtain the Apache-2.0 source + trained weights; the installed module files were read directly (`amino_acid.py`, `encodable_sequences.py`, `allele_encoding.py`, `class1_neural_network.py`, `regression_target.py`, `ensemble_centrality.py`, `class1_affinity_predictor.py`); the model `manifest.csv` (`network_json` Keras graphs), the `weights_*.npz` arrays, and `allele_sequences.csv` were inspected with numpy/pandas; oracle IC50 values were produced by re-running the `Class1AffinityPredictor` Python API in-session; the `LICENSE` (`mhcflurry-2.1.5.dist-info/licenses/LICENSE`) was opened and confirmed to be the full **Apache License 2.0**.

**Key Extracted Points:**

1. **Amino-acid alphabet + BLOSUM62 (verbatim, `amino_acid.py`):** `COMMON_AMINO_ACIDS` sorted, then `"X"` appended → order **`ACDEFGHIKLMNPQRSTVWYX`** (21 symbols). `BLOSUM62_MATRIX` is the standard BLOSUM62 with an added all-zero `X` row/column except `X,X = 1`. Both peptide and allele use the `"BLOSUM62"` vector encoding (width 21).
2. **Peptide encoding (verbatim, `encodable_sequences.py`):** the pan-allele models use `alignment_method="left_pad_centered_right_pad"`, `max_length=15`. The peptide is placed three times in a `3·max_length = 45`-position layout — left-aligned (`result[:length]`), centred (`offset = max_length + floor((max_length−length)/2)`), and right-aligned (`result[-length:]`) — with the `X` index filling the rest; then BLOSUM62-encoded → 45×21 = **945** values. Minimum length 5, maximum 15.
3. **Allele pseudosequence (verbatim, `allele_encoding.py` + `allele_sequences.csv`):** each allele maps to a fixed **37-residue** MHC-pocket pseudosequence (e.g. `HLA-A*02:01 → YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA`), BLOSUM62-encoded and flattened → 37×21 = **777** values. (14 993 alleles; 11 609 are `HLA-`.)
4. **Network (verbatim, `class1_neural_network.make_network`):** inputs `peptide` (45×21) and `allele` (index → `Embedding` holding the pseudosequence representation). Both are `Flatten`ed and `concatenate`d (peptide-flat then allele-flat → **1722**). Then `layer_sizes` `Dense` layers with `activation="tanh"`, optional `Dropout` (inactive at inference), and a final `Dense(1, activation="sigmoid")`. Two topologies appear in the ensemble: `feedforward` (plain stack) and `with-skip-connections` (each hidden layer i≥1 receives `concatenate(densenet_layers[-2:])`, i.e. `concat(prev_prev_input, prev_activation)`); the output layer has no skip.
5. **Output → IC50 (verbatim, `regression_target.to_ic50`):** `to_ic50(x) = max_ic50 ** (1 − x)` with `max_ic50 = 50000.0`.
6. **Ensemble combiner (verbatim, `class1_affinity_predictor.predict_to_dataframe` + `ensemble_centrality.py`):** default `centrality_measure="mean"`; `logs = log(per-model ic50)`, `prediction = exp(mean(logs))` — i.e. the **geometric mean** of the per-network IC50s.
7. **`models_class1_pan` 20200610 ensemble:** 10 networks; topologies/`layer_sizes`: `[1024,512]`, `[256,512,512]`*, `[1024,1024]`, `[256,512]`*, `[1024,512]`, `[1024,1024]`, `[256,256,512]`*, `[256,512,512]`*, `[512,512]`, `[1024,512]` (`*` = with-skip-connections). Total ≈ 20.0 M parameters ≈ 80 MB of float32 weights.
8. **Apache-2.0 licence (verbatim header):** "`Apache License / Version 2.0, January 2004`". MHCflurry source and the bundled `models_class1_pan` weights/pseudosequences are redistributable with attribution; the `NOTICE`/attribution is preserved in `Resources/MHCFLURRY_NOTICE.txt`.

---

## Documented Corner Cases and Failure Modes

### From MHCflurry source

1. **Peptide length out of [5, 15]:** `EncodableSequences` raises `EncodingError`; the C# port throws `ArgumentOutOfRangeException`.
2. **Non-canonical residue:** mapped to the `X` index (`allow_unsupported_amino_acids` / the `X` fallback), encoded as the all-but-`X`-zero vector.
3. **Allele not in the table:** the pan-allele predictor needs a known pseudosequence; the port throws `KeyNotFoundException` for an unknown allele.
4. **Two topologies in one ensemble:** a forward pass that ignores the `with-skip-connections` wiring produces wrong inputs to the 2nd/3rd dense layer (dimension mismatch) — the port reads `topology` per network and wires the skips accordingly.
5. **Embedding weights are external:** the `weights_*.npz` `array_0` embedding entry is shape `(0, 777)` (empty) — the allele representation is supplied at predict time from `allele_sequences.csv`, not from the npz; the port drops the empty embedding array and supplies the pseudosequence via `EncodePseudosequence`.

### From BIMAS scoring documentation / SMM

1. **Unlisted residue:** BIMAS product → coefficient 1.0 (no effect); SMM sum → contribution 0.0 (no effect).
2. **Peptide length must equal the matrix position count:** the loader stores one row per position; `Predict*` throws `ArgumentException` if `peptide.Length ≠ Rows.Count`.
3. **SMM IC50 is always finite & > 0** for any finite score (50000^x > 0), so it always satisfies the `ClassifyBindingAffinity` precondition.
4. **No embeddable trained matrix:** BIMAS coefficient files (dynamic CGI, unarchived) and Parker 1994's table (paywalled) were not retrievable; IEDB SMM is non-commercial → matrix is caller-supplied (Framework status).

### From Reynisson et al. (2020)

1. **Class I length out of range:** lengths below 8 or above 14 are not valid class I peptides; the canonical neoantigen search default is 8–11. A peptide whose length is outside the class's accepted range is not a valid binder candidate.

### From Sette et al. (1994) / IEDB

1. **Boundary semantics:** the cutoffs are stated as strict "<" inequalities (`<50 nM`, `<500 nM`). A peptide at exactly 50 nM is therefore NOT strong by the strict tier text; a peptide at exactly 500 nM is NOT a binder by the strict tier text. The tests encode the strict-inequality reading from the source.

---

## Test Datasets

### Dataset: NetMHCpan-4.1 class I %Rank decision points

**Source:** Reynisson et al. (2020), PMC7319546 — "%Rank < 0.5% and %Rank < 2% ... for SBs and WBs for class I".

| Parameter | Value |
|-----------|-------|
| Strong binder %Rank cutoff (class I) | < 0.5 |
| Weak binder %Rank cutoff (class I) | < 2.0 |
| %Rank = 0.5 (boundary) | NOT strong (strict `<`) → Weak (0.5 < 2.0) |
| %Rank = 2.0 (boundary) | NOT weak (strict `<`) → NonBinder |
| %Rank = 0.4 | Strong |
| %Rank = 1.0 | Weak |
| %Rank = 5.0 | NonBinder |

### Dataset: IEDB / Sette class I IC50 (nM) decision points

**Source:** Sette et al. (1994), PMID 7527444; IEDB threshold tiers (high <50, intermediate <500).

| Parameter | Value |
|-----------|-------|
| Strong binder IC50 cutoff | < 50 nM |
| Weak (intermediate) binder IC50 cutoff | < 500 nM |
| IC50 = 10 nM | Strong |
| IC50 = 50 nM (boundary) | NOT strong → Weak (50 < 500) |
| IC50 = 200 nM | Weak |
| IC50 = 500 nM (boundary) | NOT weak → NonBinder |
| IC50 = 1000 nM | NonBinder |

### Dataset: peptide length validity by MHC class

**Source:** Reynisson et al. (2020) (class I 8–14, default 8–11); IEDB class II tool description (13–25).

| Parameter | Value |
|-----------|-------|
| Class I accepted length range | 8–11 (canonical neoantigen search default) |
| Class II accepted length range | 13–25 |
| Class I length 7 | invalid (too short) |
| Class I length 9 | valid |
| Class I length 12 | invalid (above canonical default range) |
| Class II length 15 | valid |
| Class II length 12 | invalid |

### Dataset: SMM transform anchor points (derived from the IEDB log50k formula)

**Source:** IEDB `log50k = 1 − log(IC50)/log(50000)` inverted to IC50 = 50000^(1−score) (Peters & Sette 2005).

| score (log50k) | Predicted IC50 (nM) | `ClassifyBindingAffinity` |
|----------------|---------------------|----------------------------|
| 0.0 | 50000 | NonBinder |
| 0.5 | 223.6067977499790 (= √50000) | Weak |
| 0.9 | 2.9505093853369178 | Strong |
| 1.0 | 1 | Strong |
| round-trip IC50 = 500 | score = 0.4256251898085073 | — (boundary) |

### Dataset: BIMAS product rule (hand-built 3-position matrix; neutral coefficient = 1.0)

**Source:** BIMAS scoring rule: T½ = finalConstant · ∏ coefficients.

| Peptide | Per-position coefficients | Final constant | Predicted T½ |
|---------|---------------------------|----------------|--------------|
| `LMV` | 2.0 · 3.0 · 1.5 = 9.0 | 10.0 | 90.0 |
| `AAA` (unlisted at every row) | 1.0 · 1.0 · 1.0 = 1.0 | 10.0 | 10.0 |

### Dataset: SMM 9-mer strong-vs-non-binder ranking (hand-computed)

**Source:** SMM additive model + IC50 = 50000^(1−score). Matrix gives `GILGFVFTL` (influenza M1 58–66, the paradigm HLA-A*02:01 binder) per-position contributions summing (with intercept) to score 1.0 → IC50 = 1 nM (Strong); a poly-`W` 9-mer gets contributions summing to score 0.0 → IC50 = 50000 nM (NonBinder). Strong binder IC50 ≪ non-binder IC50.

### Dataset: MHCflurry pan-allele affinity oracle (mhcflurry 2.1.5, models_class1_pan 20200610)

**Source:** Oracle IC50 (nM) produced in-session by re-running the MHCflurry Python API. **Single-network** column = the smallest ensemble member (feedforward `[512,512]`, `PAN-CLASS1-1-3ed9fb2d2dcc9803`), which is the network embedded for the C# parity test. **Full-ensemble** column = geometric mean over all 10 members (`Class1AffinityPredictor.predict`). SIINFEKL/HLA-A\*02:01 = a self peptide (non-binder, µM range); GILGFVFTL (influenza M1 58–66), NLVPMVATV (CMV pp65), ELAGIGILTV (MART-1 26–35 A27L), AAAWYLWEV, SLYNTVATL (HIV Gag), CINGVCWTV = known HLA-A\*02:01 binders.

| Peptide | Allele | Single-net IC50 (nM) | Full-ensemble IC50 (nM) |
|---------|--------|----------------------|--------------------------|
| SIINFEKL | HLA-A\*02:01 | 11483.195201 | 11927.160413 |
| GILGFVFTL | HLA-A\*02:01 | 19.123150 | 19.955122 |
| NLVPMVATV | HLA-A\*02:01 | 17.542640 | 16.570969 |
| ELAGIGILTV | HLA-A\*02:01 | 119.054961 | 83.553826 |
| AAAWYLWEV | HLA-A\*02:01 | 16.559303 | 18.886903 |
| SIINFEKL | HLA-B\*07:02 | 28830.796646 | 29061.956599 |
| SLYNTVATL | HLA-A\*02:01 | 28.972028 | 35.353629 |
| CINGVCWTV | HLA-A\*02:01 | 92.105940 | 99.784630 |

The C# port reproduces both the single-network IC50s (tested in CI against the embedded member) and the full-ensemble IC50s (verified in-session against the live model) to within **< 0.03%** relative error.

---

## Assumptions

1. **ASSUMPTION: Full 10-network ensemble weights not embedded (size).** The MHCflurry `models_class1_pan` weights total ≈ 80 MB of float32, which are near-incompressible (the `.npz` files are already zip-format; gzip → 75 MB). Embedding 80 MB in the source repo would be the single largest data artifact and is declined for repo health. The peptide/allele encoders, the pseudosequence table (Apache-2.0, bundled, ~0.7 MB), the forward-pass engine, the `to_ic50` transform, and the geometric-mean combiner are all ported and source-verified; **one** ensemble member (~4.6 MB) is embedded as a test fixture for CI-portable forward-pass parity; the full ensemble is loaded from a caller-supplied MHCflurry weight pack via `LoadWeightPack`. This is not a correctness assumption — the algorithm is exact; it is a packaging boundary (analogous to the caller-supplied LM22 matrix in ONCO-IMMUNE-001).

1. **ASSUMPTION: Caller-supplied coefficient matrix (matrix-based prediction).** No redistributable, cross-verifiable trained HLA coefficient matrix could be retrieved this session: the public BIMAS coefficient files are served by a now-defunct dynamic CGI (not archived), the Parker (1994) 180-value table is paywalled, and the IEDB SMM matrices carry a non-commercial / no-redistribution licence. The library therefore embeds only the published *scoring rules* (BIMAS product; SMM `IC50 = 50000^(1−score)`) and a `LoadScoringMatrix` loader, and the caller supplies the matrix values under their own licence. This mirrors ONCO-IMMUNE-001 (CIBERSORT LM22 caller-supplied). The scoring rules are fully sourced and cross-verifiable; the trained weights are not embedded.

1. **ASSUMPTION: Class I canonical length range = 8–11.** The source gives 8–14 as the full class I range with 8–11 as the default. This unit adopts 8–11 as the accepted class I range to match the existing `OncologyAnalyzer.MhcClassIMin/MaxPeptideLength` constants (ONCO-NEO-001) and the pVACtools canonical neoantigen search. This is a documented default, not an invented value; callers can pass an explicit range. It affects `IsValidPeptideLength` output for lengths 12–14.

---

## Recommendations for Test Coverage

1. **MUST Test:** IC50 classification at and around both cutoffs (10/50/200/500/1000 nM) yields Strong/Weak/NonBinder per the strict-`<` tiers — Evidence: Sette 1994; IEDB tiers.
2. **MUST Test:** %Rank classification at and around both cutoffs (0.4/0.5/1.0/2.0/5.0) yields Strong/Weak/NonBinder for class I — Evidence: Reynisson 2020.
3. **MUST Test:** class II %Rank cutoffs (2.0 / 10.0) classify SB/WB — Evidence: Reynisson 2020.
4. **MUST Test:** peptide length validity for class I (8–11) and class II (13–25) at boundaries — Evidence: Reynisson 2020; IEDB class II tool description.
5. **SHOULD Test:** invalid inputs rejected — negative/NaN/infinite IC50, negative/NaN %Rank, %Rank > 100 — Rationale: %Rank is a percentile in [0,100]; IC50 is a positive concentration (invariant IC50 > 0).
6. **COULD Test:** a combined "classify candidate" helper that gates on length validity then affinity — Rationale: mirrors how a caller wires ONCO-NEO-001 windows to a supplied affinity.
7. **MUST Test (prediction):** SMM transform reproduces the exact derived anchor IC50 values (score 0→50000, 0.5→√50000, 1→1) — Evidence: IEDB log50k formula.
8. **MUST Test (prediction):** `PredictAndClassifySmm` chains predict→classify, ranking a strong-binder score (→ ≤2.95 nM, Strong) far below a non-binder score (→ 50000 nM, NonBinder) — Evidence: IEDB formula + classification cutoffs.
9. **MUST Test (prediction):** BIMAS product rule = finalConstant · ∏ coefficients with unlisted residues = 1.0 — Evidence: BIMAS scoring documentation.
10. **SHOULD Test (prediction):** loader parses `CONST=`/`RESIDUE=VALUE`, comments, blanks, upper-cases residues, round-trips into a prediction; malformed/non-numeric tokens → `FormatException`; null peptide / empty matrix / length mismatch → documented exceptions — Rationale: loader is the caller's entry point.

---

## References

1. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449–W454. https://doi.org/10.1093/nar/gkaa379 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC7319546/)
2. Sette A, Vitiello A, Reherman B, et al. (1994). The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *Journal of Immunology* 153(12):5586–5592. https://pubmed.ncbi.nlm.nih.gov/7527444/
3. Roomp K, Antes I, Lengauer T (2010). Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics* 11:90. https://doi.org/10.1186/1471-2105-11-90 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC2836306/)
4. IEDB. What thresholds (cut-offs) should I use for MHC class I and II binding predictions. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
5. IEDB. T Cell Epitopes - MHC Class II Binding Prediction Tools Description. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description
6. Parker KC, Bednarek MA, Coligan JE (1994). Scheme for ranking potential HLA-A2 binding peptides based on independent binding of individual peptide side-chains. *J. Immunol.* 152(1):163–175. https://pubmed.ncbi.nlm.nih.gov/8254189/
7. BIMAS — Information & background on the HLA peptide motif searches (NIH/CIT/CBEL; R. Taylor & K. Parker). http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html (archived: https://web.archive.org/web/20041016022153/http://www-bimas.cit.nih.gov/molbio/hla_bind/hla_motif_search_info.html)
8. Peters B, Sette A (2005). Generating quantitative models describing the sequence specificity of biological processes with the stabilized matrix method. *BMC Bioinformatics* 6:132. https://doi.org/10.1186/1471-2105-6-132
9. IEDB MHC class I log50k linearisation `log50k = 1 − log(IC50)/log(50000)` (restated in dmnfarrell, "Create an MHC-Class I binding predictor in Python"). https://dmnfarrell.github.io/bioinformatics/mhclearning
10. O'Donnell TJ, Rubinsteyn A, Laserson U (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing. *Cell Systems* 11(1):42–48.e7. https://doi.org/10.1016/j.cels.2020.06.010
11. MHCflurry source code and trained models (Apache License 2.0), OpenVax. Source modules `amino_acid.py`, `encodable_sequences.py`, `allele_encoding.py`, `class1_neural_network.py`, `regression_target.py`, `ensemble_centrality.py`, `class1_affinity_predictor.py`; `models_class1_pan` release 20200610. https://github.com/openvax/mhcflurry ; licence https://github.com/openvax/mhcflurry/blob/master/LICENSE

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-06-25**: Added matrix-based prediction — BIMAS product rule (Parker 1994 / BIMAS docs) and SMM `IC50 = 50000^(1−score)` (Peters & Sette 2005 / IEDB log50k). Documented that no redistributable trained matrix was obtainable → matrix is caller-supplied (Framework).
- **2026-06-25**: Added the ported MHCflurry 2.0 Class I pan-allele binding-AFFINITY neural network (O'Donnell et al. 2020; Apache-2.0 source + weights) — BLOSUM62 `left_pad_centered_right_pad` peptide encoding, 37-residue allele pseudosequence (bundled), feed-forward forward pass (feedforward + with-skip-connections), `IC50 = 50000^(1−x)`, geometric-mean ensemble. Oracle = `mhcflurry` 2.1.5 / `models_class1_pan` 20200610; C# port verified to <0.03%. Full 80 MB ensemble weights not embedded (size); one member embedded for CI parity, full ensemble loaded via `LoadWeightPack`.
