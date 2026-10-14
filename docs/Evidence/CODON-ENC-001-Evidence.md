# Evidence Artifact: CODON-ENC-001

**Test Unit ID:** CODON-ENC-001
**Algorithm:** Effective Number of Codons (ENC / Nc), Wright 1990
**Date Collected:** 2026-06-13

---

## Online Sources

### Fuglsang A. (2004) "The 'effective number of codons' revisited" (Biochem. Biophys. Res. Commun. 317:957–964)

**URL:** https://eclass.uoa.gr/modules/document/file.php/D473/%CE%92%CE%B9%CE%B2%CE%BB%CE%B9%CE%BF%CE%B3%CF%81%CE%B1%CF%86%CE%AF%CE%B1/DNA%20Composition/Fuglsang_2004.pdf
**Accessed:** 2026-06-13
**How retrieved:** WebSearch query `Wright 1990 "effective number of codons" Nc formula gene 87` → fetched the PDF with WebFetch (binary saved locally) → extracted full text with `pdftotext`.
**Authority rank:** 1 (peer-reviewed paper that reproduces Wright's original equations verbatim)

**Key Extracted Points:**

1. **Codon homozygosity (Eq. 1, verbatim):** `F̂ = ( n·Σ_{i=1..k} p_i² − 1 ) / ( n − 1 )`, "where n is the total count for the amino acid in the gene, and p_i is the codon frequency for the ith synonymous codon for the particular amino acid." (k = number of synonymous codons for that amino acid.)
2. **Per-amino-acid effective codons (Eq. 2):** `N̂c(aa) = 1 / F̂_aa`.
3. **Gene-level aggregation (Eq. 3, verbatim):** `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`, where `F̂₂` is the average homozygosity for amino acids with degeneracy two (histidine, glutamine, etc.) "and so on".
4. **Missing-amino-acid rule (Eq. 4):** if a gene lacks (e.g.) threonine, then `F̂₄` is the average of the codon homozygosities of the remaining four-fold amino acids: `F̂₄ = (F̂_pro + F̂_gly + F̂_ala + F̂_val) / 4` — i.e. average the estimable F values within the same degeneracy class.
5. **Upper-bound rule:** "There is a chance that N̂c, calculated through use of Eq. (3), will exceed 61. In that case, Wright recommends re-adjusting the result down to 61."
6. **Isoleucine (3-fold) fallback:** "Wright suggested … using `F̂₃ = (F̂₂ + F̂₄)/2` when the isoleucine estimator was not possible to calculate." (Eq. 5a.)
7. **Range:** "In extremely biased genes the effective number of codons can approach 20, while in unbiased genes it will approach 61."
8. **Calculability constraint:** "Nc can only be calculated if there are at least two codons for each amino acid with synonymous codons in the sequence" — homozygosity needs n ≥ 2 (the (n−1) denominator is undefined for n ≤ 1).
9. **Reference organism:** "In his original paper, Wright used Escherichia coli K12 as a reference organism."

### Fuglsang A. (2006) "Estimating the 'effective number of codons': the Wright way…" (Genetics 172(2):1301–1307)

**URL:** https://academic.oup.com/genetics/article/172/2/1301/5923091
**Accessed:** 2026-06-13
**How retrieved:** WebSearch query `effective number of codons Wright degeneracy classes "nine" twofold "five" fourfold "three" sixfold isoleucine Met Trp standard genetic code` → fetched with WebFetch.
**Authority rank:** 1 (peer-reviewed, Genetics)

**Key Extracted Points:**

1. **F formula confirmed:** `F = (n·Σ p_i² − 1)/(n−1)` with `p_i = n_i/n`; Wright uses sampling **without** replacement, which the paper concludes is the superior estimator.
2. **Constraints confirmed:** Nc is constrained between 20 and 61; absent / insufficient amino acids handled by within-class averaging; isoleucine has special 3-fold handling.

### Standard genetic code degeneracy classes (via WebSearch)

**How retrieved:** WebSearch (same query as above) returned the standard-genetic-code degeneracy partition.
**Authority rank:** 2 (NCBI standard genetic code, Table 1)

**Key Extracted Points:**

1. **Degeneracy partition of the 20 amino acids + stop:** "5 quartets (4 codons each), 9 doublets (2 codons each), 3 sextets (6 codons each), 1 triplet (3 codons) and 2 singlets (1 codon each)." The two singlets are Met (ATG) and Trp (TGG); the single triplet is isoleucine (ATT/ATC/ATA). This is exactly the partition assumed by Eq. (3): the constant `2` = the two single-codon amino acids (Met + Trp), and `9, 1, 5, 3` count the two-, three-, four- and six-fold amino acids. Stop codons are excluded from Nc.

### CodonW 1.4.4 (Peden 1999) — reference implementation (review 2026-09)

**Opened:** original source `codonw_1.4.4.orig.tar.gz` (compiled locally): `codon_us.c` `enc_out`, `codonW.h` genetic-code tables, `README_indices.txt`. Also downloaded https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf (13-page excerpt; text extraction unavailable in the sandbox, not used).
**Authority rank:** 1 (the de-facto reference implementation of Wright's Nc; cites Wright 1990)

**Key Extracted Points:**

1. `README_indices.txt`: "If amino acids are rare or missing, adjustments must be made. When there are no amino acids in a synonymous family, Nc is not calculated as the gene is either too short or has extremely skewed amino acid usage (Wright 1990). An exception to this is made for genetic codes where isoleucine is the only 3-fold synonymous amino acid, and is not used in the protein gene. The reported value of Nc is always between 20 … and 61 … If the calculated Nc is greater than 61 … it is adjusted to 61."
2. `enc_out`: per amino acid (stops skipped) `bb = (n·Σp² − 1)/(n − 1)` (0 when n ≤ 1); the amino acid enters its class sum/count only `if (bb > 0.0000001)` — F̂ = 0 is **not** an estimate; `fold[z]` counts every amino acid of the code with z codons; `enc = fold[1] + Σ_z fold[z]/mean_z`; missing class z: `z==3 && numaa[2] && numaa[4] && fold[3]==1` → mean₃ = (mean₂ + mean₄)/2, else print `*****` (not calculated); print `61.00` when enc > 61.
3. Classes come from the selected genetic code (`-code`); CodonW supports 8 codes (NCBI 1, 2, 3, 4, 5, 6, 10, 9).

### codonbias 0.5.0 (PyPI `codon-bias`) — secondary reference

`codonbias/scores.py` `EffectiveNumberOfCodons(robust=False, pseudocount=0, mean='unweighted', genetic_code=t)` reduces to Wright's F̂ (Novembre's F with uniform expectation), drops F ≤ 1e-6 and n ≤ 1, averages per degeneracy class of NCBI table t and caps at `len(P)` = number of sense codons. It imputes F = 1/deg for an empty class (differs from CodonW; not followed for Wright). Its defaults are the Sun, Yang & Xia (2013) variant.

### codonbias 0.5.0 — reference for the opt-in Novembre 2002 / Sun–Yang–Xia 2013 estimators (review 2026-10-14, F50)

Wheel `codon_bias-0.5.0-py3-none-any.whl` (PyPI `codon-bias`), read in full: `scores.py` `EffectiveNumberOfCodons.__init__/_calc_score/_calc_F/_calc_BNC/_calc_BCC` (lines 482–660) and `stats.py` `CodonCounter.count_array`, `BaseCounter.count_array`. With `k_mer = 1`:
1. Counts of the sense codons of NCBI table t (`ignore_stop=True`; upper-case, U→T, triplets with non-ACGT dropped) **+ pseudocount** on every sense codon; per amino acid N = Σn_i, P_i = n_i/N. Families are whole amino acids (a 6-fold family is **not** split into 2 + 4).
2. Background (`bg_correction=True`): base frequencies = (A, C, G, T counts over the **whole** sequence — every position, stop codons and trailing bases included — each + 1) / total (`_calc_BNC`); or of `get_score(…, background=seq)`. Expected codon frequency e_i = product of its three base frequencies renormalised within the amino acid (`_calc_BCC`); uniform 1/k otherwise.
3. χ²_aa = Σ_i (P_i − e_i)²/e_i (NaN entries — N = 0 — skipped ⇒ 0).
4. `robust=True` (Sun, Yang & Xia 2013): F = (χ² + 1)/k (= Σp² for uniform e). `robust=False` (Novembre 2002): F = (N·χ² + N − k)/((N − 1)·k) (= Wright's F̂ for uniform e).
5. Class mean over amino acids with N > 1, F > 1e-6, finite F: weighted by N (`mean='weighted'`, default) or unweighted.
6. Empty class ⇒ F̄ = 1/z; then if the 3-fold class is empty F̄₃ = (F̄₂ + F̄₄)/2 (using imputed values if needed).
7. Nc = Σ_z K_z/F̄_z, `min(len(P), Nc)` (sense-codon cap; no lower clamp).
Novembre 2002 preset = `bg_correction=True, robust=False, pseudocount=0, mean='unweighted'`; Sun–Yang–Xia 2013 preset = codonbias defaults `robust=True, pseudocount=1, mean='weighted', bg_correction=False`.

### ENCprime (Novembre's own program) — reproduced opt-in (`EncPrimeNc` / `EncPrimeNcPrime`, review 2026-10-14, A3-2)

**Source:** github.com/jnovembre/ENCprime (`git clone`, HEAD 709ace50b12fe30ce8d17a3c5d47c83e7f0ccd23, 2021-12-10), `src/ENCprime.c` (header "Date: July 17th, 2002"; 3-fold fix "JN 2/28/06"), `src/SeqCount.c`, `README`, `example/`. Compiled locally with gcc and run as `ENCprime <codcnt> <acgtfreq> <code> <out> 0 -q`.

**Rules read from the source (all reproduced):**

1. `get_code`: synonymous families = amino-acid letters of the built-in GenBank table (1–6, 9–16, 21–23); `redun` = family size; stops (`*`) excluded from Nc (`aa[j].abbrev != '*'`).
2. `calc_cod_data`: `ex_num` = Π f(base) over the codon's **variable sites** (positions where some synonymous codon differs), `ex = ex_num / Σ ex_num` within the amino acid (`ex_denom == 0` → 1). f_A, f_C, f_G, f_T are read from the nucleotide-composition file written by `SeqCount -n` = base counts / sequence length (no pseudocount) of **any** FASTA — the shipped example uses a separate `ExBackgroundSeqs.fasta`; the README workflow applies it to the genes themselves. Only ratios matter (scale cancels on renormalisation).
3. `test_null`: χ² = Σ n·(p − ex)²/ex per amino acid; χ² = −2 for `redun ≤ 1`; an `ex < eps² = 1e-6` sets χ² = −1 and the loop continues accumulating.
4. `calc_F`: "if enough observations that chi-sq is legit (i.e. >5)": for n > 5, `F_p = (χ² + n − redun)/(redun·(n − 1))`, `F = (n·Σp² − 1)/(n − 1)`; otherwise both −1 (not used).
5. `calc_Nc`: class mean over amino acids with F ≠ −1; missing class → `(F̄₂ + F̄₄)/2` for 3-fold when both exist (arithmetic mean, the 2006 fix of Fuglsang's harmonic-mean bug report), else `1/i`; then "Force F_bar to be bounded (1/k,1)": `if (F_bar < 1/k) F_bar = 1/k`; `Nc += n_fold_o/F_bar` (n_fold_o = all non-stop amino acids of the class). No cap.
6. Output `%.4f` (`Nc`, `Ncp` columns).

**Program defects not reproduced:** `get_code` sets `n_aa = i + 1` for every codon, i.e. the index of the last codon's (GGG) amino acid — under table 13 G first appears at AGA, so V, A, D, E (indices beyond G) are silently dropped (M3: 37.6786 instead of 45.9643). Verified by rebuilding with `if (i+1 > n_aa) n_aa = i+1`: the library then agrees on all 397 genes. `SeqCount` allocates `calloc(seqlength)` without the terminator (heap overflow, aborts on large inputs — README: "unresolved bug in SeqCount") and indexes non-ACGT codons out of bounds; the library counts codons with its canonical counter (U = T, ambiguous triplets skipped). `GammaP` aborts when a p-value χ² is negative (p-values only); the oracle runs used a copy returning NaN there.

**Oracle runs (2026-10-14):**

- Shipped example: `SeqCount -c ExSeqs.fasta` + `SeqCount -n ExBackgroundSeqs.fasta` → ENCprime reproduces `example/ExResults` (Nc identical; Ncp ±1e-4 from the 6-decimal composition file); C# (`Background` = each gene's background sequence) equals all 9 Nc and 9 Ncp at 4 dp (18/18).
- 397 genes (9 ENCprime example genes, with and without their backgrounds; the 7 ENC test genes M3/Det/Over/R53/R203/R15/Phe, with and without background `GGGCCCGCGATATGCGCGGCCGCAAT`; 306 CodonW-run genes of WP16; 60 random skewed-composition genes, 15–900 codons, seed 20261014) × 17 built-in codes (1–6, 9–16, 21–23), composition file written at full precision: **6352/6352 Nc and 6352/6352 Ncp identical at 4 dp** for the 16 codes other than 13; code 13: 397/397 identical to the `n_aa`-fixed program.
- Rule probes (n = 5 vs 6, missing Ile imputation, 1/k floor, e < 1e-6 sentinel, no sense codon, Met/Trp under table 2, RNA) × tables 1/2/11: 30/30 Nc and 30/30 Ncp identical at 4 dp.
- Faithful pipeline (`SeqCount -c` / `SeqCount -n` on the same FASTA, buffer-fixed SeqCount, table 1, 287 clean genes): Nc 287/287 identical; Ncp 271/287 identical, the rest within 8e-4 — from SeqCount's 6-decimal `%f` composition output (the full-precision file gives 287/287).

### Fuglsang (2006) — re-checked for A3-14 (2026-10-14)

**Opened:** WebSearch results for the title (PubMed 16299393 abstract snippet; Semantic Scholar record de3ea6bf9683c59430d63434c97929d7d625bf2f); PMC1456227 / academic.oup.com unreachable from the sandbox (DNS blocked). The abstract: five new Nc estimators derived and tested with the two published ones by resampling; F can be estimated in two closely related ways corresponding to sampling **with or without replacement**, the latter used by Wright; the methods based on sampling without replacement were much more accurate at short gene lengths ("the Wright way … leads to superior estimates"), while those with replacement correlated better with E. coli mRNA levels. Nothing in it describes a "rounding" variant: the previous algorithm-doc bullet ("Fuglsang (2006) N̂c sampling-without-replacement rounding variant — not implemented") was wrong — Wright's without-replacement F̂ is exactly what `Wright1990` (and `EncPrimeNc`) compute. Bullet corrected to list the with-replacement alternatives as the unimplemented ones.

### Numerical cross-check (2026-09-28)

- C# vs CodonW 1.4.4 binary (`-enc -nomenu -silent -machine -code c`, c = 0..7): 807 genes (random composition, codon-sampled with skewed weights, missing-family/short genes, RNA, lower case, N) × 8 codes = 6456 cases: 3136 "*****" ⇔ C# 0 (all), all others equal at CodonW's 2-dp precision; the only 5 differences are CodonW's fixed 61 cap under codes with 62/63 sense codons.
- C# vs Python port of `enc_out` over all 27 NCBI tables: 21789 cases, max |Δ| < 1e-9.
- Review 2026-10-14 (F36, `TryCalculateEnc`): CodonW 1.4.4 rebuilt from the Ubuntu `codonw_1.4.4.orig.tar.gz` and run as `codonw genes.fa -base -enc -gc -gc3s -L_aa -code c -silent -nomenu -nowarn -machine` (c = 0..7) on 304 genes (seed 20261014: 4 hand-made — terminal stop, Phe-only, internal stops + NNN, stops only — and 300 random 1–400-codon genes with N/R/Y, lower-case, RNA spelling, trailing partial codons) = 2432 cases: all 273 `*****` ⇔ `TryCalculateEnc` false / `IsEncCalculable` false / `CalculateEnc` 0; all 2159 calculated ⇔ true and equal at 2 dp, except the documented cap (CodonW prints values up to 61 where table 2 caps at 60, or 61 where tables with 62/63 sense codons do not). Locked: gene M3 41.288461538461526 / 45.0 (tables 1/2), random gene r53 42.85340746733685 / 40.37735849056604 (CodonW 42.85 / 40.38), r155 `*****` (tables 1/2).
- C# vs codonbias 0.5.0 (Wright mode) over all 27 tables, 4050 fully-populated genes (1237 at/above 61): 0 mismatches (> 1e-6).
- Review 2026-10-14 (F50, Novembre / Sun–Yang–Xia): 306 genes (M3, M10, M12, Phe-only, short, Met/Trp-only, 300 random 1–400-codon genes, seed 20261014, with N, lower case, RNA, trailing bases) × 27 NCBI tables × 7 option sets (Novembre preset; Sun–Yang–Xia preset; SYX + bg_correction; Novembre pseudocount 0.5 weighted; SYX pseudocount 0 unweighted; Novembre without bg; Novembre with `background="GGGCCCGCGATATGCGCGGCCGCAAT"`) = 57834 cases: 57351 equal to codonbias with max relative difference 2.2e-15; the 483 others are genes with no sense codon under that table (e.g. `TAG`, `AGG` under table 2), where Seqeron returns `false` (contract, as F36) and codonbias the imputed maximum. Locked to 1e-9 in `CodonUsageAnalyzer_EncVariants_Tests.cs` (6 genes × tables 1/2 × the option sets) plus the hand derivation SYX(Phe TTT×3 TTC×1) = 1028/17 = 60.470588235294116.
- Review 2026-10-14 (F51, CodonW cap): the same 306 genes through CodonW 1.4.4 (`-enc -code c -silent -nomenu -nowarn -machine`, c = 0..7 = NCBI 1, 2, 3, 4, 5, 6, 10, 9) = 2448 cases: `EncOptions { CodonWCap = true }` equals CodonW at 2 dp in 2447 (incl. all 1289 `*****` and all 167 cases where the CodonW cap differs from the sense-codon cap: table 2 values in (60, 61], tables 3/4/5/9/10 values in (61, 62], table 6 in (61, 63]); the single difference is r191 under NCBI 3, 54.814998869… vs CodonW 54.82 — a 2-dp rounding boundary of CodonW's single-precision `float enc_tot`, independent of the cap. Locked: M12 gene 61/60/62/62/62/63/62/62 → 61 under NCBI 1/2/3/4/5/6/9/10; r203 table 2 60 → 60.71285921278657 (CodonW 60.71); r15 table 3 61.89427123695478 → 61 (CodonW 61.00).

---

## Documented Corner Cases and Failure Modes

### From Fuglsang (2004)

1. **Amino acid with n ≤ 1:** F̂ is undefined (denominator n−1 → 0); such an amino acid cannot contribute its own F. Nc requires ≥ 2 codons for each represented amino acid; the within-class average (Eq. 4) is the prescribed substitute when a class member is missing.
2. **Empty degeneracy class:** if no amino acid in a class has an estimable F (e.g. no isoleucine at all), the class average is undefined. For the 3-fold class Wright gives the explicit fallback `F̂₃ = (F̂₂ + F̂₄)/2` (Eq. 5a).
3. **Overshoot past 61:** Eq. (3) can yield N̂c > 61 for nearly-uniform short genes; re-adjust the final value down to 61.
4. **Homozygosity becomes zero:** if all codons of an amino acid are used equally and counts are small, F̂ can be very low, making the apparent per-aa Nc exceed the degeneracy. The gene-level value is still re-adjusted to ≤ 61.

---

## Test Datasets

### Dataset: Fully unbiased gene (asymptotic limit)

**Source:** Fuglsang (2004), §"The effective number of codons" — unbiased gene gives F̂₂ = 0.5, F̂₄ = 0.25 (and by the same reasoning F̂₃ = 1/3, F̂₆ = 1/6).

| Class | F̂ (unbiased) | contribution to Nc |
|-------|---------------|--------------------|
| 2-fold | 0.5 | 9 / 0.5 = 18 |
| 3-fold (Ile) | 1/3 | 1 / (1/3) = 3 |
| 4-fold | 0.25 | 5 / 0.25 = 20 |
| 6-fold | 1/6 | 3 / (1/6) = 18 |
| singlets (Met,Trp) | — | 2 |
| **Total** | | **61** |

### Dataset: Fuglsang (2004) "no bias discrepancy" simulation (Nc = 40.5)

**Source:** Fuglsang (2004), Simulation results / Table 4: "all twofold degenerate aa have Nc = 1.5, isoleucine Nc = 2, the fourfold degenerate aa have Nc = 2.5, and the sixfold degenerate have Nc = 3.5 (these numbers sum up to 40.5 effective codons)."

| Class | per-aa Nc | count | class Nc total |
|-------|-----------|-------|----------------|
| 2-fold | 1.5 | 9 | 13.5 |
| 3-fold | 2.0 | 1 | 2.0 |
| 4-fold | 2.5 | 5 | 12.5 |
| 6-fold | 3.5 | 3 | 10.5 |
| singlets | 1.0 | 2 | 2.0 |
| **Total** | | | **40.5** |

### Dataset: Hand-derived two-fold example (exact F by Eq. 1)

Phe codons only, TTT × 3 and TTC × 1 (n = 4, p = (3/4, 1/4)):
Σp² = 9/16 + 1/16 = 10/16 = 0.625. F̂ = (4·0.625 − 1)/(4 − 1) = (2.5 − 1)/3 = 0.5. Nc(Phe) = 1/0.5 = 2.

Perfectly even Phe, TTT × 2 and TTC × 2 (n = 4): Σp² = 0.5, F̂ = (4·0.5 − 1)/3 = 1/3, Nc(Phe) = 3 (> 2, illustrating the per-aa overshoot corner case).

---

## Assumptions

0. **SUPERSEDED (2026-09):** an empty synonymous class no longer contributes its full codon count; per CodonW / Wright 1990 Nc is not calculated (`TryCalculateEnc` returns `false`; the legacy `CalculateEnc` returns 0 — F36 replaced the 0-only signal with the explicit API).
1. **ASSUMPTION: Lower clamp at 20.** Wright/Fuglsang state Nc *approaches* 20 in extreme bias and explicitly prescribe re-adjusting **down to 61** at the top. They do not prescribe a hard clamp at 20; 20 is the structural minimum (every degeneracy class collapses to one codon ⇒ Nc(aa)=1). Retaining `Math.Max(20, …)` is consistent with the stated range and cannot raise a legitimately-computed value, but it is not an explicit Wright instruction. Treated as a defensive bound, not an algorithmic parameter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Fully unbiased gene → Nc = 61 (use codon counts equal within every amino acid, large enough n that the asymptotic F values are reached, then the value re-adjusts to exactly 61). — Evidence: Fuglsang (2004) unbiased dataset.
2. **MUST Test:** Single-amino-acid two-fold gene with exact hand-derived F (TTT×3,TTC×1 ⇒ F=0.5 ⇒ Nc(Phe)=2; missing classes filled by Eq.4 within-class average — verify against an explicit derivation). — Evidence: Eq. 1, Eq. 2, Eq. 4.
3. **MUST Test:** Maximally biased gene (one codon per amino acid) → Nc = 20. — Evidence: range statement.
4. **MUST Test:** Invariant 20 ≤ Nc ≤ 61 for arbitrary inputs. — Evidence: range statement + Eq. 3 re-adjustment.
5. **MUST Test:** Isoleucine-absent gene uses F̂₃ = (F̂₂+F̂₄)/2 fallback. — Evidence: Eq. 5a.
6. **MUST Test:** null DnaSequence → ArgumentNullException; empty/whitespace string → 0. — Evidence: contract (no Wright rule for empty gene; degenerate input).
7. **SHOULD Test:** Lowercase input normalized; invalid (non-ACGT) codons skipped. — Rationale: parser robustness consistent with sibling methods.
8. **COULD Test:** Overshoot case (near-uniform short gene) re-adjusts to 61. — Rationale: Eq. 3 overshoot corner case.

---

## References

1. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene* 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9 (original; equations reproduced verbatim in refs 2 & 3 below, which were the retrieved sources).
2. Fuglsang, A. (2004). The 'effective number of codons' revisited. *Biochem. Biophys. Res. Commun.* 317(3):957–964. https://doi.org/10.1016/j.bbrc.2004.03.138 (retrieved PDF: https://eclass.uoa.gr/modules/document/file.php/D473/%CE%92%CE%B9%CE%B2%CE%BB%CE%B9%CE%BF%CE%B3%CF%81%CE%B1%CF%86%CE%AF%CE%B1/DNA%20Composition/Fuglsang_2004.pdf).
3. Fuglsang, A. (2006). Estimating the 'effective number of codons': the Wright way of determining codon homozygosity leads to superior estimates. *Genetics* 172(2):1301–1307. https://academic.oup.com/genetics/article/172/2/1301/5923091

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-09-28**: Review 2026-09 (B02 F15–F17): CodonW 1.4.4 source/binary and codonbias 0.5.0 added; not-calculated rule, F̂ = 0 exclusion, genetic-code-aware classes.
- **2026-10-14**: Review 2026-09 (B02 F50/F51, audit A2-7/A2-8): opt-in `EncMethod.Novembre2002` / `SunYangXia2013` (codonbias 0.5.0, 57351 cases ≤ 2.2e-15 rel.) and `EncOptions.CodonWCap` (CodonW 1.4.4, 2448 cases).
- **2026-10-14**: Review 2026-09 (B02 F36, audit A1-11): explicit not-calculable API `TryCalculateEnc`, `CodonUsageStatistics.IsEncCalculable`, MCP `calculable` / `isEncCalculable`; re-verified on the CodonW 1.4.4 binary (2432 cases).
- **2026-10-14**: Review 2026-09 (B02 FIN WP20, audit A3-2/A3-14): opt-in `EncMethod.EncPrimeNc` / `EncPrimeNcPrime` reproducing the ENCprime program (6352/6352 + 30/30 cases identical at its 4-dp output; table 13 `n_aa` bug not reproduced); Fuglsang (2006) re-checked — the "rounding variant" doc bullet corrected.
