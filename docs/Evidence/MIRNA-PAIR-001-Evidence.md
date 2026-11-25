# Evidence Artifact: MIRNA-PAIR-001

**Test Unit ID:** MIRNA-PAIR-001
**Algorithm:** MiRNA-Target Pairing Analysis (miRNA–mRNA duplex base pairing)
**Date Collected:** 2026-06-13

---

## Online Sources

### Bartel DP (2009) — MicroRNAs: Target Recognition and Regulatory Functions (Cell)

**URL:** https://www.cell.com/fulltext/S0092-8674(09)00008-7 (full text 403 to fetcher; facts taken from PMC mirror PMC4870184 and PMC4532895 below and search-result abstract)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed review, Cell)

**Key Extracted Points:**

1. **Seed region:** The miRNA seed comprises nucleotides 2–8 at the 5′ end of the miRNA; presentation of nucleotides 2–8 "prearranged in a geometry resembling an A-form helix would enhance both the affinity and specificity for matched mRNA segments, enabling 7–8 nt sites to suffice for most targeting functions." (retrieved via WebSearch query "Bartel 2009 MicroRNAs target recognition canonical seed pairing Watson-Crick Cell").
2. **Watson–Crick pairing:** miRNAs direct the silencing complex "primarily through Watson–Crick pairing between the miRNA seed (miRNA nucleotides 2–7) and complementary sites within the 3′ untranslated regions (3′ UTRs) of target RNAs."

### MicroRNA Target Recognition: Insights from Transcriptome-Wide Non-Canonical Interactions (PMC4870184)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 1 (peer-reviewed, PMC)

**Key Extracted Points:**

1. **Seed position:** "Functional miRNA-target interactions are known to majorly require as few as 6-nt matches within the seed region (position 2-8)."
2. **Pairing type:** Diagrams state "solid lines indicate Watson-Crick base pairing", distinguished from wobble pairs (G:U).
3. **Site types:** 6-mers at positions 1–6, 2–7, 3–8; 7-mers at positions 2–8 and 1–7; 8-mer at positions 1–8; offset 6-mer is a "6-mer match to position 3–8" with marginal repression.
4. **Consecutive pairing:** "local short stretches (≥ 6 nt) of consecutive base-pairing significantly contribute to target recognition."

### Predicting effective microRNA target sites in mammalian mRNAs — Agarwal et al. (PMC4532895, eLife 2015)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 1 (peer-reviewed, eLife/PMC)

**Key Extracted Points:**

1. **Canonical site types (Watson-Crick pairing pattern):**
   - 8mer = "Watson–Crick match to miRNA positions 2–8 with an A opposite position 1".
   - 7mer-m8 = "position 2–8 match".
   - 7mer-A1 = "position 2–7 match with an A opposite position 1".
   - 6mer = "position 2–7 match".
   - offset-6mer = "position 3–8 match".
2. **Pairing chemistry:** "These bases can each pair up with one specific other base—'A' pairs with 'U', and 'C' pairs with 'G'."
3. **A1 mechanism:** "the preference for an adenosine opposite position 1 is independent of the miRNA nucleotide identity" — recognised by Argonaute, not by miRNA base pairing.

### Crick FHC (1966) — Codon–anticodon pairing: the wobble hypothesis (via Wikipedia "Wobble base pair", which cites the primary)

**URL:** https://en.wikipedia.org/wiki/Wobble_base_pair
**Accessed:** 2026-06-13 (WebFetch) · re-fetched 2026-10-14 (SRC citation audit: see corrections below)
**Authority rank:** 4 (Wikipedia citing the primary: Crick FHC, J. Mol. Biol. 19(2):548–555, 1966)

**Key Extracted Points:**

1. **Canonical pairs:** the standard (Watson–Crick) base pairings are A-U and G-C (paraphrase — this sentence is not on the current page, which shows the Watson–Crick pairs in its pairing table; re-checked 2026-10-14).
2. **Wobble pair in RNA:** the principal RNA wobble pair is **G-U** (others listed involve inosine, not relevant to standard RNA duplexes). Source primary: Crick FHC, "Codon–anticodon pairing: the wobble hypothesis," J. Mol. Biol. 19(2):548–555 (1966).

### Lewis BP et al. (2005) — Conserved seed pairing… (PubMed 15652477, search-result abstract)

**URL:** https://pubmed.ncbi.nlm.nih.gov/15652477/
**Accessed:** 2026-06-13 (WebSearch query "Lewis 2005 Conserved seed pairing targetscan reverse complement miRNA 3'UTR antiparallel")
**Authority rank:** 1 (peer-reviewed, Cell)

**Key Extracted Points:**

1. **Reverse-complement seed match:** Targets are identified as "mRNAs with conserved complementarity to the seed (nucleotides 2-7) of the miRNA" — i.e. the target site is the reverse complement of the miRNA seed read antiparallel.

### NNDB — Turner 2004 RNA folding parameters (Mathews/Turner)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/index.html (and rna_2004_watson_crick_helices.html)
**Accessed:** 2026-06-13 (WebSearch); pages and parameter table opened directly 2026-11-25 (rna_2004/* URLs below)
**Authority rank:** 1/3 (primary thermodynamic database; Xia et al. 1998 / Turner 2004)

**Key Extracted Points:**

1. **Nearest-neighbor model (opened 2026-11-25, rna.urmc.rochester.edu/NNDB/rna_2004/rna_2004_watson_crick_helices.html §1.1):** "ΔG°37 Watson-Crick-Franklin = ΔG°37 intermolecular initiation + ΔG°37 AU end penalty (per AU end) + ΔG°37 symmetry (self-complementary duplexes) + Σ[ΔG°37 stacking]". Worked examples: self-complementary 5′AGCGCU3′ = 4.09 + 2×0.45 + 0.43 − 2.08 − 3.42 − 2.36 − 3.42 − 2.08 = −7.94; non-self-complementary 5′GCACG3′/3′CGUGC5′ = 4.09 − 3.42 − 2.11 − 2.24 − 2.36 = −6.04.
2. **Dangling ends (…/rna_2004_dangling_ends.html, table turner_2004_dangle_dg.txt):** "if a helix end is extended on both the 5′ and 3′ strands, then a terminal mismatch exists (not the sum of 5′ and 3′ dangling ends)"; example 3′ dangling C on GC −0.4, 5′ dangling A on GC −0.2. Table layout: 3′ dangle X sits 3′ of the top (5′→3′) base; 5′ dangle X sits 5′ of the bottom base — the key is (top base, X, bottom base).
3. **Terminal mismatches (…/rna_2004_terminal_mismatches.html):** example "UA followed by GA" = −1.1 added to the −7.94 helix → −9.0.
4. **Exterior loops (…/rna_2004_exterior_loops.html):** stabilised by terminal mismatches, dangling ends and coaxial stacks; a nucleotide participates in only one.

### ViennaRNA 2.8.0 (reference implementation, Python `import RNA`)

- `RNA.fold_compound(m + "&" + t).eval_structure(s)` (dangles = 2) evaluates the duplex structure of `AlignMiRnaToTarget`; `RNA.duplexfold(m, t)` is the MFE duplex (`RNAduplex`) reproduced by `FoldMiRnaDuplex`. Parameter file rna_turner2004: DuplexInit 4.10, terminal AU 0.50, 0.1-kcal-rounded stacks; no symmetry term in eval; no bulge degeneracy / special-C bonus.
- Cross-check (scratch fixtures, 11 hand-picked + 40 generated mutated complements): `AlignMiRnaToTarget` vs `eval_structure` max |Δ| 0.46 kcal/mol over all 51 (0.43 of it is the NNDB symmetry term on `GCGCGCGC`·`GCGCGCGC`, which `eval_structure` omits); 0.45 on the fixtures with 2×2 loops. Before B12 RNA-STRUCT-001 F3 (NNDB int11/int21/int22 tables + `vrna_E_internal` rules in the canonical kernel, commit 0409ca9) the 2×2-loop fixtures deviated by up to 2.2. `FoldMiRnaDuplex` vs `duplexfold`: identical pair set 49/50 (was 48/50; the other is a near-degenerate alternative, |Δ| 0.2), |Δ| ≤ 0.46 except one 1-nt C-bulge (1.37: NNDB special-C bonus + degeneracy, absent in Vienna).
- 1-nt C bulge re-check after B12 RNA-ENERGY-001 F14 (commit 8e059fe; B13 F31, 2026-11-25): F14 made `CalculateBulgeLoopEnergy(..., bulgeOn5PrimeStrand: null)` apply the special-C bonus only when a paired C flanks the bulge on both strands. `FoldMiRnaDuplex`/`EvaluateDuplexEnergy` did not pass the strand, so the C-bulge fixture `CGCAGUAUGCCAAGACU`·`GUCGUCCUUGGCCUACUGCGCA` (target C7 bulged between m13 A·U t8 and m14 G·C t6) moved −23.97 → −23.07. NNDB `rna_2004_bulge_loops.html` §6.1 (opened): ΔG(n=1) = initiation + special C bulge ("a bulged C adjacent to at least one paired C") + stack as though no bulge − RT ln(states); worked example 3 C states −0.616×1.099. Hand value 3.8 − 2.08 − 0.9 − 0.6163·ln 2 = 0.39; ViennaRNA 2.8 `eval_int_loop` 1.70 (`eval_internal.c` v2.7.2: `if (nl == 1) energy += P->stack[type][type_2]`, no special C / states). With the strand passed (n1 > 0 ⇒ miRNA strand) the fold is −23.97 in both strand orders; Vienna `eval_structure` of the same structure −22.60 = −23.97 + 1.31 (bulge terms) + 0.06 (rounding). The 1.37 `duplexfold` outlier is therefore the sourced NNDB term, not an error.
- Re-lock after B12 (2026-11-25): every internal loop of the affected fixtures equals ViennaRNA 2.8 `eval_int_loop` and the NNDB hand value (`turner_2004_internal_loop_dg.txt`: initiation(8) 2.3, (10) 2.5, (14) 2.9, AU/GU closure 0.7, first mismatch GA −1.0; `turner_2004_int21_dg.txt` 5′GGU/3′CAAA = 1.9; int22 5′GGUU/3′CACA = 0.5). `Vienna_4x4Loop` 2.55 → 2.35 (Vienna 2.4), `Vienna_5x5And7x7Loops` 3.00 → 2.80 (2.9), let-7a 1-nt-bulge site ungapped −10.81 → −11.11 (−11.0), `Fold_1x2InternalLoop` −28.00 → −28.50 (duplexfold −28.2). Old values used the pre-B12 generic formula for 1×2/2×2 loops and read the inner-pair first mismatch in the wrong orientation.
- Dangle orientation confirmed numerically: hairpin CGCGAAAGCGA (3′ A on the C·G closing pair) exterior term Vienna −1.1 = NNDB key "GAC" (3′ dangle sits on G).

### Biopython 1.88 — `Bio.Seq.reverse_complement_rna`

- "ACGR" → "YCGU"; "NNRYKMBDHVSW" → "WSBDHVKMRYNN"; "acgu-X" → "X-acgu" (case kept; Seqeron upper-cases recognised codes). Oracle for `GetReverseComplement` (canonical `SequenceExtensions.GetRnaComplementBase`).

---

## Documented Corner Cases and Failure Modes

### From Agarwal et al. (PMC4532895) / PMC4870184

1. **G:U wobble is not Watson-Crick:** wobble (G:U) pairs are explicitly distinguished from Watson-Crick pairs; a duplex aligner must count them separately.
2. **A opposite position 1 is not a base-pairing event:** the A1 preference is an Argonaute-mediated recognition, not miRNA–target complementarity, so it is not a "match" in a pure base-pairing aligner.

### From Crick (1966) / Wikipedia Wobble base pair

1. **Only G-U is a standard RNA wobble:** A-C, C-U, etc. do not pair; an aligner must treat them as mismatches.

---

## Test Datasets

### Dataset: Perfect Watson-Crick complement (derived, trivially verifiable)

**Source:** Watson-Crick rules A-U, G-C (Agarwal et al. PMC4532895)

| Parameter | Value |
|-----------|-------|
| miRNA | `AAAA` |
| target | `UUUU` |
| Alignment (antiparallel) | each miRNA A pairs target U → `||||` |
| Matches | 4 |
| Mismatches | 0 |
| GU wobbles | 0 |

### Dataset: G:U wobble duplex (derived)

**Source:** Crick (1966) wobble G-U; G:U ≠ Watson-Crick

| Parameter | Value |
|-----------|-------|
| miRNA | `GGGG` |
| target | `UUUU` |
| Alignment | each miRNA G : target U wobble → `::::` |
| Matches (Watson-Crick) | 0 |
| GU wobbles | 4 |
| Mismatches | 0 |

### Dataset: Non-pairing duplex (derived)

**Source:** A pairs only with U (Agarwal et al.)

| Parameter | Value |
|-----------|-------|
| miRNA | `AAAA` |
| target | `AAAA` |
| Alignment | A vs A cannot pair → no match symbols |
| Matches | 0 |
| GU wobbles | 0 |
| Mismatches | 4 |

### Dataset: hsa-let-7a-5p seed reverse complement (real miRNA)

**Source:** miRBase let-7a-5p `UGAGGUAGUAGGUUGUAUAGUU`; seed pos 2–8 `GAGGUAG`; Watson-Crick reverse complement (Lewis 2005)

| Parameter | Value |
|-----------|-------|
| miRNA seed (pos 2–8) | `GAGGUAG` |
| Reverse complement (target 5′→3′) | `CUACCUC` |

`GetReverseComplement("GAGGUAG")`: complement each base (G→C,A→U,G→C,G→C,U→A,A→U,G→C) = `CUCCAUC`, reversed = `CUACCUC`. Trivially verifiable from A-U/G-C rules.

---

## Assumptions

None. (The former ASSUMPTION 1 — Turner values not re-opened, sign-only tests — is retired: the NNDB pages were opened on 2026-11-25 and the energies are locked to NNDB worked examples and ViennaRNA 2.8.0.)

---

## Recommendations for Test Coverage

1. **MUST Test:** `CanPair` returns true exactly for A-U, U-A, G-C, C-G, G-U, U-G and false otherwise — Evidence: Agarwal et al. (A-U,C-G) + Crick (1966) (G-U).
2. **MUST Test:** `IsWobblePair` true only for G-U/U-G; false for Watson-Crick pairs — Evidence: PMC4870184 distinguishes wobble from Watson-Crick.
3. **MUST Test:** `GetReverseComplement` produces the antiparallel Watson-Crick reverse complement (RNA, T→U) including the let-7a seed `GAGGUAG`→`CUACCUC` — Evidence: A-U/G-C rules (Agarwal); Lewis (2005) reverse-complement seed.
4. **MUST Test:** `AlignMiRnaToTarget` perfect complement → all Watson-Crick matches, 0 mismatch/wobble — Evidence: Watson-Crick rules.
5. **MUST Test:** `AlignMiRnaToTarget` G:U duplex → counted as wobbles not matches — Evidence: Crick (1966); PMC4870184.
6. **MUST Test:** `AlignMiRnaToTarget` non-pairing duplex → all mismatches — Evidence: A pairs only with U.
7. **MUST Test:** `AlignMiRnaToTarget` empty input → empty duplex (defensive) — Evidence: implementation contract (no source defines behaviour for empty).
8. **MUST Test:** duplex free energy equals the NNDB Turner 2004 sum (worked examples −7.94, −6.04; dangle −6.44; terminal mismatch −9.47; AAAA·UUUU +2.20) and ViennaRNA eval within rounding; no pair ⇒ 0.
10. **MUST Test:** `FoldMiRnaDuplex` reproduces ViennaRNA `duplexfold` base pairs and energy (bulged let-7a sites).
9. **COULD Test:** DNA input (T) normalised to RNA (U) before pairing — Rationale: implementation normalises T→U.

---

## References

1. Bartel DP (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell 136(2):215–233. https://www.cell.com/fulltext/S0092-8674(09)00008-7 (DOI: 10.1016/j.cell.2009.01.002)
2. Broughton JP, Lovci MT, Huang JL, Yeo GW, Pasquinelli AE (2016). Pairing beyond the Seed Supports MicroRNA Targeting Specificity / Transcriptome-Wide Non-Canonical Interactions. PMC4870184. https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/
3. Agarwal V, Bell GW, Nam JW, Bartel DP (2015). Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/ (DOI: 10.7554/eLife.05005)
4. Crick FHC (1966). Codon–anticodon pairing: the wobble hypothesis. J. Mol. Biol. 19(2):548–555 — via Wikipedia "Wobble base pair". https://en.wikipedia.org/wiki/Wobble_base_pair
5. Lewis BP, Burge CB, Bartel DP (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120(1):15–20. https://pubmed.ncbi.nlm.nih.gov/15652477/
6. Turner DH, Mathews DH (2010). NNDB: the nearest neighbor parameter database. Nucleic Acids Res. 38:D280–D282; Turner 2004 parameters at https://rna.urmc.rochester.edu/NNDB/turner04/index.html (Xia et al. 1998, Biochemistry 37:14719–14735).
7. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL (2011). ViennaRNA Package 2.0. Algorithms Mol. Biol. 6:26 — ViennaRNA 2.8.0 Python bindings.
8. Cock PJA et al. (2009). Biopython. Bioinformatics 25:1422–1423 — Biopython 1.88 `Bio.Seq.reverse_complement_rna`.

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-11-25** (B13 review F11–F14): full NNDB Turner 2004 duplex energy (NNDB pages opened), ViennaRNA 2.8.0 cross-check, `FoldMiRnaDuplex` (RNAduplex model), IUPAC reverse complement (Biopython oracle); ASSUMPTION 1 retired.
