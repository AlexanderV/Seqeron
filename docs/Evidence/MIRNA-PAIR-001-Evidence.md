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
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 4 (Wikipedia citing the primary: Crick FHC, J. Mol. Biol. 19(2):548–555, 1966)

**Key Extracted Points:**

1. **Canonical pairs:** "The standard base pairings are A-U and G-C."
2. **Wobble pair in RNA:** the principal RNA wobble pair is **G-U** (others listed involve inosine, not relevant to standard RNA duplexes). Source primary: Crick FHC, "Codon–anticodon pairing: the wobble hypothesis," J. Mol. Biol. 19(2):548–555 (1966).

### Lewis BP et al. (2005) — Conserved seed pairing… (PubMed 15652477, search-result abstract)

**URL:** https://pubmed.ncbi.nlm.nih.gov/15652477/
**Accessed:** 2026-06-13 (WebSearch query "Lewis 2005 Conserved seed pairing targetscan reverse complement miRNA 3'UTR antiparallel")
**Authority rank:** 1 (peer-reviewed, Cell)

**Key Extracted Points:**

1. **Reverse-complement seed match:** Targets are identified as "mRNAs with conserved complementarity to the seed (nucleotides 2-7) of the miRNA" — i.e. the target site is the reverse complement of the miRNA seed read antiparallel.

### NNDB — Turner 2004 RNA folding parameters (Mathews/Turner)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/index.html (and rna_2004_watson_crick_helices.html)
**Accessed:** 2026-06-13 (WebSearch; direct WebFetch of NNDB pages returned 404/blocked — see Assumptions)
**Authority rank:** 1/3 (primary thermodynamic database; Xia et al. 1998 / Turner 2004)

**Key Extracted Points:**

1. **Nearest-neighbor model:** RNA duplex/helix folding free energy is computed by summing nearest-neighbor stacking free energies over consecutive base pairs, plus initiation and terminal-AU/GU penalties (Xia et al. 1998; Turner 2004). The 21 distinct base-pair stacks are firmly grounded in experimental melting data (search result "Analysis of RNA nearest neighbor parameters…", rnajournal.cshlp.org/content/24/11/1568.full).
2. **Sign:** all Watson-Crick nearest-neighbor stacking free energies at 37 °C are negative (stabilising). The numeric stacking table used by the implementation is the same NNDB Turner 2004 set already cited in `MiRnaAnalyzer.StackingEnergies` (committed in a prior unit). Exact per-stack values were NOT independently re-retrieved this session (NNDB pages 404 to the fetcher) — see ASSUMPTION 1.

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

1. **ASSUMPTION: Turner 2004 stacking numeric values not re-retrieved this session** — The exact per-stack free-energy values reused by the duplex free-energy estimate (`StackingEnergies`) come from the NNDB Turner 2004 set already committed in a prior unit; the NNDB pages returned 404 to the web fetcher this session, so the individual numbers could not be independently re-opened now. Consequence: tests for `AlignMiRnaToTarget` assert only base-pairing structure (matches/mismatches/wobbles, alignment symbols, orientation) and the **sign** invariant of the free energy (ΔG ≤ 0 for a fully paired duplex, ≥ 0 sign behaviour for an all-mismatch duplex), NOT exact kcal/mol values. The free-energy magnitude is documented as "Intentionally simplified".

---

## Recommendations for Test Coverage

1. **MUST Test:** `CanPair` returns true exactly for A-U, U-A, G-C, C-G, G-U, U-G and false otherwise — Evidence: Agarwal et al. (A-U,C-G) + Crick (1966) (G-U).
2. **MUST Test:** `IsWobblePair` true only for G-U/U-G; false for Watson-Crick pairs — Evidence: PMC4870184 distinguishes wobble from Watson-Crick.
3. **MUST Test:** `GetReverseComplement` produces the antiparallel Watson-Crick reverse complement (RNA, T→U) including the let-7a seed `GAGGUAG`→`CUACCUC` — Evidence: A-U/G-C rules (Agarwal); Lewis (2005) reverse-complement seed.
4. **MUST Test:** `AlignMiRnaToTarget` perfect complement → all Watson-Crick matches, 0 mismatch/wobble — Evidence: Watson-Crick rules.
5. **MUST Test:** `AlignMiRnaToTarget` G:U duplex → counted as wobbles not matches — Evidence: Crick (1966); PMC4870184.
6. **MUST Test:** `AlignMiRnaToTarget` non-pairing duplex → all mismatches — Evidence: A pairs only with U.
7. **MUST Test:** `AlignMiRnaToTarget` empty input → empty duplex (defensive) — Evidence: implementation contract (no source defines behaviour for empty).
8. **SHOULD Test:** free energy of a fully paired duplex is ≤ 0 (stabilising); all-mismatch duplex not stabilising — Rationale: NNDB Turner stacking energies are negative for paired stacks.
9. **COULD Test:** DNA input (T) normalised to RNA (U) before pairing — Rationale: implementation normalises T→U.

---

## References

1. Bartel DP (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell 136(2):215–233. https://www.cell.com/fulltext/S0092-8674(09)00008-7 (DOI: 10.1016/j.cell.2009.01.002)
2. Broughton JP, Lovci MT, Huang JL, Yeo GW, Pasquinelli AE (2016). Pairing beyond the Seed Supports MicroRNA Targeting Specificity / Transcriptome-Wide Non-Canonical Interactions. PMC4870184. https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/
3. Agarwal V, Bell GW, Nam JW, Bartel DP (2015). Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/ (DOI: 10.7554/eLife.05005)
4. Crick FHC (1966). Codon–anticodon pairing: the wobble hypothesis. J. Mol. Biol. 19(2):548–555 — via Wikipedia "Wobble base pair". https://en.wikipedia.org/wiki/Wobble_base_pair
5. Lewis BP, Burge CB, Bartel DP (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120(1):15–20. https://pubmed.ncbi.nlm.nih.gov/15652477/
6. Turner DH, Mathews DH (2010). NNDB: the nearest neighbor parameter database. Nucleic Acids Res. 38:D280–D282; Turner 2004 parameters at https://rna.urmc.rochester.edu/NNDB/turner04/index.html (Xia et al. 1998, Biochemistry 37:14719–14735).

---

## Change History

- **2026-06-13**: Initial documentation.
