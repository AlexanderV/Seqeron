# Evidence Artifact: DISORDER-LC-001

**Test Unit ID:** DISORDER-LC-001
**Algorithm:** Low-Complexity Region Detection in Protein Sequences (SEG algorithm; Wootton & Federhen)
**Date Collected:** 2026-06-14

---

## Online Sources

### NCBI BLAST `blast_seg.c` — reference implementation (NCBI C++ Toolkit)

**URL:** https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_seg.c (also doxygen reference: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html)
**Accessed:** 2026-06-14 (fetched via WebFetch of the raw source file and of the doxygen page)
**Authority rank:** 3 (reference implementation in an established library — BLAST/NCBI)

**Key Extracted Points:**

1. **Default parameters (verbatim from the macro definitions in the fetched source):** `kSegWindow = 12`, `kSegLocut = 2.2`, `kSegHicut = 2.5`. These are the SEG trigger window length W, the trigger (locut) complexity K1, and the extension (hicut) complexity K2.
2. **Complexity measure = Shannon entropy in bits:** the doxygen page documents `static double s_Entropy(Int4 *sv)` — "Calculates entropy of an integer array" — which iterates composition counts, normalizes each count by the window length, takes the logarithm, and converts to base-2 (bits) using the constant `NCBIMATH_LN2`. This is the Shannon entropy H = −Σ pᵢ·log₂(pᵢ) of the residue composition of the window, expressed in bits per residue.
3. **Two-stage scan:** `s_EntropyOn` "Calculates entropy of a sequence window" via `s_StateOn`; segments with entropy ≤ locut trigger, and are extended while entropy ≤ hicut (mirrors the GCG/manpage description below).

---

### SEG program help (GCG / Weizmann mirror) and `ncbi-seg` Ubuntu manpage

**URL:** https://bip.weizmann.ac.il/education/materials/gcg/seg.html ; https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html ; https://rothlab.ucdavis.edu/genhelp/seg.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (official program documentation for the SEG reference tool)

**Key Extracted Points:**

1. **Default parameters (verbatim):** Window `-WINdow = 12`, trigger complexity `-LOWcut (K1) = 2.2`, extension complexity `-HIGhcut (K2) = 2.5`. The window default "sets the minimum size of first stage segment."
2. **Units and bound:** complexity is measured in **bits/residue**; "If 20 different characters were distributed randomly, but with equal probability along a sequence, then each character would add 4.322 bits of information to the sequence (log(base 2) 20)." So the maximum complexity for an amino-acid alphabet is log₂(20) ≈ 4.322 bits/residue.
3. **Stage 1 (trigger):** "Seg identifies segments having a complexity equal to or less than the cutoff in bits/residue set by this parameter" (K1 = 2.2).
4. **Stage 2 (extension):** "Seg extends the low-complexity segments found in the first stage into overlapping low-complexity segments that have a complexity equal to or less than the cutoff" (K2 = 2.5).
5. **Complexity definition reference:** complexity "is defined by equation (3) of Wootton & Federhen (1993)" with a maximum of 4.322 = log₂20.

---

### Wootton J.C. & Federhen S. (1993/1996) — primary literature (search-level confirmation)

**URL:** WebSearch surfaced Semantic Scholar entry https://www.semanticscholar.org/paper/Statistics-of-Local-Complexity-in-Amino-Acid-and-Wootton-Federhen/8c865ac68cd1b5f1ad2d69d8840ffce0e0f732ea and Oxford Bioinformatics 21(2):160 review which restates the parameters.
**Accessed:** 2026-06-14 (WebSearch result blocks read)
**Authority rank:** 1 (peer-reviewed primary literature)

**Key Extracted Points:**

1. **Method identity:** Wootton & Federhen (1993, Computers & Chemistry 17(2):149–163; 1996, Methods Enzymol 266:554–571) defined "local compositional complexity" and the SEG algorithm for partitioning proteins into low- and high-complexity segments.
2. **Three user parameters:** "trigger window length (W), trigger complexity (K1), and extension complexity (K2). Each trigger window is then extended into a contig in both directions by merging with extension windows of length W and complexity less than or equal to K2." Default parameters "W = 12, K1 = 2.2 bits, K2 = 2.5 bits."
3. **Complexity measure:** the K2 / complexity measure "can be described by Shannon's Entropy" (the bits/residue entropy of the window composition), consistent with the reference implementation above.

---

## Documented Corner Cases and Failure Modes

### From NCBI/GCG SEG documentation

1. **Window longer than sequence:** the first-stage scan requires at least one full trigger window of length W; a sequence shorter than W has no window to evaluate and yields no segments.
2. **Maximal-complexity windows:** a window containing W distinct residues has entropy log₂(W); for W = 12 that is ≈ 3.585 bits > K2 (2.5), so it neither triggers nor extends.
3. **Homopolymeric / near-homopolymeric windows:** entropy 0 (single residue) up to ~1 bit (two residues) is far below K1, the canonical low-complexity case the algorithm is designed to detect.

---

## Test Datasets

### Dataset: Hand-derived window entropies (Shannon entropy, bits)

**Source:** Direct evaluation of H = −Σ pᵢ·log₂(pᵢ) for window length L = 12 (matches `s_Entropy` in NCBI `blast_seg.c`). Independently computed.

| Window composition (L = 12) | Distinct residues | H (bits) | Relation to K1=2.2 / K2=2.5 |
|------------------------------|-------------------|----------|------------------------------|
| 12 × one residue (e.g. `QQQQQQQQQQQQ`) | 1 | 0.000000 | ≤ K1 → triggers |
| 6 + 6 two residues (e.g. `AAAAAALLLLLL`) | 2 | 1.000000 | ≤ K1 → triggers |
| 3 + 3 + 3 + 3 four residues (`AAABBBCCCDDD`) | 4 | 2.000000 | ≤ K1 → triggers; > 0.5 (strict K1) → no trigger |
| 12 distinct (`ACDEFGHIKLMN`) | 12 | 3.584963 | > K2 → no trigger, no extension |

### Dataset: Maximum amino-acid complexity

**Source:** GCG/NCBI SEG help (point 2 above).

| Parameter | Value |
|-----------|-------|
| Max complexity (20-letter alphabet) | log₂(20) = 4.321928 bits/residue |

---

## Assumptions

1. **ASSUMPTION: Region-type label string (`"X-rich"`, `"X/Y-rich"`)** — The SEG specification defines *where* low-complexity segments are, not a textual composition label. The repository adds a convenience classification: the single most frequent residue when its fraction > 0.5 (`"X-rich"`), else the top two residues (`"X/Y-rich"`). This labelling is a presentation extension, not part of Wootton & Federhen; only the dominant-residue threshold (> 50 %) affects which label is produced and it does not change segment boundaries. Documented as a deviation, not source-derived.
2. **ASSUMPTION: Greedy single-residue extension** — The reference SEG extends by merging length-W extension windows with complexity ≤ K2. The repository implements the equivalent contig-growth by extending one residue at a time while the *whole growing segment's* entropy stays ≤ K2. For the homopolymer/dipeptide test cases used here the boundaries are determined by the trigger spans and are identical; documented as an implementation variant.

---

## Recommendations for Test Coverage

1. **MUST Test:** Homopolymer of length ≥ W (e.g. 26×Q) → exactly one segment spanning [0, n−1]; window entropy 0 ≤ K1. — Evidence: NCBI/GCG defaults + s_Entropy (H=0 for single residue).
2. **MUST Test:** Sequence whose every W-window has W distinct residues (entropy log₂12 ≈ 3.585 > K2) → no segments. — Evidence: GCG max-complexity statement; hand-derived H.
3. **MUST Test:** Four-types-×-3 window `AAABBBCCCDDD` (H = 2.0): triggers at default K1=2.2 (one segment) but NOT at a strict K1=0.5. — Evidence: hand-derived H = 2.0 vs the two cutoffs.
4. **MUST Test:** Two-residue block `12×A + 12×L` (each window H ≤ 1.0 ≤ K1) → single merged segment over [0,23]. — Evidence: hand-derived H = 1.0.
5. **MUST Test:** Two homopolymer runs separated by a high-complexity spacer → two distinct segments, neither overlapping the spacer. — Evidence: trigger/extension semantics (stage 1+2) + entropy of spacer windows > K2.
6. **SHOULD Test:** `minLength` filter removes a segment shorter than the threshold. — Rationale: documented minimum-length post-filter behavior.
7. **SHOULD Test:** Sequence shorter than W → empty. — Rationale: no full trigger window exists.
8. **COULD Test:** Case-insensitivity (lowercase input gives identical segments). — Rationale: implementation upper-cases input; sequences are case-insensitive by convention.

---

## References

1. Wootton J.C., Federhen S. (1993). Statistics of local complexity in amino acid sequences and sequence databases. Computers & Chemistry 17(2):149–163. https://doi.org/10.1016/0097-8485(93)85006-X (search-confirmed via https://www.semanticscholar.org/paper/8c865ac68cd1b5f1ad2d69d8840ffce0e0f732ea, accessed 2026-06-14)
2. Wootton J.C., Federhen S. (1996). Analysis of compositionally biased regions in sequence databases. Methods in Enzymology 266:554–571. https://doi.org/10.1016/S0076-6879(96)66035-2 (cited in the GCG/NCBI SEG documentation, accessed 2026-06-14)
3. NCBI C++ Toolkit. `blast_seg.c` — SEG implementation (defaults `kSegWindow=12`, `kSegLocut=2.2`, `kSegHicut=2.5`; `s_Entropy`). https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_seg.c and https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html (accessed 2026-06-14)
4. SEG program help (GCG/Weizmann mirror). https://bip.weizmann.ac.il/education/materials/gcg/seg.html ; `ncbi-seg` manpage https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html ; rothlab genhelp https://rothlab.ucdavis.edu/genhelp/seg.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
