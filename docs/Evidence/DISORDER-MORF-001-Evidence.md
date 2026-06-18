# Evidence Artifact: DISORDER-MORF-001

**Test Unit ID:** DISORDER-MORF-001
**Algorithm:** MoRF (Molecular Recognition Feature) Prediction
**Date Collected:** 2026-06-14

---

## Online Sources

### Mohan et al. (2006) — "Analysis of molecular recognition features (MoRFs)"

**URL:** https://pubmed.ncbi.nlm.nih.gov/16935303/
**Accessed:** 2026-06-14 (retrieved via WebFetch on the PubMed record; discovered via WebSearch query
"Mohan 2006 Analysis of molecular recognition features MoRFs intrinsically disordered J Mol Biol 362 1043")
**Authority rank:** 1 (peer-reviewed primary paper, J Mol Biol)

**Key Extracted Points:**

1. **Length range:** MoRFs are "relatively short (10-70 residues), loosely structured protein regions"
   (extracted verbatim from the retrieved PubMed record).
2. **Embedding context:** MoRFs are "short ... protein regions within longer, largely disordered sequences"
   that "undergo disorder-to-order transitions" when binding partners.
3. **Structural classification:** three types by bound conformation — α-MoRFs (form α-helices),
   β-MoRFs (form β-strands), and ι-MoRFs (irregular secondary structure).

### Wikipedia — "Molecular recognition feature" (citing primaries)

**URL:** https://en.wikipedia.org/wiki/Molecular_recognition_feature
**Accessed:** 2026-06-14 (retrieved via WebFetch; discovered via WebSearch query
"MoRF molecular recognition feature definition short segment order within disordered region alpha-MoRF length")
**Authority rank:** 4 (Wikipedia citing primaries — the cited primary is Mohan 2006)

**Key Extracted Points:**

1. **Definition:** MoRFs are "small (10-70 residues) intrinsically disordered regions in proteins" that
   "undergo a disorder-to-order transition upon binding to their partners".
2. **Pre-binding state:** they are "disordered prior to binding to their partners, whereas they form a
   common 3D structure after interacting".
3. **Primary citation (cited by the article):** Mohan A, Oldfield CJ, Radivojac P, Vacic V, Cortese MS,
   Dunker AK, Uversky VN (October 2006). "Analysis of molecular recognition features (MoRFs)".
   J Mol Biol 362(5):1043–59. doi:10.1016/j.jmb.2006.07.087. PMID 16935303.

### Cheng/Oldfield et al. — "Mining α-helix-forming molecular recognition features (α-MoRFs) with cross species sequence alignments"

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2570644/
**Accessed:** 2026-06-14 (retrieved via WebFetch; discovered via the same MoRF-definition WebSearch)
**Authority rank:** 1 (peer-reviewed primary paper, Biochemistry)

**Key Extracted Points:**

1. **Operational definition (the "dip"):** the heuristic "identifies short regions of order within longer
   regions of disorder – or 'dips' – in disorder prediction profiles" (extracted verbatim).
2. **Core length:** an α-MoRF is "a short (around 20 residues) structural element"; candidate regions were
   retrieved as "regions of 30 residues or less".
3. **Disorder/order threshold:** the disorder prediction profile uses "the threshold of 0.5"; values
   "above the threshold of 0.5" indicate disorder and a "dip" below it indicates relative (predicted) order.
4. **Attribution of exact dip parameters:** the paper states its heuristic is "similar to the one used
   previously" and refers to Oldfield et al. (2005) for the original algorithmic details.

### Oldfield et al. (2005) — "Coupled folding and binding with α-helix-forming molecular recognition elements"

**URL:** https://pubmed.ncbi.nlm.nih.gov/16156658/
**Accessed:** 2026-06-14 (retrieved via WebFetch on the PubMed record; discovered via WebSearch query
"Oldfield 2005 Coupled folding binding alpha-helix molecular recognition elements Biochemistry PONDR dip")
**Authority rank:** 1 (peer-reviewed primary paper, Biochemistry)

**Key Extracted Points:**

1. **Element definition:** the MoRE "consists of a short region that undergoes coupled binding and folding
   within a longer region of disorder" (extracted verbatim from the retrieved abstract).
2. **Limitation of the public record:** the exact numeric dip-detection parameters (precise flank lengths,
   the exact ordered-run length window) live in the paywalled Methods section and could NOT be retrieved in
   this session — see Assumptions. The qualitative criterion (short ordered run flanked by disorder, 0.5
   threshold, 10–70 residue total length) is fully retrievable and is what this unit implements.

---

## Documented Corner Cases and Failure Modes

### From Cheng/Oldfield (PMC2570644) and Mohan (2006)

1. **Fully ordered protein:** a sequence predicted ordered throughout has no surrounding disorder, so it
   contains no "dip within disorder" → no MoRFs.
2. **Fully disordered protein:** a uniformly disordered sequence has no ordered dip → no MoRFs.
3. **Length filter:** regions outside the 10–70 residue band are not MoRFs (Mohan 2006 length range).
4. **Edge embedding:** a dip at the very start or end of the sequence is not flanked by disorder on both
   sides, so it does not satisfy the "within a longer region of disorder" requirement.

---

## Test Datasets

### Dataset: Synthetic dip-in-disorder construct (derived from the retrieved definition)

**Source:** Cheng/Oldfield (PMC2570644) "dips" criterion + Campen et al. (2008) TOP-IDP scale used by the
repository's `PredictDisorder` to produce per-residue disorder scores.

The repository disorder score is the normalized TOP-IDP propensity (higher = more disordered, range [0,1]).
A residue is **ordered** when score < 0.5 and **disordered** when score ≥ 0.5 (PMC2570644 0.5 threshold).
Per-residue TOP-IDP normalized values (`(prop − (−0.884)) / 1.871`, Campen Table 2) needed for the cases:

| Residue | TOP-IDP raw | Normalized = (raw + 0.884)/1.871 | Class at 0.5 |
|---------|-------------|----------------------------------|--------------|
| P | 0.987 | 1.000 | disordered |
| E | 0.736 | 0.866 | disordered |
| L | −0.326 | 0.298 | ordered |
| I | −0.486 | 0.213 | ordered |
| W | −0.884 | 0.000 | ordered |

A short homopolymer window of an ordered residue (e.g. L, score 0.298) embedded inside long P/E disordered
flanks yields a single ordered dip flanked by disorder = one MoRF. (Window averaging in `PredictDisorder`
smooths near boundaries; tests use flanks long enough that interior residues reach the pure per-residue
score.)

---

## Assumptions

1. **ASSUMPTION: Exact dip flank/length detection parameters** — Oldfield et al. (2005) defines the precise
   numeric dip parameters (flank length, ordered-run window) but the Methods section is paywalled and could
   not be retrieved. This unit therefore implements the fully-retrievable qualitative criterion: an ordered
   run (per-residue disorder score < 0.5) of total length within the Mohan 10–70 residue band, flanked on
   BOTH sides by at least one disordered residue inside a predicted disordered region. This is a
   correctness-affecting modeling choice for the flank-length detail only; the threshold (0.5), the length
   band (10–70), and the "order within disorder" shape are all source-traceable and are NOT assumptions.

---

## Recommendations for Test Coverage

1. **MUST Test:** An ordered dip (run of score < 0.5) flanked by disorder, of length within 10–70, is
   reported as exactly one MoRF at the dip's coordinates — Evidence: PMC2570644 "dips"; Mohan 2006 length.
2. **MUST Test:** A fully ordered sequence returns no MoRFs — Evidence: no surrounding disorder (PMC2570644).
3. **MUST Test:** A fully disordered sequence returns no MoRFs — Evidence: no ordered dip (PMC2570644).
4. **MUST Test:** An ordered dip shorter than the 10-residue minimum is not a MoRF — Evidence: Mohan 2006.
5. **MUST Test:** An ordered dip longer than the 70-residue maximum is not a MoRF — Evidence: Mohan 2006.
6. **MUST Test:** A dip at the sequence terminus (not flanked by disorder on both sides) is not a MoRF —
   Evidence: "within a longer region of disorder" (Oldfield 2005 / Mohan 2006).
7. **MUST Test:** MoRF score lies in [0,1] and increases with dip depth (distance below 0.5) — Evidence:
   0.5 threshold (PMC2570644); bounded normalization is a derivation, documented in the algorithm doc.
8. **SHOULD Test:** Two separate dips give two non-overlapping MoRFs — Rationale: independence of regions.
9. **SHOULD Test:** Case-insensitive input — Rationale: sibling methods upper-case input.
10. **COULD Test:** Null/empty input returns empty — Rationale: standard guard.

---

## References

1. Mohan A, Oldfield CJ, Radivojac P, Vacic V, Cortese MS, Dunker AK, Uversky VN. 2006. Analysis of
   molecular recognition features (MoRFs). J Mol Biol 362(5):1043–1059.
   https://doi.org/10.1016/j.jmb.2006.07.087 (PMID 16935303, https://pubmed.ncbi.nlm.nih.gov/16935303/)
2. Cheng Y, Oldfield CJ, Meng J, Romero P, Uversky VN, Dunker AK. (Mining α-helix-forming molecular
   recognition features with cross species sequence alignments.) Biochemistry.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC2570644/
3. Oldfield CJ, Cheng Y, Cortese MS, Romero P, Uversky VN, Dunker AK. 2005. Coupled folding and binding
   with α-helix-forming molecular recognition elements. Biochemistry 44(37):12454–12470.
   https://pubmed.ncbi.nlm.nih.gov/16156658/
4. Wikipedia contributors. Molecular recognition feature. https://en.wikipedia.org/wiki/Molecular_recognition_feature (accessed 2026-06-14)
5. Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK. 2008. TOP-IDP-Scale: A New Amino Acid
   Scale Measuring Propensity for Intrinsic Disorder. Protein Pept Lett 15(9):956–963.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ (used by `PredictDisorder` for per-residue scores)

---

## Change History

- **2026-06-14**: Initial documentation. Sources retrieved live via WebSearch/WebFetch in this session.
