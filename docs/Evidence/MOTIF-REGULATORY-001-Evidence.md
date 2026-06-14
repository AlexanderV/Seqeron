# Evidence Artifact: MOTIF-REGULATORY-001

**Test Unit ID:** MOTIF-REGULATORY-001
**Algorithm:** Regulatory Elements (scan a DNA sequence for known regulatory consensus motifs)
**Date Collected:** 2026-06-14

---

## Online Sources

### TATA box — Wikipedia (citing primaries) + Bucher (1990)

**URL:** https://en.wikipedia.org/wiki/TATA_box (search: "TATA box consensus sequence TATAAA TATA-binding protein core promoter")
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing primaries) / 1 (Bucher 1990, J Mol Biol)

**Key Extracted Points:**

1. **Consensus:** "The TATA box has the consensus sequence of 5'-TATAAA-3'." A wider IUPAC consensus is `TATA(A/T)A(A/T)(A/G)`.
2. **Location:** Core promoter, ~25–35 bp upstream of the transcription start site (RNA Pol II).
3. **Bucher weight-matrix:** Bucher (1990, J Mol Biol) derived weight-matrix descriptions of four eukaryotic Pol II promoter elements from 502 unrelated promoters; the TATA box is the most prominent.

### Pribnow box (-10) and -35 box — Wikipedia + Harley & Reynolds (1987)

**URL:** https://en.wikipedia.org/wiki/Pribnow_box (search: "prokaryotic promoter consensus -10 TATAAT -35 TTGACA"); Google Scholar lookup for Harley & Reynolds 1987 (search: "-35 box consensus TTGACA sigma 70 promoter Harley Reynolds 1987")
**Accessed:** 2026-06-14
**Authority rank:** 1 (Harley & Reynolds 1987, Nucleic Acids Res) / 4 (Wikipedia)

**Key Extracted Points:**

1. **-10 consensus:** "The Pribnow box ... is a sequence of TATAAT of six nucleotides" — the -10 hexamer.
2. **-35 consensus:** Harley & Reynolds (1987): "all bases in the -35 (TTGACA) and -10 (TATAAT) hexamers were highly conserved"; 92% of promoters had inter-region spacing of 17±1 bp.
3. **Citation:** Harley C.B., Reynolds R.P. (1987). "Analysis of E. coli promoter sequences." Nucleic Acids Res 15(5):2343-2361. DOI 10.1093/nar/15.5.2343, PMID 3550697.

### CCAAT box — Bucher (1990) via NAR survey

**URL:** https://academic.oup.com/nar/article/26/5/1135/1282635 (search: "CCAAT box consensus CAAT box promoter NF-Y Bucher 1990")
**Accessed:** 2026-06-14
**Authority rank:** 1 (Bucher 1990, J Mol Biol)

**Key Extracted Points:**

1. **Pentanucleotide:** "the CCAAT pentanucleotide is present in approximately 30% of promoters" (Bucher 1990). Core element = `CCAAT`.
2. **Extended motif:** flanking-extended consensus `RRCCAATSR` (Bucher 1990); core five-nucleotide element is `CCAAT`.

### GC box (Sp1) — Wikipedia (citing Lundin et al. 1994)

**URL:** https://en.wikipedia.org/wiki/GC_box (search: "GC box Sp1 binding site consensus GGGCGG")
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia) / 1 (Lundin, Nehlin & Ronne 1994, Mol Cell Biol)

**Key Extracted Points:**

1. **Consensus:** "The GC box has a consensus sequence GGGCGG which is position-dependent and orientation-independent." Primary: Lundin M., Nehlin J.O., Ronne H. (1994) Mol Cell Biol 14(3):1979-1985.
2. **Binding:** Sp1 zinc fingers bind `(G/T)GGGCGG(G/A)(G/A)(C/T)`; the core element is `GGGCGG`.

### Kozak consensus — Wikipedia (citing Kozak 1987)

**URL:** https://en.wikipedia.org/wiki/Kozak_consensus_sequence (fetched)
**Accessed:** 2026-06-14
**Authority rank:** 1 (Kozak 1987, Nucleic Acids Res)

**Key Extracted Points:**

1. **Canonical sequence:** `5'-gccRccAUGG-3'` (R = A/G, A more frequent); determined by sequencing 699 vertebrate mRNAs. As the most-preferred-base DNA string this is `GCCGCCACCATGG` (positions -9..+4, ATG = +1..+3).
2. **Key positions:** -3 must be a purine (A/G); +4 must be G for a strong context.
3. **Citation:** Kozak M. (1987). "An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs." Nucleic Acids Res 15(20):8125-8148.

### Shine-Dalgarno (RBS) — Wikipedia / ScienceDirect

**URL:** https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence (search: "Shine-Dalgarno consensus sequence AGGAGG 16S rRNA ribosome binding site")
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia) / 1 (Shine & Dalgarno 1974, PNAS)

**Key Extracted Points:**

1. **Consensus:** "the six-base consensus sequence is AGGAGG"; located ~8 bases upstream of the AUG start codon in bacterial/archaeal mRNA.
2. **Mechanism:** complementary `CCUCCU` at the 3' end of 16S rRNA base-pairs with the SD sequence to initiate translation.

### Poly(A) signal — Proudfoot & Brownlee (1976)

**URL:** https://genesdev.cshlp.org/content/25/17/1770.full (search: "poly(A) signal AATAAA polyadenylation Proudfoot Brownlee 1976 consensus hexamer")
**Accessed:** 2026-06-14
**Authority rank:** 1 (Proudfoot & Brownlee 1976, Nature)

**Key Extracted Points:**

1. **Consensus:** "the principal cis-acting motif ... is a hexanucleotide A(A/U)UAAA sequence (Proudfoot and Brownlee, 1976)"; the canonical hexamer is `AAUAAA` (DNA `AATAAA`).
2. **Location:** ~10–30 nt upstream of the cleavage site.
3. **Citation:** Proudfoot N.J., Brownlee G.G. (1976). "3' non-coding region sequences in eukaryotic messenger RNA." Nature 263:211-214.

### E-box — Wikipedia (citing Massari & Murre 2000)

**URL:** https://en.wikipedia.org/wiki/E-box (fetched)
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia) / 1 (Massari & Murre 2000, Mol Cell Biol)

**Key Extracted Points:**

1. **Consensus:** "The E-box is a DNA sequence CANNTG (where N can be any nucleotide)"; canonical palindrome `CACGTG`. We use the IUPAC family pattern `CANNTG` so any central dinucleotide matches.
2. **Citation:** Massari M.E., Murre C. (2000). "Helix-loop-helix proteins: regulators of transcription in eucaryotic organisms." Mol Cell Biol 20(2):429-440.

### AP-1 (TRE) — Lee, Mitchell & Tjian (1987)

**URL:** https://pubmed.ncbi.nlm.nih.gov/3034433/ (fetched)
**Accessed:** 2026-06-14
**Authority rank:** 1 (Lee, Mitchell & Tjian 1987, Cell)

**Key Extracted Points:**

1. **Recognition motif:** "high-affinity AP-1-binding sites, each with a conserved recognition motif, TGACTCA." Consensus palindrome `TGA(C/G)TCA`; canonical `TGACTCA`.
2. **Citation:** Lee W., Mitchell P., Tjian R. (1987). "Purified transcription factor AP-1 interacts with TPA-inducible enhancer elements." Cell 49(6):741-752.
3. **DEFECT FOUND:** prior repository code used `TGAGTCA` (G at position 4) — wrong; corrected to `TGACTCA` per this primary.

### NF-κB κB site — Sen & Baltimore (1986)

**URL:** https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101490 (search: "NF-kB consensus binding site GGGACTTTCC GGGRNNYYCC Sen Baltimore 1986")
**Accessed:** 2026-06-14
**Authority rank:** 1 (Sen & Baltimore 1986, Cell) / 1 (PLoS One mutation-matrix study)

**Key Extracted Points:**

1. **Consensus:** κB sites have consensus `GGGRNWYYCC` (R=purine, W=A/T, Y=pyrimidine, N=any).
2. **Reference site:** "the wild-type binding site GGGACTTTCC" used in p50-binding studies — a canonical strong κB site; we scan for `GGGACTTTCC`.
3. **Citation:** Sen R., Baltimore D. (1986). "Multiple nuclear factors interact with the immunoglobulin enhancer sequences." Cell 46(5):705-716.

### CREB CRE — Montminy et al. (1986)

**URL:** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1796609/ (search: "CREB binding site CRE consensus TGACGTCA cyclic AMP response element Montminy 1986")
**Accessed:** 2026-06-14
**Authority rank:** 1 (Montminy et al. 1986, PNAS)

**Key Extracted Points:**

1. **Consensus:** "a consensus palindrome—TGACGTCA—that recognizes the CREB dimer with high affinity, as reported by Montminy et al. in 1986." First identified in the somatostatin gene.
2. **Citation:** Montminy M.R. et al. (1986). PNAS 83:6682-6686.

---

## Documented Corner Cases and Failure Modes

### From the consensus definitions

1. **Degenerate codes (E-box):** `CANNTG` contains IUPAC `N`; any base matches the two central positions, so `CACGTG`, `CAGCTG`, `CATATG`, etc. all match.
2. **Overlapping / multiple occurrences:** a sequence may contain several occurrences of one element and occurrences of several different elements; each occurrence is reported with its 0-based start position.
3. **Substring containment:** because `TATAAT` (-10 box) and `TATAAA` (TATA box) differ only in the final base, a real promoter region scanned will report whichever exact hexamer(s) are present; both are distinct entries.
4. **Empty sequence:** no occurrences → empty result (scanning window `0 <= i <= n - m`).
5. **Null sequence:** invalid input → `ArgumentNullException`.

---

## Test Datasets

### Dataset: Hand-constructed exact-match probes (one per element)

**Source:** consensus strings copied verbatim from the primary sources above.

| Element | Consensus (5'→3') | Probe sequence | Expected position(s) |
|---------|-------------------|----------------|----------------------|
| TATA Box | TATAAA | `GGGTATAAAGGG` | 3 |
| -10 Box | TATAAT | `CCTATAATCC` | 2 |
| -35 Box | TTGACA | `AATTGACAGG` | 2 |
| CAAT Box | CCAAT | `GGCCAATGG` | 2 |
| GC Box | GGGCGG | `AAGGGCGGTT` | 2 |
| Kozak | GCCGCCACCATGG | `TTGCCGCCACCATGGAA` | 2 |
| Shine-Dalgarno | AGGAGG | `TTAGGAGGTTT` | 2 |
| Poly(A) Signal | AATAAA | `CCAATAAACC` | 2 |
| E-box | CANNTG (CACGTG) | `GGCACGTGGG` | 2 |
| AP-1 | TGACTCA | `AATGACTCAGG` | 2 |
| NF-κB | GGGACTTTCC | `AAGGGACTTTCCAA` | 2 |
| CREB | TGACGTCA | `CCTGACGTCAGG` | 2 |

### Dataset: AP-1 negative control (defect regression)

**Source:** Lee, Mitchell & Tjian (1987).

| Parameter | Value |
|-----------|-------|
| Sequence | `AATGAGTCAGG` (contains the old wrong pattern TGAGTCA) |
| Expected AP-1 hits | 0 (TGAGTCA is NOT the consensus; correct consensus is TGACTCA) |

---

## Assumptions

1. **ASSUMPTION: NF-κB reference string.** The general κB consensus is `GGGRNWYYCC`. We scan for the specific strong reference site `GGGACTTTCC` (used in Sen & Baltimore 1986 and as the wild-type in p50-binding studies) rather than expanding the degenerate consensus, matching the established repository behavior. This is a representative-site choice, not a fabricated value; the exact string is source-backed.
2. **ASSUMPTION: Kozak as exact DNA string.** Kozak (1987) defines `gccRccATGG`. We scan for the single most-preferred-base string `GCCGCCACCATGG`. The position -3 purine and +4 G degeneracy is not expanded; only the canonical exact context is detected.

---

## Recommendations for Test Coverage

1. **MUST Test:** each of the 12 elements is detected at the exact 0-based position in its probe, with the correct `Name`, `Pattern`, and `Sequence`. — Evidence: primary consensus strings above.
2. **MUST Test:** AP-1 corrected consensus `TGACTCA` matches; old `TGAGTCA` does NOT produce an AP-1 hit (regression of the corrected defect). — Evidence: Lee, Mitchell & Tjian (1987).
3. **MUST Test:** E-box `CANNTG` matches a degenerate occurrence (`CACGTG`). — Evidence: Massari & Murre (2000).
4. **MUST Test:** null sequence → `ArgumentNullException`; empty sequence → no elements. — Evidence: contract / scanning window.
5. **SHOULD Test:** multiple occurrences of one element reported at all positions; multiple distinct elements in one sequence. — Rationale: documented corner case.
6. **COULD Test:** `KnownMotifs` constants equal their source consensus strings. — Rationale: locks the constants to the cited values.

---

## References

1. Bucher P. (1990). Weight matrix descriptions of four eukaryotic RNA polymerase II promoter elements derived from 502 unrelated promoter sequences. J Mol Biol 212(4):563-578. https://doi.org/10.1016/0022-2836(90)90223-9
2. Harley C.B., Reynolds R.P. (1987). Analysis of E. coli promoter sequences. Nucleic Acids Res 15(5):2343-2361. https://doi.org/10.1093/nar/15.5.2343
3. Lundin M., Nehlin J.O., Ronne H. (1994). Importance of a flanking AT-rich region in target site recognition by the GC box-binding zinc finger protein MIG1. Mol Cell Biol 14(3):1979-1985. https://doi.org/10.1128/mcb.14.3.1979
4. Kozak M. (1987). An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs. Nucleic Acids Res 15(20):8125-8148. https://doi.org/10.1093/nar/15.20.8125
5. Proudfoot N.J., Brownlee G.G. (1976). 3' non-coding region sequences in eukaryotic messenger RNA. Nature 263:211-214. https://doi.org/10.1038/263211a0
6. Massari M.E., Murre C. (2000). Helix-loop-helix proteins: regulators of transcription in eucaryotic organisms. Mol Cell Biol 20(2):429-440. https://doi.org/10.1128/MCB.20.2.429-440.2000
7. Lee W., Mitchell P., Tjian R. (1987). Purified transcription factor AP-1 interacts with TPA-inducible enhancer elements. Cell 49(6):741-752. https://doi.org/10.1016/0092-8674(87)90612-X (PMID 3034433: https://pubmed.ncbi.nlm.nih.gov/3034433/)
8. Sen R., Baltimore D. (1986). Multiple nuclear factors interact with the immunoglobulin enhancer sequences. Cell 46(5):705-716. https://doi.org/10.1016/0092-8674(86)90346-6
9. Montminy M.R., Sevarino K.A., Wagner J.A., Mandel G., Goodman R.H. (1986). Identification of a cyclic-AMP-responsive element within the rat somatostatin gene. PNAS 83(18):6682-6686. https://doi.org/10.1073/pnas.83.18.6682
10. Wikipedia: TATA box — https://en.wikipedia.org/wiki/TATA_box ; Pribnow box — https://en.wikipedia.org/wiki/Pribnow_box ; GC box — https://en.wikipedia.org/wiki/GC_box ; E-box — https://en.wikipedia.org/wiki/E-box ; Shine–Dalgarno sequence — https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence ; Kozak consensus sequence — https://en.wikipedia.org/wiki/Kozak_consensus_sequence (all accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation. Recorded AP-1 defect (TGAGTCA → TGACTCA) and addition of -10/-35 prokaryotic promoter hexamers.
