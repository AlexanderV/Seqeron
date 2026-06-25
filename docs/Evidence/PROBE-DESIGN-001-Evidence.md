# Evidence Artifact: PROBE-DESIGN-001

**Test Unit ID:** PROBE-DESIGN-001
**Algorithm:** Hybridization Probe Design — TaqMan (5'-nuclease hydrolysis probe) rules (opt-in extension)
**Date Collected:** 2026-06-24

---

## Online Sources

This Evidence file covers the **TaqMan-specific design rules** added as an opt-in mode.
The generic hybridization-probe designer (length/Tm/GC windows, homopolymer, self-complementarity,
hairpin, specificity) remains the unchanged default and is documented in the algorithm doc
(§5) and the original validation report; its thermodynamic Tm formulas were validated under
PRIMER-TM-001.

### PREMIER Biosoft — "TaqMan probe design tips"

**URL:** http://www.premierbiosoft.com/tech_notes/TaqMan.html
**Accessed:** 2026-06-24 (fetched via WebFetch)
**Authority rank:** 3 (established vendor design-tool documentation, restating Applied Biosystems guidance)

**Key Extracted Points (verbatim quotes from the retrieved page):**

1. **Probe length:** "TaqMan® probes consist of a 18-22 bp oligonucleotide probe".
2. **G+C content:** "The G+C content should ideally be 30-80%".
3. **No 5' G + more C than G:** "there should be more Cs than Gs, and not a G at the 5' end".
4. **No run of ≥4 Gs:** "The probes should not have runs of identical nucleotides (especially four or more consecutive Gs)".
5. **Probe Tm vs primer Tm:** "TaqMan® probe Tm should be 10 oC higher than the Primer Tm".

### Applied Biosystems / Thermo Fisher — TaqMan assay design guidance

**URL:** https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/real-time-pcr-learning-center/gene-expression-analysis-real-time-pcr-information/designing-taqman-gene-expression-assay.html
**Accessed:** 2026-06-24 (located via WebSearch; the page redirect-loops to WebFetch, but the
guideline text was returned in the search-result extract and corroborates PREMIER Biosoft point-for-point)
**Authority rank:** 2 (manufacturer of the TaqMan chemistry; canonical design guidance)

**Key Extracted Points:**

1. **No G at the 5' end:** "Probes do not have a G at the 5'-end, because this nucleotide
   interferes with fluorescence from the reporter dye that will be attached to this end."
2. **Probe Tm:** "the Tm of the probe should be about 10°C above that of the primer to enable
   the probe to bind to the template strand before Taq polymerase reaches it."
3. **Strand fallback:** "if it is not possible to select a probe without a guanine residue at the
   5' end, you will need to design a probe on the complement (antisense) strand."

### TaqMan chemistry mechanism (5'-nuclease quenching rationale)

**URL:** https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/taqman
**Accessed:** 2026-06-24 (page is 403 to WebFetch; the quoted sentence was returned in the
WebSearch result extract)
**Authority rank:** 1–3 (peer-reviewed reference-work topic page)

**Key Extracted Points:**

1. **5'-G quenching rationale:** "there be no G at the 5′ end, as a 'G' adjacent to the reporter
   dye quenches reporter fluorescence even after cleavage." This is why the 5'-G rule cannot be
   rescued by cleavage and is treated as a hard rule.

### IDT — "Designing PCR primers and probes" (corroborating)

**URL:** https://www.idtdna.com/pages/education/decoded/article/designing-pcr-primers-and-probes
**Accessed:** 2026-06-24 (page redirect-loops to WebFetch; guideline values returned in WebSearch extract)
**Authority rank:** 3 (oligo-manufacturer design guidance)

**Key Extracted Points:**

1. Probe Tm ≈ 10 °C higher than primer Tm; primer Tm ≈ 58–60 °C.
2. Probe length 18–30 nt; more Cs than Gs; no G at the 5' end; runs of identical nucleotides ≤ 4.
   (The tighter 18–22 nt length window from Applied Biosystems / PREMIER Biosoft is used as the default.)

---

## Documented Corner Cases and Failure Modes

### From Applied Biosystems / Thermo Fisher

1. **5' G unavoidable on the sense strand:** design the probe on the complement (antisense)
   strand instead. Implemented by `SelectTaqManStrand`.

### From PREMIER Biosoft

1. **Run of ≥4 consecutive Gs:** explicitly called out as the worst homopolymer case; flagged
   separately from the generic homopolymer cap.

---

## Test Datasets

### Dataset: Hand-derived TaqMan rule examples (this unit)

**Source:** Derived from the rules above; Tm via the repository's salt-adjusted formula
`81.5 + 16.6·log₁₀[Na⁺] + 41·GC − 600/N` ([Na⁺] = 0.05 M), validated under PRIMER-TM-001.

| Probe (5'→3') | Len | 5' base | C | G | maxGrun | GC | Tm (°C) | Expected |
|---------------|-----|---------|---|---|---------|------|---------|----------|
| `CCATCACCCTACATCACC` | 18 | C | 10 | 0 | 0 | 0.5556 | 49.3473 | passes all (primerTm ≤ 39.35) |
| `GCATCACCCTACATCACC` | 18 | G | 9 | 1 | 1 | 0.5556 | 49.3473 | fails: 5'-G |
| `ACCCCGGGGACCCTACAT` | 18 | A | 8 | 4 | 4 | 0.6667 | — | fails: GGGG run |
| `ACGGGAGGTAGGTAGGTA` | 18 | A | 1 | 9 | 3 | — | — | fails: more G than C |
| `CCATCACCCTACATCA`   | 16 | C | — | — | — | — | — | fails: length < 18 |
| `CCCGCCCCGCCCCGCCCC` | 18 | C | 15 | 3 | 1 | 1.0000 | — | fails: GC = 100% |
| sense `GTTAGGGTTAGGGTTAGG` → RC `CCTAACCCTAACCCTAAC` | 18 | G→C | 0→9 | 9→0 | 3→0 | 0.50 | — | strand selection picks antisense |

---

## Assumptions

1. **ASSUMPTION: 18–22 nt default length window.** Applied Biosystems / PREMIER Biosoft give
   18–22; IDT/Thermo elsewhere allow up to 30. We default to the tighter 18–22 and expose
   `minLength`/`maxLength` parameters, so this does not hard-code an un-citable value.
2. **ASSUMPTION: Probe-Tm gate uses the repository salt-adjusted Tm.** The "+10 °C above primer"
   rule is sourced; the exact Tm engine is the repository's existing (PRIMER-TM-001-validated)
   formula rather than a TaqMan-specific nearest-neighbor calc. The caller supplies the primer Tm.

---

## Recommendations for Test Coverage

1. **MUST Test:** a probe with a 5' G is flagged (`NoGuanineAt5Prime == false`) and `PassesAll == false`
   — Evidence: Applied Biosystems / ScienceDirect (5'-G quenches reporter even after cleavage).
2. **MUST Test:** the more-C-than-G rule (C=1, G=9 example fails) — Evidence: PREMIER Biosoft / ABI.
3. **MUST Test:** run of ≥4 Gs flagged — Evidence: PREMIER Biosoft.
4. **MUST Test:** GC outside 30–80% flagged; length outside 18–22 flagged — Evidence: PREMIER Biosoft.
5. **MUST Test:** probe-Tm gate vs supplied primer Tm (fails when Tm < primerTm + 10) — Evidence: PREMIER Biosoft / ABI.
6. **MUST Test:** a fully compliant probe is accepted (`PassesAll == true`, no violations).
7. **MUST Test:** strand selection picks the more-C-than-G / no-5'-G strand on the known example — Evidence: ABI antisense fallback.

---

## References

1. PREMIER Biosoft. "TaqMan® Probes | TaqMan® probe design tips." http://www.premierbiosoft.com/tech_notes/TaqMan.html (accessed 2026-06-24).
2. Applied Biosystems / Thermo Fisher Scientific. "Designing a TaqMan Gene Expression Assay." https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/real-time-pcr-learning-center/gene-expression-analysis-real-time-pcr-information/designing-taqman-gene-expression-assay.html (accessed 2026-06-24).
3. ScienceDirect Topics. "TaqMan — an overview." https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/taqman (accessed 2026-06-24).
4. Integrated DNA Technologies. "Rules and Tips for PCR & qPCR Primer Design." https://www.idtdna.com/pages/education/decoded/article/designing-pcr-primers-and-probes (accessed 2026-06-24).

---

## Change History

- **2026-06-24**: Initial documentation — TaqMan opt-in rules (no 5'-G, C>G strand, ≥4-G run, GC 30–80%, length 18–22, probe Tm ≥ primer Tm + 10 °C).
