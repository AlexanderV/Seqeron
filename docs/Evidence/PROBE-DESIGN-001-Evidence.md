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

---

## 2026-10-01 review (B07, F17–F21) — sources opened

- **primer3 `src/libprimer3.cc`** (raw.githubusercontent.com/primer3-org/primer3/main): internal-oligo defaults
  (`o_args`: opt/min/max size 20/18/27, Tm 60/57/63, GC 20–80, max_poly_x 5, salt 50 mM, divalent 0, dNTP 0,
  DNA 50 nM, max_self_any_th/self_end_th/hairpin_th 47), `pick_primer_range` (internal-oligo enumeration,
  no_internal_oligo_5_G / rev_comp_internal_oligo_g_c options default 0), `calc_and_check_oligo_features`
  (GC → poly-X → seqtm → thermodynamic self-any/self-end → hairpin), `p_obj_fn` OT_INTL branch, `primer_rec_comp`.
- **primer3-py 2.3.1** (installed): `design_primers(PRIMER_TASK=pick_hyb_probe_only)` — 950 random templates/settings:
  `DesignProbesPrimer3` identical (positions, order, Tm, penalty, SELF_ANY_TH, SELF_END_TH, HAIRPIN_TH; max |Δ| = 0);
  `calc_tm(…, mv 50, dv 0, dntp 0, dna 50)` for the probe Tm (CCATCACCCTACATCACC = 48.27871680775473 °C).
- **Biopython 1.88** `Bio.SeqUtils.molecular_weight` — 500 random DNA/RNA oligos identical to the routed
  `CalculateMolecularWeight` (e.g. ACGU RNA 1303.7737; old code 1239.8 with U = 308.2).
- **Applied Biosystems Primer Express TaqMan guidelines** (WebSearch extracts; Thermo/Fisher PDFs blocked):
  "the melting temperature (Tm) should be 68 to 70 °C when using Primer Express software", "Keep G-C content in the
  30-80% range", "runs of four or more Gs should be avoided", "as short as possible, without being shorter than 13 nucleotides".
- **OligoArray 2.0** (Rouillard, Zuker & Gulari 2003, NAR 31:3057; WebSearch extracts — academic.oup.com, PMC,
  ResearchGate and the umich OligoArray site blocked; audit round 3, A3-10): "This Tm is computed using the nearest
  neighbor (NN) model (10) with Na+ and DNA concentrations set to 1 M and 1 µM" (the PDF extract reads 1 mM; 1 µM taken,
  1 mM oligo is not a hybridization condition); A. thaliana design: length 45–47 nt, Tm range 82–90 °C; OligoArray 2.1
  installation example `-l 45 -L 47 -t 82 -T 88`. Microarray preset: 1 M / 1 µM, NN over the whole oligo, Tm 82–90 °C.
  Primer3 `oligotm.c` `seqtm` (raw.githubusercontent.com): `if (len > nn_max_len) long_seq_tm … else oligotm`.
  Cross-check: primer3-py 2.3.1 `calc_tm(AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAAT, mv 1000, dv 0, dntp 0,
  dna 1000, max_nn_length 60)` = 86.25468810488348 = Biopython `Tm_NN(DNA_NN3, Na 1000, dnac1 = dnac2 = 500,
  saltcorr 0)` 86.2547; old window max (60-mer, 60 % G+C, 50 mM, long_seq_tm) 74.5029020719779 < 75.
- **`DesignProbes` ranking (audit round 3, A3-11)**: Primer3 `src/libprimer3.cc` (raw.githubusercontent.com/primer3-org/primer3/main)
  `p_obj_fn` OT_INTL branch (`o_args.weights.temp_gt/temp_lt` × |Tm − opt_tm|, `length_lt/length_gt` × |len − opt_size|, …),
  `primer_rec_comp` (quality ↑, start ↓, length ↑), `_pr_data_control` ("PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE" /
  "< MIN_SIZE", "Optimum internal oligo Tm lower than minimum or higher than maximum"). The default additive score of
  `DesignProbes` has no published source (library heuristic); `ProbeRanking.Primer3Penalty` is the sourced option. Cross-check:
  primer3-py 2.3.1 `design_primers` (pick_hyb_probe_only, PRIMER_PICK_INTERNAL_OLIGO, size 18–27, Tm 57–63, GC 20–80, poly-X 5,
  PRIMER_NUM_RETURN 12) on 4 random 70-nt templates (random.seed(20261008)): 48/48 Primer3 probes are `DesignProbes` candidates
  with identical PRIMER_INTERNAL_n_PENALTY and Tm (|Δ| ≤ 1e-9, e.g. 27/22 0.007361159225695246, 14/19 0.007730506722964492)
  and order (ties 16/20 before 15/20 at 2.374406691313709, 22/24 before 21/24 at 2.6378890608341408).
- **Kane et al. 2000** (NAR 28:4552; WebSearch extracts): 50-mer probes; > 75 % identity or a ≥ 15-nt contiguous stretch to
  non-targets may cross-hybridize (specificity rule — PROBE-VALID-001 scope).
- **Molecular beacons** (Tyagi & Kramer 1996; Marras/Vet design protocol; WebSearch extracts, publisher pages blocked):
  loop 15/18–30 nt, stem 5–7 bp, stem Tm 7–10 °C above the detection temperature, probe Tm 7–10 °C above the PCR annealing temperature.
- **ε260 nearest-neighbour tables** (Cantor, Warshaw & Shapiro 1970; Warshaw & Tinoco 1966; WebSearch extracts of
  vendor/ATDBio tables — ATDBio/TriLink pages blocked): DNA ApA 27400, ApC 21200, ApG 25000, CpC 14600, CpG 18000,
  CpT 15200 …; RNA ApA 27400, ApC 21000, ApG 25000, ApU 24000, CpC 14200, CpG 17800, CpU 16200 …

## 2026-10-08 review (B07, F59, audit round 4 A4-1) — genome-index overload walks all candidates

- **Contract source opened:** `ProbeDesigner.DesignProbes(string, ISuffixTree, ProbeParameters?, int, bool)` XML
  ("Maximum number of probes to return", "If true, only return probes unique in the genome") and
  `docs/algorithms/MolTools/Hybridization_Probe_Design.md` §3.1/§3.3/§5.3 (which admitted the `maxProbes * 5`
  shortlist as an intentional simplification). No external reference defines this library overload; the contract is
  the specification, so the shortlist was a defect (unique probes that exist are not returned) and the missing
  re-rank after the `Score × specificity` scaling broke the documented score-descending order.
- **Repro fixture (auditor):** `genome = X + T10 + X + T10 + Y`, `target = X + Y` with X (400 nt) and Y (80 nt) drawn
  from `new Random(7)` as uniform ACGT; Microarray preset, `maxProbes = 2`, `requireUnique = true` → 0 probes before
  the fix, the unique 50-mers at target 402 and 417 (additive score 0.85) after it.
- **Independent cross-check (Python 3, `str.find` overlapping-occurrence count over the same two strings dumped from
  the C# fixture):** of the 4686 candidate windows (lengths 50–60) 3806 occur exactly twice in the genome (inside the
  duplicated X), 583 not at all (crossing the X|Y junction) and **297 occur exactly once** (all starting at ≥ 399, i.e.
  inside Y) — so unique probes exist although none of them is in the top `maxProbes * 5 = 10` raw-score candidates.
  The 50-mers at 402 / 417 have G+C 0.50 / 0.46 (inside the preset's 0.40–0.60) and occurrence count 1.
  In-test oracle: the exhaustive non-index `DesignProbes(target, param, int.MaxValue)` candidate list (the documented
  ranking) filtered / scaled with that naive occurrence count reproduces the overload's output for
  `maxProbes` ∈ {1, 3, 100000} and both `requireUnique` values, and for the `Primer3Penalty` ranking (F51).

## 2026-10-08 review (B07, F60, audit round 4 A4-2) — genome-index overload on both strands (opt-in)

- **Source opened:** BLAST+ 2.12.0 `blastn -help` (`-strand` default `both`); `CheckSpecificity(…, bothStrands)` contract (F24: reverse-palindromic probe counted once).
- **Oracle:** Python `str.find` overlapping-occurrence counts + Biopython `reverse_complement` over the dumped fixture strings (X, Y from C# `System.Random(7)`; index X + T10 + Y + T10 + rc(X), target X + Y): of the 4686 candidate windows (50–60 nt) 3806 have (indexed-strand, both-strand) sites (1, 2), 297 (1, 1) — all inside Y —, 561 (0, 0) and 22 (0, 1) (X|Y-crossing windows with Y's first 1–2 nt `AA`, whose reverse complement `TT` lies in the T10 before rc(X)).
- **Result:** default unchanged (indexed strand); `bothStrands = true`, Microarray, `maxProbes = 2`, `requireUnique` → 351 (X|Y-crossing window `GTCTTTGATCGCCACGACCTACCCCTAACAAACATTAGTCCAATTTAATA`, G+C 0.40, Python 0 / 1 sites — present only as its reverse complement; score 0.98) / 402 (0.85); brute force over all candidates equals the overload for `requireUnique` × `maxProbes` ∈ {3, 100000}.

## 2026-10-08 review (B07, F67, audit round 7 A7-2) — tiling covers the 3′ tail; argument guards

- **Source opened:** CATCH v1.5.2 (Metsky HC et al. 2019, Nat Biotechnol 37:160–168; git tag `v1.5.2` of
  github.com/broadinstitute/catch, cloned 2026-10-08) `catch/filter/candidate_probes.py`
  `make_candidate_probes_from_sequence(seq, probe_length, probe_stride, …)`:
  - grid: `for start in np.arange(0, len(seq), probe_stride): if start + probe_length > len(seq): break` — the same
    grid as `DesignTilingProbes` (`step = probeLength − overlap`);
  - end anchor: `if len(seq) % probe_stride != 0: # There are bases on the right that were never covered, so add
    another probe for this` → `add_probe_from_subsequence(len(seq) - probe_length, len(seq))`;
  - short sequence: `if len(seq) < probe_length:` → `ValueError("An input sequence is smaller than the probe length …")`
    unless `allow_small_seqs`.
- **Adopted:** the end-anchored final window (start `L − probeLength`) and rejection of a target shorter than the probe.
  **Not adopted (documented):** CATCH's trigger `len(seq) mod stride ≠ 0` — it does not test what its comment states:
  100 nt / 60 / stride 40 adds start 40 a second time (CATCH notes duplicates are possible), 120 nt / 60 / stride 40
  (120 mod 40 = 0) leaves 100–119 uncovered. The library adds the window exactly when the last grid window ends before
  the target end, `(L − probeLength) mod step ≠ 0`.
- **Repros (auditor):** 30-nt target, default 60 → `InvalidOperationException` ("Sequence contains no elements") also
  through MCP; `probeLength ≤ 0` or `overlap ≥ probeLength` → infinite loop (MCP guarded); 110 nt / 60 / 20 → starts
  0, 40 only, positions 100–109 uncovered while the doc claimed full coverage. After: `ArgumentOutOfRangeException`,
  `ArgumentOutOfRangeException`, starts 0, 40, 50 with Coverage 110.
- **Oracle:** tests compute the expected starts by hand from the rule above and Coverage independently (boolean
  union over `[start, start + P)`), for 9 (L, P, overlap) cases incl. both CATCH corner cases and negative overlap.
