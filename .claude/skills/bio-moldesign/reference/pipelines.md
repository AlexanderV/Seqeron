# bio-moldesign — fuller pipelines & parameter guidance

Deeper recipes for the families in `SKILL.md`. Read the per-tool doc under
[`docs/mcp/tools/moltools/`](../../../../docs/mcp/tools/moltools) for exact I/O; compute with tools
only and attach a **Provenance** block (see [`bio-rigor`](../../bio-rigor/SKILL.md)). All positions
0-based; Tm in °C; codon GC as a fraction (0..1).

---

## 1. PCR primer pair design + QC

**Steps.** `design_primers`(template, target_start, target_end[, parameters]) →
`PrimerDesigner.DesignPrimers`. Then per primer QC:
- Tm: `primer_melting_temperature_salt`(primer, na_concentration) → `…CalculateMeltingTemperatureWithSalt`.
- Cross-dimer: `primer_dimer`(fwd, rev) → `…HasPrimerDimer`.
- Hairpin: `hairpin_potential`(primer) → `…HasHairpinPotential`.
- 3′ stability: `three_prime_stability`(primer) → `…Calculate3PrimeStability`.
- Low-complexity flags: `longest_homopolymer`, `longest_dinucleotide_repeat`.

**Parameters / gotchas.**
- `target_start`/`target_end` are **0-based inclusive**, `target_start < target_end < len`.
- Defaults when `parameters` is null: 18–25 nt, 40–60% GC, 57–63 °C Tm.
- The tool scans a **200 bp window** upstream of the target for fwd candidates and 200 bp downstream
  for rev (rev evaluated on the reverse complement); it returns the single best valid candidate each side.
- Pair `isValid` requires ΔTm ≤ 5 °C **and** no 3′ dimer. `productSize = rev.Position + rev.Length − fwd.Position`.
- Salt Tm adds `16.6·log10([Na+]/1000)` (Na+ in mM, default 50). Report the [Na+] you used.
- If no valid candidate exists, `forward`/`reverse` are null — report that, don't fabricate.

**Provenance example.** `design_primers(target=[100,150))` → `primer_melting_temperature_salt(Na=50)` ×2
→ `primer_dimer` → `hairpin_potential` ×2. Report GC%, Tm (°C, 1 dp), ΔTm, product size, dimer/hairpin flags.

---

## 2. CRISPR guide design → off-targets → specificity → rank

**Steps.** `design_guide_rnas`(sequence, region_start, region_end[, system_type, parameters]) →
`CrisprDesigner.DesignGuideRnas` returns guides ≥ `MinScore` (GC%, seed GC%, polyT, self-complementarity).
For each guide: `find_off_targets`(guide, genome[, max_mismatches, system_type]) →
`CrisprDesigner.FindOffTargets`; then `crispr_specificity_score`(guide, genome[, system_type]) →
`CrisprDesigner.CalculateSpecificityScore`. Rank by (on-target score, specificity).
Optional: `find_pam_sites`, `evaluate_guide_rna`, `crispr_system_info` for context.

**Parameters / gotchas.**
- Region is **0-based inclusive**; guide length must equal the system's guide length (else error 4001).
- Default system `SpCas9`. `crispr_system_info` reports PAM + guide length per system.
- `find_off_targets`: `max_mismatches` 0..5 (default 3); the exact on-target (0 mismatches) is excluded.
- Specificity: internally runs the scan at maxMismatches = 4; `100` = no off-targets, else
  `max(0, 100 − Σ penalty)`; **seed = PAM-proximal 12 bp** weighted 5 vs 2 per mismatch.
- **Envelope/perf:** naïve O(genome × guide); keep the scanned reference **≲ 1 Mb**. This is a screen,
  not a genome-wide off-target caller — say so, and add the not-for-clinical-use caveat for real assays.

**Provenance example.** `design_guide_rnas(region=[20,45], SpCas9)` → per guide
`find_off_targets(genome, max_mismatches=3)` → `crispr_specificity_score(genome)` → ranked table.

---

## 3. Codon optimization for a target organism

**Steps.** (optional) `find_rare_codons`(cds, organism, threshold) → inspect. Then
`optimize_codons`(coding_sequence, target_organism, strategy, gc_target_min, gc_target_max,
rare_codon_threshold) → `CodonOptimizer.OptimizeSequence`. It returns original/optimized RNA, protein,
**original/optimized CAI**, **original/optimized GC**, changedCodons, and each (position, original, optimized).
Standalone CAI: `codon_adaptation_index` (`CodonUsageAnalyzer.CalculateCai`) or `cai_from_organism_table`
(`CodonOptimizer.CalculateCAI`); usage stats via `rscu`, `effective_number_of_codons`, `codon_usage_statistics`.

**Parameters / gotchas.**
- `target_organism`: preset `EColiK12` / `Yeast` / `Human`, or an inline **RNA-alphabet** custom table.
- Input is upper-cased, T→U, trimmed to whole codons; **stop codons and Met/Trp are never changed**.
- Strategies: `MaximizeCAI` (deterministic, most-frequent codon), `BalancedOptimization` (default:
  most-frequent above rare threshold, then GC-balance), `AvoidRareCodeons`, `MinimizeSecondary`,
  `HarmonizeExpression` (**non-deterministic** weighted-random — repeat / fix a seed and report it).
- GC band default 0.40–0.60 (fractions in [0,1], min ≤ max); `rare_codon_threshold` default 0.15.
- Always report **CAI before→after** and **GC before→after** so the improvement is auditable.

**Provenance example.** `optimize_codons(preset=EColiK12, strategy=MaximizeCAI, GC=[0.40,0.60])`
→ report CAI 0.42→1.00, GC 0.38→0.52, changedCodons=N.

---

## 4. Restriction analysis (sites / digest)

**Steps.** Pick enzymes: `get_enzyme`, `compatible_enzymes`, `enzymes_compatible`, `blunt_cutters`,
`sticky_cutters`, `enzymes_by_cut_length`. Then `find_restriction_sites`(sequence, enzyme_names) →
`RestrictionAnalyzer.FindSites` (or `find_all_restriction_sites` for all known enzymes).
Simulate: `restriction_digest`(sequence, enzyme_names) → `RestrictionAnalyzer.Digest`;
summarize with `digest_summary`; overview map with `restriction_map`.

**Parameters / gotchas.**
- `enzyme_names` is a non-empty array. Digest is on a **linear** molecule; `k` distinct forward-strand
  cut positions → `k + 1` fragments; zero cuts → one whole-molecule fragment.
- Each fragment: `{sequence, startPosition, length, leftEnzyme, rightEnzyme, fragmentNumber}`;
  end fragments have a null flanking enzyme. Positions are **0-based**.
- For cloning compatibility, confirm overhangs with `compatible_enzymes` / `enzymes_compatible` before digesting.

**Provenance example.** `find_restriction_sites(["EcoRI"])` → `restriction_digest(["EcoRI"])`
→ fragment table (start, length, flanking enzymes).

---

## 5. Probe design + salt-corrected Tm  (⚠ PROBE-DESIGN-001)

**Steps.** `design_probes`(target_sequence, parameters=preset, max_probes) →
`ProbeDesigner.DesignProbes` (score-descending). Validate: `validate_probe`(`…ValidateProbe`) /
`analyze_oligo`(`…AnalyzeOligo`). Tm: `primer_melting_temperature_salt`(probe, na_concentration).
Variants: `design_tiling_probes`, `design_antisense_probes`, `design_molecular_beacon`.
Oligo prep: `oligo_extinction_coefficient`, `oligo_concentration_from_absorbance`.

**Parameters / gotchas.**
- Presets: `Microarray` (default), `FISH`, `NorthernBlot`, `qPCR`, `SouthernBlot`. `max_probes` > 0 (default 10).
- A target shorter than the preset's min probe length returns an **empty** list — report that.
- **STOP rule:** this unit (MinimumMode Moderate) does **not** compute quantitative **MGB ΔTm** or model
  dual-quencher labelling (empirical/proprietary, no closed form). If the assay needs an MGB / dual-quencher
  Tm, stop and report the limitation and point to a chemistry-specific tool — do not force a number.
  See [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) row PROBE-DESIGN-001
  and `bio-rigor` [`reference/envelope.md`](../../bio-rigor/reference/envelope.md).

**Provenance example.** `design_probes(preset=qPCR, max_probes=5)` → `validate_probe` →
`primer_melting_temperature_salt(Na=50)` → ranked probe table with Tm and warnings.

---

## Cross-checks (delegate detail to `bio-rigor` rule 4)

- Validate the input is real DNA/RNA before deriving anything (Sequence/Parsers validator tools).
- Corroborate a primer/probe Tm against its GC% and length; confirm a restriction cut position by
  re-finding the site; sanity-check optimized GC against the codon table's expectation.
- Never let a design number look more certain — or more clinical — than the alpha library warrants.
