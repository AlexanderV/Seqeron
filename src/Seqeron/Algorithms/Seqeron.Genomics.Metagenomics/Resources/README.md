# Bundled universal single-copy marker HMMs (META-BIN-001)

These nine files are HMMER3/f ASCII profile-HMM files for **universal bacterial single-copy
ribosomal protein** families. They constitute the small, redistributable default marker set used
by `MetagenomicsAnalyzer.EstimateBinQualityFromMarkers` to compute CheckM-style genome
completeness and contamination (Parks et al. 2015). Each profile is scored against a bin's
predicted proteins with the `Plan7ProfileHmm` engine (Viterbi log-odds vs the Pfam gathering
threshold GA1).

Ribosomal proteins are the canonical core of universal single-copy gene (USCG) sets: "the vast
majority of USCGs in every set are annotated with the functional category J (Translation,
ribosomal structure, and biogenesis)" (Xu et al. 2022, *A revisit to universal single-copy genes
in bacterial genomes*, PMC9411617).

| File | Pfam acc | Name | LENG | Description |
|------|----------|------|------|-------------|
| `PF00318_Ribosomal_S2.hmm`  | PF00318.26 | Ribosomal_S2  | 216 | Ribosomal protein S2 |
| `PF00177_Ribosomal_S7.hmm`  | PF00177.28 | Ribosomal_S7  | 148 | Ribosomal protein S7p/S5e |
| `PF00410_Ribosomal_S8.hmm`  | PF00410.25 | Ribosomal_S8  | 125 | Ribosomal protein S8 |
| `PF00380_Ribosomal_S9.hmm`  | PF00380.26 | Ribosomal_S9  | 121 | Ribosomal protein S9/S16 |
| `PF00338_Ribosomal_S10.hmm` | PF00338.28 | Ribosomal_S10 |  98 | Ribosomal protein S10p/S20e |
| `PF00411_Ribosomal_S11.hmm` | PF00411.25 | Ribosomal_S11 | 110 | Ribosomal protein S11 |
| `PF00203_Ribosomal_S19.hmm` | PF00203.27 | Ribosomal_S19 |  81 | Ribosomal protein S19 |
| `PF00687_Ribosomal_L1.hmm`  | PF00687.27 | Ribosomal_L1  | 198 | Ribosomal protein L1p/L10e family |
| `PF00297_Ribosomal_L3.hmm`  | PF00297.29 | Ribosomal_L3  | 369 | Ribosomal protein L3 |

## Provenance

- Retrieved 2026-06-25 via the EMBL-EBI InterPro web API, e.g.
  `https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00318/?annotation=hmm`
  (gzip of the HMMER3/f ASCII profile; decompressed verbatim, no edits).
- Built by the Pfam Consortium with HMMER `hmmbuild` (HMMER 3.3, Nov 2019 format).

## Licence — CC0 (public domain)

Pfam is freely available under the Creative Commons Zero ("CC0") licence
(InterPro/Pfam documentation: https://interpro-documentation.readthedocs.io/en/latest/pfam.html —
"Pfam is freely available under the Creative Commons Zero ('CC0') licence"). CC0 places the data
in the public domain, so it is freely redistributable; these files are embedded as-is.

## Domain-level GTDB universal marker sets (bac120 / ar122 Pfam subsets)

In addition to the 9-marker ribosomal set above, the package bundles the **Pfam-defined members
of the GTDB domain-level universal single-copy marker sets** — `bac120` (Bacteria) and `ar122`
(Archaea) — from Parks et al. (2018) *Nat Biotechnol* 36:996 and the GTDB-Tk reference data.
These feed `EstimateBinQualityFromMarkers` for routine domain-level completeness/contamination
out of the box (`LoadBundledBacterialMarkerHmms` / `LoadBundledArchaealMarkerHmms`).

- **Marker list source (verbatim):** the GTDB-Tk `bac120` / `ar122` `marker_info.tsv` files,
  Ecogenomics/GTDBTk repo (`tests/data/align_dir_reference/align/intermediate_results/gtdbtk.{bac120,ar122}.marker_info.tsv`),
  retrieved 2026-06-25 via the GitHub contents API. bac120 = 120 markers (6 Pfam + 114 TIGRFAM);
  ar122 = 122 markers (35 Pfam + 87 TIGRFAM).
- **What is bundled (CC0 only):** the **6 Pfam** members of bac120 (PF00380, PF00410, PF00466,
  PF01025, PF02576, PF03726) and the **35 Pfam** members of ar122 (PF00368, PF00410, PF00466,
  PF00687, PF00827, PF00900, PF01000, PF01015, PF01090, PF01092, PF01157, PF01191, PF01194,
  PF01198, PF01200, PF01269, PF01280, PF01282, PF01496, PF01655, PF01798, PF01864, PF01866,
  PF01868, PF01984, PF01990, PF02006, PF02978, PF03874, PF04019, PF04104, PF04919, PF07541,
  PF13656, PF13685). Distinct union across both = 39 Pfam HMMs (3 of them already in the 9-marker
  ribosomal set). Each HMM retrieved 2026-06-25 from the EMBL-EBI InterPro Pfam HMM API
  (`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<ACC>/?annotation=hmm`), HMMER3/f ASCII,
  licence **CC0**.
- **What is NOT bundled (licence):** the **TIGRFAM-defined** members of bac120 (114) and ar122 (87)
  are licensed **CC BY-SA 4.0** (TIGRFAMs at NCBI; reusabledata.org/tigrfams), *not* public domain.
  Per the redistribution policy they are NOT bundled; callers who have the TIGRFAM HMMs supply them
  via `LoadMarkerHmms`. This is why the bundled domain sets are the **Pfam subsets** of bac120/ar122.

## Scope / residual

These are domain-level (Bacteria / Archaea) universal single-copy sets, not CheckM's full
**lineage-specific** marker database. The full `checkm_data` (lineage-specific Pfam/TIGRFAM
collocated marker sets + reference genome tree for tree-based lineage placement, plus operon-based
marker-set collocation) is a large gated trained DB and is NOT bundled; callers who have it (or the
CC BY-SA TIGRFAM HMMs) supply their own marker set via the caller-supplied loader
(`LoadMarkerHmms`). See the algorithm doc residual.
