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

## Scope / residual

This is a SMALL universal single-copy set, not CheckM's full lineage-specific marker database.
The full `checkm_data` (lineage-specific Pfam/TIGRFAM collocated marker sets + reference genome
tree for lineage placement) is a large gated trained DB and is NOT bundled; callers who have it
supply their own marker set via the caller-supplied loader. See the algorithm doc residual.
