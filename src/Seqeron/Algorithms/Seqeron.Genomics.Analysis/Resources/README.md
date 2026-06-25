# Bundled Pfam profile HMMs (PROTMOTIF-DOMAIN-001)

These three files are HMMER3/f ASCII profile-HMM files for protein domains that have **no
deterministic PROSITE pattern** (they are trained profile HMMs). They are consumed by the
`Plan7ProfileHmm` engine (`ProteinMotifFinder.FindDomainsByHmm`).

| File | Pfam acc | Name | LENG (match states) |
|------|----------|------|---------------------|
| `PF00018_SH3_1.hmm` | PF00018.35 | SH3_1 (SH3 domain) | 48 |
| `PF00595_PDZ.hmm`   | PF00595.30 | PDZ (PDZ domain)   | 81 |
| `PF00400_WD40.hmm`  | PF00400.39 | WD40 (WD domain, G-beta repeat) | 39 |

## Provenance

- Retrieved 2026-06-25 via the EMBL-EBI InterPro web API, e.g.
  `https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00018/?annotation=hmm`
  (gzip of the HMMER3/f ASCII profile; decompressed verbatim, no edits).
- Built by the Pfam Consortium with HMMER `hmmbuild` (HMMER 3.3, Nov 2019 format).

## Licence — CC0 (public domain)

Pfam is freely available under the Creative Commons Zero ("CC0") licence
(InterPro/Pfam documentation: https://interpro-documentation.readthedocs.io/en/latest/pfam.html).
CC0 places the data in the public domain, so it is freely redistributable; these files are
therefore embedded as-is.
