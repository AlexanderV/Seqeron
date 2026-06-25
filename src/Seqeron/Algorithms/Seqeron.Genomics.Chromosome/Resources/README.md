# Bundled alpha-satellite reference consensus monomers (CHROM-CENT-001)

`AlphaSatelliteReference.fasta` holds **three CC0 (public-domain) human alpha-satellite
consensus monomers** from the Dfam database. They are the reference set used by
`ChromosomeAnalyzer.AssignSuprachromosomalFamily(...)` to (a) confirm a query monomer is
alpha-satellite and (b) type each monomer as **A-type** or **B-type** (the sequence-defined
monomer dichotomy of human alpha satellite).

| Record | Dfam acc | LENG | Box type | Title (Dfam) |
|--------|----------|------|----------|--------------|
| `ALR`  | DF000000029 | 171 | A | Human alpha satellite DNA |
| `ALRa` | DF000000014 | 172 | A | Human alpha satellite DNA variant a |
| `ALRb` | DF000000015 | 169 | B | Human alpha satellite DNA variant b |

The **A/B box** typing of these three consensus monomers is sequence-derived and verifiable:
`ALRb` carries the 17-bp CENP-B box `YTTCGTTGGAARCGGGA` (Masumoto et al. 1989, *J Cell Biol*
109:1963) at consensus position 126 → **B-type**; `ALR` and `ALRa` lack it → **A-type**. Per
McNulty & Sullivan (2018, *Chromosome Res*; PMC6121732): "A-type monomers include J1, D2, W4,
W5, M1, and R2 monomers, while B-type consist of J2, D1, W1–W3, and R1 monomers. B-type
monomers contain CENP-B boxes; A-type contain pJα binding sites."

## Provenance

- Retrieved 2026-06-25 via the Dfam REST API, e.g.
  `https://www.dfam.org/api/families/DF000000029` (JSON `consensus_sequence` field), copied
  verbatim into FASTA, no base edits. The other two records are `DF000000014` (ALRa) and
  `DF000000015` (ALRb).
- Dfam classification of all three: `root; Tandem_Repeat; Satellite; Centromeric`.

## Licence — CC0 (public domain)

Dfam is released under the Creative Commons Zero ("CC0") public-domain dedication
(Storer et al. 2021, *Mobile DNA* 12:2, "The Dfam community resource …"; dfam.org). CC0 places
the data in the public domain, so it is freely redistributable; these consensus strings are
embedded as-is.

## Scope / residual

These three consensus monomers resolve **alpha-satellite-vs-not** and **A-type-vs-B-type
monomer**. They do NOT carry the **suprachromosomal-family-resolved** consensus monomer
sequences (J1, J2, D1, D2, W1–W5, M1, R1, R2). Those SF-resolved consensus sequences are
published only in (i) supplementary matrices that are not machine-retrievable as plain FASTA
and (ii) third-party HMM repositories (e.g. enigene/HumAS-HMMER, logsdon-lab) that ship **no
LICENSE file** and are therefore not redistributable here. Consequently
`AssignSuprachromosomalFamily` assigns the family from the two reproducible SF-determining
signals defined in the literature — the **HOR periodicity** (SF1/SF2 dimeric, SF3 pentameric,
SF4 monomeric, SF5 irregular) and the **A/B-type monomer composition** — which uniquely
resolves SF3, SF4 and SF5 and narrows dimeric arrays to the {SF1, SF2} pair. Separating SF1
from SF2 (both dimeric, same A→B alternation) would need the SF-resolved consensus monomer
library, which is not CC0 / not redistributable — see the algorithm-doc residual and
`docs/Validation/LIMITATIONS.md`.
