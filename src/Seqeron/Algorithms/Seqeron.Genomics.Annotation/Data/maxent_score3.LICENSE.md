# Provenance & Licence — `maxent_score3.txt` (MaxEntScan score3ss tables)

**File:** `maxent_score3.txt` — the precomputed maximum-entropy probability tables for the
Yeo & Burge (2004) MaxEntScan **score3ss** (3' splice-site / acceptor) model.

## What it is

A flat `index hash value` table, one record per line (tab-separated):

```
<matrixIndex 0..8>   <hashseq value>   <probability>
```

The nine matrices and their sizes (4^k entries each) reproduce the
MaxEntScan factorisation of the 21 non-consensus positions of the 23-nt acceptor window:

| index | sub-sequence (positions of the 21-nt rest) | length k | entries (4^k) |
|-------|--------------------------------------------|----------|---------------|
| 0 | rest[0:7]   | 7 | 16384 |
| 1 | rest[7:14]  | 7 | 16384 |
| 2 | rest[14:21] | 7 | 16384 |
| 3 | rest[4:11]  | 7 | 16384 |
| 4 | rest[11:18] | 7 | 16384 |
| 5 | rest[4:7]   | 3 | 64    |
| 6 | rest[7:11]  | 4 | 256   |
| 7 | rest[11:14] | 3 | 64    |
| 8 | rest[14:18] | 4 | 256   |

Total: 82 560 records.

## Source (retrieved this session, 2026-06-24)

- **Port / data file:** `kepbod/maxentpy` — a Python port of MaxEntScan.
  Data file fetched verbatim from:
  `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/data/score3_matrix.txt`
  Reference factorisation (`score3`) fetched from:
  `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/maxent.py`
- **Original algorithm:** Yeo G, Burge CB (2004). "Maximum entropy modeling of short sequence
  motifs with applications to RNA splicing signals." *J Comput Biol* 11(2–3):377–394.
  DOI 10.1089/1066527041410418.

## Licence (read carefully)

- **The `maxentpy` port — including this data file — is distributed under the MIT License**
  (`https://github.com/kepbod/maxentpy/blob/master/LICENSE`). Full MIT text:

  > MIT License
  > Permission is hereby granted, free of charge, to any person obtaining a copy of this
  > software and associated documentation files (the "Software"), to deal in the Software
  > without restriction, including without limitation the rights to use, copy, modify, merge,
  > publish, distribute, sublicense, and/or sell copies of the Software … The above copyright
  > notice and this permission notice shall be included in all copies or substantial portions
  > of the Software. THE SOFTWARE IS PROVIDED "AS IS" …

  MIT permits redistribution, so this table is bundled here under MIT.

- **FLAG — original MaxEntScan terms are academic.** maxentpy's own README states:
  *"The original algorithm and perl scripts are under license described in
  http://genes.mit.edu/burgelab/maxent/download/READTHIS. The python version of maxent is
  under the MIT License."* The Burge-lab `READTHIS` (academic terms) governs the *original*
  Perl scripts/models; the redistributable artifact embedded here is the **MIT-licensed**
  `maxentpy` table. A maintainer wanting belt-and-suspenders clearance for commercial
  redistribution should review the upstream Burge-lab terms directly.

## Cross-check (this session)

Reproduced the documented MaxEntScan worked examples exactly with this table + the
factorisation: `score3('ttccaaacgaacttttgtAGgga')` → 2.886773 (→ 2.89, the canonical doc value);
`score3('tgtctttttctgtgtggcAGtgg')` → 8.190965 (→ 8.19); `score3('ttctctcttcagacttatAGcaa')`
→ -0.080278 (→ -0.08). Expected values from the maxentpy `score3` docstring.
