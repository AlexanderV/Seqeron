# Provenance & Licence — `maxent_score5.txt` (MaxEntScan score5ss tables)

**File:** `maxent_score5.txt` — the precomputed maximum-entropy probability table for the
Yeo & Burge (2004) MaxEntScan **score5ss** (5' splice-site / donor) model.

## What it is

A flat `7mer value` table, one record per line (tab-separated):

```
<7-nt rest sequence>   <probability>
```

Unlike score3 (nine overlapping sub-matrices), the score5 model is **single-matrix**: the
9-nt donor window (3 exon + 6 intron nt) has its invariant `GT` dinucleotide (window
positions 3..4, 0-based) scored separately by the consensus/background model and removed,
leaving a 7-nt "rest" sequence (`window[0:3] + window[5:9]`). The maximum-entropy probability
of that rest sequence is looked up directly by its string key.

| sub-sequence (positions of the 9-nt window) | length k | entries (4^k) |
|----------------------------------------------|----------|---------------|
| rest = window[0:3] + window[5:9]             | 7        | 16384         |

Total: 16 384 records (one matrix, 4^7 entries).

## Source (retrieved this session, 2026-06-25)

- **Port / data file:** `kepbod/maxentpy` — a Python port of MaxEntScan.
  Data file fetched verbatim from:
  `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/data/score5_matrix.txt`
  Reference factorisation (`score5`) fetched from:
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

Reproduced the documented MaxEntScan score5 worked examples exactly with this table + the
factorisation: `score5('cagGTAAGT')` → 10.858313 (→ 10.86, the canonical doc value);
`score5('gagGTAAGT')` → 11.078494 (→ 11.08); `score5('taaATAAGT')` → -0.116791 (→ -0.12).
Expected values from the maxentpy `score5` docstring.
