# Validation Report: SPLICE-MAXENT5-001 — MaxEntScan score5 (5' Donor)

- **Validated:** 2026-06-25   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.ScoreDonorMaxEnt(string window)`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
`ScoreDonorMaxEnt` (MaxEntScan score5ss)

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs` (lines ~1229–1379)
- **Embedded table:** `src/.../Seqeron.Genomics.Annotation/Data/maxent_score5.txt` (16384 records, 4^7)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` (ME1–ME9)

## Authoritative sources opened this session
- **Yeo G, Burge CB (2004).** "Maximum entropy modeling of short sequence motifs with
  applications to RNA splicing signals." *J Comput Biol* 11(2–3):377–394. DOI 10.1089/1066527041410418.
  — defines the score5ss maximum-entropy 5' donor model: 9-nt window (3 exon + 6 intron),
  conserved GT, score = log2( P_maxent / P_background ).
- **maxentpy (`kepbod/maxentpy`, MIT)** `maxent.py::score5` + `data/score5_matrix.txt`, fetched this
  session from raw.githubusercontent.com. This is the canonical reference port of the Burge-lab
  `score5.pl`.

## Stage A — Description

**Model (Yeo & Burge 2004 / maxentpy `score5`).** The 9-nt donor window is
`(exon)XXX|GTXXXX(intron)` with the conserved GT at 0-based positions 3–4. The score factorises as:

```
score5(fa) = log2[ (cons1[G]·cons2[T]) / (bgd[G]·bgd[T])  ×  matrix[rest] ]
rest = fa[0:3] + fa[5:9]     (7-nt; the GT removed)
```

Single-matrix model (unlike score3's nine overlapping sub-matrices): the full maximum-entropy
probability of the 7-mer `rest` is stored directly, keyed by string. Constants confirmed verbatim
against maxentpy:

| const | A | C | G | T |
|-------|------|------|------|------|
| `bgd_5` | 0.27 | 0.23 | 0.23 | 0.27 |
| `cons1_5` (G pos) | 0.004 | 0.0032 | 0.9896 | 0.0032 |
| `cons2_5` (T pos) | 0.0034 | 0.0039 | 0.0042 | 0.9884 |

**Edge-case semantics (sourced).** maxentpy `score5` requires `len(fa)==9` (else `sys.exit`);
upper-cases input (case-insensitive); T/U equivalent for the rest key and GT consensus. A
non-canonical donor (e.g. `ATAAGT` at pos 3–4 → key `AT`) is *not* rejected — it scores via the
consensus table (very low, often negative), matching the published `taaATAAGT → -0.12` example.

**Independent cross-check (numbers).** Reimplemented `score5` in Python over the *embedded*
`maxent_score5.txt` (no `maxentpy` install needed) and reproduced the published docstring values
exactly:

| window | Python (this session) | round(,2) | maxentpy doc |
|--------|-----------------------|-----------|--------------|
| `cagGTAAGT` | 10.858313 | 10.86 | **10.86** |
| `gagGTAAGT` | 11.078494 | 11.08 | **11.08** |
| `taaATAAGT` | -0.116791 | -0.12 | **-0.12** |

Additional reproductions: `CAGGTAAGT` → 10.858313 (T/U-form identity); `ATGGTAAGG` → 9.331905;
`TTTGTTTTT` → -22.575296. Table loaded = 16384 entries.

**Stage A findings:** none. Description, formula, constants, log base, window layout, and
edge-case semantics match Yeo & Burge (2004) / maxentpy `score5` exactly.

## Stage B — Implementation

**Code path reviewed:** `ScoreDonorMaxEnt` (lines 1347–1377), the GT consensus dictionaries
(1259–1268), `MaxEntBackground` (reused from score3, 1031–1034), `NormalizeNucleotide` (U→T,
1219–1225), table loader (1293–1318).

**Formula realised correctly?** Yes — line-for-line equivalent to maxentpy:
- `consScore = cons1[G]·cons2[T] / (bgd[G]·bgd[T])` (1362–1364) ✓
- `rest = window[0:3] + window[5:9]`, 7 chars, U→T normalised (1368–1372) ✓
- `return Math.Log2(consScore · restScore)` (1376) ✓ — `Math.Log2` = log base 2 ✓
- Constants identical to maxentpy `bgd_5`/`cons1_5`/`cons2_5`; background dict shared with score3
  (same {0.27,0.23,0.23,0.27}, correct for score5). ✓

**Edge cases in code:** null → `ArgumentNullException` (1349); length≠9 → `ArgumentException`
(1350–1353); non-A/C/G/T(/U) → `ArgumentException` via `NormalizeNucleotide` (1223); U≡T and
case-insensitive via `ToUpperInvariant` + U→T. Non-GT donors are scored, not rejected (matches
source).

**Cross-verification table recomputed vs code:** the C# fixture asserts the exact
reference-reproduction values, which match the Python oracle bit-for-bit:

| window | oracle (Python/maxentpy) | C# assertion (ME1–ME4) |
|--------|--------------------------|-------------------------|
| `cagGTAAGT` | 10.858313 (→10.86) | `== 10.858313 ±1e-6`, `round==10.86` |
| `gagGTAAGT` | 11.078494 (→11.08) | `== 11.078494 ±1e-6`, `round==11.08` |
| `taaATAAGT` | -0.116791 (→-0.12) | `== -0.116791 ±1e-6`, `round==-0.12` |

**Variant/delegate consistency:** single public method; no `*Fast`/delegate variants. Shares the
score3 background dict and `NormalizeNucleotide` — both correct for score5.

**Test quality audit:** ME1–ME9 trace expected values to maxentpy/Yeo&Burge (not code echoes):
full-precision + 2-dp on all three published examples, strong>weak ordering (ME5), DNA/RNA
identity (ME6), case-insensitivity (ME7), null (ME8), wrong length 8 & 10 + non-ACGT 'N' (ME9).
Covers every Stage-A edge case named in the protocol (window≠9, no GT, non-ACGT, lowercase,
`cagGTAAGT→10.86`). No green-washing, no tautologies. Suite is genuine.

**Stage B findings:** none.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.** No defect found; implementation reproduces the
  Yeo & Burge (2004) / maxentpy `score5` oracle exactly (10.858313 / 11.078494 / -0.116791).
- Full unfiltered `dotnet test Seqeron.sln -c Debug` (net10.0): Failed: 0
  (Seqeron.Genomics.Tests 18762 passed), 2026-06-25.
- No follow-ups; no entry needed in FINDINGS_REGISTER.
