# Evidence Artifact: ASSEMBLY-CONSENSUS-001

**Test Unit ID:** ASSEMBLY-CONSENSUS-001
**Algorithm:** Consensus Computation (column-wise majority/threshold consensus from aligned reads)
**Date Collected:** 2026-06-13

---

## Online Sources

### Biopython `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus` (reference implementation)

**URL:** https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established bioinformatics library, Biopython v1.79)

**Retrieval:** Located via WebSearch query `Biopython SummaryInfo dumb_consensus gap_consensus threshold ambiguous source code`; the file's raw source was fetched at the `biopython-179` tag and the `dumb_consensus` method quoted verbatim.

**Key Extracted Points:**

1. **Signature / defaults:** `def dumb_consensus(self, threshold=0.7, ambiguous="X", require_multiple=False):` — default plurality `threshold = 0.7`, default ambiguous symbol `"X"`.
2. **Consensus length:** `con_len = self.alignment.get_alignment_length()` — the consensus is built over the full alignment length (the maximum column index), not the first record's length.
3. **Per-column tally (gap skipping):** for each column `n`, residues are counted into `atom_dict` only when `record.seq[n] != "-" and record.seq[n] != "."`; `num_atoms` is the number of non-gap residues contributing to that column. Sequences shorter than `con_len` are skipped at that column via `if n < len(record.seq)`.
4. **Max-count set (tie detection):** the code builds `max_atoms` — the list of residues sharing the maximum count: `if atom_dict[atom] > max_size: max_atoms=[atom]; max_size=...` / `elif atom_dict[atom] == max_size: max_atoms.append(atom)`.
5. **Decision rule (verbatim):**
   ```python
   if require_multiple and num_atoms == 1:
       consensus += ambiguous
   elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
       consensus += max_atoms[0]
   else:
       consensus += ambiguous
   ```
   A residue is emitted **only when** exactly one residue has the maximum count (`len(max_atoms) == 1`) **and** its frequency among non-gap residues meets the threshold (`max_size / num_atoms >= threshold`). Otherwise the ambiguous symbol is emitted.
6. **Tie → ambiguous:** when two or more residues tie for the maximum count (`len(max_atoms) > 1`), the first elif is False, so the `else` branch runs and the ambiguous symbol is emitted (NOT an arbitrary pick of one tied residue).
7. **All-gap column (num_atoms == 0):** `max_atoms` stays `[]`, so `len(max_atoms) == 1` is False; Python short-circuits the `and`, the division `max_size/num_atoms` is never evaluated, and the `else` branch emits the ambiguous symbol. No division by zero.

### EMBOSS `cons` — consensus from a multiple alignment (specification of the "plurality" concept)

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/cons.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation / tool documentation, EMBOSS suite)

**Retrieval:** WebSearch query `EMBOSS cons consensus sequence majority rule plurality alignment definition`; documentation describes the plurality cut-off.

**Key Extracted Points:**

1. **Purpose:** `cons` "calculates a consensus sequence from a multiple sequence alignment" — the consensus represents the most common residue at each column.
2. **Plurality cut-off:** the "plurality" qualifier "sets the cut-off for the number of positive matches (weighted) below which there is no consensus" — i.e. when support falls below the cut-off the column has no consensus residue (emitted as a non-committal symbol). This corroborates Biopython's threshold semantics: insufficient support → ambiguous/no-consensus.

### Wikipedia "Consensus sequence" (definition + cited primaries)

**URL:** https://en.wikipedia.org/wiki/Consensus_sequence
**Accessed:** 2026-06-13
**Authority rank:** 4 (encyclopedic; used for definition and to surface primary sources)

**Retrieval:** WebFetch of the article. Definition extracted; primaries surfaced (Schneider & Stephens 1990; Schneider 2002; Pierce 2002).

**Key Extracted Points:**

1. **Definition (verbatim):** a consensus sequence is "the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment."
2. **IUPAC ambiguity:** the article uses IUPAC degenerate notation (e.g. `N` = any base, `Y` = pyrimidine, `R` = purine) for positions without a single committed residue.

---

## Documented Corner Cases and Failure Modes

### From Biopython `dumb_consensus`

1. **Gap characters skipped:** `-` and `.` are not counted; a column's consensus is decided only among non-gap residues (point 3).
2. **Tie between residues:** two or more residues with equal maximum count → ambiguous symbol, not an arbitrary winner (point 6).
3. **Sub-threshold majority:** a single most-common residue whose frequency `< threshold` → ambiguous symbol (point 5).
4. **All-gap / empty column:** `num_atoms == 0` → ambiguous symbol, no division by zero (point 7).
5. **Ragged (variable-length) reads:** the consensus spans the longest read; columns beyond a shorter read contribute nothing from that read (points 2–3).

### From EMBOSS `cons`

1. **Below plurality cut-off → no consensus residue** at that column (point 2).

---

## Test Datasets

### Dataset: Sub-threshold majority (Biopython rule, threshold 0.7)

**Source:** Derived from Biopython `dumb_consensus` decision rule (threshold = 0.7).

| Parameter | Value |
|-----------|-------|
| Column residues | A, A, T (counts A=2, T=1; num_atoms=3) |
| max_size / num_atoms | 2/3 ≈ 0.667 |
| 0.667 ≥ 0.7 ? | No |
| Consensus symbol | ambiguous (`N` for DNA / `X` Biopython default) |

### Dataset: Threshold met

**Source:** Derived from Biopython rule (threshold = 0.7).

| Parameter | Value |
|-----------|-------|
| Column residues | A, A, A, T (counts A=3, T=1; num_atoms=4) |
| max_size / num_atoms | 3/4 = 0.75 |
| 0.75 ≥ 0.7 ? | Yes |
| Consensus symbol | `A` |

### Dataset: Tie → ambiguous

**Source:** Biopython rule (point 6).

| Parameter | Value |
|-----------|-------|
| Column residues | A, G (counts A=1, G=1; max_atoms = {A,G}) |
| len(max_atoms) | 2 |
| Consensus symbol | ambiguous (`N`) |

### Dataset: Gaps skipped + ragged length

**Source:** Biopython rule (points 2–3).

| Parameter | Value |
|-----------|-------|
| Reads | `A-GT`, `ACGT`, `ACG` (lengths 4,4,3) |
| Col 0 | A,A,A → A |
| Col 1 | (gap),C,C → C (gap skipped, 2/2 ≥ 0.7) |
| Col 2 | G,G,G → G |
| Col 3 | T,T,(absent) → T (2/2 ≥ 0.7) |
| Consensus | `ACGT` |

---

## Assumptions

1. **ASSUMPTION: Ambiguous symbol default = `N` (not Biopython's `X`).** Biopython's `dumb_consensus` default ambiguous symbol is `X` (point 1). For DNA/RNA assembly the IUPAC "any base" symbol is `N` (Wikipedia point 2; existing repository convention in `ComputeConsensus`). The repository default is therefore `N`; the symbol is exposed as a parameter so callers can request `X` (protein) or any other symbol. Default value choice (`N` vs `X`) is presentation-only — it does not change the decision rule — but is documented here because it changes the emitted character.
2. **ASSUMPTION: Default threshold = 0.5.** Biopython's documented default is `0.7` (point 1), but "plurality / majority" (EMBOSS plurality default is half the total weight) corresponds to a simple-majority cut-off of 0.5. To preserve a true majority-vote default while remaining configurable, the default threshold is set to `0.5`; the exact Biopython behavior is reproduced by passing `threshold: 0.7`. The threshold is a parameter, so all source-defined values are reachable. The decision rule itself (strict `>= threshold`, tie→ambiguous, gaps skipped) is fully source-backed and not assumed.

---

## Recommendations for Test Coverage

1. **MUST Test:** unanimous column returns that residue — Evidence: Biopython rule, frequency 1.0 ≥ threshold.
2. **MUST Test:** majority above threshold returns the majority residue (3/4 ≥ 0.7) — Evidence: Biopython threshold formula.
3. **MUST Test:** sub-threshold single majority returns ambiguous (2/3 < 0.7) — Evidence: Biopython point 5.
4. **MUST Test:** exact tie returns ambiguous, not an arbitrary residue — Evidence: Biopython point 6.
5. **MUST Test:** gap characters (`-`, `.`) are skipped from the tally — Evidence: Biopython point 3.
6. **MUST Test:** consensus length equals the longest read; ragged columns handled — Evidence: Biopython points 2–3.
7. **MUST Test:** all-gap / empty column returns ambiguous (no division by zero) — Evidence: Biopython point 7.
8. **MUST Test:** empty read list returns empty string — Evidence: trivial (no columns).
9. **MUST Test:** configurable threshold reproduces Biopython 0.7 behavior — Evidence: Biopython point 1.
10. **SHOULD Test:** null input throws — Rationale: standard repository contract (sibling methods use `ArgumentNullException.ThrowIfNull`).
11. **COULD Test:** custom ambiguous symbol (`X`) is emitted — Rationale: parameter exposure for protein alignments (Biopython default).

---

## References

1. Cock PJA, Antao T, Chang JT, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 — `dumb_consensus` source (v1.79): https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py
2. Rice P, Longden I, Bleasby A (2000). EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 — `cons` documentation: https://emboss.sourceforge.net/apps/cvs/emboss/apps/cons.html
3. Schneider TD, Stephens RM (1990). Sequence logos: a new way to display consensus sequences. *Nucleic Acids Research* 18(20):6097–6100. https://doi.org/10.1093/nar/18.20.6097 (cited via Wikipedia "Consensus sequence": https://en.wikipedia.org/wiki/Consensus_sequence).

---

## Change History

- **2026-06-13**: Initial documentation.
