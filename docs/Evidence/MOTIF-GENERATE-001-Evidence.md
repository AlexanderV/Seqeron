# Evidence Artifact: MOTIF-GENERATE-001

**Test Unit ID:** MOTIF-GENERATE-001
**Algorithm:** IUPAC-Degenerate Consensus Generation (`MotifFinder.GenerateConsensus`)
**Date Collected:** 2026-06-14

---

## Online Sources

### Cornish-Bowden / NC-IUB (1984) — Nomenclature for incompletely specified bases in nucleic acid sequences

**URL:** https://academic.oup.com/nar/article/13/9/3021/2381659 (DOI 10.1093/nar/13.9.3021)
**Accessed:** 2026-06-14
**Authority rank:** 2 (official IUPAC/IUB nomenclature standard)
**Retrieved how:** WebSearch "Cornish-Bowden 1985 Nomenclature for incompletely specified bases…" located the NAR article; WebFetch of the Oxford Academic article page returned the symbol table.

**Key Extracted Points:**

1. **Single-letter degenerate symbols (verbatim mapping from the fetched table):** R = A or G (purine); Y = C or T/U (pyrimidine); M = A or C (amino); K = G or T/U (keto); S = G or C (strong); W = A or T/U (weak); B = C, G, or T/U (not-A); D = A, G, or T/U (not-C); H = A, C, or T/U (not-G); V = A, C, or G (not-T/U); N = any base (A, C, G, T/U).
2. **Meaning rule:** "each symbol stands for the specific set of bases listed" — a degenerate symbol denotes exactly the set of standard bases it abbreviates. The set→symbol mapping is bijective over the 11 non-trivial subsets.

### UCSC Genome Browser — IUPAC Nucleotide Code Table

**URL:** https://genome.ucsc.edu/goldenPath/help/iupac.html
**Accessed:** 2026-06-14
**Authority rank:** 5 (well-maintained genomics database documentation; corroborates the standard)
**Retrieved how:** WebSearch "IUPAC nucleotide ambiguity codes…"; WebFetch of the UCSC page returned the full table.

**Key Extracted Points:**

1. **Full table (verbatim):** G=G, A=A, T=T, C=C; R=G or A; Y=T or C; M=A or C; K=G or T; S=G or C; W=A or T; H=A or C or T (not-G); B=G or T or C (not-A); V=G or C or A (not-T); D=G or A or T (not-C); N=G or A or T or C (any).
2. **Use context:** single-character codes represent multiple observed alleles at a single position — i.e., the same set-of-bases semantics used when summarising an aligned column.

### Wikipedia — Nucleic acid notation (citing NC-IUB 1984)

**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_notation
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the NC-IUB 1984 primary)
**Retrieved how:** WebFetch of the article; Table 1 returned with its primary reference.

**Key Extracted Points:**

1. **Table 1 mapping** matches Cornish-Bowden/UCSC exactly (W,S,M,K,R,Y two-base; B,D,H,V three-base; N four-base).
2. **Primary reference cited:** "Nomenclature Committee of the International Union of Biochemistry (NC-IUB) (1984), Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences, Nucleic Acids Research" — confirming the table is the IUPAC/NC-IUB standard, not Wikipedia-original.

### Bioconductor DECIPHER — `ConsensusSequence` (reference implementation of threshold-based degenerate consensus)

**URL:** https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (established library reference implementation)
**Retrieved how:** WebSearch "degenerate consensus sequence IUPAC threshold…"; WebFetch returned the man page text.

**Key Extracted Points:**

1. **Threshold-consensus mechanism (verbatim):** "ConsensusSequence removes the least frequent characters at each position, so long as they represent less than `threshold` fraction of the sequences in total"; remaining characters are then "represented using IUPAC degeneracy codes."
2. **Tie / equal-abundance rule (verbatim):** "Degeneracy codes are always used in cases where multiple characters are equally abundant."
3. **Establishes the family:** a degenerate consensus is parameterised by a frequency threshold; bases above the threshold at a column are combined into the IUPAC symbol for that base set. (DECIPHER's default threshold is 0.05; the threshold value is a tunable parameter, not a fixed universal constant.)

---

## Documented Corner Cases and Failure Modes

### From DECIPHER `ConsensusSequence`

1. **Frequency threshold governs inclusion:** a minority base whose frequency is below the threshold is dropped and not encoded in the ambiguity code; the threshold value is a design parameter chosen per tool.
2. **Equal abundance → degeneracy code:** when ≥2 bases pass the inclusion rule, the column emits the IUPAC symbol for that set rather than picking one base arbitrarily.

### From NC-IUB 1984 / UCSC

1. **N is the four-base symbol:** N denotes the full set {A,C,G,T}; any single missing base yields a three-base not-X symbol (B/D/H/V) instead of N.

---

## Test Datasets

### Dataset: IUPAC set→symbol mapping (NC-IUB 1984 / UCSC)

**Source:** Cornish-Bowden NAR 13(9):3021 (DOI 10.1093/nar/13.9.3021); UCSC IUPAC table.

| Base set present at column | IUPAC symbol |
|----------------------------|--------------|
| {A} | A |
| {C} | C |
| {G} | G |
| {T} | T |
| {A,G} | R |
| {C,T} | Y |
| {C,G} | S |
| {A,T} | W |
| {G,T} | K |
| {A,C} | M |
| {C,G,T} | B |
| {A,G,T} | D |
| {A,C,T} | H |
| {A,C,G} | V |
| {A,C,G,T} | N |

### Dataset: Threshold-inclusion worked examples (this implementation, threshold = n × 0.25, base included iff count > threshold)

**Source:** derivation from the DECIPHER threshold-consensus family applied to this implementation's documented 25 % design constant; the resulting *symbol* per base set is dictated by the NC-IUB table above.

| Input column (n seqs) | n | threshold = n×0.25 | bases with count > threshold | symbol |
|-----------------------|---|--------------------|------------------------------|--------|
| A,G | 2 | 0.5 | {A,G} (each 1>0.5) | R |
| C,T | 2 | 0.5 | {C,T} | Y |
| C,G | 2 | 0.5 | {C,G} | S |
| A,T | 2 | 0.5 | {A,T} | W |
| G,T | 2 | 0.5 | {G,T} | K |
| A,C | 2 | 0.5 | {A,C} | M |
| C,G,T | 3 | 0.75 | {C,G,T} (each 1>0.75) | B |
| A,G,T | 3 | 0.75 | {A,G,T} | D |
| A,C,T | 3 | 0.75 | {A,C,T} | H |
| A,C,G | 3 | 0.75 | {A,C,G} | V |
| A,A,G,G,C | 5 | 1.25 | {A(2),G(2)}; C(1) dropped | R |
| A,A,A,G | 4 | 1.0 | {A(3)}; G(1) at ≤threshold dropped | A |
| A,C,G,T | 4 | 1.0 | none (each 1, not >1.0) → fallback to most-frequent | A (tie→alphabetical) |

---

## Assumptions

1. **ASSUMPTION: 25 % inclusion threshold is a documented design constant.** This implementation includes a base in a column's IUPAC code iff its count is strictly greater than 25 % of the number of sequences (`count > total × 0.25`). The *threshold-consensus family* and the "include bases above a frequency threshold, encode the set with an IUPAC code" rule are authoritative (DECIPHER); the specific 25 % cut and the strict `>` boundary are this implementation's design choice (DECIPHER's own default is 0.05 and tools vary). It is correctness-affecting but documented and named (`threshold = total * 0.25`), not invented-untraceable. Tests pin the boundary behaviour explicitly and otherwise use inputs where the inclusion decision is unambiguous so the verified *symbol* is dictated solely by the authoritative NC-IUB table.
2. **ASSUMPTION: Fallback when no base passes the threshold.** When no base exceeds the threshold (e.g. four equally-frequent bases at exactly 25 %), the implementation falls back to the single most-frequent base (ties broken by dictionary/alphabetical order). No authoritative spec defines this corner; it is an implementation contract, verified as a documented edge case.
3. **ASSUMPTION: Column length taken from the first sequence; case-insensitive over {A,C,G,T}; non-ACGT characters at a position are ignored in the counts.** Inputs are upper-cased; only A/C/G/T are counted per the four-base alphabet.

---

## Recommendations for Test Coverage

1. **MUST Test:** every two-base set maps to the correct IUPAC code (R,Y,S,W,K,M) and every three-base set to (B,D,H,V) — Evidence: NC-IUB 1984 / UCSC table.
2. **MUST Test:** unanimous columns reproduce the input base (A/C/G/T) — Evidence: NC-IUB (singleton set → standard base).
3. **MUST Test:** strict-`>` 25 % boundary — a base at exactly 25 % is excluded; a base above 25 % is included — Evidence: implementation design constant (documented).
4. **MUST Test:** fallback to most-frequent base when no base passes the threshold (four bases each at 25 %) — Evidence: implementation contract.
5. **SHOULD Test:** minority base below threshold dropped (A,A,G,G,C → R, not a three-base code) — Rationale: confirms threshold filtering precedes IUPAC encoding.
6. **SHOULD Test:** case-insensitivity (lowercase input) and empty-collection → "" — Rationale: documented input normalisation / guard.
7. **COULD Test:** null input throws `ArgumentNullException` — Rationale: documented guard.

---

## References

1. Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. Nucleic Acids Research 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021
2. UCSC Genome Browser. IUPAC ambiguity codes. https://genome.ucsc.edu/goldenPath/help/iupac.html (accessed 2026-06-14)
3. Wikipedia. Nucleic acid notation (Table 1, citing NC-IUB 1984). https://en.wikipedia.org/wiki/Nucleic_acid_notation (accessed 2026-06-14)
4. Wright E.S. DECIPHER `ConsensusSequence` (Bioconductor). https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
