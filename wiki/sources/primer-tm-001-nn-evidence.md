---
type: source
title: "Evidence: PRIMER-TM-001-NN (Per-oligo nearest-neighbour salt-corrected design Tm)"
tags: [validation, primer]
doc_path: docs/Evidence/PRIMER-TM-001-NN-Evidence.md
sources:
  - docs/Evidence/PRIMER-TM-001-NN-Evidence.md
source_commit: 01b5cc55b34e3f57b5833e3d47e818b8acd2c7c6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PRIMER-TM-001-NN

The validation-evidence artifact for test unit **PRIMER-TM-001** — the **per-oligo
(single-strand) nearest-neighbour salt-corrected design Tm** (opt-in). This is the core
ΔH°/ΔS° → Tm path plus the salt corrections; it is the same SantaLucia unified NN model
synthesized in [[primer-dimer-thermodynamics-tm]] (which also covers the `ntthal`
dimer/hairpin DPs). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for
how units are tracked. Siblings: [[primer-tm-001-dimer-evidence]] (intermolecular dimer),
[[primer-tm-001-hairpin-evidence]] (intramolecular hairpin), [[primer-tm-001-evidence]]
(the base weighted-penalty objective).

This artifact adds the **ΔH°/ΔS° → Tm** path and salt corrections on top of the NN unified
ΔG°37 table already validated for 3'-end stability under SEQ-THERMO-001; the legacy
Wallace / Marmur-Doty Tm and the Primer3 penalty objective are unchanged (documented
elsewhere).

## What this file records

- **Online sources:**
  - **SantaLucia & Hicks (2004)** (rank 1, Annu Rev Biophys 33:415–440) — read pages 419–423
    directly. **Table 1** unified NN ΔH°/ΔS° at 1 M NaCl (10 WC stacks + initiation +0.2/−5.7,
    terminal-A·T penalty +2.2/+6.9 per end, symmetry 0.0/−1.4 self-comp only). **Tm Eq. 3:**
    `Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15`, R = 1.9872 cal/(K·mol), x = 4
    non-self-complementary / x = 1 self-complementary. **Worked example (p.419):** ΔH°=−43.5,
    ΔS°=−122.5, 0.2 mM each strand → **Tm = 35.8 °C** (reproduced). **Sodium Eq. 5:**
    `ΔS°[Na⁺] = ΔS°[1 M] + 0.368·(N/2)·ln[Na⁺]`, N = total phosphates, ΔH° independent of
    [Na⁺]; valid 0.05–1.1 M Na⁺, duplexes ≤ 16 bp (6-bp N=10 → ΔG°37[0.115 M]=−4.12 kcal/mol
    reproduced exactly). **Tables 2 & 3** (`pdftotext -layout`) cross-check the mismatch and
    dangling-end parameters below.
  - **Owczarzy et al. (2004)** (rank 1, Biochemistry 43:3537–54; coefficients via Biopython
    method 6) — **monovalent quadratic 1/Tm correction:**
    `1/Tm[Na] = 1/Tm[1 M] + (4.29·f(GC) − 3.95)·1e−5·ln[Na⁺] + 9.40e−6·(ln[Na⁺])²`, f(GC)=GC
    fraction; Tm vs ln[Na⁺] is **non-linear** over 69 mM–1.02 M.
  - **Owczarzy et al. (2008)** (rank 1, Biochemistry 47:5336–53; coefficients via Biopython
    method 7) — **divalent Mg²⁺/dNTP 1/Tm correction:**
    `corr = a + b·ln[Mg] + f(GC)·(c + d·ln[Mg]) + (1/(2(N−1)))·(e + f·ln[Mg] + g·ln[Mg]²)`,
    base coeff (×1e−5) a=3.92,b=−0.911,c=6.26,d=1.42,e=−48.2,f=52.5,g=8.31; mixed regime
    (0.22 ≤ R < 6, R=√[Mg]/[Mon]) reparameterises a,d,g; free Mg²⁺ reduced by dNTP chelation
    via Ka = 3×10⁴.
  - **Biopython `Bio.SeqUtils.MeltingTemp`** (rank 3, raw source fetched) — independent
    cross-check of the tables and equation. `DNA_NN4` (SantaLucia 1998) matches Table 1
    verbatim; `Tm_NN` uses R=1.987, `Tm = 1000·ΔH/(ΔS + R·ln(k)) − 273.15` (selfcomp ⇒ k=C_T
    + `sym`; else k uses the (dnac1 − dnac2/2) form ≡ x=4 for equal strands). `salt_correction`
    method 5 = SantaLucia Eq. 5 entropy form `0.368·(len−1)·ln(mon)`; method 6 = Owczarzy 2004;
    method 7 = Owczarzy 2008. `DNA_IMM` (Allawi & SantaLucia 1997/1998, Peyret 1999) 52
    single-mismatch keys copied verbatim into `NnInternalMismatch` (WC placeholders + inosine
    excluded); `DNA_DE` (Bommarito 2000) 32 dangling-end keys into `NnDanglingEnd`.
  - **`Tm_NN` mismatch/dangling convention (read verbatim):** `c_seq = complement` (NOT
    reverse-complement); NN key `seq[i:i+2]+"/"+cseq[i:i+2]` looked up forward then
    character-reversed; left/right dangling ends strip the outer column; `ends = seq[0]+seq[-1]`
    (terminal-AT count on the un-dotted top strand); `sym` only for self-complementary. The C#
    `CalculateNearestNeighborThermodynamicsMismatch` mirrors this exactly.
  - **Bommarito et al. (2000)** (rank 1, NAR 28:1929–34) — primary dangling-end measurement,
    32 single-nt dangling ends by UV melting in 1 M NaCl, ΔG°37 +0.48 to −0.96 kcal/mol.
  - **OligoPool tutorial** (rank 4) — corroborates the Eq. 3 form only.
- **Datasets (documented oracles):**
  - *Published worked example:* ΔH°=−43.5, ΔS°=−122.5, 0.2 mM each → **Tm 35.8 °C**.
  - *Hand-derived NN sums (Table 1), C_T=0.5 µM @ 1 M Na:* `GCGCGC` self x=1 → ΔH=−50.4,
    ΔS=−134.7, Tm 35.0473059911; `ATGCATGC` x=4 → ΔH=−57.1, ΔS=−156.5, Tm 30.4338060665;
    `CGCGAATTCGCG` self x=1 → ΔH=−100.6, ΔS=−272.1, Tm 61.1452300219. @ 50 mM Na: `GCGCGC`
    Owczarzy-2004 28.1593085080 / SantaLucia-Eq.5 24.9976652723; `ATGCATGC` Owczarzy 18.1899960529.
  - *Mismatch/dangling (DNA_IMM/DNA_DE), no-salt C_T=0.5 µM:* `CGTGAC/GCGCTG` one internal T·G
    mismatch → ΔH=−35.5, ΔS=−101.5, Tm −6.4060879279; `AGCGCGC/.CGCGCG` one 5'-dangling A →
    ΔH=−51.9, ΔS=−136.4, Tm 35.8034921829; `GCGCGC/CGCGCG` perfect → Tm 35.0473059911 (equals
    perfect-match path). SantaLucia & Hicks Table 2 example `GGACTGACG/CCTGGCTGC` (two internal
    mismatches) ΔG°37 −8.32 (paper, rounded increments) vs −8.349 (Biopython DNA_NN+DNA_IMM,
    ΔH=−58.5/ΔS=−161.7) — within rounding; all 32 DNA_DE ΔH° values reproduce Table 3 exactly.
- **Corner cases / failure modes:**
  - Self-complementary duplex: x=1 (not 4), add symmetry ΔS°=−1.4, no terminal-A·T penalty on
    G·C ends; both-A·T-ends incurs +2.2/+6.9 twice.
  - Salt-correction validity: Eq. 5 holds 0.05–1.1 M Na⁺, duplexes ≤ 16 bp.
  - Non-ACGT base: no NN parameter → thermodynamics not computable → returns null / NaN.

## Deviations and assumptions

- **ASSUMPTION: default C_T = 0.5 µM** — common PCR primer working concentration, exposed as a
  caller-overridable parameter; affects only the default operating point, not correctness.
- **ASSUMPTION: Eq. 5 phosphate count N = 2·(length − 1)** — confirmed against the paper's
  6-bp → N=10 example (no 5'-terminal phosphates); −4.12 kcal/mol reproduces exactly.

No contradictions among sources; the Biopython tables are verified faithful transcriptions of
the SantaLucia/Allawi/Peyret/Bommarito primaries. Recommended coverage — **MUST:** Table 1
ΔH°/ΔS° sums (self- and non-self-comp) term-by-term; the 35.8 °C worked example; Owczarzy 2004
lowers Tm to the hand-derived 1/Tm value. **SHOULD:** Eq. 5 entropy correction; self-comp x=1
vs non-self-comp x=4; Mg²⁺ raises Tm; divalent-with-no-Mg falls back to monovalent. **COULD:**
monotonic Tm vs [Na⁺]; default-Tm-method-unchanged guard; invalid input → NaN.
