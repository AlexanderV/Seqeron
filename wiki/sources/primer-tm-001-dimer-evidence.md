---
type: source
title: "Evidence: PRIMER-TM-001-DIMER (Self-/hetero-dimer Tm via ntthal thermodynamic alignment)"
tags: [validation, primer]
doc_path: docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
sources:
  - docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
source_commit: 6c22f2b1ee758fb6d2b4c748c139d0193a4e313a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PRIMER-TM-001-DIMER

The validation-evidence artifact for test unit **PRIMER-TM-001** — **self-dimer /
hetero-dimer melting temperature** for PCR primers, computed by a **thermodynamic
alignment** over the SantaLucia unified nearest-neighbour (NN) DNA model. It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the NN thermodynamics, bimolecular Tm equation, salt correction, symmetry
concentration factor, and the full `ntthal` dimer DP are synthesized in
[[primer-dimer-thermodynamics-tm]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **SantaLucia & Hicks (2004)**, *A unified view of ... oligonucleotide DNA NN
    thermodynamics* (rank 1, Annu Rev Biophys 33:415–440) — the primary unified parameter
    set: the 10 Watson-Crick NN ΔH°/ΔS° at 1 M NaCl (repo `NnUnifiedParams`), duplex
    initiation (ΔH°=+0.2, ΔS°=−5.7), per-end terminal A·T penalty (ΔH°=+2.2, ΔS°=+6.9),
    self-complementary symmetry ΔS°=−1.4. **Bimolecular Tm (Eq. 3):**
    `Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15`, R = 1.9872, x = 4 (non-self-complementary)
    or x = 1 (self-complementary). **Entropy salt correction (Eq. 5):**
    `ΔS°[Na⁺] = ΔS°[1 M] + 0.368·(N/2)·ln[Na⁺]`.
  - **Untergasser et al. (2012)**, *Primer3 — new capabilities* (rank 1, Nucleic Acids Res
    40(15):e115) — describes the `ntthal` **thermodynamic-alignment** engine that scores
    oligo–oligo / oligo–template secondary structures (dimers, hairpins) over the SantaLucia
    NN model, exposed as `calcHomodimer` / `calcHeterodimer`.
  - **Primer3 `thal.c`** (ntthal reference implementation, primer3-py vendored libprimer3;
    rank 3) — the verbatim constants and DP recurrence: `dplx_init_H = 200` cal,
    `dplx_init_S = −5.7`; `AT_H = 2200`, `AT_S = 6.9` per A·T-closed end;
    `saltCorrection = 0.368·ln((mv + 120·√max(0, dv − dntp))/1000)` × N stacks;
    concentration divisor `RC = R·ln(dna_conc/1e9)` when both strands symmetric else
    `/4e9`; symmetry = even-length reverse-complement palindrome; defaults `dna_conc = 50` nM,
    `temp = 310.15` K.
  - **primer3-py 2.3.0** `calc_homodimer` / `calc_heterodimer` (rank 3) — the numeric oracle
    (mv=50, dv=0, dntp=0, dna_conc=50 nM) for both contiguous-WC and non-contiguous cases.
- **Datasets (documented oracles), mv=50 mM, dna_conc=50 nM:**
  - *Contiguous-WC parity:* `GCGCGCGC` self (x=1) → ΔH=−70800, ΔS=−192.6170, Tm=40.0906 °C;
    `TGCATGCATG`/`CATGCATGCA` non-palindromic (x=4) → ΔH=−74100, ΔS=−211.8219, Tm=25.6596 °C;
    plus `ACGTACGTACGT`, `ATCGATCGATCG`/`CGATCGATCGAT`, `CGATCGATCG` (palindrome → x=1),
    `GCATGC`, `GGGGCCCC`. GCGCGCGC and the x=4 case are also **hand-derived** from Table 1
    to machine precision, independent of primer3.
  - *Non-contiguous ntthal DP:* `GCGCATGCGC` (2×2 internal loop, Tm=43.1572), 3×3 loop,
    1-base bulge, terminal overhang — the C# port reproduces every value to machine precision
    (ΔTm=0), and still reproduces all contiguous cases (regression).
- **Corner cases / failure modes:**
  - **No stable dimer:** poly-A `AAAAAAAA` self-dimer → `structure_found = False` → null / NaN.
  - **Symmetry divisor:** x = 1 only when *both* aligned oligos are reverse-complement
    palindromes; a self-dimer of a non-palindromic oligo uses x = 4.
  - **ntthal terminal-stack / overhang extension** adds `tstack2` terms beyond plain
    contiguous WC (a documented small terminal deviation for the retained contiguous scorer).
  - Invalid input (null, < 2 bases, non-ACGT) → null / NaN.

## Deviations and assumptions

**Deviations: none remaining (RESOLVED 2026-06-25)** — the dimer alignment now ships the
**complete** Primer3 `ntthal` oligo–oligo DP (mode ANY): internal single mismatch, internal
loop, single/multi-base bulge, and `tstack2` terminal overhang/dangling end, a verbatim port
of `thal.c` (`fillMatrix`, `LSH`, `RSH`, `maxTM`, `calc_bulge_internal`, `traceback`,
`calcDimer`) with the parameter tables embedded from `primer3_config/*.dh,*.ds`. Method
`PrimerDesigner.CalculateDimerThermodynamicsNtthal`; `CalculateDimerMeltingTemperature` /
`CalculateSelfDimerMeltingTemperature` delegate to it (the earlier contiguous `FindMostStableDimer`
scorer is retained for its existing contiguous-WC results). The only unmodelled capability is the
optional caller-supplied tri/tetraloop & terminal-mismatch **hairpin bonus tables** (a
hairpin/monomer feature, not a dimer one). No contradictions among sources.

Recommended coverage (MUST): GCGCGCGC self-dimer ΔH°/ΔS°/Tm exactly; the x=4 non-palindromic
hetero-dimer; the primer3-parity set; poly-A / non-complementary → no dimer (null/NaN); invalid
input → null/NaN. SHOULD: self-dimer convenience method equals the two-argument dimer with the
same sequence; monotonic salt dependence (lower [Na⁺] → lower Tm). COULD: most-stable selection
prefers the higher-Tm alignment.
