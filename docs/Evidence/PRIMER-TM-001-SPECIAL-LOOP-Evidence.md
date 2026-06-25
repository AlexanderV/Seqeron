# Evidence Artifact: PRIMER-TM-001-SPECIAL-LOOP

**Test Unit ID:** PRIMER-TM-001 (special tri/tetraloop hairpin bonus tables — bundled & applied)
**Algorithm:** Full Primer3 `ntthal` intramolecular-hairpin DP with bundled sequence-specific
special triloop / tetraloop stability bonuses
**Date Collected:** 2026-06-25

---

## Online Sources

### primer3 `primer3_config/triloop.dh` & `triloop.ds` (libprimer3, vendored in primer3-py)

**URL:** https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_config/triloop.dh
and …/triloop.ds
**Accessed:** 2026-06-25 (both files fetched verbatim with `curl`, saved locally)
**Authority rank:** 3 (reference implementation — the authoritative parameter set `ntthal` loads)

**Key Extracted Points:**

1. **Keying:** each line is `<loop> <value>`, where `<loop>` is the **5-character** loop string
   *including the closing base pair* — closing-5′ base + 3 loop nt + closing-3′ base
   (`readTLoop` copies 5 chars for a triloop; thal.c line 1161). There are exactly **16** triloop
   entries.
2. **ΔH (triloop.dh, cal/mol), verbatim:** AGAAT −1500, AGCAT −1500, AGGAT −1500, AGTAT −1500,
   CGAAG −2000, CGCAG −2000, CGGAG −2000, CGTAG −2000, GGAAC −2000, GGCAC −2000, GGGAC −2000,
   GGTAC −2000, TGAAA −1500, TGCAA −1500, TGGAA −1500, TGTAA −1500.
3. **ΔS (triloop.ds, cal/(K·mol)), verbatim:** all 16 entries = **0**.

### primer3 `primer3_config/tetraloop.dh` & `tetraloop.ds`

**URL:** https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_config/tetraloop.dh
and …/tetraloop.ds
**Accessed:** 2026-06-25 (both files fetched verbatim with `curl`, saved locally)
**Authority rank:** 3

**Key Extracted Points:**

1. **Keying:** each line is `<loop> <value>`, where `<loop>` is the **6-character** loop string
   *including the closing base pair* — closing-5′ base + 4 loop nt + closing-3′ base
   (`readTLoop` copies 6 chars for a tetraloop; thal.c line 1164). There are **76** tetraloop entries.
2. **ΔH (tetraloop.dh, cal/mol) / ΔS (tetraloop.ds, cal/(K·mol)):** transcribed verbatim into
   `NtthalHairpin.TetraloopTable` (e.g. CGAAAG ΔH=−1100 ΔS=0; GGGGAC ΔH=−1100 ΔS=0;
   AAAAAT ΔH=500 ΔS=−650; ACTTGT ΔH=0 ΔS=4190; …). The full 76-row table is embedded in the source
   with this provenance header.

### primer3 `thal.c` — `calc_hairpin` (how the bonus is applied)

**URL:** https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c
**Accessed:** 2026-06-25 (fetched verbatim; `calc_hairpin` lines 2067-2146 read)
**Authority rank:** 3

**Key Extracted Points:**

1. **Application (verbatim, lines 2090-2128):** the hairpin loop free energy starts from the
   size-keyed `hairpinLoop{Enthalpies,Entropies}[loopSize-1]`. For `loopSize > 3` a terminal-mismatch
   (`tstack2`) increment is added; for `loopSize == 3` the closing-A·T penalty
   (`atPenaltyH/atPenaltyS`) is added. Then **for `loopSize == 3` the triloop bonus** (bsearch over
   `triloopEnthalpies`/`triloopEntropies`, key `numSeq1 + i`) and **for `loopSize == 4` the tetraloop
   bonus** (bsearch over `tetraloop{Enthalpies,Entropies}`, key `numSeq1 + i`) are **added** to ΔH and ΔS.
2. **Bonus convention:** ΔH in cal/mol, ΔS in cal/(K·mol); both are *added* to the loop ΔH°/ΔS°
   (negative ΔH stabilises). The bsearch key starts at index `i` (the closing-pair 5′ base), so the
   key is the full loop string including the closing pair, exactly as the table files are written.
3. **Loop-length tables (`loops.dh`/`loops.ds`, hairpin column):** ΔH° = 0 for all sizes; ΔS° by size
   (size 3=−11.28, 4=−11.28, 5=−10.64, … size 30=−20.31). Used verbatim in `NtthalHairpin.HairpinLoopS`.
4. **Unimolecular init / Tm:** for a monomer (`type==4`) `dplx_init_H=0`, `dplx_init_S=−1e-11`,
   `RC=0` (lines 583-585). `calcHairpin` (lines 3229-3266): `Tm = mh/(ms + (N/2 − 1)·saltCorrection) − 273.15`,
   where `saltCorrection = 0.368·ln(mv/1000)` (dv=dntp=0), `mh=HEND5(len1)`, `ms=SEND5(len1)`, and N is
   the count of paired positions over `bp[0 .. len1−2]`. No strand-concentration term (intramolecular).

### SantaLucia & Hicks (2004) — special hairpin loops (primary source for the bonus tables)

**URL:** https://doi.org/10.1146/annurev.biophys.32.110601.141800
(PDF: https://users.cs.duke.edu/~reif/courses/molcomplectures/DNA.Thermodynamics&Kinetics/Annu._Rev._Biophys._Biomol._Struct._2004_SantaLucia_Jr.pdf)
**Accessed:** 2026-06-25 (reused from PRIMER-TM-001-HAIRPIN evidence; §"Hairpin Loops", Eqs. 8-9)
**Authority rank:** 1 (peer-reviewed review)

**Key Extracted Points:**

1. Length-3 hairpins add a **triloop bonus** + a closing-A·T penalty; length-4 hairpins add a
   **tetraloop bonus** + a terminal-mismatch increment (Eqs. 8-9). The triloop/tetraloop bonus values
   are tabulated in the supplementary material — the same values libprimer3 ships in
   `triloop.*`/`tetraloop.*`. The article body confirms the model; the numeric tables come from the
   libprimer3 files retrieved above.

### primer3-py 2.3.0 — `calc_hairpin` (ground-truth numbers)

**URL:** `pip3 install --user primer3-py` (version 2.3.0, installed in prior PRIMER-TM work)
**Accessed:** 2026-06-25 (`primer3.calc_hairpin(seq, mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=50)`)
**Authority rank:** 3

**Key Extracted Points (captured this session):**

| Sequence | loop | dH (cal/mol) | dS (cal/K/mol) | dG (cal/mol) | Tm (°C) |
|----------|------|-------------:|---------------:|-------------:|--------:|
| GGGGCGAAAGCCCC | tetraloop CGAAAG | −40900 | −114.1872884299936 | −5484.812493437487 | 85.03347700825856 |
| GGGGGGGACCCCC | tetraloop GGGGAC | −34000 | −94.1872884299836 | −4787.81249344059 | 87.8328944728006 |
| GGGCGAAGCCC | triloop CGAAG | −27800 | −77.68485895331574 | −3706.040995629125 | 84.7060915802943 |
| GGGGGAACCCC | triloop GGAAC | −26000 | −73.18485895331571 | −3301.715995629134 | 82.11474153055735 |
| GGGCTTTTGCCC | non-special 4-nt | −32400 | −94.58485895332572 | −3064.5059956260297 | 69.39954078842845 |
| GGGCAAAAAGCCC | non-special 5-nt | −30100 | −87.74485895332572 | −2885.9319956260306 | 69.89004085311882 |
| AAAAAAAAAAAA | (none) | structure_found = False | | | |
| GCGC | (none) | structure_found = False | | | |

---

## Documented Corner Cases and Failure Modes

### From thal.c / primer3-py

1. **Only 3-nt and 4-nt loops have special tables.** A 5-nt-or-longer loop never gets a special bonus
   (bonus applies only for `loopSize == 3` / `loopSize == 4`). Non-special 3/4-nt loops fall through
   the bsearch with no match → no bonus.
2. **Key includes the closing pair.** A loop with the right 3/4-nt sequence but a different closing
   base pair is NOT in the table → no bonus.
3. **No structure** (homopolymer / too short) → ntthal `no_structure` (`structure_found = False`).

---

## Test Datasets

### Dataset: primer3-py 2.3.0 `calc_hairpin` reference (mv=50, dv=0, dntp=0, dna_conc=50 nM)

**Source:** primer3-py 2.3.0 (this session); see the table above.

| Parameter | Value |
|-----------|-------|
| tetraloop CGAAAG dH/dS/Tm | −40900 / −114.1872884299936 / 85.03347700825856 |
| triloop CGAAG dH/dS/Tm | −27800 / −77.68485895331574 / 84.7060915802943 |
| non-special TTTT dH/dS/Tm | −32400 / −94.58485895332572 / 69.39954078842845 |

---

## Assumptions

1. **ASSUMPTION: none correctness-affecting.** The bonus tables, their keying (full loop string incl.
   closing pair), the cal/mol-ΔH / cal-K-mol-ΔS convention, and their additive application are all
   taken verbatim from the libprimer3 config files + `thal.c`, and cross-verified to machine precision
   against primer3-py `calc_hairpin`. No value is invented or recalled.

---

## Recommendations for Test Coverage

1. **MUST Test:** a recognised tetraloop (CGAAAG) reproduces primer3 `calc_hairpin` ΔH/ΔS/ΔG/Tm
   (the bonus −1100 cal/mol is present) — Evidence: tetraloop.dh + thal.c calc_hairpin + primer3-py.
2. **MUST Test:** a recognised triloop (CGAAG) reproduces primer3 (bonus −2000 cal/mol) — Evidence: triloop.dh + primer3-py.
3. **MUST Test:** a non-special 4-nt loop (TTTT) is unchanged and still matches primer3 (regression) —
   Evidence: bsearch no-match + primer3-py.
4. **MUST Test:** a homopolymer / too-short input → null (ntthal no_structure) — Evidence: primer3-py structure_found=False.
5. **SHOULD Test:** a second recognised tetraloop (GGGGAC) and triloop (GGAAC) — Rationale: a second table row each.
6. **SHOULD Test:** the legacy Table-4 `FindMostStableHairpin` and its caller-supplied
   `loopBonusDeltaG37` path are unchanged — Rationale: opt-in / no-regression guarantee.

---

## References

1. SantaLucia J Jr, Hicks D (2004). The Thermodynamics of DNA Structural Motifs. Annu Rev Biophys
   Biomol Struct 33:415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
2. Untergasser A et al. (2012). Primer3 — new capabilities and interfaces. Nucleic Acids Res 40:e115.
   https://doi.org/10.1093/nar/gks596
3. libprimer3 `primer3_config/{triloop,tetraloop,loops}.{dh,ds}` and `thal.c`, vendored in primer3-py:
   https://github.com/libnano/primer3-py/tree/master/primer3/src/libprimer3 (GPL-2.0).
4. primer3-py 2.3.0 — https://pypi.org/project/primer3-py/2.3.0/

---

## Change History

- **2026-06-25**: Initial documentation (bundled special tri/tetraloop hairpin bonus tables + full
  ntthal hairpin DP for primer3 `calc_hairpin` parity).
