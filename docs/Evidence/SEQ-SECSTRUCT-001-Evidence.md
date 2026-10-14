# Evidence Artifact: SEQ-SECSTRUCT-001

**Test Unit ID:** SEQ-SECSTRUCT-001
**Algorithm:** Protein Secondary Structure Prediction — Chou-Fasman conformational propensities (sliding-window profile)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wikipedia — Chou–Fasman method (cites primary Chou & Fasman papers)

**URL:** https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method
**Accessed:** 2026-06-13 (WebFetch of the article URL)
**Authority rank:** 4 (Wikipedia citing primaries: Chou & Fasman 1974)

**Key Extracted Points:**

1. **Method definition:** The Chou-Fasman method assigns each amino acid conformational parameters Pα, Pβ and Pt (α-helix, β-sheet, β-turn propensities) derived from relative frequencies of each residue in known structures.
2. **Helix nucleation rule:** "Four out of any six contiguous amino acids were sufficient to nucleate helix" with helix-former threshold 1.03.
3. **Sheet nucleation rule:** "Three out of any contiguous five were sufficient for a sheet" with threshold 1.00.
4. **Turn rule:** turn probability p(t) = p_t(j)·p_t(j+1)·p_t(j+2)·p_t(j+3); original cutoff 7.5e-3.
5. **Primary citations recorded by the article:** Chou PY, Fasman GD (1974) "Conformational parameters for amino acids…"; Chou PY, Fasman GD (1974) "Prediction of protein conformation", Biochemistry 13:222-245.

### Kelley bioinfo lecture — "Protein 2° Structure: Chou-Fasman Algorithm" (PDF)

**URL:** https://www.kelleybioinfo.org/algorithms/background/BCho.pdf
**Accessed:** 2026-06-13 (WebFetch saved the PDF locally; text extracted with `pdftotext -layout`)
**Authority rank:** 4 (academic lecture material; reproduces primary parameters)

**Key Extracted Points:**

1. **Parameter meaning:** P(a) = propensity in an α-helix, P(b) = propensity in a β-sheet, P(turn) = propensity in a turn, based on observed propensities in proteins of known structure; propensity = observed/expected (×100 in the integer convention).
2. **Helix algorithm (verbatim steps):** STEP 1 find 4 of 6 residues with P(a) > 100 (window 6); STEP 2 extend right until 4 contiguous residues with P(a) < 100; STEP 3 compute mean P(a) and P(b) over the region; STEP 4 region is an α-helix if length > 5 and P(a) > P(b).
3. **Worked values (verbatim, integer ×100 convention):** A: P(a)=142, P(b)=83; E: P(a)=151, P(b)=37; S: P(a)=77, P(b)=75; P: P(a)=55, P(b)=55; D: P(a)=101, P(b)=54.

### CSB|SJU (Jakubowski) — Chou-Fasman propensities table

**URL:** https://employees.csbsju.edu/hjakubowski/classes/ch331/protstructure/tablechoufas.htm
**Accessed:** 2026-06-13 (WebFetch of the table URL)
**Authority rank:** 4 (curated academic reproduction of Chou-Fasman 1978)

**Key Extracted Points (verbatim helix Pα / sheet Pβ for all 20, with former/breaker classes):**

1. Ala 1.42 / 0.83; Cys 0.70 / 1.19; Asp 1.01 / 0.54; Glu 1.51 / 0.37; Phe 1.13 / 1.38;
   Gly 0.57 / 0.75; His 1.00 / 0.87; Ile 1.08 / 1.60; **Lys 1.16 / 0.74**; Leu 1.21 / 1.30;
   Met 1.45 / 1.05; Asn 0.67 / 0.89; Pro 0.57 / 0.55; Gln 1.11 / 1.10; Arg 0.98 / 0.93;
   Ser 0.77 / 0.75; Thr 0.83 / 1.19; Val 1.06 / 1.70; Trp 1.08 / 1.37; Tyr 0.69 / 1.47.
2. **Conflict noted:** this source lists Lys Pα = 1.16; two other retrieved sources (Przytycka NCBI lecture and the reference implementation below) list Lys Pα = 1.14. See Assumptions / conflict resolution.

### Przytycka (NCBI/NLM) — "Protein secondary structure prediction" lecture (PDF)

**URL:** https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/download/lectures/CAMS_02_Prot_Sec_Str.pdf
**Accessed:** 2026-06-13 (WebFetch saved the PDF locally; text extracted with `pdftotext -layout`)
**Authority rank:** 4 (NCBI computational-biology lecture; reproduces primary parameters)

**Key Extracted Points (verbatim, "The Chow–Fasman propensities" Pa Pb Pt):**

1. Alanine 1.42 0.83 0.66
2. Asparagine 0.67 0.89 1.56
3. **Lysine 1.14 0.74 1.01**
4. Phenylalanine 1.13 1.38 0.60
5. Tyrosine 0.69 1.47 1.14

### Reference implementation — ravihansa3000/ChouFasman (ChouFasman.py)

**URL:** https://raw.githubusercontent.com/ravihansa3000/ChouFasman/master/ChouFasman.py
**Accessed:** 2026-06-13 (WebFetch of the raw source file)
**Authority rank:** 3 (reference implementation; values reproduced as the canonical Chou-Fasman 1978 set)

**Key Extracted Points (verbatim Pa / Pb / Pt, integer ×100 convention; full 20-residue set):**

| AA | Pa | Pb | Pt | AA | Pa | Pb | Pt |
|----|----|----|----|----|----|----|----|
| Ala | 142 | 83 | 66 | Leu | 121 | 130 | 59 |
| Arg | 98 | 93 | 95 | Lys | 114 | 74 | 101 |
| Asp | 101 | 54 | 146 | Met | 145 | 105 | 60 |
| Asn | 67 | 89 | 156 | Phe | 113 | 138 | 60 |
| Cys | 70 | 119 | 119 | Pro | 57 | 55 | 152 |
| Glu | 151 | 37 | 74 | Gln | 111 | 110 | 98 |
| Gly | 57 | 75 | 156 | Ser | 77 | 75 | 143 |
| His | 100 | 87 | 95 | Thr | 83 | 119 | 96 |
| Ile | 108 | 160 | 47 | Trp | 108 | 137 | 96 |
|     |     |    |    | Tyr | 69 | 147 | 114 |
|     |     |    |    | Val | 106 | 170 | 50 |

Turn position frequencies f(i)..f(i+3) are also tabulated in this source but are not used by the
sliding-window profile under test (only Pa/Pb/Pt are averaged).

### BMC Bioinformatics — "Improved Chou-Fasman method…" (PMC1780123)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1780123/
**Accessed:** 2026-06-13 (WebFetch of the PMC article)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Original CFM rules (verbatim summary):** nucleation predicted when 4 of 6 sequential residues are helix-formers (3 of 5 for strand); nucleation regions are extended in both directions until the average 4-peptide propensity drops below 1; a segment is a helix if ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩ (strand: > 1.05).
2. **Primary citation given by the paper:** Chou PY, Fasman GD. Prediction of protein conformation. Biochemistry. 1974;13:222–245.

---

### 2026-09 review — sources for the full Chou-Fasman assignment (opened 2026-09-28)

- **Chen, Gu & Huang (2006) BMC Bioinformatics 7(S4):S14**, PDF bundled in
  greg-gardner/Chou-Fasman-Protein-Prediction (git clone), text extracted with pypdf.
  "Methods" states the three CFM rules verbatim: nucleation 4 of 6 (helix) / 3 of 5 (strand);
  extension "along both directions … until the average 4-peptides propensities drops below 1";
  helix if ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩, strand if ⟨Pβ⟩ > 1.05 and ⟨Pβ⟩ > ⟨Pα⟩; in overlapped
  regions "the secondary structure conformation with higher average propensities is predicted".
- **β-turn rule** (Grinnell ExBioPy 7.5 text, quoted in hassan11196/Chou-Fasman): "If (1)
  p(t) > 0.000075; (2) the average value for P(turn) > 1.00 in the tetrapeptide; and (3) the
  averages for the tetrapeptide obey the inequality P(a-helix) < P(turn) > P(b-sheet), then a
  beta-turn is predicted"; p(t) = f(i)·f(i+1)·f(i+2)·f(i+3). Same predicate coded in
  ravihansa3000 `CF_find_turns` and kalliapap `Turn.checkConditions`.
- **1978 table incl. bend frequencies f(i)..f(i+3)**: ravihansa3000 `ChouFasman.py` (rows keyed by
  residue name, copied from prowl.rockefeller.edu/aainfo/chou.htm) and hassan11196 C++ (rows
  keyed by one-letter code) agree on all 20 × 7 values. NOTE: ravihansa3000 assigns the symbols
  'N' to the "Aspartic Acid" row and 'D' to the "Asparagine" row (swap); values follow the names
  (Asp 101/54/146, f 0.147/0.110/0.179/0.081; Asn 67/89/156, f 0.161/0.083/0.191/0.091), as in
  hassan11196. Anyone porting that script by letter gets N and D wrong.
- greg-gardner `amino_tables.py` / evanimenon use the **1974** Pα/Pβ table (15 proteins, e.g.
  A 1.45/0.97) and no turn parameters — not the 1978 set; not used.
- Blocked this session: prowl.rockefeller.edu, crdd.osdd.net, cs.princeton.edu, arxiv.org,
  biogem.org (egress 403); WebSearch snippets only.

### 2026-10 completeness audit (B03 F34) — graded 1978 rule set (opened 2026-10-14)

- **Prevelige & Fasman (1989)**, book chapter PDF hosted by the cfPred project
  (https://cfpred.sourceforge.net/preveligeandfasman1989.pdf, downloaded; scanned, pages rendered at
  300 dpi and read visually because the OCR layer is garbled):
  - **Table II** (29 proteins, "From Chou and Fasman (1977, 1978)") brackets the 1978 parameters into
    classes: Pα — Hα Glu 1.51 Met 1.45 Ala 1.42 Leu 1.21; hα Lys 1.16 Phe 1.13 Gln 1.11 Trp 1.08 Ile 1.08
    Val 1.06; Iα Asp 1.01 His 1.00; iα Arg 0.98 Thr 0.83 Ser 0.77 Cys 0.70; bα Tyr 0.69 Asn 0.67;
    Bα Pro 0.57 Gly 0.57. Pβ — Hβ Val 1.70 Ile 1.60 Tyr 1.47; hβ Phe 1.38 Trp 1.37 Leu 1.30 Cys 1.19
    Thr 1.19 Gln 1.10 Met 1.05; iβ Arg 0.93 Asn 0.89 His 0.87 Ala 0.83; bβ Ser 0.75 Gly 0.75 Lys 0.74;
    Bβ Pro 0.55 Asp 0.54 Glu 0.37 (no Iβ residue in 1978). All Pα/Pβ values equal the code table.
    Independent confirmation (WebSearch extract of the CSB|SJU table [5]): Glu H/B, Met H/h, Ala H/i,
    Leu H/h, Asp I/B, Pro B/B, Gly B/b, Ser i/b, Lys h/b.
  - **Table I** (1974b, 15 proteins) footnotes: "Hα strong α former; hα α former; Iα weak α former;
    iα α indifferent; bα α breaker; Bα strong α breaker. Iα assignments are also given to Pro and Asp
    (near the N-terminal helix) as well as Arg (near the C-terminal helix)." / "… bβ assignment is also
    given to Trp (near the C-terminal β region)."
  - **Rules 1-4 (Fasman 1985)**: "A cluster of four helical residues (Hα or hα) out of six … will initiate
    a helix … extended in both directions until sets of tetrapeptide breakers (⟨Pα⟩ < 1.00) are reached.
    Proline cannot occur in the inner helix or at the C-terminal helical end but can occur within the last
    three residues at the N-terminal end … Any segment that is at least six residues long with
    ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩ is predicted as helical." β: "three β formers out of five … terminated
    by tetrapeptide breakers (⟨Pβ⟩ < 1.00) … ⟨Pβ⟩ > 1.05 as well as ⟨Pβ⟩ > ⟨Pα⟩". Turn and overlap rules
    are the ones already implemented.
  - §V "data reduction" is a by-hand procedure (judging terminal class-i residues, multiples of four,
    overlapping turns) — explicitly not automated by the authors; not implemented.
- **CFSSP (Kumar 2013) / MIT OCW 20.320 PS2 (2009) Criteria 1-3** (WebSearch extracts; core.ac.uk,
  ocw.mit.edu, biogem.org, researchgate blocked by egress): "Weak helical residues (Iα) count as 0.5 hα
  (i.e. three hα and two Iα residues out of six could also nucleate a helix)"; "Helix formation is
  unfavorable if the segment contains ⅓ or more helix breakers (bα or Bα), or less than ½ helix formers";
  β: three of five, "unfavorable if … ⅓ or more β-sheet breakers (bβ or Bβ)"; "Pro, Asp(−), Glu(−)
  prefer the N-terminal helical end. His(+), Lys(+), Arg(+) prefer the C-terminal helical end."
- **Maclin & Shavlik UCI domain theory** (raw.githubusercontent graingert/PyPred, fetched): helix nucleus
  "at least four are helix formers and less than two are helix breakers (two helix high-indifferent amino
  acids count as a helix former)"; sheet "at least three sheet formers and less than two sheet breakers" —
  independent confirmation of the ½ count and the ≤ 1 breaker limit.
- **Runnable reference attempt**: cfPred (sourceforge `cfPred.zip`, `cfPred_New.zip`, a C port of
  Prevelige's PREDICT.c) ships Windows binaries only and prints the parameter/tetrapeptide table (the
  assignment is done by hand, by design). 23 GitHub Chou-Fasman ports cloned (greg-gardner, seongminlab,
  Alpha-Park, sau-arv-gul, moyankate, nipunlakshank, berserk-23115, vidushikapoor, 003733m, Oneyedeagle,
  grpndr, Arnav1215, hamzafarooq009, jseidel5, alickgf/ICF, ZainFormatics, ardahuseyinoglu, 745talha,
  kmn01, ahmedhesham47/PyAA, EishaMazhar, deepesh336, evanimenon): none implements graded classes
  (grep for breaker/indifferent/weak former). **No runnable graded reference exists** → expected outputs
  are hand-traced from the rules above (below) and cross-checked against an independent Python reference.

#### Hand traces (Graded1978; integer parameters ×100)

| Input | Trace | Graded | Standard |
|-------|-------|--------|----------|
| AAAHHS | A Hα×3 = 3, H Iα×2 = ½+½ → 4 formers, 0 breakers → nucleus; ΣPα 703 > 618, ΣPβ 498; sheet: no formers; turns ⟨Pt⟩ ≤ 0.9975 | HHHHHH | CCCCCC (His Pα = 1.00 is not > 1.00) |
| EEEENY | 4 Hα, N + Y = 2 bα → rejected | CCCCCC | HHHHHH |
| EEESRR | E×3 + Arg at offsets 4,5 (C-terminal zone) → Iα ½+½ = 4; ΣPα 726 > 618 | HHHHHH | CCCCCC |
| RRSEEE | Arg at N-terminal side stays iα → 3 formers; turns p(t) 5.9e-5 / 4.8e-5 | CCCCCC | CCCCCC |
| PNEEEE | Pro offset 0 → Iα; 4 Hα + ½, 1 breaker (N); ΣPα 728 > 618 | HHHHHH | HHHHHH |
| VVVSW | sheet: 3 Hβ, S bβ + W (offset 4) bβ = 2 breakers → rejected; helix: ΣPα 503 < 618 | CCCCC | EEEEE |
| WVVVS | W offset 0 = hβ → 4 formers, 1 breaker; ΣPβ 722 > 525 and > ΣPα 503 | EEEEE | EEEEE |
| EEEEEEP | nucleus 0..5; C-extension to Pro forbidden; window 1..6 has Pro at offset 5 | HHHHHHC | HHHHHHH |
| AEEPEEEEE | window 0 rejected (Pro offset 3); nuclei 1..6, 2..7, 3..8 extend to 1..8; N-extension to A would put Pro at helix position 4 | CHHHHHHHH | HHHHHHHHH |
| EEEEGG | G,G = 2 Bα → no helix; EEGG turn (p(t) 9.70e-5) | CCTTTT | HHTTTT |

Executed cross-check: independent Python reference (`cf_graded.py`, scratchpad) vs C# on 3000 random
sequences (lengths 1-60, P/R/W/E/A/V/D/H/N/Y-enriched): 3000/3000 identical; 1302/3000 differ from
`Standard`. Human haemoglobin α (P69905) locked (differs from Standard at residues 1-3 and 38).

### Stage B (2026-09-28) — executed cross-checks

1. Independent Python reference of rules 1-5 (`cf_ref.py`, scratchpad) vs the ravihansa3000
   script (py3 port, N/D symbols corrected): tables identical; turn start positions identical on
   203 sequences (3 bundled proteins + 200 random); raw nucleated+extended helix/strand regions
   identical on 303 sequences after correcting three documented quirks of that script (last
   window skipped, `> 100` instead of "drops below 1", N-side tetrapeptide omits the new residue);
   the unmodified script differs on 294/303 — only because of those quirks.
2. C# `PredictSecondaryStructureChouFasman` vs the Python reference: 3000/3000 random sequences
   (lengths 1-250, helix-/strand-biased compositions, lowercase, X/B/Z/*/- junk) identical.
3. Hand-derived boundary cases (PPPP p(t) = 7.098e-5 < 7.5e-5 → coil; NPDG turn; SGIAKQ
   ⟨Pα⟩ 1.015 ≤ 1.03 → coil; PTIQGQ ⟨Pβ⟩ 1.048 ≤ 1.05 → coil; EEEEEEVVVVVV overlap → strand)
   locked in unit tests.

## Documented Corner Cases and Failure Modes

### From Kelley bioinfo lecture / CSB|SJU

1. **Non-standard residues:** the parameter table covers exactly the 20 standard amino acids; symbols such as X, B, Z, `*` and gap characters have no defined propensity and must be excluded from the averages.
2. **Window vs sequence length:** nucleation windows (6 for helix, 5 for sheet) require the window to fit inside the sequence; a window larger than the sequence yields no scan positions.

### From Wikipedia

1. **Reliability caveat:** the original parameters were derived from a small, non-representative sample (29 proteins) and have limited accuracy (~50-60% Q3); the propensities are nonetheless the formally defined Chou-Fasman values.

---

## Test Datasets

### Dataset: Single-residue propensities (Chou-Fasman 1978)

**Source:** Przytycka NCBI lecture + ravihansa3000/ChouFasman (both retrieved 2026-06-13).

| Residue | Helix Pa | Sheet Pb | Turn Pt |
|---------|----------|----------|---------|
| A | 1.42 | 0.83 | 0.66 |
| E | 1.51 | 0.37 | 0.74 |
| V | 1.06 | 1.70 | 0.50 |
| N | 0.67 | 0.89 | 1.56 |
| K | 1.14 | 0.74 | 1.01 |

### Dataset: Window mean derivation (window = sequence length)

**Source:** Definitional — mean of the per-residue propensities above.

| Sequence | Window | Helix mean | Sheet mean | Turn mean |
|----------|--------|------------|------------|-----------|
| "AE" | 2 | (1.42+1.51)/2 = 1.465 | (0.83+0.37)/2 = 0.60 | (0.66+0.74)/2 = 0.70 |
| "AEV" | 3 | (1.42+1.51+1.06)/3 = 1.330 | (0.83+0.37+1.70)/3 = 0.9666… | (0.66+0.74+0.50)/3 = 0.6333… |

---

## Assumptions

1. **ASSUMPTION: Lysine Pα conflict resolved to 1.14, not 1.16** — Two retrieved academic
   sources disagree on the lysine α-helix parameter: CSB|SJU lists 1.16, while the Przytycka
   NCBI lecture and the ravihansa3000 reference implementation list 1.14. The value 1.14 is
   adopted because it is supported by two independent retrieved sources (one of them a
   reference implementation that lists the integer parameter 114) versus a single source for
   1.16, and 114 is the value consistent with the Chou-Fasman integer convention used by the
   reference implementation for all other residues. The remaining 19 residues are identical
   across all retrieved sources, so this is the only contested value.
2. **ASSUMPTION: default window size = 7** — Chou-Fasman defines a 6-residue helix nucleation
   window and a 5-residue sheet nucleation window, not a single 7-residue averaging window.
   The method under test is a generic sliding-window *mean-propensity profile*, and the window
   length is a caller-supplied parameter; the default of 7 is an API convenience, not a
   Chou-Fasman constant. Test expectations therefore pass the window explicitly and verify the
   arithmetic mean over that window rather than the published nucleation/extension state machine.
3. **ASSUMPTION: unknown-residue handling = skip-and-exclude** — No retrieved source specifies
   behaviour for non-standard residues inside a window. The implementation excludes them from
   the per-window count/average; a window of only unknown residues emits nothing. This is the
   documented, deterministic contract for this profile method.

---

## Recommendations for Test Coverage

1. **MUST Test:** single-residue window equals that residue's (Pa, Pb, Pt) tuple for A, E, V — Evidence: Przytycka / reference-impl parameter table.
2. **MUST Test:** Lysine window returns Pa = 1.14 (the conflict-resolved value), not 1.16 — Evidence: Przytycka + reference impl.
3. **MUST Test:** multi-residue window returns the exact arithmetic mean of member propensities ("AE", "AEV") — Evidence: definition of the profile + parameter table.
4. **MUST Test:** sliding window steps by one and yields (n − w + 1) windows in N-terminus order — Evidence: window-scan definition (Kelley lecture STEP 1/2).
5. **MUST Test:** case-insensitive (lowercase equals uppercase) — Rationale: implementation uppercases input.
6. **MUST Test:** unknown residues excluded from the average ("AXE" with window 3 averages only A and E) — Evidence: ASSUMPTION 3 contract.
7. **SHOULD Test:** null / empty input → empty result — Rationale: documented precondition.
8. **SHOULD Test:** window larger than sequence → empty result — Evidence: window-vs-length corner case.
9. **SHOULD Test:** non-positive window → empty result — Rationale: validation contract.
10. **COULD Test:** helix-favouring vs sheet-favouring peptides show mean Pa > mean Pb (and vice versa) — Rationale: qualitative biological sanity, derived from exact means.

---

## References

1. Chou PY, Fasman GD (1978). Empirical predictions of protein conformation. Annual Review of Biochemistry 47:251-276. https://pubmed.ncbi.nlm.nih.gov/354496/
2. Chou PY, Fasman GD (1974). Prediction of protein conformation. Biochemistry 13(2):222-245. (cited by Wikipedia and PMC1780123, retrieved 2026-06-13)
3. Wikipedia. Chou–Fasman method. https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method (accessed 2026-06-13)
4. Kelley bioinfo. Protein 2° Structure: Chou-Fasman Algorithm. https://www.kelleybioinfo.org/algorithms/background/BCho.pdf (accessed 2026-06-13)
5. Jakubowski H. Chou-Fasman propensities (CSB|SJU CH331). https://employees.csbsju.edu/hjakubowski/classes/ch331/protstructure/tablechoufas.htm (accessed 2026-06-13)
6. Przytycka T. Protein secondary structure prediction (NCBI/NLM lecture). https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/download/lectures/CAMS_02_Prot_Sec_Str.pdf (accessed 2026-06-13)
7. ravihansa3000. ChouFasman reference implementation (ChouFasman.py). https://raw.githubusercontent.com/ravihansa3000/ChouFasman/master/ChouFasman.py (accessed 2026-06-13)
8. Chen H, Gu F, Huang Z (2006). Improved Chou-Fasman method for protein secondary structure prediction. BMC Bioinformatics 7(Suppl 4):S14. https://pmc.ncbi.nlm.nih.gov/articles/PMC1780123/

---

## Change History

- **2026-06-13**: Initial documentation (SEQ-SECSTRUCT-001).
- **2026-09-28**: Review 2026-09 (B03): full Chou-Fasman assignment sources, 1978 bend frequencies, ravihansa3000 N/D symbol swap noted, Stage B cross-checks.
- **2026-10-14**: B03 F34: graded 1978 rule set (Prevelige & Fasman 1989 Table II classes, Table I footnotes, rules 1-4; CFSSP/MIT ½-count and breaker limit; UCI domain theory), hand traces, Python cross-check 3000/3000.
