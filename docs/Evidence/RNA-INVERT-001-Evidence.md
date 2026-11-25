# Evidence Artifact: RNA-INVERT-001

**Test Unit ID:** RNA-INVERT-001
**Algorithm:** RNA Inverted Repeats (potential stem regions)
**Date Collected:** 2026-06-14 (revised 2026-09-25)

---

## Online Sources

### IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences (Alamro, Alzamel, Iliopoulos, Pissis, Watts; BMC Bioinformatics 2021)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/
**Accessed:** 2026-06-14 (fetched the PMC full-text page; query "inverted repeat definition reverse complement EMBOSS einverted stem loop")
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Core definition:** "An inverted repeat (IR) is a single stranded sequence of nucleotides with a subsequent downstream sequence consisting of its reverse complement."
2. **Formal model (gapped IR):** "a string that can be expressed in the form WGW̄ᴿ for some pair of strings W and G where |G|≥0" — W is the left arm, G the gap/spacer/loop (length ≥ 0), W̄ᴿ the right arm equal to the reverse complement of W.
3. **Mismatch model:** "a gapped inverted repeat within k mismatches when it can be expressed in the form WGW̄ᴿ with δH(W,W̄ᴿ)≤k", where δH is the Hamming distance (count of differing positions). A perfect IR is the case k = 0.
4. **Watson-Crick basis:** "the natural choice of a complement function over the alphabet Σ={A,C,G,T} is such that A⟷T and C⟷G" (for RNA, A⟷U and C⟷G).
5. **Tool parameters:** IUPACpal accepts minimum length, maximum length, maximum gap size, and maximum mismatches allowed.
6. **Worked example:** the sequence `CT-CGCAGTCACCG-GA` is an IR with left arm W = `CT`, an internal gap, and right arm = reverse complement of the left arm, with a single mismatch toward the tail ends.

### Wikipedia: Inverted repeat (definition cites Ussery et al. 2008)

**URL:** https://en.wikipedia.org/wiki/Inverted_repeat
**Accessed:** 2026-06-14 (fetched the article; primary citation is Ussery, Wassenaar & Borini 2008)
**Authority rank:** 4 (Wikipedia citing a primary textbook source — Ussery et al. 2008)

**Key Extracted Points:**

1. **Definition:** an inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement"; the intervening sequence may be any length including zero.
2. **Worked example:** `5'---TTACGnnnnnnCGTAA---3'` is an inverted repeat, where `nnnnnn` is the (any) intervening loop. The right arm `CGTAA` is the reverse complement of the left arm `TTACG` (TTACG → complement AATGC → reverse CGTAA).
3. **Palindrome relationship:** with zero intervening nucleotides the composite `5'TTACGCGTAA3'` is a palindromic sequence (equal to its own reverse complement).

### EMBOSS einverted manual (EMBOSS application documentation)

**URL:** https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted
**Accessed:** 2026-06-14 (fetched the application help page)
**Authority rank:** 3 (reference implementation / established tool, EMBOSS suite)

**Key Extracted Points:**

1. **Purpose:** "einverted finds inverted repeats (stem loops) in nucleotide sequences. It identifies regions of local alignment of the input sequence and its reverse complement that exceed a threshold score." This confirms inverted repeats are equivalent to potential stem-loops, found by comparing a sequence to its reverse complement.
2. **Stem-loop = inverted repeat:** the two arms are complementary; the intervening region is the loop/bulge. Mismatches and gaps correspond to bulges in the stem-loop.
3. **Scoring (default, EMBOSS):** match +3, mismatch −4, gap penalty 12, minimum reportable score 50. (Not used here — the repository implementation reports perfect, ungapped arms only; see Assumptions.)

---

## Documented Corner Cases and Failure Modes

### From IUPACpal (PMC7866733)

1. **Zero gap (palindrome):** |G| ≥ 0; when |G| = 0 the IR is a palindrome (`WW̄ᴿ`). This repository's RNA variant requires a loop (minSpacing ≥ 1 by default 3), so true palindromes with no loop are out of scope of the default call.
2. **Minimum/maximum arm length:** arms shorter than the minimum length are not reported; arms are extended to their maximal complementary run.

### From Wikipedia / Ussery et al. (2008)

1. **Intervening length:** the loop may be any length including zero; for RNA hairpins a non-zero loop is biologically required (a hairpin loop typically needs ≥ 3 unpaired nucleotides).

---

## Test Datasets

### Dataset: Wikipedia inverted-repeat example (perfect arms)

**Source:** Wikipedia "Inverted repeat" (citing Ussery, Wassenaar & Borini 2008), example `5'---TTACGnnnnnnCGTAA---3'`.

Adapted to RNA with a concrete 6-nt loop (`AAAAAA`): `UUACGAAAAAACGUAA` (length 16).

| Parameter | Value |
|-----------|-------|
| Sequence | `UUACGAAAAAACGUAA` |
| Left arm W | `UUACG` (positions 0–4) |
| Loop G | `AAAAAA` (positions 5–10, length 6) |
| Right arm W̄ᴿ | `CGUAA` (positions 11–15) |
| Arm length | 5 |
| Right arm = revcomp(left arm)? | Yes: UUACG → complement AAUGC → reverse CGUAA |

### Dataset: Palindromic IR with minimal RNA loop

**Source:** Derived directly from the WGW̄ᴿ definition (IUPACpal) with a 3-nt loop.

| Parameter | Value |
|-----------|-------|
| Sequence | `GGCCAAAGGCC` (length 11) |
| Left arm W | `GGCC` (positions 0–3) |
| Loop G | `AAA` (positions 4–6, length 3) |
| Right arm W̄ᴿ | `GGCC` (positions 7–10) |
| Right arm = revcomp(left arm)? | Yes: GGCC → complement CCGG → reverse GGCC |

---

## Assumptions (revised 2026-09, B12 F22)

1. **Pairing rule = ViennaRNA `CanPair`** (A·U, G·C, G·U wobble; `T` read as `U`; N/IUPAC never pair) — Varani & McClain 2000 / Crick 1966 for G·U; RNA-PAIR-001 evidence. `allowWobble:false` restricts to A·U, G·C. *Superseded:* the 2026-06 "strict Watson-Crick + IUPAC complement" rule (it paired N·N, S·S, W·W and any non-IUPAC character with itself and had no G·U) — unsourced, removed.
2. **Maximal stems with EMBOSS `palindrome` nesting filter** (see REP-INV-001 Evidence): a stem is dropped only when contained in both arms of another stem. *Superseded:* "maximal arm + non-overlapping greedy reporting" — the EMBOSS/IUPACpal papers report every maximal stem; the greedy longest-first rule had no source and lost valid stems (e.g. `GGGGAAACCCCAAAGGGG`, `UUUUAAAAAAAA`).
3. **Mismatches optional** (`maxMismatches`): IUPACpal "k mismatches" (δH ≤ k) with EMBOSS `palindrome -nummismatches` end-trimming (interior mismatches only). *Superseded:* "perfect stems only".
4. **Gapped / scored (einverted, Gotoh-style DP; match +3, mismatch −4, gap −12, threshold 50)** is not re-implemented here: `RepeatFinder.FindInvertedRepeatsScored` (REP-INV-001) is the single implementation (no duplication).
5. **Loop bounds** `minSpacing/maxSpacing` = IUPACpal gap bounds; defaults 3/100 are repository conventions.

---

## Review 2026-09 (B12 F22) — numerical confirmation

- **Independent oracle:** Python brute force (every exact stem `(L,R,len)` whose pairs all pair, filtered by arm/loop window; reported = not contained in both arms of another stem) vs `RnaSecondaryStructure.FindInvertedRepeats` through a console harness: **8 seeds × 1 500 cases = 12 000 cases, 0 mismatches** (alphabets ACGU, ACGT, ACGUN, acguNR, ACG, GU, ACGTUacgtu; length 0–40; minLength 1–4; minSpacing 0/1/3; maxSpacing 3/5/100; both pairing modes; 6 000 non-empty results). The same brute force is re-implemented in the unit tests (C# oracle, 20 seeded cases × 300 strings).
- **Canonical engine:** wrapper == `RepeatFinder.FindInvertedRepeats` on 600 random ACGT strings × (k 0/1) × (wobble on/off), 0 differences.
- **Proof-by-example of the old rule's defects (old = greedy longest-first non-overlap):** `GGGGAAACCCCAAAGGGG` old `{(0,3,7,10,4)}` → new `{(0,3,7,10,4),(7,10,14,17,4)}`; `GGGGAAAAUUUU` old `{}` → new (wobble) `{(0,3,8,11,4),(1,4,8,11,4)}`; `NNNNAAANNNN` old `{(0,3,7,10,4)}` (N·N) → new `{}`.
- **k mismatches:** `GGAGGAAACCACC`, minLength 5: k = 0 → `{}`, k = 1 → `{(0,4,8,12,5)}` (central A·A mismatch).
- Source of the tool semantics re-opened this review: EMBOSS einverted/palindrome manuals (via REP-INV-001 Evidence), IUPACpal (PMC7866733).

## Recommendations for Test Coverage

1. **MUST Test:** A known IR (`UUACGAAAAAACGUAA`) returns arms at the exact positions with the right arm = reverse complement of the left arm (antiparallel pairing). — Evidence: Wikipedia/Ussery example; IUPACpal WGW̄ᴿ.
2. **MUST Test:** A palindromic IR with a 3-nt loop (`GGCCAAAGGCC`) returns left 0–3, right 7–10, length 4. — Evidence: IUPACpal WGW̄ᴿ definition.
3. **MUST Test:** Antiparallel (not parallel) pairing — a sequence whose arms match in parallel but not antiparallel must NOT be reported. — Evidence: IUPACpal "reverse complement" (reverse, not just complement).
4. **MUST Test:** No inverted repeat present → empty result. — Evidence: definition (no WGW̄ᴿ decomposition exists).
5. **SHOULD Test:** Arm shorter than minLength is not reported; arm exactly minLength is. — Rationale: minimum-length parameter boundary (IUPACpal).
6. **SHOULD Test:** Loop length below minSpacing / above maxSpacing excludes the IR. — Rationale: gap-size bounds (IUPACpal).
7. **COULD Test:** Maximal-arm extension — a perfect arm longer than minLength is reported at its full length, not truncated. — Rationale: einverted maximal local alignment.

---

## References

1. Alamro H, Alzamel M, Iliopoulos CS, Pissis SP, Watts S. 2021. IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/)
2. Ussery DW, Wassenaar TM, Borini S. 2008. Computing for Comparative Microbial Genomics: Bioinformatics for Microbiologists. Springer. https://doi.org/10.1007/978-1-84800-255-5 (definition accessed via Wikipedia "Inverted repeat": https://en.wikipedia.org/wiki/Inverted_repeat)
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite (einverted application). Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (einverted manual: https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-25**: B12 review F22 — stem model delegated to `RepeatFinder`; ViennaRNA pairing (G·U, no N pairing); greedy non-overlap removed; `maxMismatches`, `allowWobble`; oracle cross-check.
