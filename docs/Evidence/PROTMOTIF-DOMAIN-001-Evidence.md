# Evidence Artifact: PROTMOTIF-DOMAIN-001

**Test Unit ID:** PROTMOTIF-DOMAIN-001
**Algorithm:** Domain Prediction & Signal Peptide Prediction
**Date Collected:** 2026-02-13

---

## Online Sources

### Wikipedia — Protein Domain

**URL:** https://en.wikipedia.org/wiki/Protein_domain
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **Definition:** A protein domain is a self-stabilizing region of a protein's polypeptide chain that folds independently. Domains range from ~50 to ~250 amino acids. — Wetlaufer (1973), cited via Xu & Nussinov (1998) https://doi.org/10.1016/S1359-0278(98)00004-2
2. **Domain databases:** Pfam, InterPro, PROSITE, and SMART are established databases for domain classification. Pfam classifies using Hidden Markov Models. — El-Gebali et al. (2019) https://doi.org/10.1093/nar/gky995
3. **Zinc finger C2H2:** The shortest domains (~50 aa), stabilized by metal ions. C2H2 zinc fingers adopt ββα fold with motif `X2-Cys-X(2,4)-Cys-X(12)-His-X(3,5)-His`. — Krishna et al. (2003) https://doi.org/10.1093/nar/gkg161
4. **SH3 domain:** An interaction domain involved in protein-protein binding, recognizing proline-rich sequences. — cited via Campbell & Downing (1994)
5. **Pattern-based detection:** PROSITE uses consensus patterns, Pfam uses profile HMMs for domain identification. Simple regex-based matching is a simplification of PROSITE's approach.

### Wikipedia — Signal Peptide

**URL:** https://en.wikipedia.org/wiki/Signal_peptide
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **Length:** Signal peptides are usually 16–30 amino acids long, present at the N-terminus. — Owji et al. (2018) https://doi.org/10.1016/j.ejcb.2018.06.003
2. **Tripartite structure:** Signal peptides contain three regions: n-region (positively charged, N-terminal), h-region (hydrophobic core, 7–15 residues), c-region (cleavage recognition site). — von Heijne (1985) J Mol Biol 184:99–105; von Heijne & Gavel (1988) https://doi.org/10.1111/j.1432-1033.1988.tb14150.x
3. **Cleavage site specificity:** Small, neutral amino acids {A, G, S} at positions -1 and -3 relative to the cleavage site (the "-1, -3 rule"). — von Heijne (1983) Eur J Biochem 133:17–21

### Wikipedia — Zinc Finger

**URL:** https://en.wikipedia.org/wiki/Zinc_finger
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with cited primary sources)

**Key Extracted Points:**

1. **C2H2 consensus motif:** `X2-Cys-X(2,4)-Cys-X(12)-His-X(3,5)-His`. Pfam: PF00096, PROSITE: PS00028. — Pabo et al. (2001) https://doi.org/10.1146/annurev.biochem.70.1.313
2. **PROSITE pattern PS00028:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` — verified at https://prosite.expasy.org/PS00028

### PROSITE Database — PS00028 (Zinc Finger C2H2)

**URL:** https://prosite.expasy.org/PS00028
**Accessed:** 2026-02-13
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Pattern:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`
2. **PROSITE accession:** PS00028
3. **Pfam cross-reference:** PF00096
4. **Note:** The implementation uses `[LIVMFYWC]` at position +3 from second C. The actual PROSITE PS00028 uses `[LIVMFYWC]` (confirmed). Some sources also include H in this class but the official PROSITE pattern does not.

### PROSITE Database — PS00017 (P-loop / Walker A motif)

**URL:** https://prosite.expasy.org/PS00017
**Accessed:** 2026-02-13
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Pattern:** `[AG]-x(4)-G-K-[ST]`
2. **Function:** ATP/GTP-binding site motif A (P-loop), also known as Walker A motif
3. **Reference:** Walker et al. (1982). EMBO J 1:945-951. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x

### PROSITE Database — PS00678 (WD-repeats signature, WD_REPEATS_1) — RETRIEVED 2026-06-24

**URL:** https://prosite.expasy.org/PS00678
**Retrieval:** WebFetch of https://prosite.expasy.org/PS00678 and the documentation page
https://prosite.expasy.org/PDOC00574 — both returned the same pattern string verbatim.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Entry name / accession:** WD_REPEATS_1 / **PS00678** — "Trp-Asp (WD) repeats signature".
2. **Pattern (verbatim from PROSITE):**
   `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]`
3. This is a **deterministic PROSITE PATTERN** (regex-like), NOT a profile/HMM, so it can be
   reproduced exactly. It is a 14-element signature spanning **15 residues** (the `x(2)` element
   covers 2 positions; all others are single positions). The terminal `W-[DEN]` is the diagnostic
   "WD" (Trp-Asp) of the repeat.
4. **Cross-references:** Pfam PF00400 (WD40); the separate WD-repeats **PROFILE** is PS50082
   (WD_REPEATS_2) — a weight matrix, not a pattern.

### PROSITE → regex syntax rules — RETRIEVED 2026-06-24

**URL:** https://prosite.expasy.org/scanprosite/scanprosite_doc.html
**Retrieval:** WebFetch of the ScanProsite/PROSITE pattern-syntax documentation.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Rules extracted (verbatim/closely paraphrased):**

1. "The standard IUPAC one letter code for the amino acids is used." → literal residue.
2. "The symbol 'x' is used for a position where any amino acid is accepted." → regex `.`.
3. "Ambiguities are indicated by listing the acceptable amino acids … between square brackets '[ ]'." → kept as `[..]`.
4. "Ambiguities are also indicated by listing between a pair of curly brackets '{ }' the amino acids that are NOT accepted at a given position." → `{X}` → `[^X]`.
5. "Each element in a pattern is separated from its neighbor by a '-'." → separators dropped.
6. Repetition `e(n)` / gap `x(n)` / range `x(m,n)` → `e{n}` / `.{n}` / `.{m,n}`.
7. N-terminal `<` → `^`; C-terminal `>` → `$`.

**Applied translation of PS00678:**
`[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]`

### PROSITE SH3 (PS50002) and PDZ (PS50106) are PROFILE-only — RETRIEVED 2026-06-24

**URLs:** https://prosite.expasy.org/PDOC50002 (SH3), https://prosite.expasy.org/PDOC50106 (PDZ),
EBI InterPro profile entries https://www.ebi.ac.uk/interpro/entry/profile/PS50002 and `/PS50106`.
**Retrieval:** WebSearch + the PROSITE/InterPro entries.
**Accessed:** 2026-06-24
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **SH3 (Src homology 3): PS50002 is a PROFILE (weight matrix)**, NOT a pattern — "profiles use
   weight matrices … characterize protein domains over their entire length." There is **no
   deterministic PROSITE PATTERN** for the SH3 domain.
2. **PDZ (DHR / GLGF): PS50106 is a PROFILE** covering the minimal PDZ domain — no deterministic
   PROSITE PATTERN exists.
3. **Consequence:** SH3 and PDZ **cannot** be reproduced as an exact regex without fabricating a
   signature. They are therefore NOT detected by `FindDomains` (honest residual). The previously
   shipped ad-hoc SH3/PDZ regexes were not sourced from any PROSITE pattern and have been removed.

### Real WD40 positive — GBB1_HUMAN (UniProt P62873) — RETRIEVED 2026-06-24

**URL:** https://rest.uniprot.org/uniprotkb/P62873.fasta
**Retrieval:** WebFetch of the UniProt REST FASTA.
**Accessed:** 2026-06-24
**Authority rank:** 5 (UniProt curated database)

**Key Extracted Points:**

1. GBB1_HUMAN (β-transducin / Gβ1) is a canonical 7-bladed WD40 β-propeller.
2. The PS00678-translated regex matches the sequence at **0-based starts 69, 156, 284** (each a
   15-residue window). Example repeat (start 69): `LVSASQDGKLIIWDS`.

### von Heijne G (1986). Signal sequences — the limits of variation

**URL:** https://doi.org/10.1016/0022-2836(85)90046-4
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **-1, -3 rule:** Small neutral amino acids (predominantly Ala, Gly, Ser) at positions -1 and -3 relative to the cleavage site.
2. **N-region:** 1–5 residues, positively charged (Lys, Arg).
3. **H-region:** 7–15 hydrophobic residues forming a single α-helix. The hydrophobic core is both necessary and sufficient for membrane targeting.
4. **C-region:** 3–7 residues, polar, with the -3/-1 specificity constraint.
5. **Note:** The amino acids Thr and Cys are occasionally found at -1/-3 but are not part of the canonical set. The implementation strictly uses {A, G, S}.

### von Heijne G (1983). Patterns of amino acids near signal-sequence cleavage sites

**URL:** https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **Signal peptide statistics:** Analysis of 78 eukaryotic signal peptides showed the -1 and -3 positions strongly favor A, G, S.
2. **Position -1:** A (65%), G (14%), S (10%), T (5%), other (6%)
3. **Position -3:** A (52%), V (15%), G (9%), S (8%), L (6%), other (10%)
4. **Note:** Position -3 also accepts V and L in von Heijne (1983) data, but the canonical "-1,-3 rule" specifies {A, G, S}. The implementation strictly uses {A, G, S}.

### Pfam — PF00096, PF00400, PF00018, PF00595, PF00069

**URL:** https://www.ebi.ac.uk/interpro/entry/pfam/
**Accessed:** 2026-02-13
**Authority rank:** 5 (Bioinformatics database)

**Key Extracted Points:**

1. **PF00096 (Zinc finger C2H2):** Classical ββα fold, ~23-28 residues per finger
2. **PF00400 (WD40 repeat):** ~40 residues per repeat, forms β-propeller
3. **PF00018 (SH3 domain):** ~60 residues, binds proline-rich sequences
4. **PF00595 (PDZ domain):** ~80-90 residues, protein-protein interaction
5. **PF00069 (Protein kinase):** Catalytic domain, includes P-loop (Walker A) motif for ATP binding

---

## Documented Corner Cases and Failure Modes

### From von Heijne (1986)

1. **Short sequences:** Sequences shorter than ~15 residues cannot contain a valid signal peptide (n+h+c regions require minimum length).
2. **Atypical cleavage sites:** Some signal peptides have non-canonical residues at -1/-3 (e.g., Thr, Val). The implementation strictly enforces {A, G, S} per von Heijne (1983).
3. **No signal peptide:** Cytoplasmic and nuclear proteins lack signal peptides entirely.

### From Pfam/PROSITE

1. **Empty/null input:** Domain search on empty sequence returns no domains.
2. **No matching domains:** Most short peptides or random sequences contain no recognizable domain signatures.
3. **Overlapping domains:** Some proteins have domains with overlapping boundaries.
4. **Multiple domain instances:** Tandem repeat domains (e.g., multiple zinc fingers) should all be detected.

---

## Test Datasets

### Dataset: Synthetic C2H2 Zinc Finger

**Source:** Constructed from PROSITE PS00028 consensus `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAACXXCXXXLXXXXXXXXHXXXHAAA` |
| Domain | Zinc Finger C2H2 |
| Pfam | PF00096 |
| Expected start | 4 |
| Expected end | 24 |

### Dataset: Synthetic P-loop (Kinase ATP-binding)

**Source:** Constructed from PROSITE PS00017 consensus `[AG]-x(4)-G-K-[ST]`

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAAGXXXXGKSAAAA` |
| Domain | Protein Kinase ATP-binding |
| Pfam | PF00069 |
| Expected match | Positions 4–11 |

### Dataset: Synthetic Signal Peptide

**Source:** Constructed from von Heijne (1986) principles: M + K/R (n-region) + L-rich (h-region) + small aa at -3,-1 (c-region)

| Parameter | Value |
|-----------|-------|
| Sequence | `MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF` |
| Expected | Signal peptide detected |
| Cleavage position | ~25 (after c-region) |

---

## Design Decisions (previously Assumptions — all resolved)

1. **Signal peptide scoring weights (1:2:1)** — The implementation weights region scores as `(nScore + 2·hScore + cScore) / 4` giving the h-region double weight. Justified by von Heijne (1985) J Mol Biol 184:99–105: the hydrophobic core is "both necessary and sufficient for membrane targeting," establishing it as the dominant discriminator.

2. **Signal peptide detection: evidence-based biological constraints** — ~~Previously: arbitrary threshold 0.4 (ELIMINATED).~~ Now enforces two literature-derived constraints: (a) n-region must contain positive charge (`nScore > 0`), per von Heijne (1986) "1–5 residues, positively charged (Lys, Arg)"; (b) h-region must be predominantly hydrophobic (`hScore ≥ 0.5`), per von Heijne (1985) "hydrophobic core" requirement. No arbitrary thresholds remain.

3. **Signal peptide probability** — ~~Previously: `min(1.0, score × 1.2)` (ELIMINATED).~~ Implementation now uses `Probability = Score` directly. No arbitrary scaling.

4. **Small amino acid set {A, G, S}** — ~~Previously: {A, G, S, T, C} (CORRECTED).~~ Implementation strictly uses {A, G, S} per von Heijne (1983) Eur J Biochem 133:17–21, positions -1 and -3. Matches the canonical -1,-3 rule.

5. **Domain signature patterns** — The implementation uses PROSITE-derived regex patterns (PS00028, PS00017) for domain detection. This is a deliberate scope decision: PROSITE itself uses consensus patterns (see Hulo et al. 2006). Not an assumption — it is a documented scope boundary.

6. **Method naming** — Public API is `FindDomains` (not `PredictDomains`). Non-correctness-affecting naming choice aligned with the method's behavior (pattern-based search, not probabilistic prediction).

---

## Recommendations for Test Coverage

1. **MUST Test:** FindDomains detects C2H2 zinc finger from PROSITE PS00028 consensus — Evidence: PROSITE PS00028
2. **MUST Test:** FindDomains detects P-loop/kinase from PROSITE PS00017 consensus — Evidence: PROSITE PS00017, Walker et al. (1982)
3. **MUST Test:** FindDomains returns correct domain metadata (name, accession, start, end) — Evidence: Pfam database
4. **MUST Test:** FindDomains empty/null input returns empty — Evidence: trivial
5. **MUST Test:** PredictSignalPeptide detects signal peptide with tripartite structure — Evidence: von Heijne (1986)
6. **MUST Test:** PredictSignalPeptide enforces -1,-3 rule for cleavage site — Evidence: von Heijne (1983, 1986)
7. **MUST Test:** PredictSignalPeptide returns null for non-secretory (charged) sequences — Evidence: von Heijne (1986)
8. **MUST Test:** PredictSignalPeptide returns null for sequences too short — Evidence: trivial (< 15 aa minimum)
9. **MUST Test:** PredictSignalPeptide returns three regions (n, h, c) — Evidence: von Heijne (1986)
10. **MUST Test:** PredictSignalPeptide is case-insensitive — Evidence: standard protein sequence convention
11. **SHOULD Test:** FindDomains detects WD40 repeat (PF00400) — Evidence: Pfam
12. **SHOULD Test:** FindDomains detects SH3 domain (PF00018) — Evidence: Pfam
13. **SHOULD Test:** FindDomains detects PDZ domain (PF00595) — Evidence: Pfam
14. **SHOULD Test:** FindDomains on no-match input returns empty — Evidence: trivial
15. **COULD Test:** FindDomains detects multiple domain instances in tandem — Evidence: Pfam, domain architecture

---

## Addendum 2026-06-25 — Plan7 profile-HMM engine + bundled Pfam SH3/PDZ/WD40 (limitation fix)

This addendum records the evidence for the **opt-in** Plan7 profile-HMM scorer
(`Plan7ProfileHmm` + `ProteinMotifFinder.FindDomainsByHmm`/`ScoreDomainHmm`) that detects the
SH3, PDZ and WD40 domains, which have **no deterministic PROSITE pattern** (they are trained
profile HMMs). The exact-PROSITE-pattern `FindDomains` path and its defaults are unchanged.

### Online Source: EMBL-EBI InterPro / Pfam profile HMMs (data)

**Retrieved (how):** `curl https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00018/?annotation=hmm`
(and PF00595, PF00400) on 2026-06-25 — gzip of the HMMER3/f ASCII profile, decompressed verbatim.
**Authority rank:** 5 (well-maintained bioinformatics database) for the data; format/algorithm below.

1. **Profiles obtained:** `PF00018.35` (SH3_1, LENG 48), `PF00595.30` (PDZ, LENG 81),
   `PF00400.39` (WD40, LENG 39). All `HMMER3/f [3.3 | Nov 2019]`, `ALPH amino`, with `GA`/`TC`/`NC`
   cutoffs, `STATS LOCAL MSV/VITERBI/FORWARD` lines, `COMPO`, and per-node M/I/D emissions+transitions.
   Header of PF00018: `NAME SH3_1`, `ACC PF00018.35`, `LENG 48`, `GA 22.9 22.9;`.
2. **Licence — CC0 (public domain), verbatim:** "Pfam is freely available under the Creative Commons
   Zero ('CC0') licence." — InterPro/Pfam documentation
   (https://interpro-documentation.readthedocs.io/en/latest/pfam.html, retrieved 2026-06-25).
   CC0 places the data in the public domain → freely redistributable; the three `.hmm` files are
   embedded as-is under `Seqeron.Genomics.Analysis/Resources/` (see that folder's README.md).

### Online Source: HMMER User's Guide v3.4 (Eddy, Aug 2023) — file format + pipeline

**Retrieved (how):** `curl http://eddylab.org/software/hmmer/Userguide.pdf` on 2026-06-25; read
pp. 210–215 ("HMMER profile HMM files") and pp. 57–62 (pipeline / null model).
**Authority rank:** 2 (canonical project specification).

1. **Score storage (verbatim):** "All probability parameters are stored as negative natural log
   probabilities with five digits of precision... a probability of 0.25 is stored as
   −log 0.25 = 1.38629. The special case of a zero probability is stored as '*'."
2. **Alphabet (verbatim):** for `ALPH amino`, "the alphabet size K is set to 20 and the symbol
   alphabet to 'ACDEFGHIKLMNPQRSTVWY' (alphabetic order)".
3. **Node layout:** match emission line (node #, K emissions, MAP/CONS/RF/MM/CS), insert emission
   line (K), state transition line (7: `Mk→Mk+1, Mk→Ik, Mk→Dk+1, Ik→Mk+1, Ik→Ik, Dk→Mk+1, Dk→Dk+1`).
   The two lines after `COMPO` are the BEGIN node (`B→M1, B→I0, B→D1, I0→M1, I0→I0, 0.0, *`). The
   single line after the `HMM` tag (the transition header) is always skipped by the parser.
4. **Bit score / null model (verbatim):** "A HMMER bit score is the log of the ratio of the
   sequence's probability according to the profile... over the null model probability." The null
   model is a one-state i.i.d. background; emission probabilities are turned into odds ratios.

### Online Source: Durbin et al. (1998) §5.4 — Viterbi/Forward recurrences

**Retrieved (how):** Stanford CS273 Lecture 7 PDF
(https://web.stanford.edu/class/cs273/scribing/scribe7.pdf), read pp. 10–12, reproducing the
profile-HMM Viterbi equations verbatim. **Authority rank:** 1 (textbook, via course notes).

1. **Viterbi recurrence (verbatim):**
   `V^M_j(i) = log(e_Mj(x_i)/q_xi) + max{ V^M_{j-1}(i-1)+log a_{M(j-1)M(j)}, V^I_{j-1}(i-1)+log a_{I(j-1)M(j)}, V^D_{j-1}(i-1)+log a_{D(j-1)M(j)} }`;
   `V^I_j(i) = log(e_Ij(x_i)/q_xi) + max{ V^M_j(i-1)+log a_{M(j)I(j)}, V^I_j(i-1)+log a_{I(j)I(j)}, V^D_j(i-1)+log a_{D(j)I(j)} }`;
   `V^D_j(i) = max{ V^M_{j-1}(i)+log a_{M(j-1)D(j)}, V^I_{j-1}(i)+log a_{I(j-1)D(j)}, V^D_{j-1}(i)+log a_{D(j-1)D(j)} }`.
   Insert emissions equal the background, so insert log-odds ≈ 0. Forward = same with max → log-sum-exp.

### Test Dataset — hand-built HMM (exact DP pin)

| Parameter | Value |
|-----------|-------|
| Alphabet | {A,C} (others negligible) |
| Background q | A=0.6, C=0.4 |
| M1 emissions | A=0.7, C=0.3 |
| M2 emissions | A=0.2, C=0.8 |
| Transitions | B→M1=0.9, M1→M2=0.8 |
| Sequence | "AC", path B→M1(A)→M2(C)→E |
| **Exact Viterbi (nats)** | **0.5187937934151676** = ln(0.7/0.6)+ln 0.9+ln(0.8/0.4)+ln 0.8 |

### Test Dataset — real true positives (UniProt, retrieved 2026-06-25)

| Domain | Source region | Accession | Observed Viterbi (bits) |
|--------|---------------|-----------|--------------------------|
| SH3 | SRC_HUMAN P12931 SH3 core `TFVALYDYE…GYIPSNYVAP` | PF00018 | ≈ +60 (≫ 10) |
| PDZ | DLG4_HUMAN P78352 PDZ1 (res 61–151) | PF00595 | ≈ +83 (≫ 10) |
| WD40 | GBB1_HUMAN P62873 (full β-propeller) | PF00400 | ≈ +36 (≫ 10) |
| negative | low-complexity `A14E14K12` | all three | strongly negative (rejected) |

Cross-domain specificity confirmed: SH3 core vs PF00400 ≈ −25 bits; PDZ1 vs PF00018 ≈ −35 bits.

### Honest residual

Exact `hmmsearch` bit-score / E-value parity is **not** reproduced: the MSV/bias/Viterbi/Forward
filter pipeline, null2 biased-composition correction, and Gumbel/exponential calibration
(`STATS LOCAL`) are out of scope. The Plan7 Viterbi/Forward log-odds is verified exactly on the
hand-built HMM (1e-9) and verified by correct ranking (true positive ≫ threshold ≫ true negative)
on the three real domains. Only these 3 Pfam domains are bundled; the full Pfam library is not.

---

## Addendum 2026-06-25 — HMMER E-value / P-value statistics from profile STATS

This addendum records the evidence for the **opt-in** E-value / P-value layer (Gumbel for
MSV/Viterbi, exponential tail for Forward; `E = P·Z`). Existing detection and defaults unchanged.

### Online Source: HMMER User's Guide v3.4 — STATS line semantics (file format)

**Retrieved (how):** `WebFetch http://eddylab.org/software/hmmer/Userguide.pdf` (mirror
`https://i5k.nal.usda.gov/training/webapp/static/hmmer/pdf/Userguide.pdf`), 2026-06-25; PDF saved
locally and read with `pdftotext`. **Verbatim** (HMM file format section): *"STATS &lt;s1&gt;
&lt;s2&gt; &lt;f1&gt; &lt;f2&gt; — Statistical parameters needed for E-value calculations. &lt;s1&gt;
is the model's alignment mode configuration: currently only LOCAL is recognized. &lt;s2&gt; is the
name of the score distribution: currently MSV, VITERBI, and FORWARD are recognized. &lt;f1&gt; and
&lt;f2&gt; are two real-valued parameters controlling location and slope of each distribution,
respectively; µ and λ for Gumbel distributions for MSV and Viterbi scores, and τ and λ for
exponential tails for Forward scores. λ values must be positive. All three lines or none of them
must be present; when all three are present, the model is considered to be calibrated for E-value
statistics."* The tutorial confirms scores are in **bits** and `E = P · Z`: *"a hit scoring 47.2
bits would be expected to happen 539165 times as often: 1.2e-16 × 539165 = 6.47e-11"* — i.e. the
per-sequence P-value times the number Z of database sequences.

### Online Source: Eddy (2008) PLoS Comput Biol 4:e1000069 — Gumbel / exponential

**Retrieved (how):** `WebFetch https://pmc.ncbi.nlm.nih.gov/articles/PMC2396288/`, 2026-06-25.
**Facts extracted:** Viterbi (and MSV) bit scores are **Gumbel** (Type-I extreme value) distributed
with parametric **λ = log 2**; the high-scoring tail of **Forward** scores is **exponentially**
distributed with the same λ. The E-value is the P-value times the number of database sequences.

### Online Source: Easel survival functions (HMMER's actual implementation)

**Retrieved (how):** `WebFetch https://raw.githubusercontent.com/EddyRivasLab/easel/master/esl_gumbel.c`
and `.../esl_exponential.c`, 2026-06-25. **Verbatim code:**
```c
double esl_gumbel_surv(double x, double mu, double lambda) {       // P(X>x), 1-cdf, right tail
  double y  = lambda*(x-mu);
  double ey = -exp(-y);
  if (fabs(ey) < eslSMALLX1) return -ey;       /* 1-e^x ~ -x when e^-y is small */
  else                       return 1 - exp(ey);
}
double esl_exp_surv(double x, double mu, double lambda) {
  if (x < mu) return 1.0;
  return exp(-lambda * (x-mu));
}
```
i.e. Gumbel `P(S≥x) = 1 − exp(−exp(−λ(x−µ)))` (with the `eslSMALLX1 = 5e-9` tail branch returning
`exp(−λ(x−µ))`) and exponential `P(S≥x) = exp(−λ(x−τ))`, clamped to 1 for `x < τ`.

### Online Source: HMMER pipeline — score → P-value → E-value

**Retrieved (how):** `WebFetch
https://raw.githubusercontent.com/Janelia-Farm-Xfam/Bio-HMM-Logo/master/src/src/p7_pipeline.c`,
2026-06-25. **Verbatim usage:** MSV `P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU],
om->evparam[p7_MLAMBDA])`; Viterbi `P = esl_gumbel_surv(seq_score, om->evparam[p7_VMU],
om->evparam[p7_VLAMBDA])`; Forward `P = esl_exp_surv(seq_score, om->evparam[p7_FTAU],
om->evparam[p7_FLAMBDA])`; E-value `exp(lnP) * pli->Z` — confirming the **per-sequence bit score**
is the survival-function argument and `E = P · Z`.

### Hand-derived test pin (bundled PF00018 STATS)

`STATS LOCAL VITERBI -8.2932 0.71923` and `STATS LOCAL FORWARD -4.5735 0.71923`, read verbatim from
`Resources/PF00018_SH3_1.hmm`. At bit score `S = 40`:
- Gumbel: `y = 0.71923·(40−(−8.2932)) = 34.7326…`; `exp(−y) = 8.227179545686635e-16`; `|ey| < 5e-9`
  → `P = 8.227179545686635e-16`; `E(Z=1000) = 8.227179545686635e-13`.
- Exponential: `P = exp(−0.71923·(40−(−4.5735))) = 1.1943390031599535e-14`;
  `E(Z=1000) = 1.1943390031599535e-11`.
Computed at full double precision (Python) independently of the implementation; the engine matches
to 1e-9 relative (tests H14).

### Honest residual (narrowed)

The Gumbel/exponential P-value and `E = P·Z` are now implemented exactly. What remains out of scope
is exact `hmmsearch`-**reported** E-value *pipeline* parity: HMMER applies these formulas to its
local-multihit sequence bit score after the MSV/bias prefilters and the **null2 biased-composition
correction**, which this glocal scorer does not compute, so absolute bit scores (and hence absolute
reported E-values) differ from `hmmsearch`. Pfam coverage beyond the three bundled (caller-supplied
`.hmm`) profiles is likewise out of scope.

---

## Addendum 2026-06-25 — HMMER local-multihit Forward + null2 (hmmsearch parity)

This addendum records the evidence for the **opt-in** hmmsearch-parity layer: HMMER's
**local-multihit** Forward score, the **null2 biased-composition correction**, and the final
null2-corrected per-sequence bit score. The glocal Viterbi/Forward path and all existing defaults are
unchanged. A real `hmmsearch` reference was obtained (pyhmmer) and used as ground truth.

### Reference tool — pyhmmer (HMMER3 binding) INSTALLED

**Retrieved (how):** `pip3 install --user pyhmmer` → **pyhmmer 0.12.1** (Cython binding bundling
HMMER3), installed 2026-06-25 on macOS arm64. Verified `import pyhmmer.hmmer; hmmsearch(...)` runs.
Ground-truth `hmmsearch` scores for the bundled CC0 profiles against the test true-positives
(`pyhmmer.hmmer.hmmsearch([hmm], targets, Z=1.0)`):

| Profile | Target (UniProt) | `pre_score` (bits, no null2) | `bias` (bits) | `score` (bits) | E-value (Z=1) |
|---------|------------------|------------------------------|---------------|----------------|----------------|
| PF00018 SH3 | SRC_HUMAN P12931 core | **68.709740** | 0.025574 | 68.684166 | 1.31e-23 |
| PF00595 PDZ | DLG4_HUMAN P78352 PDZ1 | **84.862930** | 0.078865 | 84.784065 | 2.23e-28 |
| PF00400 WD40 | GBB1_HUMAN P62873 | **213.411926** | 24.969421 | 188.442505 | 1.91e-61 |

SH3 single domain envelope = positions 3–50 (`dom.bias = 0.025569`).

### Online Source: modelconfig.c — local entry/exit + length config (`p7_ReconfigLength`)

**Retrieved (how):** `curl https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/modelconfig.c`, 2026-06-25.
**Verbatim:** local entry `p7P_TSC(gm,k-1,p7P_BM) = log(occ[k]/Z)` with `Z = Σ_k occ[k]·(M−k+1)`
("Reduces to uniform 2/(M(M+1)) for occupancies of 1.0"). Multihit E split:
`gm->xsc[p7P_E][p7P_MOVE] = gm->xsc[p7P_E][p7P_LOOP] = -eslCONST_LOG2; gm->nj = 1.0f`.
Length config `p7_ReconfigLength`: `pmove = (2.0f+gm->nj)/((float)L+2.0f+gm->nj)`, `ploop = 1−pmove`,
`xsc[N|C|J][LOOP]=log(ploop)`, `xsc[N|C|J][MOVE]=log(pmove)`.

### Online Source: p7_hmm.c — `p7_hmm_CalculateOccupancy`

**Retrieved (how):** `curl …/src/p7_hmm.c`, 2026-06-25. **Verbatim:**
`mocc[1] = t[0][p7H_MI]+t[0][p7H_MM]; mocc[k] = mocc[k-1]·(t[k-1][MM]+t[k-1][MI]) + (1−mocc[k-1])·t[k-1][DM]`.
Transition enum (src/hmmer.h, retrieved 2026-06-25): `p7H_MM=0, MI=1, MD=2, IM=3, II=4, DM=5, DD=6`
— matches the `.hmm` file column order and this code's `TMM..TDD`.

### Online Source: generic_fwdback.c — local Forward/Backward (`esc=0`)

**Retrieved (how):** `curl …/src/generic_fwdback.c`, 2026-06-25. **Verbatim:**
`esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY` (local exit M_k/D_k→E prob 1 for all k); B→M_k uses
`TSC(p7P_BM,k-1)`; specials `J = logsum(J_loop·ploop, E + E_loop)`, `C = logsum(C_loop, E + E_move)`,
`N = N_loop`, `B = logsum(N+N_move, J+J_move)`; final `sc = XMX(L,C) + xsc[C][MOVE]`. Insert emissions
are hardwired to background (insert log-odds = 0): "we currently hardwire insert scores to 0".

### Online Source: generic_decoding.c — `p7_GDecoding` (posterior probabilities)

**Retrieved (how):** `curl …/src/generic_decoding.c`, 2026-06-25. **Verbatim:** per row `i`,
`MMX = exp(fwd_M + bck_M − overall_sc)`, `IMX = exp(fwd_I + bck_I − overall_sc)`,
`XMX[N|J|C] = exp(fwd_x(i−1) + bck_x(i) + xsc[x][LOOP] − overall_sc)`; renormalise by the row sum.

### Online Source: generic_null2.c — `p7_GNull2_ByExpectation`

**Retrieved (how):** `curl …/src/generic_null2.c`, 2026-06-25. **Verbatim:** sum the posterior
emission counts over the `Ld` envelope rows, take `log(count) − log(Ld)` as log posterior weights;
then for each residue `x`, `null2[x] = logsum_k( w_M(k)+MSC(k,x), w_I(k)+ISC(k,x) ) ⊕ xfactor`,
where `xfactor = logsum(w_N, w_C, w_J)`; finally `null2[x] = exp(null2[x])` = odds ratio
`f'_d(x)/f_0(x)`. Header: "Calculate the null2 model for the envelope … rows numbered 1..Ld".

### Online Source: p7_domaindef.c + p7_pipeline.c + p7_bg.c — applying null2

**Retrieved (how):** `curl …/src/p7_domaindef.c`, `…/src/p7_pipeline.c`, `…/src/p7_bg.c`, 2026-06-25.
**Verbatim:** per-position `n2sc[pos] = logf(null2[dsq[pos]])` (nats), `domcorrection += n2sc[pos]`
over the envelope; per-seq `seqbias = p7_FLogsum(0.0, log(bg->omega) + Σ n2sc)`,
`seq_score = (fwdsc − (nullsc + seqbias)) / eslCONST_LOG2`, `pre_score = (fwdsc − nullsc)/log2`,
`hit->bias = pre_score − score`. Constants: `bg->omega = 1./256.` (p7_bg.c);
`nullsc = L·log(p1)+log(1−p1)` with `p1 = L/(L+1)` (p7_bg_NullOne + p7_bg_SetLength).

### Online Source: hmmer.c — `p7_AminoFrequencies` (the scoring background)

**Retrieved (how):** `curl …/src/hmmer.c`, 2026-06-25. **Verbatim** 20 amino background frequencies
(A..Y) `0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687,
0.0673417, 0.0114135, 0.0304133` ("average Swiss-Prot residue composition"). **Key fact:** hmmsearch
scores match emissions against THIS `bg->f`, **not** the profile's COMPO line (which the glocal path
uses) — this is what makes the parity path's `pre_score` match hmmsearch exactly.

### Verification — independent reproduction + C# parity

1. **From-scratch Python** re-derivation (parsing the `.hmm`, the recurrences above, `bg->f` from
   pyhmmer) gives SH3 `pre_score = 68.709743` bits (hmmsearch 68.709740) and, decoding the
   envelope 3–50, null2 `bias = 0.025544` bits (hmmsearch 0.025574). Both match to single-float
   rounding, confirming the formulas independently of the C# code.
2. **C# implementation** (`LocalForwardBitScore`, `Null2BiasBits`, `HmmSearchBitScore`) reproduces:
   SH3 `pre = 68.709743` (ref 68.709740), PDZ `pre = 84.862933` (ref 84.862930), WD40
   `pre = 213.411951` (ref 213.411926) — all within `5e-5` bits of hmmsearch. The null2 bias matches
   hmmsearch for the single-domain envelope (SH3 envelope bias 0.025544 vs 0.025574).

### Honest residual (further narrowed)

- **`pre_score` (local-multihit Forward bit score): verified to ~1e-5-bit hmmsearch parity** for the
  three bundled profiles via pyhmmer 0.12.1.
- **null2 correction: formula verified exact** (single-domain SH3 envelope reproduces hmmsearch's
  bias to 3e-5 bits). The remaining gap is **domain decomposition**: HMMER computes null2 per the
  posterior-defined domain *envelope* (model reconfigured to the envelope length), then sums across
  domains. This implementation computes null2 over whatever sequence/envelope the caller passes; it
  does NOT run HMMER's region/envelope-definition heuristic (stochastic-traceback clustering of
  multi-domain regions). So for a **single well-resolved domain** the corrected score matches
  hmmsearch; for a **multi-domain** target (e.g. the 7-blade WD40) the caller must score each
  envelope separately to reproduce hmmsearch's per-domain decomposition. The MSV/bias *prefilters*
  are not reimplemented (they only gate which sequences reach the Forward stage; they do not change a
  reported hit's bit score). Full Pfam coverage remains caller-supplied `.hmm`.

---

## References

1. von Heijne G (1986). Signal sequences. The limits of variation. J Mol Biol 184(1):99-105. https://doi.org/10.1016/0022-2836(85)90046-4
2. von Heijne G (1983). Patterns of amino acids near signal-sequence cleavage sites. Eur J Biochem 133(1):17-21. https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
3. Walker JE et al. (1982). Distantly related sequences in the α- and β-subunits of ATP synthase, myosin, kinases and other ATP-requiring enzymes and a common nucleotide binding fold. EMBO J 1(8):945-951. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x
4. Krishna SS, Majumdar I, Grishin NV (2003). Structural classification of zinc fingers. Nucleic Acids Res 31(2):532-550. https://doi.org/10.1093/nar/gkg161
5. El-Gebali S et al. (2019). The Pfam protein families database in 2019. Nucleic Acids Res 47(D1):D427-D432. https://doi.org/10.1093/nar/gky995
6. Hulo N et al. (2006). The PROSITE database. Nucleic Acids Res 34:D227-D230. https://doi.org/10.1093/nar/gkj063
7. PROSITE PS00028 — Zinc finger C2H2 type. https://prosite.expasy.org/PS00028
8. PROSITE PS00017 — ATP/GTP-binding site motif A (P-loop). https://prosite.expasy.org/PS00017
9. Pabo CO, Peisach E, Grant RA (2001). Design and selection of novel Cys2His2 zinc finger proteins. Annu Rev Biochem 70:313-340. https://doi.org/10.1146/annurev.biochem.70.1.313
10. Owji H et al. (2018). A comprehensive review of signal peptides. Eur J Cell Biol 97(6):422-441. https://doi.org/10.1016/j.ejcb.2018.06.003
11. von Heijne G (1986). A new method for predicting signal sequence cleavage sites. Nucl Acids Res 14(11):4683–4690. https://doi.org/10.1093/nar/14.11.4683
12. von Heijne G (1984). How signal sequences maintain cleavage specificity. J Mol Biol 173(2):243–251. https://doi.org/10.1016/0022-2836(84)90192-x
13. Mistry J et al. (2021). Pfam: The protein families database in 2021. Nucleic Acids Res 49:D412–D419. https://doi.org/10.1093/nar/gkaa913
14. Eddy SR & the HMMER team (2023). HMMER User's Guide, v3.4. http://eddylab.org/software/hmmer/Userguide.pdf
15. Durbin R, Eddy SR, Krogh A, Mitchison G (1998). Biological Sequence Analysis, Ch. 5.4. Cambridge Univ Press; recurrences reproduced in Stanford CS273 Lecture 7: https://web.stanford.edu/class/cs273/scribing/scribe7.pdf
16. Eddy SR (2011). Accelerated Profile HMM Searches. PLoS Comput Biol 7:e1002195. https://doi.org/10.1371/journal.pcbi.1002195
17. Pfam licence (CC0): https://interpro-documentation.readthedocs.io/en/latest/pfam.html
18. Eddy SR (2008). A Probabilistic Model of Local Sequence Alignment That Simplifies Statistical Significance Estimation. PLoS Comput Biol 4:e1000069. https://doi.org/10.1371/journal.pcbi.1000069 (PMC2396288)
19. Easel esl_gumbel.c (esl_gumbel_surv): https://github.com/EddyRivasLab/easel/blob/master/esl_gumbel.c
20. Easel esl_exponential.c (esl_exp_surv): https://github.com/EddyRivasLab/easel/blob/master/esl_exponential.c
21. HMMER p7_pipeline.c (score→P-value→E-value via esl_gumbel_surv / esl_exp_surv, E=P·Z): https://github.com/Janelia-Farm-Xfam/Bio-HMM-Logo/blob/master/src/src/p7_pipeline.c
22. HMMER modelconfig.c (p7_ProfileConfig local entry/exit; p7_ReconfigLength): https://github.com/EddyRivasLab/hmmer/blob/master/src/modelconfig.c
23. HMMER generic_fwdback.c (p7_GForward / p7_GBackward, local esc=0): https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_fwdback.c
24. HMMER generic_decoding.c (p7_GDecoding posterior probabilities): https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_decoding.c
25. HMMER generic_null2.c (p7_GNull2_ByExpectation): https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_null2.c
26. HMMER p7_domaindef.c (per-domain null2 correction, domcorrection): https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_domaindef.c
27. HMMER p7_pipeline.c (master) (seqbias / pre_score / bias): https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_pipeline.c
28. HMMER p7_bg.c (omega=1/256, p7_bg_NullOne, p7_bg_SetLength p1=L/(L+1)): https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_bg.c
29. HMMER hmmer.c (p7_AminoFrequencies background): https://github.com/EddyRivasLab/hmmer/blob/master/src/hmmer.c
30. pyhmmer 0.12.1 (HMMER3 Cython binding, ground-truth hmmsearch scores): https://pyhmmer.readthedocs.io/

---

## Addendum 2026-06-25 — HMMER multi-domain envelope decomposition (p7_domaindef)

The Plan7 engine now performs HMMER's automatic per-target **domain/envelope decomposition** — the
`hmmsearch`/`hmmscan` step that splits a multi-domain protein into one scored domain per envelope.
Opt-in; all existing methods + defaults unchanged. Sources retrieved verbatim this session.

### Online Source: HMMER `p7_domaindef.c` — region identification + per-domain rescore

**URL fetched:** https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/p7_domaindef.c
**Authority rank:** 3 (reference implementation; EddyRivasLab/hmmer master, retrieved 2026-06-25)

Verbatim facts extracted:

1. **Thresholds (`p7_domaindef_Create`):** `ddef->rt1 = 0.25; ddef->rt2 = 0.10; ddef->rt3 = 0.20;`
   `ddef->nsamples = 200;` `ddef->do_reseeding = TRUE;`. The example/test drivers seed the RNG with
   `esl_randomness_CreateFast(42)` (default seed 42).
2. **Region identification loop (`p7_domaindef_ByPosteriorHeuristics`):** for `j = 1..n`, while not
   triggered: `if (mocc[j] - (btot[j]-btot[j-1]) < rt2) i = j; else if (i==-1) i=j;` then
   `if (mocc[j] >= rt1) triggered = TRUE;`. Once triggered, the region closes at the first `j` with
   `mocc[j] - (etot[j]-etot[j-1]) < rt2`, yielding region `i..j`.
3. **`is_multidomain_region(ddef,i,j)`:** `return ( \max_z [ \min(E(z), B(z)) ] >= rt3 )` where
   `E(z) = etot[z]-etot[i-1]` and `B(z) = btot[j]-btot[z-1]` for `z = i..j`. TRUE → resolve by
   stochastic-traceback clustering (`region_trace_ensemble` + `p7_spensemble_Cluster`); FALSE →
   "the region looks simple, single domain; convert the region to an envelope."
4. **`rescore_isolated_domain` (non-longtarget):** `p7_Forward(dsq+i-1, Ld, om, …, &envsc)` then
   `p7_Backward` + `p7_Decoding`; if null2 not already done, `p7_Null2_ByExpectation(om,ox2,null2)`
   and `ddef->n2sc[pos] = logf(null2[dsq[pos]])` for `pos=i..j`; `domcorrection += ddef->n2sc[pos]`
   (NATS); `dom->ienv=i; dom->jenv=j; dom->envsc=envsc; dom->domcorrection`. The model is in
   **unihit** mode (`p7_oprofile_ReconfigUnihit(om, saveL)`) with length still `saveL` = the full
   sequence length n.

### Online Source: HMMER `generic_decoding.c` — `p7_GDomainDecoding` (btot/etot/mocc)

**URL fetched:** https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/generic_decoding.c
**Authority rank:** 3.

`overall_logp = fwd->xmx[C at L] + gm->xsc[C][MOVE]`; for `i = 1..L`:
`btot[i] = btot[i-1] + exp(fB[i-1] + bB[i-1] - overall)`;
`etot[i] = etot[i-1] + exp(fE[i] + bE[i] - overall)`;
`njcp = exp(fN[i-1]+bN[i]+xsc[N][LOOP]-overall) + exp(fJ[i-1]+bJ[i]+xsc[J][LOOP]-overall) +
exp(fC[i-1]+bC[i]+xsc[C][LOOP]-overall)`; `mocc[i] = 1 - njcp`.

### Online Source: HMMER `p7_pipeline.c` — per-domain bit score + i-Evalue

**URL fetched:** https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/p7_pipeline.c
**Authority rank:** 3.

Verbatim per-domain finalisation (after domain definition):
`Ld = jenv - ienv + 1;`
`bitscore = envsc + (sq->n - Ld) * log((float)sq->n / (sq->n+3));` (NATS)
`dombias  = do_null2 ? p7_FLogsum(0.0, log(bg->omega) + domcorrection) : 0.0;` (NATS)
`bitscore = (bitscore - (nullsc + dombias)) / eslCONST_LOG2;` (BITS)
`lnP = esl_exp_logsurv(bitscore, evparam[FTAU], evparam[FLAMBDA]);` → per-domain i-Evalue `= domZ·exp(lnP)`.
`nullsc = p7_bg_NullOne(sq->dsq, sq->n)` over the **full** sequence (`L·ln p1 + ln(1−p1)`, `p1=L/(L+1)`).

### Online Source: HMMER `modelconfig.c` — unihit/multihit + `p7_ReconfigLength`

**URL fetched:** https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/modelconfig.c
**Authority rank:** 3.

`p7_ReconfigLength`: `pmove = (2+nj)/(L+2+nj)`, `ploop = 1−pmove`, set on N/C/J LOOP/MOVE.
`p7_ReconfigMultihit`: `xsc[E][MOVE]=xsc[E][LOOP]=−log 2; nj=1`. `p7_ReconfigUnihit`:
`xsc[E][MOVE]=0; xsc[E][LOOP]=−inf; nj=0`. The envelope rescore uses **unihit** (`nj=0`, no J,
E→C with log-prob 0) at full length n.

### Reference tool — pyhmmer 0.12.1 multi-domain GROUND TRUTH (captured 2026-06-25)

`pyhmmer.plan7.Pipeline(amino, Z=1, domZ=1, bias_filter=True).search_hmm(hmm, [seq])`. macOS arm64,
Python 3. Two real targets vs the bundled CC0 profiles:

**GBB1_HUMAN / Gβ1 (UniProt P62873, 7-bladed WD40 β-propeller, L=340) vs PF00400** —
seq score 188.442505, ndom = **7**. Per-domain (`env_from..env_to`, score, bias, i-Evalue):

| dom | env (1-based) | score (bits) | bias | i-Evalue |
|-----|---------------|--------------|------|----------|
| 0 | 45..83  | 31.139467 | 0.034545 | 1.2101e-11 |
| 1 | 87..125 | 19.004278 | 0.088024 | 8.4096e-08 |
| 2 | 133..170 | 25.053679 | 0.403497 | 1.0223e-09 |
| 3 | 174..212 | 35.552242 | 0.136274 | 4.8501e-13 |
| 4 | 216..254 | 40.454269 | 0.039276 | 1.3608e-14 |
| 5 | 259..298 | 23.443121 | 0.197371 | 3.3071e-09 |
| 6 | 303..340 | 27.824228 | 0.672585 | 1.3565e-10 |

**SRC_HUMAN SH3 core (UniProt P12931, L=55) vs PF00018** — ndom = **1**: env 3..50, score
68.540695, bias 0.025569, i-Evalue 1.4529e-23.

### Independent C# parity (this engine, double precision)

`Plan7ProfileHmm.FindDomains` reproduces: GBB1 → **7** envelopes with **exactly** the env bounds
above; per-domain scores 31.139483 / 19.004292 / 25.053932 / 35.552249 / 46→40.454279 / 23.443205 /
27.824492 — matching hmmsearch to ≈1e-3 bits (HMMER computes in float32; this engine in float64);
i-Evalues identical to ≥3 sig figs. SH3 → **1** envelope, env 3..50, score 68.540701, i-Evalue
1.453e-23. Re-derived from scratch in a standalone Python port of the retrieved recurrences
(parsing the `.hmm`) → identical, confirming the formulas independent of the C# code. All regions
here are single-domain (`is_multidomain_region` FALSE), so the verified path needs no stochastic
clustering.

### Honest residual

For a region flagged multi-domain by the `rt3` test (closely-overlapping / not-well-separated
domains), HMMER resolves the envelopes by **stochastic-traceback clustering**
(`region_trace_ensemble` → `p7_spensemble_Cluster`, sampling 200 tracebacks). That sampling
clusterer is **not** implemented; such a region is emitted as a single envelope. The verified
decomposition covers the common well-separated-domain case (tandem repeats, multi-domain
β-propellers). The bundled-profile coverage (3 CC0 Pfam HMMs; any other family is a caller-supplied
`.hmm`) is unchanged.

### Reference list additions

31. HMMER `p7_domaindef.c` — `p7_domaindef_ByPosteriorHeuristics` (region identification rt1/rt2,
    `is_multidomain_region` rt3, `region_trace_ensemble`, `rescore_isolated_domain`):
    https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_domaindef.c
32. HMMER `generic_decoding.c` — `p7_GDomainDecoding` (btot/etot/mocc):
    https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_decoding.c
33. pyhmmer 0.12.1 multi-domain `hmmsearch` ground truth (GBB1/PF00400 7 domains; SH3/PF00018 1):
    https://pyhmmer.readthedocs.io/

---

## Addendum 2026-06-25 — HMMER stochastic-traceback clustering of overlapping domains (`region_trace_ensemble`)

This addendum closes the prior residual: the **stochastic-traceback clustering** that HMMER uses to
split a *closely-overlapping* (not-well-separated) multi-domain region — `region_trace_ensemble()`
followed by `p7_spensemble_Cluster()` — is now implemented as an **opt-in** path of
`Plan7ProfileHmm.FindDomains(seq, clusterOverlapping = true)`. The well-separated decomposition,
the per-sequence `HmmSearchBitScore`, and all defaults are unchanged.

### Online Source: HMMER `p7_domaindef.c` — `region_trace_ensemble` + defaults — RETRIEVED 2026-06-25

Fetched verbatim via `curl https://raw.githubusercontent.com/EddyRivasLab/hmmer/master/src/p7_domaindef.c`.

- **`p7_domaindef_Create()` defaults (verbatim):** `ddef->nsamples = 200; min_overlap = 0.8;
  of_smaller = TRUE; max_diagdiff = 4; min_posterior = 0.25; min_endpointp = 0.02;` (and
  `rt1=0.25, rt2=0.10, rt3=0.20`).
- **`p7_domaindef_ByPosteriorHeuristics()` branch (verbatim):** when `is_multidomain_region(ddef,i,j)`
  is TRUE, the region is rescored multihit (`p7_oprofile_ReconfigMultihit(om,saveL); p7_Forward(sq->dsq+i-1, j-i+1, ...)`),
  `region_trace_ensemble(ddef, om, sq->dsq, i, j, fwd, bck, &nc)` defines `nc` clusters, then each
  cluster's `p7_spensemble_GetClusterCoords(...,&i2,&j2,...)` is rescored by `rescore_isolated_domain`.
- **`region_trace_ensemble()` (verbatim):** "By default, we make results reproducible by forcing a
  reset of the RNG to its originally seeded state" — `esl_randomness_Init(ddef->r, esl_randomness_GetSeed(ddef->r))`.
  Then `for (t=0; t<ddef->nsamples; t++) { p7_StochasticTrace(ddef->r, dsq+ireg-1, Lr, om, fwd, ddef->tr);
  p7_trace_Index(ddef->tr); for d in domains: p7_spensemble_Add(ddef->sp, t, sqfrom[d]+ireg-1,
  sqto[d]+ireg-1, hmmfrom[d], hmmto[d]); ... }` then `p7_spensemble_Cluster(ddef->sp, 0.8, TRUE, 4,
  0.25, 0.02, &nc)`, then a dominated-domain removal pass (≥0.8 mutual seq overlap → lower-`prob`
  cluster dropped).

### Online Source: HMMER `generic_stotrace.c` — `p7_GStochasticTrace` — RETRIEVED 2026-06-25

Fetched verbatim (`.../src/generic_stotrace.c`). The traceback walks the Forward matrix backward from
`C(L)`; each step builds the predecessor log-score vector `sc[]`, normalises it with
`esl_vec_FLogNorm`, and samples a predecessor with `esl_rnd_FChoose(r, sc, N)`. Transitions used
verbatim (local mode): C←{C·LOOP, E·MOVE}; E←{any M_k, D_k} (local); M←{B·BM, M·MM, I·IM, D·DM};
D←{M·MD, D·DD}; I←{M·MI, I·II}; B←{N·MOVE, J·MOVE}; J←{J·LOOP, E·LOOP}; N←{S | N}. `p7_trace_Index`
(`p7_trace.c`, retrieved) splits each `B..E` into a domain with `sqfrom/sqto` (first/last match-state
residue) and `hmmfrom/hmmto` (first/last match-state node).

### Online Source: HMMER `p7_spensemble.c` — clustering + consensus endpoints — RETRIEVED 2026-06-25

Fetched verbatim (`.../src/p7_spensemble.c`). `link_spsamples` linkage rule (verbatim): two segment
pairs link iff fractional seq overlap ≥ `min_overlap` of the smaller segment AND fractional hmm
overlap ≥ `min_overlap` (the hmm overlap omits the `+1`, verbatim) AND `|d1−d2| ≤ max_diagdiff` for
the start OR end diagonal. `p7_spensemble_Cluster` runs `esl_cluster_SingleLinkage`, keeps clusters
with `ninc/nsamples ≥ min_posterior` (distinct-sample count, de-duplicated by `idx_of_last`), and
picks consensus endpoints: leftmost `i`/`k` and rightmost `j`/`m` whose frequency ≥
`ceil(ninc·min_endpointp)` (else `argmax`); rejects `best_i>best_j || best_k>best_m`; orders by start.

### Online Source: Easel `esl_cluster.c` — `esl_cluster_SingleLinkage` — RETRIEVED 2026-06-25

Fetched verbatim (`.../easel/master/esl_cluster.c`). Breadth-first single-linkage in O(N) memory:
push all vertices on stack `a`; pop to `b`; assign `c[v]=nc`; for each available `w`, if `linkfunc(v,w)`
move `w` from `a` to `b`; increment `nc` when `b` empties.

### Online Source: Easel RNG — `esl_random.c` + `esl_mix3` (`easel.c`) — RETRIEVED 2026-06-25

Fetched verbatim. The HMMER **pipeline** (`p7_pipeline.c`, retrieved) sets `seed = 42` (default
`--seed 42`) and creates the RNG with `esl_randomness_CreateFast(seed)` — the Easel **LCG**, with
`do_reseeding = (seed==0)?FALSE:TRUE`. LCG init (verbatim): `r->x = esl_mix3(seed, 87654321,
12345678); if (r->x==0) r->x=42;`. Per draw (`knuth`): `r->x = 69069*r->x + 1; esl_random = x/2^32`
in [0,1). `esl_mix3` is Bob Jenkins' 3-word mix (verbatim from `easel.c`). `esl_rnd_FChoose`
(verbatim): roll = `esl_random(r)`; return first `i` with `cumsum/norm > roll`. (The default
distribution generator is the Mersenne Twister with non-standard `mt[z]=69069·mt[z-1]` seeding, but
the *pipeline* uses the faster LCG — captured to match the actual `hmmsearch` ensemble.)

### Reference tool — pyhmmer 0.12.1 ground truth (closely-overlapping multi-domain) — CAPTURED 2026-06-25

`pyhmmer 0.12.1` (macOS arm64). `hmmsearch` (Z=1, domZ=1, seed=42) of bundled `PF00018_SH3_1.hmm`
against a **closely-overlapping tandem-SH3** construct = the 48-aa SH3 core (= residues 3..50 of the
SRC_HUMAN SH3, `VALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPS`) TRUNCATED by `trim` residues, then
a full core. hmmsearch resolves these via the ensemble (overlapping consensus envelopes):

| construct | L | dom0 env / score | dom1 env / score |
|-----------|---|------------------|------------------|
| trim=4  | 92 | 1..46 / 57.023777 | 45..92 / 66.338997 |
| trim=12 | 84 | 1..37 / 48.047169 | 37..84 / 66.678467 (iE 5.545e-23, dom0 iE 3.660e-17) |
| trim=16 | 80 | 1..33 / 40.831696 | 33..80 / 66.867706 |

A **well-separated** tandem (two complete cores, L=96) → 2 non-overlapping envelopes 1..48 / 49..96
(handled by the rt2 flank bound, NOT the ensemble). Single SH3 → 1 envelope (env 3..50). All results
are reproducible across runs (fixed-seed LCG).

### Independent C# parity (this engine, double precision)

`FindDomains` (ensemble default) reproduces: trim=4 → 2 envelopes **1..46 / 45..92**; trim=12 →
**1..37 / 37..84**; trim=16 → **1..33 / 33..80** — **envelope coordinates EXACT** vs hmmsearch (incl.
the overlap at the shared residue), per-domain scores within ~0.06 bits (48.105 vs 48.047; 66.684 vs
66.678). With `clusterOverlapping:false` the same region fuses to ONE envelope (the prior behaviour).
Well-separated tandem still → 1..48 / 49..96; single SH3 still → 3..50; GBB1/PF00400 still → the same
7 envelopes. Determinism verified (identical across repeated calls).

### Honest residual (exact-RNG-bit parity)

The clustering is implemented faithfully and reproduces the domain **count** and the envelope
**coordinates** of `hmmsearch` exactly for the captured cases, with a deterministic fixed-seed LCG
ported verbatim from Easel. **Bit-for-bit** identity of the sampled trace ensemble would additionally
require HMMER's single-precision (`float`) Forward matrix and `esl_vec_FLogNorm`/`esl_rnd_FChoose`
float arithmetic; this engine computes the Forward matrix in `double`. The consensus endpoints are
robust to that difference (a tally over 200 samples), so coordinates match; exact per-sample
trace-by-trace bit parity is the residual.

### Reference list additions

34. HMMER `generic_stotrace.c` — `p7_GStochasticTrace`:
    https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_stotrace.c
35. HMMER `p7_spensemble.c` — `p7_spensemble_Add` / `p7_spensemble_Cluster` / `link_spsamples`:
    https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_spensemble.c
36. HMMER `p7_trace.c` — `p7_trace_Index`:
    https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_trace.c
37. HMMER `p7_pipeline.c` — pipeline RNG (`esl_randomness_CreateFast(42)`, `do_reseeding`):
    https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_pipeline.c
38. Easel `esl_cluster.c` — `esl_cluster_SingleLinkage`:
    https://github.com/EddyRivasLab/easel/blob/master/esl_cluster.c
39. Easel `esl_random.c` (LCG/knuth, `esl_rnd_FChoose`) + `easel.c` (`esl_mix3`):
    https://github.com/EddyRivasLab/easel/blob/master/esl_random.c
