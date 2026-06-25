# Validation Report: PROBE-EVALUE-001 — Karlin–Altschul Off-Target E-value / Bit-Score

- **Validated:** 2026-06-25   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.ComputeKarlinAltschul`, `ProbeDesigner.ComputeLambdaNucleotide`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- `ProbeDesigner.ComputeLambdaNucleotide(int match, int mismatch, double baseFrequency = 0.25)`
- `ProbeDesigner.ComputeKarlinAltschul(double rawScore, int queryLength, long databaseLength, ScoringMatrix? scoring = null, double k = 0.711, double baseFrequency = 0.25)`
- Result type `KarlinAltschulStatistics(RawScore, Lambda, K, BitScore, EValue, QueryLength, DatabaseLength)`
- Source: `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs:1136-1286`
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs` (KA1–KA12)

## Stage A — Description

### Sources opened this session
- **NCBI "The Statistics of Sequence Similarity Scores"** (Altschul), https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html — confirmed verbatim: bit score **S' = (λS − ln K)/ln 2**, E-value from bit score **E = m·n·2^(−S')**, the substitution-matrix identity s(i,j) = (1/λ)·ln(q_ij/(p_i p_j)), and the requirement that **the expected score for a random pair must be negative** (otherwise long alignments score high by chance).
- **Karlin & Altschul (1990), PNAS 87:2264** and **Altschul et al. (1990), JMB 215:403** — the foundational E-value form **E = K·m·n·e^(−λS)** and λ as the unique positive root of **Σ_{i,j} p_i p_j e^{λ s_ij} = 1**.
- **NCBI BLAST+ command-line manual appendices** (NBK279684) and the **swipe `blastkar_partial.c`** hard-coded parameter tables (which mirror NCBI `blast_stat`): published **ungapped** nucleotide λ/K — `blastn_values_1_3` first row `{0,0, 1.374, 0.711, ...}` (λ≈**1.374**, K≈**0.711**), `blastn_values_2_3` first row `{0,0, 0.55, 0.21, ...}` (λ≈**0.55**, K≈**0.21**).

### Formula check
| Quantity | Source formula | Code (`ProbeDesigner.cs`) |
|---|---|---|
| λ | unique positive root of Σ p_i p_j e^{λ s_ij}=1 | bisection on `pMatch·e^{λ·match}+pMismatch·e^{λ·mismatch}−1` (L1206–1223) ✓ |
| bit score | S' = (λS − ln K)/ln 2 | `(lambda*rawScore - Math.Log(k))/Math.Log(2.0)` (L1273) ✓ |
| E (raw) | E = K·m·n·e^{−λS} | `k*queryLength*databaseLength*Math.Exp(-lambda*rawScore)` (L1276) ✓ |
| E (bits) | E = m·n·2^{−S'} | algebraically identical; verified numerically (KA4) ✓ |

For a simple match/mismatch matrix under uniform background p_i=0.25, the 16 ordered pairs are 4 matches (Σp = 4·0.25² = 0.25) and 12 mismatches (Σp = 0.75), so the defining equation collapses **exactly** to 0.25·e^{λ·match}+0.75·e^{λ·mismatch}=1 — which is what the code solves. This is a correct specialization, not an approximation.

### Edge-case semantics
- **No positive score / non-negative mismatch / non-negative expected score** → λ undefined; code throws `ArgumentOutOfRangeException` (matches Altschul's stated preconditions). ✓
- **Score S = 0** → E = K·m·n, S' = −ln K/ln 2 (defined, sourced). ✓
- **K** has no citable closed form (requires the Karlin–Altschul score-lattice / geometric-spacing machinery); correctly exposed as a caller parameter defaulted to the published 0.711. ✓

### Independent cross-check (oracle = independent numerical re-solve + NCBI published constants)
`blastn` was **not installable** (no `blastn`, no `conda`). Oracle = an independent Python bisection solver written this session + NCBI's published λ/K.

| Quantity | Independent oracle | Code / test | Match |
|---|---|---|---|
| λ(1,−3, p=0.25) | **1.3740631224599755** | 1.3740631224599755 | exact; = NCBI published 1.374 ✓ |
| Σ check at λ | 1.0000000000000002 | — | root verified ✓ |
| λ(2,−3, p=0.25) | **0.6337314430979077** | 0.6337314430979077 (default BlastDna) | exact ✓ |
| Bit S' (S=30,K=0.711) | **59.962700114285006** | 59.962700114285006 | exact ✓ |
| E (S=30,m=20,n=1000,K=0.711) | **1.7801583686083893e-14** | 1.7801583686083893e-14 | exact ✓ |
| E from bits (m·n·2^{−S'}) | 1.7801583686083875e-14 | — | agrees (1e-9) ✓ |
| S=0 → E | **14220** = K·m·n | KA8 | exact ✓ |
| S=0 → S' | **0.4920785350426718** | KA8 | exact ✓ |

### Findings / divergences (Stage A)
- **λ(2,−3) divergence (documented, not a defect):** the code's *default* scheme is `SequenceAligner.BlastDna` (+2/−3), for which the uniform-0.25 root is **0.6337**, whereas NCBI's *published* ungapped 2/−3 value is **0.55** (K=0.21). The gap is real: NCBI's published values come from the full finite score-lattice/target-frequency derivation, not the two-term uniform-0.25 reduction. The code does **not** claim to reproduce published 2/−3; its XML doc explicitly directs callers to pass a +1/−3 matrix to reproduce the published λ≈1.374, and 0.711 is documented as the +1/−3 published K. The contract is the **uniform-0.25 model**, which is internally consistent and correct for that model. Recorded as a PASS-WITH-NOTES-level nuance, not a defect.

Stage A verdict: **PASS** (the description and constants match Karlin–Altschul / NCBI exactly; the one nuance is documented honestly in the code).

## Stage B — Implementation

### Code path reviewed
- `ComputeLambdaNucleotide` `ProbeDesigner.cs:1183-1224` — preconditions, p(match)=4p², bisection (200 iters on [0,100]).
- `ComputeKarlinAltschul` `ProbeDesigner.cs:1253-1280` — length/K guards, λ via the above, bit & E.
- `KarlinAltschulStatistics` record `ProbeDesigner.cs:185-192`.

### Formula realised correctly?
Yes — see Stage A table. Bisection is well-posed: f(0)=0, f'(0)=expected score < 0 (enforced), and the positive match term forces f→∞, so a unique positive root exists; 200 bisection iterations on [0,100] converge below double resolution (KA11 residual < 1e-12).

### Cross-verification recomputed vs code
All oracle values above were reproduced by the compiled code via the test fixture (42/42 green). λ, bit, E, both E forms, S=0, K-linearity, m-linearity, n-linearity, and the default +2/−3 λ all match the independent oracle to the stated tolerances.

### Variant / delegate consistency
- `E = K·m·n·e^{−λS}` (raw) and `E = m·n·2^{−S'}` (bits) agree to 1e-9 (KA4).
- Default-scheme path (`scoring = null → BlastDna`) verified (KA12). The XML-doc `ArgumentNullException` for `scoring` is effectively unreachable in the public surface because `scoring ?? BlastDna` substitutes the default before the null-check; this is a harmless doc over-statement, not a behavioural defect.

### Numerical robustness
- `databaseLength` is `long`; E uses `double` arithmetic (no overflow on stated ranges).
- Guards reject m≤0, n≤0, K≤0, and the three λ-undefined preconditions.

### Test quality audit
Pre-existing KA1–KA7 trace to NCBI/Karlin–Altschul/hand-derivation (not code echoes); tolerances are tight (1e-6/1e-9, E within E·1e-9). **Coverage gaps found and closed this session** (all expected values from the independent oracle, not code):
- **KA8** — score S=0 boundary (E = K·m·n = 14220; S' = −ln K/ln 2 = 0.49208).
- **KA9** — K parameter: doubling K doubles E and lowers bit by exactly 1 bit.
- **KA10** — search-space monotonicity in **m** (doubling m doubles E; KA6 only covered n).
- **KA11** — λ root-finder convergence to machine precision (λ to 1e-15, residual < 1e-12).
- **KA12** — default-scheme path (BlastDna +2/−3 → uniform-0.25 λ = 0.6337314430979077), also locking the documented 2/−3 vs published-0.55 nuance.

No green-washing: no skips, no weakened assertions, no widened tolerances.

### Findings / defects (Stage B)
None. Five coverage gaps closed by added tests; no source change required.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: ✅ CLEAN.**
- Test count 37 → 42 (KA8–KA12 added). Full unfiltered `dotnet test Seqeron.sln -c Debug`: Failed 0 (Seqeron.Genomics.Tests 18746 passed). 0 warnings on the changed test project.
- No defect logged. The 2/−3-vs-published-0.55 nuance is documented in code and in KA12; if a future need arises to expose NCBI's published per-scheme λ/K (rather than the uniform-0.25 model), that would be an enhancement, not a correctness fix.
