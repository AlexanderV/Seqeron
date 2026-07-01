# Test Specification: PROTMOTIF-HMM-001

**Test Unit ID:** PROTMOTIF-HMM-001
**Area:** ProteinMotif
**Algorithm:** Plan7 Profile-HMM Domain Search (HMMER3-style)
**Status:** ☑ Validated — independent Stage A/B re-validation 2026-06-25 (state CLEAN)
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Durbin, Eddy, Krogh & Mitchison (1998) *Biological Sequence Analysis* §5.4 — Plan7 Viterbi/Forward log-odds recurrences |
| 2 | Eddy (2011) *PLoS Comput Biol* 7:e1002195 — HMMER3 local N/B/M/I/D/E/J/C architecture, local entry/exit, length model, null2 |
| 3 | HMMER User's Guide v3.4 (Eddy 2023) — HMMER3/f file format, STATS LOCAL calibration lines |
| 4 | Eddy (2008) *PLoS Comput Biol* 4:e1000069 — Gumbel for Viterbi/MSV, exponential tail for Forward |
| 5 | EddyRivasLab/hmmer source: `p7_modelconfig.c`, `p7_domaindef.c`, `generic_null2.c`, `p7_bg.c`, `p7_spensemble.c`, `esl_random.c`, `esl_gumbel.c`, `esl_exponential.c`, `hmmer.c` p7_AminoFrequencies |
| 6 | Pfam PF00018 (SH3), PF00595 (PDZ), PF00400 (WD40), CC0 — bundled profiles |
| 7 | pyhmmer 0.12.1 hmmsearch ground truth (captured 2026-06-25) |

## 2. Canonical Method(s)

`Plan7ProfileHmm.{Viterbi,Forward,LocalForward}Score`, `HmmSearchBitScore`, `LocalForwardBitScore`,
`Null2BiasBits`, `FindDomains(string,bool)`, `GumbelSurvival`, `ExponentialSurvival`, `EValue`,
`{Viterbi,Msv,Forward}PValue`, `{Viterbi,Forward}EValue`; `ProteinMotifFinder.{FindDomainsByHmm,
ScoreDomainHmm, FindDomainHitsByHmm, ScoreDomainHmmEValue, FindDomainEnvelopes(×2)}`.

- **Source file:** `Plan7ProfileHmm.cs` / `ProteinMotifFinder.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`

## 3. Contract / Invariants

- R: Forward ≥ Viterbi (log-odds); strictly greater when alternative paths exist.
- R: `HmmSearchBitScore = LocalForwardBitScore − Null2BiasBits` (pipeline identity).
- R: null2 bias = `logsumexp(0, …) ≥ 0` (never raises the score).
- R: E-value `= P·Z ≥ 0`, monotone decreasing in bit score, linear in Z.
- D: deterministic given profile + sequence (ensemble reseeded to fixed pipeline seed 42 per region).
- Envelopes 1-based inclusive (HMMER env from/to). Out-of-alphabet residue → background odds 1.

## 4. Cross-check / Differential Oracle

- **Reference:** pyhmmer 0.12.1 `hmmsearch` (Z=1, domZ=1, seed=42); hand-derivation to ~1e-9.
- **Comparison:** local Forward pre-bits ±1e-4 (HMMER float32 vs double); multi-domain per-domain
  scores ≤1e-2 bit, i-Evalues ≤5%, envelope coords exact; hand DP / P-values to 1e-9 relative.
- Exact numbers recorded in `docs/Validation/reports/PROTMOTIF-HMM-001.md`.

## 5. Validation Checklist (restored ☑)

- [x] Stage A: every source opened; formulas/constants confirmed against the primary literature and
      HMMER/Easel source.
- [x] Stage B: implementation reviewed against source; pyhmmer ground truth reproduced (SH3/PDZ/WD40,
      7-domain GBB1, overlapping tandem).
- [x] Coverage gap fixed: added region H21 (MsvPValue, ExponentialSurvival, EValue positive path,
      all-X residues).
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18737 passed); 0 warnings.
- [x] `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
