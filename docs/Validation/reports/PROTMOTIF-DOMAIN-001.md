# Validation Report: PROTMOTIF-DOMAIN-001 — Protein Domain Prediction & Signal Peptide Prediction

- **Validated:** 2026-06-12   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindDomains(string)`, `ProteinMotifFinder.PredictSignalPeptide(string, int maxLength=70)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_DomainPrediction_Tests.cs` (24 tests)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

1. **PROSITE PS00028** (https://prosite.expasy.org/PS00028) — fetched. Pattern returned verbatim:
   `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`. **Exactly matches** the spec/evidence and the
   implementation regex `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H`.
2. **PROSITE PS00017** (https://prosite.expasy.org/PS00017) — fetched. Pattern returned verbatim:
   `[AG]-x(4)-G-K-[ST]`. **Exactly matches** spec/evidence and the implementation regex
   `[AG].{4}GK[ST]`.
3. **Wikipedia "Signal peptide"** — confirms tripartite structure: n-region (positively charged),
   h-region (hydrophobic core, ~5–16 residues, forms α-helix), c-region (cleaved by signal
   peptidase); overall length 16–30 aa; cites von Heijne. Matches the model in the spec.
4. **von Heijne (1983)** −1,−3 rule {A,G,S} — corroborated by evidence doc statistics
   (pos −1: A 65%, G 14%, S 10%; pos −3: A 52%, G 9%, S 8%). The implementation strictly uses
   {A,G,S}, which is the canonical rule (V/L/T occasionally appear but are not canonical).
5. **WD40 (PF00400), SH3 (PF00018), PDZ (PF00595)** — confirmed via Pfam/literature that these are
   **HMM-based Pfam families with NO authoritative PROSITE consensus pattern**. SH3/PDZ recognize
   proline-rich / C-terminal S/T-X-V partner motifs; the *domain fold* itself (~60–100 aa) is not
   captured by a short regex.

### Formula / pattern check
- C2H2 zinc finger and Walker-A/P-loop patterns are genuine PROSITE patterns reproduced **exactly**.
- Signal-peptide scoring: n-region = positive-charge density (2 K/R = 1.0), h-region = hydrophobic
  fraction {A,I,L,M,F,V,W}, c-region = small/polar fraction {A,G,S,T,N}; combined with a 1:2:1
  (n:h:c) weighting `(n + 2h + c)/4` justified by von Heijne (1985) ("hydrophobic core both
  necessary and sufficient"). This is an evidence-based heuristic, not a trained model — honestly
  declared in spec §1.4 and the Evidence "Design Decisions".

### Edge-case semantics
- Empty/null → empty (domains) / null (signal): defined.
- Short (<15 aa) signal → null: consistent with 16–30 aa minimum.
- All-charged → null (no hydrophobic core): consistent with h-region requirement.
- Case insensitivity: ToUpperInvariant + IgnoreCase: defined.

### Independent cross-check (hand-computed)
- Zinc finger `AAAACAACAAALEEEEEEEEHAAAHAAAA`: C@4, x(2)=A5A6, C@7, x(3)=A8A9A10, L@11,
  x(8)=E12..E19, H@20, x(3)=A21A22A23, H@24 → **start=4, end=24** ✓ (matches PS00028).
- P-loop `AAAAGAAEAGKSAAAA`: G@4, x(4)=A5A6E7A8, G9 K10 S11 → **start=4, end=11** ✓ (PS00017).
- Classic signal `MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF`, cleavage=25: m3=S(22), m1=G(24) ∈ {A,G,S};
  n="MKRLL" (1.0), h=15 L (1.0), c="LASAG" (0.8) → score=(1.0+2·1.0+0.8)/4 = **0.95** ✓.

### Findings / divergences (NOTES)
- **Honest-scope NOTE (not a defect):** Only PS00028 and PS00017 are authoritative PROSITE
  patterns. The WD40/SH3/PDZ "patterns" are ad-hoc simplified regexes that carry real Pfam
  accessions (PF00400/PF00018/PF00595) but are **not** the corresponding Pfam HMMs and have no
  PROSITE equivalent; they can produce false positives. This is **transparently declared** in the
  spec (§1.4.1, rows labeled "Simplified pattern"), in the tests ("Simplified pattern" comments),
  and in the Evidence doc (point 5; Design Decision 5). The code's own XML summary says "signature
  patterns", and nothing advertises HMM/Pfam-grade detection. Scope is honest → PASS-WITH-NOTES,
  not FAIL.

---

## Stage B — Implementation

### Code path reviewed
- `FindDomains` (ProteinMotifFinder.cs:986–1035) — runs 5 regexes via `FindMotifByPattern`,
  wrapping each match into a `ProteinDomain(Name, Accession, Start, End, Score, Description)`.
- `FindMotifByPattern` (:161–198) — overlapping (lookahead) matching, IgnoreCase; Start/End are
  0-based inclusive of the captured substring.
- `PredictSignalPeptide` (:383–460) — scans cleavage 15..min(35,len); enforces {A,G,S} at −1/−3,
  h-region ≥ 7, nScore > 0, hScore ≥ 0.5; keeps best 1:2:1-weighted score; Probability = Score.

### Realised correctly?
- Domain patterns are literal copies of the validated PROSITE strings (PS00028, PS00017) and the
  declared simplified patterns. Position reporting (Start ≤ End, inclusive) verified by hand above.
- Signal-peptide region slicing and scoring reproduce the spec's worked examples exactly
  (classic = 0.95, alternate = 0.825, both hand-recomputed and matched).

### Cross-verification table recomputed vs code
| Case | Expected | Code result | Match |
|------|----------|-------------|-------|
| Zinc finger M1 | start 4, end 24, PF00096 | 4 / 24 / PF00096 | ✓ |
| P-loop M2 | start 4, end 11, PF00069 | 4 / 11 / PF00069 | ✓ |
| Classic signal M7 | cleav 25, score 0.95, n="MKRLL", c="LASAG", h len 15 | identical | ✓ |
| Alternate signal S5/S6 | cleav 18, score 0.825, h len 8 | identical | ✓ |
| Threonine@−3 M15 | null (strict {A,G,S}) | null | ✓ |
| All-charged M9 | null | null | ✓ |

### Variant/delegate consistency
- `FindDomains` reuses the same `FindMotifByPattern` engine used by `FindCommonMotifs`/PROSITE
  motif search — consistent matching/score semantics across the class.

### Test quality audit
- 24 tests assert exact sourced values (positions, accessions, scores ±1e-10), null/empty edges,
  case-insensitivity, strict {A,G,S} enforcement (M15 is a real discriminating test, not a
  tautology), and multi-domain detection. No "no-throw"-only assertions.

### Findings / defects
- None. Build clean (0 warnings), all 24 domain tests pass, full suite 4486/4486 pass.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — the two real PROSITE patterns (PS00028, PS00017) match their
  authoritative sources verbatim; WD40/SH3/PDZ are honestly-declared simplified regexes (not Pfam
  HMMs). Scope is not misadvertised.
- **Stage B:** PASS — code faithfully realises the validated patterns and the von-Heijne tripartite
  heuristic; all worked examples recomputed and matched.
- **End state:** CLEAN — no defect found; no code change required.

**Sources:**
- PROSITE PS00028: https://prosite.expasy.org/PS00028
- PROSITE PS00017: https://prosite.expasy.org/PS00017
- Wikipedia "Signal peptide": https://en.wikipedia.org/wiki/Signal_peptide
- von Heijne G (1983) Eur J Biochem 133:17–21; von Heijne G (1985) J Mol Biol 184:99–105
- Krishna et al. (2003) NAR 31:532–550; Walker et al. (1982) EMBO J 1:945–951; Owji et al. (2018) Eur J Cell Biol 97:422–441
