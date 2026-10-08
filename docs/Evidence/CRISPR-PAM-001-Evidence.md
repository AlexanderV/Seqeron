# Evidence Artifact: CRISPR-PAM-001

**Test Unit ID:** CRISPR-PAM-001
**Algorithm:** PAM Site Detection (`CrisprDesigner.FindPamSites`, `CrisprDesigner.GetSystem`, `PamSite` / `CrisprSystem` records)
**Date Collected:** 2026-09-28 (campaign `docs/Validation/review-2026-09`, batch B08)

---

## Online Sources

### CRISPOR — reference implementation source (`crispor.py`, master)

**URL:** https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py
**Accessed:** 2026-09-28 (fetched with `curl`, 8956 lines; functions read verbatim)
**Authority rank:** 3 (reference implementation — CRISPOR is the standard public guide-design server, Concordet & Haeussler, NAR 2018)

**Key Extracted Points (line numbers in the fetched file):**

1. **PAM/guide-length table** (`pamDesc`, L173-200 and `setupPamInfo`, L458-520):
   - `NGG` → `20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1`; default `GUIDELEN = 20`.
   - `NNGRRT` → `21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus`; `GUIDELEN = 21`, `saCas9Mode = True`.
   - `TTTV` → `TTT(A/C/G)-23bp - Cas12a (Cpf1) - recommended, 23bp guides`; `pamIsFirst = True`, `GUIDELEN = 23`
     (CRISPOR does **not** distinguish AsCas12a from LbCas12a; `TTTV-21` is offered as the IDT variant).
   - `pamIsCasX(pam)` (L897-899) = `pam in ["TTCN"]` → `GUIDELEN = 20`, `pamIsFirst = True`.
2. **Both strands are scanned in forward coordinates** (`findAllPams`, L4627-4641): the forward sequence is
   searched for the PAM motif (strand `+`) and for `revComp(pam)` (strand `-`); the key of `startDict` is the
   **forward-strand start of the PAM match** on both strands.
3. **Boundary rule** (`findPams`, L1033-1093): a match is kept only when a full `GUIDELEN` guide fits —
   `minPosPlus = GUIDELEN` / `maxPosMinus = len(seq)-(GUIDELEN+len(pam))` for PAM-after systems, mirrored for
   `pamIsFirst` (Cpf1/CasX) systems.
4. **Coordinate & orientation convention** (`flankSeqIter`, L1582-1626) — the authoritative definition used by
   this unit:
   - `guideStart` is a **forward-strand** coordinate on both strands:
     PAM-after (`Cas9`): `+` → `pamStart-GUIDELEN`, `-` → `pamStart+pamLen`;
     PAM-first (`Cpf1`/`CasX`): `+` → `pamStart+pamLen`, `-` → `pamStart-GUIDELEN`.
   - `pamSeq` and `flankSeq` (the guide) are read **on the protospacer strand**: for `-` strand hits both are
     reverse-complemented (`pamSeq = revComp(seq[pamStart:pamStart+pamLen])`), so a reverse-strand SpCas9 PAM is
     reported as `TGG`, never as the forward-strand `CCA`.
5. **IUPAC semantics** (`patMatch`, L2347-2394): `N` = any, `R` ∈ {A,G}, `V` ∈ {A,C,G} — identical to
   `IupacHelper.MatchesIupac`.
6. **Cas12a cut site** (L1565, CRISPOR UI text citing Zetsche et al. 2015 Fig. 3): "cleavage occurs usually … after
   the 18th base on the non-targeted strand which has the TTTV PAM motif … Cleavage mostly occurs after the 23rd
   base on the targeted strand".

### Primary literature (PAM definitions, carried forward and re-confirmed against CRISPOR's tables)

| Source | Confirms |
|--------|----------|
| Jinek et al. (2012) Science 337:816 | SpCas9 `NGG` PAM 3' of a 20-nt protospacer; blunt cut 3 bp 5' of the PAM |
| Hsu et al. (2013) Nat Biotechnol 31:827 | `NAG` as the secondary, lower-activity SpCas9 PAM |
| Ran et al. (2015) Nature 520:186 | SaCas9 `NNGRRT` PAM, 21-nt spacer (21-23 functional) |
| Zetsche et al. (2015) Cell 163:759 | Cas12a/Cpf1 `TTTV` PAM 5' of the spacer, canonical 23-nt spacer; staggered cut |
| Liu et al. (2019) Nature 566:218 | CasX/Cas12e `TTCN` PAM 5' of a 20-nt spacer |

---

## Numerical Cross-Check vs CRISPOR (2026-09-28)

The CRISPOR functions `revComp`, `patMatch`, `findPat`, `findPams`, `flankSeqIter` were copied **verbatim** from
the fetched `crispor.py` into a standalone Python module and driven directly (no re-implementation), then compared
field-by-field with `CrisprDesigner.FindPamSites`.

**Scope:** 4 hand-built sequences (SpCas9 reverse-only, SpCas9 mixed, SaCas9, Cas12a) + 60 pseudo-random sequences
(seed 20260928, lengths 10-200 nt) across all 7 systems → **136 reference sites**.

Mapping: CRISPOR `pamStart` ↔ `PamSite.Position`, `guideStart` ↔ `PamSite.TargetStart`, `pamSeq` ↔
`PamSite.PamSequence`, `flankSeq` ↔ `PamSite.TargetSequence`, strand ↔ `IsForwardStrand`.

**Result before the fix:** every forward-strand site agreed exactly; **every reverse-strand site disagreed** in two
fields, e.g. for `CCAACGTACGTACGTACGTACGTACGTACGT` (SpCas9):

| Field | CRISPOR | Seqeron (before) | Seqeron (after) |
|-------|---------|------------------|-----------------|
| `pamStart` / `Position` | 0 | 0 | 0 |
| `guideStart` / `TargetStart` | 3 | **8** (index into the reverse-complement string) | 3 |
| `pamSeq` / `PamSequence` | `TGG` | **`CCA`** (reverse-complemented back to forward) | `TGG` |
| `flankSeq` / `TargetSequence` | `ACGTACGTACGTACGTACGT` | `ACGTACGTACGTACGTACGT` | same |

Further diverging reference rows used as locked test values (`ACGTACGTACGTACGTACGTAGGCCATTGACCAGTTGCAGTCGGATCC
GTAAGCTTGGCCTAGCTAGCTAGGTACCGGATCCAAGCTT`, SpCas9, reverse strand):

| `pamStart` | `guideStart` | `pamSeq` | `flankSeq` |
|---|---|---|---|
| 23 | 26 | TGG | ATCCGACTGCAACTGGTCAA |
| 30 | 33 | TGG | CTTACGGATCCGACTGCAAC |
| 46 | 49 | CGG | AGCTAGCTAGGCCAAGCTTA |
| 58 | 61 | AGG | GATCCGGTACCTAGCTAGCT |

and for Cas12a (`TTTAACGTACGTACGTACGTACGTACGTACGTACGTAAACCATTGACCAGTTGCAGTCGGATCCGTAAGCTTGG`):
`+` strand `pamStart=0, guideStart=4, pamSeq=TTTA`; `-` strand `pamStart=35, guideStart=12, pamSeq=TTTA,
flankSeq=CGTACGTACGTACGTACGTACGT`.

**Result after the fix:** 0 mismatches over all 136 sites (positions, guide starts, PAM strings, guide strings,
strand and site counts identical to CRISPOR).

---

## Divergences Recorded (not defects)

- **`LbCas12a` guide length 24 vs CRISPOR's 23.** CRISPOR uses `GUIDELEN = 23` for every Cpf1 PAM and does not
  model As/Lb separately; the canonical Zetsche 2015 spacer is 23 nt (20-24 functional, IDT recommends 21).
  Seqeron's 24 lies inside the published functional range but is not supported by any source opened in this
  session — flagged for the ledger rather than changed (campaign rule 1: no unsourced changes).
