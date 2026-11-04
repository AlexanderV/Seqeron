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
| Kleinstiver et al. (2015) Nature 523:481, "Engineered CRISPR-Cas9 nucleases with altered PAM specificities" (PMC4540238, opened 2026-10-21, neutral User-Agent `Seqeron-review/1.0`) | VQR variant "strongly depleted sites bearing NGAN and NGCG PAMs"; the human-cell assay shows VQR "robustly cleaving sites bearing NGAN PAMs (with relative efficiencies NGAG > NGAT = NGAA > NGAC)"; the same paper's abstract names NAG and NGA as non-canonical SpCas9 PAMs. Supports `NGA` as a SpCas9 PAM motif (`SpCas9_NGA`). It does **not** contain a SaCas9 `NNGRRN` PAM (its SaCas9 site-depletion result lists only `NNGGGT`, `NNGAAT`, `NNGAGT`, all `NNGRRT`) |

### `SpCas9_NGA` and `SaCas9_NNGRRN` (added 2026-10-21, B08 round 3 G8)

**URL:** https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py (fetched again 2026-10-21 with a
neutral User-Agent; same 8956-line file). Verbatim lines:

- L188: `('NGA','20bp-NGA - Cas9 S. Pyogenes mutant VQR'),` — `NGA` is an own CRISPOR PAM entry (SpCas9 VQR mutant), 20-bp guide.
- L243-249 (`offtargetPams`): `"NGG" : ["NAG","NGA"],` … `"NGA" : ["NGG"],` `"NNGRRT" : ["NNGRRN"],` — `NNGRRN` exists in CRISPOR only
  as the alternative off-target pattern of `NNGRRT`; there is no `pamDesc` entry for it. Comment above the table: "The NGA -> NGG rule
  was described by Kleinstiver...Young 2015 'Improved Cas9 Specificity...'", "the length of the alternate PAM has to be the same as the
  original PAM".
- `setupPamInfo` L458-520: `NGA` falls in the final `else: GUIDELEN = 20`; `NNGRRT`/`NNNRRT` set `GUIDELEN = 21`. An alternative
  PAM is searched with the main system's guide length, so `NNGRRN` is 21 nt (PAM after the target — `pamIsFirst = False`).
- `patMatch` L2347-2394: pattern `N` matches any character (including a genome `N`); `R` matches `A`/`G` only; any other pattern letter
  matches only itself (so a genome `N` never satisfies `G`, `A`, `R`, `T`).

**Library mapping:** `SpCas9_NGA` = (`NGA`, 20, PAM after) and `SaCas9_NNGRRN` = (`NNGRRN`, 21, PAM after), i.e. exactly the two rows
above (`CrisprDesigner.GetSystem`, added by B08 F66 so that the alternative-PAM off-target scan can run them).

**Oracle for the new tests:** a standalone Python script (kept in the review scratchpad, not in the repository) with
(a) an independent two-strand brute-force scan — forward strand: every `i` where `seq[i:i+k]` matches the IUPAC motif and `i-guide >= 0`;
reverse strand: scan `revcomp(seq)`, map the match at `j` back to forward `Position = L - j - k`, `TargetStart = Position + k`,
PAM/target read on the protospacer strand — and (b) CRISPOR's `patMatch`, `findPat`, `findPams`, `findAllPams`, `flankSeqIter` and `revComp`
exec'd verbatim from the fetched file (`pamIsFirst = False`, `GUIDELEN` 20 / 21). The two agreed on every sequence used by the tests
(sites listed in `CrisprDesigner_PAM_Tests.cs`, region "FindPamSites - SpCas9_NGA (VQR PAM) and SaCas9_NNGRRN"): NGA 74-nt sequence 3 forward + 6 reverse
sites; NNGRRN 87-nt sequence 4 forward + 6 reverse sites; sequence-edge pairs (overrun excluded / exact fit kept) on both strands for both
systems; all 16 `NNG{A,G}{A,G}{A,C,G,T}` PAMs accepted and 6 non-R / non-G decoys rejected; `TGG` is not an NGA site; `CCGAAC` is an
`NNGRRN` but not an `NNGRRT` site; `N` in the target kept on both strands; `N` at a concrete PAM position rejected.

**Recorded divergence (not pinned by a test):** a genome `N` against a PAM `N`. CRISPOR `patMatch` accepts it; `IupacHelper.MatchesIupac('N','N')`
accepts only `A/C/G/T`, so Seqeron does not report that site. The brute-force scan in the oracle script follows Seqeron's rule and therefore
disagrees with CRISPOR on exactly those inputs (two probe sequences, `…TGN`/`CCGAAN`); no test uses them. CRISPOR's optional `doFilterNs` is not applied either.

---

## Numerical Cross-Check vs CRISPOR (2026-09-28)

The CRISPOR functions `revComp`, `patMatch`, `findPat`, `findPams`, `flankSeqIter` were copied **verbatim** from
the fetched `crispor.py` into a standalone Python module and driven directly (no re-implementation), then compared
field-by-field with `CrisprDesigner.FindPamSites`.

**Scope:** 4 hand-built sequences (SpCas9 reverse-only, SpCas9 mixed, SaCas9, Cas12a) + 60 pseudo-random sequences
(seed 20260928, lengths 10-200 nt) across all 7 systems (the then-existing ones; `SpCas9_NGA` and `SaCas9_NNGRRN` came later, see above) → **136 reference sites**.

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

- **`LbCas12a` guide length (resolved, B08 round 4 H11).** It was 24 while `Cas12a`/`AsCas12a` were 23, a split no source
  opened in the review supports. CRISPOR uses `GUIDELEN = 23` for every Cpf1 PAM (`setupPamInfo` L492-496) and does not
  model As/Lb separately; the canonical Zetsche 2015 spacer is 23 nt (IDT recommends 21). CHOPCHOP's 24
  (`CPF1_DEFAULT["GUIDE_SIZE"] = 24`) is for Cpf1 in general, again with no ortholog distinction. All three aliases now
  share the PAM and the 23-nt spacer, so one genomic site yields one protospacer, `TargetStart` and cut site whichever alias is
  chosen (`FindPamSites_AllCas12aAliases_ReportTheSameSites`). The ranges once quoted for Zetsche 2015 in tests and docs
  ("23-25", "20-24") are not in the opened source and have been removed.
