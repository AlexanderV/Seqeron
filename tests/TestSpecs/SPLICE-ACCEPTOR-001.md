# Test Specification: SPLICE-ACCEPTOR-001

**Test Unit ID:** SPLICE-ACCEPTOR-001
**Area:** Splicing
**Algorithm:** Acceptor Site Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Shapiro & Senapathy (1987). Nucleic Acids Res 15(17):7155–7174 | 1 | — | 2026-02-12 |
| 2 | Burge, Tuschl & Sharp (1999). The RNA World, 2nd ed. CSHL Press | 1 | — | 2026-02-12 |
| 3 | Yeo & Burge (2004). J Comput Biol 11(2–3):377–394 | 1 | — | 2026-02-12 |
| 4 | Patel & Steitz (2003). Nat Rev Mol Cell Biol 4(12):960–970 | 1 | — | 2026-02-12 |
| 5 | Hall & Padgett (1994). J Mol Biol 239(3):357–365 | 1 | — | 2026-03-16 |
| 6 | Jackson (1991). Nucleic Acids Res 19(14):3795–3798 | 1 | — | 2026-03-16 |
| 7 | Dietrich, Incorvaia & Padgett (1997). Molecular Cell 1(1):151–160 | 1 | — | 2026-03-16 |
| 8 | Gao, Masuda, Matsuura & Ohno (2008). Nucleic Acids Res 36(7):2257–2267 (yUnAy) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2367711/ | 2026-06-24 |
| 9 | Mercer et al. (2015). Genome Res 25(2):290–303 (branch-point distribution) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4315302/ | 2026-06-24 |
| 10 | Wikipedia — RNA splicing | 4 | https://en.wikipedia.org/wiki/RNA_splicing | 2026-02-12 |
| 11 | Wikipedia — Polypyrimidine tract | 4 | https://en.wikipedia.org/wiki/Polypyrimidine_tract | 2026-02-12 |
| 12 | Wikipedia — Minor spliceosome | 4 | https://en.wikipedia.org/wiki/Minor_spliceosome | 2026-03-16 |
| 13 | Shapiro & Senapathy (1987) scanned PDF — Table 1 (3' ss), Methods p. 7158, Table 5 | 1 | europepmc getPdf PMC306199 | 2026-11-25 |
| 14 | Gao et al. (2008) PDF — Table 4, yeast UACUAAC | 1 | europepmc getPdf PMC2367711 | 2026-11-25 |
| 15 | Brock, Dietrich & Padgett (2008) RNA 14:2430 (U12 BPS UUCCUUAAC) | 1 | PMC2578861 (snippets) | 2026-11-25 |
| 16 | Levine & Durbin (2001) NAR 29:4006 (U12 BPS–acceptor 8–21 bp) — Levine thesis ch. 4 | 1 | ftp.sanger.ac.uk | 2026-11-25 |
| 17 | Turunen et al. (2013) WIREs RNA 4:61 (U12: no PPT, BPS–3'ss 11–13 nt) | 2 | europepmc getPdf PMC3584512 | 2026-11-25 |

### 1.2 Key Evidence Points

1. 3' splice site consensus is (Y)nNCAG|G — AG dinucleotide is almost invariant — Shapiro & Senapathy (1987)
2. PPT is 15–20 nt pyrimidine-rich region upstream of AG — Lodish et al. (2004) via Wikipedia
3. PPT quality (continuous pyrimidines) determines splice site strength — Burge et al. (1999)
4. U12-type introns end in AC (AT-AC) or AG (GT-AG), lack a PPT and are recognised through the branch site UUCCUUAAC (branch = second A) 8–21 nt upstream — Brock et al. (2008), Levine & Durbin (2001), Turunen et al. (2013). (The former `YCCAC` 3'ss consensus has no source.)
4a. **S&S 3' ss score** (Methods p. 7158): `((t1−l1)/(h1−l1)+(t2−l2)/(h2−l2))/2` over −14..−5 (best 8 of 10) and −3,−2,−1,+1; AG in −14..−5 eliminates the candidate; Table 5 reproduced 40/41 at 1 dp
5. Position -2 = A (100%), position -1 = G (100%), position -3 = C (~70%) — Shapiro & Senapathy (1987)
6. First exonic nucleotide (position 0) favors G (~50%) — Shapiro & Senapathy (1987)
7. Terminal dinucleotides alone do not distinguish U2- from U12-type introns — Dietrich et al. (1997)
8. Human branch-point consensus is **yUnAy** at positions −3..+1, branch adenosine at position 0; conservation y@−3 = 79.0%, U@−2 = 74.6%, A@0 = 92.3%, y@+1 = 75.1% — Gao et al. (2008)
9. Branch point lies 18–40 nt upstream of the 3' AG (83% at −34..−21; median −26; Mercer "19–35 nt"); the PPT spans 4–24 nt between branch point and AG — Gao et al. (2008), Mercer et al. (2015)
10. The Yeo & Burge (2004) MaxEntScan **score3ss** 23-nt maximum-entropy 3' model: `log2(P_maxent/P_background)` over a 23-nt window (20 intron + 3 exon, conserved AG at 0-based 18–19); the AG is removed (consensus/background term) and the 21-nt rest is factorised over 9 overlapping sub-sequences (5 multiplied, 4 divided) via precomputed probability tables. Documented worked examples: `score3('ttccaaacgaacttttgtAGgga')=2.89`, `score3('tgtctttttctgtgtggcAGtgg')=8.19`, `score3('ttctctcttcagacttatAGcaa')=-0.08` — Yeo & Burge (2004); factorisation + tables from the MIT-licensed maxentpy port (now bundled, opt-in `ScoreAcceptorMaxEnt`)

### 1.3 Documented Corner Cases

1. **U12 acceptor**: AC or AG terminus with a UUCCUUAAC branch site 8–21 nt upstream — Brock et al. (2008), Levine & Durbin (2001)
2. **Weak PPT**: Interrupted PPT reduces score; may be skipped in alternative splicing — Burge et al. (1999)
3. **Cryptic AG sites**: AG dinucleotides in intronic sequence can mimic acceptor sites; context and PPT quality distinguish real from cryptic — Shapiro & Senapathy (1987)

### 1.4 Known Failure Modes / Pitfalls

1. **Insufficient context**: a candidate whose −14..+1 (MaxEnt −20..+3) window runs off the sequence or contains a non-ACGTU character is not scored
2. **T/U ambiguity**: DNA (T) and RNA (U) representations must be handled equivalently — implementation converts T→U

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindAcceptorSites(sequence, minScore, includeNonCanonical)` | SpliceSitePredictor | Canonical | S&S 1987 primate |
| `FindAcceptorSites(sequence, minScore, includeNonCanonical, AcceptorScoringModel, ShapiroSenapathyCategory)` | SpliceSitePredictor | Canonical | S&S subtable / MaxEntScan score3 |
| `FindBranchPoints(sequence, searchStart, searchEnd, minScore[, BranchPointModel])` | SpliceSitePredictor | Canonical | yUnAy (default) / UACUAAC |
| `FindAcceptorBranchPoint(sequence, acceptorAgPosition, minScore)` | SpliceSitePredictor | Canonical | best yUnAy 18–40 nt upstream + PPT fraction |
| `ScoreAcceptorMaxEnt(window)` | SpliceSitePredictor | Canonical | MaxEntScan score3ss |
| `ScoreAcceptorSiteShapiroSenapathy`, `ScoreAcceptorSiteMaxEnt`, `ScoreU12AcceptorSite`, `ScoreBranchPointConsensus` | SpliceSitePredictor | Internal | tested via the public scanners |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | S&S / U12 Score ∈ [0, 1] (MaxEnt in bits) | Yes | (t−l)/(h−l) halves; matches/9 |
| INV-2 | Confidence ∈ [0, 1] | Yes | clamp / logistic |
| INV-3 | AG sites have Type = Acceptor | Yes | |
| INV-4 | U12 candidates (AC/AG) have Type = U12Acceptor | Yes | |
| INV-5 | Position = i + 1 (last intron nucleotide) | Yes | |
| INV-6 | Motif = −20..+3 context (truncated at ends) | Yes | |
| INV-7 | Empty/null input → empty result; incomplete window → no site | Yes | full-window rule |
| INV-8 | Higher minScore → subset | Yes | `score >= minScore` |
| INV-9 | An AG inside −14..−5 eliminates the candidate | Yes | S&S 1987 Methods |
| INV-10 | Branch-point score ∈ [0, 1]; perfect yUnAy ⇒ 1.0 | Yes | Gao (2008) Table 4 weights |
| INV-11 | `FindAcceptorBranchPoint` window 18–40 nt | Yes | Gao (2008), Mercer (2015) |
| INV-12 | U12 branch site counted only 8–21 nt upstream | Yes | Levine & Durbin (2001) |
| INV-13 | `ScoreAcceptorMaxEnt` = maxentpy score3 (2.89 / 8.19 / −0.08) | Yes | Yeo & Burge (2004) |
| INV-14 | DNA ≡ RNA, case-insensitive | Yes | |
| INV-15 | `ScoreAcceptorMaxEnt` rejects windows ≠ 23 nt / non-ACGTU | Yes | |

---

## 4. Test Cases (file `Unit/Annotation/SpliceSitePredictor_AcceptorSite_Tests.cs`)

| ID | Test Case | Expected | Evidence |
|----|-----------|----------|----------|
| M1 | `UUUUUUUUUUUUUUUUCAGGG` | 1 site, pos 18, score 1.0, motif = input | S&S Table 1 primate |
| M1-T5 | 41 published Table 5 3' ss (6 categories) | 40 identical at 1 dp + rounding outlier exact fraction | S&S Table 5 |
| M1-AG | AG in −14..−5 eliminated; AA instead → 0.978125 | | S&S Methods |
| M1-2P | two purines at the two lowest U rows → 1.0 | | S&S Methods ("two purines can score as well") |
| M1-ME | MaxEnt mode = score3 (2.8867730651152104), partial window skipped, bits threshold | | Yeo & Burge 2004 |
| M2/M3 | no AG / empty / null → empty | | |
| M4 | 14 nt → empty; 15 nt → site; no +1 → empty | | full window |
| M5 | all-U 1.0 > all-A 0.5171875 | | S&S |
| M6–M9 | ranges, DNA≡RNA, case | | |
| M10 | two AG 1.0/1.0; shadowed AG eliminated, first 0.9423076923… | | S&S |
| S1 | U12 branch site at distance 8 → 1.0 at pos 16; window 8/21 in, 7/22 out; GT-AG U12 reported as both; 8/9 partial | | Brock 2008, Levine & Durbin 2001 |
| S2–S5, C1 | defaults, position, motif, threshold, < 12 upstream nt | | |
| C2 | independent best-8-of-10 recomputation (363, 57, 377) | 0.978125 | S&S |
| C3 | N inside window → no site; N outside → 1.0 | | full window |
| BP1–BP8 | `FindAcceptorBranchPoint` (unchanged) + N placement skipped | | Gao 2008 |
| FBP | `FindBranchPoints` yUnAy (CUAAC 1.0, UUCAG 2.459/3.210), yeast UACUAAC 7/7, 6/7, N skipped, unknown model throws | | Gao 2008, Langford & Gallwitz 1983 |
| ME1–ME9 | `ScoreAcceptorMaxEnt` references | | maxentpy |

Mutation/fuzz/metamorphic/snapshot tests re-based on the S&S formula (independent Table 1 oracle), U12 branch-site
window bounds, yUnAy oracle; snapshots oracle-checked (acceptor 0.5362179487179487 / 1.0; gene structure 0.8767814236229544).

---

## 5. Audit of Existing Tests

Rewritten 2026-11-25 (B11 F27–F31): tests that locked the unsourced PWM/PPT score (0.8393), the `YCCAC` U12 weights
(x/3.5), the YNYURAC branch PWM and the `< 20 nt` guard were re-based on sourced values (S&S Table 1/5, Brock/Levine,
Gao Table 4) confirmed by the independent Python oracle (Evidence Dataset 6).

---

## 6. Assumption Register

**Total assumptions:** 0

All previously documented assumptions have been resolved:

| # | Former Assumption | Resolution | Evidence |
|---|-------------------|------------|----------|
| 1 | PWM weights approximate | **Superseded (B11 F27)** — the 8-offset PWM had no table source; replaced by the S&S (1987) Table 1 3' ss score | S&S (1987) Table 1, Methods |
| 2 | Normalization formula implementation-specific | **Superseded (B11 F27)** — published S&S formula, score/100 | S&S (1987) Methods |
| 3 | U12 acceptor scoring | **Superseded (B11 F29)** — the YCCAC (1, 1, 0.5) + PPT /3.5 weights had no source; replaced by the UUCCUUAAC branch-site match fraction 8–21 nt upstream | Brock (2008), Levine & Durbin (2001), Turunen (2013) |
| 4 | Branch-point detection (consensus + window + weights) | **Resolved** — `FindAcceptorBranchPoint` uses the yUnAy consensus, 18–40 nt window, and per-position conservation weights all from Gao et al. (2008), corroborated by Mercer et al. (2015). No invented constants | Gao et al. (2008), Mercer et al. (2015) |
| 5 | MaxEntScan score3ss model (tables + factorisation) | **Resolved** — `ScoreAcceptorMaxEnt` embeds the precomputed probability tables and applies the published factorisation from the MIT-licensed maxentpy port; reproduces the documented 2.89/8.19/-0.08 worked examples exactly. No invented or guessed values | Yeo & Burge (2004); maxentpy (MIT) |

---

## 7. Bug Fixes

| Date | Bug | Fix | Impact |
|------|-----|-----|--------|
| 2026-03-16 | AcceptorPwm alignment off-by-2 in `ScoreAcceptorSite` — PWM offsets are splice-site-relative but were applied to A-position | Changed `int pos = position + offset` → `int pos = position + 2 + offset` | Canonical site with strong PPT: score 0.27 → 0.84. Gene structure prediction now detects proper introns. |

---

## 8. Open Questions / Decisions

**MaxEntScan score3ss — now bundled (2026-06-24).** The Yeo & Burge (2004) MaxEntScan
23-nt maximum-entropy 3' acceptor model is implemented as the opt-in `ScoreAcceptorMaxEnt`.
The precomputed probability tables (`Data/maxent_score3.txt`, 82 560 records) and the
factorisation were retrieved this session from the **MIT-licensed maxentpy port**
(redistribution permitted; provenance + licence in `Data/maxent_score3.LICENSE.md`) and
embedded as a resource. The implementation reproduces the documented worked examples
exactly (2.89 / 8.19 / -0.08), so a wrong table/factorisation fails the 2.89 cross-check.
The default `FindAcceptorSites` PWM+PPT scorer is unchanged.

**Residual (honest, stated):** only the **3' acceptor** `score3ss` model is bundled; the
corresponding Yeo & Burge (2004) **5' donor** `score5ss` maximum-entropy model is not yet
included (it needs the analogous `me2x5` / `score5_matrix` tables — a separate import). That
narrowed residual is recorded in `docs/Validation/LIMITATIONS.md`.

**B11 review (2026-11-25):** default = Shapiro & Senapathy (1987) 3' ss score (F27); full-window scoring (F28);
U12 acceptors by branch site (F29); `FindBranchPoints` yUnAy / UACUAAC (F30); MaxEntScan scan mode + MCP
`model`/`category` (F31). The 5' MaxEnt residual above is resolved (`ScoreDonorMaxEnt`, SPLICE-MAXENT5-001).
