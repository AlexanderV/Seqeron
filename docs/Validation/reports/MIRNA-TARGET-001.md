# Validation Report: MIRNA-TARGET-001 — microRNA Target Site Prediction (incl. TA_3UTR)

- **Validated:** 2026-06-26 (fresh re-validation after limitation-fix 5f2fbd40; supersedes the 2026-06-25 entry)   **Area:** MiRNA
- **Scope:** the OWN canonical surface of MIRNA-TARGET-001 — seed-match site detection (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer), the base/default scorer, the context++ wiring — PLUS the newly-added TA_3UTR capability `ComputeTa3Utr` / `CountSeedSites3Utr` (commit 5f2fbd40). This session FOCUSES on the new TA surface and confirms the pre-existing surface still holds.
- **Canonical method(s):** `MiRnaAnalyzer.FindTargetSites` (site detection); `ComputeTa3Utr(MiRna, IEnumerable<string>)` and `CountSeedSites3Utr(...)` (new — TA_3UTR feature); supporting `GetSeedSequence`, `GetReverseComplement`, private `CountSeedSitesInUtr`, and the context++ wiring `ScoreTargetSiteContextPlusPlus` / `TaContribution`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs` (M-/S-/E-/CTX-/CTX-TA- regions).
- **Stage A verdict:** PASS (✅)
- **Stage B verdict:** PASS (✅) — one divergence from the literal source ("non-overlapping") found and fully fixed this session.
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened THIS session (external, first-source)
- **Garcia D, Baek D, Shin C, Bell GW, Grimson A, Bartel DP (2011) *Nat Struct Mol Biol* 18:1139–1146**, "Weak seed-pairing stability and high target-site abundance decrease the proficiency of lsy-6 and other microRNAs." Open-access mirror **PMC3190056** (Online Methods) quoted **verbatim** this session:
  - **TA definition:** *"TA in the human transcriptome was calculated as the number of **non-overlapping** 3′UTR **8mer, 7mer-m8, and 7mer-A1** sites in the reference mRNAs."* → site set = exactly the three high-confidence canonical types (bare 6mer / offset-6mer EXCLUDED); counting is **non-overlapping**.
  - **log10 scale:** *"each additional CG dinucleotide imparted an additional **log10** reduction in TA."* → TA is carried on a log10 scale.
- **Agarwal V, Bell GW, Nam JW, Bartel DP (2015) *eLife* 4:e05005**, "Predicting effective microRNA target sites in mammalian mRNAs." PMC4532895 + Bartel-lab reprint. **Table 1**: the context++ feature **`TA_3UTR` = "Number of sites in all annotated 3′ UTRs"**, attributed to *Arvey et al., 2010; Garcia et al., 2011*. Feature selected in 100% of site-type models; its contribution is `coeff × min-max-scaled(value)`.
- **TargetScan 7.0 reference implementation** (`targetscan_70_context_scores.pl`, nsoranzo mirror) — the decisive artifact for the log10 convention and the scaling form:
  - `read_TA_SPS`: `($seedRegion, $SPS_1, $SPS_2, $TA) = split(/\t/, $_); $garcia{$seedRegion}{"TA"} = $TA;` — TA is **column 4** of `TA_SPS_by_seed_region.txt`, stored and used **as-is**.
  - `getAgarwalContribution`: `scaled = (raw - min)/(max - min); contribution = coeff × scaled` — exactly the min-max form.
  - `$TA_contribution = getAgarwalContribution($siteType, "TA_3UTR", $garcia{$seedRegion}{"TA"});` — the file value is fed **without any log transform in the script**, so the stored value IS `log10(count)` (Garcia computes the log10; the script consumes it). Min/max ≈ 3.07–3.89 ⇒ raw counts ≈ 10^3.1…10^3.9 ≈ 1.3k…7.9k sites genome-wide — biologically sensible.
- **Grimson et al. (2007)** (PMC3800283, opened) — confirms the canonical site-type definitions used to define which sites TA counts: **8mer** = seed match flanked by the m8 match AND the A1; **7mer-m8** = seed match + m8; **7mer-A1** = seed match + A at target pos 1; **6mer** = perfect match to miRNA nt 2–7. (Bartel 2009 / Lewis 2005 reconfirm seed = nt 2–7 and reverse-complement targeting; re-litigated in the prior report, unchanged.)

### Formula validated
`TA_3UTR = log10(N)`, where **N = total number of non-overlapping 8mer + 7mer-m8 + 7mer-A1 sites of the miRNA seed across the supplied 3′UTR set**. This matches the Garcia (2011) Online-Methods definition exactly (site-type membership, non-overlapping, log10 scale) and feeds the Agarwal (2015) TA coefficient via the published min-max scaling.

### Edge-case semantics (defined & sourced)
- **Empty UTR set / only empty-or-null UTRs / no sites** ⇒ N = 0 ⇒ TA = 0 (log10(1) floor; log10(0) is undefined and TargetScan never emits a seed with zero sites — the 0 floor is the safe convention).
- **DNA UTRs** (T) normalised to RNA (U) before counting.
- **Seed < 7 nt / degenerate** ⇒ no defined site set ⇒ 0.
- **Non-ACGU** never matches the core (the RC of any ambiguous base is `N`, which still participates only as a literal char; spurious matches impossible for the canonical ACGU core).
- **Overlapping sites** — counted **non-overlapping** per Garcia (see Stage-B finding/fix).
- **null enumerable** ⇒ `ArgumentNullException` (contract).

### Independent cross-check (hand-derived, reproduced by an independent Python reimplementation written from the Garcia definition, NOT from the C# code)
- **CTX-TA-001** — synthetic seed `ACGUACG` (pos 2-8), seedRC `CGUACGU`, pos8Rc `C`, core `GUACGU`. UTRs `{CGUACGUA, GGUACGUG, CGUACGUGAAAGUACGUC}`:
  - `CGUACGUA`: core@1, upstream `C`=m8, downstream `A`=A1 → **8mer** (1).
  - `GGUACGUG`: core@1, upstream `G` (no m8), downstream `G` (no A1) → bare 6mer → **0**.
  - `CGUACGUGAAAGUACGUC`: core@1 m8-only → 7mer-m8 (1); the second core@11 is bare → 0 → **1**.
  - **N = 2 ⇒ TA = log10(2) = 0.301029995663981.**
- **CTX-TA-002** — let-7a seed `GAGGUAG`, seedRC `CUACCUC`, pos8Rc `C`, core `UACCUC`. 8mer + 7mer-m8 + 7mer-A1 + (bare 6mer excluded) + (8mer + 7mer-m8) = **N = 5 ⇒ TA = log10(5) = 0.698970004336019.**
- **TA→context++:** for an 8mer site, TA enters as `0.222 × ((log10(5) − 3.113)/(3.865 − 3.113))` (Agarwal 8mer TA row); for a 7mer-m8 site, `0.139 × ((log10(2) − 3.067)/(3.887 − 3.067))` — both reproduced exactly by `ScoreTargetSiteContextPlusPlus`/`TaContribution`.

**Stage A finding:** No biological or mathematical error. The formula, the site-type membership, the log10 scale and the scaling-into-context++ are all confirmed against Garcia (2011), Agarwal (2015) and the TargetScan reference implementation. **PASS.** (The original code comment's claim that "non-overlapping" is automatically satisfied was the one inaccuracy — see Stage B.)

---

## Stage B — Implementation

### Code path reviewed
`ComputeTa3Utr` (MiRnaAnalyzer.cs:~888), `CountSeedSites3Utr` (~906), private `CountSeedSitesInUtr` (~938), `GetSeedSequence` (95), `GetReverseComplement` (296), `TaContribution` (1482) + the `ScoreTargetSiteContextPlusPlus` wiring (806/835); pre-existing `FindTargetSites` (156) reconfirmed unchanged by the limitation-fix.

### Realises the validated description? (evidence)
- **TA = log10(N), N over the three site types** — `CountSeedSites3Utr` derives the seed (SeedSequence, else `GetSeedSequence` = positions 2-8), takes `seedRC`, `pos8Rc = seedRC[0]`, `sixmerCore = seedRC[1..7]`, sums `CountSeedSitesInUtr` over each UTR (T→U, upper); `ComputeTa3Utr` returns `Math.Log10(total)` (0 when total = 0). A site counts iff the 6mer core matches AND (m8 upstream OR A1 downstream) — i.e. 8mer/7mer-m8/7mer-A1, bare 6mer excluded. Matches Garcia's membership exactly. ✓
- **log10 / floor** — `total > 0 ? Math.Log10(total) : 0.0` matches the TargetScan convention. ✓
- **Feeds context++ with the bundled Agarwal TA coefficient + min-max scaling** — `TaContribution` selects per-site-type `(coeff,min,max)` (8mer 0.222 / 3.113–3.865; 7mer-m8 0.139 / 3.067–3.887; 7mer-A1 0.117 / 3.145–3.887; 6mer 0.058 / 3.113–3.887) and applies `coeff × (ta−min)/(max−min)` — exactly the perl `getAgarwalContribution`. ✓
- **Edge cases** — null ⇒ `ArgumentNullException.ThrowIfNull`; empty/short seed ⇒ 0; empty/null UTRs skipped; DNA normalised; all reproduced by tests CTX-TA-005…008. ✓

### Finding & fix (the one divergence)
- **DEFECT (fixed this session):** Garcia's definition requires **non-overlapping** sites. The original `CountSeedSitesInUtr` incremented once per qualifying 6mer-core anchor with a plain `for (i++)` scan, and the code comment asserted overlap "is avoided because each qualifying anchor yields exactly one mutually-exclusive site type." That reasoning is **false**: site-TYPE mutual exclusivity *per anchor* does not make *distinct anchor positions* non-overlapping in sequence space. For a self-similar / periodic 6mer core the same physical region matches at several offsets and was over-counted:
  - core `AAAAAA` (seed pos 2-7 = `UUUUUU`) vs `AAAAAAAAAA` → counted **5** (should be 1).
  - core `ACACAC` (period-2) vs `ACACACACACA` → counted **3** (should be 1).
  - **Fix:** `CountSeedSitesInUtr` now scans **left-to-right greedily** — on counting a site at `i`, advance to `i+6` so the next candidate core cannot overlap the just-counted site's 6mer-core footprint; non-qualifying positions still advance by 1. On non-self-overlapping cores (real miRNA seeds vs typical 3′UTRs) this is **identical** to the prior per-anchor count (all existing CTX-TA fixtures unchanged: 2, 5, 2, 1); on periodic cores it now yields the non-overlapping count (5→1, 3→1). The misleading comment + XML `<remarks>` were corrected to describe the greedy non-overlapping scan and cite Garcia explicitly.
  - **Lock test added — CTX-TA-009:** periodic seed `GUGUGUU` (core `ACACAC`): `ACACACACACA` ⇒ 1 (overlapping anchors collapse); `ACACACAGGACACACA` (two cores ≥6 apart) ⇒ 2.

### Pre-existing surface re-confirmed (still holds after the fix)
`FindTargetSites` reverse-complement targeting, exact-match seed scan, nt8/A1 classification, antiparallel geometry, 0-based-inclusive coordinates, monotone Grimson-proportional base scorer, and the context++ wiring all remain green and unchanged (M-/S-/E-/CTX- tests pass; the limitation-fix touched only the new TA methods + comments).

### Test quality audit
The 9 CTX-TA tests assert exact hand-derived values (`log10(2)`, `log10(5)`, exact per-site-type scaled contributions, the non-overlapping counts 1 and 2), not code echoes; they cover both public methods, both overloads' contract (null throws), and the edge cases (empty set, no sites, DNA, short seed, overlapping cores). Expected numbers trace to Garcia (2011) / Agarwal (2015) / the TargetScan parameter rows.
- Target-prediction class: **72 tests pass**.
- Full unfiltered `dotnet test Seqeron.sln -c Debug`: **Seqeron.Genomics.Tests = 18861 passed, 0 failed**; whole solution green; **0 warnings** on the changed project.

### Honest residual (open boundary — NOT a defect of this unit)
TA_3UTR is now computable from a caller-supplied 3′UTR set. The remaining context++ residual is honest and unchanged: **PCT per-family parameters, SPS, Len_ORF and ORF-8mer count are caller-supplied** (listed in `OmittedFeatures` when absent), and **no default human transcriptome / 3′UTR set is bundled** — the caller must supply the 3′UTR set over which abundance is counted. This is a declared scope boundary (LIMITATIONS.md §3, trimmed by the fix), not an error.

---

## Verdict & follow-ups
- **Stage A: PASS** — `TA_3UTR = log10(N)`, N = non-overlapping 8mer+7mer-m8+7mer-A1 sites, log10 scale, Agarwal min-max scaling — all confirmed verbatim against Garcia (2011) Online Methods (PMC3190056), Agarwal (2015) Table 1, and the TargetScan 7.0 reference perl.
- **Stage B: PASS** — the code realises the formula; the one divergence from the literal "non-overlapping" requirement (over-counting on periodic cores) was found and **completely fixed** (greedy non-overlapping scan + corrected comments + lock test CTX-TA-009); all existing fixtures unchanged; full suite green.
- **End-state: ✅ CLEAN** — defect found and fully fixed in-session; the residual (PCT params / SPS / Len_ORF / ORF8m caller-supplied; no bundled transcriptome) is an honest, declared boundary.
- Defect logged in FINDINGS_REGISTER.md.
