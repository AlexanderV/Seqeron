# Validation Report: SPLICE-PREDICT-001 — Gene Structure Prediction (intron/exon + spliced sequence)

- **Validated:** 2026-06-24   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.PredictGeneStructure`, `PredictIntrons`, and internal helpers `SelectNonOverlappingIntrons`, `DeriveExons`, `GenerateSplicedSequence` (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Context

A prior fix (commit `8b8c3cfe`) changed `GenerateSplicedSequence` to concatenate the **same
reported-exon set** that `DeriveExons` produces (consistent-filter), so that
`splicedSequence == concat(reportedExons)` and the coverage/length invariants (INV-3 / INV-4 /
INV-5) hold *by construction* for every `minExonLength`. This session re-confirms that
independently against external sources, by hand, and by code trace.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — RNA splicing** (fetched 2026-06-24): the **GT-AG rule** — the 5' donor site is
  the "almost invariant sequence GU at the 5' end of the intron"; the 3' acceptor site
  "terminates with an almost invariant AG sequence". mRNA assembly: the second
  transesterification "joins the exons and releases the intron lariat", yielding "a continuous
  mature mRNA molecule composed of spliced exons, with all intervening introns removed". This
  confirms: donor = GU, acceptor = AG, and **mature mRNA = exon concatenation, introns excised**.
- TestSpec / Evidence doc sources (Breathnach & Chambon 1981; Burge et al. 1999; Gilbert 1978;
  Alberts et al. 2002) are authoritative rank-1 references and state what is attributed to them:
  >99% GT-AG introns, donor consensus MAG|GURAGU, acceptor (Y)nNCAG|G, U2 = GT-AG, exon types
  Initial/Internal/Terminal/Single, phase = (Σ preceding exon lengths) mod 3.

### Formula / model check
- Intron model: donor (GU) paired with downstream acceptor (AG); intron spans
  [donor.Position .. acceptor.Position]; exons are the complementary regions. Matches GT-AG rule.
- Spliced sequence = concatenation of exon sequences — matches the splicing definition (S1/Gilbert).
- Phase = cumulative preceding-exon length mod 3 — standard (Alberts 2002); first exon → 0 trivially.

### Edge-case semantics
- Empty/null → empty structure (0 exons, 0 introns, "" spliced, 0 score). Defined & trivial.
- No splice sites → single Single-type exon spanning [0, len-1]. Sourced (single-exon genes, S6).
- Overlapping intron candidates → greedy non-overlapping selection (design decision; biology
  prescribes only the non-overlap invariant, not the algorithm — correctly documented as such).
- Sub-`minExonLength` inter-intron regions are dropped from the exon list (parameter contract).

### Invariants
INV-1..INV-10 are genuine properties. The two load-bearing ones for this re-confirmation:
- **INV-4** spliced = concat(reported exons), **INV-3** reported-exon + intron coverage internally
  consistent, **INV-5** spliced length = total − Σ intron length. The "consistent-filter" design
  makes these hold for all `minExonLength` (including when regions are dropped → spliced excludes
  them too). All are true under the current code (see Stage B).

### Findings / divergences
None. The description is biologically and mathematically sound.

## Stage B — Implementation

### Code path reviewed
- `PredictGeneStructure` (`SpliceSitePredictor.cs:518-562`) — orchestrator. Line 544 derives exons,
  line 551 calls `GenerateSplicedSequence(exons)` on the **same** list (this is the prior fix).
- `DeriveExons` (`:592-644`) — Initial/Internal/Terminal/Single typing, `minExonLength` filter,
  phase via `CalculatePhase` (`:646-650`).
- `GenerateSplicedSequence(List<Exon>)` (`:659-669`) — concatenates `exon.Sequence` in start order.
- `PredictIntrons` (`:439-489`), `SelectNonOverlappingIntrons` (`:564-590`).

### Formula realised correctly? (evidence)
- Donor pairing uses GU dinucleotide; acceptor uses AG (`FindDonorSites`/`FindAcceptorSites`).
  Acceptor `Position = i+1` is the G of AG (last intron nt); intron `[donor..acceptor]`,
  `Length = acceptor.Position − donor.Position + 1`. Confirms intron flanked by GU…AG.
- `GenerateSplicedSequence` now consumes the reported `exons` list, so spliced ≡ concat(exons)
  by construction — INV-4 holds independent of any filtering.

### Hand cross-check (designed pre-mRNA, 2 introns, GU…AG each)
Sequence: `flank50 + intron83 + mid40 + intron83 + flank50` = 306 nt, introns at
[50..132] and [173..255], with `intron83 = GUAAGU + 60·A + 14·U + CAG` (starts GU, ends AG).
- `minExonLength = 50`: DeriveExons → Initial[0..49] (50), mid[133..172] (40 < 50, **dropped**),
  Terminal[256..305] (50). Reported exons = {flank, flank}. `GenerateSplicedSequence` →
  flank+flank (100 nt); mid excluded from BOTH exon list and spliced product. INV-4: spliced ==
  concat(reported exons) ✓. INV-3 internally consistent ✓ (dropped region not double-counted).
- `minExonLength = 10`: mid kept → 3 exons; spliced = flank+mid+flank (140 nt); exon+intron
  coverage = 306 = full length (nothing dropped) ✓.
- Single-intron M3 case (153 nt): intron [35..117] len 83 (GU…AG), exons [0..34]+[118..152],
  spliced = Exon1+Exon2 = 70 = 153 − 83 ✓ (INV-5).
All three match the locked test assertions exactly.

### Cross-verification vs code (tests run)
`dotnet test --filter SpliceSitePredictor_GeneStructure_Tests` → **24 passed, 0 failed**. These
include the three dedicated drop/keep/all-dropped cases that assert
`SplicedSequence == concat(reported exons)`, length == Σ reported-exon lengths, and
`Does.Not.Contain(MidRegion40)` when dropped. No code was changed, so the full suite was not re-run.

### Variant/delegate consistency
No-intron path: `DeriveExons` returns one Single exon over the whole (uppercased, T→U) sequence,
so `GenerateSplicedSequence` returns the full sequence — equivalent to the old `introns==0 →
return sequence` branch; C2 identity test confirms. DNA (T) and lowercase inputs normalize to the
same RNA result (S2/S5 tests pass).

### Test quality audit
Assertions check exact sourced values (intron boundaries, GU/AG termini, exact spliced content
Exon1+Exon2, phase = 35 mod 3 = 2, drop/keep coverage), are deterministic, and the INV-3/INV-4
consistency helper (`AssertSplicedConsistency`) is applied to the dropped-region cases. Non-vacuous
(prerequisite count assertions guard the meaningful branches). Strong.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.** The prior consistent-filter fix is
  independently re-confirmed: donor (GU) + acceptor (AG) define introns, exons are the
  complementary regions, the spliced sequence is exactly the concatenation of the reported exons,
  and INV-3/INV-4/INV-5 hold by construction for all `minExonLength` (including dropped
  sub-threshold regions, which are excluded from both the exon list and the spliced product).
- No code changed; no follow-ups.
