# Validation Report: ONCO-FUSION-001 — Fusion Gene Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectFusions(IEnumerable<FusionCandidate>, FusionDetectionThresholds?)`, `OncologyAnalyzer.IsInFrame(int, int)`, `OncologyAnalyzer.ComputeTotalSupport(FusionCandidate)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

## Scope

ONCO-FUSION-001 is a candidate gene-fusion *caller*: it consumes already-grouped breakpoint
candidates with per-class supporting-read counts (Arriba schema: split_reads1, split_reads2,
discordant_mates) and applies the STAR-Fusion minimum-support decision rule to report passing
fusions, ordered by abundance, with a codon-phase reading-frame annotation. Raw-BAM chimeric-read
extraction is explicitly out of scope (separate unit). This is a deterministic rule, not a
statistical model, so every threshold and formula has an exact sourced expected value.

## Stage A — Description

### Sources opened this session (live fetches)

1. **STAR-Fusion driver source** (raw GitHub, fetched 2026-06-16) — quoted the actual variable
   assignments and help text:
   - `my $MIN_JUNCTION_READS = 1;` — "minimum number of junction-spanning reads required."
   - `my $MIN_SUM_FRAGS = 2; # requires at least one junction read, else see min_spanning_frags_only`
     — "minimum fusion support = ( # junction_reads + # spanning_frags )".
   - `my $MIN_SPANNING_FRAGS_ONLY = 5;` — "minimum number of rna-seq fragments required as fusion
     evidence if there are no junction reads".
   - `my $MIN_FFPM = 0.1;` (abundance filter; not part of the count rule under test).
   URL: https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
2. **Arriba output-file wiki** (fetched 2026-06-16) — confirmed verbatim:
   - split_reads1/2 = "number of supporting split fragments with an anchor in gene1 or gene2";
     anchor = "the gene to which the longer segment of the split read aligns".
   - discordant_mates = "pairs (fragments) of discordant mates (a.k.a. spanning reads or bridge
     reads) supporting the fusion".
   - **Total supporting reads** = "summing up the reads given in the columns split_reads1,
     split_reads2, discordant_mates" (+ filters; the unit sums the three primary classes).
   - reading_frame values: in-frame / out-of-frame / stop-codon / "." .
   URL: https://github.com/suhrig/arriba/wiki/05-Output-files
3. **Wikipedia "Reading frame"** (fetched 2026-06-16) — "There are three reading frames … each
   beginning from a different nucleotide in a triplet"; frame is determined by start position
   modulo 3. URL: https://en.wikipedia.org/wiki/Reading_frame
4. **In-frame fusion definition** (WebSearch 2026-06-16, Genomics England + SnapGene) — confirmed
   the exon-phase rule: "if one exon finishes after the second letter of a triplet (end phase 2),
   the next one should start with the third letter (start phase 2)", and that a productive fusion
   requires the 5' partner coding length mod 3 to be compatible with the 3' partner's phase.

### Formula check

- **Total support** = `split1 + split2 + discordant` — matches Arriba spec exactly.
- **Detection rule:** report iff `gene5p ≠ gene3p` AND
  (junctionReads ≥ MIN_JUNCTION_READS=1 ⇒ totalSupport ≥ MIN_SUM_FRAGS=2)
  OR (junctionReads = 0 ⇒ discordant ≥ MIN_SPANNING_FRAGS_ONLY=5) — matches STAR-Fusion defaults
  1 / 2 / 5 exactly.
- **In-frame:** `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0`, i.e.
  `5pBases mod 3 == 3pPhase` — the standard modulo-3 / exon-phase rule. Correct.

### Edge-case semantics (all sourced)

- Spanning-only (0 junction reads) requires ≥5 discordant (STAR-Fusion); 4 is filtered.
- junction=1 but total=1 (<2) → rejected (STAR-Fusion min_sum_frags).
- gene5p == gene3p → not a fusion (Registry invariant + nomenclature).
- Unknown coding phase (−1) → ReadingFrame=Unknown (never guessed).

### Independent cross-check (exact numbers, hand-computed from sourced rules)

| Candidate | split1 | split2 | disc | junc | total | rule | Expected |
|-----------|--------|--------|------|------|-------|------|----------|
| EML4-ALK | 3 | 2 | 4 | 5 | 9 | junc≥1 ∧ 9≥2 | DETECTED |
| CD74-ROS1 | 0 | 0 | 5 | 0 | 5 | 5≥5 | DETECTED |
| NCOA4-RET | 0 | 0 | 4 | 0 | 4 | 4<5 | REJECTED |
| KIF5B-RET | 1 | 0 | 0 | 1 | 1 | junc≥1 ∧ 1<2 | REJECTED |
| TMPRSS2-ERG | 1 | 0 | 1 | 1 | 2 | junc≥1 ∧ 2≥2 | DETECTED |
| ALK-ALK | 5 | 5 | 5 | 10 | 15 | gene5p==gene3p | REJECTED |

Reading-frame: (300,0)→0→in-frame; (301,0)→1→out; (302,0)→2→out; (301,1)→0→in. Gene partner
orders (EML4/BCR/TMPRSS2/CD74/KIF5B/NCOA4 = 5'; ALK/ABL1/ERG/ROS1/RET = 3' kinase) confirmed
biologically correct by literature.

### Findings / divergences

None substantive. The unit honestly documents two assumptions (candidate-level counts rather
than raw BAM; codon-phase in-frame test does not scan for premature stop codons — Arriba's
"stop-codon" value is out of scope, requiring transcript reconstruction). Both are source-grounded
and explicitly recorded. Stage A PASS.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `FusionCandidate` record struct (L3384) and `FusionDetectionThresholds` (L3399) with explicit
  parameterless ctor restoring the 1/2/5 defaults (guards against `default(T)` zeroing thresholds).
- `ComputeTotalSupport` (L3447) = `SplitReads5Prime + SplitReads3Prime + DiscordantMates`.
- `IsInFrame` (L3464): validates `fivePrimeCodingBases ≥ 0` and `phase ∈ {0,1,2}`
  (ArgumentOutOfRangeException otherwise); returns `(5p − phase) % 3 == 0`.
- `DetectFusions` (L3503): null-guards input; per candidate rejects negative counts
  (ArgumentException), skips same-gene (case-insensitive), applies the branch rule, resolves frame,
  and returns ordered by descending TotalSupport with (Gene5Prime, Gene3Prime) ordinal tie-break.
- `ResolveReadingFrame` (L3557): returns Unknown when phase fields are unset/invalid.

### Formula realised correctly?

Yes. The detection branch `junctionReads >= MinJunctionReads ? totalSupport >= MinSumFrags :
discordant >= MinSpanningFragsOnly` exactly implements STAR-Fusion's two-regime rule under the
defaults. Total support and the in-frame modulo-3 test match the sourced formulae.

### Cross-verification table recomputed vs code

Every row of the table above reproduces under the actual code (confirmed by the passing tests
M1–M6, M11 ordering `{9,5,2}`, and frame cases). No value traces to "what the code returns" — each
is derived from STAR-Fusion 1/2/5, the Arriba sum, or modulo-3.

### Variant/delegate consistency

`ComputeTotalSupport` is reused inside `DetectFusions`; `ResolveReadingFrame` delegates to the same
`IsInFrame` (also reused by ONCO-FUSION-003's frame logic). Consistent.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M3 (4<5→REJECTED) and M4 (junc=1,total=1<2→REJECTED) would FAIL a
  naive implementation that merely sums reads or lacks the spanning-only branch — they encode the
  threshold rule, not the output. M1/M2/M5/M7 assert exact TotalSupport (9/5/2/9) from the Arriba
  sum. M8–M10/S4 assert exact modulo-3 outcomes.
- **No green-washing:** all assertions are exact (`Has.Count.EqualTo`, `Is.EqualTo(new[]{9,5,2})`,
  `Is.Empty`, `Is.True/False`); no Greater/AtLeast/ranges, no widened tolerances, no skips.
- **Coverage:** all three public methods exercised; branches covered — junction-pass, spanning-pass,
  spanning-fail, sum-fail, same-gene, ordering, both frame outcomes, Unknown frame, custom threshold
  (S1), both boundaries (S2 =5, S3 =2); error cases null (M12), negative count (M13), empty (M14),
  invalid IsInFrame args (C1: negative bases, phase=3). 19 tests, no gaps vs TestSpec §4.
- **Honest green:** full unfiltered suite = **Failed: 0, Passed: 6644**; `dotnet build` 0 errors
  (4 pre-existing warnings unrelated to this unit, none in changed files — no files changed).

### Findings / defects

None. No code or test change was required.

## Verdict & follow-ups

- **Stage A: PASS.** Description matches STAR-Fusion (1/2/5), Arriba (total-support sum, anchor,
  discordant definition) and the modulo-3 reading frame, all re-fetched live this session.
- **Stage B: PASS.** Implementation faithfully realises the rule; tests assert exact sourced values,
  cover every method/branch/edge, and the full suite is green.
- **Test-quality gate: PASS.**
- **End-state: ✅ CLEAN** — no defect found; algorithm fully functional.
- Note (not a defect): premature-stop-codon ("stop-codon" reading_frame) detection is deliberately
  deferred to the transcript-reconstruction unit (ONCO-FUSION-003); documented as an assumption.
