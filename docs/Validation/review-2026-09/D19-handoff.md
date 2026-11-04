# D19 hand-off package — gap/bulge-tolerant off-target search (Cas-OFFinder)

**Why this document exists.** D19 is the last open DOABLE item of B08 round 1. Three attempts to do it in this
environment were halted by a safety classifier — two subagents in resume-4 (one of them during read-only recon of the
existing code) and the lead's own context in resume-5. Nothing was written: `CrisprDesigner.cs` is untouched by D19.
Per the campaign rule of 2026-10-21 ("all DOABLE before the duplication sweep, without exception; a halted item becomes
a hand-off package"), the item stays **DOABLE**, B08 is **not finished**, and the duplication sweep must not precede it.

**What is asked of you:** run this in another environment (your own machine, or any session that is not halted), then
return the four artefacts listed in "What to send back". The next B08 session integrates them under the normal campaign
rules — sourced, oracle-checked, tests, fast tier — and no result is accepted without its oracle numbers.

---

## 1. The algorithm and its source

Primary source: **Bae S, Park J, Kim J-S. "Cas-OFFinder: a fast and versatile algorithm that searches for potential
off-target sites of Cas9 RNA-guided endonucleases." *Bioinformatics* 30(10):1473–1475 (2014)**,
doi:10.1093/bioinformatics/btu048 — plus the implementation at <https://github.com/snugel/cas-offinder> and the bulge
variant `cas-offinder-bulge`. Use the source and the paper, not a blog restatement.

What it adds over the substitution-only scan we ship today:
- a **DNA bulge** of size *b*: the genomic target carries *b* extra bases that the guide has no partner for, so the
  protospacer window is **longer** than the guide by *b*;
- an **RNA bulge** of size *b*: the guide carries *b* extra bases, so the window is **shorter** than the guide by *b*;
- both are bounded by `maxDnaBulge` / `maxRnaBulge`, and a hit is reported when the mismatch count outside the bulge is
  within `maxMismatches`. The reference enumerates window lengths from `guideLength + maxDnaBulge` down to
  `guideLength − maxRnaBulge` and aligns with the bulge at each admissible offset.

## 2. Where it goes in this repository

File: `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs` (line numbers are from branch
`claude/stoic-maxwell-0olr2z` at the time of writing and will drift — anchor on the member names):

| Member | Line | Role |
|---|---|---|
| `FindOffTargets` | 587 | public entry, guard wrapper |
| `ValidateScanArguments` | ~600 | the guard order every entry point shares (null/empty guide, null genome, `maxMismatches` 0–5, guide length, guide alphabet) |
| `FindOffTargetsCore` | 651 | projection of the PAM-annotated scan |
| `FindOffTargetsWithPamCore` | 662 | `(OffTargetSite, Pam)` pairs |
| `ScanOffTargetsCore` | 676 | **the scan itself**; fixed-length window, `CountMismatches`, CRISPOR's on-target skip |
| `FindPamSitesCore` | 60 | PAM enumeration; yields `PamSite` with a **fixed-length** `TargetSequence` |
| `CountMismatches` | 830 | equal-length comparison, `int.MaxValue` sentinel on length mismatch |
| `OffTargetSite` | 2694 | `(Position, Sequence, Mismatches, MismatchPositions, IsForwardStrand, OffTargetScore)` |
| `PamSite` | 2486 | `(Position, PamSequence, TargetSequence, TargetStart, IsForwardStrand, System)` |

**The structural problem:** `FindPamSitesCore` compares a window of exactly `system.GuideLength`, and `CountMismatches`
returns `int.MaxValue` for unequal lengths. Bulges need a **variable-length** protospacer window. So the work is a new
enumerator, not an edit of the existing one.

## 3. The shape the result must have

1. **A new core enumerator** (e.g. `ScanOffTargetsWithBulgesCore`) plus a public entry point
   (e.g. `FindOffTargetsWithBulges(string guideSequence, DnaSequence genome, int maxMismatches, int maxDnaBulge,
   int maxRnaBulge, CrisprSystemType systemType)`). **Leave the existing chain literally untouched** —
   `FindOffTargets`, `FindOffTargetsCore`, `FindOffTargetsWithPamCore`, `ScanOffTargetsCore` keep their current
   behaviour and output order. Reuse `ValidateScanArguments` for the shared guards and add the two bulge guards in the
   same style (`0 ≤ maxDnaBulge, maxRnaBulge ≤ 2`, the reference's practical ceiling; state the bound in the XML doc).
2. **A new sibling record**, not a change to `OffTargetSite` — that record is serialised through MCP and widening it
   breaks consumers. Something like
   `BulgedOffTargetSite(int Position, string AlignedGuide, string AlignedTarget, int Mismatches, BulgeType Bulge,
   int BulgeSize, int BulgePosition, bool IsForwardStrand)` with `-` for the gap in the aligned strings, and
   `enum BulgeType { None, Dna, Rna }`. **No value tuples in anything MCP-serialised**: F62 is the precedent — they
   serialise as `{}`.
3. **The three specificity scorers stay substitution-only** — `CalculateHeuristicSpecificityScore`,
   `CalculateMitSpecificityScore`, `CalculateCfdSpecificityScore`. Neither the MIT nor the CFD model defines a weight
   for a bulged position, so feeding bulged hits into them would be unsourced. Document that explicitly in the new
   method's XML doc and in `docs/algorithms/MolTools/Off_Target_Analysis.md`.
4. **Amend the one stale doc line:** `FindOffTargets`'s XML doc says "Substitutions only — bulges/indels are out of
   scope"; it must now point at the new entry point instead of claiming the library cannot do it.
5. Optional, only if it is natural: an MCP tool wrapper. If you add one, the catalog must be regenerated
   (`python3 scripts/skills/gen-catalog.py`, and `--check` must pass).

## 4. The oracle — this is the part that cannot be skipped

Every number in this project comes from running the reference, never from reasoning about it. Please produce the
comparison yourself and send the raw output.

```bash
# option A — the reference implementation (needs OpenCL; conda is the easy route)
conda install -c bioconda cas-offinder          # provides cas-offinder and cas-offinder-bulge
# input file (Cas-OFFinder format): line 1 = genome file or directory, line 2 = PAM pattern,
# then one line per guide: <guide+Ns> <mismatches>
printf '%s\n' /path/to/toy_genome.fa NNNNNNNNNNNNNNNNNNNNNRG 'GGCCGAGATGTCTCGCTCCGNNN 3' > in.txt
cas-offinder-bulge in.txt C out.txt             # C = CPU; add the DNA/RNA bulge sizes per its README
```

```bash
# option B — if OpenCL is unavailable, use any maintained Python port / the paper's worked examples,
# and say in your reply exactly which implementation and version produced the numbers.
```

Use a **small toy genome you can paste into the Evidence file** (a few hundred bp with hand-placed bulged targets) so
the comparison is reproducible here, plus one realistic case if you can. Record: the exact tool version, the input
file, the full output, and the command line.

## 5. Tests that must exist and pass

- New unit tests in `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_OffTarget_Tests.cs` (or a new
  `..._OffTargetBulge_Tests.cs`) locking the oracle rows **bit-exactly**: per case the guide, the target, the bulge type,
  size and position, and the mismatch count.
- Guard tests mirroring the existing ones: null/empty guide, null genome, `maxMismatches` out of range, bulge sizes out
  of range, bad alphabet — same exception types and order as `ValidateScanArguments`.
- A test asserting the **existing** `FindOffTargets` output is unchanged for a case that has bulged hits nearby
  (regression guard for point 3.1).
- `tests/TestSpecs/CRISPR-OFF-001.md` gains the new rows; `docs/Evidence/CRISPR-OFF-001-Evidence.md` gains the oracle
  section (tool version, command, input, output, and the comparison table).
- Fast tier must be green: `dotnet test tests/Seqeron/Seqeron.Genomics.Tests` filtered as the campaign does
  (Fuzzing/Properties/Metamorphic excluded), plus `Seqeron.Mcp.MolTools.Tests` if you touched the MCP layer. The current
  baseline is **18639 passed / 0 failed** and MolTools MCP **441/0** — report your numbers.

## 6. What to send back

1. A **patch** (`git format-patch` or a plain diff) against branch `claude/stoic-maxwell-0olr2z`, or a pushed branch name.
2. The **raw oracle output** and the exact command line + tool version that produced it.
3. The **test run output** (counts, and the new tests' names).
4. One paragraph on anything you judged differently from this document and why — that goes into the Evidence as a
   decision record.

On integration I will re-run your oracle command here if the tool installs, check every test expectation against your
recorded numbers rather than against the code, re-run the fast tier, and only then tick D19 and let the duplication
sweep proceed.
