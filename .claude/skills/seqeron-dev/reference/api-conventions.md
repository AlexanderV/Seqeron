# API conventions — construction, Method-ID lookup, worked pipelines

Code-side conventions for the `Seqeron.Genomics` C# API. Rigor (provenance, cross-check, disclaimer)
is delegated to [`bio-rigor`](../../bio-rigor/SKILL.md); envelope mechanics are in
[`limitation-policy.md`](limitation-policy.md). All symbols below are verified against source
(paths relative to repo root `../../../../`).

## Construction: `TryCreate` vs constructor

`DnaSequence` (and `RnaSequence`, `ProteinSequence`) in `Seqeron.Genomics.Core` offer both:

```csharp
using Seqeron.Genomics.Core;

var dna = new DnaSequence("AAAGAATTCAAA");            // DnaSequence.cs:22 — throws on invalid input
if (!DnaSequence.TryCreate("ACGTNN", out var seq))   // DnaSequence.cs:129 — false on invalid, no throw
    return;                                            // seq is DnaSequence? here
double gc = seq!.GcContent();                          // DnaSequence.cs:82
```

Rule of thumb: **`TryCreate` for untrusted input** (parsers, user data, file records); the throwing
ctor only where a failure is a programming error you *want* to surface.

## Finding the C# entry point for any MCP tool (via `Method ID`)

Each per-tool doc `docs/mcp/tools/<server>/<tool>.md` has a **`Method ID`** row (`Type.Method`) and a
**source link** to the exact `.cs#Lnnn`. That pair is the MCP→C# bridge. Example
([`primer_melting_temperature.md`](../../../../docs/mcp/tools/moltools/primer_melting_temperature.md)):

| Row | Value |
|---|---|
| Method ID | `PrimerDesigner.CalculateMeltingTemperature` |
| Source | `PrimerDesigner.cs#L197` |

So the C# call is `PrimerDesigner.CalculateMeltingTemperature(primer)` in
`Seqeron.Genomics.MolTools`. When you don't know the tool name, ask
[`seqeron-discovery`](../../seqeron-discovery/SKILL.md) (it searches the tool + algorithm docs
without loading all 427 schemas). Algorithm contracts / invariants / formulas live under
[`docs/algorithms/`](../../../../docs/algorithms/) — read the contract before relying on an edge case.

## Worked pipeline 1 — global align + statistics

`SequenceAligner` is a `static class` in `Seqeron.Genomics.Alignment`; results are records in
`Seqeron.Genomics.Infrastructure`.

```csharp
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

var a = new DnaSequence("ACGTACGTACGT");
var b = new DnaSequence("ACGTTCGTACGT");

// GlobalAlign (Needleman–Wunsch). Overloads take DnaSequence or raw string; optional ScoringMatrix
// defaults to SequenceAligner.SimpleDna.  (SequenceAligner.cs:58, 72)
AlignmentResult aln = SequenceAligner.GlobalAlign(a, b);
//                    AlignmentResult(AlignedSequence1, AlignedSequence2, Score,
//                                    AlignmentType, StartPosition1/2, EndPosition1/2)  (AlignmentTypes.cs:27)

// Derive identity/similarity/gaps from the alignment. (SequenceAligner.cs:570)
AlignmentStatistics stats = SequenceAligner.CalculateStatistics(aln);
Console.WriteLine($"Score={aln.Score}  Identity={stats.Identity:F1}%  Gaps={stats.Gaps}");
//                 AlignmentStatistics(Matches, Mismatches, Gaps, AlignmentLength,
//                                     Identity, Similarity, GapPercent)  (AlignmentTypes.cs:43)
```

Positions in `AlignmentResult` are **0-based** (a `bio-rigor` unit/coordinate check). For MSA use
`SequenceAligner.MultipleAlign` / `MultipleAlignProgressive` (`SequenceAligner.cs:702, 1170`),
which return `MultipleAlignmentResult(AlignedSequences, Consensus, TotalScore)`.

## Worked pipeline 2 — primer melting temperature (Tm)

```csharp
using Seqeron.Genomics.MolTools;

string primer = "GATCGATCGGGATCCAA";

// Wallace (<14 nt) / Marmur–Doty (≥14 nt); non-ACGT ignored. (PrimerDesigner.cs:197)
double tm = PrimerDesigner.CalculateMeltingTemperature(primer);            // °C

// Salt-corrected variant. (PrimerDesigner.cs:227)
double tmSalt = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, naConcentration: 50);

// Nearest-neighbour thermodynamic Tm (more accurate). (PrimerDesigner.cs:664)
double tmNN = PrimerDesigner.CalculateMeltingTemperatureNN(primer);
Console.WriteLine($"Tm(Wallace/MD)={tm:F1}°C  Tm(+salt)={tmSalt:F1}°C  Tm(NN)={tmNN:F1}°C");
```

These back the MCP tools `primer_melting_temperature` (Method ID `PrimerDesigner.CalculateMeltingTemperature`)
and `primer_melting_temperature_salt`.

## Worked pipeline 3 — validated construction + composition

Mirrors the README "validation-friendly path": validate first, then compute — never derive from a
sequence you haven't validated (a `bio-rigor` unit/validity check).

```csharp
using Seqeron.Genomics.Core;

foreach (var raw in userSuppliedSequences)
{
    if (!DnaSequence.TryCreate(raw, out var dna))   // DnaSequence.cs:129
    {
        Console.WriteLine($"skip invalid: {raw}");
        continue;
    }
    Console.WriteLine($"len={dna.Length}  GC%={dna.GcContent():F2}  revcomp={dna.ReverseComplement()}");
}
```

## Provenance

Every result you report from these calls carries a `bio-rigor` provenance block: the Method IDs /
C# calls in order, their parameters, and — if you touched a guarded unit — the `Envelope:` line with
the chosen `LimitationMode`. See [`bio-rigor`](../../bio-rigor/SKILL.md).
