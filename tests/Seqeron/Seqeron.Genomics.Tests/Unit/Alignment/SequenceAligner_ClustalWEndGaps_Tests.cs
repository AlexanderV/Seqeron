using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10d part 2, B09 F67): the end-gap options of ClustalW 2.1's profile stage, <c>UserParameters::useEndGaps</c> (<see cref="ClustalWOptions.UseEndGaps"/>)
/// and <c>endGapPenalties</c> (<see cref="ClustalWOptions.EndGapPenalties"/>), and the program's abort when a first-pass profile alignment scores below zero
/// (<see cref="ClustalWAbortedException"/>).
///
/// Oracle (TestData/Alignment/clustalw_end_gaps_oracle.tsv, scripts/alignment/build_clustalw_end_gaps_oracle.py): clustalw 2.1+lgpl-7 built from the Debian source with the single edit
/// that the <c>-ENDGAPS</c> handler of <c>CommandLineParser.cpp</c> reads <c>CW_USEENDGAPS</c> / <c>CW_ENDGAPPEN</c> (Docker package docs/Validation/review-2026-09/handoff/A1-10d, run
/// by the user): 340 sets (31 DNA + 3 protein sets are among the 34 aborting ones; the sets are the committed DNA / protein ClustalW oracle sets with their option variants) x
/// 4 modes; the control mode equals the installed binary on 340 / 340 sets. In the modes endgappen / both the program writes NO output (exit status 0) for 34 sets: the score of a
/// first-pass profile alignment is negative, <c>MSA::multiSeqAlign</c> returns -1 and <c>Clustal::align</c> returns (MSA.cpp 150-155, Clustal.cpp 205-208).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWEndGaps_Tests
{
    private sealed record Row(string Id, string Type, string Mode, string[] Seqs, string Options, string[]? Rows);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_end_gaps_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), f[4], f[5] == "ABORT" ? null : f[5].Split(';')));
        }

        return list;
    });

    private static (bool UseEndGaps, bool EndGapPenalties) Flags(string mode) => mode switch
    {
        "control" => (false, false),
        "useendgaps" => (true, false),
        "endgappen" => (false, true),
        "both" => (true, true),
        _ => throw new ArgumentException(mode),
    };

    private static string[] AlignOrThrow(Row r, bool useEndGaps, bool endGapPenalties)
    {
        if (r.Type == "DNA")
        {
            var o = SequenceAligner_ClustalWMultiple_Tests.Parse(r.Options) with { UseEndGaps = useEndGaps, EndGapPenalties = endGapPenalties };
            return SequenceAligner.MultipleAlignClustalW(r.Seqs, o).AlignedSequences.ToArray();
        }

        var p = SequenceAligner_ClustalWProteinMultiple_Tests.Parse(r.Options) with { UseEndGaps = useEndGaps, EndGapPenalties = endGapPenalties };
        return SequenceAligner.MultipleAlignClustalW(r.Seqs, p).AlignedSequences.ToArray();
    }

    /// <summary>Runs the library on an oracle row; <c>null</c> = ClustalW's abort.</summary>
    private static string[]? Align(Row r)
    {
        var (useEndGaps, endGapPenalties) = Flags(r.Mode);
        try
        {
            return AlignOrThrow(r, useEndGaps, endGapPenalties);
        }
        catch (ClustalWAbortedException)
        {
            return null;
        }
    }

    [Test]
    public void OracleHasEveryModeForEverySet()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(1360));
            foreach (string mode in new[] { "control", "useendgaps", "endgappen", "both" })
                Assert.That(rows.Count(r => r.Mode == mode), Is.EqualTo(340), mode);
            Assert.That(rows.Where(r => r.Mode == "control").Select(r => r.Id).Distinct().Count(), Is.EqualTo(340));
            Assert.That(rows.Count(r => r.Mode == "control" && r.Type == "DNA"), Is.EqualTo(120));
            Assert.That(rows.Count(r => r.Mode == "control" && r.Type == "PROTEIN"), Is.EqualTo(220));
            Assert.That(rows.Count(r => r.Rows is null), Is.EqualTo(68), "34 aborting sets x 2 modes");
        });
    }

    [TestCase("control")]
    [TestCase("useendgaps")]
    public void MultipleAlignClustalW_EndGapModes_RowsEqualClustalW_OnAll340Sets(string mode)
    {
        int bad = 0;
        string? first = null;
        foreach (var r in Oracle.Value.Where(r => r.Mode == mode))
        {
            var got = Align(r);
            Assert.That(r.Rows, Is.Not.Null, $"{mode} never aborts in the oracle");
            if ((got is null || !got.SequenceEqual(r.Rows!)) && ++bad == 1)
                first = $"set {r.Id} {r.Type} ({r.Options}) {mode}: expected {string.Join('/', r.Rows!)}, got {(got is null ? "abort" : string.Join('/', got))}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    [TestCase("endgappen")]
    [TestCase("both")]
    public void MultipleAlignClustalW_EndGapPenaltyModes_RowsEqualClustalW_OnAll306SetsAndAbortOn34(string mode)
    {
        var rows = Oracle.Value.Where(r => r.Mode == mode).ToList();
        int aligned = 0, aborted = 0, bad = 0;
        string? first = null;
        foreach (var r in rows)
        {
            var got = Align(r);
            if (r.Rows is null)
            {
                aborted++;
                if (got is not null && ++bad == 1)
                    first = $"set {r.Id} {r.Type} ({r.Options}) {mode}: ClustalW writes no alignment, got {string.Join('/', got)}";
            }
            else
            {
                aligned++;
                if ((got is null || !got.SequenceEqual(r.Rows)) && ++bad == 1)
                    first = $"set {r.Id} {r.Type} ({r.Options}) {mode}: expected {string.Join('/', r.Rows)}, got {(got is null ? "abort" : string.Join('/', got))}";
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(aligned, Is.EqualTo(306));
            Assert.That(aborted, Is.EqualTo(34));
            Assert.That(bad, Is.Zero, first);
        });
    }

    [Test]
    public void ControlMode_EqualsTheDefaultOptions()
    {
        // the default options never abort and carry both flags false
        Assert.That(new ClustalWOptions().UseEndGaps, Is.False);
        Assert.That(new ClustalWOptions().EndGapPenalties, Is.False);
        Assert.That(new ClustalWProteinOptions().UseEndGaps, Is.False);
        Assert.That(new ClustalWProteinOptions().EndGapPenalties, Is.False);
        foreach (var r in Oracle.Value.Where(r => r.Mode == "control").Take(40))
            Assert.That(AlignOrThrow(r, false, false), Is.EqualTo(r.Rows), r.Id);
    }

    /// <summary>The options really change the alignment (the oracle is not an accident of easy inputs): many sets differ from the control rows in each mode.</summary>
    [Test]
    public void OracleDiscriminatesTheOptions()
    {
        var byId = Oracle.Value.GroupBy(r => r.Id).ToDictionary(g => g.Key, g => g.ToDictionary(r => r.Mode));
        int Differing(string mode) => byId.Values.Count(m => m[mode].Rows is { } rows && !rows.SequenceEqual(m["control"].Rows!));
        Assert.Multiple(() =>
        {
            Assert.That(Differing("useendgaps"), Is.GreaterThan(100));
            Assert.That(Differing("endgappen"), Is.GreaterThan(150));
            Assert.That(Differing("both"), Is.GreaterThan(150));
            Assert.That(byId.Values.Count(m => m["useendgaps"].Rows is { } a && m["both"].Rows is { } b && !a.SequenceEqual(b)), Is.GreaterThan(100));
        });
    }

    /// <summary>The abort is the documented contract: a first-pass score below zero throws, carrying the group number and the score; never raised without end-gap penalties.</summary>
    [Test]
    public void Abort_IsRaisedExactlyForTheProgramsNoOutputSets_AndNeverWithoutEndGapPenalties()
    {
        var failing = Oracle.Value.Where(r => r.Mode == "endgappen" && r.Rows is null).ToList();
        Assert.That(failing, Has.Count.EqualTo(34));
        Assert.That(failing.Count(r => r.Type == "DNA"), Is.EqualTo(31));
        Assert.That(failing.Count(r => r.Type == "PROTEIN"), Is.EqualTo(3));
        var r0 = failing[0];
        var ex = Assert.Throws<ClustalWAbortedException>(() => AlignOrThrow(r0, false, true))!;
        Assert.Multiple(() =>
        {
            Assert.That(ex.Score, Is.LessThan(0));
            Assert.That(ex.Group, Is.GreaterThanOrEqualTo(1));
            Assert.That(ex, Is.InstanceOf<InvalidOperationException>());
        });

        // the same set without end-gap penalties (control, useendgaps) aligns
        foreach (var r in failing)
        {
            Assert.That(Oracle.Value.First(c => c.Id == r.Id && c.Mode == "control").Rows, Is.Not.Null);
            Assert.That(Oracle.Value.First(c => c.Id == r.Id && c.Mode == "useendgaps").Rows, Is.Not.Null);
        }
    }

    /// <summary>The profile-alignment score is truncated toward zero after the division by 100: scores in (-100, 0) do not abort (checked on the abort cases: |score| >= 1).</summary>
    [Test]
    public void Abort_ScoresOfAbortingSetsAreNegativeIntegers()
    {
        foreach (var r in Oracle.Value.Where(r => r.Mode == "both" && r.Rows is null).Take(12))
        {
            var ex = Assert.Throws<ClustalWAbortedException>(() => AlignOrThrow(r, true, true), r.Id)!;
            Assert.That(ex.Score, Is.LessThanOrEqualTo(-1), r.Id);
            Assert.That(ex.Message, Does.Contain("scores " + ex.Score));
        }
    }
}
