using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit R14-1, B09 F288): ONE convention for <see cref="AlignmentResult.EndPosition1"/> / <see cref="AlignmentResult.EndPosition2"/> across every producer:
/// 0-based INCLUSIVE (length − 1 for a global result, the last consumed index of a local / fitting result, <c>Start − 1</c> for an empty consumed interval, −1 for an empty local result).
/// Before the fix NeedleAlign / WaterAlign / StretcherAlign / Global- and LocalAlignDualAffine returned the exclusive end (R14-1: "TTACGTTT" vs "ACGT": Needle 8,4, Water 2..6, Dual 2..6 against GlobalAlign 7,3 and LocalAlign 2..5).
/// The contract is checked over a corpus of 60 seeded pairs for every aligner of the unit.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class AlignmentResult_EndConvention_Tests
{
    private static readonly List<(string A, string B)> Corpus = BuildCorpus();

    private static List<(string, string)> BuildCorpus()
    {
        var list = new List<(string, string)> { ("TTACGTTT", "ACGT"), ("ACGT", "TTACGTTT"), ("ACGT", "ACGT"), ("AAAA", "CCCC"), ("A", "A"), ("A", "C"), ("ACGTACGT", "ACGTTACGT") };
        var rnd = new Random(14001);
        string R(int n) => new([.. Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(4)])]);
        for (int i = 0; i < 60; i++)
            list.Add((R(rnd.Next(1, 26)), R(rnd.Next(1, 26))));
        return list;
    }

    private static readonly DualAffineScoring Dual = DualAffineScoring.Minimap2;

    private static IEnumerable<TestCaseData> Producers()
    {
        yield return new TestCaseData("GlobalAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlign(a, b))).SetName("End_GlobalAlign");
        yield return new TestCaseData("GlobalAlignAffine", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlignAffine(a, b))).SetName("End_GlobalAlignAffine");
        yield return new TestCaseData("GlobalAlignLinearSpace", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlignLinearSpace(a, b, null, 4))).SetName("End_GlobalAlignLinearSpace");
        yield return new TestCaseData("GlobalAlignGeneralGap", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlignGeneralGap(a, b, k => -2.0 - k))).SetName("End_GlobalAlignGeneralGap");
        yield return new TestCaseData("GlobalAlignCoOptimal", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlignCoOptimal(a, b)[0])).SetName("End_GlobalAlignCoOptimal");
        yield return new TestCaseData("LocalAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlign(a, b))).SetName("End_LocalAlign");
        yield return new TestCaseData("LocalAlignAffine", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignAffine(a, b))).SetName("End_LocalAlignAffine");
        yield return new TestCaseData("LocalAlignLinearSpace", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignLinearSpace(a, b, null, 4))).SetName("End_LocalAlignLinearSpace");
        yield return new TestCaseData("LocalAlignGeneralGap", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignGeneralGap(a, b, k => -2.0 - k))).SetName("End_LocalAlignGeneralGap");
        yield return new TestCaseData("LocalAlignTopN", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignTopN(a, b, null, 1) is { Count: > 0 } l ? l[0] : SequenceAligner.LocalAlign(a, b))).SetName("End_LocalAlignTopN");
        yield return new TestCaseData("LocalAlignAffineTopN(Sim)", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignAffineTopN(a, b, null, 1, 0, TopNTieOrder.Sim) is { Count: > 0 } l ? l[0] : SequenceAligner.LocalAlign(a, b))).SetName("End_LocalAlignAffineTopN_Sim");
        yield return new TestCaseData("SemiGlobalAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.SemiGlobalAlign(a, b))).SetName("End_SemiGlobalAlign");
        yield return new TestCaseData("SemiGlobalAlignAffine", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.SemiGlobalAlignAffine(a, b))).SetName("End_SemiGlobalAlignAffine");
        yield return new TestCaseData("SemiGlobalAlign(EndGaps.All)", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.SemiGlobalAlign(a, b, (ScoringMatrix?)null, GapScoreProfile.FromEndGaps(EndGaps.All)))).SetName("End_SemiGlobalAlign_EndsFree");
        yield return new TestCaseData("NeedleAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.NeedleAlign(a, b))).SetName("End_NeedleAlign");
        yield return new TestCaseData("WaterAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.WaterAlign(a, b))).SetName("End_WaterAlign");
        yield return new TestCaseData("StretcherAlign", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.StretcherAlign(a, b))).SetName("End_StretcherAlign");
        yield return new TestCaseData("GlobalAlignDualAffine", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.GlobalAlignDualAffine(a, b, Dual))).SetName("End_GlobalAlignDualAffine");
        yield return new TestCaseData("LocalAlignDualAffine", (Func<string, string, AlignmentResult>)((a, b) => SequenceAligner.LocalAlignDualAffine(a, b, Dual))).SetName("End_LocalAlignDualAffine");
    }

    private static string Ungapped(string row) => row.Replace("-", "");

    /// <summary>The single End rule: Start − 1 ≤ End ≤ length − 1 and the consumed interval [Start, End] holds exactly the residues the rows align (rows that carry the whole sequence list the free flanks too).</summary>
    private static void AssertContract(string name, string s, string row, int start, int end, AlignmentType type, string label)
    {
        string ctx = $"{name} {label}";
        string residues = Ungapped(row);
        Assert.That(end, Is.LessThanOrEqualTo(s.Length - 1), $"{ctx}: End is the INCLUSIVE last index (≤ length − 1), got {end} for length {s.Length}");
        Assert.That(end, Is.GreaterThanOrEqualTo(start - 1), $"{ctx}: End ≥ Start − 1");
        if (residues.Length == 0 && type == AlignmentType.Local)
        {
            Assert.That((start, end), Is.EqualTo((-1, -1)), $"{ctx}: an empty local result is −1/−1");
            return;
        }

        Assert.That(start, Is.GreaterThanOrEqualTo(0), ctx);
        int consumed = end - start + 1;
        if (residues.Length == s.Length)
        {
            Assert.That(residues, Is.EqualTo(s), $"{ctx}: a full row is the sequence");
            Assert.That(consumed, Is.LessThanOrEqualTo(residues.Length), ctx);
            if (type == AlignmentType.Global)
                Assert.That((start, end), Is.EqualTo((0, s.Length - 1)), $"{ctx}: global = 0 .. length − 1");
        }
        else
        {
            Assert.That(consumed, Is.EqualTo(residues.Length), $"{ctx}: a trimmed row holds exactly [Start, End]");
            Assert.That(s.Substring(start, consumed), Is.EqualTo(residues), $"{ctx}: [Start, End] is the row");
        }
    }

    [TestCaseSource(nameof(Producers))]
    public void EveryProducer_EndIsZeroBasedInclusive(string name, Func<string, string, AlignmentResult> run)
    {
        foreach (var (a, b) in Corpus)
        {
            var r = run(a, b);
            string label = $"({a} / {b})";
            AssertContract(name, a, r.AlignedSequence1, r.StartPosition1, r.EndPosition1, r.AlignmentType, label + " seq1");
            AssertContract(name, b, r.AlignedSequence2, r.StartPosition2, r.EndPosition2, r.AlignmentType, label + " seq2");
        }
    }

    [Test]
    public void TheMeasuredPair_TTACGTTT_ACGT_HasTheSameEndConventionInAllProducers()
    {
        // R14-1 measurement: GlobalAlign end 7,3; LocalAlign 2..5 / 0..3; Needle was 8,4, Water 2..6, Dual 2..6 (exclusive)
        const string a = "TTACGTTT", b = "ACGT";
        Assert.That(SequenceAligner.GlobalAlign(a, b), Has.Property("EndPosition1").EqualTo(7).And.Property("EndPosition2").EqualTo(3));
        var needle = SequenceAligner.NeedleAlign(a, b);
        Assert.That((needle.StartPosition1, needle.StartPosition2, needle.EndPosition1, needle.EndPosition2), Is.EqualTo((0, 0, 7, 3)));
        var stretcher = SequenceAligner.StretcherAlign(a, b);
        Assert.That((stretcher.EndPosition1, stretcher.EndPosition2), Is.EqualTo((7, 3)));
        var dualGlobal = SequenceAligner.GlobalAlignDualAffine(a, b, Dual);
        Assert.That((dualGlobal.EndPosition1, dualGlobal.EndPosition2), Is.EqualTo((7, 3)));
        var local = SequenceAligner.LocalAlign(a, b);
        Assert.That((local.StartPosition1, local.EndPosition1, local.StartPosition2, local.EndPosition2), Is.EqualTo((2, 5, 0, 3)));
        var water = SequenceAligner.WaterAlign(a, b);
        Assert.That((water.StartPosition1, water.EndPosition1, water.StartPosition2, water.EndPosition2), Is.EqualTo((2, 5, 0, 3)));
        var dualLocal = SequenceAligner.LocalAlignDualAffine(a, b, Dual);
        Assert.That((dualLocal.StartPosition1, dualLocal.EndPosition1, dualLocal.StartPosition2, dualLocal.EndPosition2), Is.EqualTo((2, 5, 0, 3)));
    }

    [Test]
    public void EmbossPairAlignment_Flanks_AreDerivedFromTheRows_NotFromTheEndConvention()
    {
        // local: rows hold the aligned part, offset = Start, offend = the rest
        var local = SequenceAligner.LocalAlign("TTACGTTT", "ACGT");
        var pair = EmbossPairAlignment.From(local, "a", 8, "b", 4);
        Assert.That((pair.First.Offset, pair.First.Offend, pair.Second.Offset, pair.Second.Offend), Is.EqualTo((2, 2, 0, 0)));
        // fitting: the rows carry the whole sequences (free flanks included), so there is no offset although Start2 = 2
        var fitting = SequenceAligner.SemiGlobalAlign("ACGT", "TTACGTTT");
        Assert.That(fitting.StartPosition2, Is.EqualTo(2));
        var fitted = EmbossPairAlignment.From(fitting, "q", 4, "ref", 8);
        Assert.That((fitted.First.Offset, fitted.First.Offend, fitted.Second.Offset, fitted.Second.Offend), Is.EqualTo((0, 0, 0, 0)));
        // empty local result: the whole sequence lies after the (empty) row
        var none = SequenceAligner.LocalAlign("AAAA", "CCCC");
        var empty = EmbossPairAlignment.From(none, "a", 4, "b", 4);
        Assert.That((empty.First.Offset, empty.First.Offend, empty.Second.Offset, empty.Second.Offend), Is.EqualTo((0, 4, 0, 4)));
    }

    [Test]
    public void EmbossPrograms_OffsetsOfNeedleWaterStretcher_AreUnchanged()
    {
        // needle / stretcher: full rows; water: the local part (the oracle of the programs is EmbossPrograms_Oracle_Tests and the committed needle / water / stretcher oracles)
        const string a = "TTACGTTT", b = "ACGT";
        foreach (var r in new[] { SequenceAligner.NeedleAlign(a, b), SequenceAligner.StretcherAlign(a, b) })
        {
            var p = EmbossPairAlignment.From(r, "a", 8, "b", 4);
            Assert.That((p.First.Offset, p.First.Offend, p.Second.Offset, p.Second.Offend), Is.EqualTo((0, 0, 0, 0)));
        }

        var w = EmbossPairAlignment.From(SequenceAligner.WaterAlign(a, b), "a", 8, "b", 4);
        Assert.That((w.First.Offset, w.First.Offend, w.Second.Offset, w.Second.Offend), Is.EqualTo((2, 2, 0, 0)));
    }
}
