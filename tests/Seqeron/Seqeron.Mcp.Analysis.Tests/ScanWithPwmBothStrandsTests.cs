using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// scan_with_pwm_both_strands — delegation to <c>MotifFinder.ScanWithPwmBothStrands</c>.
/// Expected values: Biopython 1.88 search(both=True) on the Wikipedia PWM example, locked in
/// MotifFinder_PwmStrandsAndThresholds_Tests.
/// </summary>
[TestFixture]
public class ScanWithPwmBothStrandsTests
{
    internal static readonly string[] WikipediaSequences =
    {
        "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
        "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
    };

    private const string Target = "CCTAGGTAAGTAACAGGTCAGTGG";

    internal static PwmInput WikipediaPwm()
    {
        var pwm = AnalysisTools.CreatePwm(WikipediaSequences, 0.25);
        return new PwmInput(pwm.Matrix, pwm.Length);
    }

    [Test]
    public void ScanWithPwmBothStrands_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.ScanWithPwmBothStrands(Target, WikipediaPwm()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwmBothStrands("", WikipediaPwm()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwmBothStrands("ACGU", WikipediaPwm()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwmBothStrands(Target, null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwmBothStrands(Target, new PwmInput(new double[3][], 2)));
    }

    [Test]
    public void ScanWithPwmBothStrands_WikipediaExample_EqualsBiopython()
    {
        var items = AnalysisTools.ScanWithPwmBothStrands(Target, WikipediaPwm(), 0.0).Items;
        Assert.Multiple(() =>
        {
            Assert.That(items.Select(i => i.BiopythonPosition), Is.EqualTo(new[] { 2, 6, -16, 13 }));
            Assert.That(items.Select(i => i.Position), Is.EqualTo(new[] { 2, 6, 8, 13 }));
            Assert.That(items.Select(i => i.Strand), Is.EqualTo(new[] { "+", "+", "-", "+" }));
            Assert.That(items.Select(i => i.Score),
                Is.EqualTo(new[] { 11.70753002166748, 4.779160022735596, 2.387691020965576, 9.316061019897461 }).Within(1e-5));
            Assert.That(items[2].MatchedSequence, Is.EqualTo("CCTGTTACT"));
            Assert.That(items.All(i => i.Pattern == "TAGGTAAGT"), Is.True);
        });
    }

    [Test]
    public void ScanWithPwmBothStrands_EqualsLibrary()
    {
        var lib = global::Seqeron.Genomics.Analysis.MotifFinder.ScanWithPwmBothStrands(
            new global::Seqeron.Genomics.Core.DnaSequence(Target),
            global::Seqeron.Genomics.Analysis.MotifFinder.CreatePwm(WikipediaSequences, 0.25), -3.0).ToArray();
        var items = AnalysisTools.ScanWithPwmBothStrands(Target, WikipediaPwm(), -3.0).Items;
        Assert.That(items.Select(i => (i.Position, i.BiopythonPosition, i.Strand, i.MatchedSequence, i.Pattern, i.Score)),
            Is.EqualTo(lib.Select(m => (m.Position, m.BiopythonPosition, m.Strand.ToString(), m.MatchedSequence, m.Pattern, m.Score))));
    }
}
