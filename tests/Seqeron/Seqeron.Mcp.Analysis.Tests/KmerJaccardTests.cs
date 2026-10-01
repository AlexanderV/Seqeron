using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_jaccard</c> MCP tool. Expected values: Python set Jaccard, sourmash 4.9.4
/// MinHash(scaled=1) and the Mash 2.3 binary (<c>mash dist -s 100000 [-n]</c>), not the wrapper output.
/// </summary>
[TestFixture]
public class KmerJaccardTests
{
    [Test]
    public void KmerJaccard_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerJaccard("ACGT", "ACGA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("", "ACGT", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("ACGT", null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("ACGT", "ACGT", 0));
    }

    [Test]
    public void KmerJaccard_Binding_MatchesMashAndSourmash()
    {
        // ACGTTGCAACGGT / ACGTAGCATCGGTA k=2: literal J = 6/13; canonical J = 0.5 (mash 5/10, sourmash 0.5),
        // Mash distance 0.20273255405408222 (mash dist 0.202733); mash -n (acgtOnly) 0.229766.
        var literal = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2);
        var canon = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2, canonical: true);
        var acgt = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2, acgtOnly: true);
        Assert.Multiple(() =>
        {
            Assert.That(literal.Jaccard, Is.EqualTo(6.0 / 13).Within(1e-12));
            Assert.That(canon.Jaccard, Is.EqualTo(0.5));
            Assert.That(canon.MashDistance, Is.EqualTo(0.20273255405408222).Within(1e-12));
            Assert.That(acgt.MashDistance, Is.EqualTo(0.2297661646892201).Within(1e-12));
        });
    }
}
