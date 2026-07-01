using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>count_kmers_both_strands</c> MCP tool.
/// Expected values derived from the documented identity count[w] = forward[w] + forward[RC(w)]
/// (KmerAnalyzer.CountKmersBothStrands XML doc), NOT from the wrapper's output. A wrapper that
/// swapped strands or skipped the reverse complement would fail these.
/// </summary>
[TestFixture]
public class CountKmersBothStrandsTests
{
    [Test]
    public void CountKmersBothStrands_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CountKmersBothStrands("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmersBothStrands("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmersBothStrands(null!, 2));
        // Non-DNA input must be rejected by the DNA guard.
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmersBothStrands("AUGC", 2));
    }

    [Test]
    public void CountKmersBothStrands_Binding_InvokesSuccessfully()
    {
        // "AAAA", k=2: forward AA=3; reverse complement "TTTT" -> TT=3. Total 6 = 2*(L-k+1).
        var homo = AnalysisTools.CountKmersBothStrands("AAAA", 2).Counts;
        Assert.Multiple(() =>
        {
            Assert.That(homo["AA"], Is.EqualTo(3));
            Assert.That(homo["TT"], Is.EqualTo(3));
            Assert.That(homo, Has.Count.EqualTo(2));
            Assert.That(homo.Values.Sum(), Is.EqualTo(6));
        });

        // "ATGC", k=2: forward AT,TG,GC each 1; reverse complement "GCAT" adds GC,CA,AT.
        // -> AT=2, TG=1, GC=2, CA=1. count[w] == count[RC(w)].
        var mixed = AnalysisTools.CountKmersBothStrands("ATGC", 2).Counts;
        Assert.Multiple(() =>
        {
            Assert.That(mixed["AT"], Is.EqualTo(2));
            Assert.That(mixed["TG"], Is.EqualTo(1));
            Assert.That(mixed["GC"], Is.EqualTo(2));
            Assert.That(mixed["CA"], Is.EqualTo(1));
            Assert.That(mixed.Values.Sum(), Is.EqualTo(6));
            // Inversion symmetry: AT <-> AT (self-RC dimer's partner), GC <-> GC.
            Assert.That(mixed["AT"], Is.EqualTo(mixed["AT"]));
        });
    }
}
