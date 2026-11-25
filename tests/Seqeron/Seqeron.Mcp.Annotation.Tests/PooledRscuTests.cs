using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PooledRscuTests
{
    private static readonly string[] Cds =
    {
        "ATGCTTCTTCTGTTAAAAAAGTAA",
        "augcugcuuaaaaaaggcgcc",       // RNA, partial-free
        "ATGNNNCTTGGCGGCGGAGGG",       // ambiguous triplet skipped, frame kept
    };

    [Test]
    public void PooledRscu_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PooledRscu(Cds));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PooledRscu(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PooledRscu(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PooledRscu(Cds, geneticCode: 99));
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(11)]
    public void PooledRscu_EqualsCalculateRscuOnPooledCounts_WithoutStopCodons(int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        var pooled = new Dictionary<string, int>();
        foreach (var s in Cds)
            foreach (var (codon, n) in AnnotationTools.CodonUsage(s).Usage)
                pooled[codon] = pooled.GetValueOrDefault(codon) + n;
        var expected = CodonUsageAnalyzer.CalculateRscu(pooled, code);
        foreach (var (rna, aa) in code.CodonTable)
            if (aa == '*') expected.Remove(rna.Replace('U', 'T'));

        var actual = AnnotationTools.PooledRscu(Cds, table).Rscu;

        Assert.That(actual.Keys, Is.EquivalentTo(expected.Keys));
        foreach (var kv in expected)
            Assert.That(actual[kv.Key], Is.EqualTo(kv.Value).Within(1e-12), kv.Key);
    }

    [Test]
    public void PooledRscu_LeucineFamily_MatchesCodonWWorkedValues()
    {
        // CTT x3, CTG x1, TTA x1 pooled over two sequences: Leu n = 6, total 5 -> n*x/sum with n = 6, sum = 5
        var r = AnnotationTools.PooledRscu(new[] { "CTTCTTCTG", "CTTTTA" }).Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(r["CTT"], Is.EqualTo(6 * 3 / 5.0).Within(1e-12));
            Assert.That(r["CTG"], Is.EqualTo(6 * 1 / 5.0).Within(1e-12));
            Assert.That(r["TTA"], Is.EqualTo(6 * 1 / 5.0).Within(1e-12));
            Assert.That(r["CTC"], Is.EqualTo(0.0));
            Assert.That(r.ContainsKey("TAA"), Is.False);
        });
    }
}
