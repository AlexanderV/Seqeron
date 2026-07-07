using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.DifferentialAbundance (Welch's t-test).
// The reported FoldChange field is the log2 fold-change log2(mean2/mean1); a taxon is
// Significant when p < threshold AND |log2FC| > 1. Reference values hand-derived from the
// algorithm: zero-variance groups with differing means give an exact p = 0.
[TestFixture]
public class DifferentialAbundanceTests
{
    private static AbundanceSample Sample(params (string Name, double Fraction)[] items)
        => new(items.Select(t => new AbundanceItem(t.Name, t.Fraction)).ToArray());

    [Test]
    public void DifferentialAbundance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.DifferentialAbundance(
            new[] { Sample(("T", 1)), Sample(("T", 1)) },
            new[] { Sample(("T", 4)), Sample(("T", 4)) }));

        // Empty condition -> defined (no items yielded), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.DifferentialAbundance(
            Array.Empty<AbundanceSample>(),
            new[] { Sample(("T", 4)) }));
    }

    [Test]
    public void DifferentialAbundance_Binding_InvokesSuccessfully()
    {
        // Taxon T: condition1 mean = 1, condition2 mean = 4 across 3 zero-variance replicates.
        // foldChange = 4; log2FC = 2; Welch p = 0 (differing means, zero variance) -> significant.
        var result = MetagenomicsTools.DifferentialAbundance(
            new[] { Sample(("T", 1)), Sample(("T", 1)), Sample(("T", 1)) },
            new[] { Sample(("T", 4)), Sample(("T", 4)), Sample(("T", 4)) },
            pValueThreshold: 0.05);

        var t = result.Items.Single(i => i.Taxon == "T");
        Assert.Multiple(() =>
        {
            Assert.That(t.FoldChange, Is.EqualTo(2.0).Within(1e-10),
                "log2(mean2/mean1) = log2(4/1) = 2.");
            Assert.That(t.PValue, Is.EqualTo(0.0).Within(1e-12),
                "Zero within-group variance with differing means -> Welch p = 0.");
            Assert.That(t.Significant, Is.True, "p < 0.05 AND |log2FC| = 2 > 1 -> significant.");
        });

        // Identical samples: log2FC = 0, p = 1.0 -> not significant.
        var same = MetagenomicsTools.DifferentialAbundance(
            new[] { Sample(("T", 2)), Sample(("T", 2)) },
            new[] { Sample(("T", 2)), Sample(("T", 2)) });
        var st = same.Items.Single(i => i.Taxon == "T");
        Assert.Multiple(() =>
        {
            Assert.That(st.FoldChange, Is.EqualTo(0.0).Within(1e-10), "Equal means -> log2FC 0.");
            Assert.That(st.PValue, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(st.Significant, Is.False);
        });
    }
}
