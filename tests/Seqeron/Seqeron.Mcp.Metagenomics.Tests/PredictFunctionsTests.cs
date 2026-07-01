using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.PredictFunctions (motif containment + BLAST bit-score/E-value,
// best hit = lowest E-value). Reference values from Seqeron.Genomics.Tests
// MetagenomicsAnalyzer_PredictFunctions_Tests (BLOSUM62 lambda=0.3176, K=0.134, W=11).
[TestFixture]
public class PredictFunctionsTests
{
    [Test]
    public void PredictFunctions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.PredictFunctions(
            new[] { new ProteinInput("g1", "WWW") },
            new[] { new FunctionDatabaseEntry("WWW", "tryptophanase", "Amino acid metabolism", "K01667") }));

        // Empty protein set is defined (no annotations), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.PredictFunctions(
            Array.Empty<ProteinInput>(),
            new[] { new FunctionDatabaseEntry("WWW", "f", "p", "k") }));
    }

    [Test]
    public void PredictFunctions_Binding_InvokesSuccessfully()
    {
        // Exact self-match "WWW": BitScore = 18.0202932787533, EValue = 3.3852730346546e-5.
        var exact = MetagenomicsTools.PredictFunctions(
            new[] { new ProteinInput("g1", "WWW") },
            new[] { new FunctionDatabaseEntry("WWW", "tryptophanase", "Amino acid metabolism", "K01667") })
            .Items.Single();

        Assert.Multiple(() =>
        {
            Assert.That(exact.Function, Is.EqualTo("tryptophanase"));
            Assert.That(exact.Pathway, Is.EqualTo("Amino acid metabolism"));
            Assert.That(exact.KoNumber, Is.EqualTo("K01667"));
            Assert.That(exact.BitScore, Is.EqualTo(18.0202932787533).Within(1e-12),
                "BLAST bit score for raw score S = 33 (WWW, BLOSUM62 W = 11).");
            Assert.That(exact.EValue, Is.EqualTo(3.3852730346546e-5).Within(1e-18),
                "E = K*m*n*e^(-lambda*S), S=33, m=n=3.");
        });

        // Best-hit selection: "AAAAWW" contains AAAA (S=16) and WW (S=22); WW has the lower E-value.
        var best = MetagenomicsTools.PredictFunctions(
            new[] { new ProteinInput("g1", "AAAAWW") },
            new[]
            {
                new FunctionDatabaseEntry("AAAA", "weakHit", "PathA", "K00001"),
                new FunctionDatabaseEntry("WW", "strongHit", "PathB", "K00002"),
            })
            .Items.Single();

        Assert.Multiple(() =>
        {
            Assert.That(best.Function, Is.EqualTo("strongHit"),
                "Best hit is the higher-scoring 'WW' (lower E-value).");
            Assert.That(best.EValue, Is.EqualTo(0.0014851955539388528).Within(1e-15));
        });
    }
}
