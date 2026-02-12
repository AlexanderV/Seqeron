using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Residual tests for DisorderPredictor methods not yet covered by a canonical Test Unit.
/// MoRF and LowComplexity tests remain here until their own Test Units are processed.
///
/// Region detection and classification tests → DisorderPredictor_DisorderedRegion_Tests.cs (DISORDER-REGION-001)
/// Per-residue prediction tests → DisorderPredictor_DisorderPrediction_Tests.cs (DISORDER-PRED-001)
/// </summary>
[TestFixture]
public class DisorderPredictorTests
{
    #region MoRF Prediction Tests

    [Test]
    public void PredictMoRFs_OrderedProtein_FewMoRFs()
    {
        string ordered = "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL";

        var morfs = DisorderPredictor.PredictMoRFs(ordered).ToList();

        // Ordered proteins shouldn't have many MoRFs
        Assert.That(morfs.Count, Is.LessThanOrEqualTo(2));
    }

    [Test]
    public void PredictMoRFs_ReturnsValidCoordinates()
    {
        string sequence = "LLLLLLPPPPPPPPPPPPPPPPPPLLLLLL";

        var morfs = DisorderPredictor.PredictMoRFs(sequence).ToList();

        foreach (var morf in morfs)
        {
            Assert.That(morf.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(morf.End, Is.LessThan(sequence.Length));
            Assert.That(morf.Score, Is.InRange(0.0, 1.0));
        }
    }

    #endregion

    #region Low Complexity Region Tests

    [Test]
    public void PredictLowComplexityRegions_PolyQ_FindsRegion()
    {
        string polyQ = "AAAAAQQQQQQQQQQQQQQQQQQAAAAAA";

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            polyQ, windowSize: 10, minLength: 8).ToList();

        Assert.That(regions.Any(r => r.Type.Contains("Q")), Is.True);
    }

    [Test]
    public void PredictLowComplexityRegions_ComplexSequence_NoRegions()
    {
        string complex = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            complex, minLength: 10).ToList();

        Assert.That(regions.Count, Is.LessThanOrEqualTo(1));
    }

    [Test]
    public void PredictLowComplexityRegions_ValidBoundaries()
    {
        string sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        foreach (var region in regions)
        {
            Assert.That(region.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(region.End, Is.LessThan(sequence.Length));
        }
    }

    #endregion
}
