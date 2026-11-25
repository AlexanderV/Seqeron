namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for ThermodynamicProbeDesigner.DesignProbes, DesignTilingProbes, and probe scoring.
/// Test Unit: PROBE-DESIGN-001
/// </summary>
[TestFixture]
public class ProbeDesigner_ProbeDesign_Tests
{

    // Good sequence for microarray probes (moderate GC, no extreme features)
    private const string MicroarrayTargetSequence =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    // Longer sequence for various tests
    private static readonly string LongSequence =
        new string('A', 30) + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + new string('T', 30);

    private static int IndependentCoverage(int length, int[] starts, int probeLength)
    {
        var covered = new bool[length];
        foreach (int s in starts)
            for (int i = s; i < s + probeLength; i++)
                covered[i] = true;
        return covered.Count(c => c);
    }



    private const string BeaconTarget = "ATGCGTACGTTAGCCGATCGATCGGCTAGCTAGGATCCGATCGTAGCTAGCATCGACTGAC";

    [TestCase(null)]
    [TestCase("")]
    public void CalculateExtinctionCoefficient_NullOrEmpty_ReturnsZero(string? sequence)
    {
        // A7-4: null was NullReferenceException; the siblings CalculateExtinctionCoefficientNearestNeighbor and
        // CalculateMolecularWeight return 0 for null / empty, and so does the mononucleotide sum of no bases.
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficient(sequence!), Is.EqualTo(0));
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor(sequence!), Is.EqualTo(0));
    }

    [Test]
    public void CalculateExtinctionCoefficient_NonDegenerate_Unchanged()
    {
        // Mononucleotide sum: A 15400 + C 7400 + G 11500 + T 8700 + U 9900 + N 10000.
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficient("acgtUN"), Is.EqualTo(62900));
    }



    private const string GuardTarget = "ATGCGTACGTTAGCCGATCGATCGGCTAGCTAGGATCCGATCGTAGCTAGCATCGACTGACTTAGGCA";

}
