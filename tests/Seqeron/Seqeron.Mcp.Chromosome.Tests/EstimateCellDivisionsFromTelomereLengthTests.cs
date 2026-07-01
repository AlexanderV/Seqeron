using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>estimate_cell_divisions_from_telomere_length</c>.
/// ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength = max(0, (birth - current) / loss).
/// </summary>
[TestFixture]
public class EstimateCellDivisionsFromTelomereLengthTests
{
    [Test]
    public void EstimateCellDivisions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.EstimateCellDivisionsFromTelomereLength(10000));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCellDivisionsFromTelomereLength(-1));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCellDivisionsFromTelomereLength(10000, birthLength: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCellDivisionsFromTelomereLength(10000, lossPerDivision: 0));
    }

    [Test]
    public void EstimateCellDivisions_Binding_InvokesSuccessfully()
    {
        // (15000 - 10000) / 50 = 100
        Assert.That(ChromosomeTools.EstimateCellDivisionsFromTelomereLength(10000).CellDivisions, Is.EqualTo(100.0).Within(1e-9));
        // current >= birth -> lost <= 0 -> clamped to 0
        Assert.That(ChromosomeTools.EstimateCellDivisionsFromTelomereLength(15000).CellDivisions, Is.EqualTo(0.0).Within(1e-9));
        Assert.That(ChromosomeTools.EstimateCellDivisionsFromTelomereLength(20000).CellDivisions, Is.EqualTo(0.0).Within(1e-9));
        // (12000 - 2000) / 100 = 100
        Assert.That(ChromosomeTools.EstimateCellDivisionsFromTelomereLength(2000, birthLength: 12000, lossPerDivision: 100).CellDivisions,
            Is.EqualTo(100.0).Within(1e-9));
    }
}
