namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for RNA secondary structure prediction.
/// Verifies energy, dot-bracket, and stem-loop invariants.
///
/// Test Units: RNA-STRUCT-001, RNA-STEMLOOP-001, RNA-ENERGY-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class RnaStructureProperties
{
    /// <summary>
    /// Minimum free energy is ≤ 0 (folding is energetically favorable or neutral).
    /// </summary>
    [Test]
    [Category("Property")]
    public void MinimumFreeEnergy_IsNonPositive()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);
        Assert.That(mfe, Is.LessThanOrEqualTo(0.001),
            $"MFE={mfe} should be ≤ 0");
    }

    /// <summary>
    /// Predicted structure dot-bracket length equals sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictStructure_DotBracketLength_EqualsSequenceLength()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var structure = RnaSecondaryStructure.PredictStructure(rna);
        Assert.That(structure.DotBracket.Length, Is.EqualTo(rna.Length),
            $"DotBracket length={structure.DotBracket.Length}, sequence length={rna.Length}");
    }

    /// <summary>
    /// Dot-bracket notation only contains valid characters: '.', '(', ')'.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictStructure_DotBracket_HasValidCharacters()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var structure = RnaSecondaryStructure.PredictStructure(rna);
        Assert.That(structure.DotBracket.All(c => c == '.' || c == '(' || c == ')'), Is.True,
            $"Invalid characters in dot-bracket: {structure.DotBracket}");
    }

    /// <summary>
    /// Dot-bracket notation has balanced parentheses.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictStructure_DotBracket_IsBalanced()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var structure = RnaSecondaryStructure.PredictStructure(rna);
        int count = 0;
        foreach (char c in structure.DotBracket)
        {
            if (c == '(') count++;
            else if (c == ')') count--;
            Assert.That(count, Is.GreaterThanOrEqualTo(0), "Unbalanced ')' in dot-bracket");
        }
        Assert.That(count, Is.EqualTo(0), "Unbalanced '(' in dot-bracket");
    }

    /// <summary>
    /// Stem-loop total free energy is finite.
    /// </summary>
    [Test]
    [Category("Property")]
    public void StemLoop_Energy_IsFinite()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var stemLoops = RnaSecondaryStructure.FindStemLoops(rna).ToList();

        foreach (var sl in stemLoops)
            Assert.That(double.IsFinite(sl.TotalFreeEnergy), Is.True,
                $"Stem-loop energy={sl.TotalFreeEnergy} should be finite");
    }

    /// <summary>
    /// Stem-loop positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void StemLoop_Positions_WithinBounds()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var stemLoops = RnaSecondaryStructure.FindStemLoops(rna).ToList();

        foreach (var sl in stemLoops)
        {
            Assert.That(sl.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(sl.End, Is.LessThanOrEqualTo(rna.Length));
        }
    }

    /// <summary>
    /// Base pairs in predicted structure are valid (A-U, G-C, or G-U wobble).
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictStructure_BasePairs_AreValid()
    {
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        var structure = RnaSecondaryStructure.PredictStructure(rna);

        foreach (var bp in structure.BasePairs)
        {
            Assert.That(RnaSecondaryStructure.CanPair(bp.Base1, bp.Base2), Is.True,
                $"Invalid base pair: {bp.Base1}-{bp.Base2} at positions {bp.Position1}-{bp.Position2}");
        }
    }

    /// <summary>
    /// Poly-A sequence has no stem-loops (no pairing partners).
    /// </summary>
    [Test]
    [Category("Property")]
    public void PolyA_HasNoStemLoops()
    {
        string rna = "AAAAAAAAAAAAAAAAAAA";
        var stemLoops = RnaSecondaryStructure.FindStemLoops(rna).ToList();
        Assert.That(stemLoops, Is.Empty,
            "Poly-A should have no stem-loops");
    }

    /// <summary>
    /// ValidateDotBracket returns true for balanced and false for unbalanced.
    /// </summary>
    [TestCase("(((...)))", true)]
    [TestCase("....", true)]
    [TestCase("((..)", false)]
    [TestCase("(..))", false)]
    [Category("Property")]
    public void ValidateDotBracket_CorrectResult(string dotBracket, bool expected)
    {
        Assert.That(RnaSecondaryStructure.ValidateDotBracket(dotBracket), Is.EqualTo(expected));
    }
}
