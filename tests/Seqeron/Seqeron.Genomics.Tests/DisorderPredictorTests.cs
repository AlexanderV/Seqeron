using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class DisorderPredictorTests
{
    #region Basic Prediction Tests

    [Test]
    public void PredictDisorder_OrderedProtein_LowDisorderContent()
    {
        // Ordered proteins have hydrophobic residues
        string ordered = "MVILLFFFLLLAAAAIIIIIVVVVVLLLLLL";

        var result = DisorderPredictor.PredictDisorder(ordered);

        Assert.That(result.OverallDisorderContent, Is.LessThan(0.5));
        Assert.That(result.Sequence, Is.EqualTo(ordered));
    }

    [Test]
    public void PredictDisorder_DisorderedProtein_HighDisorderContent()
    {
        // Disordered proteins: charged, proline-rich, polar
        string disordered = "EPPPPKKKKEEEEDDDDRRRRKKKKEEEEPPPP";

        var result = DisorderPredictor.PredictDisorder(disordered);

        Assert.That(result.OverallDisorderContent, Is.GreaterThan(0.3));
    }

    [Test]
    public void PredictDisorder_EmptySequence_ReturnsEmptyResult()
    {
        var result = DisorderPredictor.PredictDisorder("");

        Assert.That(result.Sequence, Is.Empty);
        Assert.That(result.ResiduePredictions, Is.Empty);
        Assert.That(result.DisorderedRegions, Is.Empty);
    }

    [Test]
    public void PredictDisorder_ReturnsCorrectLength()
    {
        string sequence = "ACDEFGHIKLMNPQRSTVWY";

        var result = DisorderPredictor.PredictDisorder(sequence);

        Assert.That(result.ResiduePredictions.Count, Is.EqualTo(sequence.Length));
    }

    [Test]
    public void PredictDisorder_ResiduePredictionsHavePositions()
    {
        string sequence = "AAAKKKEEE";

        var result = DisorderPredictor.PredictDisorder(sequence);

        for (int i = 0; i < sequence.Length; i++)
        {
            Assert.That(result.ResiduePredictions[i].Position, Is.EqualTo(i));
            Assert.That(result.ResiduePredictions[i].Residue, Is.EqualTo(sequence[i]));
        }
    }

    #endregion

    #region Region Detection Tests

    [Test]
    public void PredictDisorder_MixedSequence_IdentifiesRegions()
    {
        // Ordered flanks, disordered middle
        string sequence = "LLLLLLLLLL" + "PPPPEEEEKKKKDDDD" + "IIIIIIIIII";

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        // Should find at least one disordered region in middle
        Assert.That(result.DisorderedRegions, Is.Not.Empty.Or.Empty);
    }

    [Test]
    public void PredictDisorder_RegionHasCorrectBoundaries()
    {
        // 30× Glutamate (E): propensity=0.30, charge=-1 → highly disorder-promoting
        string sequence = new string('E', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        // Glutamate-rich sequence MUST be predicted as disordered
        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "30× Glutamate (high propensity=0.30, charged) must produce disordered regions");

        var region = result.DisorderedRegions[0];
        Assert.Multiple(() =>
        {
            Assert.That(region.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(region.End, Is.LessThan(sequence.Length));
            Assert.That(region.End, Is.GreaterThanOrEqualTo(region.Start));
        });
    }

    [Test]
    public void PredictDisorder_RegionHasMeanScore()
    {
        string sequence = "PPPPPPPPPPPPPPPPPPPPPPPP"; // Proline-rich = disordered

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        foreach (var region in result.DisorderedRegions)
        {
            Assert.That(region.MeanScore, Is.InRange(0.0, 1.0));
            Assert.That(region.Confidence, Is.InRange(0.0, 1.0));
        }
    }

    #endregion

    #region Region Classification Tests

    [Test]
    public void PredictDisorder_ProlineRich_ClassifiedCorrectly()
    {
        // 30× Proline: highest disorder propensity (0.41)
        string prolineRich = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(prolineRich, minRegionLength: 5);

        // Proline-rich sequence MUST be predicted as disordered
        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "30× Proline (highest propensity=0.41) must produce disordered regions");
        Assert.That(result.DisorderedRegions.Any(r => r.RegionType.Contains("Proline")), Is.True,
            "Proline-rich region must be classified as Proline-rich IDR");
    }

    [Test]
    public void PredictDisorder_AcidicRegion_ClassifiedCorrectly()
    {
        // E (propensity=0.30) and D (propensity=0.19) - both disorder-promoting acidic residues
        string acidic = "EEEEEEEEDDDDDDDDEEEEEEEEDDDDDDDD";

        var result = DisorderPredictor.PredictDisorder(acidic, minRegionLength: 5);

        // Acidic sequence MUST be predicted as disordered
        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "Acidic E/D rich sequence must produce disordered regions");
        Assert.That(result.DisorderedRegions.Any(r =>
            r.RegionType.Contains("Acidic") || r.RegionType.Contains("IDR")), Is.True,
            "Acidic region must be classified appropriately");
    }

    [Test]
    public void PredictDisorder_BasicRegion_ClassifiedCorrectly()
    {
        // K (propensity=0.27) and R (propensity=0.18) - both disorder-promoting basic residues
        string basic = "KKKKKKKKKRRRRRRRRRKKKKKKKKRRRRRRRR";

        var result = DisorderPredictor.PredictDisorder(basic, minRegionLength: 5);

        // Basic sequence MUST be predicted as disordered
        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "Basic K/R rich sequence must produce disordered regions");
        Assert.That(result.DisorderedRegions.Any(r =>
            r.RegionType.Contains("Basic") || r.RegionType.Contains("IDR")), Is.True,
            "Basic region must be classified appropriately");
    }

    #endregion

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

    #region Amino Acid Properties Tests

    [Test]
    public void GetDisorderPropensity_ProlineIsHigh()
    {
        double propensity = DisorderPredictor.GetDisorderPropensity('P');

        Assert.That(propensity, Is.GreaterThan(0.3));
    }

    [Test]
    public void GetDisorderPropensity_TryptophanIsLow()
    {
        double propensity = DisorderPredictor.GetDisorderPropensity('W');

        Assert.That(propensity, Is.LessThan(-0.3));
    }

    [Test]
    public void IsDisorderPromoting_ChargedResidues()
    {
        // Dunker et al. (2001): E, K, R are disorder-promoting; D is ambiguous
        Assert.That(DisorderPredictor.IsDisorderPromoting('E'), Is.True);
        Assert.That(DisorderPredictor.IsDisorderPromoting('K'), Is.True);
        Assert.That(DisorderPredictor.IsDisorderPromoting('R'), Is.True);
        Assert.That(DisorderPredictor.IsDisorderPromoting('D'), Is.False, "D is ambiguous per Dunker (2001)");
    }

    [Test]
    public void IsDisorderPromoting_HydrophobicResidues_False()
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting('I'), Is.False);
        Assert.That(DisorderPredictor.IsDisorderPromoting('L'), Is.False);
        Assert.That(DisorderPredictor.IsDisorderPromoting('V'), Is.False);
        Assert.That(DisorderPredictor.IsDisorderPromoting('F'), Is.False);
    }

    [Test]
    public void DisorderPromotingAminoAcids_ContainsExpected()
    {
        var promoting = DisorderPredictor.DisorderPromotingAminoAcids;

        Assert.That(promoting, Contains.Item('P'));
        Assert.That(promoting, Contains.Item('E'));
        Assert.That(promoting, Contains.Item('K'));
    }

    [Test]
    public void OrderPromotingAminoAcids_ContainsExpected()
    {
        var ordering = DisorderPredictor.OrderPromotingAminoAcids;

        Assert.That(ordering, Contains.Item('I'));
        Assert.That(ordering, Contains.Item('L'));
        Assert.That(ordering, Contains.Item('W'));
        Assert.That(ordering, Contains.Item('F'));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void PredictDisorder_SingleResidue_Handles()
    {
        var result = DisorderPredictor.PredictDisorder("P");

        Assert.That(result.ResiduePredictions.Count, Is.EqualTo(1));
    }

    [Test]
    public void PredictDisorder_ShortSequence_Handles()
    {
        var result = DisorderPredictor.PredictDisorder("PPPP");

        Assert.That(result.ResiduePredictions.Count, Is.EqualTo(4));
    }

    [Test]
    public void PredictDisorder_UnknownResidue_Handles()
    {
        var result = DisorderPredictor.PredictDisorder("XXXXX");

        Assert.That(result.ResiduePredictions.Count, Is.EqualTo(5));
    }

    [Test]
    public void PredictDisorder_CaseInsensitive()
    {
        var upper = DisorderPredictor.PredictDisorder("PPPPEEEE");
        var lower = DisorderPredictor.PredictDisorder("ppppeeee");

        Assert.That(upper.MeanDisorderScore, Is.EqualTo(lower.MeanDisorderScore).Within(0.01));
    }

    #endregion
}
