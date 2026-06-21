namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the RnaStructure area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("RnaStructure")]
public class RnaStructureCombinatorialTests
{
    /// <summary>A guaranteed-foldable hairpin of the requested length: G-run stem, AAAA loop, C-run stem.</summary>
    private static string Hairpin(int length)
    {
        int arm = (length - 4) / 2;
        int extra = (length - 4) - 2 * arm; // pad the loop if length is odd
        return new string('G', arm) + new string('A', 4 + extra) + new string('C', arm);
    }

    private static bool IsBalanced(string dotBracket)
    {
        int depth = 0;
        foreach (char c in dotBracket)
        {
            if (c == '(') depth++;
            else if (c == ')') { depth--; if (depth < 0) return false; }
        }
        return depth == 0;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-STRUCT-001 — RNA secondary-structure prediction (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 71.
    // Spec: tests/TestSpecs/RNA-STRUCT-001.md (canonical PredictStructure /
    //       CalculateMinimumFreeEnergy).
    // Dimensions: algorithm(2: greedy/MFE) × seqLen(3) × temperature(3). Grid 2×3×3 = 18.
    //
    // Model (Nussinov 1978; Zuker 1981): secondary-structure prediction yields a set of nested
    // Watson-Crick/wobble base pairs minimising free energy (≤ 0). Structure prediction here uses
    // fixed Turner-2004 (37 °C) parameters, so it is temperature-independent; temperature enters
    // only the Boltzmann stability of a structure within its ensemble (CalculateStructureProbability).
    //
    // The combinatorial point: algorithm, length and temperature interact — a foldable hairpin
    // yields a valid balanced structure with non-positive energy under both algorithms at every
    // length, and the structure's Boltzmann probability is a proper temperature-dependent value.
    // ═══════════════════════════════════════════════════════════════════════

    public enum FoldAlgorithm { GreedyStemLoop, MinimumFreeEnergy }

    [Test, Combinatorial]
    public void RnaStruct_ValidStructureAndStability_AcrossAlgorithmLengthTemperature(
        [Values(FoldAlgorithm.GreedyStemLoop, FoldAlgorithm.MinimumFreeEnergy)] FoldAlgorithm algorithm,
        [Values(20, 40, 60)] int seqLen,
        [Values(283.15, 310.15, 337.15)] double temperatureK)
    {
        string rna = Hairpin(seqLen);

        double energy;
        if (algorithm == FoldAlgorithm.GreedyStemLoop)
        {
            var s = RnaSecondaryStructure.PredictStructure(rna);
            s.DotBracket.Length.Should().Be(seqLen, "dot-bracket spans the sequence");
            IsBalanced(s.DotBracket).Should().BeTrue("brackets are balanced");
            s.BasePairs.Should().OnlyContain(bp => RnaSecondaryStructure.CanPair(bp.Base1, bp.Base2),
                "every reported pair is a valid base pair");
            energy = s.MinimumFreeEnergy;
        }
        else
        {
            energy = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);
        }

        energy.Should().BeLessThanOrEqualTo(0, "folding does not raise free energy");

        // Temperature enters the Boltzmann stability relative to a 1 kcal/mol more-stable ensemble.
        double p = RnaSecondaryStructure.CalculateStructureProbability(energy, energy - 1.0, temperatureK);
        p.Should().BeInRange(0.0, 1.0 + 1e-9, "a structure probability is in (0,1]");
    }

    /// <summary>
    /// Interaction witness: a foldable hairpin forms base pairs (negative energy) while a
    /// non-pairing poly-A RNA has no pairs and zero free energy.
    /// </summary>
    [Test]
    public void RnaStruct_HairpinFolds_PolyADoesNot()
    {
        var hairpin = RnaSecondaryStructure.PredictStructure(Hairpin(30));
        hairpin.BasePairs.Should().NotBeEmpty("a G/C hairpin folds");
        hairpin.MinimumFreeEnergy.Should().BeLessThan(0);

        var polyA = RnaSecondaryStructure.PredictStructure(new string('A', 30));
        polyA.BasePairs.Should().BeEmpty("poly-A cannot pair");
        polyA.MinimumFreeEnergy.Should().Be(0);
    }

    /// <summary>
    /// Interaction witness: for a fixed structure below its ensemble, the Boltzmann structure
    /// probability rises with temperature (the free-energy gap is discounted by RT).
    /// </summary>
    [Test]
    public void RnaStruct_StructureProbability_RisesWithTemperature()
    {
        double pCold = RnaSecondaryStructure.CalculateStructureProbability(-10, -11, 283.15);
        double pWarm = RnaSecondaryStructure.CalculateStructureProbability(-10, -11, 337.15);
        pWarm.Should().BeGreaterThan(pCold, "a higher RT shrinks the energy gap's penalty");
    }
}
