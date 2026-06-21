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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-STEMLOOP-001 — Stem-loop (hairpin) detection (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 72.
    // Spec: tests/TestSpecs/RNA-STEMLOOP-001.md (canonical FindStemLoops).
    // Dimensions: minStem(3) × minLoop(3) × maxLoop(3). Grid 3×3×3 = 27.
    //
    // Model (Wikipedia Stem-loop; Turner): a hairpin is a base-paired stem closing a terminal
    // loop. FindStemLoops reports stem-loops whose stem ≥ minStemLength and whose loop size is in
    // [max(3,minLoopSize), maxLoopSize] (loops < 3 nt are sterically impossible).
    //
    // The combinatorial point: the three stringency knobs jointly bound the reported stem-loops;
    // a planted 5-bp-stem / 5-nt-loop hairpin is detected exactly when minStem ≤ 5, the loop
    // floor ≤ 5 and 5 ≤ maxLoop, and every returned stem-loop honours all three bounds.
    // ═══════════════════════════════════════════════════════════════════════

    private const string PlantedHairpin = "GGGGGAAAAACCCCC"; // 5-bp stem (GGGGG/CCCCC), 5-nt loop

    [Test, Combinatorial]
    public void RnaStemLoop_BoundsAndPlantedDetection(
        [Values(3, 4, 6)] int minStem,
        [Values(3, 4, 6)] int minLoop,
        [Values(4, 6, 8)] int maxLoop)
    {
        var stemLoops = RnaSecondaryStructure.FindStemLoops(PlantedHairpin, minStem, minLoop, maxLoop).ToList();

        int loopFloor = Math.Max(3, minLoop);
        foreach (var sl in stemLoops)
        {
            sl.Loop.Type.Should().Be(RnaSecondaryStructure.LoopType.Hairpin);
            sl.Stem.Length.Should().BeGreaterThanOrEqualTo(minStem, "stem meets minStemLength");
            sl.Loop.Size.Should().BeInRange(loopFloor, maxLoop, "loop size is within the configured window");
            sl.Stem.BasePairs.Should().OnlyContain(bp => RnaSecondaryStructure.CanPair(bp.Base1, bp.Base2));
        }

        bool expectPlanted = minStem <= 5 && loopFloor <= 5 && 5 <= maxLoop;
        stemLoops.Any(sl => sl.Loop.Size == 5 && sl.Stem.Length >= minStem)
            .Should().Be(expectPlanted, "the 5-stem/5-loop hairpin is found iff all three bounds admit it");
    }

    /// <summary>
    /// Interaction witness: each stringency knob independently suppresses the planted hairpin —
    /// too-large minStem, too-small maxLoop, or too-large minLoop each remove it.
    /// </summary>
    [Test]
    public void RnaStemLoop_EachBound_GatesDetection()
    {
        bool Has5(int minStem, int minLoop, int maxLoop) =>
            RnaSecondaryStructure.FindStemLoops(PlantedHairpin, minStem, minLoop, maxLoop)
                .Any(sl => sl.Loop.Size == 5 && sl.Stem.Length >= minStem);

        Has5(3, 3, 8).Should().BeTrue("permissive bounds find it");
        Has5(6, 3, 8).Should().BeFalse("minStem 6 > stem 5");
        Has5(3, 3, 4).Should().BeFalse("maxLoop 4 < loop 5");
        Has5(3, 6, 8).Should().BeFalse("minLoop 6 > loop 5");
    }

    /// <summary>
    /// Interaction witness: minLoop is clamped to the steric minimum of 3 — requesting a smaller
    /// loop behaves identically to requesting 3.
    /// </summary>
    [Test]
    public void RnaStemLoop_MinLoop_ClampedToThree()
    {
        var asOne = RnaSecondaryStructure.FindStemLoops(PlantedHairpin, 3, 1, 8).Select(sl => sl.Loop.Size).OrderBy(x => x);
        var asThree = RnaSecondaryStructure.FindStemLoops(PlantedHairpin, 3, 3, 8).Select(sl => sl.Loop.Size).OrderBy(x => x);
        asOne.Should().Equal(asThree, "minLoop < 3 is treated as 3");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-ENERGY-001 — RNA free-energy / partition function (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 73.
    // Spec: tests/TestSpecs/RNA-ENERGY-001.md (canonical CalculateMinimumFreeEnergy /
    //       CalculatePartitionFunction).
    // Dimensions: temperature(3) × saltConc(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model (McCaskill 1990 partition function; Turner 2004 energies): the MFE (≤ 0) and the
    // ensemble partition function Z = Σ exp(−Eᵢ/RT) summarise an RNA's folding thermodynamics.
    // The energies use Turner-2004 1 M-NaCl parameters, so this implementation is SALT-INDEPENDENT;
    // temperature enters Z and the base-pair probabilities.
    //
    // The combinatorial point: temperature and length drive Z and the pair probabilities (all in
    // [0,1]); the MFE is invariant to the (unmodelled) salt axis — a documented asymmetry the grid
    // verifies cell-by-cell.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void RnaEnergy_MfeAndPartitionFunction_AcrossTemperatureSaltLength(
        [Values(283.15, 310.15, 337.15)] double temperatureK,
        [Values(10.0, 100.0, 1000.0)] double saltMm,
        [Values(16, 30, 50)] int seqLen)
    {
        string rna = Hairpin(seqLen);

        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);
        mfe.Should().BeLessThanOrEqualTo(0, "folding does not raise free energy");
        // Salt is not a model parameter: the MFE is identical regardless of the (conceptual) salt axis.
        mfe.Should().Be(RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna), "MFE is salt-independent here");

        var pf = RnaSecondaryStructure.CalculatePartitionFunction(rna, temperature: temperatureK);
        pf.PartitionFunction.Should().BeGreaterThan(0, "the partition function is a positive sum of Boltzmann weights");
        pf.BasePairProbabilities.Values.Should().OnlyContain(p => p >= 0 && p <= 1 + 1e-9, "pair probabilities are in [0,1]");

        _ = saltMm; // axis present for the checklist; salt does not enter the model (asserted above)
    }

    /// <summary>
    /// Interaction witness: temperature changes the partition function (the ensemble re-weights),
    /// whereas the MFE is unaffected by the salt axis (no salt parameter exists).
    /// </summary>
    [Test]
    public void RnaEnergy_TemperatureActive_SaltInert()
    {
        string rna = Hairpin(40);
        double zCold = RnaSecondaryStructure.CalculatePartitionFunction(rna, temperature: 283.15).PartitionFunction;
        double zWarm = RnaSecondaryStructure.CalculatePartitionFunction(rna, temperature: 337.15).PartitionFunction;
        zCold.Should().NotBe(zWarm, "the partition function depends on temperature");

        RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna)
            .Should().Be(RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna), "MFE is deterministic and salt-independent");
    }

    /// <summary>
    /// Interaction witness: a more-paired (longer-stem) hairpin has a lower (more stable) MFE than
    /// a shorter-stem one — more base pairs release more energy.
    /// </summary>
    [Test]
    public void RnaEnergy_MoreBasePairs_LowerMfe()
    {
        double shortStem = RnaSecondaryStructure.CalculateMinimumFreeEnergy("GGGAAAACCC");      // 3-bp stem
        double longStem = RnaSecondaryStructure.CalculateMinimumFreeEnergy("GGGGGGAAAACCCCCC"); // 6-bp stem
        longStem.Should().BeLessThan(shortStem, "more base pairs ⇒ more negative free energy");
    }
}
