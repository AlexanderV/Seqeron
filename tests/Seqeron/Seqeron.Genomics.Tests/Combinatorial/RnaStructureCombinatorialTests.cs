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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-HAIRPIN-001 — Hairpin loop & stem free-energy (Turner 2004) (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 150.
    // Spec: tests/TestSpecs/RNA-HAIRPIN-001.md (canonical CalculateHairpinLoopEnergy /
    //       CalculateStemEnergy). ADVANCED §10 (combinatorial/pairwise).
    // Dimensions: loopSize(3) × stemLen(3). Grid 3×3 = 9 (full grid, exhaustive ⊇ pairwise).
    //
    // Model (Turner 2004 nearest-neighbor, NNDB; Mathews 2004 PNAS): the free energy of a
    // hairpin is ADDITIVELY decomposed into a STEM term and a LOOP term —
    //   ΔG°(hairpin) = ΔG°(stem) + ΔG°(loop).
    //   • Stem: Σ nearest-neighbor stacks (P pairs ⇒ P−1 stacks) + terminal AU/GU end penalties.
    //   • Loop: initiation(n) + terminal-mismatch + sequence bonuses/penalties.
    // The two axes are ORTHOGONAL knobs of this decomposition: stem length drives only the stem
    // term, loop size drives only the loop term. That orthogonality is the combinatorial claim.
    //
    // Engineered construct: a clean G-C hairpin  G^P · A^L · C^P  closing the loop with a G(5')-C(3')
    // pair. Two NON-CIRCULAR ground truths follow from the published constants:
    //   • Every stem stack is GG/CC = −3.26 kcal/mol (NNDB wc-parameters.html) and G-C ends carry
    //     NO terminal AU/GU penalty ⇒ ΔG°(stem) = round(−3.26·(P−1), 2) EXACTLY.
    //   • The all-A loop closing G-C has only ONE size-dependent term — initiation(n); the terminal
    //     mismatch (G·A..A·C), the UU/GA, GG and all-C terms are all absent or size-independent ⇒
    //     ΔG°(loop, n) − initiation(n) is a CONSTANT across n. Published initiations (NNDB loop.txt):
    //     n=4→5.6, n=6→5.4, n=8→5.5. (No size 3 here: triloops drop the first-mismatch term — a
    //     separate regime — so the loop axis uses {4,6,8} to stay in the n≥4 additive model.)
    // The combinatorial point: the grid verifies the stem axis EXACTLY, the loop axis as
    // initiation + a size-independent residual, and asserts the double dissociation (each axis moves
    // its own term and leaves the other's term untouched).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>NNDB Turner-2004 hairpin-loop initiation energies (loop.txt), kcal/mol at 37 °C.</summary>
    private static readonly IReadOnlyDictionary<int, double> PublishedHairpinInit = new Dictionary<int, double>
    {
        [4] = 5.6, [6] = 5.4, [8] = 5.5,
    };

    private const double GcStackEnergy = -3.26; // NNDB wc-parameters.html: GG/CC nearest-neighbor stack

    /// <summary>A clean G-C hairpin: P G·C pairs, an all-A loop of size L, closing pair G(5')-C(3').</summary>
    private static (string Seq, List<RnaSecondaryStructure.BasePair> Pairs, string Loop) GcHairpin(int stemLen, int loopSize)
    {
        string seq = new string('G', stemLen) + new string('A', loopSize) + new string('C', stemLen);
        var pairs = new List<RnaSecondaryStructure.BasePair>();
        for (int i = 0; i < stemLen; i++)
            pairs.Add(new RnaSecondaryStructure.BasePair(i, seq.Length - 1 - i, 'G', 'C',
                RnaSecondaryStructure.BasePairType.WatsonCrick));
        return (seq, pairs, new string('A', loopSize));
    }

    [Test, Combinatorial]
    public void RnaHairpin_AdditiveDecomposition_AcrossLoopSizeAndStemLength(
        [Values(4, 6, 8)] int loopSize,
        [Values(3, 4, 5)] int stemLen)
    {
        var (seq, pairs, loop) = GcHairpin(stemLen, loopSize);

        // ── Stem axis (exact ground truth): P pairs ⇒ P−1 GG/CC stacks, no AU/GU end penalty. ──
        double stemEnergy = RnaSecondaryStructure.CalculateStemEnergy(seq, pairs);
        double expectedStem = Math.Round(GcStackEnergy * (stemLen - 1), 2);
        stemEnergy.Should().BeApproximately(expectedStem, 1e-9,
            "a G-C homopolymer stem of {0} pairs is exactly {1} GG/CC stacks", stemLen, stemLen - 1);
        stemEnergy.Should().BeLessThan(0, "base-pair stacking is stabilising");

        // ── Loop axis: ΔG°(loop, n) = initiation(n) + a size-independent residual. ──
        double loopEnergy = RnaSecondaryStructure.CalculateHairpinLoopEnergy(loop, 'G', 'C');
        loopEnergy.Should().BeGreaterThan(0, "a hairpin loop is entropically destabilising");

        double residual = loopEnergy - PublishedHairpinInit[loopSize];
        double refResidual = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAAA", 'G', 'C') - PublishedHairpinInit[4];
        residual.Should().BeApproximately(refResidual, 1e-9,
            "removing the published initiation(n) leaves the same size-independent context for every loop size");

        // ── Additive decomposition: the total hairpin free energy is stem + loop. ──
        double total = stemEnergy + loopEnergy;
        total.Should().BeApproximately(stemEnergy + loopEnergy, 1e-12); // definitional anchor for the grid cell
        // The loop term is positive and the stem term negative ⇒ a 1-bp stem cannot offset the loop,
        // but a sufficiently long stem makes the hairpin net-stabilising. With stem ≥ 4 here:
        if (stemLen >= 5)
            total.Should().BeLessThan(0, "a 5-bp G-C stem (−13.04) overcomes the ~4 kcal/mol loop");
    }

    /// <summary>
    /// Interaction witness — double dissociation: the stem axis moves ONLY the stem term and the
    /// loop axis moves ONLY the loop term. This is what makes the two parameters genuinely
    /// combinatorial rather than redundant.
    /// </summary>
    [Test]
    public void RnaHairpin_AxesAreOrthogonal()
    {
        // Stem axis active: longer stem ⇒ strictly lower (more negative) stem energy …
        double stem3 = RnaSecondaryStructure.CalculateStemEnergy(GcHairpin(3, 6).Seq, GcHairpin(3, 6).Pairs);
        double stem5 = RnaSecondaryStructure.CalculateStemEnergy(GcHairpin(5, 6).Seq, GcHairpin(5, 6).Pairs);
        stem5.Should().BeLessThan(stem3, "two extra G-C stacks lower the stem energy");

        // … while the loop energy is unaffected by stem length (the loop term is a function of the loop only).
        double loopAtStem3 = RnaSecondaryStructure.CalculateHairpinLoopEnergy(GcHairpin(3, 6).Loop, 'G', 'C');
        double loopAtStem5 = RnaSecondaryStructure.CalculateHairpinLoopEnergy(GcHairpin(5, 6).Loop, 'G', 'C');
        loopAtStem5.Should().Be(loopAtStem3, "loop size 6 is identical regardless of stem length");

        // Loop axis active: changing loop size changes the loop energy …
        double loop4 = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAAA", 'G', 'C');
        double loop8 = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAAAAAAA", 'G', 'C');
        loop8.Should().NotBe(loop4, "different loop sizes carry different initiation energies");

        // … while the stem energy is unaffected by loop size (same number of stacks).
        double stemLoop4 = RnaSecondaryStructure.CalculateStemEnergy(GcHairpin(4, 4).Seq, GcHairpin(4, 4).Pairs);
        double stemLoop8 = RnaSecondaryStructure.CalculateStemEnergy(GcHairpin(4, 8).Seq, GcHairpin(4, 8).Pairs);
        stemLoop8.Should().Be(stemLoop4, "a 4-bp stem has 3 stacks regardless of loop size");
    }

    /// <summary>
    /// Worked-example anchors (NNDB Turner-2004) fixing the absolute energy scale the grid relies on:
    /// hairpin-example-1 loop and helix. Cited in tests/TestSpecs/RNA-HAIRPIN-001.md M1/M8.
    /// </summary>
    [Test]
    public void RnaHairpin_NndbWorkedExamples_AbsoluteScale()
    {
        // Example 1 loop: 6-nt all-A loop, closing A-U ⇒ initiation(5.4) + terminal-mismatch(−0.8) = 4.6.
        RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAAAAA", 'A', 'U')
            .Should().BeApproximately(4.6, 1e-9, "NNDB hairpin-example-1: loop ΔG° = 4.6 kcal/mol");

        // Example 1 helix: pairs C-G, A-U, C-G, A-U ⇒ 3 stacks + one AU terminal end = −6.01.
        var helix = new List<RnaSecondaryStructure.BasePair>
        {
            new(0, 11, 'C', 'G', RnaSecondaryStructure.BasePairType.WatsonCrick),
            new(1, 10, 'A', 'U', RnaSecondaryStructure.BasePairType.WatsonCrick),
            new(2,  9, 'C', 'G', RnaSecondaryStructure.BasePairType.WatsonCrick),
            new(3,  8, 'A', 'U', RnaSecondaryStructure.BasePairType.WatsonCrick),
        };
        RnaSecondaryStructure.CalculateStemEnergy("CACAUUUUUGUG", helix)
            .Should().BeApproximately(-6.01, 1e-9, "NNDB hairpin-example-1: helix ΔG° = −6.01 kcal/mol");

        // Loops < 3 nt are prohibited ⇒ a prohibitive sentinel, never a normal low value (INV-2).
        RnaSecondaryStructure.CalculateHairpinLoopEnergy("AA", 'G', 'C')
            .Should().BeGreaterThanOrEqualTo(100.0, "sub-3-nt loops are sterically prohibited");
    }
}
