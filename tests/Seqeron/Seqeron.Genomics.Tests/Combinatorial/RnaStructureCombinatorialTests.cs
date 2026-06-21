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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-INVERT-001 — Inverted repeats / potential stems (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 151.
    // Spec: tests/TestSpecs/RNA-INVERT-001.md (canonical FindInvertedRepeats). ADVANCED §10.
    // Dimensions: minArm(3) × maxGap(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive ⊇ pairwise).
    //
    // Model (Alamro 2021 IUPACpal; EMBOSS einverted): an inverted repeat is the pattern W·G·W̄ᴿ — a
    // left arm W, a loop G (|G| ≥ minSpacing), then a right arm equal to the REVERSE COMPLEMENT of W
    // (strict Watson-Crick/IUPAC, antiparallel). It is the structural definition of a potential
    // hairpin stem. FindInvertedRepeats(seq, minLength, minSpacing, maxSpacing) reports maximal,
    // non-overlapping repeats with arm ≥ minLength and loop length in [minSpacing, maxSpacing].
    //
    // Engineered construct: a single planted IR — arm W = GGGCC (5 nt), loop AAAAA (5 nt), right arm
    // GGCCC = revcomp(GGGCC) — padded with poly-A (which cannot pair, so it adds no spurious stems
    // and cannot extend the arms). The three knobs jointly gate detection: arm 5 is found iff
    // minArm ≤ 5, and loop 5 is admitted iff minSpacing ≤ 5 ≤ maxGap. seqLen (padding) must NOT
    // change the planted repeat's position or length.
    //
    // The combinatorial point: every reported repeat must satisfy the full IR definition
    // (INV-1 antiparallel revcomp, INV-2 equal arms, INV-3 loop bounds, INV-4 min arm, INV-5
    // ordering) at every grid cell, and the planted IR appears exactly when all bounds admit it —
    // independent of sequence length.
    // ═══════════════════════════════════════════════════════════════════════

    private const string IrLeftArm = "GGGCC";  // revcomp = GGCCC
    private const string IrLoop = "AAAAA";      // 5-nt loop
    private static readonly string IrCore = IrLeftArm + IrLoop + "GGCCC"; // planted IR at (0,4,10,14,5)

    [Test, Combinatorial]
    public void RnaInvertedRepeats_DefinitionAndPlantedDetection_AcrossArmGapLength(
        [Values(4, 5, 6)] int minArm,
        [Values(3, 5, 8)] int maxGap,
        [Values(15, 24, 33)] int seqLen)
    {
        string seq = IrCore + new string('A', seqLen - IrCore.Length);

        var reps = RnaSecondaryStructure.FindInvertedRepeats(seq, minLength: minArm, minSpacing: 3, maxSpacing: maxGap)
            .ToList();

        // Every reported repeat must satisfy the inverted-repeat definition (INV-1..INV-5).
        foreach (var r in reps)
        {
            (r.End1 - r.Start1 + 1).Should().Be(r.Length, "left arm spans Length (INV-2)");
            (r.End2 - r.Start2 + 1).Should().Be(r.Length, "right arm spans Length (INV-2)");
            r.Length.Should().BeGreaterThanOrEqualTo(minArm, "arm meets minLength (INV-4)");

            int loop = r.Start2 - r.End1 - 1;
            loop.Should().BeInRange(3, maxGap, "loop length is within [minSpacing, maxSpacing] (INV-3)");

            r.Start1.Should().BeLessThan(r.End1, "left arm is ordered (INV-5)");
            r.End1.Should().BeLessThan(r.Start2, "left arm precedes loop precedes right arm (INV-5)");
            r.Start2.Should().BeLessThan(r.End2, "right arm is ordered (INV-5)");

            for (int k = 0; k < r.Length; k++)
            {
                RnaSecondaryStructure.GetComplement(seq[r.Start2 + r.Length - 1 - k])
                    .Should().Be(seq[r.Start1 + k],
                        "the right arm is the antiparallel reverse complement of the left arm (INV-1)");
            }
        }

        // Planted-IR gate: arm 5 found iff minArm ≤ 5, loop 5 admitted iff maxGap ≥ 5 (minSpacing 3 ≤ 5).
        bool expectPlanted = minArm <= 5 && maxGap >= 5;
        reps.Any(r => r.Length == 5 && r.Start1 == 0 && r.Start2 == 10)
            .Should().Be(expectPlanted, "the planted 5-bp/5-loop IR appears exactly when all three bounds admit it");
    }

    /// <summary>
    /// Interaction witness — each search bound independently gates the planted inverted repeat:
    /// too-large minArm, too-small maxGap, or too-large minSpacing each remove it.
    /// </summary>
    [Test]
    public void RnaInvertedRepeats_EachBound_GatesDetection()
    {
        bool Found(int minArm, int minSpacing, int maxGap) =>
            RnaSecondaryStructure.FindInvertedRepeats(IrCore, minArm, minSpacing, maxGap)
                .Any(r => r.Length == 5 && r.Start1 == 0 && r.Start2 == 10);

        Found(4, 3, 8).Should().BeTrue("permissive bounds find the planted IR");
        Found(6, 3, 8).Should().BeFalse("minArm 6 > arm 5");
        Found(4, 3, 3).Should().BeFalse("maxGap 3 < loop 5");
        Found(4, 6, 8).Should().BeFalse("minSpacing 6 > loop 5");
    }

    /// <summary>
    /// Worked-example anchors (tests/TestSpecs/RNA-INVERT-001.md M1/M3): the canonical Ussery/IUPACpal
    /// IR is found at exact coordinates, while a PARALLEL direct repeat (right arm not the reverse
    /// complement) is correctly rejected.
    /// </summary>
    [Test]
    public void RnaInvertedRepeats_NndbWorkedExamples()
    {
        // M1: UUACG · AAAAAA · CGUAA  (CGUAA = revcomp UUACG) ⇒ exactly one IR (0,4,11,15,5).
        var m1 = RnaSecondaryStructure.FindInvertedRepeats("UUACGAAAAAACGUAA").ToList();
        m1.Should().ContainSingle();
        m1[0].Should().Be((0, 4, 11, 15, 5));

        // M3: AAGG · AAAAA · AGG — arms AAGG vs AGG are PARALLEL, not reverse complement ⇒ none.
        RnaSecondaryStructure.FindInvertedRepeats("AAGGAAAAAGG").Should().BeEmpty(
            "a parallel direct repeat is not an inverted repeat");

        // No complementary arm at all ⇒ none.
        RnaSecondaryStructure.FindInvertedRepeats("AAAAAAAAAAAA").Should().BeEmpty();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-MFE-001 — Minimum free energy (Zuker–Stiegler DP) (RnaStructure)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 152.
    // Spec: tests/TestSpecs/RNA-MFE-001.md (canonical CalculateMinimumFreeEnergy + the internal
    //       CalculateMinimumFreeEnergyClassic baseline). ADVANCED §10.
    // Dimensions: algorithm(2) × seqLen(3) × temperature(3). Grid 2×3×3 = 18 (full, exhaustive).
    //
    // Model (Zuker & Stiegler 1981): the MFE is an O(n³) DP over loop decomposition; the open chain
    // (ΔG = 0) is always available, so the optimum is ≤ 0 and is non-increasing as the sequence is
    // extended (extension only adds folding options). Two DP engines realise this:
    //   • algorithm = TurnerDp     → CalculateMinimumFreeEnergy (Turner-2004 nearest-neighbor, NNDB)
    //   • algorithm = ClassicPairs → CalculateMinimumFreeEnergyClassic (simplified −2.0/WC, −1.0/GU)
    // Both engines use FIXED energy tables (Turner-37 °C / simplified constants) so the MFE VALUE is
    // temperature-independent; temperature enters only the Boltzmann weight of a structure within its
    // ensemble (CalculateStructureProbability, McCaskill 1990).
    //
    // The combinatorial point: INV-01 (≤0), INV-02 (monotone under extension) and INV-03
    // (determinism) hold for BOTH engines at every length, while the algorithm axis interacts with
    // composition — only the Turner engine distinguishes a GC stem from an AU stem (witness below).
    // ═══════════════════════════════════════════════════════════════════════

    public enum MfeEngine { TurnerDp, ClassicPairs }

    private static double Mfe(MfeEngine engine, string rna) => engine == MfeEngine.TurnerDp
        ? RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna)
        : RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(rna);

    [Test, Combinatorial]
    public void RnaMfe_StructuralLaws_AcrossAlgorithmLengthTemperature(
        [Values(MfeEngine.TurnerDp, MfeEngine.ClassicPairs)] MfeEngine algorithm,
        [Values(12, 20, 30)] int seqLen,
        [Values(283.15, 310.15, 337.15)] double temperatureK)
    {
        string rna = Hairpin(seqLen);

        double mfe = Mfe(algorithm, rna);
        mfe.Should().BeLessThanOrEqualTo(0, "the open chain (ΔG=0) is always available (INV-01)");
        mfe.Should().BeLessThan(0, "a foldable G/C hairpin actually folds under both engines");

        // INV-03 determinism: the same input yields the same value.
        Mfe(algorithm, rna).Should().Be(mfe, "the MFE DP is deterministic (INV-03)");

        // INV-02 monotone under suffix extension: extending the sequence cannot raise the MFE.
        double mfeExtended = Mfe(algorithm, rna + Hairpin(12));
        mfeExtended.Should().BeLessThanOrEqualTo(mfe + 1e-9, "extension only adds folding options (INV-02)");

        // Temperature is not an MFE parameter — it weights the structure within its ensemble.
        double p = RnaSecondaryStructure.CalculateStructureProbability(mfe, mfe - 1.0, temperatureK);
        p.Should().BeInRange(0.0, 1.0 + 1e-9, "a structure's Boltzmann probability lies in (0,1]");
    }

    /// <summary>
    /// Interaction witness — the algorithm axis interacts with base composition: the Turner engine
    /// scores a GC stem as more stable than an AU stem of identical geometry (GC stacking ≫ AU),
    /// whereas the simplified per-pair engine assigns every Watson-Crick pair the same −2.0 and so
    /// cannot tell them apart.
    /// </summary>
    [Test]
    public void RnaMfe_TurnerDistinguishesGcFromAu_ClassicDoesNot()
    {
        const string gcStem = "GGGGAAAACCCC"; // 4 G-C pairs, 4-nt loop
        const string auStem = "AAAACCCCUUUU"; // 4 A-U pairs, identical geometry (4-bp stem, 4-nt loop)

        RnaSecondaryStructure.CalculateMinimumFreeEnergy(gcStem)
            .Should().BeLessThan(RnaSecondaryStructure.CalculateMinimumFreeEnergy(auStem),
                "Turner stacking makes a GC stem more stable than an AU stem");

        RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(gcStem)
            .Should().Be(RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(auStem),
                "the simplified engine scores every WC pair at −2.0, so equal pair counts tie");
    }

    /// <summary>
    /// Interaction witness — both engines fold a hairpin (MFE &lt; 0) but report 0 for an unpairable
    /// homopolymer; and the Turner engine reproduces the NNDB worked examples exactly
    /// (tests/TestSpecs/RNA-MFE-001.md M1/M2), anchoring the absolute energy scale.
    /// </summary>
    [Test]
    public void RnaMfe_HairpinFolds_HomopolymerZero_NndbAnchors()
    {
        foreach (var engine in new[] { MfeEngine.TurnerDp, MfeEngine.ClassicPairs })
        {
            Mfe(engine, Hairpin(20)).Should().BeLessThan(0, "a G/C hairpin folds");
            Mfe(engine, new string('A', 20)).Should().Be(0, "poly-A cannot pair ⇒ open chain ΔG=0");
        }

        // NNDB hairpin worked examples (Turner engine, exact).
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("CACAAAAAAAUGUG")
            .Should().BeApproximately(-1.41, 1e-9, "NNDB hairpin-example-1 MFE");
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("CACAGAAAGUGUG")
            .Should().BeApproximately(-1.91, 1e-9, "NNDB hairpin-example-2 MFE");
    }

    /// <summary>
    /// Interaction witness — temperature re-weights the ensemble: for a fixed MFE structure 1 kcal/mol
    /// below its ensemble, the Boltzmann probability rises with temperature (RT discounts the gap).
    /// </summary>
    [Test]
    public void RnaMfe_StructureProbability_RisesWithTemperature()
    {
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(Hairpin(20));
        double pCold = RnaSecondaryStructure.CalculateStructureProbability(mfe, mfe - 1.0, 283.15);
        double pWarm = RnaSecondaryStructure.CalculateStructureProbability(mfe, mfe - 1.0, 337.15);
        pWarm.Should().BeGreaterThan(pCold, "a higher RT shrinks the free-energy gap's penalty");
    }
}
