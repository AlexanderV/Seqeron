namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for CODON-OPT-001: windowed GC constraint of the GC-balancing pass
/// (review 2026-09, B02 F43, audit item A1-7).
///
/// Canonical methods: <c>CodonOptimizer.EvaluateGcContent(...)</c> and the additive
/// <c>gcWindow</c> parameter of <c>CodonOptimizer.OptimizeSequence(...)</c>.
///
/// Reference: DNA Chisel 3.2.16 <c>EnforceGCContent(mini, maxi, window).evaluate(problem)</c>
/// (<c>builtin_specifications/EnforceGCContent.py</c>, <c>biotools/gc_content.py</c>,
/// <c>biotools/indices_operations.group_nearby_segments</c>). Every expected score, breach-window
/// set and breach segment below is the value that library returns for the same arguments, and the
/// post-repair scores were obtained by feeding our optimizer output back into dnachisel's own
/// <c>evaluate</c> (docs/Evidence/CODON-OPT-001-Evidence.md, review 2026-09 F43).
/// </summary>
[TestFixture]
[Category("Codon")]
[Category("CODON-OPT-001")]
public class CodonOptimizer_WindowedGcContent_Tests
{
    // AT-rich head + GC-rich block + AT-rich tail (54 nt), protein MKEKILKGGGAAAMKKF*.
    private const string MixedGene = "ATGAAAGAAAAAATTCTTAAAGGCGGCGGCGCGGCGGCCATGAAAAAATTTTAA";

    // GC-rich (51 nt), protein MAAAAAALLLRRGGAA*.
    private const string GcRichGene = "ATGGCCGCCGCCGCGGCGGCCCTGCTGCTGCGCCGCGGCGGCGCCGCCTGA";

    // AT-rich (54 nt), protein MKFKFKFKFKFKFKFKF* — Lys/Phe only, so no synonymous choice reaches
    // 40 % GC in any 20-nt window (dnachisel's own solver also reports NoSolutionError).
    private const string AtRichGene = "ATGAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTTAA";

    // M-(KE)x7-stop, 48 nt: starts far below 40 % GC and IS repairable (Lys AAA/AAG, Glu GAA/GAG).
    private const string RepairableLowGcGene = "ATGAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAATAA";

    // M-(GS)x7-stop, 48 nt: starts far above 60 % GC and IS repairable (Gly GGC/GGU, Ser AGC/UCU).
    private const string RepairableHighGcGene = "ATGGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCTAA";

    // AT-rich head (36 nt) + balanced middle (105 nt) + AT-rich tail (36 nt) = 177 nt; the two
    // breaching stretches are far enough apart to be reported as two separate segments.
    private const string TwoSegmentGene =
        "AAATTTAAATTTAAATTTAAATTTAAATTTAAATTT" +
        "ATGGCAACTGCAACTGCAACTATGGCAACTGCAACTGCAACTATGGCAACTGCAACTGCAACT" +
        "ATGGCAACTGCAACTGCAACTATGGCAACTGCAACTGCAACT" +
        "AAATTTAAATTTAAATTTAAATTTAAATTTAAATTT";

    #region Equivalence with DNA Chisel EnforceGCContent.evaluate

    [Test]
    [Description("R1: windowed evaluation score equals dnachisel EnforceGCContent(mini, maxi, window).evaluate(problem).score")]
    public void EvaluateGcContent_Windowed_ScoreEqualsDnaChisel()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CodonOptimizer.EvaluateGcContent(MixedGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-4.25).Within(1e-12), "dnachisel: -4.25");
            Assert.That(CodonOptimizer.EvaluateGcContent(GcRichGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-9.45).Within(1e-12), "dnachisel: -9.45");
            Assert.That(CodonOptimizer.EvaluateGcContent(AtRichGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-13.85).Within(1e-12), "dnachisel: -13.85");
            Assert.That(CodonOptimizer.EvaluateGcContent(TwoSegmentGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-20.10).Within(1e-12), "dnachisel: -20.1");
            Assert.That(CodonOptimizer.EvaluateGcContent(MixedGene, 0.30, 0.70, window: 9).Score,
                Is.EqualTo(-7.655555555556).Within(1e-11), "dnachisel: -7.655555555556");
            Assert.That(CodonOptimizer.EvaluateGcContent(RepairableLowGcGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-6.70).Within(1e-12), "dnachisel: -6.7");
            Assert.That(CodonOptimizer.EvaluateGcContent(RepairableLowGcGene, 0.40, 0.60, window: 30).Score,
                Is.EqualTo(-4.366666666667).Within(1e-11), "dnachisel: -4.366667");
            Assert.That(CodonOptimizer.EvaluateGcContent(RepairableHighGcGene, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-6.40).Within(1e-12), "dnachisel: -6.4");
            Assert.That(CodonOptimizer.EvaluateGcContent(RepairableHighGcGene, 0.40, 0.60, window: 30).Score,
                Is.EqualTo(-4.166666666667).Within(1e-11), "dnachisel: -4.166667");
        });
    }

    [Test]
    [Description("R2: the set of breaching window starts equals dnachisel's (breaches > 0).nonzero()")]
    public void EvaluateGcContent_Window20_BreachWindowStartsEqualDnaChisel()
    {
        var mixed = CodonOptimizer.EvaluateGcContent(MixedGene, 0.40, 0.60, window: 20);
        var twoSegment = CodonOptimizer.EvaluateGcContent(TwoSegmentGene, 0.40, 0.60, window: 20);

        Assert.Multiple(() =>
        {
            // dnachisel: [0..7, 13..27, 33, 34] — windows 8..12 and 28..32 pass.
            Assert.That(mixed.BreachWindowStarts, Is.EqualTo(new[]
            {
                0, 1, 2, 3, 4, 5, 6, 7,
                13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                33, 34
            }));
            // dnachisel: [0..31, 125..157]
            Assert.That(twoSegment.BreachWindowStarts, Is.EqualTo(
                Enumerable.Range(0, 32).Concat(Enumerable.Range(125, 33)).ToArray()));
        });
    }

    [Test]
    [Description("R3: breach segments equal dnachisel's SpecEvaluation.locations (group_nearby_segments, max_start_spread = 50)")]
    public void EvaluateGcContent_BreachLocations_EqualDnaChisel()
    {
        var mixed = CodonOptimizer.EvaluateGcContent(MixedGene, 0.40, 0.60, window: 20);
        var twoSegment = CodonOptimizer.EvaluateGcContent(TwoSegmentGene, 0.40, 0.60, window: 20);

        Assert.Multiple(() =>
        {
            // dnachisel message: "Out of bound on segments 0-54"
            Assert.That(mixed.BreachLocations.Select(b => (b.Start, b.End)),
                Is.EqualTo(new[] { (0, 54) }));
            // dnachisel message: "Out of bound on segments 0-51, 125-177"
            Assert.That(twoSegment.BreachLocations.Select(b => (b.Start, b.End)),
                Is.EqualTo(new[] { (0, 51), (125, 177) }));
        });
    }

    [Test]
    [Description("R4: global (window = null) evaluation equals dnachisel EnforceGCContent(0.40, 0.60) — whole sequence as the breach location")]
    public void EvaluateGcContent_Global_EqualsDnaChisel()
    {
        var mixed = CodonOptimizer.EvaluateGcContent(MixedGene, 0.40, 0.60);
        var gcRich = CodonOptimizer.EvaluateGcContent(GcRichGene, 0.40, 0.60);
        var atRich = CodonOptimizer.EvaluateGcContent(AtRichGene, 0.40, 0.60);

        Assert.Multiple(() =>
        {
            // dnachisel: score -0.0, "Passed !"
            Assert.That(mixed.Score, Is.EqualTo(0).Within(1e-12));
            Assert.That(mixed.BreachLocations, Is.Empty);
            // dnachisel: -0.262745098039, segment 0-51 (GC fraction 0.862745… > 0.60)
            Assert.That(gcRich.Score, Is.EqualTo(-0.262745098039).Within(1e-11));
            Assert.That(gcRich.BreachLocations.Select(b => (b.Start, b.End)), Is.EqualTo(new[] { (0, 51) }));
            // dnachisel: -0.381481481481, segment 0-54 (GC fraction 0.018518… < 0.40)
            Assert.That(atRich.Score, Is.EqualTo(-0.381481481481).Within(1e-11));
            Assert.That(atRich.BreachLocations.Select(b => (b.Start, b.End)), Is.EqualTo(new[] { (0, 54) }));
        });
    }

    [Test]
    [Description("R5: a window longer than the sequence yields no windows and score 0 (dnachisel evaluates an empty window array); bounds are inclusive; invalid arguments throw")]
    public void EvaluateGcContent_WindowLongerThanSequence_AndInclusiveBounds()
    {
        Assert.Multiple(() =>
        {
            var tooLong = CodonOptimizer.EvaluateGcContent(AtRichGene, 0.40, 0.60, window: AtRichGene.Length + 1);
            Assert.That(tooLong.Score, Is.EqualTo(0));
            Assert.That(tooLong.BreachWindowStarts, Is.Empty);

            // "GCAT" windows of 4 have GC = 0.50 exactly: inclusive bounds → no breach.
            var exactlyAtBounds = CodonOptimizer.EvaluateGcContent("GCATGCATGCAT", 0.50, 0.50, window: 4);
            Assert.That(exactlyAtBounds.Score, Is.EqualTo(0));
            Assert.That(exactlyAtBounds.BreachLocations, Is.Empty);

            Assert.That(CodonOptimizer.EvaluateGcContent("").Score, Is.EqualTo(0));

            Assert.Throws<ArgumentOutOfRangeException>(
                () => CodonOptimizer.EvaluateGcContent(AtRichGene, 0.40, 0.60, window: 0));
            Assert.Throws<ArgumentNullException>(() => CodonOptimizer.EvaluateGcContent(null!));
        });
    }

    #endregion

    #region OptimizeSequence(gcWindow: …) — windowed repair

    [Test]
    [Description("R6: default (gcWindow = null) GC behaviour is unchanged — global pass only, windows may still breach")]
    public void OptimizeSequence_WithoutGcWindow_KeepsGlobalBehaviour()
    {
        var global = CodonOptimizer.OptimizeSequence(
            RepairableHighGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var explicitNull = CodonOptimizer.OptimizeSequence(
            RepairableHighGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: null);

        Assert.Multiple(() =>
        {
            Assert.That(explicitNull.OptimizedSequence, Is.EqualTo(global.OptimizedSequence));
            Assert.That(global.OptimizedSequence, Is.EqualTo("AUGGGUUCUGGUUCUGGUUCUGGUUCUGGCAGCGGCAGCGGCAGCUAA"));
            // Global GC inside [0.40, 0.60] …
            Assert.That(CodonOptimizer.EvaluateGcContent(global.OptimizedSequence, 0.40, 0.60).Score,
                Is.EqualTo(0).Within(1e-12));
            // … while 20-nt windows still breach (dnachisel on this output: -1.5): exactly the
            // limitation that gcWindow removes.
            Assert.That(CodonOptimizer.EvaluateGcContent(global.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-1.5).Within(1e-12));
        });
    }

    [Test]
    [Description("R7: with gcWindow every sliding window satisfies [gcTargetMin, gcTargetMax] when that is achievable; verified by the dnachisel-equivalent evaluation; protein preserved")]
    public void OptimizeSequence_WithGcWindow_SatisfiesEveryWindow_WhenAchievable()
    {
        foreach (string input in new[] { RepairableLowGcGene, RepairableHighGcGene })
        {
            string protein = CodonOptimizer.OptimizeSequence(
                    input, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization)
                .ProteinSequence;

            foreach (int window in new[] { 20, 30 })
            {
                var result = CodonOptimizer.OptimizeSequence(
                    input, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                    gcWindow: window);

                Assert.Multiple(() =>
                {
                    Assert.That(CodonOptimizer.EvaluateGcContent(result.OptimizedSequence, 0.40, 0.60, window: window).Score,
                        Is.EqualTo(0).Within(1e-12),
                        $"every {window}-nt window must satisfy [0.40, 0.60] for {input}");
                    Assert.That(result.ProteinSequence, Is.EqualTo(protein), "protein preserved");
                    Assert.That(result.OptimizedSequence, Has.Length.EqualTo(result.OriginalSequence.Length));
                });
            }
        }
    }

    [Test]
    [Description("R8: locked deterministic windowed outputs (dnachisel scores them 0 on the requested window)")]
    public void OptimizeSequence_WithGcWindow_IsDeterministic_LockedOutputs()
    {
        var low = CodonOptimizer.OptimizeSequence(
            RepairableLowGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);
        var lowAgain = CodonOptimizer.OptimizeSequence(
            RepairableLowGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);
        var high = CodonOptimizer.OptimizeSequence(
            RepairableHighGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);

        Assert.Multiple(() =>
        {
            Assert.That(lowAgain.OptimizedSequence, Is.EqualTo(low.OptimizedSequence), "deterministic");
            Assert.That(low.OptimizedSequence, Is.EqualTo("AUGAAAGAGAAGGAGAAGGAGAAAGAGAAGGAGAAGGAGAAAGAGUAA"));
            Assert.That(low.ProteinSequence, Is.EqualTo("MKEKEKEKEKEKEKE*"));
            Assert.That(low.ChangedCodons, Is.EqualTo(11));
            Assert.That(high.OptimizedSequence, Is.EqualTo("AUGGGCUCUGGUUCUGGUUCUGGCUCUGGUUCUGGUUCUGGCAGCUAA"));
            Assert.That(high.ProteinSequence, Is.EqualTo("MGSGSGSGSGSGSGS*"));
            Assert.That(high.ChangedCodons, Is.EqualTo(10));
        });
    }

    [Test]
    [Description("R9: where the band is unreachable with synonymous codons the pass reduces the breach as far as possible and keeps the protein (dnachisel's own solver reports NoSolutionError for these genes)")]
    public void OptimizeSequence_WithUnreachableGcWindow_ReducesBreach_AndPreservesProtein()
    {
        var atRich = CodonOptimizer.OptimizeSequence(
            AtRichGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);
        var gcRich = CodonOptimizer.OptimizeSequence(
            GcRichGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);
        var twoSegment = CodonOptimizer.OptimizeSequence(
            TwoSegmentGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);
        var mixed = CodonOptimizer.OptimizeSequence(
            MixedGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 20);

        Assert.Multiple(() =>
        {
            // dnachisel on our outputs: -2.4 / -2.0 / -2.6 / -0.2 (inputs: -13.85 / -9.45 / -20.1 / -4.25).
            Assert.That(CodonOptimizer.EvaluateGcContent(atRich.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-2.4).Within(1e-12));
            Assert.That(atRich.ProteinSequence, Is.EqualTo("MKFKFKFKFKFKFKFKF*"));
            Assert.That(CodonOptimizer.EvaluateGcContent(gcRich.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-2.0).Within(1e-12));
            Assert.That(gcRich.ProteinSequence, Is.EqualTo("MAAAAAALLLRRGGAA*"));
            Assert.That(CodonOptimizer.EvaluateGcContent(twoSegment.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-2.6).Within(1e-12));
            Assert.That(CodonOptimizer.EvaluateGcContent(mixed.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-0.2).Within(1e-12));
            Assert.That(mixed.ProteinSequence, Is.EqualTo("MKEKILKGGGAAAMKKF*"));
            Assert.That(mixed.OptimizedSequence,
                Is.EqualTo("AUGAAGGAGAAGAUCCUGAAAGGUGGUGGUGCAGCAGCGAUGAAGAAGUUCUAA"));
        });
    }

    [Test]
    [Description("R10: a single-codon amino acid run cannot escape an impossible band — no exception, protein and frame intact")]
    public void OptimizeSequence_SingleCodonAminoAcids_UnchangedUnderImpossibleBand()
    {
        const string metGene = "ATGATGATGATGATGATGATGATG";
        var result = CodonOptimizer.OptimizeSequence(
            metGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcTargetMin: 0.0, gcTargetMax: 0.10, rareCodonThreshold: 0.15, gcWindow: 9);

        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("MMMMMMMM"));
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGAUGAUGAUGAUGAUGAUGAUG"));
        });
    }

    [Test]
    [Description("R11: gcWindow must be positive")]
    public void OptimizeSequence_NonPositiveGcWindow_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CodonOptimizer.OptimizeSequence(
            AtRichGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            gcWindow: 0));
    }

    [Test]
    [Description("R12: gcWindow also applies to MinimizeSecondary and is ignored by the strategies without a GC pass")]
    public void OptimizeSequence_GcWindow_AppliesToGcPassStrategiesOnly()
    {
        var maximizeCai = CodonOptimizer.OptimizeSequence(
            RepairableLowGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            gcWindow: 20);
        var maximizeCaiPlain = CodonOptimizer.OptimizeSequence(
            RepairableLowGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var minimizeSecondary = CodonOptimizer.OptimizeSequence(
            RepairableLowGcGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
            gcWindow: 20);

        Assert.Multiple(() =>
        {
            Assert.That(maximizeCai.OptimizedSequence, Is.EqualTo(maximizeCaiPlain.OptimizedSequence),
                "MaximizeCAI has no GC pass, so gcWindow changes nothing");
            Assert.That(minimizeSecondary.ProteinSequence, Is.EqualTo("MKEKEKEKEKEKEKE*"));
        });
    }

    #endregion
}
