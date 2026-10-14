using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// CODON-OPT-001 — thermodynamic (Turner 2004) window scoring for
/// <see cref="CodonOptimizer.ReduceSecondaryStructure"/> and the MinimizeSecondary strategy
/// (review 2026-09, B02 F41 / F42; audit items A1-6, A1-8).
///
/// Sources
/// - Energy model: NNDB Turner 2004 nearest-neighbour parameters, evaluated by the canonical
///   folding engine <see cref="RnaSecondaryStructure.CalculateMinimumFreeEnergy"/> (Zuker-style
///   O(n³) DP) — the same model ViennaRNA RNAfold uses by default.
/// - Reference implementation for the cross-check: ViennaRNA 2.7.2 (Python RNA module,
///   <c>RNA.fold_compound(window, md).mfe()</c> with <c>md.dangles = 2</c>, the RNAfold default).
/// - DNA Chisel has NO free-energy objective (its <c>AvoidHairpins</c> is the IDT
///   reverse-complement stem rule, verified on dnachisel 3.2.16: no module references ViennaRNA,
///   RNA.fold or a free energy), so the energy model is taken from ViennaRNA / NNDB, not from
///   DNA Chisel.
/// </summary>
[TestFixture]
public class CodonOptimizer_SecondaryStructureScoring_Tests
{
    // Window width used by the structure pass: windowSize 40 → 40/3 + 1 = 14 codons = 42 nt …
    // in whole codons the pass concatenates 14 codons from the codon list, i.e. 42 nt; the
    // locked windows below are the 45-nt (15-codon) spans used for the independent
    // before/after profile, which is scored the same way (Turner 2004 MFE of the window).
    private const int ProfileWindow = 45;

    // Test gene A — a designed CDS with two long complementary blocks (15×Gly GGC, 5×Lys AAG,
    // 15×Ala GCC): the GGC block pairs with the GCC block, so the gene is strongly structured.
    private const string GeneA =
        "AUGGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCAAGAAGAAGAAGAAG" +
        "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCUAA";

    // Test gene B — a 33-residue CDS of mixed composition (MKLADLAQLAFREDLRQALDRFMAQLDKALAE*).
    private const string GeneB =
        "AUGAAACUGGCAGAUCUGGCGCAGCUGGCAUUUCGUGAAGAUCUGCGUCAGGCGCUGGAUCGU" +
        "UUCAUGGCGCAGCUGGAUAAAGCGCUGGCGGAAUAA";

    private const string GeneAProtein = "MGGGGGGGGGGGGGGGKKKKKAAAAAAAAAAAAAAA*";
    private const string GeneBProtein = "MKLADLAQLAFREDLRQALDRFMAQLDKALAE*";

    private static string Translate(string rna) =>
        string.Concat(Enumerable.Range(0, rna.Length / 3)
            .Select(i => GeneticCode.Standard.Translate(rna.Substring(i * 3, 3))));

    private static double SumWindowMfe(string rna)
    {
        double total = 0;
        for (int i = 0; i + ProfileWindow <= rna.Length; i += 3)
            total += RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna.Substring(i, ProfileWindow));
        return total;
    }

    private static double MinWindowMfe(string rna)
    {
        double min = 0;
        for (int i = 0; i + ProfileWindow <= rna.Length; i += 3)
            min = Math.Min(min, RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna.Substring(i, ProfileWindow)));
        return min;
    }

    // ── Lock 1: the window energies the Turner mode scores with, vs ViennaRNA RNAfold ──

    [Test]
    [Description("A1-6 lock 1: the Turner-2004 window energies used for scoring are the canonical "
        + "engine's values and agree with ViennaRNA 2.7.2 RNAfold (dangles=2) on six 45-nt windows "
        + "of the two test genes, within the deviation documented for the engine (missing NNDB "
        + "int21/int22 lookup tables, RNA-STRUCT-001 Evidence D3/D4).")]
    public void TurnerWindowEnergies_MatchViennaRnaRnafold()
    {
        // ours (CalculateMinimumFreeEnergy) vs ViennaRNA 2.7.2 md.dangles = 2, kcal/mol:
        (string Window, double Ours, double ViennaRna)[] windows =
        {
            (GeneA.Substring(0, 45), -25.82, -25.80),
            (GeneA.Substring(3, 45), -29.34, -29.00),
            (GeneA.Substring(6, 45), -25.82, -26.50),
            (GeneB.Substring(0, 45),  -9.22,  -9.20),
            (GeneB.Substring(3, 45), -13.44, -14.00),
            (GeneB.Substring(6, 45), -13.82, -14.00),
        };

        Assert.Multiple(() =>
        {
            foreach (var (window, ours, viennaRna) in windows)
            {
                double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(window);
                Assert.That(mfe, Is.EqualTo(ours).Within(1e-9),
                    $"locked Turner 2004 MFE of window {window}");
                Assert.That(mfe, Is.EqualTo(viennaRna).Within(1.0),
                    $"ViennaRNA 2.7.2 RNAfold (dangles=2) gives {viennaRna:F2} kcal/mol for {window}; "
                    + "max observed deviation over these windows is 0.68 kcal/mol");
            }
        });
    }

    // ── Lock 2/3: Turner mode improves the window-MFE profile and preserves the protein ──

    [Test]
    [Description("A1-6 locks 2+3 (gene A): Turner2004Mfe scoring raises (makes less negative) both "
        + "the summed and the minimum 45-nt window MFE, preserves the protein, and the output is "
        + "locked exactly for determinism.")]
    public void ReduceSecondaryStructure_TurnerMode_GeneA_RaisesWindowMfe_AndPreservesProtein()
    {
        string optimized = CodonOptimizer.ReduceSecondaryStructure(
            GeneA, CodonOptimizer.EColiK12,
            scoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe);

        Assert.Multiple(() =>
        {
            Assert.That(optimized, Is.EqualTo(
                "AUGGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGGGGAGGAGGAGGAAAGAAGAAGAAGAAA" +
                "GCCGCCGCAGCCGCAGCAGCAGCCGCAGCAGCAGCAGCCGCAGCCUAA"),
                "deterministic Turner-mode output");
            Assert.That(Translate(optimized), Is.EqualTo(GeneAProtein), "protein preserved");
            // Summed window MFE −477.92 → −144.26 kcal/mol; worst window −37.34 → −16.04.
            Assert.That(SumWindowMfe(optimized), Is.EqualTo(-144.26).Within(0.01));
            Assert.That(SumWindowMfe(GeneA), Is.EqualTo(-477.92).Within(0.01));
            Assert.That(SumWindowMfe(optimized), Is.GreaterThanOrEqualTo(SumWindowMfe(GeneA)),
                "the summed window MFE must not get worse");
            Assert.That(MinWindowMfe(optimized), Is.EqualTo(-16.04).Within(0.01));
            Assert.That(MinWindowMfe(GeneA), Is.EqualTo(-37.34).Within(0.01));
            Assert.That(MinWindowMfe(optimized), Is.GreaterThanOrEqualTo(MinWindowMfe(GeneA)),
                "the most stable window must not get worse");
        });
    }

    [Test]
    [Description("A1-6 locks 2+3 (gene B): same contract on a mixed-composition CDS.")]
    public void ReduceSecondaryStructure_TurnerMode_GeneB_RaisesWindowMfe_AndPreservesProtein()
    {
        string optimized = CodonOptimizer.ReduceSecondaryStructure(
            GeneB, CodonOptimizer.EColiK12,
            scoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe);

        Assert.Multiple(() =>
        {
            Assert.That(optimized, Is.EqualTo(
                "AUGAAGCUGGCAGAUCUUGCUCAACUUGCAUUUCGUGAAGACCUUCGGCAAGCUUUAGAUCGU" +
                "UUCAUGGCACAAUUAGAUAAAGCAUUAGCAGAAUAA"),
                "deterministic Turner-mode output");
            Assert.That(Translate(optimized), Is.EqualTo(GeneBProtein), "protein preserved");
            // Summed window MFE −252.88 → −96.45 kcal/mol; worst window −17.00 → −7.18.
            Assert.That(SumWindowMfe(optimized), Is.EqualTo(-96.45).Within(0.01));
            Assert.That(SumWindowMfe(GeneB), Is.EqualTo(-252.88).Within(0.01));
            Assert.That(SumWindowMfe(optimized), Is.GreaterThanOrEqualTo(SumWindowMfe(GeneB)));
            Assert.That(MinWindowMfe(optimized), Is.EqualTo(-7.18).Within(0.01));
            Assert.That(MinWindowMfe(GeneB), Is.EqualTo(-17.00).Within(0.01));
            Assert.That(MinWindowMfe(optimized), Is.GreaterThanOrEqualTo(MinWindowMfe(GeneB)));
        });
    }

    [Test]
    [Description("A1-6: the default scoring is unchanged (PairCount) — the opt-in Turner mode does "
        + "not alter the historical default behaviour.")]
    public void ReduceSecondaryStructure_DefaultScoring_IsPairCount()
    {
        string byDefault = CodonOptimizer.ReduceSecondaryStructure(GeneB, CodonOptimizer.EColiK12);
        string explicitPairCount = CodonOptimizer.ReduceSecondaryStructure(
            GeneB, CodonOptimizer.EColiK12,
            scoring: CodonOptimizer.SecondaryStructureScoring.PairCount);

        Assert.Multiple(() =>
        {
            Assert.That(byDefault, Is.EqualTo(explicitPairCount));
            // No 42-nt window of gene B exceeds the 0.5 canonical-pair fraction, so the
            // pair-count pass leaves the gene as it is — while the Turner pass rewrites 19 codons.
            Assert.That(byDefault, Is.EqualTo(GeneB));
        });
    }

    [Test]
    [Description("A1-6: the Turner-mode fold budget bounds the runtime and stops the pass "
        + "deterministically; with a budget of one fold only the first window is rewritten.")]
    public void ReduceSecondaryStructure_TurnerMode_FoldBudget_StopsDeterministically()
    {
        string budgeted = CodonOptimizer.ReduceSecondaryStructure(
            GeneA, CodonOptimizer.EColiK12,
            scoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe,
            maxFoldEvaluations: 1);
        string again = CodonOptimizer.ReduceSecondaryStructure(
            GeneA, CodonOptimizer.EColiK12,
            scoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe,
            maxFoldEvaluations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(budgeted, Is.EqualTo(again), "deterministic under a fold budget");
            Assert.That(Translate(budgeted), Is.EqualTo(GeneAProtein), "protein preserved");
            Assert.That(budgeted, Is.Not.EqualTo(GeneA), "the first window is still rewritten");
            Assert.That(budgeted.Length, Is.EqualTo(GeneA.Length));
        });
    }

    // ── A1-8: MinimizeSecondary actually runs the structure pass ──

    [Test]
    [Description("A1-8: OptimizationStrategy.MinimizeSecondary = balanced codon selection + GC "
        + "pass + the structure pass. On a CDS of two complementary blocks (8×Glu, 8×Phe) the "
        + "structure pass changes the result, which BalancedOptimization alone does not.")]
    public void OptimizeSequence_MinimizeSecondary_RunsStructurePass()
    {
        // 8 × GAA (Glu) followed by 8 × UUU (Phe): GAG-repeat and UUC-repeat blocks are
        // complementary, so the window pair-count score exceeds 0.5 and the pass engages.
        const string gene = "GAAGAAGAAGAAGAAGAAGAAGAAUUUUUUUUUUUUUUUUUUUUUUUU";

        var balanced = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var minimizeSecondary = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MinimizeSecondary);

        Assert.Multiple(() =>
        {
            Assert.That(balanced.OptimizedSequence,
                Is.EqualTo("GAGGAGGAGGAGGAGGAGGAGGAGUUCUUCUUCUUCUUUUUUUUUUUU"),
                "BalancedOptimization: codon selection + GC pass only");
            Assert.That(minimizeSecondary.OptimizedSequence,
                Is.EqualTo("GAGGAAGAAGAAGAAGAAGAAGAAUUCUUCUUCUUCUUCUUCUUCUUU"),
                "MinimizeSecondary additionally runs ReduceSecondaryStructure");
            Assert.That(minimizeSecondary.OptimizedSequence, Is.Not.EqualTo(balanced.OptimizedSequence),
                "review 2026-09 F42: the strategy used to be an alias of BalancedOptimization");
            Assert.That(minimizeSecondary.ProteinSequence, Is.EqualTo("EEEEEEEEFFFFFFFF"),
                "protein preserved");
            Assert.That(Translate(minimizeSecondary.OptimizedSequence), Is.EqualTo("EEEEEEEEFFFFFFFF"));
            Assert.That(minimizeSecondary.OptimizedSequence.Length,
                Is.EqualTo(minimizeSecondary.OriginalSequence.Length));
            // The reported CAI / GC describe the sequence that is actually returned.
            Assert.That(minimizeSecondary.OptimizedCAI, Is.EqualTo(
                CodonOptimizer.CalculateCAI(minimizeSecondary.OptimizedSequence, CodonOptimizer.EColiK12))
                .Within(1e-12));
            Assert.That(minimizeSecondary.ChangedCodons, Is.EqualTo(
                Enumerable.Range(0, 16).Count(i =>
                    minimizeSecondary.OriginalSequence.Substring(i * 3, 3)
                        != minimizeSecondary.OptimizedSequence.Substring(i * 3, 3))));
        });
    }

    [Test]
    [Description("A1-8: MinimizeSecondary can be scored with the Turner 2004 energy model through "
        + "the optional secondaryScoring parameter; on gene B that lowers the folding propensity "
        + "of the balanced sequence substantially (summed window MFE −329.72 → −99.86 kcal/mol).")]
    public void OptimizeSequence_MinimizeSecondary_Turner2004Mfe_OptIn()
    {
        var balanced = CodonOptimizer.OptimizeSequence(
            GeneB, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var turner = CodonOptimizer.OptimizeSequence(
            GeneB, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
            secondaryScoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe);

        Assert.Multiple(() =>
        {
            Assert.That(turner.OptimizedSequence, Is.EqualTo(
                "AUGAAACUUGCCGACUUAGCUCAAUUAGCAUUUCGUGAAGAUUUGCGUCAAGCACUCGAUCGU" +
                "UUUAUGGCACAAUUGGACAAGGCUCUGGCAGAAUAA"),
                "deterministic Turner-scored MinimizeSecondary output");
            Assert.That(turner.ProteinSequence, Is.EqualTo(GeneBProtein));
            Assert.That(Translate(turner.OptimizedSequence), Is.EqualTo(GeneBProtein));
            Assert.That(SumWindowMfe(balanced.OptimizedSequence), Is.EqualTo(-329.72).Within(0.01));
            Assert.That(SumWindowMfe(turner.OptimizedSequence), Is.EqualTo(-99.86).Within(0.01));
            Assert.That(SumWindowMfe(turner.OptimizedSequence),
                Is.GreaterThan(SumWindowMfe(balanced.OptimizedSequence)),
                "the Turner-scored pass must leave the sequence less structured than the balanced one");
            Assert.That(MinWindowMfe(turner.OptimizedSequence),
                Is.GreaterThan(MinWindowMfe(balanced.OptimizedSequence)));
        });
    }
}
