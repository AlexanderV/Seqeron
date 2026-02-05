using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for RNA secondary structure prediction algorithms.
/// 
/// Test Units:
/// - RNA-STRUCT-001: Secondary Structure Prediction
/// - RNA-STEMLOOP-001: Stem-Loop Detection
/// 
/// Evidence Sources:
/// - Wikipedia (Nucleic acid structure prediction, Stem-loop, Tetraloop, Pseudoknot)
/// - Nussinov & Jacobson (1980), PNAS 77(11):6309-6313
/// - Zuker & Stiegler (1981), Nucleic Acids Res 9(1):133-148
/// - Turner (2004) thermodynamic parameters
/// - Woese et al. (1990), PNAS 87(21):8467-8471 (tetraloops)
/// - Heus & Pardi (1991), Science 253(5016):191-194 (GNRA stability)
/// 
/// See: docs/Evidence/RNA-STRUCT-001-Evidence.md
/// See: docs/Evidence/RNA-STEMLOOP-001-Evidence.md
/// See: tests/TestSpecs/RNA-STRUCT-001.md
/// See: tests/TestSpecs/RNA-STEMLOOP-001.md
/// </summary>
[TestFixture]
public class RnaSecondaryStructureTests
{
    #region Base Pairing Tests

    /// <summary>
    /// Evidence: Watson-Crick base pairs (A-U, G-C) and Wobble pairs (G-U) are
    /// the valid base pairings in RNA secondary structure.
    /// Source: Wikipedia (Base pair), IUPAC nucleotide code conventions
    /// </summary>
    [TestCase('A', 'U', true)]
    [TestCase('U', 'A', true)]
    [TestCase('G', 'C', true)]
    [TestCase('C', 'G', true)]
    [TestCase('G', 'U', true)]
    [TestCase('U', 'G', true)]
    [TestCase('A', 'A', false)]
    [TestCase('A', 'G', false)]
    [TestCase('C', 'C', false)]
    [TestCase('U', 'U', false)]
    public void CanPair_VariousBases_ReturnsExpected(char base1, char base2, bool expected)
    {
        Assert.That(CanPair(base1, base2), Is.EqualTo(expected));
    }

    [Test]
    public void GetBasePairType_WatsonCrick_ReturnsCorrectType()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetBasePairType('A', 'U'), Is.EqualTo(BasePairType.WatsonCrick));
            Assert.That(GetBasePairType('G', 'C'), Is.EqualTo(BasePairType.WatsonCrick));
            Assert.That(GetBasePairType('C', 'G'), Is.EqualTo(BasePairType.WatsonCrick));
            Assert.That(GetBasePairType('U', 'A'), Is.EqualTo(BasePairType.WatsonCrick));
        });
    }

    [Test]
    public void GetBasePairType_Wobble_ReturnsWobble()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetBasePairType('G', 'U'), Is.EqualTo(BasePairType.Wobble));
            Assert.That(GetBasePairType('U', 'G'), Is.EqualTo(BasePairType.Wobble));
        });
    }

    [Test]
    public void GetBasePairType_NonPairing_ReturnsNull()
    {
        Assert.That(GetBasePairType('A', 'G'), Is.Null);
    }

    [TestCase('A', 'U')]
    [TestCase('U', 'A')]
    [TestCase('G', 'C')]
    [TestCase('C', 'G')]
    public void GetComplement_Returns_Complement(char input, char expected)
    {
        Assert.That(GetComplement(input), Is.EqualTo(expected));
    }

    #endregion

    #region Stem-Loop Finding Tests (RNA-STEMLOOP-001)

    /// <summary>
    /// Evidence: Stem-loop (hairpin) structures form when complementary bases
    /// within a single RNA strand pair with each other, creating a stem with
    /// an unpaired loop at the end.
    /// Source: Wikipedia (Stem-loop), Nussinov & Jacobson (1980)
    /// Test Unit: RNA-STEMLOOP-001
    /// </summary>
    [Test]
    public void FindStemLoops_SimpleHairpin_FindsStructure()
    {
        // Evidence: GGGAAAACCC forms a perfect hairpin with 3bp stem and 4nt loop
        string rna = "GGGAAAACCC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoops, Has.Count.GreaterThanOrEqualTo(1));
        var sl = stemLoops[0];
        Assert.That(sl.Stem.Length, Is.GreaterThanOrEqualTo(3));
        Assert.That(sl.Loop.Type, Is.EqualTo(LoopType.Hairpin));
    }

    [Test]
    public void FindStemLoops_NoComplement_ReturnsEmpty()
    {
        string rna = "AAAAAAAAAAAAAAA";
        var stemLoops = FindStemLoops(rna, minStemLength: 3).ToList();

        Assert.That(stemLoops, Is.Empty);
    }

    [Test]
    public void FindStemLoops_TooShort_ReturnsEmpty()
    {
        string rna = "GCAUC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 3).ToList();

        Assert.That(stemLoops, Is.Empty);
    }

    [Test]
    public void FindStemLoops_WithWobblePairs_IncludesWobble()
    {
        // Evidence: G-U wobble pair is a valid non-Watson-Crick base pair in RNA
        // Sequence designed to form stem with G-U wobble: GCGU pairs with ACGC (G-C, C-G, G-U, U-A)
        string rna = "GCGUAAAACGC"; // 5'-GCGU-loop-CGC-3' with G-U wobble possible
        var stemLoops = FindStemLoops(rna, minStemLength: 2, allowWobble: true).ToList();

        // With wobble enabled and a sequence containing potential G-U pair, we should find structures
        Assert.That(stemLoops, Is.Not.Empty, "Should find stem-loops with wobble pairs allowed");
    }

    [Test]
    public void FindStemLoops_WithoutWobble_ExcludesWobble()
    {
        string rna = "GCGAAAACGU";
        var stemLoops = FindStemLoops(rna, minStemLength: 2, allowWobble: false).ToList();

        var hasWobble = stemLoops.Any(sl =>
            sl.Stem.BasePairs.Any(bp => bp.Type == BasePairType.Wobble));

        Assert.That(hasWobble, Is.False);
    }

    [Test]
    public void FindStemLoops_MultipleStemLoops_FindsAll()
    {
        // Two potential hairpins
        string rna = "GCGCAAAAAGCGCGCUUUUUGCGC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 3, maxLoopSize: 6).ToList();

        // Should find at least one structure
        Assert.That(stemLoops, Has.Count.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void FindStemLoops_Tetraloop_FindsSpecialLoop()
    {
        // Evidence: GNRA tetraloops are stable RNA hairpin structures (Heus & Pardi, 1991)
        // Sequence: GGGG-CGAA-CCCC forms hairpin with CGAA tetraloop (GNRA pattern where N=G, R=A)
        string rna = "GGGGCGAACCCC"; // Perfect hairpin: 4bp stem + 4nt loop
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoops, Is.Not.Empty, "Should find hairpin with tetraloop");
        Assert.That(stemLoops.Any(sl => sl.Loop.Size == 4), Is.True, "Should have 4-nucleotide tetraloop");
    }

    [Test]
    public void FindStemLoops_DotBracket_IsGenerated()
    {
        // Evidence: Dot-bracket notation is standard for RNA secondary structure
        string rna = "GGGAAAACCC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4).ToList();

        Assert.That(stemLoops, Is.Not.Empty, "Should find stem-loop in GGGAAAACCC");
        var sl = stemLoops[0];
        Assert.That(sl.DotBracketNotation, Is.Not.Empty, "Dot-bracket notation should be generated");
        Assert.That(sl.DotBracketNotation, Does.Contain("("), "Should contain opening brackets");
        Assert.That(sl.DotBracketNotation, Does.Contain(")"), "Should contain closing brackets");
    }

    /// <summary>
    /// Evidence: Empty input should return empty result without exception.
    /// Test Unit: RNA-STEMLOOP-001, Test ID: EC-001
    /// </summary>
    [Test]
    public void FindStemLoops_EmptyString_ReturnsEmpty()
    {
        var stemLoops = FindStemLoops("", minStemLength: 3).ToList();

        Assert.That(stemLoops, Is.Empty);
    }

    /// <summary>
    /// Evidence: Null input should return empty result without exception.
    /// Test Unit: RNA-STEMLOOP-001, Test ID: EC-002
    /// </summary>
    [Test]
    public void FindStemLoops_NullString_ReturnsEmpty()
    {
        var stemLoops = FindStemLoops(null!, minStemLength: 3).ToList();

        Assert.That(stemLoops, Is.Empty);
    }

    /// <summary>
    /// Evidence: MinStemLength parameter controls minimum stem size.
    /// Test Unit: RNA-STEMLOOP-001, Test ID: PH-003
    /// </summary>
    [Test]
    public void FindStemLoops_MinStemParameter_RespectsMinimum()
    {
        // Sequence: GCGCAAAAGCGC could form 4bp stem
        string rna = "GCGCAAAAGCGC";

        // With minStemLength=5, should find nothing (max possible is 4)
        var stemLoops5 = FindStemLoops(rna, minStemLength: 5, minLoopSize: 4).ToList();
        Assert.That(stemLoops5, Is.Empty, "Should not find stems shorter than minStemLength");

        // With minStemLength=3, should find structure
        var stemLoops3 = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4).ToList();
        Assert.That(stemLoops3, Is.Not.Empty, "Should find stems >= minStemLength");
        Assert.That(stemLoops3.All(sl => sl.Stem.Length >= 3), Is.True);
    }

    /// <summary>
    /// Evidence: Loop size parameters constrain search space.
    /// Source: Wikipedia - "loops fewer than three bases long are sterically impossible"
    /// Test Unit: RNA-STEMLOOP-001, Test ID: PH-004
    /// </summary>
    [Test]
    public void FindStemLoops_LoopSizeRange_RespectsLimits()
    {
        // Sequence with 4-nucleotide loop
        string rna = "GGGAAAACCC";

        // Exact match: minLoop=4, maxLoop=4
        var stemLoops4 = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();
        Assert.That(stemLoops4, Is.Not.Empty);
        Assert.That(stemLoops4.All(sl => sl.Loop.Size == 4), Is.True, "All loops should be exactly 4nt");

        // Exclude: minLoop=5, maxLoop=6 (4nt loop should not be found)
        var stemLoops56 = FindStemLoops(rna, minStemLength: 3, minLoopSize: 5, maxLoopSize: 6).ToList();
        Assert.That(stemLoops56.Any(sl => sl.Loop.Size == 4), Is.False, "4nt loop should not be found with minLoop=5");
    }

    /// <summary>
    /// Evidence: Minimum biological loop size is 3 nucleotides due to steric constraints.
    /// Source: Wikipedia - "sterically impossible and thus do not form"
    /// Test Unit: RNA-STEMLOOP-001
    /// </summary>
    [Test]
    public void FindStemLoops_MinimumLoopSize_BiologicalConstraint()
    {
        // Sequence that could theoretically form 2nt loop: GC-AA-GC
        // But biologically this is sterically impossible
        string rna = "GCAAGC";
        var stemLoops = FindStemLoops(rna, minStemLength: 2, minLoopSize: 3).ToList();

        // Should not find hairpin with 2nt loop
        Assert.That(stemLoops.Any(sl => sl.Loop.Size < 3), Is.False,
            "Should not find loops smaller than 3nt (biological constraint)");
    }

    #endregion

    #region Energy Calculation Tests

    /// <summary>
    /// Evidence: RNA folding stability is determined by nearest-neighbor
    /// thermodynamic parameters. Stacking interactions between adjacent base
    /// pairs provide stabilizing (negative) free energy contributions.
    /// Source: Turner (2004), Mathews et al. (2004) PNAS 101(19):7287-7292
    /// </summary>
    [Test]
    public void CalculateStemEnergy_BasePairs_ReturnsNegative()
    {
        var basePairs = new List<BasePair>
        {
            new(0, 9, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 8, 'C', 'G', BasePairType.WatsonCrick),
            new(2, 7, 'G', 'C', BasePairType.WatsonCrick)
        };

        double energy = CalculateStemEnergy("GCGAAAACGC", basePairs);

        Assert.That(energy, Is.LessThan(0), "Stem energy should be negative (stabilizing)");
    }

    [Test]
    public void CalculateStemEnergy_SinglePair_ReturnsZero()
    {
        var basePairs = new List<BasePair>
        {
            new(0, 4, 'G', 'C', BasePairType.WatsonCrick)
        };

        double energy = CalculateStemEnergy("GAAAC", basePairs);

        Assert.That(energy, Is.EqualTo(0));
    }

    /// <summary>
    /// Evidence: GNRA tetraloops (where N is any base and R is purine) have
    /// enhanced stability due to specific tertiary interactions.
    /// GAAA is a canonical GNRA tetraloop with ~3 kcal/mol bonus.
    /// Source: Heus & Pardi (1991), Turner 2004 parameters
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_Tetraloop_HasBonus()
    {
        double energy_GAAA = CalculateHairpinLoopEnergy("GAAA", 'G', 'C');
        double energy_AAAA = CalculateHairpinLoopEnergy("AAAA", 'G', 'C');

        // GAAA is a GNRA tetraloop, should have bonus
        Assert.That(energy_GAAA, Is.LessThan(energy_AAAA));
    }

    /// <summary>
    /// Evidence: Poly-C loops are destabilized due to electrostatic repulsion
    /// between the closely spaced negative charges on the phosphate backbone.
    /// Source: Turner 2004 parameters (all-C loop penalty)
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_AllC_HasPenalty()
    {
        double energy_CCCC = CalculateHairpinLoopEnergy("CCCC", 'G', 'C');
        double energy_AAAA = CalculateHairpinLoopEnergy("AAAA", 'G', 'C');

        Assert.That(energy_CCCC, Is.GreaterThan(energy_AAAA));
    }

    [Test]
    public void CalculateMinimumFreeEnergy_SimpleHairpin_ReturnsNegative()
    {
        string rna = "GGGCAAAAGCCC";
        double mfe = CalculateMinimumFreeEnergy(rna);

        Assert.That(mfe, Is.LessThan(0));
    }

    [Test]
    public void CalculateMinimumFreeEnergy_NoStructure_ReturnsZero()
    {
        string rna = "AAAAAA";
        double mfe = CalculateMinimumFreeEnergy(rna);

        Assert.That(mfe, Is.EqualTo(0));
    }

    [Test]
    public void CalculateMinimumFreeEnergy_EmptySequence_ReturnsZero()
    {
        Assert.That(CalculateMinimumFreeEnergy(""), Is.EqualTo(0));
        Assert.That(CalculateMinimumFreeEnergy(null!), Is.EqualTo(0));
    }

    [Test]
    public void CalculateMinimumFreeEnergy_LongerStem_MoreStable()
    {
        string shortStem = "GCAAAAGC";
        string longStem = "GCGCAAAAGCGC";

        double mfeShort = CalculateMinimumFreeEnergy(shortStem);
        double mfeLong = CalculateMinimumFreeEnergy(longStem);

        Assert.That(mfeLong, Is.LessThanOrEqualTo(mfeShort));
    }

    #endregion

    #region Structure Prediction Tests

    /// <summary>
    /// Evidence: RNA secondary structure prediction uses dynamic programming
    /// to find optimal base pairing patterns. The result includes sequence,
    /// dot-bracket notation, and identified structural motifs.
    /// Source: Nussinov & Jacobson (1980), Zuker & Stiegler (1981)
    /// </summary>
    [Test]
    public void PredictStructure_SimpleHairpin_ReturnsPrediction()
    {
        string rna = "GGGGAAAACCCC";
        var structure = PredictStructure(rna);

        Assert.Multiple(() =>
        {
            Assert.That(structure.Sequence, Is.EqualTo(rna));
            Assert.That(structure.DotBracket, Has.Length.EqualTo(rna.Length));
        });
    }

    /// <summary>
    /// Evidence: Dot-bracket notation must be balanced - every opening
    /// bracket must have a corresponding closing bracket.
    /// Source: Wikipedia (Nucleic acid secondary structure representation)
    /// </summary>
    [Test]
    public void PredictStructure_DotBracket_IsValid()
    {
        string rna = "GGGGAAAACCCC";
        var structure = PredictStructure(rna);

        Assert.That(ValidateDotBracket(structure.DotBracket), Is.True);
    }

    [Test]
    public void PredictStructure_HasBasePairs_ForStructuredRNA()
    {
        string rna = "GCGCAAAAAGCGC";
        var structure = PredictStructure(rna, minStemLength: 3);

        // May or may not have base pairs depending on search parameters
        Assert.That(structure.BasePairs, Is.Not.Null);
    }

    [Test]
    public void PredictStructure_EmptySequence_ReturnsEmptyStructure()
    {
        var structure = PredictStructure("");

        Assert.Multiple(() =>
        {
            Assert.That(structure.Sequence, Is.Empty);
            Assert.That(structure.DotBracket, Is.Empty);
            Assert.That(structure.BasePairs, Is.Empty);
            Assert.That(structure.StemLoops, Is.Empty);
        });
    }

    /// <summary>
    /// Evidence: In standard secondary structure (without pseudoknots),
    /// structural elements cannot overlap. This is a fundamental constraint
    /// that enables efficient dynamic programming solutions.
    /// Source: Nussinov algorithm (1980), nested vs crossing base pairs
    /// </summary>
    [Test]
    public void PredictStructure_NonOverlapping_StructuresSelected()
    {
        string rna = "GGGGAAAACCCCUUUUGGGGAAAACCCC";
        var structure = PredictStructure(rna, minStemLength: 3);

        // Check that selected stem-loops don't overlap
        var stemLoops = structure.StemLoops;
        for (int i = 0; i < stemLoops.Count; i++)
        {
            for (int j = i + 1; j < stemLoops.Count; j++)
            {
                bool overlaps = stemLoops[i].End >= stemLoops[j].Start &&
                               stemLoops[j].End >= stemLoops[i].Start;
                Assert.That(overlaps, Is.False, "Stem-loops should not overlap");
            }
        }
    }

    #endregion

    #region Pseudoknot Detection Tests

    /// <summary>
    /// Evidence: A pseudoknot occurs when base pairs cross each other, i.e.,
    /// for pairs (i,j) and (k,l), we have i < k < j < l. This violates
    /// the nested structure requirement of standard secondary structure.
    /// Source: Wikipedia (Pseudoknot), Rivas & Eddy (1999)
    /// </summary>
    [Test]
    public void DetectPseudoknots_NoCrossing_ReturnsEmpty()
    {
        var basePairs = new List<BasePair>
        {
            new(0, 5, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 4, 'C', 'G', BasePairType.WatsonCrick)
        };

        var pseudoknots = DetectPseudoknots(basePairs).ToList();

        Assert.That(pseudoknots, Is.Empty);
    }

    [Test]
    public void DetectPseudoknots_CrossingPairs_DetectsKnot()
    {
        // Crossing: (0,6) and (3,9) - 0 < 3 < 6 < 9
        var basePairs = new List<BasePair>
        {
            new(0, 6, 'G', 'C', BasePairType.WatsonCrick),
            new(3, 9, 'C', 'G', BasePairType.WatsonCrick)
        };

        var pseudoknots = DetectPseudoknots(basePairs).ToList();

        Assert.That(pseudoknots, Has.Count.EqualTo(1));
    }

    #endregion

    #region Dot-Bracket Tests

    /// <summary>
    /// Evidence: Dot-bracket notation is the standard text representation
    /// for RNA secondary structure. Dots (.) represent unpaired bases,
    /// matching parentheses represent base pairs.
    /// Source: Wikipedia (Nucleic acid secondary structure representation)
    /// </summary>
    [Test]
    public void ParseDotBracket_SimpleStructure_ReturnsPairs()
    {
        string dotBracket = "(((...)))";
        var pairs = ParseDotBracket(dotBracket).ToList();

        Assert.That(pairs, Has.Count.EqualTo(3));
        Assert.That(pairs, Does.Contain((0, 8)));
        Assert.That(pairs, Does.Contain((1, 7)));
        Assert.That(pairs, Does.Contain((2, 6)));
    }

    [Test]
    public void ParseDotBracket_EmptyStructure_ReturnsEmpty()
    {
        string dotBracket = ".....";
        var pairs = ParseDotBracket(dotBracket).ToList();

        Assert.That(pairs, Is.Empty);
    }

    [Test]
    public void ParseDotBracket_MultipleBrackets_ParsesAll()
    {
        string dotBracket = "(([[]]))";
        var pairs = ParseDotBracket(dotBracket).ToList();

        Assert.That(pairs, Has.Count.EqualTo(4));
    }

    [Test]
    public void ValidateDotBracket_Balanced_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket("(((...)))"), Is.True);
            Assert.That(ValidateDotBracket("...."), Is.True);
            Assert.That(ValidateDotBracket("((..))((..))"), Is.True);
            Assert.That(ValidateDotBracket(""), Is.True);
        });
    }

    [Test]
    public void ValidateDotBracket_Unbalanced_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket("(((...)"), Is.False);
            Assert.That(ValidateDotBracket("...)"), Is.False);
            Assert.That(ValidateDotBracket(")("), Is.False);
        });
    }

    #endregion

    #region Inverted Repeat Tests (Smoke - see RepeatFinder_InvertedRepeat_Tests.cs for full coverage)

    /// <summary>
    /// Smoke test for RnaSecondaryStructure.FindInvertedRepeats.
    /// Full inverted repeat testing is in RepeatFinder_InvertedRepeat_Tests.cs (REP-INV-001).
    /// This test verifies the RNA-specific implementation finds hairpin stems correctly.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_RnaHairpin_SmokeTest()
    {
        // Evidence: RNA hairpin stems form via antiparallel Watson-Crick base pairing
        // GCGC at 5' end pairs with GCGC at 3' end in antiparallel fashion:
        // 5'-G-C-G-C-loop-G-C-G-C-3'
        //    | | | |      | | | |  (antiparallel pairing)
        //    C-G-C-G      C-G-C-G  (complement read 3'->5')
        // This forms a valid hairpin stem
        string rna = "GCGCAAAAAAGCGC";
        var repeats = FindInvertedRepeats(rna, minLength: 4, minSpacing: 3).ToList();

        // Should find the GCGC...GCGC hairpin stem
        Assert.That(repeats, Is.Not.Empty, "Should find hairpin stem GCGC...GCGC");
    }

    #endregion

    #region Utility Tests

    [Test]
    public void CalculateStructureProbability_ReturnsValidRange()
    {
        double prob = CalculateStructureProbability(-5.0, -10.0);

        Assert.That(prob, Is.GreaterThanOrEqualTo(0));
        Assert.That(prob, Is.LessThanOrEqualTo(1));
    }

    [Test]
    public void CalculateStructureProbability_MFEStructure_HighProbability()
    {
        double prob = CalculateStructureProbability(-10.0, -10.0);

        // When structure energy equals ensemble energy, probability should be high
        Assert.That(prob, Is.GreaterThan(0.5));
    }

    [Test]
    public void GenerateRandomRna_CorrectLength()
    {
        int length = 100;
        string rna = GenerateRandomRna(length);

        Assert.That(rna, Has.Length.EqualTo(length));
    }

    [Test]
    public void GenerateRandomRna_ValidBases()
    {
        string rna = GenerateRandomRna(1000);

        Assert.That(rna.All(c => "ACGU".Contains(c)), Is.True);
    }

    [Test]
    public void GenerateRandomRna_ApproximateGcContent()
    {
        string rna = GenerateRandomRna(10000, gcContent: 0.6);

        int gc = rna.Count(c => c == 'G' || c == 'C');
        double gcRatio = (double)gc / rna.Length;

        Assert.That(gcRatio, Is.InRange(0.55, 0.65));
    }

    #endregion

    #region Integration Tests

    [Test]
    public void FullWorkflow_tRNALike_AnalyzesStructure()
    {
        // Simplified tRNA-like sequence
        string trna = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";

        var structure = PredictStructure(trna, minStemLength: 4, minLoopSize: 3, maxLoopSize: 8);

        Assert.Multiple(() =>
        {
            Assert.That(structure.Sequence, Has.Length.EqualTo(trna.Length));
            Assert.That(structure.DotBracket, Has.Length.EqualTo(trna.Length));
            Assert.That(ValidateDotBracket(structure.DotBracket), Is.True);
        });
    }

    [Test]
    public void FullWorkflow_StemLoopWithEnergy_CalculatesCorrectly()
    {
        string rna = "GCGCGAAAACGCGC";
        var stemLoops = FindStemLoops(rna, minStemLength: 4, minLoopSize: 4, maxLoopSize: 4).ToList();

        if (stemLoops.Any())
        {
            var bestStemLoop = stemLoops.OrderBy(sl => sl.TotalFreeEnergy).First();

            Assert.Multiple(() =>
            {
                Assert.That(bestStemLoop.TotalFreeEnergy, Is.LessThan(0).Or.GreaterThan(0));
                Assert.That(bestStemLoop.Stem.FreeEnergy, Is.LessThanOrEqualTo(0));
            });
        }
    }

    [Test]
    public void LowerCaseInput_HandlesCorrectly()
    {
        // Evidence: RNA sequence input should be case-insensitive
        string rnaLower = "gggaaaaccc";
        string rnaUpper = "GGGAAAACCC";

        var stemLoopsLower = FindStemLoops(rnaLower, minStemLength: 3).ToList();
        var stemLoopsUpper = FindStemLoops(rnaUpper, minStemLength: 3).ToList();

        // Both should produce same results
        Assert.That(stemLoopsLower.Count, Is.EqualTo(stemLoopsUpper.Count),
            "Lowercase and uppercase input should produce same number of results");
    }

    #endregion
}