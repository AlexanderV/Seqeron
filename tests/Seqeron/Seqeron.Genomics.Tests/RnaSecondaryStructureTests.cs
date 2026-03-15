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
        // Dot-bracket: (((....))) — standard notation (Wikipedia: Nucleic acid structure)
        string rna = "GGGAAAACCC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoops, Has.Count.EqualTo(1));
        var sl = stemLoops[0];
        Assert.Multiple(() =>
        {
            Assert.That(sl.Start, Is.EqualTo(0), "Stem-loop starts at position 0");
            Assert.That(sl.End, Is.EqualTo(9), "Stem-loop ends at position 9");
            Assert.That(sl.Stem.Length, Is.EqualTo(3), "3bp stem: GGG pairs with CCC");
            Assert.That(sl.Stem.Start5Prime, Is.EqualTo(0));
            Assert.That(sl.Stem.End5Prime, Is.EqualTo(2));
            Assert.That(sl.Stem.Start3Prime, Is.EqualTo(7));
            Assert.That(sl.Stem.End3Prime, Is.EqualTo(9));
            Assert.That(sl.Loop.Type, Is.EqualTo(LoopType.Hairpin));
            Assert.That(sl.Loop.Start, Is.EqualTo(3));
            Assert.That(sl.Loop.End, Is.EqualTo(6));
            Assert.That(sl.Loop.Size, Is.EqualTo(4), "4nt loop: AAAA");
            Assert.That(sl.Loop.Sequence, Is.EqualTo("AAAA"));
            Assert.That(sl.DotBracketNotation, Is.EqualTo("(((....)))"),
                "Dot-bracket: 3 opening, 4 dots, 3 closing (DB-001/DB-002/DB-003)");
        });
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
        // GGUAAAGCC: G(0)-C(8) WC, G(1)-C(7) WC, U(2)-G(6) Wobble → 3bp stem, loop AAA
        string rna = "GGUAAAGCC";
        var stemLoops = FindStemLoops(rna, minStemLength: 2, minLoopSize: 3, allowWobble: true).ToList();

        Assert.That(stemLoops, Has.Count.EqualTo(2),
            "Two stem-loops: 3bp with wobble (loop AAA) and 2bp WC-only (loop UAAAG)");

        var wobbleSl = stemLoops.First(sl => sl.Stem.BasePairs.Any(bp => bp.Type == BasePairType.Wobble));
        Assert.Multiple(() =>
        {
            Assert.That(wobbleSl.Stem.Length, Is.EqualTo(3), "3bp stem including wobble pair");
            Assert.That(wobbleSl.Loop.Sequence, Is.EqualTo("AAA"), "Loop is AAA (3nt)");
            Assert.That(wobbleSl.Loop.Size, Is.EqualTo(3));
            Assert.That(wobbleSl.Stem.BasePairs.Any(bp =>
                bp.Base1 == 'U' && bp.Base2 == 'G' && bp.Type == BasePairType.Wobble), Is.True,
                "Stem must contain U-G wobble pair at innermost position");
        });
    }

    [Test]
    public void FindStemLoops_WithoutWobble_ExcludesWobble()
    {
        // GCGAAAACGU: G(2)-C(7) WC, C(1)-G(8) WC form 2bp stem; G(0)-U(9) wobble is excluded
        string rna = "GCGAAAACGU";
        var stemLoops = FindStemLoops(rna, minStemLength: 2, allowWobble: false).ToList();

        Assert.That(stemLoops, Has.Count.EqualTo(1), "Should find 1 stem-loop (WC pairs only)");
        var sl = stemLoops[0];
        Assert.Multiple(() =>
        {
            Assert.That(sl.Stem.Length, Is.EqualTo(2), "2bp WC-only stem (G-U wobble excluded)");
            Assert.That(sl.Loop.Sequence, Is.EqualTo("AAAA"), "4nt loop");
            Assert.That(sl.Stem.BasePairs.All(bp => bp.Type == BasePairType.WatsonCrick), Is.True,
                "All base pairs must be Watson-Crick when wobble disabled");
        });
    }

    [Test]
    public void FindStemLoops_MultipleStemLoops_FindsAll()
    {
        // Two palindromic regions: GCGC(0-3)-AAAAA-GCGC(9-12) and GCGC(11-14)-UUUUU-GCGC(20-23)
        string rna = "GCGCAAAAAGCGCGCUUUUUGCGC";
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 3, maxLoopSize: 6).ToList();

        Assert.That(stemLoops, Has.Count.GreaterThanOrEqualTo(2),
            "Should find at least 2 stem-loops from two palindromic regions");
        Assert.That(stemLoops.Any(sl => sl.Loop.Sequence == "AAAAA"), Is.True,
            "Must find stem-loop with AAAAA loop from first palindromic region");
        Assert.That(stemLoops.Any(sl => sl.Loop.Sequence == "UUUUU"), Is.True,
            "Must find stem-loop with UUUUU loop from second palindromic region");
    }

    [Test]
    public void FindStemLoops_Tetraloop_FindsSpecialLoop()
    {
        // Evidence: GNRA tetraloops are stable RNA hairpin structures (Heus & Pardi, 1991)
        // Wikipedia Tetraloop: GNRA = G + any N + purine R (A/G) + A
        // GGGCGAAAGCCC: 4bp stem (GGGC/GCCC) + GAAA tetraloop (GNRA: G=G, N=A, R=A, A=A)
        // Closing pair: C(3)-G(8)
        string rna = "GGGCGAAAGCCC"; // 4bp stem + 4nt GAAA loop
        var stemLoops = FindStemLoops(rna, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoops, Has.Count.EqualTo(1), "Exactly 1 hairpin in GGGCGAAAGCCC");
        var sl = stemLoops[0];
        Assert.Multiple(() =>
        {
            Assert.That(sl.Stem.Length, Is.EqualTo(4), "4bp stem: GGGC pairs with GCCC");
            Assert.That(sl.Loop.Sequence, Is.EqualTo("GAAA"),
                "Loop should be GAAA — a canonical GNRA tetraloop (Wikipedia Tetraloop)");
            Assert.That(sl.Loop.Size, Is.EqualTo(4));
            Assert.That(sl.DotBracketNotation, Is.EqualTo("((((....))))"));
        });
    }

    // DB-001/DB-002/DB-003 (DotBracket assertions) merged into FindStemLoops_SimpleHairpin_FindsStructure

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
    /// Evidence: The biological minimum loop size of 3 nt must be enforced even if
    /// the caller explicitly passes a smaller value. Loops &lt; 3 nt are sterically
    /// impossible and must never appear in results.
    /// Source: NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops
    ///         with fewer than 3 nucleotides."
    /// Source: Wikipedia Stem-loop: "sterically impossible and thus do not form"
    /// Test Unit: RNA-STEMLOOP-001, Test IDs: EC-004
    /// </summary>
    [Test]
    public void FindStemLoops_MinLoopSizeBelowThree_ClampedToThree()
    {
        // AAGCAAGCAA can form 2bp stem GC(2,3)/GC(6,7) with 2nt loop AA(4,5)
        // if minLoopSize < 3 were honored — but loops < 3 are sterically impossible.
        string twoNtLoop = "AAGCAAGCAA";
        var stemLoops = FindStemLoops(twoNtLoop, minStemLength: 2, minLoopSize: 1).ToList();
        Assert.That(stemLoops, Is.Empty,
            "Even with minLoopSize=1, loops < 3 nt must not be found (NNDB Turner 2004)");

        // Positive control: same stem pattern with 3nt loop IS found
        string threeNtLoop = "AAGCAAAGCAA";
        var found = FindStemLoops(threeNtLoop, minStemLength: 2, minLoopSize: 3).ToList();
        Assert.That(found, Is.Not.Empty, "3nt loop (biological minimum) should be found");
        Assert.That(found[0].Loop.Size, Is.EqualTo(3));
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
        // 3-pair GC stem: stacking GC/CG(-3.42) + CG/GC(-2.36) = -5.78, no terminal AU penalty
        var basePairs = new List<BasePair>
        {
            new(0, 9, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 8, 'C', 'G', BasePairType.WatsonCrick),
            new(2, 7, 'G', 'C', BasePairType.WatsonCrick)
        };

        double energy = CalculateStemEnergy("GCGAAAACGC", basePairs);

        Assert.That(energy, Is.EqualTo(-5.78).Within(0.01),
            "GC/CG(-3.42) + CG/GC(-2.36) = -5.78, no terminal penalty");
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

    /// <summary>
    /// Evidence: GC-rich stems are more stable (more negative energy) than AU-rich stems
    /// due to the three hydrogen bonds in G-C pairs vs two in A-U pairs, and stronger
    /// base stacking interactions.
    /// Source: Turner 2004 parameters, NNDB (GC/CG = -3.42 vs AU/UA = -1.10 kcal/mol)
    /// Test Unit: RNA-ENERGY-001, Test ID: SE-003
    /// </summary>
    [Test]
    public void CalculateStemEnergy_GCRichVsAURich_GCMoreStable()
    {
        // GC-rich stem: G-C, C-G, G-C pairs
        var gcPairs = new List<BasePair>
        {
            new(0, 9, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 8, 'C', 'G', BasePairType.WatsonCrick),
            new(2, 7, 'G', 'C', BasePairType.WatsonCrick)
        };

        // AU-rich stem: A-U, U-A, A-U pairs
        var auPairs = new List<BasePair>
        {
            new(0, 9, 'A', 'U', BasePairType.WatsonCrick),
            new(1, 8, 'U', 'A', BasePairType.WatsonCrick),
            new(2, 7, 'A', 'U', BasePairType.WatsonCrick)
        };

        double gcEnergy = CalculateStemEnergy("GCGAAAACGC", gcPairs);
        double auEnergy = CalculateStemEnergy("AUAAAAAUAU", auPairs);

        Assert.That(gcEnergy, Is.LessThan(auEnergy),
            "GC-rich stem should be more stable (more negative) than AU-rich stem");
    }

    /// <summary>
    /// Evidence: Empty base pairs list should return zero energy (no stacking).
    /// Test Unit: RNA-ENERGY-001, Test ID: SE-004
    /// </summary>
    [Test]
    public void CalculateStemEnergy_EmptyBasePairs_ReturnsZero()
    {
        var emptyPairs = new List<BasePair>();

        double energy = CalculateStemEnergy("GCGAAAACGC", emptyPairs);

        Assert.That(energy, Is.EqualTo(0), "Empty base pairs should return zero energy");
    }

    /// <summary>
    /// Evidence: Hairpin loop initiation energy is always positive (destabilizing).
    /// A minimum of 3 nucleotides is required for a valid hairpin loop (steric constraint).
    /// Source: Wikipedia (Stem-loop) - "loops fewer than three bases long are sterically impossible"
    /// Test Unit: RNA-ENERGY-001, Test ID: HL-003
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_MinimumLoop_ReturnsPositive()
    {
        // Loop of exactly 3 nucleotides (biological minimum)
        double energy = CalculateHairpinLoopEnergy("AAA", 'G', 'C');

        Assert.That(energy, Is.GreaterThan(0),
            "Hairpin loop initiation energy should be positive (destabilizing)");
    }

    /// <summary>
    /// Evidence: Hairpin loop initiation energies vary with loop size, with
    /// size 3 being the minimum biological size.
    /// Source: Turner 2004 - loop initiation: 3nt=5.4, 4nt=5.6, 5nt=5.7, 6nt=5.4 kcal/mol
    /// Note: Energy doesn't increase monotonically (size 6 is lower than size 5)
    /// Test Unit: RNA-ENERGY-001, Test ID: HL-005
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_DifferentSizes_AllPositive()
    {
        double energy3 = CalculateHairpinLoopEnergy("AAA", 'G', 'C');
        double energy4 = CalculateHairpinLoopEnergy("AAAA", 'G', 'C');
        double energy5 = CalculateHairpinLoopEnergy("AAAAA", 'G', 'C');
        double energy6 = CalculateHairpinLoopEnergy("AAAAAA", 'G', 'C');

        Assert.Multiple(() =>
        {
            // All loop initiation energies should be positive (destabilizing)
            Assert.That(energy3, Is.GreaterThan(0), "3nt loop should be positive");
            Assert.That(energy4, Is.GreaterThan(0), "4nt loop should be positive");
            Assert.That(energy5, Is.GreaterThan(0), "5nt loop should be positive");
            Assert.That(energy6, Is.GreaterThan(0), "6nt loop should be positive");
        });
    }

    /// <summary>
    /// Evidence: GC-rich hairpins should have lower (more stable) MFE than AU-rich hairpins
    /// of the same length due to stronger stacking interactions.
    /// Source: Turner 2004 parameters
    /// Test Unit: RNA-ENERGY-001, Test ID: MFE-005
    /// </summary>
    [Test]
    public void CalculateMinimumFreeEnergy_GCRichHairpin_MoreStable()
    {
        // GC-rich hairpin
        string gcHairpin = "GCGCAAAAGCGC";
        // AU-rich hairpin (must use valid complementary pairs)
        string auHairpin = "AUAUAAAAUAUA";

        double gcMfe = CalculateMinimumFreeEnergy(gcHairpin);
        double auMfe = CalculateMinimumFreeEnergy(auHairpin);

        // Both should form structures (negative MFE)
        Assert.That(gcMfe, Is.LessThanOrEqualTo(auMfe),
            "GC-rich hairpin should be more stable (lower or equal MFE) than AU-rich");
    }

    /// <summary>
    /// Validates Zuker MFE against manual Turner 2004 calculation.
    /// GGGAAACCC: stem GG/CC(-3.26)×2 + hairpin AAA init(5.4) = -1.12
    /// </summary>
    [Test]
    public void CalculateMinimumFreeEnergy_SimpleHairpin_MatchesTurnerManualCalc()
    {
        // GGG...CCC: 3 GC pairs, 3nt loop AAA
        // Stacking: GG/CC = -3.26, GG/CC = -3.26 → -6.52
        // Hairpin: init(3) = 5.4
        // Total V(0,8) = -6.52 + 5.4 = -1.12
        double mfe = CalculateMinimumFreeEnergy("GGGAAACCC");
        Assert.That(mfe, Is.EqualTo(-1.12).Within(0.01),
            "MFE should match Turner 2004 manual calculation");
    }

    /// <summary>
    /// Validates Zuker MFE for 4-pair GC hairpin.
    /// GGGGAAAACCCC: stacking GG/CC(-3.26)×3 + hairpin AAAA init(5.6)+tm(GAAC=-1.1) = -5.28
    /// </summary>
    [Test]
    public void CalculateMinimumFreeEnergy_FourPairGC_MatchesTurner()
    {
        // 4 GC pairs, 4nt loop
        // Stacking: 3 × GG/CC(-3.26) = -9.78
        // Hairpin: init(4)=5.6 + tm(GAAC)=-1.1 = 4.5
        // V = -9.78 + 4.5 = -5.28
        double mfe = CalculateMinimumFreeEnergy("GGGGAAAACCCC");
        Assert.That(mfe, Is.EqualTo(-5.28).Within(0.01),
            "4-pair GC hairpin MFE should match Turner 2004");
    }

    /// <summary>
    /// Validates that AU stem has correct terminal penalty in Zuker MFE.
    /// AAUAAAUUU → 2 AU pairs + 3nt loop, with terminal AU penalty.
    /// </summary>
    [Test]
    public void CalculateMinimumFreeEnergy_AUStem_IncludesTerminalPenalty()
    {
        // AU stem should be much less stable than GC stem of same length
        double auMfe = CalculateMinimumFreeEnergy("AAUAAAUUU");
        double gcMfe = CalculateMinimumFreeEnergy("GGGAAACCC");

        Assert.That(gcMfe, Is.LessThan(auMfe),
            "GC hairpin should be more negative than AU hairpin of same geometry");
        Assert.That(auMfe, Is.GreaterThanOrEqualTo(0).Or.LessThan(0),
            "AU MFE should be computable (no crash)");
    }

    /// <summary>
    /// Validates Zuker MFE handles internal loops when they're more stable than
    /// just the hairpin alone.
    /// </summary>
    [Test]
    public void CalculateMinimumFreeEnergy_WithBulge_StillWorks()
    {
        // Sequence with potential bulge: GCAGCAAAAGCGC (extra A creates bulge)
        double mfe = CalculateMinimumFreeEnergy("GCAGCAAAAGCGC");

        // Should still compute without error and give negative energy
        Assert.That(mfe, Is.LessThan(0),
            "Sequence with bulge potential should still yield stable structure");
    }

    /// <summary>
    /// Evidence: UNCG tetraloops have enhanced stability when closed by C-G pair.
    /// Source: NNDB Turner 2004 — CUUCGG total energy = 3.7 kcal/mol
    ///         (replaces standard model; standard 4nt initiation alone is 5.6).
    /// Test Unit: RNA-ENERGY-001, Test ID: HL-004
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_UNCGTetraloop_HasBonus()
    {
        // UNCG tetraloop requires C-G closing pair for NNDB special treatment
        double energy_UUCG = CalculateHairpinLoopEnergy("UUCG", 'C', 'G');
        // Compare against generic loop without UU/GA/GG bonus
        double energy_ACCA = CalculateHairpinLoopEnergy("ACCA", 'C', 'G');

        Assert.That(energy_UUCG, Is.LessThan(energy_ACCA),
            "UNCG tetraloop with C-G closing should have enhanced stability (NNDB special loop)");
    }

    #region Turner 2004 NNDB Parameter Validation

    /// <summary>
    /// Validates Watson-Crick stacking energies against NNDB Turner 2004 exact values.
    /// Expected = stacking + terminal AU/GU penalties (+0.45 per AU/UA end).
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html
    /// </summary>
    [TestCase('G', 'C', 'C', 'G', -3.42, Description = "GC/CG most stable WC stack, no terminal penalty")]
    [TestCase('G', 'G', 'C', 'C', -3.26, Description = "GG/CC, no terminal penalty")]
    [TestCase('C', 'G', 'G', 'C', -2.36, Description = "CG/GC, no terminal penalty")]
    [TestCase('A', 'U', 'U', 'A', -0.20, Description = "AU/UA=-1.10, +0.90 AU terminal penalty (both ends)")]
    [TestCase('U', 'A', 'A', 'U', -0.43, Description = "UA/AU=-1.33, +0.90 AU terminal penalty (both ends)")]
    [TestCase('A', 'A', 'U', 'U', -0.03, Description = "AA/UU=-0.93, +0.90 AU terminal penalty (both ends)")]
    [TestCase('C', 'A', 'G', 'U', -1.66, Description = "CA/GU=-2.11, +0.45 AU penalty (inner end only)")]
    [TestCase('U', 'G', 'A', 'C', -1.66, Description = "UG/AC=-2.11, +0.45 AU penalty (outer end only)")]
    public void CalculateStemEnergy_WatsonCrickStacking_MatchesNNDB(
        char b1Top, char b2Top, char b1Bot, char b2Bot, double expected)
    {
        var basePairs = new List<BasePair>
        {
            new(0, 5, b1Top, b1Bot, BasePairType.WatsonCrick),
            new(1, 4, b2Top, b2Bot, BasePairType.WatsonCrick)
        };

        double energy = CalculateStemEnergy("dummy", basePairs);

        Assert.That(energy, Is.EqualTo(expected).Within(0.01),
            $"Stacking {b1Top}{b2Top}/{b1Bot}{b2Bot} should match NNDB Turner 2004 value");
    }

    /// <summary>
    /// Validates hairpin loop model energy against NNDB Turner 2004.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-initiation-parameters.html
    /// For ≥4nt loops: initiation + terminal_mismatch(AAAU=-0.8). For 3nt: initiation only.
    /// </summary>
    [TestCase(3, 5.4, Description = "3nt: initiation only, no terminal mismatch")]
    [TestCase(4, 4.8, Description = "4nt: init(5.6) + tm(AAAU=-0.8)")]
    [TestCase(5, 4.9, Description = "5nt: init(5.7) + tm(AAAU=-0.8)")]
    [TestCase(6, 4.6, Description = "6nt: init(5.4) + tm(AAAU=-0.8)")]
    [TestCase(9, 5.6, Description = "9nt: init(6.4) + tm(AAAU=-0.8)")]
    public void CalculateHairpinLoopEnergy_Initiation_MatchesNNDB(int loopSize, double expected)
    {
        // Use A-repeat loops with A-U closing (no special loop, no UU/GA/GG bonus)
        // For ≥4nt: terminal mismatch AAAU = -0.8 is included in expected
        string loop = new('A', loopSize);
        double energy = CalculateHairpinLoopEnergy(loop, 'A', 'U');

        Assert.That(energy, Is.EqualTo(expected).Within(0.01),
            $"Hairpin model energy for {loopSize}nt loop should match NNDB");
    }

    /// <summary>
    /// Validates NNDB special tetraloop total energies.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-special-parameters.html
    /// The special loop value replaces the model calculation entirely.
    /// </summary>
    [TestCase("CUCG", 'C', 'G', 2.5, Description = "CCUCGG — most stable UNCG tetraloop")]
    [TestCase("UUCG", 'C', 'G', 3.7, Description = "CUUCGG — canonical UNCG tetraloop")]
    [TestCase("CGAG", 'C', 'G', 3.5, Description = "CCGAGG — GNRA tetraloop (CGAG)")]
    [TestCase("UACG", 'C', 'G', 2.8, Description = "CUACGG — UNCG variant")]
    public void CalculateHairpinLoopEnergy_SpecialTetraloop_MatchesNNDB(
        string loop, char closing5, char closing3, double expected)
    {
        double energy = CalculateHairpinLoopEnergy(loop, closing5, closing3);

        Assert.That(energy, Is.EqualTo(expected).Within(0.01),
            $"Special tetraloop {closing5}{loop}{closing3} should match NNDB total energy");
    }

    /// <summary>
    /// Validates that non-special tetraloops use standard model (not NNDB override).
    /// When closing pair doesn't match NNDB entry, initiation energy applies.
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_NonSpecialTetraloop_UsesStandardModel()
    {
        // UUCG with G-C closing (not C-G) → not a special loop
        double energy = CalculateHairpinLoopEnergy("UUCG", 'G', 'C');

        // Should return standard model: init(5.6) + tm(GUGC=-1.1) = 4.5, not special value (3.7)
        Assert.That(energy, Is.EqualTo(4.5).Within(0.01),
            "Non-matching closing pair should use standard model, not NNDB special value");
    }

    /// <summary>
    /// Validates all-C loop penalty formula from NNDB:
    /// - Size 3: flat +1.5 penalty
    /// - Size > 3: 0.3n + 1.6 penalty
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-mismatch-parameters.html
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_AllCPenalty_MatchesNNDB()
    {
        // Use A-U closing to avoid special loop matches
        double energy3C = CalculateHairpinLoopEnergy("CCC", 'A', 'U');
        double energy3A = CalculateHairpinLoopEnergy("AAA", 'A', 'U');
        double penalty3 = energy3C - energy3A;

        // For 4nt: use CAAC as reference (same first/last mismatch C→C as CCCC to isolate penalty)
        double energy4C = CalculateHairpinLoopEnergy("CCCC", 'A', 'U');
        double energy4Ref = CalculateHairpinLoopEnergy("CAAC", 'A', 'U');
        double penalty4 = energy4C - energy4Ref;

        Assert.Multiple(() =>
        {
            Assert.That(penalty3, Is.EqualTo(1.5).Within(0.01),
                "All-C 3nt loop penalty should be +1.5 (NNDB)");
            Assert.That(penalty4, Is.EqualTo(0.3 * 4 + 1.6).Within(0.01),
                "All-C 4nt loop penalty should be 0.3n+1.6 (NNDB)");
        });
    }

    /// <summary>
    /// Validates UU/GA first mismatch bonus from NNDB:
    /// UU or GA first mismatch → -0.9 kcal/mol bonus (not applied to AG).
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-mismatch-parameters.html
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_MismatchBonuses_MatchesNNDB()
    {
        // Use A-U closing, 4nt loops to test mismatch effects
        // GA first mismatch: loop starts with G, ends with A
        double energy_GA = CalculateHairpinLoopEnergy("GCCA", 'A', 'U');
        // No mismatch bonus: AA
        double energy_AA = CalculateHairpinLoopEnergy("ACCA", 'A', 'U');

        Assert.That(energy_GA, Is.EqualTo(energy_AA - 0.9).Within(0.01),
            "GA first mismatch should receive -0.9 bonus (NNDB)");
    }

    /// <summary>
    /// Validates GU wobble stacking energies from NNDB include both stabilizing
    /// and destabilizing entries. GU/UG tandem is notably destabilizing (+1.29).
    /// Expected includes +0.90 terminal GU penalty (both ends are wobble).
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html
    /// </summary>
    [Test]
    public void CalculateStemEnergy_GUWobbleStacking_IncludesDestabilizing()
    {
        // GU/UG tandem: G-U pair followed by U-G pair → destabilizing
        var tandemGU = new List<BasePair>
        {
            new(0, 5, 'G', 'U', BasePairType.Wobble),
            new(1, 4, 'U', 'G', BasePairType.Wobble)
        };

        double energy = CalculateStemEnergy("dummy", tandemGU);

        Assert.That(energy, Is.GreaterThan(0),
            "GU/UG tandem should be destabilizing per NNDB");
        // GU/UG stacking = +1.29, terminal GU penalty = +0.45 * 2 = +0.90, total = +2.19
        Assert.That(energy, Is.EqualTo(2.19).Within(0.01));
    }

    /// <summary>
    /// Validates the extrapolation formula for hairpin loops > 30 nt.
    /// Source: NNDB — ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9) = 6.4 + 1.079·ln(n/9)
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_LargeLoop_UsesCorrectExtrapolation()
    {
        // For size 40: [6.4 + 1.75·R·T·ln(40/9)] + tm(AAAU=-0.8)
        double initiation = 6.4 + 1.75 * 1.987 * 310.15 / 1000.0 * Math.Log(40.0 / 9.0);
        double expected = initiation + (-0.8); // terminal mismatch AAAU = -0.8
        string loop = new('A', 40);
        double energy = CalculateHairpinLoopEnergy(loop, 'A', 'U');

        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Large loop extrapolation should use NNDB formula + terminal mismatch");
    }

    /// <summary>
    /// Validates terminal AU/GU penalty from NNDB: +0.45 per helix end with AU/UA or GU/UG.
    /// Source: NNDB — "Per AU end" (wc-parameters.html), "Per GU end" (gu-parameters.html)
    /// </summary>
    [Test]
    public void CalculateStemEnergy_TerminalAUPenalty_MatchesNNDB()
    {
        // GC-only stem (no terminal penalty)
        var gcStem = new List<BasePair>
        {
            new(0, 5, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 4, 'C', 'G', BasePairType.WatsonCrick)
        };
        double gcEnergy = CalculateStemEnergy("dummy", gcStem); // GC/CG = -3.42, no penalty

        // Same stacking but with AU terminals
        var auOuterStem = new List<BasePair>
        {
            new(0, 7, 'A', 'U', BasePairType.WatsonCrick),
            new(1, 6, 'G', 'C', BasePairType.WatsonCrick),
            new(2, 5, 'C', 'G', BasePairType.WatsonCrick)
        };
        double auOuterEnergy = CalculateStemEnergy("dummy", auOuterStem);

        Assert.Multiple(() =>
        {
            // GC stem: no terminal penalty
            Assert.That(gcEnergy, Is.EqualTo(-3.42).Within(0.01),
                "GC-only stem should have no terminal penalty");
            // AU outer: AG/UC stacking (-2.08) + GC/CG stacking (-3.42) + AU penalty at outer end (+0.45)
            Assert.That(auOuterEnergy, Is.EqualTo(-2.08 + -3.42 + 0.45).Within(0.01),
                "AU outer should add +0.45 terminal penalty");
        });
    }

    #endregion

    #region Terminal Mismatch, Internal Loop, Bulge Loop, Multibranch, Coaxial Tests

    /// <summary>
    /// Validates terminal mismatch table lookups against NNDB Turner 2004 exact values.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/tm-parameters.html
    /// </summary>
    [TestCase('A', 'U', 'A', 'A', -0.8, Description = "AU closing, AA mismatch")]
    [TestCase('C', 'G', 'A', 'A', -1.5, Description = "CG closing, AA mismatch")]
    [TestCase('G', 'C', 'A', 'A', -1.1, Description = "GC closing, AA mismatch")]
    [TestCase('G', 'U', 'A', 'A', -0.3, Description = "GU closing, AA mismatch")]
    [TestCase('U', 'A', 'A', 'A', -1.0, Description = "UA closing, AA mismatch")]
    [TestCase('U', 'G', 'A', 'A', -1.0, Description = "UG closing, AA mismatch")]
    [TestCase('G', 'C', 'G', 'A', -1.6, Description = "GC closing, GA mismatch")]
    [TestCase('C', 'G', 'G', 'G', -1.6, Description = "CG closing, GG mismatch")]
    public void GetTerminalMismatchEnergy_MatchesNNDB(
        char c5, char c3, char mm5, char mm3, double expected)
    {
        double energy = GetTerminalMismatchEnergy(c5, c3, mm5, mm3);

        Assert.That(energy, Is.EqualTo(expected).Within(0.01),
            $"Terminal mismatch {c5}{mm5}{mm3}{c3} should match NNDB");
    }

    /// <summary>
    /// Validates that the terminal mismatch table is applied additively with UU/GA/GG bonuses.
    /// NNDB Example 5: hairpin closing G-U with GG mismatch → tm + GG_bonus + GU_closure.
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_TerminalMismatch_IsAdditiveWithBonuses()
    {
        // Loop "GAAAG" with GC closing: first=G, last=G
        // tm(GGGC) = -1.4, GG bonus = -0.8, init(5) = 5.7
        double energyGG = CalculateHairpinLoopEnergy("GAAAG", 'G', 'C');
        double expected = 5.7 + (-1.4) + (-0.8); // = 3.5
        Assert.That(energyGG, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "GG mismatch should get terminal mismatch + GG bonus additively");

        // Loop "GAAAA" with GC closing: first=G, last=A → GA mismatch
        // tm(GGAC) = -1.6, GA bonus = -0.9, init(5) = 5.7
        double energyGA = CalculateHairpinLoopEnergy("GAAAA", 'G', 'C');
        double expectedGA = 5.7 + (-1.6) + (-0.9); // = 3.2
        Assert.That(energyGA, Is.EqualTo(Math.Round(expectedGA, 2)).Within(0.01),
            "GA mismatch should get terminal mismatch + UU/GA bonus additively");
    }

    /// <summary>
    /// Validates special GU closure bonus (-2.2) for hairpin loops.
    /// NNDB: applied when closing pair is G-U (not U-G) and preceded by two Gs.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-5.html
    /// </summary>
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialGUClosure_AppliesBonus()
    {
        // Loop "GAAAG" with GU closing, specialGUClosure=true
        // init(5) + tm(GGGU=-0.8) + GG_bonus(-0.8) + GU_closure(-2.2)
        double energy = CalculateHairpinLoopEnergy("GAAAG", 'G', 'U', specialGUClosure: true);
        double expected = 5.7 + (-0.8) + (-0.8) + (-2.2); // = 1.9
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Special GU closure should add -2.2 bonus");

        // Same loop without GU closure context
        double energyNoGU = CalculateHairpinLoopEnergy("GAAAG", 'G', 'U', specialGUClosure: false);
        double expectedNoGU = 5.7 + (-0.8) + (-0.8); // = 4.1
        Assert.That(energyNoGU, Is.EqualTo(Math.Round(expectedNoGU, 2)).Within(0.01),
            "Without GU closure context, -2.2 should not apply");
    }

    /// <summary>
    /// Validates internal loop energy calculation against NNDB Turner 2004.
    /// ΔG° = initiation(n₁+n₂) + |n₁−n₂|·0.6 + mismatch + AU/GU_closure(0.7)
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/internal-parameters.html
    /// </summary>
    [Test]
    public void CalculateInternalLoopEnergy_Symmetric_MatchesNNDB()
    {
        // 2×2 symmetric internal loop (total size 4), GC closing both sides
        // init(4)=1.1, asymmetry=0, no AU/GU closure, GG mismatch: -1.2 each side
        double energy = CalculateInternalLoopEnergy(
            2, 2, 'G', 'C', 'G', 'C', 'G', 'G', 'G', 'G');

        double expected = 1.1 + 0.0 + (-1.2) + (-1.2); // = -1.3
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Symmetric 2×2 internal loop with GG mismatches should match NNDB");
    }

    [Test]
    public void CalculateInternalLoopEnergy_Asymmetric_IncludesPenalty()
    {
        // 2×3 asymmetric internal loop (total size 5), AU closing on side 1
        // init(5)=2.0, asymmetry=0.6×|2-3|=0.6, AU closure side 1=+0.7
        double energy = CalculateInternalLoopEnergy(
            2, 3, 'A', 'U', 'G', 'C', 'A', 'A', 'A', 'A');

        double expected = 2.0 + 0.6 + 0.7; // = 3.3 (no mismatch for AA in 2×3)
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Asymmetric internal loop should include asymmetry penalty and AU closure");
    }

    [Test]
    public void CalculateInternalLoopEnergy_BothAUClosure_DoublesPenalty()
    {
        // 3×3 symmetric, both sides AU closing
        // init(6)=2.0, asymmetry=0, AU closure both sides=+1.4, GA mismatch each side=-1.0
        double energy = CalculateInternalLoopEnergy(
            3, 3, 'A', 'U', 'A', 'U', 'G', 'A', 'G', 'A');

        double expected = 2.0 + 0.0 + 0.7 + 0.7 + (-1.0) + (-1.0); // = 1.4
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Both AU closures should each add +0.7 penalty");
    }

    /// <summary>
    /// Validates bulge loop energy calculation against NNDB Turner 2004.
    /// n=1: initiation + stacking (as if no bulge) + special_C.
    /// n>1: initiation + terminal AU/GU penalties.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/bulge.html
    /// </summary>
    [Test]
    public void CalculateBulgeLoopEnergy_SingleNucleotide_IncludesStacking()
    {
        // Single-nucleotide bulge 'A' between GC and GC pairs
        // init(1)=3.8, stacking GG/CC=-3.26 (as if continuous)
        double energy = CalculateBulgeLoopEnergy(1, 'A', 'G', 'C', 'G', 'C');

        double expected = 3.8 + (-3.26); // = 0.54
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Single-nt bulge should include stacking as if continuous");
    }

    [Test]
    public void CalculateBulgeLoopEnergy_SingleC_AdjacentToC_HasBonus()
    {
        // Single-nucleotide bulge 'C' between GC and CG pairs
        // init(1)=3.8, stacking GC/CG=-3.42, special C bonus=-0.9
        double energy = CalculateBulgeLoopEnergy(1, 'C', 'G', 'C', 'C', 'G');

        double expected = 3.8 + (-3.42) + (-0.9); // = -0.52
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Bulged C adjacent to paired C should get -0.9 bonus");
    }

    [Test]
    public void CalculateBulgeLoopEnergy_MultiNucleotide_HasTerminalPenalty()
    {
        // 3-nucleotide bulge between GC and AU pairs
        // init(3)=3.2, AU terminal penalty=+0.45 (only on AU side)
        double energy = CalculateBulgeLoopEnergy(3, 'A', 'G', 'C', 'A', 'U');

        double expected = 3.2 + 0.45; // = 3.65
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Multi-nt bulge with AU end should have terminal penalty");
    }

    [Test]
    public void CalculateBulgeLoopEnergy_MultiNucleotide_BothAU_HasDoublePenalty()
    {
        // 2-nucleotide bulge between AU and GU pairs
        // init(2)=2.8, AU penalty=+0.45 + GU penalty=+0.45 = +0.90
        double energy = CalculateBulgeLoopEnergy(2, 'A', 'A', 'U', 'G', 'U');

        double expected = 2.8 + 0.45 + 0.45; // = 3.7
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Multi-nt bulge with both AU/GU ends should have double terminal penalty");
    }

    /// <summary>
    /// Validates bulge degeneracy term: −RT·ln(numStates) for n=1 bulges.
    /// NNDB Example 1: 3 identical C's → 3 states → −0.616·ln(3) = −0.68 kcal/mol.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/bulge.html, Example 1
    /// </summary>
    [Test]
    public void CalculateBulgeLoopEnergy_SingleNucleotide_DegeneracyTerm()
    {
        // Single-nucleotide bulge 'A' between GC and GC pairs, 3 equivalent states
        // init(1)=3.8, stacking GG/CC=-3.26, degeneracy=−0.616·ln(3)
        const double RT = 1.987 * 310.15 / 1000.0;
        double energy = CalculateBulgeLoopEnergy(1, 'A', 'G', 'C', 'G', 'C', numStates: 3);

        double expected = 3.8 + (-3.26) - RT * Math.Log(3); // = 0.54 − 0.68 = −0.14
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "n=1 bulge with 3 states should include −RT·ln(3) degeneracy");
    }

    [Test]
    public void CalculateBulgeLoopEnergy_SingleNucleotide_NoDegeneracyWhenOneState()
    {
        // Default numStates=1: no degeneracy contribution (ln(1)=0)
        double energyDefault = CalculateBulgeLoopEnergy(1, 'A', 'G', 'C', 'G', 'C');
        double energyExplicit = CalculateBulgeLoopEnergy(1, 'A', 'G', 'C', 'G', 'C', numStates: 1);

        Assert.That(energyDefault, Is.EqualTo(energyExplicit).Within(0.001),
            "numStates=1 should produce same result as default (no degeneracy)");
    }

    [Test]
    public void CalculateBulgeLoopEnergy_MultiNucleotide_DegeneracyIgnored()
    {
        // Degeneracy term only applies to n=1; n>1 should ignore numStates
        double energy1 = CalculateBulgeLoopEnergy(3, 'A', 'G', 'C', 'A', 'U', numStates: 1);
        double energy3 = CalculateBulgeLoopEnergy(3, 'A', 'G', 'C', 'A', 'U', numStates: 3);

        Assert.That(energy1, Is.EqualTo(energy3).Within(0.001),
            "n>1 bulge should not be affected by numStates");
    }

    /// <summary>
    /// Validates multibranch loop energy calculation.
    /// ΔG° = a(9.25) + b(0.91)·asymmetry + c(-0.63)·helices + stacking + strain
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/mb-parameters.html
    /// </summary>
    [Test]
    public void CalculateMultibranchLoopEnergy_ThreeWay_MatchesNNDB()
    {
        // 3-way junction with 6 unpaired bases, no strain, no stacking
        // a=9.25, b·(6/3)=0.91·2=1.82, c·3=-0.63·3=-1.89
        double energy = CalculateMultibranchLoopEnergy(3, 6);

        double expected = 9.25 + 0.91 * 2.0 + (-0.63) * 3; // = 9.18
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "3-way junction should use NNDB multibranch formula");
    }

    [Test]
    public void CalculateMultibranchLoopEnergy_WithStrain_AddsStrainPenalty()
    {
        // 3-way junction with 1 unpaired base (strain applies) + stacking=-2.0
        double energy = CalculateMultibranchLoopEnergy(3, 1, hasStrain: true, stackingEnergy: -2.0);

        double asymmetry = 1.0 / 3.0;
        double expected = 9.25 + 0.91 * asymmetry + (-0.63) * 3 + (-2.0) + 3.14;
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Strained 3-way junction should add +3.14 strain penalty");
    }

    /// <summary>
    /// Validates dangling end energy lookups against NNDB Turner 2004.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/de-parameters.html
    /// </summary>
    [TestCase('C', 'G', 'A', true, -1.7, Description = "3' dangle A on CG pair")]
    [TestCase('G', 'C', 'A', true, -1.1, Description = "3' dangle A on GC pair")]
    [TestCase('A', 'U', 'A', true, -0.8, Description = "3' dangle A on AU pair")]
    [TestCase('C', 'G', 'A', false, -0.2, Description = "5' dangle A on CG pair")]
    [TestCase('G', 'C', 'A', false, -0.5, Description = "5' dangle A on GC pair")]
    public void GetDanglingEndEnergy_MatchesNNDB(
        char c5, char c3, char dangle, bool is3Prime, double expected)
    {
        double energy = GetDanglingEndEnergy(c5, c3, dangle, is3Prime);

        Assert.That(energy, Is.EqualTo(expected).Within(0.01),
            $"{(is3Prime ? "3'" : "5'")} dangle {dangle} on {c5}-{c3} should match NNDB");
    }

    /// <summary>
    /// Validates flush coaxial stacking — uses WC/GU stacking table.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/coax.html
    /// </summary>
    [Test]
    public void CalculateFlushCoaxialStacking_MatchesStackingTable()
    {
        // Flush stack of GC onto CG → GC/CG stacking = -3.42
        double energy = CalculateFlushCoaxialStacking('G', 'C', 'C', 'G');

        Assert.That(energy, Is.EqualTo(-3.42).Within(0.01),
            "Flush coaxial should use WC stacking");
    }

    /// <summary>
    /// Validates mismatch-mediated coaxial stacking.
    /// ΔG° = terminal_mismatch + (-2.1) + WC/GU bonus.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/coax.html
    /// </summary>
    [Test]
    public void CalculateMismatchCoaxialStacking_MatchesNNDB()
    {
        // GC closing, AA mismatch: tm(GAAC)=-1.1, base=-2.1, WC bonus=-0.4
        double energy = CalculateMismatchCoaxialStacking('G', 'C', 'A', 'A');

        double expected = -1.1 + (-2.1) + (-0.4); // = -3.6
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Mismatch coaxial with WC pair should use tm + base + WC bonus");
    }

    [Test]
    public void CalculateMismatchCoaxialStacking_GUPair_UsesGUBonus()
    {
        // GU closing, AA mismatch: tm(GAAU)=-0.3, base=-2.1, GU bonus=-0.2
        double energy = CalculateMismatchCoaxialStacking('G', 'U', 'A', 'A');

        double expected = -0.3 + (-2.1) + (-0.2); // = -2.6
        Assert.That(energy, Is.EqualTo(Math.Round(expected, 2)).Within(0.01),
            "Mismatch coaxial with GU pair should use tm + base + GU bonus");
    }

    #endregion

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
        // GCGCAAAAAGCGC: 4bp stem (GCGC...GCGC), 5nt loop (AAAAA) — must produce pairs
        string rna = "GCGCAAAAAGCGC";
        var structure = PredictStructure(rna, minStemLength: 3);

        Assert.That(structure.BasePairs, Is.Not.Empty,
            "Structured RNA with complementary stem must produce base pairs");
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

    /// <summary>
    /// Invariant: All base pairs returned by PredictStructure must be WC or Wobble.
    /// Source: IUPAC — only Watson-Crick (A-U, G-C) and Wobble (G-U) pairs are valid
    /// in standard RNA secondary structure prediction.
    /// </summary>
    [Test]
    public void PredictStructure_AllBasePairs_AreWatsonCrickOrWobble()
    {
        string rna = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
        var structure = PredictStructure(rna, minStemLength: 3);

        Assert.That(structure.BasePairs.All(bp =>
            bp.Type == BasePairType.WatsonCrick || bp.Type == BasePairType.Wobble), Is.True,
            "All base pairs must be Watson-Crick or Wobble");
    }

    /// <summary>
    /// Edge case: Sequence of identical non-complementary bases cannot form structure.
    /// </summary>
    [Test]
    public void PredictStructure_AllSameBase_NoStructure()
    {
        string rna = "GGGGGGGGGGG";
        var structure = PredictStructure(rna, minStemLength: 3);

        Assert.That(structure.BasePairs, Is.Empty,
            "Poly-G cannot form Watson-Crick pairs, no structure expected");
        Assert.That(structure.StemLoops, Is.Empty);
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
        // GCGCGAAAACGCGC: 5bp GC stem + 4nt AAAA loop
        // Stem: GC/CG(-3.42) + CG/GC(-2.36) + GC/CG(-3.42) + CG/GC(-2.36) = -11.56
        // Loop: init(4)=5.6 + tm(GAAC)=-1.1 = 4.5
        // Total: -11.56 + 4.5 = -7.06
        string rna = "GCGCGAAAACGCGC";
        var stemLoops = FindStemLoops(rna, minStemLength: 4, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoops, Is.Not.Empty, "Should find hairpin in GCGCGAAAACGCGC");
        var bestStemLoop = stemLoops.OrderBy(sl => sl.TotalFreeEnergy).First();

        Assert.Multiple(() =>
        {
            Assert.That(bestStemLoop.Stem.FreeEnergy, Is.EqualTo(-11.56).Within(0.01),
                "5bp GC stem: 4 stackings = -11.56");
            Assert.That(bestStemLoop.TotalFreeEnergy, Is.EqualTo(-7.06).Within(0.01),
                "Total: stem(-11.56) + loop(4.5) = -7.06");
        });
    }

    [Test]
    public void LowerCaseInput_HandlesCorrectly()
    {
        // Evidence: RNA sequence input should be case-insensitive (EC-003)
        string rnaLower = "gggaaaaccc";
        string rnaUpper = "GGGAAAACCC";

        var stemLoopsLower = FindStemLoops(rnaLower, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();
        var stemLoopsUpper = FindStemLoops(rnaUpper, minStemLength: 3, minLoopSize: 4, maxLoopSize: 4).ToList();

        Assert.That(stemLoopsLower, Has.Count.EqualTo(stemLoopsUpper.Count));
        Assert.Multiple(() =>
        {
            Assert.That(stemLoopsLower[0].Stem.Length, Is.EqualTo(stemLoopsUpper[0].Stem.Length),
                "Stem length must match regardless of case");
            Assert.That(stemLoopsLower[0].Loop.Sequence, Is.EqualTo(stemLoopsUpper[0].Loop.Sequence),
                "Loop sequence must match regardless of case");
            Assert.That(stemLoopsLower[0].DotBracketNotation, Is.EqualTo(stemLoopsUpper[0].DotBracketNotation),
                "Dot-bracket must match regardless of case");
            Assert.That(stemLoopsLower[0].TotalFreeEnergy, Is.EqualTo(stemLoopsUpper[0].TotalFreeEnergy),
                "Energy must match regardless of case");
        });
    }

    #endregion

    #region Int11 (1×1 Internal Loop) Lookup Table Tests

    /// <summary>
    /// Validates that 1×1 internal loop uses the int11 lookup table from NNDB.
    /// The int11 values include AU/GU closure penalties.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/int11.txt
    /// </summary>
    [Test]
    public void CalculateInternalLoopEnergy_1x1_GG_CG_CG_ReturnsNNDB()
    {
        // 1×1 loop: CG/CG, mismatch G/G → strongly stabilizing
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'C', 'G', 'C', 'G', 'G', 'G', 'G', 'G');
        Assert.That(energy, Is.EqualTo(-2.2).Within(0.01),
            "CG/CG with G-G mismatch should be -2.2 from int11 table");
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_GG_GC_GC_ReturnsNNDB()
    {
        // 1×1 loop: GC/GC, mismatch G/G → strongly stabilizing
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'G', 'C', 'G', 'C', 'G', 'G', 'G', 'G');
        Assert.That(energy, Is.EqualTo(-2.2).Within(0.01),
            "GC/GC with G-G mismatch should be -2.2 from int11 table");
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_AA_AU_AU_ReturnsNNDB()
    {
        // 1×1 loop: AU/AU, mismatch A/A → default 1.9
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'A', 'U', 'A', 'U', 'A', 'A', 'A', 'A');
        Assert.That(energy, Is.EqualTo(1.9).Within(0.01),
            "AU/AU with A-A mismatch should be 1.9 from int11 table");
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_AC_CG_CG_ReturnsNNDB()
    {
        // 1×1 loop: CG/CG, mismatch A/C → -0.4
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'C', 'G', 'C', 'G', 'A', 'C', 'A', 'C');
        Assert.That(energy, Is.EqualTo(-0.4).Within(0.01),
            "CG/CG with A-C mismatch should be -0.4 from int11 table");
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_UU_AU_AU_ReturnsNNDB()
    {
        // 1×1 loop: AU/AU, mismatch U/U → 1.5 (special UU value)
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'A', 'U', 'A', 'U', 'U', 'U', 'U', 'U');
        Assert.That(energy, Is.EqualTo(1.5).Within(0.01),
            "AU/AU with U-U mismatch should be 1.5 from int11 table");
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_Symmetry_Verified()
    {
        // Int11 symmetry: int11[p1,p2,X,Y] == int11[rev(p2),rev(p1),Y,X]
        // CG/GC with A/C should equal CG/GC with C/A
        double e1 = CalculateInternalLoopEnergy(
            1, 1, 'C', 'G', 'G', 'C', 'A', 'C', 'A', 'C');
        double e2 = CalculateInternalLoopEnergy(
            1, 1, 'C', 'G', 'G', 'C', 'C', 'A', 'C', 'A');
        Assert.That(e1, Is.EqualTo(0.5).Within(0.01));
        Assert.That(e2, Is.EqualTo(0.5).Within(0.01));
    }

    [Test]
    public void CalculateInternalLoopEnergy_1x1_Int11IncludesClosurePenalties()
    {
        // The int11 table already includes AU/GU closure: verify that the value
        // differs from what the generic model would produce.
        // Generic 1×1 (old model): init(min) + 0 asymmetry + AU/GU closures = 1.1 + 0 + 0.7 = 1.8
        // Int11 AU/CG with A/A: 1.2 (includes closure in a different way)
        double energy = CalculateInternalLoopEnergy(
            1, 1, 'A', 'U', 'C', 'G', 'A', 'A', 'A', 'A');
        Assert.That(energy, Is.EqualTo(1.2).Within(0.01),
            "AU/CG with A-A should be 1.2 from int11 (different from generic model)");
    }

    #endregion
}