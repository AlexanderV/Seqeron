using NUnit.Framework;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ASSEMBLY-COMPARE-001 (B10 sweep U4): <c>CompareAssemblies</c> = Mash/sourmash exact canonical k-mer
/// comparison; <c>FindSyntenicBlocks</c> = nucmer front end (MUMs + MUMmer 4 mgaps clustering).
/// Fixtures are produced by <c>tools/review/oracles/B10/compare_oracle.py</c> (independent brute-force MUM
/// definition + literal mgaps port + Mash formulas); the generator agreed with the C# on 600/600 random cases
/// (1037 blocks, 222 inverted).
/// </summary>
[TestFixture]
[Category("Unit")]
public class GenomeAssemblyAnalyzer_AssemblyCompare_Tests
{
    private const double Tol = 1e-12;

    private const string A = "CTGTTTCTGTCCATTAGGTGCAGGTATAGGACACCTCACTTCTCCATCGTTTTGATTGCGTGAATCTCTGGGGAAGTTTTCTGGTTACGTGCTGGCGGATACACCCGCCCGGGCCTGAAAACGTTGTGACATTCCATTATCGGACTATCT";
    private const string B = "TAACCTGAACGATGAGCACTAGGTATCGTGCTTCGATATTGCCGTATGGAGGGGGTGGGGCTCGGCATAGATACAAAGAACCTGAGCCGATAACCTATTCTTAGCCTCCACGGTAGTCCTGTCCCTTTGGGCGAACTACTGGGCGCTAAT";
    private const string C = "TACACTTCCGATGAACATAAGCAAAGGCCAGTGAGAGCCCCGACCCTGGGATATCGCCTAAAGATGGAGTAGGGGTAAACTTGGAAGGACGTCGGTCCTCGGACTCTGTGTGGTGATAAAGCGGGGTGCGTGCTCCGTGTCACGTTGCCG";
    private const string D = "CCTAGCTCGACAGGGAATCGTCGTAGCTGATAAGGGGGACTGTCTTTCTGAACTCTTTGGTGTAATAAAGCGACCGACGGTAGGCAGCGACGCAATTCGCTAGTAACTTATAACGCTAGGGCTCTTTCGCTGAGCGGCATCTGATGAAAC";

    private static string Rc(string s) => Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(s);

    public static readonly object[][] OracleCases =
    {
        new object[] { new[] { ("r0", "TTTGCCGCCTGACAAGTCAATGCGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATTAACTG"), ("r1", "TAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAACGGCGA") }, new[] { ("q0", "CAGTTAATCGTTTGTGATGCGACAGTGCCTATAGTCTTGGCATACTGCGCGCTTTACATTTGCTGTGAGAGGTACTGCCCCTACGGATCGCATTGACTTGTCAGGCGGCAAA"), ("q1", "TAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAACGGCGA") }, 8, 12, 30, 133, 158, 127, 0.774390243902439, 0.01699862501911949, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 37, "q0", 74, 111, true), ("r0", 37, 86, "q0", 0, 49, true), ("r1", 0, 59, "q1", 0, 59, false) } },
        new object[] { new[] { ("r0", "ACACGCTCCCGGCCACTGGTGGAGAGAGCCCCTACGAGTGAAATTTAGCTGTTGTGAATAGCACATAGAGTACT"), ("r1", "AAGCAAGCTCCCTTGGACTAAGTTCCGTTCCCTAGCAGTCGGCGCTAACGAGAAGCGGGGGGTTGACATCACCGGGTTGCCGAG") }, new[] { ("q0", "AGTACTCTATGTGCTATTCACAACAGCTAAATTTCACTCGTAGGGGCTCTCGCATGTTCGGCAAAGAACGAATACCTCCACCAGTGGCCGGGAGCGTGT"), ("q1", "AAGCAAGCTCCCTTGGACTAAGTTCCGTTCCCTAGCAGTCGGCGCTAACGAGAAGCGGGGGGTTGACATCACCGGGTTGCCGAG") }, 11, 12, 30, 138, 163, 129, 0.75, 0.014013698166114396, new (string, int, int, string, int, int, bool)[] { ("r0", 23, 73, "q0", 0, 50, true), ("r1", 0, 83, "q1", 0, 83, false) } },
        new object[] { new[] { ("r0", "TCGTATTTGCTTCACTCAAGTAAGTCCTCGTCCTAGATTGCGACAAGAGGCAAAGAGCTTAATGTTTATCTCGTTTGAATGCCTTGGCCTCGCAATAATGTAAATGATGCTAAACCAACACGTTGCGAATGAAATACGTGCTAGTGGGAATGCGAGGGGCTGCTT"), ("r1", "CCCAAGCGGCTTCAGACTTACTTTCGGTTTCTCGTAACACGGTTGGGCCCACCTGACCCGGGAGCTATCTTATTAACTGCAATTACTGCAGAAATCTCTGGTCCAGTCGGAGAAGGGGTTTTTGA") }, new[] { ("q0", "AAGCAGCCCCTCGCATTCCCACTAGCACGTATTTCATTCGCAACGTGTTGCACCCCCTGCGTTACACTAATAATTGTTTAGCATCATTTACATTATTGCGAGGCCAAGGCATTCAAACGAGATAAACATTAAGCTCTTTGCCTCTTGTCGCAATCTAGGACGAGGACTTACTTGAGTGAAGCAAATACGA"), ("q1", "CCCAAGCGGCTTCAGACTTACTTTCGGTTTCTCGTAACACGGTTGGGCCCACCTGACCCGGGAGCTATCTTATTAACTGCAATTACTGCAGAAATCTCTGGTCCAGTCGGAGAAGGGGTTTTTGA") }, 4, 15, 15, 117, 118, 117, 0.9915254237288136, 0.001066099696614385, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 114, "q0", 75, 189, true), ("r0", 115, 164, "q0", 0, 49, true), ("r1", 0, 124, "q1", 0, 124, false) } },
        new object[] { new[] { ("r0", "TGGTCCTGCGGAGATCCGTCnCGATCTCGAACACGCGACTTATGTGACCAACCTAAAGAAATCTACCCAGTAGCCAGCAGGAACATGGAGATGGTGTTGTTCTTTCACGTCCAAAATGTGTATTGTCTGATGGACGGTGTCCA") }, new[] { ("q0", "TGGACACCGTCCATCAGACAATACACATTTTGGACGTGAAAGAACAACACGCCGCCCTCAGTGTATCGTAGGGTACATCTCCATGTTCCTGCTGGCTACTGGGTAGATTTCTTTAGGTTGGTCACATAAGTCGCGTGTTCGAGATCGTGACGGATCTCCGCAGGACCA") }, 8, 15, 65, 128, 161, 121, 0.7202380952380952, 0.022186120244468227, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 92, "q0", 75, 167, true) } },
        new object[] { new[] { ("r0", "TTCAGAAGCAACCGGACTCGnCCTCTTCCGTCGCTGAGTAAGACGGTAAACTGGACGAGGGCTTAGGGAGAGTGGTGCAGACTAAGCTACCACT"), ("r1", "ACACACCTCCTTGACGGTAGTCTCGATCAGTTGATAATAATGCGTATTGGTCTATAGCTCCCCCGATGGAATGTGCTTTGTAATGCATCCGGAGAGGTAGGGGCCAATGCAAGCTGGGAAGGATGAGTAGGAGAACTAGAGGACATTCCGGTGTCAAACTGCTTGTCAACCGTCAAGGAATGCCATCACACCATAGTGTCTTCGTTCAATTAACGCATT") }, new[] { ("q0", "AGTGGTAGCTTAGTCTGCACCACTCTCCCTAAGCCCTCGTCCAGTTTACCTTCTTCTGACGGCCCTTTTCCCGGAGTCTTACTCAGCGACGGAAGAGGCCGAGTCCGGTTGCTTCTGAA"), ("q1", "ACACACCTCCTTGACGGTAGTCTCGATCAGTTGATAATAATGCGTATTGGTCTATAGCTCCCCCGATGGAATGTGCTTTGTAATGCATCCGGAGAGGTAGGGGCCAATGCAAGCTGGGAAGGATGAGTAGGAGAACTAGAGGACATTCCGGTGTCAAACTGCTTGTCAACCGTCAAGGAATGCCATCACACCATAGTGTCTTCGTTCAATTAACGCATT") }, 4, 15, 15, 115, 119, 115, 0.9663865546218487, 0.004310451608626526, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 43, "q0", 75, 118, true), ("r0", 44, 93, "q0", 0, 49, true), ("r1", 0, 218, "q1", 0, 218, false) } },
        new object[] { new[] { ("r0", "CCTATCCATTAAATAGTGGGCTGTCGGGCGTAGCTTTGGTTTGCGCAACGGCTTCTCCGAGGACGGCTCAACAAGTCACCCCCAAACCCAAGCACCATGAAGGAAACCTGCACCATGCACGATGTACGCTT"), ("r1", "ACTTCGTACGCTCCACATTCTAGAACTGCCCCCAGGTGTAGAAGAGTAAAGCCCCTCGCTTAATAAACCAGGCAACCTAATGACAAATACGGATGTGTATATCATGTATACCCACCGGAAAAGATAACGGCAAATTCGCGCGTTTACAGCTGTTTCAGCATGGTCGTCGCTGTGACCTAACTCTG") }, new[] { ("q0", "AAGCGTACATCGTGCATGGTGCAGGTTTCCTTCATGGTGCTTGGGTTTGGAGCCCGAATTGAGTTGCGCCGTGTAGGGTGACTTGTTGAGCCGTCCTCGGAGAAGCCGTTGCGCAAACCAAAGCTACGCCCGACAGCCCACTATTTAATGGATAGG"), ("q1", "ACTTCGTACGCTCCACATTCTAGAACTGCCCCCAGGTGTAGAAGAGTAAAGCCCCTCGCTTAATAAACCAGGCAACCTAATGACAAATACGGATGTGTATATCATGTATACCCACCGGAAAAGATAACGGCAAATTCGCGCGTTTACAGCTGTTTCAGCATGGTCGTCGCTGTGACCTAACTCTG") }, 4, 15, 15, 125, 126, 125, 0.9920634920634921, 0.0009980053173843624, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 80, "q0", 75, 155, true), ("r0", 81, 130, "q0", 0, 49, true), ("r1", 0, 184, "q1", 0, 184, false) } },
        new object[] { new[] { ("r0", "GCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAAATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTAGCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGATACCTCTCCATCTGACCCAAGATTGTGCTTGTTCAATTC") }, new[] { ("q0", "TTCTTAACGTGATAACAGAATCAAACCTGCCAGGCGGTCGTCGCGGACCTCGGTCGAAGTAGTGGTGCGGATCCAGGGGAACCGTTGACTCAAAAGGAGCNNNGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGC") }, 6, 10, 30, 169, 131, 35, 0.1320754716981132, 0.24254787210114034, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 39, "q0", 103, 142, false), ("r0", 0, 39, "q0", 143, 182, false) } },
        new object[] { new[] { ("r0", "GAAAAAATGGCTAGGTTCCAGCTTTTGGGGAGACGTCTTTCTGAGGGTCAGCCGTGATTCCGATTCGATTAGACTGGTCCCCACGGGTCCATGAGTACGAGGAAACTCGGTATCGAGCCTAAAAGTTATAAGGCATCTCGCCCAGGAAAGTAACGACGTATGGGTAGTTCTCCATCA"), ("r1", "CAGCTATAATGGCTAGCGCACTCTCGTTCCAGGGCGTAGTTACACTGAGCGTGCCATGTCAGCATGCTAGCGTATCGCCCCCCAATGCCCCGCAATAGGG") }, new[] { ("q0", "TAATTCGCCGACGAGTAAGCGTAGATTACACACCCAGGAAACGATCTAGACAGATTGAAATCCCCTTCATTATAGGTCGTGTAGCGCTAGACAGTCACCTNNNGAAAAAATGGCTAGGTTCCAGCTTTTGGGGAGACGTCTTTGAAAAAATGGCTAGGTTCCAGCTTTTGGGGAGACGTCTTT") }, 11, 8, 15, 257, 130, 30, 0.08403361344537816, 0.16946183007342552, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 39, "q0", 103, 142, false), ("r0", 0, 39, "q0", 143, 182, false) } },
        new object[] { new[] { ("r0", "GAAGAATCAGAGGCAAGATCnACGTGGCAGTCTCGTGTTGACGCCTTAGCCGGTGGCGAACAGTATTGACCTGGCCGATGCTAATATTCTGATTTGGGGTTGATTTGCGCTTCAGGCGCTAAAGTGGTTTTGAGTAACATGTCCTTTTGACGGGAGCAGGTCGCCTCAAGATAAGAGTAAACCTGCCTACCAAAACTTTAAGCCGGCAGAAGCTTAACT"), ("r1", "ATACCCACCGATGTGTACTCTGTTACACCGTCAGTGAGTGTAATGCTCTGGCTAGAGCCCACGCTTCCGGCTTCGTCCTCGTGCTCCAAGTACGATACCGCAAGGCAGACGCTGGTTCGCAGGTATCTGACGAGCATACTCGCTAGCCTGTGAAGAACAAGCGATTCGAGTTGTACTCTCAGCCCGCACGGTACGCCTTCCATCGGCCC") }, new[] { ("q0", "GAAGAATCAGATGCAAGATCTACGTGGCAGTCTCATGTTGACGCCTTAGCCGGTGGCGAACAGTATTGACCTGGCCGATGCTAATTTTCTGACTTGGGGTTGATTTGCGCTTCAGGCGCTAAAGTGGTTTTGAGTAACATGTCCTTTTGACGGGAGCAGGTCGCCTCAAGATACGAGTAAACCTGCCTGCCAAAACTTTAAGCCGGCAGAAGCTTAACT"), ("q1", "ATACCCATCGATGTGTACTCTGTTACACCGTCAGTGAGTGTAATGCTCTGGCTAGGGCCCACGCTTCCGGCTTCGTCCTCGTGCTCCAAGTACGATACCGCAAAGCAGACGCTGGGTCGCAGGTATCTGACGAGCAGACGCCCTAGCCTGTGAAGAACAAGCGATTTGAGTTGTACTCTCAGCCCGCACGGTACGCCTATCATCGGCCA") }, 11, 10, 15, 397, 408, 258, 0.4716636197440585, 0.040430501085181876, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 218, "q0", 0, 218, false), ("r1", 8, 197, "q1", 8, 197, false) } },
        new object[] { new[] { ("r0", "TAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG"), ("r1", "TGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTG") }, new[] { ("q0", "TAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG"), ("q1", "TGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTAACGAAGAGGCGCGCCCTCCTGAAGAG") }, 6, 15, 15, 177, 176, 165, 0.8776595744680851, 0.011229233745461815, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 97, "q0", 0, 97, false), ("r1", 0, 92, "q1", 0, 92, false) } },
        new object[] { new[] { ("r0", "AGTGAGAAGCCGTGCGTATCAATTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGCTTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTAGCTGCTAGTGTCAGACT"), ("r1", "GCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCAAAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAGCCGGTGACT") }, new[] { ("q0", "AGCTTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTAGCTGCTAGTGTCAGACTAGTGAGAAGCCGTGCGTATCAATTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCC"), ("q1", "GCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCAAAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAGCCGGTGACT") }, 6, 15, 15, 249, 249, 245, 0.9683794466403162, 0.0026991143199967677, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 78, "q0", 80, 158, false), ("r0", 79, 158, "q0", 0, 79, false), ("r1", 0, 110, "q1", 0, 110, false) } },
        new object[] { new[] { ("r0", "CACCTAACGTGAAGTTCCAAAATCCCAAACCTCTCGAGATATTTATCCAGCAAGGAGTGGCAACGCCCGCTGCTTTAATCGCTACCAAAACGCAAACAAA"), ("r1", "GCATACCCAAAAGTACACGGGTGAGGGAGGTGATATAGTACAGCTACGAAGTATCTGGCGCCTCAATAGGATTATA") }, new[] { ("q0", "CACCTAACGTGAAGTTCCAAAATCCCAAACCTCTCGAGATATTTATCCAGCAAGGAGTGGCAACGCCCGCTGCTTTAATCGCTACCAAAACGCAAACAAA"), ("q1", "GCATACCCAAAAGTACACGGGTGAGGGAGGTGATATAGTACAGCTACGAAGTATCTGGCGCCTCAATAGGATTATA") }, 6, 12, 65, 161, 161, 161, 1.0, 0.0, new (string, int, int, string, int, int, bool)[] { ("r0", 0, 99, "q0", 0, 99, false), ("r1", 0, 75, "q1", 0, 75, false) } },
    };

    [TestCaseSource(nameof(OracleCases))]
    public void CompareAssemblies_AndSyntenicBlocks_MatchOracle(
        (string, string)[] asm1, (string, string)[] asm2, int k, int minMatch, int minBlock,
        int k1, int k2, int shared, double jaccard, double mash,
        (string, int, int, string, int, int, bool)[] expectedBlocks)
    {
        var cmp = CompareAssemblies(asm1, asm2, "a", "b", k);
        Assert.That(cmp.DistinctKmers1, Is.EqualTo(k1));
        Assert.That(cmp.DistinctKmers2, Is.EqualTo(k2));
        Assert.That(cmp.SharedKmers, Is.EqualTo(shared));
        Assert.That(cmp.JaccardIndex, Is.EqualTo(jaccard).Within(Tol));
        Assert.That(cmp.MashDistance, Is.EqualTo(mash).Within(Tol));
        Assert.That(cmp.Containment1, Is.EqualTo(k1 == 0 ? 0.0 : (double)shared / k1).Within(Tol));
        Assert.That(cmp.Containment2, Is.EqualTo(k2 == 0 ? 0.0 : (double)shared / k2).Within(Tol));

        var blocks = FindSyntenicBlocks(asm1, asm2, minBlock, minMatch).ToList();
        Assert.That(blocks, Is.EqualTo(expectedBlocks));
    }

    #region CompareAssemblies (Mash, Ondov 2016)

    [Test]
    public void CompareAssemblies_MashFormula_PartialOverlap()
    {
        // canonical 4-mers {ACGT, CGTA, GTAC} vs {ACGT, CGTA, GTAA, TAAA, AAAA}: 2 shared, 6 in union.
        var cmp = CompareAssemblies(new[] { ("s", "ACGTACGT") }, new[] { ("s", "ACGTAAAA") }, "x", "y", 4);
        Assert.That(cmp, Is.EqualTo(new AssemblyComparison("x", "y", 4, 3, 5, 2, 2.0 / 3, 2.0 / 5, 1.0 / 3, Math.Log(2) / 4)));
    }

    [Test]
    public void CompareAssemblies_IsStrandInsensitive()
    {
        string s = A + B;
        var cmp = CompareAssemblies(new[] { ("s", s) }, new[] { ("s", Rc(s)) }, kmerSize: 21);
        Assert.That(cmp.JaccardIndex, Is.EqualTo(1.0));
        Assert.That(cmp.MashDistance, Is.EqualTo(0.0));
    }

    [Test]
    public void CompareAssemblies_ContigsArePooled_AndKmersDoNotSpanContigs()
    {
        // 2 contigs of one genome = 1 pooled set; the junction k-mers of the concatenation must not appear.
        var split = CompareAssemblies(new[] { ("1", A), ("2", B) }, new[] { ("g", A + B) }, kmerSize: 21);
        Assert.That(split.DistinctKmers1, Is.EqualTo(split.DistinctKmers2 - 20));
        Assert.That(split.Containment1, Is.EqualTo(1.0));
        Assert.That(split.Containment2, Is.EqualTo((double)split.SharedKmers / split.DistinctKmers2).Within(Tol));
    }

    [Test]
    public void CompareAssemblies_NonAcgtWindowsAreSkipped()
    {
        // Mash/Jellyfish: a window containing N (or any IUPAC code) is not a k-mer; case is ignored.
        var cmp = CompareAssemblies(new[] { ("s", "ACGTNACGT") }, new[] { ("s", "acgtRacgt") }, kmerSize: 4);
        Assert.That(cmp.DistinctKmers1, Is.EqualTo(1));
        Assert.That(cmp.DistinctKmers2, Is.EqualTo(1));
        Assert.That(cmp.JaccardIndex, Is.EqualTo(1.0));
    }

    [Test]
    public void CompareAssemblies_EmptySets_AreIdenticalForMash()
    {
        var cmp = CompareAssemblies(new[] { ("s", "ACG") }, new[] { ("s", "") }, kmerSize: 4);
        Assert.That(cmp.JaccardIndex, Is.EqualTo(0.0));
        Assert.That(cmp.MashDistance, Is.EqualTo(0.0)); // Mash CommandDistance: shared == union
        Assert.That(cmp.Containment1, Is.EqualTo(0.0));
    }

    [Test]
    public void CompareAssemblies_Validation()
    {
        var a = new[] { ("s", "ACGTACGT") };
        Assert.Throws<ArgumentNullException>(() => CompareAssemblies(null!, a));
        Assert.Throws<ArgumentNullException>(() => CompareAssemblies(a, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => CompareAssemblies(a, a, kmerSize: 0));
    }

    #endregion

    #region FindSyntenicBlocks (nucmer: MUMs + mgaps)

    private const string Ref = A + B + C + D;
    private const string E = "TTCCCCCAGTATCTCGTCCTCGAATGTAGATCGATCTAGCCCTCCAAACTTATACGATGCTCACTGTCGTACCTAAACGCCTCCGTCGAGCAGAAGCTTGTTTGACAGTTCGCGAGCCTGTATTGAGGTCGTGTCGTTCTGCGAGGCCGTCTGTAACAGCTCGTTGCAAGAAATGGGCATTCGTGCCTTTCGGCGTTCTTCACTAAGTAGAGAGTGCCATAAGCGCTATCCATCGGGACCGAACCCAATCGCCGAGTCGCACTACTACTACCTGGGATAGCCAGGGATTAAGCGGGATGCCAGTTCCTACCGCATAACCCGAGCTCAATTACCGTCTAGCCTCATCCCTCCGCACCCTTACGGGTCAAATAGCAGCGGAGCGACGTCGGCTAAGCCATAGGTCCAAACCCGGAATATAGTTTGGGAGCAACGGAGGCCTTAGATAAGTTTGCAGTCAGTTGAAAACTCTAAAGCCAATGCCGCTAGCTTAAGGGCGAAAGTTCGTCGTGGCTGCCGCCGACTCAGTAATATTAAACACGCACGGCCAAGACGTGATGTGCGCAGCCTTCGCGAAGCCTAAGGGCGAGAGTCCACATTTTG"; // unrelated 600-mer (oracle: no 20-mer MUM block)

    private static List<(string, int, int, string, int, int, bool)> Blocks(string r, string q, bool both = false)
        => FindSyntenicBlocks(new[] { ("r", r) }, new[] { ("q", q) }, 65, 20, uniqueInBoth: both).ToList();

    [Test]
    public void FindSyntenicBlocks_PlantedInversion_ReportsInvertedMiddleBlock()
    {
        // oracle: [0,149]/[0,149] +, [150,299]/[150,299] − (query carries RC(B)), [300,599]/[300,599] +
        Assert.That(Blocks(Ref, A + Rc(B) + C + D), Is.EqualTo(new[]
        {
            ("r", 0, 149, "q", 0, 149, false),
            ("r", 150, 299, "q", 150, 299, true),
            ("r", 300, 599, "q", 300, 599, false),
        }));
    }

    [Test]
    public void FindSyntenicBlocks_Rearrangement_TwoColinearBlocks()
    {
        Assert.That(Blocks(Ref, C + D + A + B), Is.EqualTo(new[]
        {
            ("r", 0, 299, "q", 300, 599, false),
            ("r", 300, 599, "q", 0, 299, false),
        }));
    }

    [Test]
    public void FindSyntenicBlocks_Unrelated_NoBlocks()
        => Assert.That(Blocks(Ref, E), Is.Empty);

    [Test]
    public void FindSyntenicBlocks_ReferenceRepeatIsNotAnAnchor()
    {
        // A occurs twice in the reference: only matches unique in the reference (nucmer --mumreference) anchor;
        // the unique A+B junction string still anchors [0,299], and the second copy is reached via C.
        Assert.That(Blocks(A + B + A + C, A + B + C), Is.EqualTo(new[]
        {
            ("r", 0, 299, "q", 0, 299, false),
            ("r", 449, 599, "q", 299, 449, false),
        }));
    }

    [Test]
    public void FindSyntenicBlocks_SingleMismatchSymbol_DoesNotBreakBlock()
    {
        // an N never matches, but the two MUMs around it are clustered (diagonal difference 0, gap 1 ≤ 90).
        string q = string.Concat(A.AsSpan(0, 60), "N", A.AsSpan(61), B + C + D);
        Assert.That(Blocks(Ref, q), Is.EqualTo(new[] { ("r", 0, 599, "q", 0, 599, false) }));
    }

    [Test]
    public void FindSyntenicBlocks_GapBeyondMaxGap_SplitsBlock()
    {
        // 614 unrelated bases inserted into the query: gap > nucmer -g 90 (and diagonal shift > 5) ⇒ two blocks
        // (oracle; the second starts at 299 because one base of the insert flank coincidentally extends the MUM).
        string ins = "AAAAAAA" + E + "TTTTTTT";
        var split = Blocks(Ref, A + B + ins + C + D);
        Assert.That(split, Has.Count.EqualTo(2));
        Assert.That(split[0], Is.EqualTo(("r", 0, 299, "q", 0, 299, false)));
        Assert.That(split[1], Is.EqualTo(("r", 299, 599, "q", 913, 1213, false)));
    }

    [Test]
    public void FindSyntenicBlocks_MinBlockSize_FiltersShortClusters()
    {
        var q = new[] { ("q", string.Concat(A, "GGGGG", B.AsSpan(0, 30))) };
        var all = FindSyntenicBlocks(new[] { ("r", Ref) }, q, 65, 20).ToList();
        var strict = FindSyntenicBlocks(new[] { ("r", Ref) }, q, 400, 20).ToList();
        Assert.That(all, Is.Not.Empty);
        Assert.That(strict, Is.Empty);
    }

    [Test]
    public void FindSyntenicBlocks_Validation()
    {
        var a = new[] { ("s", "ACGTACGT") };
        Assert.Throws<ArgumentNullException>(() => FindSyntenicBlocks(null!, a));
        Assert.Throws<ArgumentNullException>(() => FindSyntenicBlocks(a, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSyntenicBlocks(a, a, minBlockSize: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSyntenicBlocks(a, a, kmerSize: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSyntenicBlocks(a, a, maxGap: -1));
    }

    #endregion
}
