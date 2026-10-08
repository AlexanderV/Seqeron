// KMER-BOTH-001 / KMER-FIND-001 / KMER-FREQ-001 / KMER-DIST-001 (audit round 2, WP8) — option-aware both-strand counting
// (kPAL ACGT-only), most-frequent k-mers and frequency profiles (Jellyfish -C / ACGT-only), and the two spaced
// conventions: words with a non-ACGT symbol at a match position are dropped, and the default reverse-complement mode.
// Evidence: docs/Evidence/KMER-BOTH-001-Evidence.md, KMER-FIND-001-Evidence.md, KMER-FREQ-001-Evidence.md,
// KMER-DIST-001-Evidence.md; algorithm docs Both_Strand_Kmer_Counting.md, K-mer_Search.md, K-mer_Frequency_Analysis.md,
// K-mer_Euclidean_Distance.md §7.7.
// Reference values (executed):
//   - kPAL (LUMC/kPAL master, kpal/klib.py + metrics.py run from source): Profile.from_sequences([s], k) splits s on
//     [^AaCcGgTt] and counts each part; balance() adds each k-mer's count to its reverse complement's (palindromes doubled).
//   - Jellyfish 2.3.1: jellyfish count -m k -s 10000 [-C] + jellyfish dump -c (frequencies = count / sum of counts;
//     most frequent = arg-max). The kPAL forward profile equals the Jellyfish ACGT-only dump on every input.
//   - spaced 1.2.0 (Ubuntu archive binary; src/sort.h spacedDNA): spaced [-r] -t 1 -f <patterns> -d JS|EU, 12 printed
//     digits. Locked values are the Python replica (rep.py: non-ACGT letters -> N, word dropped when a match position
//     reads N, frequency = count / (L - l + 1); both-strand mode = first record forward + reverse-complement counts,
//     total 2 W, vs the second record forward), which equals all 36 printed spaced values.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests
{
    private const double Tol = 1e-12;

    private static readonly KmerCountingOptions AcgtOnly = new(AcgtOnly: true);
    private static readonly KmerCountingOptions Canonical = new(Canonical: true);

    private const string X = "GAATTCNNACGTTGCAGGATCCATGCRYacgtgcaNTTGCA";
    private const string Ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";

    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
    private const string A = "CAGAGCAGACAACTAAGTGCTATCAACTAGGCGAAAGCCGCCTGAGGTGCTACTACAGTGTCGGGCTTCGAGGTTCCGACAAATAACTACGTTTCCCTTAAAACGCACCGTGGGGATTTTCCAGAAAACAAAACATGACCCATTCCAACCTAATGAGAGTTGCCATGTCATTTTTGAGCAACCGTCTTTCAAAAGTTCGTAGTTTCATCTTCCTGAAATTATTGCGAATGGTCTTACCAGCCATCCGTCTGGGACAATACGGCGGTAGACTGGTTCGGGGCGCTGTTCGTTAACAAAGTT";
    private const string B = "ATTATTATCGTTAGGCACTAGGACTGATTTAAATAACTTAATGTAATCATTTATATGCGTAAGAAAAACCTCGACAAATGAGTATGGTCGTTCTTTCCAAATAAATTTAACGAGATTGATCCACGGAGAAACAACATGTTTTCCTATTATTTTTAGGTTTTTTTAAACATGTTTGATTGGTATAAATATTGCGTTTTAAATCGAAATTTCGGGATATTAGTTGAATATAGTATGTTATGAAGGCACCTAT";
    // S1 / S2 with N, IUPAC (R, Y) and lower case inserted.
    private const string N1 = "AGGTAAGGTGNGTTGAGATctggacTTTTGACGCCTRGAGCCCGCAGTGCTCCTCGAAAAGTAGCNNATGCCTTGGGCTGCT";
    private const string N2 = "CAAAGGCCCTACCTTCTTATAGTCCTTYCAACATACAAGTAtagttgGAAGTTCTAAGTTCAGNTTAATC";

    private static readonly string[] P5 = ["11011", "10111", "11101"];
    private static readonly string[] P7 = ["1101011", "1011101", "1110011"];
    private static readonly string[] P4 = ["1111"];

    private static Dictionary<string, double> Parse(string table) =>
        table.Split(';').Select(e => e.Split('=')).ToDictionary(p => p[0], p => double.Parse(p[1], System.Globalization.CultureInfo.InvariantCulture));

    private static void AssertTable(IReadOnlyDictionary<string, int> actual, string expected)
    {
        var e = Parse(expected);
        Assert.That(actual.Keys, Is.EquivalentTo(e.Keys));
        foreach (var (k, v) in e)
            Assert.That(actual[k], Is.EqualTo((int)v), k);
    }

    #region CountKmersBothStrands — kPAL ACGT-only balance

    [TestCase(X, 3, "AAC=1;AAT=2;ACG=4;AGG=1;ATC=2;ATG=2;ATT=2;CAA=2;CAC=1;CAG=1;CAT=2;CCA=1;CCT=1;CGT=4;CTG=1;GAA=2;GAT=2;GCA=7;GGA=2;GTG=1;GTT=1;TCC=2;TGC=7;TGG=1;TTC=2;TTG=2", 56)]
    [TestCase(X, 4, "AACG=1;AATT=2;ACGT=4;AGGA=1;ATCC=2;ATGC=1;ATGG=1;ATTC=2;CAAC=1;CACG=1;CAGG=1;CATG=2;CCAT=1;CCTG=1;CGTG=1;CGTT=1;CTGC=1;GAAT=2;GATC=2;GCAA=2;GCAC=1;GCAG=1;GCAT=1;GGAT=2;GTGC=1;GTTG=1;TCCA=1;TCCT=1;TGCA=6;TGGA=1;TTGC=2", 48)]
    [TestCase("ACGTNACGTAAcgtRTT", 3, "AAC=1;ACG=6;CGT=6;GTA=1;GTT=1;TAA=1;TAC=1;TTA=1", 18)]
    [TestCase("AAAANTTTTGGGGuCCCC", 2, "AA=6;CA=1;CC=6;GG=6;TG=1;TT=6", 26)]
    public void BothStrands_AcgtOnly_EqualsKpalBalancedProfile(string sequence, int k, string expected, int total)
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(sequence, k, AcgtOnly);
        AssertTable(counts, expected);
        Assert.That(counts.Values.Sum(), Is.EqualTo(total));
    }

    [Test]
    public void BothStrands_AcgtOnly_PalindromesDoubled_TotalTwiceAcgtWindows()
    {
        var counts = KmerAnalyzer.CountKmersBothStrands(X, 4, AcgtOnly);
        int acgtWindows = KmerAnalyzer.CountKmers(X, 4, AcgtOnly).Values.Sum();
        Assert.Multiple(() =>
        {
            Assert.That(counts.Values.Sum(), Is.EqualTo(2 * acgtWindows));
            Assert.That(counts["AATT"], Is.EqualTo(2)); // kPAL counts[i] += counts[i] (Jellyfish -C: 1)
            Assert.That(KmerAnalyzer.CountKmers(X, 4, Canonical)["AATT"], Is.EqualTo(1));
            foreach (var (w, c) in counts)
                Assert.That(counts[DnaSequence.GetReverseComplementString(w)], Is.EqualTo(c), w);
        });
    }

    [Test]
    public void BothStrands_DefaultOptions_EqualsLegacyOverload()
    {
        foreach (var s in new[] { X, Ba1b, "acgtNNacgtacgRtTTGCA", "" })
            Assert.That(KmerAnalyzer.CountKmersBothStrands(s, 3, KmerCountingOptions.Default),
                Is.EquivalentTo(KmerAnalyzer.CountKmersBothStrands(s, 3)), s);
    }

    [Test]
    public void BothStrands_AcgtOnlyOnAcgtInput_EqualsLiteral()
        => Assert.That(KmerAnalyzer.CountKmersBothStrands(Ba1b, 4, AcgtOnly), Is.EquivalentTo(KmerAnalyzer.CountKmersBothStrands(Ba1b, 4)));

    [Test]
    public void BothStrands_Canonical_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() => KmerAnalyzer.CountKmersBothStrands(X, 3, Canonical));
        Assert.That(ex!.ParamName, Is.EqualTo("options"));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.CountKmersBothStrands(X, 3, new KmerCountingOptions(true, true)));
    }

    [Test]
    public void BothStrands_Options_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.CountKmersBothStrands(null!, 3, AcgtOnly), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmersBothStrands("NNNN", 2, AcgtOnly), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmersBothStrands("ACG", 4, AcgtOnly), Is.Empty);
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmersBothStrands("ACGT", 0, AcgtOnly));
        });
    }

    #endregion

    #region FindMostFrequentKmers — Jellyfish -C / ACGT-only

    [TestCase(Ba1b, 4, true, false, new[] { "ATGC" })]
    [TestCase(Ba1b, 4, false, true, new[] { "CATG", "GCAT" })]
    [TestCase(X, 2, true, false, new[] { "CA" })]
    [TestCase(X, 2, false, true, new[] { "CA", "GC", "TG" })]
    [TestCase(X, 3, true, false, new[] { "GCA" })]
    [TestCase(X, 3, false, true, new[] { "TGC" })]
    [TestCase(X, 4, true, false, new[] { "TGCA" })]
    [TestCase("AAAANTTTTGGGGuCCCC", 2, true, false, new[] { "AA", "CC" })]
    [TestCase("AAAANTTTTGGGGuCCCC", 2, false, true, new[] { "AA", "CC", "GG", "TT" })]
    [TestCase("ACGTNACGTAAcgtRTT", 3, true, false, new[] { "ACG" })]
    public void MostFrequent_Options_EqualJellyfishDumpArgMax(string sequence, int k, bool canonical, bool acgtOnly, string[] expected)
        => Assert.That(KmerAnalyzer.FindMostFrequentKmers(sequence, k, new KmerCountingOptions(canonical, acgtOnly)), Is.EquivalentTo(expected));

    [Test]
    public void MostFrequent_DefaultOptions_EqualsLegacy_Ba1bSample()
    {
        Assert.That(KmerAnalyzer.FindMostFrequentKmers(Ba1b, 4, KmerCountingOptions.Default), Is.EquivalentTo(new[] { "CATG", "GCAT" }));
        Assert.That(KmerAnalyzer.FindMostFrequentKmers(X, 2, KmerCountingOptions.Default), Is.EquivalentTo(KmerAnalyzer.FindMostFrequentKmers(X, 2)));
    }

    [Test]
    public void MostFrequent_Options_EdgeCases()
    {
        Assert.That(KmerAnalyzer.FindMostFrequentKmers("NNNN", 2, Canonical), Is.Empty);
        Assert.That(KmerAnalyzer.FindMostFrequentKmers("", 2, Canonical), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.FindMostFrequentKmers("ACGT", 0, Canonical).ToList());
    }

    #endregion

    #region GetKmerFrequencies — kPAL / Jellyfish profiles

    [TestCase(X, 2, true, false, "AA=0.125;AC=0.125;AG=0.03125;AT=0.09375;CA=0.25;CC=0.0625;CG=0.0625;GA=0.125;GC=0.125")]
    [TestCase(X, 2, false, true, "AA=0.03125;AC=0.0625;AG=0.03125;AT=0.09375;CA=0.125;CC=0.03125;CG=0.0625;GA=0.0625;GC=0.125;GG=0.03125;GT=0.0625;TC=0.0625;TG=0.125;TT=0.09375")]
    [TestCase(X, 4, true, false, "AACG=0.041666666666666664;AATT=0.041666666666666664;ACGT=0.08333333333333333;AGGA=0.041666666666666664;ATCC=0.08333333333333333;ATGC=0.041666666666666664;ATGG=0.041666666666666664;ATTC=0.08333333333333333;CAAC=0.041666666666666664;CACG=0.041666666666666664;CAGG=0.041666666666666664;CATG=0.041666666666666664;CTGC=0.041666666666666664;GATC=0.041666666666666664;GCAA=0.08333333333333333;GCAC=0.041666666666666664;TCCA=0.041666666666666664;TGCA=0.125")]
    [TestCase(X, 4, false, true, "AATT=0.041666666666666664;ACGT=0.08333333333333333;AGGA=0.041666666666666664;ATCC=0.041666666666666664;ATGC=0.041666666666666664;ATTC=0.041666666666666664;CAGG=0.041666666666666664;CATG=0.041666666666666664;CCAT=0.041666666666666664;CGTG=0.041666666666666664;CGTT=0.041666666666666664;GAAT=0.041666666666666664;GATC=0.041666666666666664;GCAG=0.041666666666666664;GGAT=0.041666666666666664;GTGC=0.041666666666666664;GTTG=0.041666666666666664;TCCA=0.041666666666666664;TGCA=0.125;TTGC=0.08333333333333333")]
    [TestCase("ACGTNACGTAAcgtRTT", 3, true, false, "AAC=0.1111111111111111;ACG=0.6666666666666666;GTA=0.1111111111111111;TAA=0.1111111111111111")]
    [TestCase("AAAANTTTTGGGGuCCCC", 3, false, true, "AAA=0.2;CCC=0.2;GGG=0.2;TGG=0.1;TTG=0.1;TTT=0.2")]
    public void Frequencies_Options_EqualJellyfishAndKpalProfiles(string sequence, int k, bool canonical, bool acgtOnly, string expected)
    {
        var freq = KmerAnalyzer.GetKmerFrequencies(sequence, k, new KmerCountingOptions(canonical, acgtOnly));
        var e = Parse(expected);
        Assert.That(freq.Keys, Is.EquivalentTo(e.Keys));
        foreach (var (w, f) in e)
            Assert.That(freq[w], Is.EqualTo(f).Within(Tol), w);
        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(1e-12));
    }

    [Test]
    public void Frequencies_DefaultOptions_EqualsLegacy()
        => Assert.That(KmerAnalyzer.GetKmerFrequencies(X, 3, KmerCountingOptions.Default), Is.EquivalentTo(KmerAnalyzer.GetKmerFrequencies(X, 3)));

    [Test]
    public void Frequencies_Options_EdgeCases()
    {
        Assert.That(KmerAnalyzer.GetKmerFrequencies("NNNNN", 2, AcgtOnly), Is.Empty);
        Assert.That(KmerAnalyzer.GetKmerFrequencies(null!, 2, Canonical), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.GetKmerFrequencies("ACGT", -1, AcgtOnly));
    }

    #endregion

    #region CountSpacedWords — spaced word rule

    [Test]
    public void SpacedWords_AcgtOnly_DropsOnlyMatchPositionNonAcgt()
    {
        // Pattern 1011: windows ANGT (N at the don't-care position -> word AGT kept), NGTC (dropped), GTCN (dropped).
        var counts = KmerAnalyzer.CountSpacedWords("ANGTCN", "1011", AcgtOnly);
        Assert.That(counts, Is.EquivalentTo(new Dictionary<string, int> { ["AGT"] = 1 }));
        var literal = KmerAnalyzer.CountSpacedWords("ANGTCN", "1011", KmerCountingOptions.Default);
        Assert.That(literal, Is.EquivalentTo(new Dictionary<string, int> { ["AGT"] = 1, ["NTC"] = 1, ["GCN"] = 1 }));
    }

    [Test]
    public void SpacedWords_AcgtOnlyAllOnesPattern_EqualsKmerAcgtOnly()
    {
        foreach (var s in new[] { X, N1, N2, "acgtRYacgu" })
            Assert.That(KmerAnalyzer.CountSpacedWords(s, "1111", AcgtOnly), Is.EquivalentTo(KmerAnalyzer.CountKmers(s, 4, AcgtOnly)), s);
    }

    [Test]
    public void SpacedWords_DefaultOptions_EqualsLegacy()
        => Assert.That(KmerAnalyzer.CountSpacedWords(N1, "11011", KmerCountingOptions.Default), Is.EquivalentTo(KmerAnalyzer.CountSpacedWords(N1, "11011")));

    [Test]
    public void SpacedWords_Canonical_Throws()
    {
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.CountSpacedWords(S1, "1101", Canonical));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, P5, KmerDistanceMetric.JensenShannon, Canonical));
    }

    #endregion

    #region SpacedWordDistance — spaced 1.2.0 with / without -r

    private static IEnumerable<TestCaseData> SpacedCases()
    {
        yield return new TestCaseData(S1, S2, P5, KmerDistanceMetric.JensenShannon, false, 0.8163228541607374, 0.816322854161).SetName("Spaced_S1S2_P5_JS_r");
        yield return new TestCaseData(S1, S2, P5, KmerDistanceMetric.JensenShannon, true, 0.7417483114521343, 0.741748311452).SetName("Spaced_S1S2_P5_JS_both");
        yield return new TestCaseData(S1, S2, P5, KmerDistanceMetric.EuclideanCounts, false, 12.40897581662776, 12.4089758166).SetName("Spaced_S1S2_P5_EU_r");
        yield return new TestCaseData(S1, S2, P5, KmerDistanceMetric.EuclideanCounts, true, 16.165781203024796, 16.165781203).SetName("Spaced_S1S2_P5_EU_both");
        yield return new TestCaseData(S1, S2, P7, KmerDistanceMetric.JensenShannon, false, 0.9467819156191135, 0.946781915619).SetName("Spaced_S1S2_P7_JS_r");
        yield return new TestCaseData(S1, S2, P7, KmerDistanceMetric.JensenShannon, true, 0.8999302969511375, 0.899930296951).SetName("Spaced_S1S2_P7_JS_both");
        yield return new TestCaseData(S1, S2, P7, KmerDistanceMetric.EuclideanCounts, false, 11.660417980553781, 11.6604179806).SetName("Spaced_S1S2_P7_EU_r");
        yield return new TestCaseData(S1, S2, P7, KmerDistanceMetric.EuclideanCounts, true, 14.557157018207116, 14.5571570182).SetName("Spaced_S1S2_P7_EU_both");
        yield return new TestCaseData(S1, S2, P4, KmerDistanceMetric.JensenShannon, false, 0.8003989533664194, 0.800398953366).SetName("Spaced_S1S2_P4_JS_r");
        yield return new TestCaseData(S1, S2, P4, KmerDistanceMetric.JensenShannon, true, 0.7053854348362533, 0.705385434836).SetName("Spaced_S1S2_P4_JS_both");
        yield return new TestCaseData(S1, S2, P4, KmerDistanceMetric.EuclideanCounts, false, 12.489995996796797, 12.4899959968).SetName("Spaced_S1S2_P4_EU_r");
        yield return new TestCaseData(S1, S2, P4, KmerDistanceMetric.EuclideanCounts, true, 15.7797338380595, 15.7797338381).SetName("Spaced_S1S2_P4_EU_both");
        yield return new TestCaseData(A, B, P5, KmerDistanceMetric.JensenShannon, false, 0.4145172554123738, 0.414517255412).SetName("Spaced_AB_P5_JS_r");
        yield return new TestCaseData(A, B, P5, KmerDistanceMetric.JensenShannon, true, 0.3475325682040406, 0.347532568204).SetName("Spaced_AB_P5_JS_both");
        yield return new TestCaseData(A, B, P5, KmerDistanceMetric.EuclideanCounts, false, 25.704410975393987, 25.7044109754).SetName("Spaced_AB_P5_EU_r");
        yield return new TestCaseData(A, B, P5, KmerDistanceMetric.EuclideanCounts, true, 35.63491330667269, 35.6349133067).SetName("Spaced_AB_P5_EU_both");
        yield return new TestCaseData(A, B, P7, KmerDistanceMetric.JensenShannon, false, 0.7579794284965504, 0.757979428497).SetName("Spaced_AB_P7_JS_r");
        yield return new TestCaseData(A, B, P7, KmerDistanceMetric.JensenShannon, true, 0.6723119547618994, 0.672311954762).SetName("Spaced_AB_P7_JS_both");
        yield return new TestCaseData(A, B, P7, KmerDistanceMetric.EuclideanCounts, false, 23.845599914433468, 23.8455999144).SetName("Spaced_AB_P7_EU_r");
        yield return new TestCaseData(A, B, P7, KmerDistanceMetric.EuclideanCounts, true, 29.602587870482225, 29.6025878705).SetName("Spaced_AB_P7_EU_both");
        yield return new TestCaseData(A, B, P4, KmerDistanceMetric.JensenShannon, false, 0.418942116524784, 0.418942116525).SetName("Spaced_AB_P4_JS_r");
        yield return new TestCaseData(A, B, P4, KmerDistanceMetric.JensenShannon, true, 0.3556317547108073, 0.355631754711).SetName("Spaced_AB_P4_JS_both");
        yield return new TestCaseData(A, B, P4, KmerDistanceMetric.EuclideanCounts, false, 25.80697580112788, 25.8069758011).SetName("Spaced_AB_P4_EU_r");
        yield return new TestCaseData(A, B, P4, KmerDistanceMetric.EuclideanCounts, true, 36.61966684720111, 36.6196668472).SetName("Spaced_AB_P4_EU_both");
        yield return new TestCaseData(N1, N2, P5, KmerDistanceMetric.JensenShannon, false, 0.7061373678517501, 0.706137367852).SetName("Spaced_N1N2_P5_JS_r");
        yield return new TestCaseData(N1, N2, P5, KmerDistanceMetric.JensenShannon, true, 0.6447798019943919, 0.644779801994).SetName("Spaced_N1N2_P5_JS_both");
        yield return new TestCaseData(N1, N2, P5, KmerDistanceMetric.EuclideanCounts, false, 11.282731603082693, 11.2827316031).SetName("Spaced_N1N2_P5_EU_r");
        yield return new TestCaseData(N1, N2, P5, KmerDistanceMetric.EuclideanCounts, true, 14.444165136064031, 14.4441651361).SetName("Spaced_N1N2_P5_EU_both");
        yield return new TestCaseData(N1, N2, P7, KmerDistanceMetric.JensenShannon, false, 0.7674733740815286, 0.767473374082).SetName("Spaced_N1N2_P7_JS_r");
        yield return new TestCaseData(N1, N2, P7, KmerDistanceMetric.JensenShannon, true, 0.7249059828308658, 0.724905982831).SetName("Spaced_N1N2_P7_JS_both");
        yield return new TestCaseData(N1, N2, P7, KmerDistanceMetric.EuclideanCounts, false, 10.59658835000654, 10.59658835).SetName("Spaced_N1N2_P7_EU_r");
        yield return new TestCaseData(N1, N2, P7, KmerDistanceMetric.EuclideanCounts, true, 12.829657257375308, 12.8296572574).SetName("Spaced_N1N2_P7_EU_both");
        yield return new TestCaseData(N1, N2, P4, KmerDistanceMetric.JensenShannon, false, 0.6939755479919848, 0.693975547992).SetName("Spaced_N1N2_P4_JS_r");
        yield return new TestCaseData(N1, N2, P4, KmerDistanceMetric.JensenShannon, true, 0.6181605928211007, 0.618160592821).SetName("Spaced_N1N2_P4_JS_both");
        yield return new TestCaseData(N1, N2, P4, KmerDistanceMetric.EuclideanCounts, false, 11.357816691600547, 11.3578166916).SetName("Spaced_N1N2_P4_EU_r");
        yield return new TestCaseData(N1, N2, P4, KmerDistanceMetric.EuclideanCounts, true, 14.247806848775006, 14.2478068488).SetName("Spaced_N1N2_P4_EU_both");
    }

    [TestCaseSource(nameof(SpacedCases))]
    public void SpacedWordDistance_AcgtOnly_ReproducesSpacedBinary(string seq1, string seq2, string[] patterns, KmerDistanceMetric metric, bool bothStrands, double replica, double spacedPrinted)
    {
        double d = KmerAnalyzer.SpacedWordDistance(seq1, seq2, patterns, metric, AcgtOnly, bothStrands);
        Assert.That(d, Is.EqualTo(replica).Within(metric == KmerDistanceMetric.EuclideanCounts ? 1e-10 : Tol));
        Assert.That(d, Is.EqualTo(spacedPrinted).Within(Math.Abs(spacedPrinted) * 1e-11)); // 12 significant digits
    }

    [Test]
    public void SpacedWordDistance_BothStrands_DependsOnArgumentOrder_AsSpacedInputOrder()
    {
        // spaced (no -r) on the file N2, N1: JS 0.632325490347, EU 14.604611969 (vs 0.644779801994 / 14.4441651361 for N1, N2).
        Assert.That(KmerAnalyzer.SpacedWordDistance(N2, N1, P5, KmerDistanceMetric.JensenShannon, AcgtOnly, true),
            Is.EqualTo(0.6323254903466718).Within(Tol));
        Assert.That(KmerAnalyzer.SpacedWordDistance(N2, N1, P5, KmerDistanceMetric.EuclideanCounts, AcgtOnly, true),
            Is.EqualTo(14.604611968975481).Within(1e-10));
    }

    [Test]
    public void SpacedWordDistance_DefaultOptionsSingleStrand_EqualsLegacyOverload()
    {
        foreach (var metric in new[] { KmerDistanceMetric.Euclidean, KmerDistanceMetric.JensenShannon, KmerDistanceMetric.EuclideanCounts, KmerDistanceMetric.Cosine })
        {
            Assert.That(KmerAnalyzer.SpacedWordDistance(N1, N2, P5, metric, KmerCountingOptions.Default),
                Is.EqualTo(KmerAnalyzer.SpacedWordDistance(N1, N2, P5, metric)), metric.ToString());
            Assert.That(KmerAnalyzer.SpacedWordDistance(S1, S2, P7, metric, AcgtOnly),
                Is.EqualTo(KmerAnalyzer.SpacedWordDistance(S1, S2, P7, metric)), metric.ToString()); // ACGT input
        }
    }

    [Test]
    public void SpacedWordDistance_LiteralWordsOnNInput_DifferFromSpacedRule()
        => Assert.That(KmerAnalyzer.SpacedWordDistance(N1, N2, P5, KmerDistanceMetric.JensenShannon),
            Is.Not.EqualTo(0.7061373678517503).Within(1e-6));

    [Test]
    public void SpacedWordDistance_BothStrands_Validation()
    {
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, P5, KmerDistanceMetric.D2Star, AcgtOnly, true));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, ["1101", "11"], KmerDistanceMetric.Euclidean, AcgtOnly, true));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, ["0110"], KmerDistanceMetric.Euclidean, AcgtOnly, true));
        Assert.That(KmerAnalyzer.SpacedWordDistance(null!, null!, P5, KmerDistanceMetric.JensenShannon, AcgtOnly, true), Is.EqualTo(0));
    }

    #endregion
}
