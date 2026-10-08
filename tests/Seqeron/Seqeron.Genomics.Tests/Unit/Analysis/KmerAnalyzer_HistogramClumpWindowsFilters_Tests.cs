// KMER-FREQ-001 / KMER-FIND-001 / KMER-UNIQUE-001 — B06 audit round 1, WP3 (F12)
// (1) GetKmerHistogram = jellyfish histo -l/-h/-i/-f; (2) FindClumpWindows (qualifying (L,t) windows);
// (3) option-aware FindUniqueKmers / FindKmersWithMinCount = jellyfish count [-C] + dump -L/-U.
// Sources: Jellyfish sub_commands/histo_main.cc + histo_main_cmdline.yaggo, dump_main.cc (gmarcais/Jellyfish master);
//          Compeau & Pevzner, Bioinformatics Algorithms ch. 1 (clump definition), Rosalind BA1E.
// Evidence: docs/algorithms/K-mer/K-mer_Frequency_Analysis.md §7.3, K-mer_Search.md §7.3,
//           Unique_And_MinCount_Kmers.md §7.3.
// Every histogram / dump expectation below is the verbatim output of the real Jellyfish 2.3.1 binary
// (`jellyfish count -m k -s 10000 -t 1 [-C]` then `histo ...` or `dump -c -L -U | sort -k2,2nr -k1,1`);
// a Python replica of histo_main.cc reproduces all 72 histo rows. Clump-window expectations come from a
// Python brute force that counts every window (agrees with an independent occurrence-interval method).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_HistogramClumpWindowsFilters_Tests
{
    private const string Ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    private const string Mixed = "acgtNNacgtacgRtTTGCAnA";
    private const string HomoMix = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGTTTGCA";
    private const string Ba1e =
        "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC";

    private const string RosalindKmer =
        "CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG" +
        "CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT" +
        "TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA" +
        "AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG" +
        "GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA" +
        "CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA" +
        "CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG";

    private static KmerCountingOptions Jellyfish(bool canonical) =>
        canonical ? new KmerCountingOptions(Canonical: true) : new KmerCountingOptions(AcgtOnly: true);

    private static string Format(IEnumerable<KmerHistogramBin> rows) =>
        string.Join(" ", rows.Select(r => $"{r.Bin} {r.Frequency}"));

    #region GetKmerHistogram — jellyfish histo (72 binary runs)

    private static IEnumerable<TestCaseData> JellyfishHistoCases()
    {
        yield return new TestCaseData(RosalindKmer, 4, true, 1L, 10000L, 1L, false, "1 23 2 34 3 27 4 25 5 7 6 5 7 4 9 3 10 2")
            .SetName("Histo_Ros_k4_C_defaults");
        yield return new TestCaseData(Ba1b, 4, true, 1L, 10000L, 1L, false, "1 16 2 2 3 1 4 1")
            .SetName("Histo_Ba1b_k4_C_defaults");
        yield return new TestCaseData(HomoMix, 3, true, 1L, 10000L, 1L, false, "1 1 2 2 6 1 8 1 30 1")
            .SetName("Histo_Homo_k3_C_defaults");
        yield return new TestCaseData(RosalindKmer, 4, true, 1L, 5L, 1L, false, "1 23 2 34 3 27 4 25 5 7 6 14")
            .SetName("Histo_Ros_k4_C_h_5");
        yield return new TestCaseData(Ba1b, 4, true, 1L, 5L, 1L, false, "1 16 2 2 3 1 4 1")
            .SetName("Histo_Ba1b_k4_C_h_5");
        yield return new TestCaseData(HomoMix, 3, true, 1L, 5L, 1L, false, "1 1 2 2 6 3")
            .SetName("Histo_Homo_k3_C_h_5");
        yield return new TestCaseData(RosalindKmer, 4, true, 1L, 10000L, 2L, false, "0 23 2 61 4 32 6 9 8 3 10 2")
            .SetName("Histo_Ros_k4_C_i_2");
        yield return new TestCaseData(Ba1b, 4, true, 1L, 10000L, 2L, false, "0 16 2 3 4 1")
            .SetName("Histo_Ba1b_k4_C_i_2");
        yield return new TestCaseData(HomoMix, 3, true, 1L, 10000L, 2L, false, "0 1 2 2 6 1 8 1 30 1")
            .SetName("Histo_Homo_k3_C_i_2");
        yield return new TestCaseData(RosalindKmer, 4, true, 3L, 8L, 2L, false, "1 57 3 52 5 12 7 4 9 5")
            .SetName("Histo_Ros_k4_C_l_3_h_8_i_2");
        yield return new TestCaseData(Ba1b, 4, true, 3L, 8L, 2L, false, "1 18 3 2")
            .SetName("Histo_Ba1b_k4_C_l_3_h_8_i_2");
        yield return new TestCaseData(HomoMix, 3, true, 3L, 8L, 2L, false, "1 3 5 1 7 1 9 1")
            .SetName("Histo_Homo_k3_C_l_3_h_8_i_2");
        yield return new TestCaseData(RosalindKmer, 4, true, 2L, 6L, 1L, false, "1 23 2 34 3 27 4 25 5 7 6 5 7 9")
            .SetName("Histo_Ros_k4_C_l_2_h_6");
        yield return new TestCaseData(Ba1b, 4, true, 2L, 6L, 1L, false, "1 16 2 2 3 1 4 1")
            .SetName("Histo_Ba1b_k4_C_l_2_h_6");
        yield return new TestCaseData(HomoMix, 3, true, 2L, 6L, 1L, false, "1 1 2 2 6 1 7 2")
            .SetName("Histo_Homo_k3_C_l_2_h_6");
        yield return new TestCaseData(RosalindKmer, 4, true, 1L, 4L, 3L, false, "0 57 3 59 6 14")
            .SetName("Histo_Ros_k4_C_l_1_h_4_i_3");
        yield return new TestCaseData(Ba1b, 4, true, 1L, 4L, 3L, false, "0 18 3 2")
            .SetName("Histo_Ba1b_k4_C_l_1_h_4_i_3");
        yield return new TestCaseData(HomoMix, 3, true, 1L, 4L, 3L, false, "0 3 6 3")
            .SetName("Histo_Homo_k3_C_l_1_h_4_i_3");
        yield return new TestCaseData(RosalindKmer, 4, true, 4L, 4L, 1L, false, "3 84 4 25 5 21")
            .SetName("Histo_Ros_k4_C_l_4_h_4");
        yield return new TestCaseData(Ba1b, 4, true, 4L, 4L, 1L, false, "3 19 4 1")
            .SetName("Histo_Ba1b_k4_C_l_4_h_4");
        yield return new TestCaseData(HomoMix, 3, true, 4L, 4L, 1L, false, "3 3 5 3")
            .SetName("Histo_Homo_k3_C_l_4_h_4");
        yield return new TestCaseData(RosalindKmer, 4, true, 1L, 5L, 1L, true, "0 0 1 23 2 34 3 27 4 25 5 7 6 14")
            .SetName("Histo_Ros_k4_C_f_h_5");
        yield return new TestCaseData(Ba1b, 4, true, 1L, 5L, 1L, true, "0 0 1 16 2 2 3 1 4 1 5 0 6 0")
            .SetName("Histo_Ba1b_k4_C_f_h_5");
        yield return new TestCaseData(HomoMix, 3, true, 1L, 5L, 1L, true, "0 0 1 1 2 2 3 0 4 0 5 0 6 3")
            .SetName("Histo_Homo_k3_C_f_h_5");
        yield return new TestCaseData(RosalindKmer, 4, true, 3L, 8L, 2L, true, "1 57 3 52 5 12 7 4 9 5")
            .SetName("Histo_Ros_k4_C_f_l_3_h_8_i_2");
        yield return new TestCaseData(Ba1b, 4, true, 3L, 8L, 2L, true, "1 18 3 2 5 0 7 0 9 0")
            .SetName("Histo_Ba1b_k4_C_f_l_3_h_8_i_2");
        yield return new TestCaseData(HomoMix, 3, true, 3L, 8L, 2L, true, "1 3 3 0 5 1 7 1 9 1")
            .SetName("Histo_Homo_k3_C_f_l_3_h_8_i_2");
        yield return new TestCaseData(RosalindKmer, 4, true, 2L, 6L, 1L, true, "1 23 2 34 3 27 4 25 5 7 6 5 7 9")
            .SetName("Histo_Ros_k4_C_f_l_2_h_6");
        yield return new TestCaseData(Ba1b, 4, true, 2L, 6L, 1L, true, "1 16 2 2 3 1 4 1 5 0 6 0 7 0")
            .SetName("Histo_Ba1b_k4_C_f_l_2_h_6");
        yield return new TestCaseData(HomoMix, 3, true, 2L, 6L, 1L, true, "1 1 2 2 3 0 4 0 5 0 6 1 7 2")
            .SetName("Histo_Homo_k3_C_f_l_2_h_6");
        yield return new TestCaseData(RosalindKmer, 4, true, 5L, 7L, 1L, false, "4 109 5 7 6 5 7 4 8 5")
            .SetName("Histo_Ros_k4_C_l_5_h_7_i_1");
        yield return new TestCaseData(Ba1b, 4, true, 5L, 7L, 1L, false, "4 20")
            .SetName("Histo_Ba1b_k4_C_l_5_h_7_i_1");
        yield return new TestCaseData(HomoMix, 3, true, 5L, 7L, 1L, false, "4 3 6 1 8 2")
            .SetName("Histo_Homo_k3_C_l_5_h_7_i_1");
        yield return new TestCaseData(RosalindKmer, 4, true, 6L, 9L, 4L, false, "2 116 6 12 10 2")
            .SetName("Histo_Ros_k4_C_l_6_h_9_i_4");
        yield return new TestCaseData(Ba1b, 4, true, 6L, 9L, 4L, false, "2 20")
            .SetName("Histo_Ba1b_k4_C_l_6_h_9_i_4");
        yield return new TestCaseData(HomoMix, 3, true, 6L, 9L, 4L, false, "2 3 6 2 10 1")
            .SetName("Histo_Homo_k3_C_l_6_h_9_i_4");
        yield return new TestCaseData(RosalindKmer, 4, false, 1L, 10000L, 1L, false, "1 91 2 65 3 36 4 6 5 9 6 1 8 1")
            .SetName("Histo_Ros_k4_plain_defaults");
        yield return new TestCaseData(Ba1b, 4, false, 1L, 10000L, 1L, false, "1 17 2 2 3 2")
            .SetName("Histo_Ba1b_k4_plain_defaults");
        yield return new TestCaseData(HomoMix, 3, false, 1L, 10000L, 1L, false, "1 6 3 2 4 2 29 1")
            .SetName("Histo_Homo_k3_plain_defaults");
        yield return new TestCaseData(RosalindKmer, 4, false, 1L, 5L, 1L, false, "1 91 2 65 3 36 4 6 5 9 6 2")
            .SetName("Histo_Ros_k4_plain_h_5");
        yield return new TestCaseData(Ba1b, 4, false, 1L, 5L, 1L, false, "1 17 2 2 3 2")
            .SetName("Histo_Ba1b_k4_plain_h_5");
        yield return new TestCaseData(HomoMix, 3, false, 1L, 5L, 1L, false, "1 6 3 2 4 2 6 1")
            .SetName("Histo_Homo_k3_plain_h_5");
        yield return new TestCaseData(RosalindKmer, 4, false, 1L, 10000L, 2L, false, "0 91 2 101 4 15 6 1 8 1")
            .SetName("Histo_Ros_k4_plain_i_2");
        yield return new TestCaseData(Ba1b, 4, false, 1L, 10000L, 2L, false, "0 17 2 4")
            .SetName("Histo_Ba1b_k4_plain_i_2");
        yield return new TestCaseData(HomoMix, 3, false, 1L, 10000L, 2L, false, "0 6 2 2 4 2 28 1")
            .SetName("Histo_Homo_k3_plain_i_2");
        yield return new TestCaseData(RosalindKmer, 4, false, 3L, 8L, 2L, false, "1 156 3 42 5 10 7 1")
            .SetName("Histo_Ros_k4_plain_l_3_h_8_i_2");
        yield return new TestCaseData(Ba1b, 4, false, 3L, 8L, 2L, false, "1 19 3 2")
            .SetName("Histo_Ba1b_k4_plain_l_3_h_8_i_2");
        yield return new TestCaseData(HomoMix, 3, false, 3L, 8L, 2L, false, "1 6 3 4 9 1")
            .SetName("Histo_Homo_k3_plain_l_3_h_8_i_2");
        yield return new TestCaseData(RosalindKmer, 4, false, 2L, 6L, 1L, false, "1 91 2 65 3 36 4 6 5 9 6 1 7 1")
            .SetName("Histo_Ros_k4_plain_l_2_h_6");
        yield return new TestCaseData(Ba1b, 4, false, 2L, 6L, 1L, false, "1 17 2 2 3 2")
            .SetName("Histo_Ba1b_k4_plain_l_2_h_6");
        yield return new TestCaseData(HomoMix, 3, false, 2L, 6L, 1L, false, "1 6 3 2 4 2 7 1")
            .SetName("Histo_Homo_k3_plain_l_2_h_6");
        yield return new TestCaseData(RosalindKmer, 4, false, 1L, 4L, 3L, false, "0 156 3 51 6 2")
            .SetName("Histo_Ros_k4_plain_l_1_h_4_i_3");
        yield return new TestCaseData(Ba1b, 4, false, 1L, 4L, 3L, false, "0 19 3 2")
            .SetName("Histo_Ba1b_k4_plain_l_1_h_4_i_3");
        yield return new TestCaseData(HomoMix, 3, false, 1L, 4L, 3L, false, "0 6 3 4 6 1")
            .SetName("Histo_Homo_k3_plain_l_1_h_4_i_3");
        yield return new TestCaseData(RosalindKmer, 4, false, 4L, 4L, 1L, false, "3 192 4 6 5 11")
            .SetName("Histo_Ros_k4_plain_l_4_h_4");
        yield return new TestCaseData(Ba1b, 4, false, 4L, 4L, 1L, false, "3 21")
            .SetName("Histo_Ba1b_k4_plain_l_4_h_4");
        yield return new TestCaseData(HomoMix, 3, false, 4L, 4L, 1L, false, "3 8 4 2 5 1")
            .SetName("Histo_Homo_k3_plain_l_4_h_4");
        yield return new TestCaseData(RosalindKmer, 4, false, 1L, 5L, 1L, true, "0 0 1 91 2 65 3 36 4 6 5 9 6 2")
            .SetName("Histo_Ros_k4_plain_f_h_5");
        yield return new TestCaseData(Ba1b, 4, false, 1L, 5L, 1L, true, "0 0 1 17 2 2 3 2 4 0 5 0 6 0")
            .SetName("Histo_Ba1b_k4_plain_f_h_5");
        yield return new TestCaseData(HomoMix, 3, false, 1L, 5L, 1L, true, "0 0 1 6 2 0 3 2 4 2 5 0 6 1")
            .SetName("Histo_Homo_k3_plain_f_h_5");
        yield return new TestCaseData(RosalindKmer, 4, false, 3L, 8L, 2L, true, "1 156 3 42 5 10 7 1 9 0")
            .SetName("Histo_Ros_k4_plain_f_l_3_h_8_i_2");
        yield return new TestCaseData(Ba1b, 4, false, 3L, 8L, 2L, true, "1 19 3 2 5 0 7 0 9 0")
            .SetName("Histo_Ba1b_k4_plain_f_l_3_h_8_i_2");
        yield return new TestCaseData(HomoMix, 3, false, 3L, 8L, 2L, true, "1 6 3 4 5 0 7 0 9 1")
            .SetName("Histo_Homo_k3_plain_f_l_3_h_8_i_2");
        yield return new TestCaseData(RosalindKmer, 4, false, 2L, 6L, 1L, true, "1 91 2 65 3 36 4 6 5 9 6 1 7 1")
            .SetName("Histo_Ros_k4_plain_f_l_2_h_6");
        yield return new TestCaseData(Ba1b, 4, false, 2L, 6L, 1L, true, "1 17 2 2 3 2 4 0 5 0 6 0 7 0")
            .SetName("Histo_Ba1b_k4_plain_f_l_2_h_6");
        yield return new TestCaseData(HomoMix, 3, false, 2L, 6L, 1L, true, "1 6 2 0 3 2 4 2 5 0 6 0 7 1")
            .SetName("Histo_Homo_k3_plain_f_l_2_h_6");
        yield return new TestCaseData(RosalindKmer, 4, false, 5L, 7L, 1L, false, "4 198 5 9 6 1 8 1")
            .SetName("Histo_Ros_k4_plain_l_5_h_7_i_1");
        yield return new TestCaseData(Ba1b, 4, false, 5L, 7L, 1L, false, "4 21")
            .SetName("Histo_Ba1b_k4_plain_l_5_h_7_i_1");
        yield return new TestCaseData(HomoMix, 3, false, 5L, 7L, 1L, false, "4 10 8 1")
            .SetName("Histo_Homo_k3_plain_l_5_h_7_i_1");
        yield return new TestCaseData(RosalindKmer, 4, false, 6L, 9L, 4L, false, "2 207 6 2")
            .SetName("Histo_Ros_k4_plain_l_6_h_9_i_4");
        yield return new TestCaseData(Ba1b, 4, false, 6L, 9L, 4L, false, "2 21")
            .SetName("Histo_Ba1b_k4_plain_l_6_h_9_i_4");
        yield return new TestCaseData(HomoMix, 3, false, 6L, 9L, 4L, false, "2 10 10 1")
            .SetName("Histo_Homo_k3_plain_l_6_h_9_i_4");
    }

    [TestCaseSource(nameof(JellyfishHistoCases))]
    public void GetKmerHistogram_MatchesJellyfishHisto(
        string sequence, int k, bool canonical, long low, long high, long increment, bool full, string expected)
    {
        var rows = KmerAnalyzer.GetKmerHistogram(sequence, k, Jellyfish(canonical), low, high, increment, full);
        Assert.That(Format(rows), Is.EqualTo(expected));
    }

    [Test]
    public void GetKmerHistogram_CountTableOverload_EqualsSequenceOverload()
    {
        var counts = KmerAnalyzer.CountKmers(RosalindKmer, 4, new KmerCountingOptions(Canonical: true));
        var fromTable = KmerAnalyzer.GetKmerHistogram(counts.Values, low: 3, high: 8, increment: 2);
        var fromSequence = KmerAnalyzer.GetKmerHistogram(RosalindKmer, 4, new KmerCountingOptions(Canonical: true), 3, 8, 2);
        Assert.That(fromTable, Is.EqualTo(fromSequence));
        Assert.That(Format(fromTable), Is.EqualTo("1 57 3 52 5 12 7 4 9 5"));
    }

    [Test]
    public void GetKmerHistogram_DefaultsWithoutCap_EqualSortedSpectrum()
    {
        // With high above every count the default histo rows are the non-zero spectrum bins (ascending).
        foreach (var (sequence, k) in new[] { (RosalindKmer, 4), (Ba1b, 3), (HomoMix, 2), (Mixed, 2) })
        {
            var spectrum = KmerAnalyzer.GetKmerSpectrum(sequence, k);
            var rows = KmerAnalyzer.GetKmerHistogram(sequence, k, high: int.MaxValue);
            Assert.That(rows.Select(r => ((int)r.Bin, (int)r.Frequency)),
                Is.EqualTo(spectrum.OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value))), $"{sequence} k={k}");
        }
    }

    [Test]
    public void GetKmerHistogram_DefaultHigh_PoolsCountsAboveTenThousandInCapBin10001()
    {
        // A^10010, k=1: one k-mer with count 10010 > high 10000 -> cap bin labelled high + increment.
        var rows = KmerAnalyzer.GetKmerHistogram(new string('A', 10010), 1);
        Assert.That(Format(rows), Is.EqualTo("10001 1"));
    }

    [Test]
    public void GetKmerHistogram_FullDefaults_Has10002BinsFromZero()
    {
        // base 0, ceil 10001, nb_buckets = 10002 (labels 0..10001); count-0 bin always empty.
        var rows = KmerAnalyzer.GetKmerHistogram("ATGATG", 3, full: true);
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(10002));
            Assert.That(rows[0], Is.EqualTo(new KmerHistogramBin(0, 0)));
            Assert.That(rows[1], Is.EqualTo(new KmerHistogramBin(1, 2)));
            Assert.That(rows[2], Is.EqualTo(new KmerHistogramBin(2, 1)));
            Assert.That(rows[^1], Is.EqualTo(new KmerHistogramBin(10001, 0)));
            Assert.That(rows.Sum(r => r.Frequency), Is.EqualTo(3));
        });
    }

    [Test]
    public void GetKmerHistogram_FrequenciesSumToDistinctKmers()
    {
        var options = new KmerCountingOptions(Canonical: true);
        int distinct = KmerAnalyzer.DistinctKmers(RosalindKmer, 5, options).Count;
        foreach (var (low, high, inc) in new[] { (1L, 10000L, 1L), (1L, 3L, 1L), (4L, 9L, 3L), (0L, 0L, 5L) })
        {
            var rows = KmerAnalyzer.GetKmerHistogram(RosalindKmer, 5, options, low, high, inc);
            Assert.That(rows.Sum(r => r.Frequency), Is.EqualTo(distinct), $"-l {low} -h {high} -i {inc}");
        }
    }

    [Test]
    public void GetKmerHistogram_EmptyInput_NoRows_FullStillListsBins()
    {
        Assert.That(KmerAnalyzer.GetKmerHistogram("", 3), Is.Empty);
        Assert.That(KmerAnalyzer.GetKmerHistogram(null!, 3), Is.Empty);
        // -l 3 -h 8 -i 2 -f with no k-mers: buckets 1,3,5,7,9 all zero (as Jellyfish -f on an empty database).
        Assert.That(Format(KmerAnalyzer.GetKmerHistogram("AC", 3, default, 3, 8, 2, full: true)), Is.EqualTo("1 0 3 0 5 0 7 0 9 0"));
    }

    [Test]
    public void GetKmerHistogram_InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 2, high: 0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>()
                .With.Property("ParamName").EqualTo("high"));
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 2, increment: 0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>()
                .With.Property("ParamName").EqualTo("increment"));
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 2, low: -1), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>()
                .With.Property("ParamName").EqualTo("low"));
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 2, high: long.MaxValue), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 2, high: long.MaxValue / 2, full: true), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>()
                .With.Property("ParamName").EqualTo("full"));
            Assert.That(() => KmerAnalyzer.GetKmerHistogram("ACGT", 0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => KmerAnalyzer.GetKmerHistogram((IEnumerable<int>)null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => KmerAnalyzer.GetKmerHistogram(new[] { 2, -1 }), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void GetKmerSpectrum_Unchanged_NoCapBin()
    {
        // The existing spectrum keeps every multiplicity (no 10000 cap), unlike histo defaults.
        var spectrum = KmerAnalyzer.GetKmerSpectrum(new string('A', 10010), 1);
        Assert.That(spectrum, Is.EquivalentTo(new Dictionary<int, int> { [10010] = 1 }));
    }

    #endregion

    #region FindClumpWindows — qualifying (L, t) windows

    private static string Format(IEnumerable<KmerClump> clumps) =>
        string.Join(" ", clumps.Select(c => c.Kmer + ":" +
            string.Join(",", c.WindowRuns.Select(r => $"{r.FirstWindowStart}-{r.LastWindowStart}"))));

    private static IEnumerable<TestCaseData> ClumpWindowCases()
    {
        // Expected = Python brute force over every window Genome[i..i+L-1] (runs of qualifying i, inclusive).
        yield return new TestCaseData(Ba1e, 5, 75, 4, "CGACA:0-6 GAAGA:0-16 AATGT:16-21").SetName("ClumpWindows_BA1E_Sample");
        yield return new TestCaseData(Ba1b, 4, 30, 3, "CATG:0-0 GCAT:0-0").SetName("ClumpWindows_BA1B_WholeSequenceWindow");
        yield return new TestCaseData(Ba1b, 4, 12, 2, "GCAT:4-5,11-12 CATG:5-6,12-13 ATGA:13-14").SetName("ClumpWindows_BA1B_L12_t2_SplitRuns");
        yield return new TestCaseData("AAAAA", 3, 5, 3, "AAA:0-0").SetName("ClumpWindows_AAAAA");
        yield return new TestCaseData("ACACGTTTTTTTTTTACACGTGGGGGGGGGGACACGT", 4, 15, 2, "TTTT:0-10 GGGG:11-22").SetName("ClumpWindows_TwoHomopolymers");
        yield return new TestCaseData("AAAAAAAAAA", 2, 4, 3, "AA:0-6").SetName("ClumpWindows_Homopolymer_EveryWindow");
        yield return new TestCaseData("acgtacgt", 4, 8, 2, "ACGT:0-0").SetName("ClumpWindows_LowerCase");
    }

    [TestCaseSource(nameof(ClumpWindowCases))]
    public void FindClumpWindows_MatchesBruteForce(string sequence, int k, int windowSize, int t, string expected)
    {
        var clumps = KmerAnalyzer.FindClumpWindows(sequence, k, windowSize, t);
        Assert.That(Format(clumps), Is.EqualTo(expected));
    }

    [TestCaseSource(nameof(ClumpWindowCases))]
    public void FindClumpWindows_SameKmerSetAsFindClumps(string sequence, int k, int windowSize, int t, string expected)
    {
        var windows = KmerAnalyzer.FindClumpWindows(sequence, k, windowSize, t).Select(c => c.Kmer);
        Assert.That(windows, Is.EquivalentTo(KmerAnalyzer.FindClumps(sequence, k, windowSize, t)));
    }

    [Test]
    public void FindClumpWindows_FirstWindowStart_IsLeftmostQualifyingWindow()
    {
        // BA1E: AATGT occurs at 21, 73, 81, 86 (0-based). All four fit in a 75-window Genome[i..i+74] iff
        // 86 + 5 - 75 <= i <= 21, i.e. i in [16, 21] (Compeau & Pevzner clump definition).
        var aatgt = KmerAnalyzer.FindClumpWindows(Ba1e, 5, 75, 4).Single(c => c.Kmer == "AATGT");
        Assert.That(KmerAnalyzer.FindKmerPositions(Ba1e, "AATGT"), Is.EqualTo(new[] { 21, 73, 81, 86 }));
        Assert.That(aatgt.FirstWindowStart, Is.EqualTo(16));
        Assert.That(aatgt.WindowRuns, Is.EqualTo(new[] { new ClumpWindowRun(16, 21) }));
    }

    [Test]
    public void FindClumpWindows_SameLeavingAndEnteringKmer_DoesNotSplitRun()
    {
        // AAAAAAAAAA, k=2, L=4, t=3: every slide removes and adds AA; one run 0..6.
        var clumps = KmerAnalyzer.FindClumpWindows("AAAAAAAAAA", 2, 4, 3);
        Assert.That(clumps.Single().WindowRuns, Is.EqualTo(new[] { new ClumpWindowRun(0, 6) }));
    }

    [TestCase(null, 3, 5, 2)]
    [TestCase("", 3, 5, 2)]
    [TestCase("ACGTACGT", 0, 5, 2)]
    [TestCase("ACGTACGT", 4, 3, 2)]
    [TestCase("ACGTACGT", 2, 9, 2)]
    [TestCase("ACGTACGT", 2, 5, 0)]
    public void FindClumpWindows_DegenerateInput_Empty(string? sequence, int k, int windowSize, int t)
    {
        Assert.That(KmerAnalyzer.FindClumpWindows(sequence!, k, windowSize, t), Is.Empty);
        Assert.That(KmerAnalyzer.FindClumps(sequence!, k, windowSize, t), Is.Empty);
    }

    #endregion

    #region Option-aware FindUniqueKmers / FindKmersWithMinCount — jellyfish count [-C] + dump -L/-U

    private static IEnumerable<TestCaseData> JellyfishDumpCases()
    {
        yield return new TestCaseData(Ba1b, 4, true, 1, 1, "AACG:1 ACAT:1 ACGT:1 AGAG:1 AGCT:1 ATCA:1 CAAC:1 CATC:1 CGAC:1 CGCA:1 CTCA:1 GACA:1 GAGA:1 GAGC:1 GCAA:1 GCGA:1")
            .SetName("Dump_Ba1b_k4_C_L1_U1");
        yield return new TestCaseData(Ba1b, 4, true, 2, int.MaxValue, "ATGC:4 CATG:3 ATGA:2 TGCA:2")
            .SetName("Dump_Ba1b_k4_C_L2");
        yield return new TestCaseData(Ba1b, 4, true, 2, 3, "CATG:3 ATGA:2 TGCA:2")
            .SetName("Dump_Ba1b_k4_C_L2_U3");
        yield return new TestCaseData(Mixed, 3, true, 1, 1, "AAA:1 CAA:1")
            .SetName("Dump_Mixed_k3_C_L1_U1");
        yield return new TestCaseData(Mixed, 3, true, 2, int.MaxValue, "ACG:5 GCA:2 GTA:2")
            .SetName("Dump_Mixed_k3_C_L2");
        yield return new TestCaseData("ATGATG", 3, true, 1, 1, "ATC:1 TCA:1")
            .SetName("Dump_ATGATG_k3_C_L1_U1");
        yield return new TestCaseData("ATGATG", 3, true, 2, int.MaxValue, "ATG:2")
            .SetName("Dump_ATGATG_k3_C_L2");
        yield return new TestCaseData(RosalindKmer, 4, true, 7, int.MaxValue, "AACT:10 ACTC:10 ACTG:9 AGTC:9 CTCA:9 AGAC:7 AGTA:7 CAGC:7 GTGA:7")
            .SetName("Dump_Ros_k4_C_L7");
        yield return new TestCaseData(RosalindKmer, 4, true, 5, 6, "ACAG:6 CACA:6 CCAG:6 CCGA:6 CGAA:6 AAAA:5 AAGT:5 CCGC:5 CGGC:5 CTGC:5 GAAA:5 GTAA:5")
            .SetName("Dump_Ros_k4_C_L5_U6");
        yield return new TestCaseData(RosalindKmer, 5, true, 4, int.MaxValue, "CTCAC:6 AACTC:5 ACTCA:5 AGACT:4 AGTAC:4 CAGTC:4")
            .SetName("Dump_RosU_k5_C_L4");
        yield return new TestCaseData(Ba1b, 4, false, 1, 1, "ACGT:1 AGAG:1 AGCT:1 ATGC:1 ATGT:1 CGCA:1 CGTT:1 GAGA:1 GAGC:1 GATG:1 GTCG:1 GTTG:1 TCGC:1 TGAG:1 TGAT:1 TGTC:1 TTGC:1")
            .SetName("Dump_Ba1b_k4_plain_L1_U1");
        yield return new TestCaseData(Ba1b, 4, false, 2, int.MaxValue, "CATG:3 GCAT:3 ATGA:2 TGCA:2")
            .SetName("Dump_Ba1b_k4_plain_L2");
        yield return new TestCaseData(Ba1b, 4, false, 2, 3, "CATG:3 GCAT:3 ATGA:2 TGCA:2")
            .SetName("Dump_Ba1b_k4_plain_L2_U3");
        yield return new TestCaseData(Mixed, 3, false, 1, 1, "GCA:1 GTA:1 TAC:1 TGC:1 TTG:1 TTT:1")
            .SetName("Dump_Mixed_k3_plain_L1_U1");
        yield return new TestCaseData(Mixed, 3, false, 2, int.MaxValue, "ACG:3 CGT:2")
            .SetName("Dump_Mixed_k3_plain_L2");
        yield return new TestCaseData("ATGATG", 3, false, 1, 1, "GAT:1 TGA:1")
            .SetName("Dump_ATGATG_k3_plain_L1_U1");
        yield return new TestCaseData("ATGATG", 3, false, 2, int.MaxValue, "ATG:2")
            .SetName("Dump_ATGATG_k3_plain_L2");
        yield return new TestCaseData(RosalindKmer, 4, false, 7, int.MaxValue, "CAGT:8")
            .SetName("Dump_Ros_k4_plain_L7");
        yield return new TestCaseData(RosalindKmer, 4, false, 5, 6, "CTCA:6 AACT:5 ACTC:5 AGTA:5 AGTC:5 AGTT:5 GAGT:5 GCTG:5 GTCT:5 TCAC:5")
            .SetName("Dump_Ros_k4_plain_L5_U6");
        yield return new TestCaseData(RosalindKmer, 5, false, 4, int.MaxValue, "CAGTC:4 CTCAC:4")
            .SetName("Dump_RosU_k5_plain_L4");
    }

    [TestCaseSource(nameof(JellyfishDumpCases))]
    public void FindKmersWithMinCount_Options_MatchesJellyfishDump(
        string sequence, int k, bool canonical, int lower, int upper, string expected)
    {
        var rows = KmerAnalyzer.FindKmersWithMinCount(sequence, k, lower, upper, Jellyfish(canonical));
        Assert.That(string.Join(" ", rows.Select(r => $"{r.Kmer}:{r.Count}")), Is.EqualTo(expected));
    }

    [TestCase(Ba1b, 4, true, "AACG ACAT ACGT AGAG AGCT ATCA CAAC CATC CGAC CGCA CTCA GACA GAGA GAGC GCAA GCGA")]
    [TestCase(Mixed, 3, true, "AAA CAA")]
    [TestCase(Mixed, 3, false, "GCA GTA TAC TGC TTG TTT")]
    [TestCase("ATGATG", 3, true, "ATC TCA")]
    public void FindUniqueKmers_Options_MatchesJellyfishDumpL1U1(string sequence, int k, bool canonical, string expected)
    {
        Assert.That(string.Join(" ", KmerAnalyzer.FindUniqueKmers(sequence, k, Jellyfish(canonical))), Is.EqualTo(expected));
    }

    [Test]
    public void FindKmersWithMinCount_ThreeArgOptions_EqualsUnboundedRange()
    {
        var options = new KmerCountingOptions(Canonical: true);
        Assert.That(KmerAnalyzer.FindKmersWithMinCount(RosalindKmer, 4, 7, options),
            Is.EqualTo(KmerAnalyzer.FindKmersWithMinCount(RosalindKmer, 4, 7, int.MaxValue, options)));
    }

    [TestCase("GTAGAGCTGTNNacgt", 2)]
    [TestCase(Ba1b, 4)]
    [TestCase(Mixed, 3)]
    public void OptionOverloads_DefaultOptions_EqualLiteralOverloads(string sequence, int k)
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindUniqueKmers(sequence, k, KmerCountingOptions.Default),
                Is.EqualTo(KmerAnalyzer.FindUniqueKmers(sequence, k)));
            Assert.That(KmerAnalyzer.FindKmersWithMinCount(sequence, k, 2, KmerCountingOptions.Default),
                Is.EqualTo(KmerAnalyzer.FindKmersWithMinCount(sequence, k, 2)));
            Assert.That(KmerAnalyzer.FindKmersWithMinCount(sequence, k, 1, 2, KmerCountingOptions.Default),
                Is.EqualTo(KmerAnalyzer.FindKmersWithMinCount(sequence, k, 1, 2)));
        });
    }

    [Test]
    public void FindUniqueKmers_Canonical_CountEqualsAnalyzeKmersSingletons()
    {
        var options = new KmerCountingOptions(Canonical: true);
        // Jellyfish 2.3.1 stats -C on the Rosalind KMER sample k=4: Unique 23 (F10).
        Assert.That(KmerAnalyzer.FindUniqueKmers(RosalindKmer, 4, options).Count(), Is.EqualTo(23));
        Assert.That(KmerAnalyzer.AnalyzeKmers(RosalindKmer, 4, options).SingletonKmers, Is.EqualTo(23));
    }

    [Test]
    public void FindKmersWithMinCount_Options_Contracts()
    {
        var options = new KmerCountingOptions(Canonical: true);
        Assert.Multiple(() =>
        {
            Assert.That(() => KmerAnalyzer.FindKmersWithMinCount("ACGT", 2, 1, -1, options),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo("maxCount"));
            Assert.That(() => KmerAnalyzer.FindUniqueKmers("ACGT", 0, options).ToList(), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(KmerAnalyzer.FindUniqueKmers(null!, 3, options), Is.Empty);
            Assert.That(KmerAnalyzer.FindKmersWithMinCount("", 3, 1, options), Is.Empty);
            Assert.That(KmerAnalyzer.FindKmersWithMinCount(Ba1b, 4, 3, 2, options), Is.Empty);
        });
    }

    #endregion
}
