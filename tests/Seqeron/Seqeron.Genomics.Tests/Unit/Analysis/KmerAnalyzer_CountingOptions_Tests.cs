// KMER-COUNT-001 / KMER-BOTH-001 / KMER-STATS-001 — option-aware counting (B06 audit round 1, WP1)
// Evidence: docs/algorithms/K-mer/K-mer_Counting.md §7.3 (Jellyfish 2.3.1 reference table)
// TestSpec: tests/TestSpecs/KMER-COUNT-001.md (ACGT-only / canonical section)
// Source: Jellyfish include/jellyfish/mer_iterator.hpp (window reset on non-ACGT code; canonical m_ < rcm_),
//         include/jellyfish/mer_dna.hpp (codes[256], get_canonical), sub_commands/count_main_cmdline.yaggo (-C);
//         Marçais G, Kingsford C (2011). Bioinformatics 27(6):764-770.
// Every expected value below was produced by the real Jellyfish 2.3.1 binary
// (`jellyfish count -m k [-C] -s 10000 -t 1` then `dump -c` / `stats` / `histo`) and agrees with an
// independent Python replica of mer_iterator; entropies by scipy.stats.entropy(counts, base=2).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_CountingOptions_Tests
{
    private static readonly KmerCountingOptions Canonical = new(Canonical: true);
    private static readonly KmerCountingOptions AcgtOnly = new(AcgtOnly: true);

    private const string Ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    private const string BothStrandRef = "GAATTCACGTTGCAGGATCCATGC";
    private const string Mixed = "acgtNNacgtacgRtTTGCAnA";

    private const string RosalindKmer =
        "CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG" +
        "CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT" +
        "TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA" +
        "AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG" +
        "GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA" +
        "CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA" +
        "CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG";

    private static Dictionary<string, int> Parse(string dump) =>
        dump.Split(' ', StringSplitOptions.RemoveEmptyEntries)
            .Select(t => t.Split(':'))
            .ToDictionary(p => p[0], p => int.Parse(p[1]));

    #region Canonical (jellyfish count -C) — full dump tables

    private static IEnumerable<TestCaseData> CanonicalDumpCases()
    {
        yield return new TestCaseData(BothStrandRef, 3,
            "AAC:1 AAT:2 ACG:2 AGG:1 ATC:2 ATG:2 CAA:1 CAC:1 CAG:1 CCA:1 GAA:2 GCA:3 GGA:2 TCA:1")
            .SetName("Canonical_BothStrandRef_k3");
        yield return new TestCaseData(BothStrandRef, 4,
            "AACG:1 AATT:1 ACGT:1 AGGA:1 ATCC:2 ATGC:1 ATGG:1 ATTC:2 CAAC:1 CACG:1 CAGG:1 CATG:1 CTGC:1 GATC:1 GCAA:1 GTGA:1 TCCA:1 TGAA:1 TGCA:1")
            .SetName("Canonical_BothStrandRef_k4");
        yield return new TestCaseData(Ba1b, 4,
            "AACG:1 ACAT:1 ACGT:1 AGAG:1 AGCT:1 ATCA:1 ATGA:2 ATGC:4 CAAC:1 CATC:1 CATG:3 CGAC:1 CGCA:1 CTCA:1 GACA:1 GAGA:1 GAGC:1 GCAA:1 GCGA:1 TGCA:2")
            .SetName("Canonical_Ba1b_k4");
        yield return new TestCaseData(Mixed, 3, "AAA:1 ACG:5 CAA:1 GCA:2 GTA:2")
            .SetName("Canonical_LowerCaseNIupac_k3");
        yield return new TestCaseData("ATGATG", 3, "ATC:1 ATG:2 TCA:1").SetName("Canonical_ATGATG_k3");
        yield return new TestCaseData("ACGTNACGT", 4, "ACGT:2").SetName("Canonical_NResetsWindow_k4");
        yield return new TestCaseData("AAUUAA", 2, "AA:2").SetName("Canonical_UIsNotAcgt_k2");
        yield return new TestCaseData("AAAA", 2, "AA:3").SetName("Canonical_Homopolymer_k2");
    }

    [TestCaseSource(nameof(CanonicalDumpCases))]
    public void CountKmers_Canonical_MatchesJellyfishDump(string sequence, int k, string expectedDump)
    {
        Assert.That(KmerAnalyzer.CountKmers(sequence, k, Canonical), Is.EquivalentTo(Parse(expectedDump)));
    }

    #endregion

    #region ACGT-only (Jellyfish window rule, no -C) — full dump tables

    private static IEnumerable<TestCaseData> AcgtOnlyDumpCases()
    {
        yield return new TestCaseData(Mixed, 3, "ACG:3 CGT:2 GCA:1 GTA:1 TAC:1 TGC:1 TTG:1 TTT:1")
            .SetName("AcgtOnly_LowerCaseNIupac_k3");
        yield return new TestCaseData("ACGTNACGT", 4, "ACGT:2").SetName("AcgtOnly_NResetsWindow_k4");
        yield return new TestCaseData("AAUUAA", 2, "AA:2").SetName("AcgtOnly_UIsNotAcgt_k2");
        yield return new TestCaseData(Ba1b, 4,
            "ACGT:1 AGAG:1 AGCT:1 ATGA:2 ATGC:1 ATGT:1 CATG:3 CGCA:1 CGTT:1 GAGA:1 GAGC:1 GATG:1 GCAT:3 GTCG:1 GTTG:1 TCGC:1 TGAG:1 TGAT:1 TGCA:2 TGTC:1 TTGC:1")
            .SetName("AcgtOnly_Ba1b_k4_EqualsLiteral");
    }

    [TestCaseSource(nameof(AcgtOnlyDumpCases))]
    public void CountKmers_AcgtOnly_MatchesJellyfishDump(string sequence, int k, string expectedDump)
    {
        Assert.That(KmerAnalyzer.CountKmers(sequence, k, AcgtOnly), Is.EquivalentTo(Parse(expectedDump)));
    }

    [Test]
    public void CountKmers_AcgtOnly_DiffersFromLiteralOnlyByNonAcgtWindows()
    {
        // Literal counting keeps CGTN, GTNA, TNAC, NACG as keys; Jellyfish drops them.
        var literal = KmerAnalyzer.CountKmers("ACGTNACGT", 4);
        Assert.That(literal.Keys, Is.EquivalentTo(new[] { "ACGT", "CGTN", "GTNA", "TNAC", "NACG" }));
        Assert.That(literal.Values.Sum(), Is.EqualTo(6));
        Assert.That(KmerAnalyzer.CountKmers("ACGTNACGT", 4, AcgtOnly).Values.Sum(), Is.EqualTo(2));
    }

    #endregion

    #region Jellyfish stats / histo rows (count [-C] + stats + histo)

    // seq, k, canonical, Unique, Distinct, Total, Max_count, histo "m:n ..."
    private static IEnumerable<TestCaseData> StatsRows()
    {
        yield return new TestCaseData(RosalindKmer, 4, false, 91, 209, 412, 8, "1:91 2:65 3:36 4:6 5:9 6:1 8:1").SetName("Stats_Rosalind_k4");
        yield return new TestCaseData(RosalindKmer, 4, true, 23, 130, 412, 10, "1:23 2:34 3:27 4:25 5:7 6:5 7:4 9:3 10:2").SetName("Stats_Rosalind_k4_C");
        yield return new TestCaseData(RosalindKmer, 5, false, 292, 348, 411, 4, "1:292 2:51 3:3 4:2").SetName("Stats_Rosalind_k5");
        yield return new TestCaseData(RosalindKmer, 5, true, 181, 279, 411, 6, "1:181 2:74 3:18 4:3 5:2 6:1").SetName("Stats_Rosalind_k5_C");
        yield return new TestCaseData(Ba1b, 4, true, 16, 20, 27, 4, "1:16 2:2 3:1 4:1").SetName("Stats_Ba1b_k4_C");
        yield return new TestCaseData(BothStrandRef, 3, true, 7, 14, 22, 3, "1:7 2:6 3:1").SetName("Stats_BothStrandRef_k3_C");
        yield return new TestCaseData(BothStrandRef, 4, true, 17, 19, 21, 2, "1:17 2:2").SetName("Stats_BothStrandRef_k4_C");
        yield return new TestCaseData(Mixed, 3, false, 6, 8, 11, 3, "1:6 2:1 3:1").SetName("Stats_Mixed_k3_AcgtOnly");
        yield return new TestCaseData(Mixed, 3, true, 2, 5, 11, 5, "1:2 2:2 5:1").SetName("Stats_Mixed_k3_C");
    }

    [TestCaseSource(nameof(StatsRows))]
    public void AnalyzeKmersAndSpectrum_WithOptions_MatchJellyfishStatsAndHisto(
        string sequence, int k, bool canonical, int unique, int distinct, int total, int maxCount, string histo)
    {
        // Non-canonical rows are run ACGT-only (Jellyfish never counts non-ACGT windows).
        var options = new KmerCountingOptions(Canonical: canonical, AcgtOnly: true);
        var stats = KmerAnalyzer.AnalyzeKmers(sequence, k, options);
        var spectrum = KmerAnalyzer.GetKmerSpectrum(sequence, k, options);
        var expectedHisto = histo.Split(' ').Select(t => t.Split(':'))
            .ToDictionary(p => int.Parse(p[0]), p => int.Parse(p[1]));

        Assert.Multiple(() =>
        {
            Assert.That(stats.SingletonKmers, Is.EqualTo(unique), "Jellyfish Unique");
            Assert.That(stats.DistinctKmers, Is.EqualTo(distinct), "Jellyfish Distinct");
            Assert.That(stats.TotalKmers, Is.EqualTo(total), "Jellyfish Total");
            Assert.That(stats.MaxCount, Is.EqualTo(maxCount), "Jellyfish Max_count");
            Assert.That(spectrum, Is.EquivalentTo(expectedHisto), "Jellyfish histo");
            Assert.That(KmerAnalyzer.DistinctKmers(sequence, k, options), Has.Count.EqualTo(distinct));
        });
    }

    [TestCase(RosalindKmer, 4, 3.169230769230769, 6.779144227048732)]
    [TestCase(Ba1b, 4, 1.35, 4.1343361131944505)]
    [TestCase(Mixed, 3, 2.2, 2.0403733936884962)]
    public void AnalyzeKmers_Canonical_MeanAndEntropy_MatchScipy(string sequence, int k, double mean, double entropyBits)
    {
        var stats = KmerAnalyzer.AnalyzeKmers(sequence, k, Canonical);
        Assert.Multiple(() =>
        {
            Assert.That(stats.AverageCount, Is.EqualTo(mean).Within(1e-15));
            Assert.That(stats.Entropy, Is.EqualTo(entropyBits).Within(1e-12));
            Assert.That(stats.MinCount, Is.EqualTo(1));
        });
    }

    [Test]
    public void AnalyzeKmers_Canonical_WithCountFilter_AppliesJellyfishLU()
    {
        // BA1B -C: retained counts with -L 2 are ATGA:2 ATGC:4 CATG:3 TGCA:2 -> Distinct 4, Total 11, Max 4, Unique 0.
        var s = KmerAnalyzer.AnalyzeKmers(Ba1b, 4, Canonical, lowerCount: 2);
        Assert.Multiple(() =>
        {
            Assert.That(s.DistinctKmers, Is.EqualTo(4));
            Assert.That(s.TotalKmers, Is.EqualTo(11));
            Assert.That(s.MaxCount, Is.EqualTo(4));
            Assert.That(s.SingletonKmers, Is.EqualTo(0));
            Assert.That(s.MinCount, Is.EqualTo(2));
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.AnalyzeKmers(Ba1b, 4, Canonical, -1));
    }

    #endregion

    #region Invariants and contracts

    [TestCase(Ba1b, 4)]
    [TestCase(Mixed, 3)]
    [TestCase(RosalindKmer, 5)]
    [TestCase("ACGTNACGT", 4)]
    public void DefaultOptions_AreIdenticalToLiteralOverloads(string sequence, int k)
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.CountKmers(sequence, k, KmerCountingOptions.Default), Is.EquivalentTo(KmerAnalyzer.CountKmers(sequence, k)));
            Assert.That(KmerAnalyzer.GetKmerSpectrum(sequence, k, default), Is.EquivalentTo(KmerAnalyzer.GetKmerSpectrum(sequence, k)));
            Assert.That(KmerAnalyzer.AnalyzeKmers(sequence, k, KmerCountingOptions.Default), Is.EqualTo(KmerAnalyzer.AnalyzeKmers(sequence, k)));
            Assert.That(KmerAnalyzer.DistinctKmers(sequence, k), Is.EquivalentTo(KmerAnalyzer.CountKmers(sequence, k).Keys));
        });
    }

    [Test]
    public void Canonical_ImpliesAcgtOnly_FlagCombinationsAgree()
    {
        Assert.That(KmerAnalyzer.CountKmers(Mixed, 3, new KmerCountingOptions(Canonical: true, AcgtOnly: true)),
            Is.EquivalentTo(KmerAnalyzer.CountKmers(Mixed, 3, Canonical)));
        Assert.That(Canonical.SkipsNonAcgt, Is.True);
        Assert.That(KmerCountingOptions.Default.SkipsNonAcgt, Is.False);
    }

    [Test]
    public void Canonical_KeysAreOrdinalMinOfReverseComplementPair_AndEqualForwardFold()
    {
        var canonical = KmerAnalyzer.CountKmers(RosalindKmer, 6, Canonical);
        var forward = KmerAnalyzer.CountKmers(RosalindKmer, 6);
        var rcCounts = KmerAnalyzer.CountKmers(DnaSequence.GetReverseComplementString(RosalindKmer), 6);
        Assert.Multiple(() =>
        {
            foreach (var (key, count) in canonical)
            {
                var rc = DnaSequence.GetReverseComplementString(key);
                Assert.That(string.CompareOrdinal(key, rc), Is.LessThanOrEqualTo(0), key);
                int expected = forward.GetValueOrDefault(key) + (rc == key ? 0 : forward.GetValueOrDefault(rc));
                Assert.That(count, Is.EqualTo(expected), key);
                // Strand symmetry: the reverse-complement strand has the same canonical table.
                Assert.That(count, Is.EqualTo(rcCounts.GetValueOrDefault(key) + (rc == key ? 0 : rcCounts.GetValueOrDefault(rc))), key);
            }
            Assert.That(canonical.Values.Sum(), Is.EqualTo(RosalindKmer.Length - 6 + 1));
        });
    }

    [Test]
    public void Canonical_IsInvariantUnderReverseComplementOfInput()
    {
        var rc = DnaSequence.GetReverseComplementString(Ba1b);
        Assert.That(KmerAnalyzer.CountKmers(rc, 4, Canonical), Is.EquivalentTo(KmerAnalyzer.CountKmers(Ba1b, 4, Canonical)));
    }

    [Test]
    public void OptionOverloads_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.CountKmers(null!, 0, Canonical), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmers("", 3, AcgtOnly), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmers("ACG", 4, Canonical), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmers("NNNNNN", 2, AcgtOnly), Is.Empty);
            Assert.That(KmerAnalyzer.DistinctKmers(null!, 3), Is.Empty);
            Assert.That(KmerAnalyzer.AnalyzeKmers("NNNN", 2, Canonical), Is.EqualTo(new KmerStatistics(0, 0, 0, 0, 0, 0)));
            Assert.That(KmerAnalyzer.GetKmerSpectrum("NNNN", 2, AcgtOnly), Is.Empty);
            var ex = Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmers("ACGT", 0, Canonical));
            Assert.That(ex!.ParamName, Is.EqualTo("k"));
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.DistinctKmers("ACGT", -1, AcgtOnly));
        });
    }

    [Test]
    public void OptionOverload_HonoursCancellationAndProgress()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => KmerAnalyzer.CountKmers(RosalindKmer, 4, Canonical, cts.Token));

        var reports = new List<double>();
        var progress = new SynchronousProgress(reports.Add);
        KmerAnalyzer.CountKmers(RosalindKmer, 4, AcgtOnly, CancellationToken.None, progress);
        Assert.That(reports, Is.EqualTo(new[] { 0.0, 1.0 }));
    }

    [Test]
    public void DistinctKmers_ReturnsOrdinalCallerOwnedSet()
    {
        var set = KmerAnalyzer.DistinctKmers("atgatg", 3, Canonical);
        Assert.That(set, Is.EquivalentTo(new[] { "ATC", "ATG", "TCA" }));
        Assert.That(set.Comparer, Is.EqualTo(StringComparer.Ordinal));
        set.Add("XXX");
        Assert.That(KmerAnalyzer.DistinctKmers("atgatg", 3, Canonical), Has.Count.EqualTo(3));
    }

    private sealed class SynchronousProgress(Action<double> onReport) : IProgress<double>
    {
        public void Report(double value) => onReport(value);
    }

    #endregion
}
