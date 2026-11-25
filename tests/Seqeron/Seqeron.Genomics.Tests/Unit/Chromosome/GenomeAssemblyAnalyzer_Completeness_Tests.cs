using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-COMPLETE-001: BUSCO-like marker completeness (compleasm classification, BUSCO percentages)
/// and k-mer completeness (Merqury solid threshold / completeness / QV, GenomeScope error rate and size).
/// Expected values come from tools/review/oracles/B10/completeness_oracle.py (Biopython 1.88 local
/// alignment + transcribed compleasm v0.2.6 / BUSCO hmmer.py / Merqury / GenomeScope code); see
/// docs/Evidence/ASSEMBLY-COMPLETE-001-Evidence.md.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_Completeness_Tests
{
    private const double Tol = 1e-12;

    // Oracle scenario components (random, seed 41 in the oracle).
    private const string G = "TGCCTGGTACATCCGCGAAATGCAGTCAAACCAGGCACTAGTCAATATTCAACATTTGAATGCCAAACTATACCTGTTTAGCCGTGTGTACCCACCTCACCTGACTATCTTGGGACCCCTAAGCTGGCGTCAGTGTTCATTCTCTGACTTTAACTTGGCTACCTCATACGACTCAGCTAACGCGTATAGAGTAGTGCAAGTGCTCGCAAAGGGGACATAACGGGCTTTGGCCAGACGGGGGGGAAGGCTAGCGTTTCACAGGACTTATAATTAAAAGTTTGACATTGTCGACGCCGGTTA";
    private const string L = "CCCAAAGGCCCGCCGATAAAAGCACGGCACGACACTCCGTCGAGTTGTTGTGCAATCGTTTCTGCAGAAGGGCGCATTATCATCTCAATCCTGTCTTGGT";
    private const string R = "ATGCGCCAGTATATTGCTCTTACTAACCGAAGAGTGGGGCACAATTAATCACTTAAACTCCCATCCAACAAATGAGGACGGTATTGCGTAAATTTGCGAG";
    private const string F50a = "ATCGCCTCCTCAACGCCTAGGCATCGTGTAATCGAAGAGTGGCCGAGCCA";
    private const string F50b = "CCTGGCGCTAGTTGATAGATTAATAAGACAGCCCTTTAGAGCCCTGTCTC";
    private const string F80 = "CACGATGATCGAAAATCTGCGCTTTCTAGTGACCTGTATCTCATTGACCCACGGTTCGCCTATGCGTAGTCTTGACACCC";
    private const string F500 = "TTCTCGCAGTCTTCAAAGGGCCTTACTTTTAGTGACTTAATTCCATCAAATTGAGACCAAGACCGGCTCGTTGTGGGCTGTTGATATAGAAAAACCACAGTATTGCAACGGGTGACGTAATGGGAATGTCAATGACGGAGATTACAGCAAGGCAACAAGAGCGCATCGCCTAACATTTCTTTTACCACTCTTTGGACCCTGTCTTATAGGGGCTGAATAGTCCGGGGCGAGTAGGTGCTCCGTTCCGGATCAACTCAAGTCTCTCGCGAGCAGTTCTGGAAACCAATAATTCACCAAGAACCGCGACACAACATCTGCGTTAGCGTATTTGATGTGGCACGTAGAAGTTCAGTTCGAGGGTTGACGTCGGTAAGAACTTTCCTGAAGACTCGCGTATTCCTATGCGAGAGATGTCCCCGGTTCTTCGTCAGGATTAACTTTTAACGCTTCTGGGAGGGCGAGCAGCACGAACTCTCGGCGGGTCAAGGTTTATCTTTATA";
    private const string F40a = "TCCGGAAAGGGACACCCAGCAGACAGCTAGCACCGCATGC";
    private const string F40b = "CCTAGTGCCACTTTCGTAAGTTTGATTAAGGACGCGCGAA";

    private static string Rc(string s) => Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(s);

    private static string Complement(string s, params int[] positions)
    {
        var chars = s.ToCharArray();
        foreach (int i in positions)
            chars[i] = Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(chars[i]);
        return new string(chars);
    }

    private static (string, string)[] Asm(params string[] seqs) =>
        seqs.Select((s, i) => ($"c{i + 1}", s)).ToArray();

    private static string Label(CompletenessResult r) =>
        r.CompleteSingleCopy == 1 ? "Single" : r.CompleteDuplicated == 1 ? "Duplicated" : r.Fragmented == 1 ? "Fragmented" : "Missing";

    private static string ClassifyOne(string marker, (string, string)[] asm,
        double identity = 0.9, double coverage = 0.9, double minLength = DefaultMarkerMinLengthFraction)
        => Label(AssessCompleteness(asm, new[] { ("g", marker) }, identity, coverage, minLength));

    #region FindMarkerHits / classification vs oracle (M1–M4)

    public static IEnumerable<TestCaseData> OracleScenarios()
    {
        // name, marker, assembly, expected label, expected hits (id, strand, start, end, identity, coverage)
        yield return new TestCaseData("exact_forward", G, Asm(L + G + R), "Single",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0) });
        yield return new TestCaseData("reverse_strand", G, Asm(L + Rc(G) + R), "Single",
            new[] { ("c1", '-', 100, 399, 1.0, 1.0) });
        yield return new TestCaseData("two_contigs", G, Asm(L + G + R, F50a + G + F50b), "Duplicated",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0), ("c2", '+', 50, 349, 1.0, 1.0) });
        // compleasm get_region_clusters merges abutting half-open regions [100,400) and [400,700).
        yield return new TestCaseData("tandem_adjacent", G, Asm(L + G + G + R), "Single",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0), ("c1", '+', 400, 699, 1.0, 1.0) });
        yield return new TestCaseData("tandem_spaced", G, Asm(L + G + F80 + G + R), "Duplicated",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0), ("c1", '+', 480, 779, 1.0, 1.0) });
        yield return new TestCaseData("plus_and_minus_same_contig", G, Asm(L + G + F80 + Rc(G) + R), "Duplicated",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0), ("c1", '-', 480, 779, 1.0, 1.0) });
        yield return new TestCaseData("fragment_70pct", G, Asm(L + G[..210] + R), "Fragmented",
            new[] { ("c1", '+', 100, 309, 1.0, 0.7) });
        yield return new TestCaseData("fragment_40pct", G, Asm(L + G[..120] + R), "Fragmented",
            new[] { ("c1", '+', 100, 223, 123.0 / 124, 124.0 / 300) });
        yield return new TestCaseData("absent", G, Asm(F500), "Missing", Array.Empty<(string, char, int, int, double, double)>());
        yield return new TestCaseData("one_bp_deletion", G, Asm(L + G[..150] + G[151..] + R), "Single",
            new[] { ("c1", '+', 100, 398, 299.0 / 300, 1.0) });
        yield return new TestCaseData("complete_plus_fragment", G, Asm(L + G + R, F40a + G[..120] + F40b), "Single",
            new[] { ("c1", '+', 100, 399, 1.0, 1.0), ("c2", '+', 40, 159, 1.0, 0.4) });
        yield return new TestCaseData("three_substitutions", G, Asm(L + Complement(G, 50, 150, 250) + R), "Single",
            new[] { ("c1", '+', 100, 399, 0.99, 1.0) });
        yield return new TestCaseData("lowercase_with_N", G.ToLowerInvariant(),
            Asm((L + G[..100] + "N" + G[101..] + R).ToLowerInvariant()), "Single",
            new[] { ("c1", '+', 100, 399, 299.0 / 300, 1.0) });
    }

    [TestCaseSource(nameof(OracleScenarios))]
    public void FindMarkerHits_AndClassification_MatchOracle(string name, string marker, (string, string)[] assembly,
        string expectedLabel, (string Id, char Strand, int Start, int End, double Identity, double Coverage)[] expectedHits)
    {
        var hits = FindMarkerHits(marker, assembly);
        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(expectedHits.Length), name);
            for (int i = 0; i < Math.Min(hits.Count, expectedHits.Length); i++)
            {
                var (id, strand, start, end, identity, coverage) = expectedHits[i];
                Assert.That(hits[i].SequenceId, Is.EqualTo(id), $"{name} hit {i}");
                Assert.That(hits[i].Strand, Is.EqualTo(strand), $"{name} hit {i}");
                Assert.That(hits[i].Start, Is.EqualTo(start), $"{name} hit {i}");
                Assert.That(hits[i].End, Is.EqualTo(end), $"{name} hit {i}");
                Assert.That(hits[i].Identity, Is.EqualTo(identity).Within(Tol), $"{name} hit {i}");
                Assert.That(hits[i].Coverage, Is.EqualTo(coverage).Within(Tol), $"{name} hit {i}");
            }
            Assert.That(ClassifyOne(marker, assembly), Is.EqualTo(expectedLabel), name);
        });
    }

    [Test]
    public void FindMarkerHits_ReverseStrand_ReportsForwardCoordinates()
    {
        // RC(G[100..199]) at c1[100..199]; R starts with 3 bases that extend the alignment (oracle: c1 100–202,
        // identity 102/103, coverage 103/300). On '-' the contig end 202 pairs with marker position 97, so the
        // marker interval is reported in forward marker coordinates: 97–199.
        var hit = FindMarkerHits(G, Asm(L + Rc(G[100..200]) + R)).Single();
        Assert.That((hit.Strand, hit.Start, hit.End, hit.MarkerStart, hit.MarkerEnd), Is.EqualTo(('-', 100, 202, 97, 199)));
        Assert.That(hit.Identity, Is.EqualTo(102.0 / 103).Within(Tol));
        Assert.That(hit.Coverage, Is.EqualTo(103.0 / 300).Within(Tol));
    }

    [Test]
    public void AssessCompleteness_SeedRequirement_DivergedMarkerWithoutExactWordIsMissing()
    {
        // A substitution every 20 bp (identity 0.95) leaves no exact 28-bp word: the megablast seed
        // (word_size 28) cannot anchor it — Missing (oracle with the same seed rule agrees).
        var positions = Enumerable.Range(0, 15).Select(i => 10 + 20 * i).ToArray();
        Assert.That(ClassifyOne(G, Asm(L + Complement(G, positions) + R)), Is.EqualTo("Missing"));
    }

    [Test]
    public void AssessCompleteness_ShortMarker_IsSeededByItsWholeLength()
    {
        // m = 20 < 28: the seed is the whole marker (exact occurrence).
        var hit = FindMarkerHits(G[..20], Asm(L + G[..20] + R)).Single();
        Assert.That((hit.Start, hit.End, hit.Identity, hit.Coverage), Is.EqualTo((100, 119, 1.0, 1.0)));
        Assert.That(ClassifyOne(G[..20], Asm(L + G[..20] + R)), Is.EqualTo("Single"));
    }

    #endregion

    #region Thresholds and BUSCO percentages (M5–M7)

    [Test]
    public void AssessCompleteness_Thresholds_FollowCompleasmPassFilter()
    {
        // fragment_40pct hit: identity 123/124 ≈ 0.992, coverage 124/300 ≈ 0.413.
        var asm = Asm(L + G[..120] + R);
        Assert.Multiple(() =>
        {
            Assert.That(ClassifyOne(G, asm), Is.EqualTo("Fragmented"));                       // identity ≥ 0.9
            Assert.That(ClassifyOne(G, asm, identity: 1.0), Is.EqualTo("Missing"));          // neither clause
            Assert.That(ClassifyOne(G, asm, identity: 1.0, minLength: 0.4), Is.EqualTo("Fragmented")); // coverage ≥ 0.4
            // coverageThreshold is inclusive (Protein_mapped_rate >= min_complete): 0.7 ≥ 0.7 ⇒ Single.
            Assert.That(ClassifyOne(G, Asm(L + G[..210] + R), coverage: 0.7), Is.EqualTo("Single"));
        });
    }

    [Test]
    public void AssessCompleteness_MultiMarker_BuscoPercentagesOverN()
    {
        // g1 single, g2 on c1 and c2 (duplicated), g3 = reversed G (absent), g4 first 120 of 250 bp on c3.
        string g2 = F500[..200], g3 = new string(G.Reverse().ToArray()), g4 = F500[250..];
        var asm = Asm(L + G + F80 + g2 + R, F40a + g2 + F40b, F40b + g4[..120] + F40a);
        var r = AssessCompleteness(asm, new[] { ("g1", G), ("g2", g2), ("g3", g3), ("g4", g4) });

        Assert.Multiple(() =>
        {
            Assert.That((r.TotalGenes, r.Complete, r.CompleteSingleCopy, r.CompleteDuplicated, r.Fragmented, r.Missing),
                Is.EqualTo((4, 2, 1, 1, 1, 1)));
            Assert.That(r.CompletenessPercent, Is.EqualTo(50.0).Within(Tol));
            Assert.That(r.DuplicationPercent, Is.EqualTo(25.0).Within(Tol)); // BUSCO D/n (old D/C gave 50)
            Assert.That(r.SingleCopyPercent, Is.EqualTo(25.0).Within(Tol));
            Assert.That(r.FragmentedPercent, Is.EqualTo(25.0).Within(Tol));
            Assert.That(r.MissingPercent, Is.EqualTo(25.0).Within(Tol));
            Assert.That(r.BuscoNotation, Is.EqualTo("C:50.0%[S:25.0%,D:25.0%],F:25.0%,M:25.0%,n:4"));
        });
    }

    [TestCase(4, 2, 1, 0, 1, "C:75.0%[S:50.0%,D:25.0%],F:0.0%,M:25.0%,n:4")]
    [TestCase(3, 1, 0, 2, 0, "C:33.3%[S:33.3%,D:0.0%],F:66.7%,M:0.0%,n:3")]
    [TestCase(2000, 199, 0, 0, 1801, "C:10.0%[S:10.0%,D:0.0%],F:0.0%,M:90.0%,n:2000")]
    [TestCase(2000, 0, 1, 0, 1999, "C:0.1%[S:0.0%,D:0.1%],F:0.0%,M:100.0%,n:2000")] // Python round(0.05, 1) = 0.1
    public void BuscoNotation_MatchesBuscoOneLineSummary(int n, int s, int d, int f, int m, string expected)
    {
        var r = new CompletenessResult(n, s + d, s, d, f, m, (s + d) * 100.0 / n, d * 100.0 / n);
        Assert.That(r.BuscoNotation, Is.EqualTo(expected));
    }

    #endregion

    #region Validation (M8)

    [Test]
    public void AssessCompleteness_InvalidArguments_Throw()
    {
        var asm = Asm(L + G + R);
        var markers = new[] { ("g", G) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => AssessCompleteness(null!, markers));
            Assert.Throws<ArgumentNullException>(() => AssessCompleteness(asm, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => AssessCompleteness(asm, markers, identityThreshold: 1.1));
            Assert.Throws<ArgumentOutOfRangeException>(() => AssessCompleteness(asm, markers, coverageThreshold: -0.1));
            Assert.Throws<ArgumentOutOfRangeException>(() => AssessCompleteness(asm, markers, minLengthFraction: double.NaN));
            Assert.Throws<ArgumentException>(() => AssessCompleteness(asm, new[] { ("g", "") }));
            Assert.Throws<ArgumentException>(() => AssessCompleteness(asm, new[] { ("g", (string)null!) }));
            Assert.Throws<ArgumentException>(() => AssessCompleteness(new[] { ("c", (string)null!) }, markers));
            Assert.Throws<ArgumentNullException>(() => FindMarkerHits(null!, asm));
            Assert.Throws<ArgumentException>(() => FindMarkerHits("", asm));
        });
        var empty = AssessCompleteness(asm, Array.Empty<(string, string)>());
        Assert.That((empty.TotalGenes, empty.BuscoNotation), Is.EqualTo((0, "C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:0.0%,n:0")));
    }

    #endregion

    #region EstimateCompletenessFromKmers (M9–M14)

    // 20 canonical 3-mers: histogram 1:5, 2:2, 3:1, 4:3, 5:6, 6:3 (errors, valley at 3, peak at 5).
    private static readonly (string, int)[] ValleySpectrum =
    {
        ("AAA", 1), ("AAC", 1), ("AAG", 1), ("AAT", 1), ("ACA", 1), ("ACC", 2), ("ACG", 2), ("ACT", 3),
        ("AGA", 4), ("AGC", 4), ("AGG", 4), ("ATA", 5), ("ATC", 5), ("ATG", 5), ("CAA", 5), ("CAC", 5),
        ("CAG", 5), ("CCA", 6), ("CCC", 6), ("CCG", 6),
    };

    private static readonly (string, string)[] ValleyAssembly = { ("a", "ACGTTGCAAC"), ("b", "GGNNCATTAC") };

    [Test]
    public void EstimateCompletenessFromKmers_MerqurySolidThreshold_GenomeScopeErrorAndSize_MerquryCompletenessAndQv()
    {
        // Oracle: filt = 3 (slope −3, −1, +2 ⇒ depth 4 − 1); T = 72, E = 1·5+2·2+3·1 = 12;
        // error = 1 − (1 − 12/72)^(1/3) = 0.058963971118971 (GenomeScope uniroot agrees to 1e-16);
        // size = round(60/5) = 12; solid = 12; canonical assembly 3-mers ACG,AAC,CAA,GCA ×2, ATG,AAT,TAA,GTA ×1:
        // solid ∩ asm = {ATG, CAA} ⇒ completeness 2/12; ASM_ONLY = GCA(2)+GTA+TAA = 4 of TOTAL 12 ⇒
        // E = 1 − (1 − 4/12)^(1/3) = 0.126419535263701, QV = 8.981858105339732.
        var e = EstimateCompletenessFromKmers(ValleySpectrum, ValleyAssembly);
        Assert.Multiple(() =>
        {
            Assert.That((e.K, e.SolidThreshold, e.PeakCoverage, e.TotalKmers, e.ErrorKmers), Is.EqualTo((3, 3, 5, 72L, 12L)));
            Assert.That(e.ErrorRate, Is.EqualTo(0.05896397111897156).Within(Tol));
            Assert.That(e.EstimatedGenomeSize, Is.EqualTo(12L));
            Assert.That((e.SolidKmers, e.AssemblySolidKmers), Is.EqualTo((12L, 2L)));
            Assert.That(e.Completeness, Is.EqualTo(2.0 / 12).Within(Tol));
            Assert.That((e.AssemblyKmers, e.AssemblyOnlyKmers), Is.EqualTo((12L, 4L)));
            Assert.That(e.AssemblyErrorRate, Is.EqualTo(0.1264195352637011).Within(Tol));
            Assert.That(e.QualityValue, Is.EqualTo(8.981858105339732).Within(1e-9));
        });
    }

    [Test]
    public void EstimateCompletenessFromKmers_NonCanonical_CountsAssemblyForwardOnly()
    {
        // canonical = false: forward assembly 3-mers; only CAA is solid ⇒ 1/12; ASM_ONLY 9 ⇒ QV 4.317519438166913.
        var e = EstimateCompletenessFromKmers(ValleySpectrum, ValleyAssembly, canonical: false);
        Assert.Multiple(() =>
        {
            Assert.That(e.Completeness, Is.EqualTo(1.0 / 12).Within(Tol));
            Assert.That((e.AssemblyKmers, e.AssemblyOnlyKmers), Is.EqualTo((12L, 9L)));
            Assert.That(e.AssemblyErrorRate, Is.EqualTo(0.3700394750525634).Within(Tol));
            Assert.That(e.QualityValue, Is.EqualTo(4.317519438166913).Within(1e-9));
        });
    }

    [Test]
    public void EstimateCompletenessFromKmers_ExpectedCoverage_OverridesPeak()
    {
        // (72 − 12) / 4 = 15.
        var (completeness, errorRate, size) = EstimateCompletenessFromKmers(ValleySpectrum, expectedCoverage: 4);
        Assert.That(size, Is.EqualTo(15L));
        Assert.That(errorRate, Is.EqualTo(0.05896397111897156).Within(Tol));
        Assert.That(completeness, Is.NaN); // no assembly ⇒ no k-mer completeness
    }

    [Test]
    public void EstimateCompletenessFromKmers_MerquryQuirks_PlateauAndRisingStart()
    {
        // Plateau 1:3, 2:1, 3:1, 4:3 — slopes −2, 0, +2: the product with a zero slope is not < 0 ⇒ no
        // minimum, filt = 0 (kmerHistToPloidyDepth); peak = smallest count among the tied modes = 1 ⇒ size 20.
        var plateau = EstimateCompletenessFromKmers(new[]
        {
            ("AAA", 1), ("AAC", 1), ("AAG", 1), ("AAT", 2), ("ACA", 3), ("ACC", 4), ("ACG", 4), ("ACT", 4),
        }, null);
        Assert.That((plateau.SolidThreshold, plateau.PeakCoverage, plateau.EstimatedGenomeSize, plateau.ErrorRate),
            Is.EqualTo((0, 1, 20L, 0.0)));

        // Rising start 2:1, 3:2 — first slope vs Double.MAX_VALUE is negative, then +1 ⇒ filt = 3 − 1 = 2.
        var rising = EstimateCompletenessFromKmers(new[] { ("AAA", 2), ("AAC", 3), ("AAG", 3) }, null);
        Assert.Multiple(() =>
        {
            Assert.That((rising.SolidThreshold, rising.ErrorKmers, rising.PeakCoverage, rising.EstimatedGenomeSize, rising.SolidKmers),
                Is.EqualTo((2, 2L, 3, 2L, 2L)));
            Assert.That(rising.ErrorRate, Is.EqualTo(0.09143970358393017).Within(Tol));
        });
    }

    [Test]
    public void EstimateCompletenessFromKmers_NoAssemblyOnlyKmers_QvIsInfinite()
    {
        // "AAACC" → AAA, AAC, ACC: all in the read spectrum (none solid) ⇒ ASM_ONLY 0 ⇒ E = 0, QV = +∞,
        // completeness 0/12.
        var e = EstimateCompletenessFromKmers(ValleySpectrum, new[] { ("a", "AAACC") });
        Assert.Multiple(() =>
        {
            Assert.That((e.AssemblyKmers, e.AssemblyOnlyKmers, e.AssemblySolidKmers), Is.EqualTo((3L, 0L, 0L)));
            Assert.That(e.AssemblyErrorRate, Is.EqualTo(0.0));
            Assert.That(double.IsPositiveInfinity(e.QualityValue), Is.True);
            Assert.That(e.Completeness, Is.EqualTo(0.0));
        });

        // An assembly shorter than k has no k-mers: E and QV undefined (Merqury "asm has no kmers").
        var none = EstimateCompletenessFromKmers(ValleySpectrum, new[] { ("a", "AC") });
        Assert.That((none.AssemblyKmers, double.IsNaN(none.QualityValue), double.IsNaN(none.AssemblyErrorRate)),
            Is.EqualTo((0L, true, true)));
    }

    [Test]
    public void EstimateCompletenessFromKmers_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => EstimateCompletenessFromKmers(null!));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(Array.Empty<(string, int)>()));
            Assert.Throws<ArgumentOutOfRangeException>(() => EstimateCompletenessFromKmers(ValleySpectrum, expectedCoverage: -1));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ("AAA", 2), ("AA", 2) }));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ("AAA", 2), ("aaa", 3) }));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ("AAA", 0) }));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ((string)null!, 1) }));
            // Canonical evaluation needs canonical ACGT keys: TTT is the reverse complement of AAA.
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ("TTT", 2) }, ValleyAssembly));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(new[] { ("NAA", 2) }, ValleyAssembly));
            Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(ValleySpectrum, new[] { ("a", (string)null!) }));
        });
        // Non-canonical keys are accepted when canonical = false or without assembly.
        Assert.DoesNotThrow(() => EstimateCompletenessFromKmers(new[] { ("TTT", 2) }, ValleyAssembly, canonical: false));
        Assert.DoesNotThrow(() => EstimateCompletenessFromKmers(new[] { ("TTT", 2) }));
    }

    #endregion
}
