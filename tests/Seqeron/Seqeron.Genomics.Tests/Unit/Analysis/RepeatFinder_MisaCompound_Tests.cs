namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// MISA per-unit-size thresholds, MISA compound microsatellites, canonical motif classes (MISA / Krait) and
/// progress reporting of <see cref="RepeatFinder.FindMicrosatellites(DnaSequence,int,int,int,CancellationToken,IProgress{double})"/>
/// (REP-STR-001 / REP-TANDEM-001).
///
/// Sources (opened): MISA <c>misa.pl</c> v1.0 (Thiel et al. 2003, TAG 106:411; raw GitHub mirror
/// cfljam/SSR_marker_design) — <c>misa.ini</c> default <c>1-10 2-6 3-5 4-5 5-5 6-5</c>, interruptions 100, compound
/// assembly and the <c>.statistics</c> table "Frequency of classified repeat types (considering sequence
/// complementary)"; Krait (Du et al. 2018, lmdu/krait <c>src/motif.py</c> <c>StandardMotif</c>). Every expected
/// value below is copied from a real <c>perl misa.pl</c> run (default misa.ini) or from Krait's <c>standard()</c>
/// executed in Python — not from this implementation.
/// </summary>
[TestFixture]
public class RepeatFinder_MisaCompound_Tests
{
    /// <summary>101-bp spacer without any perfect repeat of ≥ 3 copies (units 1–6).</summary>
    private const string Spacer =
        "CTAAGCCAACTGCATTGCTAGAGCGAAGTCTTCGTAATGGACCGACCGTTCTGTCCGGACTAGTGAATCGCTGTACAAGTCCGAGGCATCAAGGACTAGTA";

    private const string StatSequence =
        "ACACACACACACT" + Spacer + "CACACACACACAT" + Spacer + "GTGTGTGTGTGTA" + Spacer + "TGTGTGTGTGTGA" + Spacer +
        "ACATACATACATACATACATG" + Spacer + "AAAAAAAAAAAAG" + Spacer + "TTTTTTTTTTTTTG";

    #region MISA per-unit-size thresholds

    [Test]
    public void MisaDefaultMinRepeats_IsMisaIniDefinition()
    {
        // misa.ini: "definition(unit_size,min_repeats): 1-10 2-6 3-5 4-5 5-5 6-5"
        Assert.That(RepeatFinder.MisaDefaultMinRepeats, Is.EquivalentTo(new Dictionary<int, int>
        {
            [1] = 10, [2] = 6, [3] = 5, [4] = 5, [5] = 5, [6] = 5,
        }));
        Assert.That(RepeatFinder.MisaDefaultMaxInterruption, Is.EqualTo(100));
    }

    [Test]
    public void FindMicrosatellites_MisaThresholds_MononucleotideNeedsTenCopies()
    {
        // misa.pl (default ini) on CCCCCCCCCGCCCCCCCCCC: only "p1 (C)10 10 11 20" — C×9 is below 1-10.
        const string seq = "CCCCCCCCCGCCCCCCCCCC";

        var misa = RepeatFinder.FindMicrosatellites(seq, RepeatFinder.MisaDefaultMinRepeats).ToList();
        var uniform = RepeatFinder.FindMicrosatellites(seq, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(misa, Is.EqualTo(new[]
            {
                new MicrosatelliteResult(10, "C", 10, 10, RepeatType.Mononucleotide),
            }));
            Assert.That(uniform.Select(m => (m.Position, m.RepeatUnit, m.RepeatCount)),
                Is.EqualTo(new[] { (0, "C", 9), (10, "C", 10) }));
        });
    }

    [Test]
    public void FindMicrosatellites_MisaThresholds_StatSequence_MatchesMisaSsrList()
    {
        // misa.pl raw SSR list: (AC)6 1-12, (AC)6 114-125, (GT)6 229-240, (TG)6 343-354, (ACAT)5 457-476,
        // (A)13 578-590, (T)13 693-705 (1-based inclusive).
        var ssrs = RepeatFinder.FindMicrosatellites(new DnaSequence(StatSequence), RepeatFinder.MisaDefaultMinRepeats)
            .OrderBy(m => m.Position)
            .Select(m => (m.RepeatUnit, m.RepeatCount, Start: m.Position + 1, End: m.Position + m.TotalLength))
            .ToList();

        Assert.That(ssrs, Is.EqualTo(new[]
        {
            ("AC", 6, 1, 12), ("AC", 6, 114, 125), ("GT", 6, 229, 240), ("TG", 6, 343, 354),
            ("ACAT", 5, 457, 476), ("A", 13, 578, 590), ("T", 13, 693, 705),
        }));
    }

    [Test]
    public void FindMicrosatellites_UniformMap_EqualsMinRepeatsOverload()
    {
        var rng = new Random(20260930);
        for (int t = 0; t < 300; t++)
        {
            var chars = new char[rng.Next(0, 300)];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = "ACGT"[rng.Next(rng.Next(1, 5))];
            string seq = new(chars);
            int k = rng.Next(2, 6);
            var map = Enumerable.Range(1, 6).ToDictionary(p => p, _ => k);

            Assert.That(RepeatFinder.FindMicrosatellites(seq, map).ToList(),
                Is.EqualTo(RepeatFinder.FindMicrosatellites(seq, 1, 6, k).ToList()), seq);
        }
    }

    [Test]
    public void FindMicrosatellites_Map_OnlyListedUnitLengthsSearched_OrderedByUnitLength()
    {
        // Unit lengths 3 and 1 only; the dinucleotide run is not searched.
        const string seq = "AAAAAGTCACACACACAGCAGCAGCAG";
        var map = new Dictionary<int, int> { [3] = 4, [1] = 5 };

        var result = RepeatFinder.FindMicrosatellites(seq, map).ToList();

        Assert.That(result.Select(m => (m.Position, m.RepeatUnit, m.RepeatCount)),
            Is.EqualTo(new[] { (0, "A", 5), (15, "CAG", 4) }));
    }

    [Test]
    public void FindMicrosatellites_Map_InvalidArguments_ThrowEagerly()
    {
        var dna = new DnaSequence("ACGT");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindMicrosatellites(dna, (IReadOnlyDictionary<int, int>)null!));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindMicrosatellites((DnaSequence)null!, RepeatFinder.MisaDefaultMinRepeats));
            Assert.Throws<ArgumentException>(() => RepeatFinder.FindMicrosatellites(dna, new Dictionary<int, int>()));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindMicrosatellites(dna, new Dictionary<int, int> { [0] = 3 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindMicrosatellites("ACGT", new Dictionary<int, int> { [2] = 1 }));
            // Unit sizes above 6 are valid (misa.pl accepts any size in its def line; B04 F61); size 0 is not.
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.GetTandemRepeatSummary(dna, new Dictionary<int, int> { [0] = 3 }));
            Assert.That(RepeatFinder.FindMicrosatellites((string)null!, RepeatFinder.MisaDefaultMinRepeats), Is.Empty);
        });
    }

    [Test]
    public void GetTandemRepeatSummary_MisaThresholds_StatSequence()
    {
        var s = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(StatSequence), RepeatFinder.MisaDefaultMinRepeats);

        // misa.pl .statistics: 7 SSRs; unit size 1: 2, 2: 4, 4: 1.
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(7));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(2));
            Assert.That(s.DinucleotideRepeats, Is.EqualTo(4));
            Assert.That(s.TrinucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.TetranucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(12 * 4 + 20 + 13 * 2));
            Assert.That(s.LongestRepeat?.RepeatUnit, Is.EqualTo("ACAT"));
            Assert.That(s.MostFrequentUnit, Is.EqualTo("AC"));
        });
    }

    #endregion

    #region Compound microsatellites (MISA c / c*)

    /// <summary>Rows copied from <c>perl misa.pl</c> (default misa.ini): type, SSR notation, size, start, end (1-based).</summary>
    private static IEnumerable<TestCaseData> MisaCompoundRows()
    {
        yield return new TestCaseData("ACGTATATATATATATccgtGAGAGAGAGAGAGAtttttAAAAAAAAAAAAT",
            "c", "(TA)6tccgt(GA)7ttttt(A)12", 48, 4, 51).SetName("MISA c: two interruptions");
        yield return new TestCaseData("nnnACACACACACACACAGCAGCAGCAGCAGCAG",
            "c*", "(AC)7(CAG)6*", 31, 4, 34).SetName("MISA c*: overlapping di/tri");
        yield return new TestCaseData("ATATATATATAT" + Spacer[..100] + "GAGAGAGAGAGA",
            "c", "(AT)6" + Spacer[..100].ToLowerInvariant() + "(GA)6", 124, 1, 124).SetName("MISA c: interruption of exactly 100 bases");
        yield return new TestCaseData("ATATATATATATGAGAGAGAGAGA",
            "c", "(AT)6(GA)6", 24, 1, 24).SetName("MISA c: adjacent SSRs, empty interruption");
        yield return new TestCaseData(
            "ACCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCGGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAAAAAAAAAAAAAAAGAGGGAGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
            "c*", "(CCCAG)7cccg(GCCCA)6(A)15*gagggaggga(G)36", 129, 2, 130).SetName("MISA c*: four components, mixed joins");
        yield return new TestCaseData(
            "CATCATCATCATCATCATCATCATCATCATCATGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATAATAATAATAATAATAATAATAATAATAA",
            "c*", "(CAT)11(GTCG)10gaaaa(T)34(TAA)11*", 144, 1, 144).SetName("MISA c*: adjacent + interrupted + overlapping");
    }

    [TestCaseSource(nameof(MisaCompoundRows))]
    public void FindCompoundMicrosatellites_MisaDefaults_MatchMisaRow(
        string sequence, string type, string notation, int size, int start1, int end1)
    {
        var compounds = RepeatFinder.FindCompoundMicrosatellites(sequence);

        Assert.That(compounds, Has.Count.EqualTo(1));
        var c = compounds[0];
        Assert.Multiple(() =>
        {
            Assert.That(c.MisaType, Is.EqualTo(type));
            Assert.That(c.Notation, Is.EqualTo(notation));
            Assert.That(c.Length, Is.EqualTo(size));
            Assert.That(c.Start + 1, Is.EqualTo(start1));
            Assert.That(c.End, Is.EqualTo(end1));
            Assert.That(c.Interruptions, Has.Count.EqualTo(c.Components.Count - 1));
        });
    }

    [Test]
    public void FindCompoundMicrosatellites_InterruptionOf101_TwoSingleSsrs()
    {
        // misa.pl: "p2 (AT)6 12 1 12" and "p2 (GA)6 12 114 125" — 101 interrupting bases > 100.
        string seq = "ATATATATATATC" + Spacer[..100] + "GAGAGAGAGAGA";

        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindCompoundMicrosatellites(seq), Is.Empty);
            Assert.That(RepeatFinder.FindCompoundMicrosatellites(seq, maxInterruption: 101), Has.Count.EqualTo(1));
        });
    }

    [Test]
    public void FindCompoundMicrosatellites_MaxInterruptionZero_OnlyAdjacentOrOverlapping()
    {
        // misa.pl with "interruptions 0": c5 → "c (AT)6(GA)6 24 1 24"; g0 → two p2 rows; with "interruptions 1":
        // g0 → "c (AT)6c(GA)6 25 1 25". Adjacent (0 bases) and overlapping SSRs always form compounds.
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindCompoundMicrosatellites("ATATATATATATGAGAGAGAGAGA", maxInterruption: 0)
                .Single().Notation, Is.EqualTo("(AT)6(GA)6"));
            Assert.That(RepeatFinder.FindCompoundMicrosatellites("ATATATATATATCGAGAGAGAGAGA", maxInterruption: 0), Is.Empty);
            Assert.That(RepeatFinder.FindCompoundMicrosatellites("ATATATATATATCGAGAGAGAGAGA", maxInterruption: 1)
                .Single().Notation, Is.EqualTo("(AT)6c(GA)6"));
        });
    }

    [Test]
    public void FindCompoundMicrosatellites_UniformMinRepeats_UsesUnitLengthsOneToSix()
    {
        // misa.pl with "definition 1-3 2-3 3-3 4-3 5-3 6-3, interruptions 10": "c (AT)3cc(GA)3 14 1 14".
        var c = RepeatFinder.FindCompoundMicrosatellites(new DnaSequence("ATATATCCGAGAGA"), 3, 10).Single();
        Assert.That(c.Notation, Is.EqualTo("(AT)3cc(GA)3"));
    }

    [Test]
    public void AssembleCompoundMicrosatellites_ChainComparesWithPreviousComponentEnd_NotMaxEnd()
    {
        // misa.pl compares start(i+1) with end(i) of the PREVIOUS SSR and sets end to the last component's end:
        // X = [0,30), Y = [5,15) nested (c*), Z starts at 20: 20 − 15 = 5 bases ≤ 5 → joined; compound end = Z end.
        string seq = new('A', 60);
        var ssrs = new[]
        {
            new MicrosatelliteResult(0, "A", 30, 30, RepeatType.Mononucleotide),
            new MicrosatelliteResult(5, "AAAAA", 2, 10, RepeatType.Pentanucleotide),
            new MicrosatelliteResult(20, "AA", 5, 10, RepeatType.Dinucleotide),
        };

        var c = RepeatFinder.AssembleCompoundMicrosatellites(seq, ssrs, maxInterruption: 5).Single();

        Assert.Multiple(() =>
        {
            Assert.That(c.Notation, Is.EqualTo("(A)30(AAAAA)2*aaaaa(AA)5"));
            Assert.That(c.IsOverlapping, Is.True);
            Assert.That((c.Start, c.End), Is.EqualTo((0, 30)));
            Assert.That(c.Interruptions, Is.EqualTo(new[] { "", "aaaaa" }));
            Assert.That(RepeatFinder.AssembleCompoundMicrosatellites(seq, ssrs, maxInterruption: 4).Single().Components,
                Has.Count.EqualTo(2));
        });
    }

    [Test]
    public void AssembleCompoundMicrosatellites_EqualStarts_KeepInputOrder()
    {
        // MISA orders equal starts by Perl hash order; the stable sort keeps the caller's order.
        string seq = "ACACACACAC";
        var a = new MicrosatelliteResult(0, "AC", 2, 4, RepeatType.Dinucleotide);
        var b = new MicrosatelliteResult(0, "ACACA", 2, 10, RepeatType.Pentanucleotide);

        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.AssembleCompoundMicrosatellites(seq, new[] { a, b }).Single().Notation,
                Is.EqualTo("(AC)2(ACACA)2*"));
            Assert.That(RepeatFinder.AssembleCompoundMicrosatellites(seq, new[] { b, a }).Single().Notation,
                Is.EqualTo("(ACACA)2(AC)2*"));
        });
    }

    [Test]
    public void FindCompoundMicrosatellites_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindCompoundMicrosatellites("ACGT", maxInterruption: -1));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindCompoundMicrosatellites((DnaSequence)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindCompoundMicrosatellites("ACGT", 1));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.AssembleCompoundMicrosatellites(null!, []));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.AssembleCompoundMicrosatellites("ACGT",
                new[] { new MicrosatelliteResult(2, "A", 3, 3, RepeatType.Mononucleotide) }));
            Assert.That(RepeatFinder.FindCompoundMicrosatellites(string.Empty), Is.Empty);
        });
    }

    /// <summary>
    /// Documented convention differences (REP-STR-001 maximal primitive runs vs MISA's resumed regex scan), found in
    /// the 6 048-sequence misa.pl cross-check; expected values from misa.pl and from the brute-force maximal-run reference.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MisaConventionDifferences_AreTheDocumentedOnes()
    {
        // (a) same-size run overlapping the previous match by < p bases: MISA resumes after the (AAAAG)9 match and
        // reports "(AGAAA)8 46-85"; the true maximal run starts one base earlier: (AAGAA)8 at 45-84 (1-based).
        const string a = "AAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGACTTACGTAGATAAGACTTACGTAGATAAG";
        var pentaA = RepeatFinder.FindMicrosatellites(a, RepeatFinder.MisaDefaultMinRepeats)
            .Where(m => m.RepeatUnit.Length == 5).Select(m => (m.RepeatUnit, m.RepeatCount, m.Position + 1)).ToList();

        // (b) greedy consumption by a rejected non-primitive match: MISA's hexamer regex matches (CTCTCT)5 at
        // 101-130, rejects it as redundant, resumes at 131 and reports "(TAAACT)6 131-166"; the maximal run is
        // (CTTAAA)7 at 129-170.
        const string b = "GCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGGCACAAAAAAGAAAACTATCAGGAATAGAGTATAGAGTAGGGGGGGGGGGGAAACAACTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTAAACTTAAACTTAAACTTAAACTTAAACTTAAACTTAAAATTAAACT";
        var hexaB = RepeatFinder.FindMicrosatellites(b, RepeatFinder.MisaDefaultMinRepeats)
            .Where(m => m.RepeatUnit == "CTTAAA").Select(m => (m.RepeatUnit, m.RepeatCount, m.Position + 1)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(pentaA, Does.Contain(("AAGAA", 8, 45)));
            Assert.That(pentaA, Does.Not.Contain(("AGAAA", 8, 46)));
            Assert.That(hexaB, Is.EqualTo(new[] { ("CTTAAA", 7, 129) }));
        });
    }

    #endregion

    #region Canonical motifs (MISA classes, Krait standard motifs)

    [TestCase("AC", "AC/GT")]
    [TestCase("CA", "AC/GT")]
    [TestCase("GT", "AC/GT")]
    [TestCase("tg", "AC/GT")]
    [TestCase("A", "A/T")]
    [TestCase("T", "A/T")]
    [TestCase("AT", "AT/AT")]
    [TestCase("ACAT", "ACAT/ATGT")]
    [TestCase("CAG", "AGC/CTG")]
    [TestCase("GCCCA", "AGCCC/CTGGG")]
    public void GetCanonicalMotifClass_MatchesMisaStatisticsRowName(string motif, string expected) =>
        Assert.That(RepeatFinder.GetCanonicalMotifClass(motif), Is.EqualTo(expected));

    // Krait motif.py StandardMotif(level).standard(motif) with a fresh cache, run in Python.
    [TestCase("GT", new[] { "GT", "TG", "AC", "AC", "AC" })]
    [TestCase("CA", new[] { "CA", "AC", "AC", "AC", "AC" })]
    [TestCase("T", new[] { "T", "T", "A", "A", "A" })]
    [TestCase("ACAT", new[] { "ACAT", "ATAC", "ATAC", "ATAC", "ATAC" })]
    [TestCase("AGCCC", new[] { "AGCCC", "AGCCC", "AGCCC", "AGCCC", "ACCCG" })]
    [TestCase("CTGGG", new[] { "CTGGG", "TGGGC", "AGCCC", "ACCCG", "ACCCG" })]
    [TestCase("ACT", new[] { "ACT", "ACT", "ACT", "ATG", "ATC" })]
    [TestCase("CTG", new[] { "CTG", "TGC", "AGC", "ACG", "ACG" })]
    public void GetStandardMotif_MatchesKraitAllLevels(string motif, string[] expectedByLevel)
    {
        for (int level = 0; level <= 4; level++)
            Assert.That(RepeatFinder.GetStandardMotif(motif, level), Is.EqualTo(expectedByLevel[level]), $"level {level}");
    }

    [Test]
    public void GetCanonicalMotif_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RepeatFinder.GetCanonicalMotifClass(""));
            Assert.Throws<ArgumentException>(() => RepeatFinder.GetCanonicalMotifClass("ACN"));
            Assert.Throws<ArgumentException>(() => RepeatFinder.GetStandardMotif(null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.GetStandardMotif("AC", 5));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.GetStandardMotifFrequencies([], -1));
        });
    }

    [Test]
    public void GetCanonicalMotifFrequencies_StatSequence_MatchesMisaClassifiedTable()
    {
        // misa.pl .statistics "Frequency of classified repeat types (considering sequence complementary)", total column:
        // A/T 2, AC/GT 4, ACAT/ATGT 1.
        var ssrs = RepeatFinder.FindMicrosatellites(StatSequence, RepeatFinder.MisaDefaultMinRepeats).ToList();

        var classes = RepeatFinder.GetCanonicalMotifFrequencies(ssrs);
        var krait = RepeatFinder.GetStandardMotifFrequencies(ssrs);

        Assert.Multiple(() =>
        {
            Assert.That(classes.ToList(), Is.EqualTo(new[]
            {
                new KeyValuePair<string, int>("A/T", 2),
                new KeyValuePair<string, int>("AC/GT", 4),
                new KeyValuePair<string, int>("ACAT/ATGT", 1),
            }));
            Assert.That(krait, Is.EquivalentTo(new Dictionary<string, int> { ["A"] = 2, ["AC"] = 4, ["ATAC"] = 1 }));
            Assert.That(RepeatFinder.GetCanonicalMotifFrequencies([]), Is.Empty);
        });
    }

    #endregion

    #region Progress reporting (cancellable overloads)

    private sealed class RecordingProgress(Action<double>? onReport = null) : IProgress<double>
    {
        public List<double> Values { get; } = [];

        public void Report(double value)
        {
            Values.Add(value);
            onReport?.Invoke(value);
        }
    }

    private static string RandomDna(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < chars.Length; i++)
            chars[i] = "ACGT"[rng.Next(4)];
        return new string(chars);
    }

    [Test]
    public void FindMicrosatellites_Progress_MonotoneInUnitIntervalEndingAtOne()
    {
        var dna = new DnaSequence(RandomDna(50_000, 17));
        var progress = new RecordingProgress();

        var withProgress = RepeatFinder.FindMicrosatellites(dna, 1, 6, 3, CancellationToken.None, progress).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(withProgress, Is.EqualTo(RepeatFinder.FindMicrosatellites(dna, 1, 6, 3).ToList()));
            Assert.That(progress.Values, Has.Count.GreaterThan(10));
            Assert.That(progress.Values, Is.Ordered.Ascending);
            Assert.That(progress.Values, Has.All.InRange(0.0, 1.0));
            Assert.That(progress.Values[^1], Is.EqualTo(1.0));
            Assert.That(progress.Values.Take(progress.Values.Count - 1), Has.All.LessThan(1.0));
        });
    }

    [Test]
    public void FindMicrosatellites_Progress_ShortSequence_ReportsOnlyFinalOne()
    {
        var progress = new RecordingProgress();
        _ = RepeatFinder.FindMicrosatellites("CACACA", 1, 6, 3, CancellationToken.None, progress).ToList();
        Assert.That(progress.Values, Is.EqualTo(new[] { 1.0 }));
    }

    [Test]
    public void FindMicrosatellites_Progress_MapOverload_MonotoneEndingAtOne()
    {
        var progress = new RecordingProgress();
        _ = RepeatFinder.FindMicrosatellites(RandomDna(30_000, 5), RepeatFinder.MisaDefaultMinRepeats,
            CancellationToken.None, progress).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(progress.Values, Has.Count.GreaterThan(10));
            Assert.That(progress.Values, Is.Ordered.Ascending);
            Assert.That(progress.Values[^1], Is.EqualTo(1.0));
        });
    }

    [Test]
    public void FindMicrosatellites_PreCancelledToken_ThrowsOnEnumeration()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.Multiple(() =>
        {
            Assert.Throws<OperationCanceledException>(() =>
                RepeatFinder.FindMicrosatellites(new DnaSequence("CACACA"), 1, 6, 3, cts.Token).ToList());
            Assert.Throws<OperationCanceledException>(() =>
                RepeatFinder.FindMicrosatellites("CACACA", 1, 6, 3, cts.Token).ToList());
            Assert.Throws<OperationCanceledException>(() =>
                RepeatFinder.FindMicrosatellites("CACACA", RepeatFinder.MisaDefaultMinRepeats, cts.Token).ToList());
        });
    }

    [Test]
    public void FindMicrosatellites_CancelledDuringScan_ThrowsAndStopsReporting()
    {
        using var cts = new CancellationTokenSource();
        var progress = new RecordingProgress(_ => cts.Cancel());

        Assert.Throws<OperationCanceledException>(() =>
            RepeatFinder.FindMicrosatellites(RandomDna(100_000, 3), 1, 6, 3, cts.Token, progress).ToList());

        // Cancelled at the first report; the next check (≤ 1000 run starts later) throws before any further report.
        Assert.That(progress.Values, Has.Count.EqualTo(1));
        Assert.That(progress.Values[0], Is.LessThan(1.0));
    }

    #endregion
}
