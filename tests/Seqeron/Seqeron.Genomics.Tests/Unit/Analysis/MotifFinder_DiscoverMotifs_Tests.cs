// MOTIF-DISCOVER-001 — Motif Discovery via Overrepresented k-mers (observed/expected enrichment)
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md
// Source: Compeau P, Pevzner P (2015). Bioinformatics Algorithms: An Active Learning Approach, 2nd ed., Ch. 2.
//         Expected occurrences of a specific k-mer: E = (N - k + 1) / 4^k under the i.i.d. uniform background.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical test class for MOTIF-DISCOVER-001: overrepresented k-mer discovery in a single
/// DNA sequence. Verifies <see cref="MotifFinder.DiscoverMotifs(DnaSequence, int, int)"/>
/// against the Compeau &amp; Pevzner expected-count formula E = (N - k + 1) / 4^k and the
/// observed/expected (O/E) enrichment ratio.
/// </summary>
[TestFixture]
public class MotifFinder_DiscoverMotifs_Tests
{
    #region DiscoverMotifs — MUST

    // M1 — Repeated k-mer in a tandem string. "ATGCATGCATGC" (N=12) contains "ATGC" at
    // windows 0,4,8 -> Count = 3. Source: deterministic window enumeration.
    [Test]
    public void DiscoverMotifs_TandemRepeat_ReturnsKmerWithExactCount()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        var atgc = motifs.Single(m => m.Sequence == "ATGC");
        Assert.That(atgc.Count, Is.EqualTo(3),
            "\"ATGC\" starts at windows 0, 4 and 8 in ATGCATGCATGC, so its observed count is exactly 3.");
    }

    // M2 — Exact occurrence positions (0-based window starts). INV-01.
    [Test]
    public void DiscoverMotifs_TandemRepeat_ReturnsExactPositions()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        var atgc = motifs.Single(m => m.Sequence == "ATGC");
        Assert.That(atgc.Positions, Is.EqualTo(new[] { 0, 4, 8 }),
            "The 0-based start positions of \"ATGC\" in ATGCATGCATGC are exactly {0,4,8} (INV-01).");
    }

    // M3 — Exact O/E enrichment (tandem). N=12, k=4 -> windows = 9, E = 9/4^4 = 9/256.
    // "ATGC" count = 3 -> enrichment = 3 / (9/256) = 768/9. Source: Compeau & Pevzner.
    [Test]
    public void DiscoverMotifs_TandemRepeat_ComputesExactObservedOverExpectedEnrichment()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        const double expectedEnrichment = 768.0 / 9.0; // 3 / (9/256), E = (12-4+1)/4^4

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        var atgc = motifs.Single(m => m.Sequence == "ATGC");
        Assert.That(atgc.Enrichment, Is.EqualTo(expectedEnrichment).Within(1e-10),
            "Enrichment = Count / E with E = (N-k+1)/4^k = 9/256, so 3/(9/256) = 768/9 (INV-02).");
    }

    // M4 — Exact O/E enrichment (homopolymer). N=10, k=3 -> windows = 8, E = 8/4^3 = 0.125.
    // "AAA" count = 8 -> enrichment = 8/0.125 = 64.0. Source: Compeau & Pevzner.
    [Test]
    public void DiscoverMotifs_Homopolymer_ComputesExactCountAndEnrichment()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        const double expectedEnrichment = 64.0; // 8 / (8/64), E = (10-3+1)/4^3 = 0.125

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 3, minCount: 1).ToList();

        var aaa = motifs.Single(m => m.Sequence == "AAA");
        Assert.Multiple(() =>
        {
            Assert.That(aaa.Count, Is.EqualTo(8),
                "\"AAA\" occurs at windows 0..7 in a 10-base homopolymer, so its count is exactly 8.");
            Assert.That(aaa.Enrichment, Is.EqualTo(expectedEnrichment).Within(1e-10),
                "Enrichment = 8 / ((10-3+1)/4^3) = 8 / 0.125 = 64.0 (INV-02).");
        });
    }

    // M5 — minCount filter both INCLUDES count>=minCount k-mers and EXCLUDES count<minCount ones.
    // In "ACGTACGTAA" (N=10), k=4 the 4-mer multiset is:
    //   ACGT -> {0,4} (count 2), CGTA -> {1,5} (count 2),
    //   GTAC -> {2}, TACG -> {3}, GTAA -> {6} (each count 1).
    // With minCount=2 exactly {ACGT, CGTA} must be returned and the three singletons excluded.
    // (Verified by hand window enumeration.) This is a non-vacuous filter test: it would fail
    // both an implementation that returns nothing and one that ignores the minCount filter.
    [Test]
    public void DiscoverMotifs_MinCountTwo_KeepsAtOrAboveThresholdAndDropsBelow()
    {
        var sequence = new DnaSequence("ACGTACGTAA");

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        var returned = motifs.Select(m => m.Sequence).OrderBy(s => s).ToArray();
        Assert.Multiple(() =>
        {
            // INV-03: every returned motif satisfies Count >= minCount.
            Assert.That(motifs.All(m => m.Count >= 2), Is.True,
                "With minCount=2 no k-mer occurring only once may be returned (INV-03).");
            // Exactly the two count-2 k-mers are returned (inclusion side).
            Assert.That(returned, Is.EqualTo(new[] { "ACGT", "CGTA" }),
                "Only ACGT (count 2) and CGTA (count 2) reach minCount=2 in ACGTACGTAA.");
            // The three count-1 k-mers are excluded (exclusion side).
            Assert.That(returned, Does.Not.Contain("GTAC")
                                    .And.Not.Contain("TACG")
                                    .And.Not.Contain("GTAA"),
                "Singletons GTAC, TACG, GTAA occur once and must be filtered out.");
        });
    }

    #endregion

    #region DiscoverMotifs — SHOULD (edge cases)

    // S1 — Null sequence throws ArgumentNullException.
    [Test]
    public void DiscoverMotifs_NullSequence_ThrowsArgumentNullException()
    {
        Assert.That(() => MotifFinder.DiscoverMotifs(null!).ToList(),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null sequence is invalid input.");
    }

    // S2 — k < 1 throws ArgumentOutOfRangeException.
    [Test]
    public void DiscoverMotifs_ZeroK_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.That(() => MotifFinder.DiscoverMotifs(sequence, k: 0).ToList(),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
            "k must be at least 1; k=0 is out of range.");
    }

    // S3 — k > N: no length-k windows exist, so the result is empty.
    [Test]
    public void DiscoverMotifs_KGreaterThanSequenceLength_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAA");

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 5, minCount: 1).ToList();

        Assert.That(motifs, Is.Empty,
            "With k=5 > N=3 there are no length-5 windows, so no motifs can be discovered.");
    }

    #endregion

    #region DiscoverMotifs — COULD

    // C1 — No floor on the expected count. N=12, k=6 -> windows = 7, E = 7/4^6 = 7/4096.
    // "AACCGG" occurs at 0 and 6 -> count 2 -> enrichment = 2 / (7/4096) = 8192/7 ~= 1170.286.
    // A reintroduced max(E, 0.1) clamp would instead give 2/0.1 = 20, so this value pins the
    // exact O/E ratio and guards against the previous defect. Source: Compeau & Pevzner.
    [Test]
    public void DiscoverMotifs_SmallExpectedCount_UsesUnclampedExpectedDenominator()
    {
        var sequence = new DnaSequence("AACCGGAACCGG");
        const double expectedEnrichment = 8192.0 / 7.0; // 2 / ((12-6+1)/4^6) = 2 / (7/4096)

        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 6, minCount: 2).ToList();

        var aaccgg = motifs.Single(m => m.Sequence == "AACCGG");
        Assert.Multiple(() =>
        {
            Assert.That(aaccgg.Count, Is.EqualTo(2),
                "\"AACCGG\" starts at windows 0 and 6, so its count is exactly 2.");
            Assert.That(aaccgg.Enrichment, Is.EqualTo(expectedEnrichment).Within(1e-9),
                "Enrichment uses the unclamped E = (12-6+1)/4^6 = 7/4096, giving 2/(7/4096) = 8192/7.");
        });
    }

    #endregion

    #region DiscoverMotifs — review 2026-09 (RSAT oligo-analysis Bernoulli model, numerics, canonical counting)

    // Reference values: RSAT oligo-analysis (van Helden et al. 1998; perl-scripts/oligo-analysis,
    // raw.githubusercontent.com/rsa-tools/rsat-code): exp_freq = prod residue_proba (Bernoulli) or
    // 1/4^k (equiprobable), exp_occ = exp_freq * sum_occurrences (overlapping windows, N-k+1),
    // ratio = occ / exp_occ — recomputed with exact Python Fractions (scratch ref_discover.py).

    private static readonly double[] RsatBackground = { 0.3, 0.2, 0.2, 0.3 };

    // Deterministic 512-mer (LCG x = (1103515245x + 12345) mod 2^31, base = "ACGT"[(x >> 16) & 3], x0 = 1).
    private static string Lcg512Mer()
    {
        var sb = new System.Text.StringBuilder(512);
        long x = 1;
        for (int i = 0; i < 512; i++)
        {
            x = (x * 1103515245 + 12345) % 2147483648L;
            sb.Append("ACGT"[(int)((x >> 16) & 3)]);
        }
        return sb.ToString();
    }

    [Test]
    public void DiscoverMotifs_BernoulliBackground_EqualsRsatRatio()
    {
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), 4, 2, RsatBackground)
            .ToDictionary(m => m.Sequence);

        Assert.Multiple(() =>
        {
            Assert.That(motifs.Keys, Is.EquivalentTo(new[] { "ATGC", "TGCA", "GCAT", "CATG" }));
            // ATGC: p = .3*.3*.2*.2 = 0.0036, E = 9*0.0036 = 0.0324, 3/0.0324
            Assert.That(motifs["ATGC"].Enrichment, Is.EqualTo(92.5925925925926).Within(1e-12));
            Assert.That(motifs["ATGC"].Positions, Is.EqualTo(new[] { 0, 4, 8 }));
            // TGCA/GCAT/CATG: count 2, E = 0.0324 -> 2/0.0324
            Assert.That(motifs["TGCA"].Enrichment, Is.EqualTo(61.72839506172839).Within(1e-12));
            Assert.That(motifs["GCAT"].Enrichment, Is.EqualTo(61.72839506172839).Within(1e-12));
            Assert.That(motifs["CATG"].Enrichment, Is.EqualTo(61.72839506172839).Within(1e-12));
        });
    }

    [Test]
    public void DiscoverMotifs_BernoulliBackground_IsNormalisedAndHomopolymerMatchesRsat()
    {
        // Unnormalised (3,2,2,3) == (0.3,0.2,0.2,0.3); AAA in A10: 8 / (8 * 0.3^3) = 37.037037...
        var raw = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), 4, 2, new double[] { 3, 2, 2, 3 })
            .Single(m => m.Sequence == "ATGC");
        var aaa = MotifFinder.DiscoverMotifs(new DnaSequence("AAAAAAAAAA"), 3, 1, RsatBackground).Single();

        Assert.Multiple(() =>
        {
            Assert.That(raw.Enrichment, Is.EqualTo(92.5925925925926).Within(1e-12));
            Assert.That(aaa.Enrichment, Is.EqualTo(37.03703703703704).Within(1e-12));
        });
    }

    [Test]
    public void DiscoverMotifs_UniformBackgroundOverload_EqualsDefaultExactly()
    {
        var rng = new Random(20260929);
        for (int trial = 0; trial < 50; trial++)
        {
            var chars = new char[rng.Next(1, 200)];
            for (int i = 0; i < chars.Length; i++) chars[i] = "ACGT"[rng.Next(4)];
            var dna = new DnaSequence(new string(chars));
            int k = rng.Next(1, 7);

            var a = MotifFinder.DiscoverMotifs(dna, k, 1).ToList();
            var b = MotifFinder.DiscoverMotifs(dna, k, 1, new[] { 1.0, 1.0, 1.0, 1.0 }).ToList();

            Assert.That(b.Select(m => (m.Sequence, m.Count, m.Enrichment)),
                Is.EqualTo(a.Select(m => (m.Sequence, m.Count, m.Enrichment))));
        }
    }

    [Test]
    public void DiscoverMotifs_InvalidBackground_Throws()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(() => MotifFinder.DiscoverMotifs(dna, 2, 1, null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => MotifFinder.DiscoverMotifs(dna, 2, 1, new[] { 0.5, 0.5, 0.0 }),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.DiscoverMotifs(dna, 2, 1, new[] { 0.5, 0.5, 0.0, 0.0 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.DiscoverMotifs(dna, 2, 1, new[] { 0.25, double.NaN, 0.25, 0.25 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.DiscoverMotifs(dna, 0, 1, new[] { 0.25, 0.25, 0.25, 0.25 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // 4^512 = 2^1024 overflows a double, so the former E = (N-k+1)/Math.Pow(4,k) became 0 and the
    // enrichment +Infinity, although the true ratio 2*4^512/513 = 7.008550233381348e+305 is finite.
    [Test]
    public void DiscoverMotifs_LongK_FiniteRatioDoesNotOverflow()
    {
        string x = Lcg512Mer();
        var dna = new DnaSequence(x + x);

        var uniform = MotifFinder.DiscoverMotifs(dna, 512, 2).Single();
        var bernoulli = MotifFinder.DiscoverMotifs(dna, 512, 2, new double[] { 26, 24, 24, 26 }).Single();

        Assert.Multiple(() =>
        {
            Assert.That(uniform.Sequence, Is.EqualTo(x));
            Assert.That(uniform.Positions, Is.EqualTo(new[] { 0, 512 }));
            Assert.That(double.IsFinite(uniform.Enrichment), Is.True);
            Assert.That(uniform.Enrichment / 7.008550233381348e+305, Is.EqualTo(1.0).Within(1e-14));
            Assert.That(double.IsFinite(bernoulli.Enrichment), Is.True);
            Assert.That(bernoulli.Enrichment / 3.2383751592973165e+306, Is.EqualTo(1.0).Within(1e-12));
        });
    }

    // Counts equal the canonical k-mer counter; positions ascending; order = first occurrence.
    [Test]
    public void DiscoverMotifs_CountsEqualCanonicalCounter_PositionsAscending_FirstOccurrenceOrder()
    {
        var ordered = MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGTAA"), 4, 1).Select(m => m.Sequence);
        Assert.That(ordered, Is.EqualTo(new[] { "ACGT", "CGTA", "GTAC", "TACG", "GTAA" }));

        var rng = new Random(7);
        for (int trial = 0; trial < 30; trial++)
        {
            var chars = new char[rng.Next(1, 300)];
            for (int i = 0; i < chars.Length; i++) chars[i] = "ACGT"[rng.Next(4)];
            string s = new(chars);
            int k = rng.Next(1, 6);

            var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(s), k, 1).ToList();
            var canonical = KmerAnalyzer.CountKmers(s, k);

            Assert.That(motifs.ToDictionary(m => m.Sequence, m => m.Count), Is.EquivalentTo(canonical));
            foreach (var m in motifs)
            {
                Assert.That(m.Positions, Is.Ordered.Ascending);
                Assert.That(m.Positions.All(p => s.Substring(p, k) == m.Sequence), Is.True);
            }
        }
    }

    #endregion
}
