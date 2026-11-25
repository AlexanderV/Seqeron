// RNA-ACCESS-001 — region accessibility: P(a stretch of u nucleotides is entirely unpaired) = Z_open / Z
//   on the whole sequence (McCaskill 1990; Bernhart et al. 2006 RNAplfold with W, L >= n; Kertesz et al. 2007).
// Evidence: docs/Evidence/RNA-ACCESS-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-ACCESS-001.md
// Oracles: (1) ViennaRNA 2.7.2 exact ratio exp(-(G_constrained - G_ensemble)/RT) with every position of the region
//          hard-constrained unpaired (RNA.fold_compound(...).hc_add_up(p, CONSTRAINT_CONTEXT_ALL_LOOPS)), NNDB-precision
//          parameter file, md passed explicitly; for 37 C it equals RNA.pfl_fold_up(seq, 14, n, n, n) (RNAplfold, W = L = n)
//          to <= 2e-6 on the cases below (600 random queries: <= 1.0e-3 = Vienna's SMOOTH() quirk at the sequence ends).
//          (2) exhaustive enumeration of every pseudoknot-free structure scored with CalculateFreeEnergy.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_RegionAccessibility_Tests
{
    private static double Rt(double kelvin) => 1.98717 * kelvin / 1000.0;

    // sequences without NNDB-only terms ViennaRNA lacks (no repeated neighbours), generated with a fixed seed
    // (value = ViennaRNA 2.7.2 ratio, see header). Columns: sequence, Celsius, 0-based window end, window length u, P(unpaired).
    private static readonly (string Seq, double Celsius, int End, int Len, double P)[] ViennaCases =
    {
        ("GGGCGAAAGCCCAAGGGCGAAAGCCC", 37, 25, 1, 0.018548),
        ("GGGCGAAAGCCCAAGGGCGAAAGCCC", 60, 13, 4, 0.001941),
        ("GCGGAAACGCCGAAAGGCGC", 60, 19, 1, 0.986769),
        ("GCGGAAACGCCGAAAGGCGC", 60, 11, 4, 0.008907),
        ("GCGGAAACGCCGAAAGGCGC", 60, 14, 14, 0.005938),
        ("GCGGAAACGCCGAAAGGCGC", 60, 19, 14, 0.005934),
        ("AGUGUACUGUCUGUCUCACUACGAGAUCUCAGCGCGUGUA", 37, 18, 4, 0.004353),
        ("AGUGUACUGUCUGUCUCACUACGAGAUCUCAGCGCGUGUA", 60, 39, 1, 0.949283),
        ("AGUGUACUGUCUGUCUCACUACGAGAUCUCAGCGCGUGUA", 60, 18, 4, 0.059001),
        ("AGUGUACUGUCUGUCUCACUACGAGAUCUCAGCGCGUGUA", 60, 29, 14, 0.005321),
        ("AGUGUACUGUCUGUCUCACUACGAGAUCUCAGCGCGUGUA", 60, 39, 14, 0.018726),
        ("GAUGCGCAUGUGCAUCGCGAGUGCAGUGCGUCGUCUGAUACAUGUAUGCGAGACUGCUCU", 37, 30, 1, 0.263641),
        ("GAUGCGCAUGUGCAUCGCGAGUGCAGUGCGUCGUCUGAUACAUGUAUGCGAGACUGCUCU", 37, 33, 8, 0.003516),
        ("GAUGCGCAUGUGCAUCGCGAGUGCAGUGCGUCGUCUGAUACAUGUAUGCGAGACUGCUCU", 37, 59, 14, 0.008010),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 37, 17, 1, 0.143155),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 37, 34, 14, 0.030532),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 20, 17, 1, 0.100328),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 20, 34, 14, 0.001816),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 60, 17, 1, 0.294539),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 60, 20, 8, 0.014357),
        ("GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA", 60, 13, 14, 0.018147),
        ("GCGAGUGUCGAGUCAUCGCACUCAGACAGUGAUACAUAUCUGCACGUGUCUGAUAGUAUGUCAGUACUCUCGAUAGACGA", 37, 79, 1, 0.988150),
        ("GCGAGUGUCGAGUCAUCGCACUCAGACAGUGAUACAUAUCUGCACGUGUCUGAUAGUAUGUCAGUACUCUCGAUAGACGA", 37, 31, 4, 0.006842),
    };

    public static IEnumerable<TestCaseData> ViennaData => ViennaCases.Select(c =>
        new TestCaseData(c.Seq, c.Celsius, c.End, c.Len, c.P).SetName($"Vienna_{c.Seq.Length}nt_{c.Celsius:0}C_end{c.End}_u{c.Len}"));

    // ViennaRNA's pf() weights sequence-end dangles with SMOOTH() (params.c) and our NNDB parameters differ from the
    // Vienna file by rounding, so the agreement is 0.0009 absolute (max over 46 locked-style cases; 0.0056 over 1440 random queries).
    [TestCaseSource(nameof(ViennaData))]
    public void CalculateRegionUnpairedProbability_MatchesViennaRnaConstrainedEnsemble(string seq, double celsius, int end, int len, double expected)
    {
        double got = CalculateRegionUnpairedProbability(seq, end, len, temperature: celsius + 273.15);
        Assert.That(got, Is.EqualTo(expected).Within(Math.Max(0.0015, 0.02 * expected)),
            $"P({len} nt ending at {end} unpaired) of {seq} at {celsius} C");
    }

    // Independent exact oracle: exhaustive enumeration of every pseudoknot-free structure (hairpin loop >= minLoop),
    // scored by CalculateFreeEnergy; P(region unpaired) = Σ_{S: no pair touches the region} w(S) / Σ_S w(S).
    private static IEnumerable<List<(int, int)>> Structures(string seq, int i, int j, int minLoop)
    {
        if (i >= j) { yield return new List<(int, int)>(); yield break; }
        foreach (var s in Structures(seq, i + 1, j, minLoop)) yield return s;
        for (int k = i + minLoop + 1; k < j; k++)
        {
            if (!CanPair(seq[i], seq[k])) continue;
            foreach (var inner in Structures(seq, i + 1, k, minLoop))
                foreach (var rest in Structures(seq, k + 1, j, minLoop))
                {
                    var l = new List<(int, int)> { (i, k) };
                    l.AddRange(inner);
                    l.AddRange(rest);
                    yield return l;
                }
        }
    }

    private static double EnumeratedRegionProbability(string seq, int start, int end, double kelvin, int minLoop)
    {
        int n = seq.Length;
        var all = new List<(List<(int, int)> S, double E)>();
        foreach (var s in Structures(seq, 0, n, minLoop))
        {
            var db = new char[n];
            Array.Fill(db, '.');
            foreach (var (a, b) in s) { db[a] = '('; db[b] = ')'; }
            all.Add((s, CalculateFreeEnergy(seq, new string(db), kelvin)));
        }
        double emin = all.Min(e => e.E), rt = Rt(kelvin), z = 0, zOpen = 0;
        foreach (var (s, e) in all)
        {
            double w = Math.Exp(-(e - emin) / rt);
            z += w;
            if (s.All(p => (p.Item1 < start || p.Item1 > end) && (p.Item2 < start || p.Item2 > end))) zOpen += w;
        }
        return zOpen / z;
    }

    // Z_open/Z for u > 1 equals exhaustive enumeration (incl. repeated neighbours: NNDB-only terms, all in CalculateFreeEnergy).
    [TestCase(310.15, 3)]
    [TestCase(293.15, 3)]
    [TestCase(343.15, 3)]
    [TestCase(310.15, 5)]
    public void CalculateRegionUnpairedProbability_EqualsExhaustiveEnumeration(double kelvin, int minLoop)
    {
        var rng = new Random(20260925);
        var seqs = Enumerable.Range(0, 14).Select(_ => new string(Enumerable.Range(0, rng.Next(8, 15)).Select(__ => "ACGU"[rng.Next(4)]).ToArray()))
            .Concat(new[] { "GGGAAACCC", "GGGGAAAACCCCUUUU", "GGCCGAAAGGCC" });
        foreach (string seq in seqs)
        {
            int n = seq.Length;
            for (int len = 1; len <= Math.Min(n, 6); len++)
            {
                int end = rng.Next(len - 1, n);
                double expected = EnumeratedRegionProbability(seq, end - len + 1, end, kelvin, minLoop);
                double got = CalculateRegionUnpairedProbability(seq, end, len, minLoop, kelvin);
                Assert.That(got, Is.EqualTo(expected).Within(1e-9), $"{seq} [{end - len + 1}..{end}] T={kelvin} minLoop={minLoop}");
            }
        }
    }

    // Region probability is NOT the product of per-base probabilities (positions of one helix are strongly correlated):
    // joint probability of a hairpin's whole 5' arm exceeds the product; always bounded by the minimum and by Bonferroni.
    [Test]
    public void CalculateRegionUnpairedProbability_IsJointProbability_NotProductOfPerBase()
    {
        const string seq = "GGGGAAAACCCCUUUUGGGG";
        var pu = CalculateUnpairedProbabilities(seq).UnpairedProbabilities;
        for (int end = 3; end < seq.Length; end++)
        {
            double joint = CalculateRegionUnpairedProbability(seq, end, 4);
            double[] each = Enumerable.Range(end - 3, 4).Select(i => pu[i]).ToArray();
            Assert.That(joint, Is.LessThanOrEqualTo(each.Min() + 1e-12), $"joint <= min per-base @ {end}");
            Assert.That(joint, Is.GreaterThanOrEqualTo(1.0 - each.Sum(p => 1.0 - p) - 1e-12), $"Bonferroni lower bound @ {end}");
        }
        double armJoint = CalculateRegionUnpairedProbability(seq, 3, 4);
        double armProduct = Enumerable.Range(0, 4).Aggregate(1.0, (a, i) => a * pu[i]);
        Assert.That(armJoint, Is.GreaterThan(armProduct * 1.2), "correlated positions: joint > product of marginals");
    }

    #region Edge cases

    [Test]
    public void CalculateRegionUnpairedProbability_WholeSequence_IsInverseZ_ForSingleStructureCase()
    {
        // u = n: only the open chain qualifies → 1/Z for every sequence.
        foreach (string seq in new[] { "GAAAC", "CAAAAG", "GGGCGAAAGCCC", "GUCUCUGUGCGCAGCACAUG" })
        {
            double z = CalculatePartitionFunction(seq).PartitionFunction;
            Assert.That(CalculateRegionUnpairedProbability(seq, seq.Length - 1, seq.Length), Is.EqualTo(1.0 / z).Within(1e-12), seq);
        }
    }

    [Test]
    public void CalculateRegionUnpairedProbability_EndsOfSequence_AreValidAndConsistent()
    {
        const string seq = "GGGCGAAAGCCCAAGGGCGAAAGCCC";
        var pu = CalculateUnpairedProbabilities(seq).UnpairedProbabilities;
        Assert.Multiple(() =>
        {
            Assert.That(CalculateRegionUnpairedProbability(seq, 0, 1), Is.EqualTo(pu[0]).Within(1e-9), "u=1 at the first base");
            Assert.That(CalculateRegionUnpairedProbability(seq, seq.Length - 1, 1), Is.EqualTo(pu[^1]).Within(1e-9), "u=1 at the last base");
            Assert.That(CalculateRegionUnpairedProbability(seq, 5, 6), Is.InRange(0.0, pu[0]), "region at the 5' end");
            Assert.That(CalculateRegionUnpairedProbability(seq, seq.Length - 1, 6), Is.InRange(0.0, 1.0), "region at the 3' end");
        });
    }

    [Test]
    public void CalculateRegionUnpairedProbability_InvalidWindows_Throw()
    {
        const string seq = "GGGGAAAACCCC";
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, 11, 13), "longer than the sequence");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, 12, 1), "end beyond the last index");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, -1, 1), "negative end");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, 5, -2), "negative length");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, 5, 0), "zero length");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(string.Empty, 0, 1), "empty sequence has no window");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(null!, 0, 1), "null sequence has no window");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionUnpairedProbability(seq, 5, 3, temperature: 0.0), "non-positive temperature");
        });
    }

    [Test]
    public void CalculateRegionUnpairedProbability_DnaAndLowerCase_EqualRna()
    {
        const string rna = "GGGCGAAAGCCCAAGGGCGAAAGCCC";
        double expected = CalculateRegionUnpairedProbability(rna, 13, 8);
        Assert.Multiple(() =>
        {
            Assert.That(CalculateRegionUnpairedProbability(rna.Replace('U', 'T'), 13, 8), Is.EqualTo(expected), "T read as U");
            Assert.That(CalculateRegionUnpairedProbability(rna.ToLowerInvariant(), 13, 8), Is.EqualTo(expected), "case-insensitive");
        });
    }

    [Test]
    public void CalculateRegionUnpairedProbability_NonAcguBases_NeverPair()
    {
        // a region made of N never takes part in a pair → probability exactly 1 (Z_open = Z)
        Assert.That(CalculateRegionUnpairedProbability("GGGCNNNNGCCC", 7, 4), Is.EqualTo(1.0).Within(1e-12));
        // N's elsewhere only remove pairing options: stays a valid probability, deterministic
        double a = CalculateRegionUnpairedProbability("GGGCGAAAGCNC", 8, 3);
        Assert.That(a, Is.InRange(0.0, 1.0));
        Assert.That(CalculateRegionUnpairedProbability("GGGCGAAAGCNC", 8, 3), Is.EqualTo(a));
    }

    [Test]
    public void CalculateRegionUnpairedProbability_LargerMinLoop_NeverLessAccessible()
    {
        // a larger minimum hairpin size removes structures but never the open chain; for a pure-hairpin region it can only help
        const string seq = "GGGGAAAACCCC";
        double p3 = CalculateRegionUnpairedProbability(seq, 5, 4, minLoopSize: 3);
        double p5 = CalculateRegionUnpairedProbability(seq, 5, 4, minLoopSize: 5);
        double p11 = CalculateRegionUnpairedProbability(seq, 5, 4, minLoopSize: 11);
        Assert.Multiple(() =>
        {
            Assert.That(p5, Is.GreaterThan(p3));
            Assert.That(p11, Is.EqualTo(1.0).Within(1e-12), "no hairpin loop of >= 11 nt fits in 12 nt");
            Assert.That(CalculateRegionUnpairedProbability(seq, 5, 4, minLoopSize: 1), Is.EqualTo(p3), "minLoopSize below the NNDB minimum is clamped to 3");
        });
    }

    // Kertesz 2007 opening energy ΔG_open = -RT ln P(unpaired) = G_constrained - G_ensemble (ViennaRNA GASCONST 1.98717).
    [TestCase(310.15)]
    [TestCase(293.15)]
    public void CalculateRegionOpeningEnergy_IsMinusRtLogOfProbability(double kelvin)
    {
        const string seq = "GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA";
        double p = CalculateRegionUnpairedProbability(seq, 34, 14, temperature: kelvin);
        double g = CalculateRegionOpeningEnergy(seq, 34, 14, temperature: kelvin);
        Assert.That(g, Is.EqualTo(-Rt(kelvin) * Math.Log(p)).Within(1e-12));
        if (kelvin == 310.15)
            Assert.That(g, Is.EqualTo(-Rt(kelvin) * Math.Log(0.030532)).Within(0.03), "ViennaRNA constrained ensemble P = 0.030532");
    }

    [Test]
    public void CalculateRegionOpeningEnergy_NoPairPossible_IsZero()
    {
        Assert.That(CalculateRegionOpeningEnergy("GCGC", 3, 4), Is.EqualTo(0.0).Within(1e-12));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRegionOpeningEnergy("GCGC", 4, 1));
    }

    #endregion
}
