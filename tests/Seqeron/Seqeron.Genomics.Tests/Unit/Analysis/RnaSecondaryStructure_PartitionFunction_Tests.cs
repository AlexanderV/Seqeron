// RNA-PARTITION-001 — RNA Partition Function (McCaskill) & Boltzmann Structure Probability
// Evidence: docs/Evidence/RNA-PARTITION-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PARTITION-001.md
// Source: McCaskill JS (1990). Biopolymers 29(6-7):1105-1119. DOI 10.1002/bip.360290621 (PMID 1695107).
//         Lorenz R et al. (2011) ViennaRNA Package 2.0, Algorithms Mol Biol 6:26; ViennaRNA 2.7.2 sdist
//         (partfunc/, probabilities/, params/params.c RESCALE_dG); Mathews et al. (2004) PNAS 101:7287 (Turner 2004).

using Seqeron.Genomics.Core;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_PartitionFunction_Tests
{
    // ViennaRNA GASCONST 1.98717 cal/(mol·K) (params/constants.h); RT in kcal/mol.
    private static double Rt(double kelvin) => 1.98717 * kelvin / 1000.0;
    private const double Rt37 = 1.98717 * 310.15 / 1000.0;

    private static IEnumerable<string> RandomSequences(int count, int minLen, int maxLen, int seed, string alphabet = "ACGU")
    {
        var rng = new Random(seed);
        for (int k = 0; k < count; k++)
            yield return new string(Enumerable.Range(0, rng.Next(minLen, maxLen + 1)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
    }

    // ───────────────────────────────────────────────────────────────────────────────────────────
    // Independent oracle: exhaustive enumeration of every pseudoknot-free structure with hairpins ≥ 3,
    // each scored by the public structure evaluator CalculateFreeEnergy (loop-by-loop Turner 2004 sum, a
    // code path that shares nothing with the inside/outside recursions). Z = Σ exp(−E/RT) and
    // P(i,j) = Σ_{S∋(i,j)} exp(−E/RT)/Z are then exact by definition.
    // ───────────────────────────────────────────────────────────────────────────────────────────
    private static IEnumerable<List<(int, int)>> Structures(string seq, int i, int j)
    {
        if (i >= j) { yield return new List<(int, int)>(); yield break; }
        foreach (var s in Structures(seq, i + 1, j)) yield return s; // i unpaired
        for (int k = i + 4; k < j; k++)
        {
            if (!CanPair(seq[i], seq[k])) continue;
            foreach (var inner in Structures(seq, i + 1, k))
                foreach (var rest in Structures(seq, k + 1, j))
                {
                    var l = new List<(int, int)> { (i, k) };
                    l.AddRange(inner);
                    l.AddRange(rest);
                    yield return l;
                }
        }
    }

    private static (double Z, double[,] P, double Mfe) BruteForce(string seq, double kelvin)
    {
        int n = seq.Length;
        double rt = Rt(kelvin);
        var energies = new List<(List<(int, int)> S, double E)>();
        foreach (var s in Structures(seq, 0, n))
        {
            var db = new char[n];
            Array.Fill(db, '.');
            foreach (var (a, b) in s) { db[a] = '('; db[b] = ')'; }
            energies.Add((s, CalculateFreeEnergy(seq, new string(db), kelvin)));
        }
        double emin = energies.Min(e => e.E);
        double z = 0;
        var p = new double[n, n];
        foreach (var (s, e) in energies)
        {
            double w = Math.Exp(-(e - emin) / rt);
            z += w;
            foreach (var (a, b) in s) p[a, b] += w;
        }
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++) p[a, b] /= z;
        return (z * Math.Exp(-emin / rt), p, emin);
    }

    #region Exact ensemble: brute-force enumeration

    // M3 — Z, every base-pair probability, MFE and P(MFE) equal exhaustive enumeration over all structures.
    [TestCase(310.15)]
    [TestCase(293.15)]
    [TestCase(343.15)]
    public void CalculatePartitionFunction_EqualsExhaustiveEnumeration(double kelvin)
    {
        // short random sequences incl. repeated neighbours (NNDB-only terms: bulge degeneracy, special C, GU closure, all-C loops)
        foreach (string seq in RandomSequences(24, 6, 15, 1234).Concat(new[] { "GGGAAACCC", "GGGGCCCC", "GGGAAAACCCUUUGGGAAACCC", "GGCCGAAAGGCCAGGCCGAAAGGCC" }))
        {
            var (z, p, mfe) = BruteForce(seq, kelvin);
            var r = CalculatePartitionFunction(seq, temperature: kelvin);

            Assert.That(r.PartitionFunction / z, Is.EqualTo(1.0).Within(1e-9), $"Z of {seq}");
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-Rt(kelvin) * Math.Log(z)).Within(1e-9), $"G of {seq}");
            Assert.That(r.MfeFreeEnergy, Is.EqualTo(mfe).Within(1e-9), $"MFE of {seq} = min over all structures");
            for (int i = 0; i < seq.Length; i++)
                for (int j = i + 1; j < seq.Length; j++)
                {
                    double got = r.BasePairProbabilities.TryGetValue((i, j), out double v) ? v : 0.0;
                    Assert.That(got, Is.EqualTo(p[i, j]).Within(1e-9), $"P({i},{j}) of {seq}");
                }
        }
    }

    #endregion

    #region ViennaRNA 2.7.2 locked values

    // ViennaRNA 2.7.2, RNA.fold_compound(seq, md).mfe()/pf()/bpp()/pr_structure()/mean_bp_distance()/centroid(),
    // md.dangles = 2, md.temperature = Celsius, NNDB-precision parameter file (RNA.params_save() with the NNDB stacks, TerminalAU 0.45
    // and the 27 NNDB-vs-ViennaRNA table differences; see Evidence §12) loaded with RNA.params_load — the parameter set of this engine.
    // Sequences contain no NNDB-only terms that ViennaRNA lacks. Unpaired = P(unpaired) at positions 0 / n/2 / n−1; Top = the two most
    // probable pairs. Vienna's pf() deviates from the exact Boltzmann sum over its own eval_structure energies by ≤ 0.017 kcal/mol
    // (SMOOTH() of exterior dangle weights, params.c); at ≠ 37 °C its MFE uses int-truncated parameters (TRUNC_MAYBE).
    public sealed record ViennaCase(string Seq, double Celsius, string Structure, double Mfe, double G, double PMfe,
        double Diversity, string Centroid, double[] Unpaired, (int I, int J, double P)[] Top);

    private static readonly ViennaCase[] ViennaCases =
    {
        new("GGGAAAACCC", 37, "(((....)))", -2.02, -2.1529, 0.8060, 1.3195, "(((....)))", new[] { 0.0593, 1.0000, 0.1762 }, new[] { (1, 8, 0.8376), (2, 7, 0.8267) }),
        new("GGGCGAAAGCCC", 37, "((((....))))", -6.64, -6.6562, 0.9740, 0.0520, "((((....))))", new[] { 0.0251, 1.0000, 0.0252 }, new[] { (2, 9, 0.9999), (1, 10, 0.9997) }),
        new("GGGGAAAACCCC", 37, "((((....))))", -5.28, -5.3941, 0.8310, 1.7016, "((((....))))", new[] { 0.0300, 1.0000, 0.1506 }, new[] { (2, 9, 0.8638), (1, 10, 0.8637) }),
        new("GCGGAAACGCCGAAAGGCGC", 37, ".......((((....)))).", -6.84, -6.8736, 0.9470, 0.1876, ".......((((....)))).", new[] { 0.9932, 0.0052, 0.9971 }, new[] { (9, 16, 0.9939), (8, 17, 0.9936) }),
        new("GGCGAUCAGGCCAAAGGCCUGAUCGCC", 37, "((((((((((((...))))))))))))", -23.57, -23.5964, 0.9580, 0.0842, "((((((((((((...))))))))))))", new[] { 0.0250, 1.0000, 0.0251 }, new[] { (9, 17, 1.0000), (8, 18, 1.0000) }),
        new("GGGCGAAAGCCCAAGGGCGAAAGCCCAAGGGCGAAAGCCC", 37, "((((....((((..((((....))))..))))....))))", -25.02, -25.0805, 0.9065, 2.1116, "((((....((((..((((....))))..))))....))))", new[] { 0.0235, 1.0000, 0.0247 }, new[] { (16, 23, 0.9932), (15, 24, 0.9932) }),
        new("GCACACGCUGACUGUGAGAGACUCAGCACUCUCACACACUGCUGUACUGAUCACUCAUAGCUGAGUAGUGUACUCG", 37, "((....)).....(((((((........)))))))((((((((((..(((....)))..))..)))))))).....", -21.71, -23.1022, 0.1045, 21.4880, ".............(((((((........))))))).(((((((....................)))))))......", new[] { 0.4246, 0.2726, 0.9827 }, new[] { (17, 30, 0.8327), (15, 32, 0.8327) }),
        new("AUAUAAAUAUAUUUAUAUAAAUUUAUAUAU", 37, "((((((((.(((....))).))))))))..", -3.95, -4.2444, 0.6203, 3.1347, "((((((((.(((....))).))))))))..", new[] { 0.1711, 0.9142, 0.9917 }, new[] { (4, 23, 0.9406), (3, 24, 0.9406) }),
        new("GGGAAAACCC", 20, "(((....)))", -3.39, -3.4935, 0.8373, 1.1602, "(((....)))", new[] { 0.0165, 1.0000, 0.1389 }, new[] { (1, 8, 0.8686), (2, 7, 0.8626) }),
        new("GGGGAAAACCCC", 20, "((((....))))", -7.20, -7.3078, 0.8310, 1.6038, "((((....))))", new[] { 0.0141, 1.0000, 0.1368 }, new[] { (1, 10, 0.8708), (2, 9, 0.8707) }),
        new("GGGCGAAAGCCCAAGGGCGAAAGCCCAAGGGCGAAAGCCC", 20, "((((....((((..((((....))))..))))....))))", -32.11, -32.3368, 0.6775, 6.0337, "((((....((((..((((....))))..))))....))))", new[] { 0.0093, 1.0000, 0.0111 }, new[] { (16, 23, 0.9928), (15, 24, 0.9928) }),
        new("GCACACGCUGACUGUGAGAGACUCAGCACUCUCACACACUGCUGUACUGAUCACUCAUAGCUGAGUAGUGUACUCG", 20, "(.(((((((.((((((((.((.((((.((..............)).)))))).))))))).).))).)))).)...", -31.01, -32.3158, 0.1063, 26.0678, ".......((..(((((((.((.((((.((..............)).)))))).)))))))...))...........", new[] { 0.2628, 0.6201, 0.9953 }, new[] { (15, 55, 0.5866), (16, 54, 0.5866) }),
        new("GGGAAAACCC", 60, "(((....)))", -0.14, -0.6425, 0.4681, 1.8597, ".((....)).", new[] { 0.4304, 1.0000, 0.4902 }, new[] { (1, 8, 0.5207), (2, 7, 0.5070) }),
        new("GGGCGAAAGCCC", 60, "((((....))))", -3.34, -3.4054, 0.9059, 0.1837, "((((....))))", new[] { 0.0642, 1.0000, 0.0652 }, new[] { (2, 9, 0.9924), (1, 10, 0.9908) }),
        new("GGGGAAAACCCC", 60, "((((....))))", -2.64, -2.8454, 0.7333, 2.0637, "((((....))))", new[] { 0.0959, 1.0000, 0.1909 }, new[] { (2, 9, 0.8271), (1, 10, 0.8267) }),
        new("GGGCGAAAGCCCAAGGGCGAAAGCCCAAGGGCGAAAGCCC", 60, "((((....((((..((((....))))..))))....))))", -15.29, -15.4272, 0.8129, 0.7398, "((((....((((..((((....))))..))))....))))", new[] { 0.0584, 1.0000, 0.0607 }, new[] { (15, 24, 0.9944), (16, 23, 0.9944) }),
        new("GCACACGCUGACUGUGAGAGACUCAGCACUCUCACACACUGCUGUACUGAUCACUCAUAGCUGAGUAGUGUACUCG", 60, "............((((((((........))))))))(((((((....(((....)))......)))))))......", -10.79, -12.8344, 0.0456, 15.9708, "............((((((((........))))))))(((((((....(((....)))......)))))))......", new[] { 0.8195, 0.2423, 0.8670 }, new[] { (15, 32, 0.9223), (17, 30, 0.9220) }),
    };

    public static IEnumerable<TestCaseData> ViennaData => ViennaCases.Select(c => new TestCaseData(c).SetName($"Vienna_{c.Seq.Length}nt_{c.Celsius:0}C_{c.Seq[..Math.Min(8, c.Seq.Length)]}"));

    [TestCaseSource(nameof(ViennaData))]
    public void CalculatePartitionFunction_MatchesViennaRna(ViennaCase c)
    {
        double kelvin = c.Celsius + 273.15;
        var r = CalculatePartitionFunction(c.Seq, temperature: kelvin);
        int n = c.Seq.Length;
        bool ref37 = Math.Abs(c.Celsius - 37) < 1e-9;

        Assert.Multiple(() =>
        {
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(c.G).Within(0.02), "ensemble free energy (pf)");
            Assert.That(r.MfeFreeEnergy, Is.EqualTo(c.Mfe).Within(ref37 ? 1e-9 : 0.15), "MFE");
            Assert.That(r.MfeProbability, Is.EqualTo(c.PMfe).Within(0.03 * Math.Max(c.PMfe, 0.1) + (ref37 ? 0.0 : 0.15 * c.PMfe)), "P(MFE structure)");
            Assert.That(r.EnsembleDiversity, Is.EqualTo(c.Diversity).Within(0.03 + 0.02 * c.Diversity), "mean base-pair distance");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(c.Unpaired[0]).Within(0.01));
            Assert.That(r.UnpairedProbabilities[n / 2], Is.EqualTo(c.Unpaired[1]).Within(0.01));
            Assert.That(r.UnpairedProbabilities[n - 1], Is.EqualTo(c.Unpaired[2]).Within(0.01));
            foreach (var (i, j, p) in c.Top)
                Assert.That(r.BasePairProbabilities[(i, j)], Is.EqualTo(p).Within(0.01), $"P({i},{j})");
            if (ref37)
            {
                Assert.That(r.MfeStructure, Is.EqualTo(c.Structure), "MFE structure");
                Assert.That(r.CentroidStructure, Is.EqualTo(c.Centroid), "centroid structure");
            }
        });
    }

    // ViennaRNA GGGAAAACCC sanity value: MFE −2.02 kcal/mol for (((....))) (NNDB-precision parameters; the default file gives −2.10).
    [Test]
    public void CalculatePartitionFunction_SanityHairpin_GGGAAAACCC()
    {
        var r = CalculatePartitionFunction("GGGAAAACCC");
        Assert.Multiple(() =>
        {
            Assert.That(r.MfeStructure, Is.EqualTo("(((....)))"));
            Assert.That(r.MfeFreeEnergy, Is.EqualTo(-2.02).Within(1e-9));
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-2.1529).Within(0.005), "ViennaRNA pf() −2.1529");
            Assert.That(r.PartitionFunction, Is.EqualTo(Math.Exp(2.1529 / Rt37)).Within(0.01 * Math.Exp(2.1529 / Rt37)));
        });
    }

    #endregion

    #region Normalisation and internal consistency

    // M6d / M7 — probabilities in [0,1]; rows of the base-pair-probability matrix sum to ≤ 1;
    // P(unpaired i) + Σ_j P(i,j) = 1; G ≤ MFE; Z = exp(ln Z) = exp(−G/RT).
    [Test]
    public void CalculatePartitionFunction_Normalisation_HoldsForRandomSequences()
    {
        foreach (string seq in RandomSequences(40, 10, 120, 99))
        {
            var r = CalculatePartitionFunction(seq);
            int n = seq.Length;
            var rows = new double[n];
            foreach (var ((i, j), p) in r.BasePairProbabilities)
            {
                Assert.That(p, Is.InRange(0.0, 1.0), $"P({i},{j}) of {seq}");
                rows[i] += p;
                rows[j] += p;
            }
            for (int i = 0; i < n; i++)
            {
                Assert.That(rows[i], Is.LessThanOrEqualTo(1.0 + 1e-9), $"row sum {i} of {seq}");
                Assert.That(rows[i] + r.UnpairedProbabilities[i], Is.EqualTo(1.0).Within(1e-9), $"P(unpaired)+ΣP(pairs) at {i} of {seq}");
            }
            Assert.That(r.EnsembleFreeEnergy, Is.LessThanOrEqualTo(r.MfeFreeEnergy + 1e-9), "G ≤ MFE");
            Assert.That(r.PartitionFunction, Is.GreaterThanOrEqualTo(1.0 - 1e-12), "INV-1: Z ≥ 1");
            Assert.That(r.LogPartitionFunction, Is.EqualTo(Math.Log(r.PartitionFunction)).Within(1e-9));
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-Rt37 * r.LogPartitionFunction).Within(1e-9));
        }
    }

    // P(MFE structure) = exp((G − E_mfe)/RT) (ViennaRNA pr_structure) and = CalculateStructureProbability(E_mfe, G).
    [Test]
    public void CalculatePartitionFunction_MfeProbability_IsBoltzmannFactorOfMfeStructure()
    {
        foreach (string seq in RandomSequences(20, 12, 90, 5))
        {
            var r = CalculatePartitionFunction(seq);
            double expected = Math.Exp((r.EnsembleFreeEnergy - r.MfeFreeEnergy) / Rt37);
            Assert.That(r.MfeProbability, Is.EqualTo(expected).Within(1e-12), seq);
            Assert.That(r.MfeProbability, Is.EqualTo(CalculateStructureProbability(r.MfeFreeEnergy, r.EnsembleFreeEnergy)).Within(1e-12));
            Assert.That(CalculateFreeEnergy(seq, r.MfeStructure), Is.EqualTo(r.MfeFreeEnergy).Within(1e-9), "MfeStructure evaluates to MfeFreeEnergy");
            Assert.That(r.MfeProbability, Is.InRange(0.0, 1.0));
        }
    }

    // Ensemble diversity = 2·Σ_{i<j} P(1−P); centroid = pairs with P > 0.5 (ViennaRNA mean_bp_distance / centroid).
    [Test]
    public void CalculatePartitionFunction_DiversityAndCentroid_FollowTheirDefinitions()
    {
        foreach (string seq in RandomSequences(15, 15, 100, 321))
        {
            var r = CalculatePartitionFunction(seq);
            double d = 2 * r.BasePairProbabilities.Values.Sum(p => p * (1 - p));
            Assert.That(r.EnsembleDiversity, Is.EqualTo(d).Within(1e-9));
            var expected = new string('.', seq.Length).ToCharArray();
            foreach (var ((i, j), p) in r.BasePairProbabilities.Where(kv => kv.Value > 0.5)) { expected[i] = '('; expected[j] = ')'; }
            Assert.That(r.CentroidStructure, Is.EqualTo(new string(expected)));
        }
    }

    // CalculateUnpairedProbabilities is the same engine reduced to the legacy record.
    [Test]
    public void CalculateUnpairedProbabilities_EqualsPartitionFunctionFields()
    {
        const string seq = "GGGCGAAAGCCCAAGGGCGAAAGCCC";
        var a = CalculatePartitionFunction(seq);
        var b = CalculateUnpairedProbabilities(seq);
        Assert.Multiple(() =>
        {
            Assert.That(b.PartitionFunction, Is.EqualTo(a.PartitionFunction));
            Assert.That(b.EnsembleFreeEnergy, Is.EqualTo(a.EnsembleFreeEnergy));
            Assert.That(b.UnpairedProbabilities, Is.EqualTo(a.UnpairedProbabilities));
            Assert.That(b.BasePairProbabilities, Is.EqualTo(a.BasePairProbabilities));
        });
    }

    #endregion

    #region Input contracts

    // M10
    [Test]
    public void CalculatePartitionFunction_NullSequence_Throws()
        => Assert.Throws<ArgumentNullException>(() => CalculatePartitionFunction(null!));

    // M11 — empty: only the empty structure, Z = 1, G = 0.
    [Test]
    public void CalculatePartitionFunction_EmptySequence_ReturnsOneNoPairs()
    {
        var r = CalculatePartitionFunction("");
        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0));
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(0.0));
            Assert.That(r.BasePairProbabilities, Is.Empty);
            Assert.That(r.UnpairedProbabilities, Is.Empty);
            Assert.That(r.MfeStructure, Is.Empty);
            Assert.That(r.MfeProbability, Is.EqualTo(1.0));
        });
    }

    // M1/M2 — no admissible pair (AAAA) or only a pair of span ≤ 3 (GC, GAAC): Z = 1, every base unpaired with probability 1.
    [TestCase("AAAA")]
    [TestCase("GC")]
    [TestCase("GAAC")]
    [TestCase("G")]
    public void CalculatePartitionFunction_NoAdmissiblePair_ReturnsOne(string seq)
    {
        var r = CalculatePartitionFunction(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(r.BasePairProbabilities, Is.Empty);
            Assert.That(r.UnpairedProbabilities, Is.All.EqualTo(1.0));
            Assert.That(r.MfeStructure, Is.EqualTo(new string('.', seq.Length)));
        });
    }

    // M5 — GAAAAC: exactly two structures (open chain, (0,5)), Z = 1 + exp(−E/RT) with E = hairpin(GC, AAAA) = 4.5 (+ nothing else).
    [Test]
    public void CalculatePartitionFunction_SingleAdmissiblePair_ProbabilityIsBoltzmann()
    {
        var r = CalculatePartitionFunction("GAAAAC");
        double e = CalculateFreeEnergy("GAAAAC", "(....)");
        double w = Math.Exp(-e / Rt37);
        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0 + w).Within(1e-12));
            Assert.That(r.BasePairProbabilities, Has.Count.EqualTo(1));
            Assert.That(r.BasePairProbabilities[(0, 5)], Is.EqualTo(w / (1.0 + w)).Within(1e-12));
        });
    }

    // DNA T, lower case: same ensemble as the upper-case RNA sequence.
    [Test]
    public void CalculatePartitionFunction_TAndLowerCase_AreReadAsUpperCaseRna()
    {
        var u = CalculatePartitionFunction("GGGCGAAAGCCCAAGGGCGAAAGCCC");
        var t = CalculatePartitionFunction("gggcgaaagcccaagggcgaaagccc".Replace('u', 't'));
        var d = CalculatePartitionFunction("GGGCGAAAGCCCAAGGGCGAAAGCCC".Replace('U', 'T'));
        Assert.Multiple(() =>
        {
            Assert.That(t.PartitionFunction, Is.EqualTo(u.PartitionFunction));
            Assert.That(d.PartitionFunction, Is.EqualTo(u.PartitionFunction));
            Assert.That(CalculatePartitionFunction("GGGAAAUCCC").PartitionFunction, Is.EqualTo(CalculatePartitionFunction("GGGAAATCCC").PartitionFunction).Within(1e-12));
        });
    }

    // N / IUPAC / gaps never pair (ViennaRNA alphabet): an N-substituted base removes its pairs and matches the ensemble with that base unpairable.
    [Test]
    public void CalculatePartitionFunction_UnknownBases_NeverPair()
    {
        var r = CalculatePartitionFunction("GGGNAAAACCC");
        foreach (var ((i, j), _) in r.BasePairProbabilities)
        {
            Assert.That("GGGNAAAACCC"[i], Is.Not.EqualTo('N'));
            Assert.That("GGGNAAAACCC"[j], Is.Not.EqualTo('N'));
        }
        Assert.That(CalculatePartitionFunction("NNNNNNNNNN").PartitionFunction, Is.EqualTo(1.0));
        Assert.That(CalculatePartitionFunction("RYKMSWRYKM").PartitionFunction, Is.EqualTo(1.0), "IUPAC ambiguity codes do not pair");
        Assert.DoesNotThrow(() => CalculatePartitionFunction("GGGAA-ACCC.R"));
    }

    // minLoopSize: values < 3 are clamped to 3 (NNDB); a larger value forbids short hairpins.
    [Test]
    public void CalculatePartitionFunction_MinLoopSize_ClampedAndHonoured()
    {
        Assert.That(CalculatePartitionFunction("GGGAAAACCC", minLoopSize: 0).PartitionFunction,
            Is.EqualTo(CalculatePartitionFunction("GGGAAAACCC", minLoopSize: 3).PartitionFunction));
        var wide = CalculatePartitionFunction("GGGAAAACCC", minLoopSize: 5);
        Assert.Multiple(() =>
        {
            Assert.That(wide.BasePairProbabilities.Keys.All(k => k.J - k.I - 1 >= 5), Is.True, "no hairpin loop shorter than 5 nt");
            Assert.That(wide.EnsembleFreeEnergy, Is.GreaterThan(CalculatePartitionFunction("GGGAAAACCC").EnsembleFreeEnergy));
        });
    }

    // S2 — non-positive / non-finite temperature is rejected, never a NaN.
    [TestCase(0.0)]
    [TestCase(-5.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void CalculatePartitionFunction_InvalidTemperature_Throws(double t)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculatePartitionFunction("GGGAAACCC", temperature: t));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateStructureProbability(-5, -6, t));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMinimumFreeEnergy("GGGAAACCC", temperature: t));
    }

    #endregion

    #region Temperature (Turner 2004 enthalpy rescaling)

    // ViennaRNA params.c RESCALE_dG: dG(T) = dH − (dH − dG37)·T/310.15 with the Turner 2004 enthalpies (stackdH).
    // C·G→C·G stack: dG37 = −3.26 (NNDB; ViennaRNA −3.30), dH = −13.40 kcal/mol.
    [TestCase(293.15)]
    [TestCase(333.15)]
    public void TurnerModel_Stack_FollowsEnthalpyRescaling(double kelvin)
    {
        var m = Seqeron.Genomics.Analysis.TurnerModel.ForTemperature(kelvin);
        double f = kelvin / 310.15;
        double expected = -13.40 - (-13.40 - -3.26) * f;
        Assert.That(m.Stack(1, 1), Is.EqualTo(expected).Within(1e-12));
        Assert.That(Seqeron.Genomics.Analysis.TurnerModel.At37.Stack(1, 1), Is.EqualTo(-3.26));
        // RT follows R·T (not only RT): rescaling is not a rescaling of 1/RT alone.
        Assert.That(m.RT, Is.EqualTo(Rt(kelvin)).Within(1e-15));
    }

    // ViennaRNA fc.eval_structure at 20 / 60 °C (md.temperature, NNDB-precision file): hairpin −3.39 / −0.14 (int-truncated parameters;
    // ours keeps the unrounded rescaled values).
    [Test]
    public void CalculateFreeEnergy_AtOtherTemperatures_MatchesViennaRna()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFreeEnergy("GGGAAAACCC", "(((....)))", 293.15), Is.EqualTo(-3.39).Within(0.02));
            Assert.That(CalculateFreeEnergy("GGGAAAACCC", "(((....)))", 333.15), Is.EqualTo(-0.14).Within(0.02));
            Assert.That(CalculateMinimumFreeEnergy("GGGAAAACCC", temperature: 293.15), Is.EqualTo(CalculateMfeStructure("GGGAAAACCC", temperature: 293.15).FreeEnergy));
        });
    }

    // A stable hairpin ensemble is more stable when cold; the default temperature is 37 °C.
    [Test]
    public void CalculatePartitionFunction_FreeEnergy_IncreasesWithTemperature_AndDefaultIs37C()
    {
        const string seq = "GGGGAAAACCCC";
        double g10 = CalculatePartitionFunction(seq, temperature: 283.15).EnsembleFreeEnergy;
        double g37 = CalculatePartitionFunction(seq, temperature: 310.15).EnsembleFreeEnergy;
        double g60 = CalculatePartitionFunction(seq, temperature: 333.15).EnsembleFreeEnergy;
        Assert.Multiple(() =>
        {
            Assert.That(g10, Is.LessThan(g37));
            Assert.That(g37, Is.LessThan(g60));
            Assert.That(CalculatePartitionFunction(seq).EnsembleFreeEnergy, Is.EqualTo(g37));
        });
    }

    #endregion

    #region Numerical stability and speed

    // Long random sequences: no overflow / NaN, normalisation holds, ln Z finite (ViennaRNA pf_scale: the matrices are scaled per nucleotide).
    [Test]
    public void CalculatePartitionFunction_Length500_IsStableAndNormalised()
    {
        string seq = GenerateRandomRna(500, new Random(7));
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = CalculatePartitionFunction(seq);
        sw.Stop();

        var rows = new double[seq.Length];
        foreach (var ((i, j), p) in r.BasePairProbabilities) { rows[i] += p; rows[j] += p; }
        Assert.Multiple(() =>
        {
            Assert.That(double.IsFinite(r.EnsembleFreeEnergy), Is.True);
            Assert.That(double.IsFinite(r.LogPartitionFunction), Is.True);
            Assert.That(r.EnsembleFreeEnergy, Is.LessThanOrEqualTo(r.MfeFreeEnergy + 1e-9));
            Assert.That(rows.Max(), Is.LessThanOrEqualTo(1.0 + 1e-9));
            for (int i = 0; i < rows.Length; i++)
                Assert.That(rows[i] + r.UnpairedProbabilities[i], Is.EqualTo(1.0).Within(1e-9));
            Assert.That(r.BasePairProbabilities.Values.All(double.IsFinite), Is.True);
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(120.0), "O(n³): n = 500 completes well within the budget even in Debug");
        });
    }

    // Extreme temperatures: ensembles stay well-formed; at a few Kelvin the ensemble is the ground state (no NaN).
    [Test]
    public void CalculatePartitionFunction_ExtremeTemperatures_AreFinite()
    {
        foreach (double t in new[] { 0.001, 1.0, 5.0, 100.0, 1000.0, 1e6 })
        {
            var r = CalculatePartitionFunction("GGGGAAAACCCC", temperature: t);
            Assert.That(double.IsFinite(r.EnsembleFreeEnergy), Is.True, $"G at {t} K");
            Assert.That(double.IsNaN(r.PartitionFunction), Is.False, $"Z at {t} K");
            Assert.That(r.BasePairProbabilities.Values.All(p => double.IsFinite(p) && p >= 0 && p <= 1.0 + 1e-9), Is.True, $"P at {t} K");
        }
    }

    #endregion

    #region CalculateStructureProbability

    // M8 — E_struct = E_ensemble → p = 1.
    [Test]
    public void CalculateStructureProbability_EnergyEqualsEnsemble_ReturnsOne()
        => Assert.That(CalculateStructureProbability(structureEnergy: -5.0, ensembleEnergy: -5.0), Is.EqualTo(1.0).Within(1e-12));

    // M9 — p = exp(−(E_s − E_ens)/RT), RT = 1.98717e-3·310.15 (ViennaRNA kT): E_s = −5, E_ens = −6 → exp(−1/RT) = 0.19735...
    [Test]
    public void CalculateStructureProbability_KnownEnergies_ReturnsBoltzmannRatio()
    {
        double p = CalculateStructureProbability(structureEnergy: -5.0, ensembleEnergy: -6.0);
        Assert.That(p, Is.EqualTo(Math.Exp(-1.0 / Rt37)).Within(1e-12));
        Assert.That(p, Is.EqualTo(0.19735).Within(5e-5));
    }

    // Evaluated in the exponent: huge |E| neither overflows nor yields NaN (exp(−E/RT)/exp(−G/RT) would be ∞/∞).
    [Test]
    public void CalculateStructureProbability_LargeEnergies_NoOverflow()
    {
        Assert.That(CalculateStructureProbability(-2000.0, -2000.5), Is.EqualTo(Math.Exp(-0.5 / Rt37)).Within(1e-12));
        Assert.That(CalculateStructureProbability(-2000.0, -2000.5), Is.InRange(0.0, 1.0));
    }

    // Temperature enters through RT: a higher T flattens the Boltzmann factor.
    [Test]
    public void CalculateStructureProbability_TemperatureDependence()
        => Assert.That(CalculateStructureProbability(-5, -6, 340.0), Is.GreaterThan(CalculateStructureProbability(-5, -6, 280.0)));

    #endregion

    #region GenerateRandomRna

    // M12
    [Test]
    public void GenerateRandomRna_SameSeed_IsDeterministicAndValidAlphabet()
    {
        string a = GenerateRandomRna(50, new Random(42));
        string b = GenerateRandomRna(50, new Random(42));

        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(b));
            Assert.That(a, Has.Length.EqualTo(50));
            Assert.That(a.All(c => "ACGU".Contains(c)), Is.True);
        });
    }

    // S3
    [Test]
    public void GenerateRandomRna_BoundaryLengths_ProduceCorrectOutput()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GenerateRandomRna(0, new Random(1)), Is.Empty);
            Assert.That(GenerateRandomRna(100, new Random(1)), Has.Length.EqualTo(100));
            Assert.That(GenerateRandomRna(100, gcContent: 0.5), Has.Length.EqualTo(100));
        });
    }

    [Test]
    public void GenerateRandomRna_NullRandom_Throws()
        => Assert.Throws<ArgumentNullException>(() => GenerateRandomRna(10, (Random)null!));

    // Delegate smoke test: parameterless-seed overload returns a valid sequence of the right length.
    [Test]
    public void GenerateRandomRna_DefaultOverload_ReturnsValidSequence()
    {
        string s = GenerateRandomRna(20);
        Assert.Multiple(() =>
        {
            Assert.That(s, Has.Length.EqualTo(20));
            Assert.That(s.All(c => "ACGU".Contains(c)), Is.True);
        });
    }

    #endregion
}
