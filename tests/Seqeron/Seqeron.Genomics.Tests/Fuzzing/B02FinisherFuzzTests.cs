using static Seqeron.Genomics.MolTools.CodonOptimizer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier fuzz tests for the entry points added by the B02 finisher
/// (docs/Validation/review-2026-09/B02.md, F31–F61, WP24/WP25). Random UTF-16 garbage (incl. U+017F 'ſ',
/// U+0131 'ı', U+212A Kelvin, surrogate halves, '\0'), random IUPAC/X nucleotide text and random option
/// values are fed to every new surface; each call must return a value satisfying its documented contract or
/// throw the documented exception (ArgumentException family) — never hang and never leak another runtime
/// exception (IndexOutOfRange, KeyNotFound, NullReference, InvalidOperation, …).
///
/// Test Units: TRANS-PROT-001, TRANS-CODON-001, TRANS-SIXFRAME-001, CODON-ENC-001, CODON-RSCU-001,
/// CODON-STATS-001, CODON-OPT-001, CODON-RARE-001, B02-SWEEP.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Translation")]
public class B02FinisherFuzzTests
{
    private const int Seed = 20261014;

    private static string Garbage(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
        {
            chars[i] = rng.Next(6) switch
            {
                0 => (char)rng.Next(0, 0x10000),
                1 => (char)rng.Next(0, 128),
                2 => "ſKı\u0000￿𐀀"[rng.Next(7)],
                _ => "ACGTUNacgtunRYSWKMBDHVXx-.*"[rng.Next(27)],
            };
        }
        return new string(chars);
    }

    private static string Nucleotides(Random rng, int length, string alphabet = "ACGTUacgtuRYSWKMBDHVNXx") =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static GeneticCode RandomCode(Random rng) =>
        GeneticCode.GetByTableNumber(GeneticCode.SupportedTableNumbers[rng.Next(GeneticCode.SupportedTableNumbers.Count)]);

    // Runs the call; only the ArgumentException family is a documented failure.
    private static void Documented(Action call, string context)
    {
        try { call(); }
        catch (ArgumentException) { }
        catch (Exception e)
        {
            Assert.Fail($"{context}: undocumented {e.GetType().Name}: {e.Message}");
        }
    }

    /// <summary>F31/F32/F33 + WP24/25: protein, codon and SequenceCase entry points on garbage.</summary>
    [Test]
    public void ProteinAndCodonSurfaces_Garbage_DocumentedOutcomesOnly()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 400; trial++)
        {
            string s = Garbage(rng, rng.Next(0, 40));
            var code = RandomCode(rng);
            Documented(() =>
            {
                bool ok = ProteinSequence.TryCreate(s, out var p);
                if (ok)
                {
                    Assert.That(p!.Sequence.All(c => c < 128), Is.True, "accepted protein holds a non-ASCII char");
                    _ = p.ToThreeLetterCode();
                    _ = p.MolecularWeight();
                    _ = p.FindMotif("A").ToList();
                }
            }, $"protein {trial}");
            Documented(() => _ = ProteinPhysicochemistry.MolecularWeight(s), $"MW {trial}");
            Documented(() => _ = ProteinPhysicochemistry.Gravy(s), $"GRAVY {trial}");
            Documented(() => _ = ProteinPhysicochemistry.IsoelectricPoint(s), $"pI {trial}");
            Documented(() => _ = Translator.Translate(s, code, rng.Next(0, 3), false, rng.Next(2) == 0), $"translate {trial}");
            if (s.Length >= 3)
            {
                string codon = s[..3];
                Documented(() => _ = code.Translate(codon, allowNucleotideX: true), $"codon {trial}");
                Documented(() => _ = code.IsStartCodon(codon, true), $"start {trial}");
                Documented(() => _ = code.IsStopCodon(codon, true), $"stop {trial}");
            }
            Documented(() => _ = code.GetCodonsForAminoAcid(s.Length > 0 ? s[0] : 'A').ToList(), $"codons-for {trial}");

            string upper = SequenceCase.ToUpperAscii(s);
            Assert.That(upper.Length, Is.EqualTo(s.Length));
            Assert.That(upper.Where((c, i) => c != s[i]).All(c => c is >= 'A' and <= 'Z'), Is.True, $"case {trial}");
        }

        for (int n = -10; n < 50; n++)
            Documented(() => _ = GeneticCode.GetByTableNumberOrThrow(n, "genetic_code"), $"table {n}");
    }

    /// <summary>F44–F46, F53, F58: getorf modes / circular / transeq options on IUPAC + garbage input and random options.</summary>
    [Test]
    public void SixFrameAndGetOrf_RandomInputsAndOptions_DocumentedOutcomesOnly()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = rng.Next(4) == 0 ? Garbage(rng, rng.Next(0, 30)) : Nucleotides(rng, rng.Next(0, 60));
            var code = RandomCode(rng);
            var mode = (OrfFindMode)rng.Next(-1, 8);
            int minSize = rng.Next(-2, 40), maxSize = rng.Next(-2, 200), flanking = rng.Next(-2, 15);
            bool circular = rng.Next(2) == 0, both = rng.Next(2) == 0, met = rng.Next(2) == 0;
            Documented(() =>
            {
                var records = Translator.FindOrfs(s, mode, code, minSize, maxSize, flanking, both, met, circular).ToList();
                int len = s.Length;
                foreach (var r in records)
                {
                    Assert.That(r.StartPosition, Is.InRange(0, Math.Max(0, len - 1)), $"getorf {trial}: start");
                    Assert.That(r.EndPosition, Is.GreaterThanOrEqualTo(r.StartPosition - 1), $"getorf {trial}: end");
                    if (!circular) Assert.That(r.EndPosition, Is.LessThan(len), $"getorf {trial}: linear end");
                    _ = r.ToGetOrfRange(len);
                }
            }, $"getorf {trial} mode {(int)mode} table {code.TableNumber}");
            Documented(() =>
            {
                foreach (var o in Translator.FindOrfs(s, code, rng.Next(-1, 10), both))
                    _ = o.ToGetOrfRange(s.Length);
            }, $"legacy orfs {trial}");

            var convention = (SixFrameConvention)rng.Next(-1, 4);
            Documented(() =>
            {
                var frames = Translator.TranslateSixFrames(s, code, convention);
                Assert.That(frames.Keys, Is.EquivalentTo(new[] { 1, 2, 3, -1, -2, -3 }));
            }, $"six-frame {trial}");

            var regions = rng.Next(3) == 0
                ? null
                : Enumerable.Range(0, rng.Next(0, 3)).Select(_ => (rng.Next(-1, 40), rng.Next(-1, 80))).ToList();
            var options = new TranseqOptions { Trim = rng.Next(2) == 0, Clean = rng.Next(2) == 0, Methionine = rng.Next(2) == 0, Regions = regions };
            Documented(() =>
            {
                var frames = Translator.TranslateSixFrames(s, code, convention, options);
                foreach (var p in frames.Values)
                {
                    if (options.Clean) Assert.That(p.Sequence, Does.Not.Contain("*"), $"transeq {trial}: clean");
                    if (options.Trim) Assert.That(p.Sequence.TrimEnd('X', '*'), Is.EqualTo(p.Sequence), $"transeq {trial}: trim");
                }
            }, $"transeq {trial}");
        }
    }

    /// <summary>F36/F37/F50–F52/F57/F59/F61: codon-usage statistics on garbage and random option values.</summary>
    [Test]
    public void CodonUsageSurfaces_RandomInputsAndOptions_DocumentedOutcomesOnly()
    {
        var rng = new Random(Seed + 2);
        var methods = Enum.GetValues<EncMethod>();
        for (int trial = 0; trial < 300; trial++)
        {
            string s = rng.Next(3) == 0 ? Garbage(rng, rng.Next(0, 60)) : Nucleotides(rng, rng.Next(0, 240), "ACGTACGTUacgtN");
            var code = RandomCode(rng);
            var options = new EncOptions
            {
                Method = rng.Next(8) == 0 ? (EncMethod)99 : methods[rng.Next(methods.Length)],
                Pseudocount = rng.Next(4) == 0 ? new[] { -1.0, double.NaN, 0.0, 0.5, 2.0 }[rng.Next(5)] : null,
                BackgroundCorrection = rng.Next(4) == 0 ? rng.Next(2) == 0 : null,
                WeightedAverage = rng.Next(4) == 0 ? rng.Next(2) == 0 : null,
                Background = rng.Next(5) == 0 ? Nucleotides(rng, rng.Next(0, 30), "ACGTN") : null,
                CodonWCap = rng.Next(4) == 0,
            };
            Documented(() =>
            {
                bool ok = CodonUsageAnalyzer.TryCalculateEnc(s, code, options, out double enc);
                Assert.That(ok ? enc > 0 && double.IsFinite(enc) : enc == 0, Is.True, $"enc {trial}: {enc}");
            }, $"enc {trial} {options}");
            Documented(() =>
            {
                bool ok = CodonUsageAnalyzer.TryCalculateEnc(s, code, out double enc);
                Assert.That(ok ? enc > 0 : enc == 0, Is.True);
            }, $"enc-plain {trial}");

            double pseudocount = new[] { -0.5, double.NaN, double.PositiveInfinity, 0.0, 1e-9, 1.0 }[rng.Next(6)];
            var level = (RscuFamilyLevel)rng.Next(-1, 3);
            Documented(() =>
            {
                var rscu = CodonUsageAnalyzer.CalculateRscu(s, code, pseudocount, level);
                Assert.That(rscu.Values.All(v => v >= 0 && double.IsFinite(v)), Is.True, $"rscu {trial}");
            }, $"rscu {trial}");

            var convention = (CodonGcConvention)rng.Next(-1, 3);
            var set = CodonWOptimalCodonSets.All[rng.Next(CodonWOptimalCodonSets.All.Count)];
            Documented(() =>
            {
                var st = CodonUsageAnalyzer.GetStatistics(s, code, convention, rng.Next(2) == 0 ? set.OptimalCodons : null);
                foreach (double v in new[] { st.Gc1, st.Gc2, st.Gc3, st.Gc3s })
                    Assert.That(v, Is.InRange(0.0, 100.0), $"stats {trial}");
                if (st.Fop is double fop) Assert.That(fop, Is.InRange(0.0, 1.0));
                if (st.Cbi is double cbi) Assert.That(cbi, Is.LessThanOrEqualTo(1.0 + 1e-12));
            }, $"stats {trial}");

            var badCodons = rng.Next(4) == 0 ? new HashSet<string> { Garbage(rng, rng.Next(0, 5)) } : null;
            Documented(() => _ = CodonUsageAnalyzer.TryCalculateFop(s, out _, badCodons ?? set.OptimalCodons, code,
                rng.Next(2) == 0 ? set.RareCodons : null), $"fop {trial}");
            Documented(() => _ = CodonUsageAnalyzer.TryCalculateCbi(s, out _, badCodons ?? set.OptimalCodons, code), $"cbi {trial}");
            Documented(() => _ = CodonWOptimalCodonSets.TryGetByName(rng.Next(2) == 0 ? null : Garbage(rng, 8), out _), $"set {trial}");
        }
    }

    /// <summary>F38–F40, F43, F48, F54: optimiser / rare-codon entry points on garbage and random parameters.</summary>
    [Test]
    public void OptimizerSurfaces_RandomInputsAndOptions_DocumentedOutcomesOnly()
    {
        var rng = new Random(Seed + 3);
        var tables = new[] { EColiK12, Yeast, Human };
        for (int trial = 0; trial < 150; trial++)
        {
            string s = rng.Next(3) == 0 ? Garbage(rng, rng.Next(0, 45)) : Nucleotides(rng, rng.Next(0, 90), "ACGTACGTUacgN");
            var table = tables[rng.Next(tables.Length)];
            var code = rng.Next(2) == 0 ? null : RandomCode(rng);
            int? window = rng.Next(3) == 0 ? null : rng.Next(-2, 40);
            var sites = rng.Next(2) == 0 ? null
                : Enumerable.Range(0, rng.Next(0, 3)).Select(_ => rng.Next(4) == 0 ? Garbage(rng, 4) : Nucleotides(rng, rng.Next(0, 7), "ACGTNRY")).ToList();
            var strategy = (OptimizationStrategy)rng.Next(0, 5);
            if (strategy == OptimizationStrategy.MinimizeSecondary) strategy = OptimizationStrategy.BalancedOptimization; // MFE cost

            Documented(() =>
            {
                var r = OptimizeSequence(s, table, strategy, gcWindow: window, geneticCode: code, avoidSites: sites);
                Assert.That(r.OptimizedSequence.Length, Is.EqualTo(r.OriginalSequence.Length), $"optimize {trial}");
                Assert.That(r.ChangedCodons, Is.EqualTo(r.Changes.Count));
            }, $"optimize {trial}");
            Documented(() =>
            {
                var e = EvaluateGcContent(s, rng.NextDouble() * 0.5, 0.5 + rng.NextDouble() * 0.5, window, rng.Next(-5, 60));
                Assert.That(e.Score, Is.LessThanOrEqualTo(0));
                Assert.That(e.Score == 0, Is.EqualTo(e.BreachWindowStarts.Count == 0));
            }, $"evaluate {trial}");
            Documented(() => _ = FindRareCodonClusters(s, table, rng.NextDouble() * 0.4, rng.Next(-1, 10), rng.Next(-1, 6), rng.Next(2) == 0), $"clusters {trial}");
            Documented(() => _ = FindRareCodons(s, table, rng.NextDouble() * 0.4, rng.Next(2) == 0, code).ToList(), $"rare {trial}");
            Documented(() => _ = CalculateMinMaxProfile(s, table, rng.Next(-1, 25), (MinMaxUsageBasis)rng.Next(-1, 3), code), $"minmax {trial}");

            int rows = rng.Next(0, 4), columns = rng.Next(0, 12);
            var alignment = Enumerable.Range(0, rows)
                .Select(_ => rng.Next(6) == 0 ? Garbage(rng, 3 * columns) : Nucleotides(rng, 3 * columns + (rng.Next(8) == 0 ? 1 : 0), "ACGT-.N"))
                .ToList();
            var opts = new SherloccOptions
            {
                AveragingWindow = rng.Next(-1, 9), ClusterWindow = rng.Next(-1, 9), MinPauses = rng.Next(-1, 6),
                MinCoverage = rng.NextDouble() * 1.2 - 0.1, UseCodonFrequencies = rng.Next(2) == 0,
            };
            Documented(() =>
            {
                var r = FindConservedRareCodonClusters(alignment, table, rng.NextDouble() * 20, opts);
                Assert.That(r.ColumnAverages.Count, Is.EqualTo(r.AlignmentLength), $"sherlocc {trial}");
                Assert.That(r.ClusterFraction, Is.InRange(0.0, 1.0));
            }, $"sherlocc {trial}");
        }
    }
}
