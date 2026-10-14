using static Seqeron.Genomics.MolTools.CodonOptimizer;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the codon-optimiser / rare-codon behaviours added by the B02
/// finisher (docs/Validation/review-2026-09/B02.md): Sherlocc "≤" pauses and conserved clusters
/// (F38/F40), %MinMax usage basis (F39), Turner 2004 structure scoring and the MinimizeSecondary
/// structure pass (F41/F42), the windowed GC constraint (F43), the genetic-code-aware optimiser (F48)
/// and restriction-site avoidance chained in OptimizeSequence, GC-aware since F60 (F54/F60).
/// Fixed-seed random inputs; MFE-based relations use tiny sequences (the fold is the cost).
///
/// Test Units: CODON-OPT-001, CODON-RARE-001, CODON-CAI-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Translation")]
public class B02FinisherMetamorphicTests
{
    private const int Seed = 20261014;

    private static string RandomDna(Random rng, int codons) =>
        new(Enumerable.Range(0, 3 * codons).Select(_ => "ACGT"[rng.Next(4)]).ToArray());

    private static string Protein(string cds, GeneticCode? code = null) =>
        Translator.Translate(cds, code ?? GeneticCode.Standard).Sequence;

    #region F38 / F40 — Sherlocc

    // A codon alignment of homologues: one random ancestor, per-row point mutations and gap codons.
    private static List<string> RandomAlignment(Random rng, int rows, int codons)
    {
        string ancestor = RandomDna(rng, codons);
        var alignment = new List<string>();
        for (int r = 0; r < rows; r++)
        {
            var chars = ancestor.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < 0.15) chars[i] = "ACGT"[rng.Next(4)];
            for (int c = 0; c < codons; c++)
                if (rng.NextDouble() < 0.08) { chars[3 * c] = chars[3 * c + 1] = chars[3 * c + 2] = '-'; }
            alignment.Add(new string(chars));
        }
        return alignment;
    }

    /// <summary>
    /// F40: the Sherlocc analysis is invariant to the order of the aligned rows and to the gap spelling
    /// ('-' vs '.'); raising the threshold only adds pause columns.
    /// </summary>
    [Test]
    public void Sherlocc_RowOrderAndGapSpelling_Invariant_ThresholdMonotone()
    {
        var rng = new Random(Seed);
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 60; trial++)
            {
                var alignment = RandomAlignment(rng, rng.Next(2, 6), rng.Next(8, 40));
                double threshold = 8 + rng.Next(0, 15);
                var baseResult = FindConservedRareCodonClusters(alignment, Yeast, threshold);

                var shuffled = alignment.OrderBy(_ => rng.Next()).ToList();
                var shuffledResult = FindConservedRareCodonClusters(shuffled, Yeast, threshold);
                Assert.That(shuffledResult.PauseColumns, Is.EqualTo(baseResult.PauseColumns), $"trial {trial}: pauses vs row order");
                Assert.That(shuffledResult.Clusters.Select(c => (c.StartColumn, c.EndColumn, c.PauseCount, c.MiddlePoint)),
                    Is.EqualTo(baseResult.Clusters.Select(c => (c.StartColumn, c.EndColumn, c.PauseCount, c.MiddlePoint))),
                    $"trial {trial}: clusters vs row order");
                Assert.That(shuffledResult.ColumnUsageCounts, Is.EqualTo(baseResult.ColumnUsageCounts), $"trial {trial}");
                for (int i = 0; i < baseResult.ColumnAverages.Count; i++)
                {
                    double? a = baseResult.ColumnAverages[i], b = shuffledResult.ColumnAverages[i];
                    Assert.That(a.HasValue, Is.EqualTo(b.HasValue), $"trial {trial} column {i}");
                    if (a.HasValue) Assert.That(b!.Value, Is.EqualTo(a.Value).Within(1e-9), $"trial {trial} column {i}");
                }

                var dotted = alignment.Select(s => s.Replace('-', '.')).ToList();
                var dottedResult = FindConservedRareCodonClusters(dotted, Yeast, threshold);
                Assert.That(dottedResult.ColumnAverages, Is.EqualTo(baseResult.ColumnAverages), $"trial {trial}: gap spelling");
                Assert.That(dottedResult.Clusters, Is.EqualTo(baseResult.Clusters), $"trial {trial}: gap spelling");

                var higher = FindConservedRareCodonClusters(alignment, Yeast, threshold + 5);
                Assert.That(higher.PauseColumns, Is.SupersetOf(baseResult.PauseColumns), $"trial {trial}: threshold monotone");
            }
        });
    }

    /// <summary>
    /// F38: the Sherlocc "≤" pause test only adds pauses to the strict "&lt;" test — every strict cluster lies
    /// inside an inclusive one and every strict rare codon is an inclusive rare codon; a threshold equal to a
    /// tabulated frequency is where the two differ.
    /// </summary>
    [Test]
    public void RareCodons_InclusiveThreshold_IsMonotoneSuperset()
    {
        var rng = new Random(Seed + 1);
        var frequencies = EColiK12.CodonFrequencies.Values.Distinct().OrderBy(v => v).ToArray();
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 150; trial++)
            {
                string cds = RandomDna(rng, rng.Next(5, 80));
                // Thresholds at tabulated frequencies exercise the ≤ / < boundary.
                double threshold = frequencies[rng.Next(frequencies.Length)];
                int window = rng.Next(1, 9), minRare = rng.Next(1, window + 1);

                var strict = FindRareCodonClusters(cds, EColiK12, threshold, window, minRare, inclusiveThreshold: false);
                var inclusive = FindRareCodonClusters(cds, EColiK12, threshold, window, minRare, inclusiveThreshold: true);
                foreach (var c in strict)
                    Assert.That(inclusive.Any(i => i.StartCodon <= c.StartCodon && c.EndCodon <= i.EndCodon && i.RareCount >= 0),
                        Is.True, $"trial {trial}: strict cluster {c} not covered");
                Assert.That(inclusive.Sum(c => c.RareCount), Is.GreaterThanOrEqualTo(0));

                var strictRare = FindRareCodons(cds, EColiK12, threshold, inclusiveThreshold: false).Select(r => r.Position).ToHashSet();
                var inclusiveRare = FindRareCodons(cds, EColiK12, threshold, inclusiveThreshold: true).Select(r => r.Position).ToHashSet();
                Assert.That(inclusiveRare.IsSupersetOf(strictRare), Is.True, $"trial {trial}: rare codons");
                Assert.That(FindRareCodons(cds, EColiK12, threshold).Select(r => r.Position), Is.EquivalentTo(strictRare),
                    $"trial {trial}: FindRareCodons default stays strict");
            }
        });
    }

    #endregion

    #region F39 — %MinMax usage basis

    /// <summary>
    /// F39: without per-thousand usage (E. coli K-12 preset) both bases give the same profile; every %MinMax
    /// value lies in [−100, 100] for every preset and basis; the default overload is PerThousandWhenAvailable.
    /// </summary>
    [Test]
    public void MinMax_Basis_EqualWithoutPerThousand_AndBounded()
    {
        var rng = new Random(Seed + 2);
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 60; trial++)
            {
                string cds = RandomDna(rng, rng.Next(5, 60));
                int window = rng.Next(1, 20);
                Assert.That(CalculateMinMaxProfile(cds, EColiK12, window, MinMaxUsageBasis.PerThousandWhenAvailable),
                    Is.EqualTo(CalculateMinMaxProfile(cds, EColiK12, window, MinMaxUsageBasis.CodonFrequencies)), $"trial {trial}");
                foreach (var table in new[] { EColiK12, Yeast, Human })
                {
                    var dflt = CalculateMinMaxProfile(cds, table, window);
                    Assert.That(dflt, Is.EqualTo(CalculateMinMaxProfile(cds, table, window, MinMaxUsageBasis.PerThousandWhenAvailable)));
                    foreach (var basis in new[] { MinMaxUsageBasis.PerThousandWhenAvailable, MinMaxUsageBasis.CodonFrequencies })
                        foreach (var w in CalculateMinMaxProfile(cds, table, window, basis))
                            Assert.That(Math.Abs(w.PercentMinMax), Is.LessThanOrEqualTo(100.0 + 1e-9), $"trial {trial} {table.OrganismName} {basis}");
                }
            }
        });
    }

    #endregion

    #region F41 / F42 — Turner 2004 structure scoring

    /// <summary>
    /// F41/F42: Turner-scored ReduceSecondaryStructure and the MinimizeSecondary strategy (which now runs the
    /// structure pass) preserve the protein and are deterministic. Tiny inputs — every window is an MFE fold.
    /// </summary>
    [Test]
    public void TurnerScoring_PreservesProtein_AndIsDeterministic()
    {
        var rng = new Random(Seed + 3);
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 6; trial++)
            {
                string cds = RandomDna(rng, rng.Next(10, 15));
                string a = ReduceSecondaryStructure(cds, EColiK12, windowSize: 24,
                    scoring: SecondaryStructureScoring.Turner2004Mfe, maxFoldEvaluations: 300);
                string b = ReduceSecondaryStructure(cds, EColiK12, windowSize: 24,
                    scoring: SecondaryStructureScoring.Turner2004Mfe, maxFoldEvaluations: 300);
                Assert.That(b, Is.EqualTo(a), $"trial {trial}: deterministic");
                Assert.That(Protein(a), Is.EqualTo(Protein(cds)), $"trial {trial}: protein");

                if (trial < 3)
                {
                    var r1 = OptimizeSequence(cds, EColiK12, OptimizationStrategy.MinimizeSecondary,
                        secondaryScoring: SecondaryStructureScoring.Turner2004Mfe);
                    var r2 = OptimizeSequence(cds, EColiK12, OptimizationStrategy.MinimizeSecondary,
                        secondaryScoring: SecondaryStructureScoring.Turner2004Mfe);
                    Assert.That(r2.OptimizedSequence, Is.EqualTo(r1.OptimizedSequence), $"trial {trial}: optimise deterministic");
                    Assert.That(Protein(r1.OptimizedSequence), Is.EqualTo(Protein(cds)), $"trial {trial}: optimise protein");
                    Assert.That(r1.ChangedCodons, Is.EqualTo(r1.Changes.Count));
                }
            }
        });
    }

    #endregion

    #region F43 — windowed GC

    /// <summary>
    /// F43: the windowed GC repair starts from the same codon choice as MaximizeCAI and only accepts
    /// breach-reducing rounds, so its DNA Chisel breach score is never worse; the protein is preserved; the
    /// evaluation's breach windows are exactly the windows outside the bounds (direct count).
    /// </summary>
    [Test]
    public void WindowedGc_NeverIncreasesBreach_PreservesProtein()
    {
        var rng = new Random(Seed + 4);
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 80; trial++)
            {
                string cds = RandomDna(rng, rng.Next(10, 50));
                int window = new[] { 9, 12, 18, 30 }[rng.Next(4)];
                double mini = 0.35 + 0.05 * rng.Next(0, 3), maxi = mini + 0.1 + 0.05 * rng.Next(0, 3);

                string selected = OptimizeSequence(cds, EColiK12, OptimizationStrategy.MaximizeCAI).OptimizedSequence;
                var repaired = OptimizeSequence(cds, EColiK12, OptimizationStrategy.BalancedOptimization,
                    gcTargetMin: mini, gcTargetMax: maxi, gcWindow: window);

                var before = EvaluateGcContent(selected, mini, maxi, window);
                var after = EvaluateGcContent(repaired.OptimizedSequence, mini, maxi, window);
                Assert.That(after.Score, Is.GreaterThanOrEqualTo(before.Score - 1e-12), $"trial {trial}");
                Assert.That(Protein(repaired.OptimizedSequence), Is.EqualTo(Protein(cds)), $"trial {trial}: protein");

                // Breach windows = windows whose G+C fraction is outside [mini, maxi] (direct count).
                string s = repaired.OptimizedSequence;
                var expected = Enumerable.Range(0, Math.Max(0, s.Length - window + 1))
                    .Where(i => { double gc = (double)s.Substring(i, window).Count(c => c is 'G' or 'C') / window; return gc < mini || gc > maxi; })
                    .ToList();
                Assert.That(after.BreachWindowStarts, Is.EqualTo(expected), $"trial {trial}: breach windows");
                Assert.That(after.Score == 0, Is.EqualTo(expected.Count == 0), $"trial {trial}: score sign");
            }
        });
    }

    #endregion

    #region F48 — genetic-code-aware optimiser

    /// <summary>F48: every strategy preserves the protein under the given genetic code, for every NCBI table.</summary>
    [Test]
    public void Optimizer_PreservesProtein_UnderEveryGeneticCode()
    {
        var rng = new Random(Seed + 5);
        var strategies = new[]
        {
            OptimizationStrategy.MaximizeCAI, OptimizationStrategy.BalancedOptimization,
            OptimizationStrategy.HarmonizeExpression, OptimizationStrategy.AvoidRareCodons,
        };
        Assert.Multiple(() =>
        {
            foreach (int table in GeneticCode.SupportedTableNumbers)
            {
                var code = GeneticCode.GetByTableNumber(table);
                for (int trial = 0; trial < 4; trial++)
                {
                    string cds = RandomDna(rng, rng.Next(5, 40));
                    var strategy = strategies[rng.Next(strategies.Length)];
                    var result = OptimizeSequence(cds, EColiK12, strategy, geneticCode: code);
                    Assert.That(Protein(result.OptimizedSequence, code), Is.EqualTo(Protein(cds, code)),
                        $"table {table} {strategy}: protein");
                    Assert.That(result.ProteinSequence, Is.EqualTo(Protein(cds, code)), $"table {table}: reported protein");

                    string reduced = RemoveRestrictionSites(cds, new[] { "GAATTC", "GGATCC" }, EColiK12, code);
                    Assert.That(Protein(reduced, code), Is.EqualTo(Protein(cds, code)), $"table {table}: site removal");
                }
            }
        });
    }

    #endregion

    #region F54 / F60 — avoidSites in OptimizeSequence

    private static readonly string[] Sites = { "GAATTC", "GGATCC", "AAGCTT", "GCNGC" };

    // Occurrences (forward strand, both site orientations) of the IUPAC sites in a DNA-spelled sequence.
    private static List<(int Position, string Pattern)> Occurrences(string sequence, IEnumerable<string> sites)
    {
        string dna = SequenceCase.ToUpperDna(sequence);
        var patterns = sites.SelectMany(s => new[] { s, DnaSequence.GetReverseComplementString(s) }).Distinct();
        var hits = new List<(int, string)>();
        foreach (var p in patterns)
            for (int i = 0; i + p.Length <= dna.Length; i++)
                if (Enumerable.Range(0, p.Length).All(k => IupacHelper.MatchesIupac(dna[i + k], p[k])))
                    hits.Add((i, p));
        return hits;
    }

    /// <summary>
    /// F54: with a GC-free strategy the chained avoidance equals the public RemoveRestrictionSites applied to
    /// the optimised sequence; the protein is preserved; every occurrence left is not removable by any single
    /// synonymous codon replacement (the pass's definition of removable).
    /// </summary>
    [Test]
    public void AvoidSites_EqualsPublicPass_RemovesRemovable_PreservesProtein()
    {
        var rng = new Random(Seed + 6);
        var code = GeneticCode.Standard;
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 80; trial++)
            {
                string cds = RandomDna(rng, rng.Next(6, 40));
                if (trial % 2 == 0) cds = cds.Insert(3 * rng.Next(0, cds.Length / 3), "GAATTCGGATCC");
                var strategy = trial % 3 == 0 ? OptimizationStrategy.HarmonizeExpression : OptimizationStrategy.MaximizeCAI;

                var plain = OptimizeSequence(cds, EColiK12, strategy);
                var chained = OptimizeSequence(cds, EColiK12, strategy, avoidSites: Sites);
                Assert.That(chained.OptimizedSequence,
                    Is.EqualTo(RemoveRestrictionSites(plain.OptimizedSequence, Sites, EColiK12)), $"trial {trial}: chain");
                Assert.That(Protein(chained.OptimizedSequence), Is.EqualTo(Protein(cds)), $"trial {trial}: protein");

                string result = chained.OptimizedSequence;
                var codons = Enumerable.Range(0, result.Length / 3).Select(i => result.Substring(3 * i, 3)).ToList();
                foreach (var (position, pattern) in Occurrences(result, Sites))
                {
                    for (int ci = position / 3; ci <= Math.Min((position + pattern.Length - 1) / 3, codons.Count - 1); ci++)
                    {
                        char aa = code.Translate(codons[ci]);
                        if (aa == '*') continue;
                        foreach (var alt in code.GetCodonsForAminoAcid(aa).Select(c => c.Replace('T', 'U')))
                        {
                            if (alt == codons[ci]) continue;
                            var copy = codons.ToList();
                            copy[ci] = alt;
                            bool stillThere = Occurrences(string.Concat(copy), new[] { pattern }).Any(o => o.Position == position && o.Pattern == pattern);
                            Assert.That(stillThere, Is.True,
                                $"trial {trial}: site {pattern}@{position} removable by codon {ci} → {alt}");
                        }
                    }
                }
            }
        });
    }

    /// <summary>
    /// F60: from the same GC-balanced sequence holding one removable site occurrence, a single GC-aware
    /// removal step never ends with a larger DNA Chisel GC breach than the GC-blind step (the GC-aware ranking
    /// prefers acceptable — or least-breach — candidates and agrees with the blind choice when that one is
    /// acceptable). Inputs where either pass needed more than one rewrite are skipped.
    /// </summary>
    [Test]
    public void GcAwareSiteRemoval_NeverWorseBreachThanGcBlind()
    {
        var rng = new Random(Seed + 7);
        var site = new[] { "GAATTC" };
        int compared = 0;
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 500; trial++)
            {
                string cds = RandomDna(rng, rng.Next(8, 30));
                cds = cds.Insert(3 * rng.Next(0, cds.Length / 3 + 1), rng.Next(2) == 0 ? "GAATTC" : "AGAATTCAA");
                int window = new[] { 12, 18, 24 }[rng.Next(3)];
                double mini = 0.45, maxi = 0.55;

                string balanced = OptimizeSequence(cds, EColiK12, OptimizationStrategy.BalancedOptimization,
                    gcTargetMin: mini, gcTargetMax: maxi, gcWindow: window).OptimizedSequence;
                if (Occurrences(balanced, site).Count != 1) continue;

                string aware = OptimizeSequence(cds, EColiK12, OptimizationStrategy.BalancedOptimization,
                    gcTargetMin: mini, gcTargetMax: maxi, gcWindow: window, avoidSites: site).OptimizedSequence;
                string blind = RemoveRestrictionSites(balanced, site, EColiK12);
                if (CodonDiff(balanced, aware) != 1 || CodonDiff(balanced, blind) != 1
                    || Occurrences(aware, site).Count != 0 || Occurrences(blind, site).Count != 0)
                    continue;

                compared++;
                double awareBreach = -EvaluateGcContent(aware, mini, maxi, window).Score;
                double blindBreach = -EvaluateGcContent(blind, mini, maxi, window).Score;
                Assert.That(awareBreach, Is.LessThanOrEqualTo(blindBreach + 1e-9), $"trial {trial}");
                Assert.That(Protein(aware), Is.EqualTo(Protein(cds)), $"trial {trial}: protein");
            }
        });
        Assert.That(compared, Is.GreaterThan(10), "too few single-step comparisons generated");
    }

    private static int CodonDiff(string a, string b) =>
        Enumerable.Range(0, a.Length / 3).Count(i => a.Substring(3 * i, 3) != b.Substring(3 * i, 3));

    #endregion
}
