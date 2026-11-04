using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class OptimizeCodonsTests
{
    private static readonly CodonUsageTableInput EColi = new(Preset: "EColiK12");

    [Test]
    public void OptimizeCodons_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.optimize_codons("ATGTTATAA", EColi,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI));
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons("", EColi));
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(null!, EColi));
        // Unknown preset resolves to an error.
        Assert.Throws<ArgumentException>(() =>
            MolToolsTools.optimize_codons("ATGTAA", new CodonUsageTableInput(Preset: "Martian")));
        // GC bounds out of range / inverted.
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons("ATGTAA", EColi,
            gc_target_min: -0.1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons("ATGTAA", EColi,
            gc_target_min: 0.7, gc_target_max: 0.3));
    }

    [Test]
    public void OptimizeCodons_MaximizeCai_ReplacesRareLeucineWithCug()
    {
        // ATG(M) + three rare-in-E.coli Leucine codons (UUA 0.13, CUU 0.10, CUA 0.04) + UAA stop.
        // MaximizeCAI picks the most frequent synonymous codon: every Leu -> CUG (freq 0.50).
        // Met and stop are single/kept, so exactly 3 codons change and optimized CAI = 1.0.
        var r = MolToolsTools.optimize_codons("ATGTTACTTCTATAA", EColi,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.Multiple(() =>
        {
            Assert.That(r.OriginalSequence, Is.EqualTo("AUGUUACUUCUAUAA"));
            Assert.That(r.OptimizedSequence, Is.EqualTo("AUGCUGCUGCUGUAA"));
            Assert.That(r.ProteinSequence, Is.EqualTo("MLLL*"));
            Assert.That(r.OptimizedCAI, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(r.OriginalCAI, Is.LessThan(r.OptimizedCAI));
            Assert.That(r.GcContentOptimized, Is.EqualTo(7.0 / 15.0).Within(1e-9));
            Assert.That(r.ChangedCodons, Is.EqualTo(3));

            var changes = r.Changes.OrderBy(c => c.Position).ToList();
            Assert.That(changes.Select(c => c.Position), Is.EqualTo(new[] { 3, 6, 9 }));
            Assert.That(changes.All(c => c.Optimized == "CUG"), Is.True);
            Assert.That(changes.Select(c => c.Original), Is.EqualTo(new[] { "UUA", "CUU", "CUA" }));
        });
    }

    [Test]
    public void OptimizeCodons_MinimizeSecondary_RunsStructurePass_AndExposesSecondaryScoring()
    {
        // review 2026-09, CODON-OPT-001 F42 (B02 A1-8): MinimizeSecondary = balanced selection +
        // GC pass + the structure pass; it used to be an alias of BalancedOptimization.
        // 8 x GAA (Glu) + 8 x UUU (Phe): the balanced GAG/UUC blocks are complementary, so the
        // pair-count structure pass engages and changes the result.
        const string gene = "GAAGAAGAAGAAGAAGAAGAAGAAUUUUUUUUUUUUUUUUUUUUUUUU";

        var balanced = MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var minimizeSecondary = MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.MinimizeSecondary);

        Assert.Multiple(() =>
        {
            Assert.That(balanced.OptimizedSequence,
                Is.EqualTo("GAGGAGGAGGAGGAGGAGGAGGAGUUCUUCUUCUUCUUUUUUUUUUUU"));
            // F63: the structure pass carries the GC constraint of the same call, so the output
            // stays inside [0.40, 0.60] (GC 0.479167). Before F63 it was
            // GAGGAAGAAGAAGAAGAAGAAGAAUUCUUCUUCUUCUUCUUCUUCUUU (GC 0.3333, dnachisel
            // EnforceGCContent(0.40, 0.60).evaluate = -0.0666667).
            Assert.That(minimizeSecondary.OptimizedSequence,
                Is.EqualTo("GAGGAGGAGGAGGAGGAGGAGGAGUUCUUCUUCUUCUUCUUCUUCUUU"));
            Assert.That(minimizeSecondary.ProteinSequence, Is.EqualTo("EEEEEEEEFFFFFFFF"));
            Assert.That(CodonOptimizer.EvaluateGcContent(minimizeSecondary.OptimizedSequence, 0.40, 0.60).Score,
                Is.EqualTo(0).Within(1e-12));
        });

        // The Turner 2004 energy model is reachable through the tool's secondary_scoring parameter.
        var turner = MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
            secondary_scoring: CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe);
        Assert.That(turner.ProteinSequence, Is.EqualTo("EEEEEEEEFFFFFFFF"));
        Assert.That(turner.OptimizedSequence.Length, Is.EqualTo(gene.Length));
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
            secondary_scoring: (CodonOptimizer.SecondaryStructureScoring)42));
    }

    [Test]
    [Description("F43: gc_window constrains every sliding window (DNA Chisel EnforceGCContent(mini, maxi, window)); default stays global")]
    public void OptimizeCodons_GcWindow_ConstrainsEverySlidingWindow()
    {
        // M-(GS)x7-stop: starts above 60 % GC; every 20-nt window is reachable inside [0.40, 0.60].
        const string gene = "ATGGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCGGCAGCTAA";

        var global = MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var windowed = MolToolsTools.optimize_codons(gene, EColi,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization, gc_window: 20);

        Assert.Multiple(() =>
        {
            // dnachisel EnforceGCContent(0.40, 0.60, window=20).evaluate: -1.5 for the global
            // default output, 0 for the windowed one (input -6.4).
            Assert.That(CodonOptimizer.EvaluateGcContent(global.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(-1.5).Within(1e-12));
            Assert.That(CodonOptimizer.EvaluateGcContent(windowed.OptimizedSequence, 0.40, 0.60, window: 20).Score,
                Is.EqualTo(0).Within(1e-12));
            Assert.That(windowed.OptimizedSequence,
                Is.EqualTo("AUGGGCUCUGGUUCUGGUUCUGGCUCUGGUUCUGGUUCUGGCAGCUAA"));
            Assert.That(windowed.ProteinSequence, Is.EqualTo("MGSGSGSGSGSGSGS*"));
        });
    }

    [Test]
    [Description("F54 (A2-13): avoid_sites clears restriction sites created by codon selection; null/empty = historical output")]
    public void OptimizeCodons_AvoidSites_ClearsSitesAfterTheOtherPasses()
    {
        // MLQHMTMAK*: MaximizeCAI creates PstI (CUG+CAG), NdeI (CAU+AUG) and NcoI (ACC+AUG).
        const string gene = "ATGCTACAACACATGACAATGGCAAAGTAA";
        string[] sites = { "CTGCAG", "CATATG", "CCATGG" };

        var plain = MolToolsTools.optimize_codons(gene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var avoided = MolToolsTools.optimize_codons(gene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            avoid_sites: sites);
        var empty = MolToolsTools.optimize_codons(gene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            avoid_sites: Array.Empty<string>());

        Assert.Multiple(() =>
        {
            foreach (string site in sites)
                Assert.That(plain.OptimizedSequence.Replace('U', 'T'), Does.Contain(site), site);
            // dnachisel 3.2.16: AvoidPattern(site).evaluate scores 0 (passes) on this output for all
            // three sites, EnforceTranslation passes, and dnachisel's own solver returns exactly
            // this sequence (CodonOptimize objective -1.3892) - review 2026-09 F54.
            Assert.That(avoided.OptimizedSequence, Is.EqualTo("AUGCUGCAACACAUGACGAUGGCGAAAUAA"));
            Assert.That(avoided.ProteinSequence, Is.EqualTo("MLQHMTMAK*"));
            Assert.That(avoided.OptimizedCAI, Is.EqualTo(0.793311).Within(1e-6));
            Assert.That(avoided.ChangedCodons, Is.EqualTo(4));
            Assert.That(empty.OptimizedSequence, Is.EqualTo(plain.OptimizedSequence));
            // Site strings are validated like remove_restriction_sites.
            Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(
                gene, EColi, avoid_sites: new[] { "GAATTZ" }));
        });
    }

    [Test]
    [Description("F55 (A2-14): the strategy parses under both spellings, case-insensitively")]
    public void OptimizeCodons_AvoidRareCodons_BothSpellingsAccepted()
    {
        const string gene = "ATGCTGCTAAGAGAGTTTTAA";
        var legacy = MolToolsTools.optimize_codons(gene, EColi, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons);
        var fixedSpelling = MolToolsTools.optimize_codons(gene, EColi, CodonOptimizer.OptimizationStrategy.AvoidRareCodons);

        Assert.Multiple(() =>
        {
            Assert.That(fixedSpelling.OptimizedSequence, Is.EqualTo(legacy.OptimizedSequence));
            foreach (string name in new[] { "AvoidRareCodons", "avoidrarecodons", "AvoidRareCodeons", "avoidRareCodeons" })
                Assert.That(
                    System.Text.Json.JsonSerializer.Deserialize<CodonOptimizer.OptimizationStrategy>(
                        $"\"{name}\"",
                        new System.Text.Json.JsonSerializerOptions
                        {
                            Converters = { new System.Text.Json.Serialization.JsonStringEnumConverter() },
                        }),
                    Is.EqualTo(CodonOptimizer.OptimizationStrategy.AvoidRareCodeons), name);
        });
    }

    [Test]
    [Description("F43: gc_window must be positive")]
    public void OptimizeCodons_GcWindow_MustBePositive()
    {
        Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(
            "ATGGGCAGCTAA", EColi, CodonOptimizer.OptimizationStrategy.BalancedOptimization, gc_window: 0));
    }
}
