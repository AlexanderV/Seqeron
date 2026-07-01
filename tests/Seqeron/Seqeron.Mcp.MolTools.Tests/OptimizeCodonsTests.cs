using System;
using System.Linq;
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
}
