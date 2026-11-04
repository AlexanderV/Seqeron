using System.Diagnostics;
using System.Reflection;
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// The MCP-visible side of B08 R14: remove_restriction_sites / optimize_codons(avoid_sites) cost time linear in the sequence length (the tool
/// needed 8.2 s for 48 000 nt with the site N), a site longer than 64 characters or with a non-IUPAC character is refused whatever the sequence
/// is (an empty or sub-codon one used to skip the check), and the 64-character limit is stated in the descriptions the caller reads.
/// </summary>
[TestFixture]
public class RemoveRestrictionSitesLimitTests
{
    private static CodonUsageTableInput EColi() => new(Preset: "EColiK12");

    private static readonly string TooLong = new('N', RestrictionAnalyzer.MaxRecognitionLength + 1);

    private static string RandomCodingDna(int codons, int seed)
    {
        var rnd = new Random(seed);
        var chars = new char[codons * 3];
        for (int i = 0; i < codons; i++)
        {
            // GCA GCC GCG GCT .. any sense codon: avoid stops so the sequence stays a plain open reading frame.
            string codon;
            do
            {
                codon = new string(new[] { "ACGT"[rnd.Next(4)], "ACGT"[rnd.Next(4)], "ACGT"[rnd.Next(4)] });
            } while (codon is "TAA" or "TAG" or "TGA");
            codon.CopyTo(0, chars, 3 * i, 3);
        }

        return new string(chars);
    }

    [Test]
    [Description("R14: 48 000 nt with the site N took 8.2 s through this tool")]
    public void RemoveRestrictionSites_48000nt_SiteN_FinishesQuickly()
    {
        string dna = RandomCodingDna(16_000, 1);

        var sw = Stopwatch.StartNew();
        var result = MolToolsTools.remove_restriction_sites(dna, ["N"], EColi());
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(3)), "elapsed");
            Assert.That(result.OptimizedSequence, Is.EqualTo(dna.Replace('T', 'U')), "N matches every position: nothing can be removed");
        });
    }

    [Test]
    public void RemoveRestrictionSites_AnEcoRiSiteInALongSequence_IsRemoved_AndTheProteinIsKept()
    {
        string dna = RandomCodingDna(16_000, 2).Remove(300, 6).Insert(300, "GAATTC");

        var sw = Stopwatch.StartNew();
        var result = MolToolsTools.remove_restriction_sites(dna, ["GAATTC"], EColi());
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(3)), "elapsed");
            Assert.That(result.OptimizedSequence, Has.Length.EqualTo(dna.Length));
            Assert.That(result.OptimizedSequence, Does.Not.Contain("GAAUUC"));
        });
    }

    // O9-2: the check no longer depends on the data.
    [TestCase("AUG")]
    [TestCase("AUGGCUAAAGAAUUCGGAUCC")]
    public void RemoveRestrictionSites_ASiteOver64Characters_IsRefused_ForAnySequence(string sequence)
    {
        var e = Assert.Throws<ArgumentException>(() => MolToolsTools.remove_restriction_sites(sequence, [TooLong], EColi()));

        Assert.That(e!.Message, Does.Contain("64").And.Contain("65"));
    }

    [TestCase("AU")]
    [TestCase("AUG")]
    [TestCase("AUGGCUAAAGAAUUCGGAUCC")]
    public void RemoveRestrictionSites_ANonIupacSite_IsRefused_ForAnySequence(string sequence)
    {
        var e = Assert.Throws<ArgumentException>(() => MolToolsTools.remove_restriction_sites(sequence, ["GAAXTC"], EColi()));

        Assert.That(e!.Message, Does.Contain("GAAXTC"));
    }

    [TestCase("AU")]
    [TestCase("AUG")]
    public void RemoveRestrictionSites_TheLongestPermittedSite_IsAccepted_ForAnySequence(string sequence)
    {
        var result = MolToolsTools.remove_restriction_sites(sequence, [new string('N', 64), "GAAUUC", ""], EColi());

        Assert.That(result.OptimizedSequence, Is.EqualTo(sequence == "AU" ? "" : "AUG"));
    }

    [TestCase("AU")]
    [TestCase("AUG")]
    [TestCase("AUGGCUAAAGAAUUCGGAUCC")]
    public void OptimizeCodons_AnInvalidAvoidSite_IsRefused_ForAnySequence(string sequence)
    {
        var tooLong = Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(sequence, EColi(), avoid_sites: [TooLong]));
        var nonIupac = Assert.Throws<ArgumentException>(() => MolToolsTools.optimize_codons(sequence, EColi(), avoid_sites: ["GAAXTC"]));

        Assert.Multiple(() =>
        {
            Assert.That(tooLong!.Message, Does.Contain("64"));
            Assert.That(tooLong.ParamName, Is.EqualTo("avoidSites"));
            Assert.That(nonIupac!.ParamName, Is.EqualTo("avoidSites"));
        });
    }

    // The limit has to be readable where a caller (or a model) picks the arguments.
    [Test]
    public void TheDescriptions_StateTheSiteLengthLimit()
    {
        var tools = typeof(MolToolsTools);
        var remove = tools.GetMethod(nameof(MolToolsTools.remove_restriction_sites))!;
        var optimize = tools.GetMethod(nameof(MolToolsTools.optimize_codons))!;

        string removeTool = remove.GetCustomAttribute<System.ComponentModel.DescriptionAttribute>()!.Description;
        string removeSites = remove.GetParameters().Single(p => p.Name == "restriction_sites").GetCustomAttribute<System.ComponentModel.DescriptionAttribute>()!.Description;
        string avoidSites = optimize.GetParameters().Single(p => p.Name == "avoid_sites").GetCustomAttribute<System.ComponentModel.DescriptionAttribute>()!.Description;

        Assert.Multiple(() =>
        {
            Assert.That(removeTool, Does.Contain("at most 64 characters"));
            Assert.That(removeSites, Does.Contain("at most 64 characters"));
            Assert.That(avoidSites, Does.Contain("at most 64 characters"));
        });
    }
}
