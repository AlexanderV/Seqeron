using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// B08 H27/H28: the restriction MCP surface over an IUPAC-bearing sequence, the ranged cut-length
/// listing, and caller-supplied enzyme definitions.
/// </summary>
[TestFixture]
public sealed class RestrictionIupacAndCustomEnzymeTests
{
    /// <summary>
    /// A soft-masked / N-containing molecule: `nnn` + `GAATTC` (EcoRI at 0-based 3) + `nnn` +
    /// `GGATCC` (BamHI at 12) + `NNN`. Every restriction tool used to reject this outright, because it
    /// built a <c>DnaSequence</c>, which accepts only A/C/G/T.
    /// </summary>
    private const string Masked = "nnnGAATTCnnnGGATCCNNN";

    /// <summary>
    /// H27: the defect itself — <c>DnaSequence</c> still rejects the molecule, which is why the
    /// wrappers had to stop constructing one.
    /// </summary>
    [Test]
    public void DnaSequence_StillRejectsTheMaskedMolecule()
    {
        Assert.Throws<ArgumentException>(() => _ = new DnaSequence(Masked));
    }

    /// <summary>
    /// H27: all six sequence-taking restriction tools now answer for it, and the answer is the
    /// library's own raw-string answer — the wrappers delegate, they do not re-derive anything.
    /// </summary>
    [Test]
    public void EveryRestrictionTool_AcceptsAnIupacSequence()
    {
        var sites = MolToolsTools.find_restriction_sites(Masked, ["EcoRI", "BamHI"]).Sites;
        var all = MolToolsTools.find_all_restriction_sites(Masked).Sites;
        var digest = MolToolsTools.restriction_digest(Masked, ["EcoRI"]).Fragments;
        var summary = MolToolsTools.digest_summary(Masked, ["EcoRI"]);
        var partial = MolToolsTools.partial_digest(Masked, ["EcoRI", "BamHI"]).Fragments;
        var map = MolToolsTools.restriction_map(Masked, ["EcoRI", "BamHI"]);

        Assert.Multiple(() =>
        {
            // Two windows, each as a forward and a reverse record (both enzymes are palindromic).
            Assert.That(sites, Has.Count.EqualTo(4));
            Assert.That(sites.Where(s => s.IsForwardStrand).Select(s => s.Position), Is.EquivalentTo(new[] { 3, 12 }));
            Assert.That(sites, Is.EqualTo(RestrictionAnalyzer
                .FindSites(Masked, MoleculeTopology.Linear, null, "EcoRI", "BamHI").ToList()),
                "the wrapper delegates to the raw-string library entry point");

            Assert.That(all.Where(s => s.IsForwardStrand).Select(s => s.Enzyme.Name),
                Does.Contain("EcoRI").And.Contains("BamHI"));

            // EcoRI cuts G^AATTC, so the 0-based top-strand bond is 3 + 1 = 4.
            Assert.That(digest.Select(f => f.Length), Is.EqualTo(new[] { 4, 17 }));
            Assert.That(digest.Sum(f => f.Length), Is.EqualTo(Masked.Length));
            Assert.That(summary.TotalFragments, Is.EqualTo(2));

            // (k+1)(k+2)/2 for k = 2 linear cuts.
            Assert.That(partial, Has.Count.EqualTo(6));
            Assert.That(map.TotalSites, Is.EqualTo(2), "one cleavage per palindromic window");
            Assert.That(map.NonCutters, Is.Empty);
        });
    }

    /// <summary>
    /// H27: the masked bases do not become wildcards. An ambiguity code in the sequence blocks a
    /// concrete site code rather than widening it, so the <c>nnn</c> runs create no extra site —
    /// Biopython's own rule, which <c>find_pattern_occurrences</c> exposes directly.
    /// </summary>
    [Test]
    public void IupacCodesInTheSequence_DoNotCreateExtraSites()
    {
        Assert.Multiple(() =>
        {
            // GAANTC does not match GAATTC: the sequence's N blocks the concrete T.
            Assert.That(MolToolsTools.find_pattern_occurrences("GAANTC", "GAATTC").Positions, Is.Empty);
            // But a pattern N matches a sequence N, exactly as SfiI.search(Seq("GGCCNNNNNGGCC")) does.
            Assert.That(MolToolsTools.find_pattern_occurrences("GGCCNNNNNGGCC", "GGCCNNNNNGGCC").Positions,
                Is.EqualTo(new[] { 0 }));
            // Overlaps are all reported, start_index skips, null yields nothing, bad pattern throws.
            Assert.That(MolToolsTools.find_pattern_occurrences("AAAA", "AA").Positions, Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(MolToolsTools.find_pattern_occurrences("AAAA", "AA", start_index: 2).Positions,
                Is.EqualTo(new[] { 2 }));
            Assert.That(MolToolsTools.find_pattern_occurrences(null, "GAATTC").Positions, Is.Empty);
            Assert.That(MolToolsTools.find_pattern_occurrences("", "GAATTC").Positions, Is.Empty);
            Assert.Throws<ArgumentException>(() => _ = MolToolsTools.find_pattern_occurrences(null, "GAZTC"));
            Assert.Throws<ArgumentException>(() => _ = MolToolsTools.find_pattern_occurrences("ACGT", ""));
            Assert.Throws<ArgumentOutOfRangeException>(
                () => _ = MolToolsTools.find_pattern_occurrences("ACGT", "ACGT", start_index: -1));
        });
    }

    /// <summary>
    /// H28: the ranged <c>GetEnzymesByCutLength(min, max)</c> had no MCP expression, so "all 4-6 bp
    /// cutters" took three calls. The common table has 9 4-cutters and 24 6-cutters and no 5-cutter,
    /// so the range returns 33.
    /// </summary>
    [Test]
    public void EnzymesByCutLength_RangedForm_IsOneCall()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.enzymes_by_cut_length(4).Enzymes, Has.Count.EqualTo(9));
            Assert.That(MolToolsTools.enzymes_by_cut_length(6).Enzymes, Has.Count.EqualTo(24));
            Assert.That(MolToolsTools.enzymes_by_cut_length(5).Enzymes, Is.Empty);
            Assert.That(MolToolsTools.enzymes_by_cut_length(4, max_length: 6).Enzymes, Has.Count.EqualTo(33),
                "9 + 24 + 0 in one call");
            Assert.That(MolToolsTools.enzymes_by_cut_length(4, max_length: 6).Enzymes
                .All(e => e.RecognitionSequence.Length is >= 4 and <= 6), Is.True);
            Assert.That(MolToolsTools.enzymes_by_cut_length(4, "rebase", 6).Enzymes, Has.Count.EqualTo(597));
            Assert.That(MolToolsTools.enzymes_by_cut_length(6, max_length: 4).Enzymes, Is.Empty,
                "an inverted range is empty, as the C# invariant INV-05 states");
            Assert.That(MolToolsTools.enzymes_by_cut_length(4, max_length: 6).Enzymes,
                Is.EqualTo(RestrictionAnalyzer.GetEnzymesByCutLength(4, 6, EnzymeSet.Common).ToList()),
                "the wrapper delegates to the ranged overload");
            Assert.Throws<ArgumentException>(() => _ = MolToolsTools.enzymes_by_cut_length(0));
            Assert.Throws<ArgumentException>(() => _ = MolToolsTools.enzymes_by_cut_length(4, max_length: 0));
        });
    }

    /// <summary>
    /// H28: a caller-supplied enzyme definition had no MCP expression at all. Written out by hand,
    /// EcoRI must give byte-identical records to the named tool.
    /// </summary>
    [Test]
    public void FindSitesCustomEnzyme_EcoRIByHand_EqualsTheNamedTool()
    {
        const string molecule = "AAAGAATTCAAA";
        var custom = MolToolsTools.find_sites_custom_enzyme(molecule, "GAATTC", 1, 5, "MyEcoRI").Sites;
        var named = MolToolsTools.find_restriction_sites(molecule, ["EcoRI"]).Sites;

        Assert.Multiple(() =>
        {
            Assert.That(custom.Select(s => (s.Position, s.IsForwardStrand, s.CutPosition)),
                Is.EqualTo(named.Select(s => (s.Position, s.IsForwardStrand, s.CutPosition))));
            Assert.That(custom.Select(s => s.Enzyme.Name), Is.All.EqualTo("MyEcoRI"));
            Assert.That(custom.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(4),
                "G^AATTC at position 3, so the 0-based top-strand bond is 4");
            // Default name, and a distal Type IIS definition.
            Assert.That(MolToolsTools.find_sites_custom_enzyme(molecule, "GAATTC", 1, 5).Sites[0].Enzyme.Name,
                Is.EqualTo("CustomEnzyme"));
            Assert.That(
                MolToolsTools.find_sites_custom_enzyme("AAAAAAAAAAGGTCTCAAAAAAAAAA", "GGTCTC", 7, 11, "MyBsaI").Sites
                    .Select(s => (s.Position, s.IsForwardStrand, s.CutPosition)),
                Is.EqualTo(new[] { (10, true, 17) }), "GGTCTCN^NNNN cuts 7 past the site start");
            // And an IUPAC-bearing sequence, which is the H27 half of the same story.
            Assert.That(MolToolsTools.find_sites_custom_enzyme("NNNGAATTCNNN", "GAATTC", 1, 5).Sites,
                Has.Count.EqualTo(2));
        });
    }

    /// <summary>
    /// H28: the digest form, and the two-sided (Type IIB) and validation paths.
    /// </summary>
    [Test]
    public void RestrictionDigestCustomEnzyme_MatchesTheNamedToolAndValidates()
    {
        const string molecule = "AAAGAATTCAAA";
        var custom = MolToolsTools.restriction_digest_custom_enzyme(molecule, "GAATTC", 1, 5, "MyEcoRI").Fragments;
        var named = MolToolsTools.restriction_digest(molecule, ["EcoRI"]).Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(custom.Select(f => f.Length), Is.EqualTo(named.Select(f => f.Length)));
            Assert.That(custom.Select(f => f.Length), Is.EqualTo(new[] { 4, 8 }));
            Assert.That(
                MolToolsTools.restriction_digest_custom_enzyme(
                    molecule, "GAATTC", 1, 5, topology: MoleculeTopology.Circular).Fragments.Select(f => f.Length),
                Is.EqualTo(new[] { 12 }), "one cut on a circle linearises it rather than splitting it");

            // A two-sided (Type IIB) definition: BaeI's own geometry, (-10, -15) and (23, 18).
            var baeI = MolToolsTools.find_sites_custom_enzyme(
                new string('A', 30) + "ACTTTTGTACC" + new string('A', 30),
                "ACNNNNGTAYC", -10, -15, "MyBaeI", second_cut_forward: 23, second_cut_reverse: 18).Sites;
            Assert.That(baeI.Count(s => s.IsForwardStrand), Is.EqualTo(2),
                "a two-sided cutter reports one record per cut pair");
            Assert.That(baeI.Count(s => s.IsSecondCut), Is.EqualTo(1));

            // Half a pair, and an inconsistent stagger, are errors.
            Assert.Throws<ArgumentException>(
                () => _ = MolToolsTools.find_sites_custom_enzyme("ACGT", "ACGT", 1, 3, second_cut_forward: 9));
            Assert.Throws<ArgumentException>(
                () => _ = MolToolsTools.find_sites_custom_enzyme("ACGT", "ACGT", 1, 3, second_cut_reverse: 9));
            Assert.Throws<ArgumentException>(
                () => _ = MolToolsTools.restriction_digest_custom_enzyme("ACGT", "ACGT", 1, 3, second_cut_forward: 9),
                "the digest form validates the same way");
            Assert.Throws<ArgumentException>(
                () => _ = MolToolsTools.find_sites_custom_enzyme("ACGT", "ACGT", 1, 3,
                    second_cut_forward: 9, second_cut_reverse: 4).ToString(),
                "the two pairs must leave the same stagger");
            Assert.Throws<ArgumentException>(
                () => _ = MolToolsTools.find_sites_custom_enzyme("ACGT", "GAZTC", 1, 3),
                "a non-IUPAC recognition sequence is rejected");
        });
    }
}
