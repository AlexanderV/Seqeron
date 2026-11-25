using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// MIRNA-PCT-001 (WP-3, F35/F36): TargetScan 7.0 conservation from a multi-species 3'UTR alignment —
/// site finding and cross-species grouping (<c>targetscan_70.pl</c>), the UTR bin
/// (<c>targetscan_70_BL_bins.pl</c>) and the end-to-end branch length / PCT (<c>targetscan_70_BL_PCT.pl</c>).
/// Oracles: the official kit outputs (perl 5.38.2 + BioPerl + Statistics::Lite re-run reproduces
/// <c>targetscan_70_output.txt</c> and <c>UTRs_median_BLs_bins.txt</c> byte for byte) and constructed
/// inputs whose expected rows are what the perl printed.
/// </summary>
[TestFixture]
[Category("MiRNA")]
public class MiRnaAnalyzer_TargetScanAlignment_Tests
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // Gene → (species → aligned UTR), in file order.
    private static List<(string Gene, Dictionary<string, string> Utrs)> KitGenes()
    {
        var genes = new List<(string, Dictionary<string, string>)>();
        foreach (var line in MiRnaTargetScanAlignmentKitData.Utrs)
        {
            var f = line.Split('\t');
            if (genes.Count == 0 || genes[^1].Item1 != f[0])
                genes.Add((f[0], new Dictionary<string, string>(StringComparer.Ordinal)));
            genes[^1].Item2[f[1]] = f[2];
        }
        return genes;
    }

    // Family → (seed, species with the family), iterated in perl `sort keys` (ordinal) order.
    private static SortedDictionary<string, (string Seed, string[] Species)> KitFamilies()
    {
        var fams = new SortedDictionary<string, (string, string[])>(StringComparer.Ordinal);
        foreach (var line in MiRnaTargetScanAlignmentKitData.MiRnaFamilies)
        {
            var f = line.Split('\t');
            fams[f[0]] = (f[1], f[2].Split(';'));
        }
        return fams;
    }

    private static string TypeName(TargetSiteType t) => t switch
    {
        TargetSiteType.Seed8mer => "8mer-1a",
        TargetSiteType.Seed7merM8 => "7mer-m8",
        TargetSiteType.Seed7merA1 => "7mer-1a",
        _ => "6mer"
    };

    // targetscan_70.pl prints one row per site; Group_num continues over genes and families.
    [Test]
    public void KitSample_AllSiteRows_EqualTargetScan70Output()
    {
        var produced = new List<string>();
        int offset = 0;
        foreach (var (gene, utrs) in KitGenes())
        {
            foreach (var (family, (seed, species)) in KitFamilies())
            {
                var groups = FindTargetScanSiteGroups(seed, utrs, species);
                for (int gi = 0; gi < groups.Count; gi++)
                {
                    var g = groups[gi];
                    Assert.That(g.GroupNumber, Is.EqualTo(gi + 1));
                    foreach (var s in g.Sites)
                    {
                        string type = TypeName(s.SiteType);
                        string perType = g.GroupType == type ? "" : string.Join(" ", g.GetSpeciesWithSiteType(s.SiteType));
                        produced.Add(string.Join("\t", gene, family, s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd,
                            offset + g.GroupNumber, type, "x", g.GroupType, string.Join(" ", g.Species), perType, s.OrfOverlap ? 1 : 0));
                    }
                }
                offset += groups.Count;
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(produced, Has.Count.EqualTo(213));
            Assert.That(produced, Is.EqualTo(MiRnaTargetScanAlignmentKitData.SiteRows));
            Assert.That(offset, Is.EqualTo(121), "last Group_num of the kit output");
        });
    }

    // targetscan_70_BL_bins.pl prints "$gene\t$median\t$bin" (perl stringifies the median with %.15g).
    [Test]
    public void KitSample_UtrBins_EqualBinsScript()
    {
        var produced = KitGenes().Select(g =>
        {
            var bin = AssignTargetScanUtrBin(g.Utrs);
            return $"{g.Gene}\t{bin.MedianBranchLength.ToString("G15", Inv)}\t{bin.Bin}";
        }).ToList();

        Assert.That(produced, Is.EqualTo(MiRnaTargetScanAlignmentKitData.Bins));
    }

    // From the raw alignment to targetscan_70_output.BL_PCT.txt: every row's Branch length score, Pct and
    // Conserved columns (rows are matched by gene, family, species and MSA columns because the script
    // re-sorts its input with `sort -k8,8n`).
    [Test]
    public void KitSample_EndToEnd_BranchLengthPctConserved_EqualTargetScan()
    {
        var expected = new Dictionary<string, string[]>(StringComparer.Ordinal);
        foreach (var line in MiRnaTargetScanAlignmentKitData.BlPctRows)
        {
            var f = line.Split('\t');
            expected.Add(string.Join("|", f[0], f[1], f[2], f[3], f[4]), f);
        }

        var mismatches = new List<string>();
        int rows = 0, pctRows = 0, nonZeroPct = 0, conserved = 0, offset = 0;
        foreach (var (gene, utrs) in KitGenes())
        {
            foreach (var (family, (seed, species)) in KitFamilies())
            {
                var result = ScoreTargetScanConservation(seed, utrs, species);
                foreach (var s in result.Sites)
                {
                    rows++;
                    string key = string.Join("|", gene, family, s.Site.Species, s.Site.MsaStart, s.Site.MsaEnd);
                    if (!expected.TryGetValue(key, out var f)) { mismatches.Add("unexpected " + key); continue; }
                    if (f[5] != s.Site.UtrStart.ToString(Inv) || f[6] != s.Site.UtrEnd.ToString(Inv)) mismatches.Add("UTR coords " + key);
                    if (f[7] != (offset + s.GroupNumber).ToString(Inv)) mismatches.Add("group " + key);
                    if (f[8] != TypeName(s.Site.SiteType) || f[10] != f[8]) mismatches.Add("type " + key);
                    if (double.Parse(f[11], Inv) != s.BranchLength) mismatches.Add($"BL {key}: {s.BranchLength} vs {f[11]}");
                    if (f[12] == "NA") { if (s.Pct is not null) mismatches.Add("PCT not NA " + key); }
                    else
                    {
                        pctRows++;
                        if (s.Pct is not double p || p != double.Parse(f[12], Inv)) mismatches.Add($"PCT {key}: {s.Pct} vs {f[12]}");
                        if (double.Parse(f[12], Inv) != 0.0) nonZeroPct++;
                    }
                    string cons = f.Length > 13 ? f[13] : "";
                    if ((cons == "x") != s.Conserved) mismatches.Add("Conserved " + key);
                    if (s.Conserved) conserved++;
                }
                offset += result.Groups.Count;
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(mismatches, Is.Empty);
            Assert.That(rows, Is.EqualTo(213));
            Assert.That(pctRows, Is.EqualTo(161));
            Assert.That(nonZeroPct, Is.EqualTo(28));
            Assert.That(conserved, Is.EqualTo(MiRnaTargetScanAlignmentKitData.BlPctRows.Count(l => l.Split('\t') is var f && f.Length > 13 && f[13] == "x")));
        });
    }

    // The bins script's per-gene bins (UTRs_median_BLs_bins.txt) feed the BL step; the computed bins equal them.
    [Test]
    public void KitSample_ComputedBins_SelectTheKitTrees()
    {
        var bins = MiRnaTargetScanAlignmentKitData.Bins.Select(l => l.Split('\t')).ToDictionary(f => f[0], f => int.Parse(f[2], Inv));
        foreach (var (gene, utrs) in KitGenes())
            Assert.That(ScoreTargetScanConservation("GAGGUAG", utrs, utrBin: null).UtrBin.Bin, Is.EqualTo(bins[gene]), gene);
    }

    [Test]
    public void BundledTrees_EqualKitFiles_SingleLine_84Leaves()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetTargetScanBinTreeNewick(1), Is.EqualTo(MiRnaPctTargetScanKitData.TreeBin01));
            Assert.That(GetTargetScanBinTreeNewick(2), Is.EqualTo(MiRnaPctTargetScanKitData.TreeBin02));
            for (int b = 1; b <= 10; b++)
            {
                string n = GetTargetScanBinTreeNewick(b);
                Assert.That(n.Any(char.IsWhiteSpace), Is.False, $"bin {b} single line");
                Assert.That(CountLeaves(GetTargetScanBinTree(b)), Is.EqualTo(84), $"bin {b}");
            }
            Assert.That(CountLeaves(GetTargetScanGenericTree()), Is.EqualTo(84));
            Assert.That(TargetScanGenericTreeNewick, Does.StartWith("((((((((((((((((9606:0.00499585,9598:0.00495799):0.00169118,"), "Tree.generic.txt first lines");
            Assert.That(TargetScanUtrBinThresholds, Is.EqualTo(new[] { 0, 1.21207417, 2.17396073, 2.80215414, 3.26272822, 3.65499277, 4.01461968, 4.40729032, 4.90457274, 5.78196252 }));
            Assert.That(() => GetTargetScanBinTreeNewick(0), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => GetTargetScanBinTreeNewick(11), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // Rows as targetscan_70.pl prints them (gene, family, "x" and the global group offset fixed).
    private static List<string> PerlRows(string gene, string family, IReadOnlyList<TargetScanSiteGroup> groups) =>
        groups.SelectMany(g => g.Sites.Select(s =>
        {
            string type = TypeName(s.SiteType);
            string perType = g.GroupType == type ? "" : string.Join(" ", g.GetSpeciesWithSiteType(s.SiteType));
            return string.Join("\t", gene, family, s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd,
                g.GroupNumber, type, "x", g.GroupType, string.Join(" ", g.Species), perType, s.OrfOverlap ? 1 : 0);
        })).ToList();

    // Constructed alignment run through the official targetscan_70.pl (perl 5.38.2), output verbatim.
    // 9606: lower-case/t/N/gap inside the 14-letter mask; a gapped 8mer (cols 23–33, its 7mers and 6mer
    // dropped); 8mers "cUACCUCA" (lower-case first letter ⇒ ORF_overlap 0) and "CUACCUCAg" (lower-case
    // letter right after the site ⇒ 1): the script tests the window one column right of the site; a final
    // "tUACCUCA" read as UUACCUCA (7mer-A1). Same-start/end and ≥ 2-column overlaps group with 10090;
    // the two unpaired 9606 sites get the last group numbers.
    [Test]
    public void ConstructedAlignment_MaskGapsSubsetsOrfWindow_EqualPerl()
    {
        var utrs = new Dictionary<string, string>
        {
            ["9606"] = "acgtNACGTACGT-ACGUACGGCU-ACC--UCAGGUACCUCGGGCUACCUCUGGAAUACCUCAAGGcUACCUCAGGGCUACCUCAgGGtUACCUCA",
            ["10090"] = "ACGUACGUACGUAC--GUACGGCUACCUCA---GGUACCUCGGG-CUACCUCGGAAUACCUCAAGGGGGGGGGGGGGGGGGGGGGGGGGUACCUCA",
        };
        string[] perl =
        {
            "GENE1\tlet-7/98\t10090\t23\t30\t21\t28\t1\t8mer-1a\tx\t8mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t23\t33\t22\t29\t1\t8mer-1a\tx\t8mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t10090\t36\t41\t31\t36\t2\t6mer\tx\t6mer\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t36\t41\t32\t37\t2\t6mer\tx\t6mer\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t10090\t46\t52\t40\t46\t3\t7mer-m8\tx\t7mer-m8\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t45\t51\t41\t47\t3\t7mer-m8\tx\t7mer-m8\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t10090\t57\t63\t51\t57\t4\t7mer-1a\tx\t7mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t57\t63\t53\t59\t4\t7mer-1a\tx\t7mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t10090\t90\t96\t84\t90\t5\t7mer-1a\tx\t7mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t90\t96\t86\t92\t5\t7mer-1a\tx\t7mer-1a\t9606 10090\t\t0",
            "GENE1\tlet-7/98\t9606\t67\t74\t63\t70\t6\t8mer-1a\tx\t8mer-1a\t9606\t\t0",
            "GENE1\tlet-7/98\t9606\t78\t85\t74\t81\t7\t8mer-1a\tx\t8mer-1a\t9606\t\t1",
        };
        Assert.That(PerlRows("GENE1", "let-7/98", FindTargetScanSiteGroups("GAGGUAG", utrs)), Is.EqualTo(perl));
        // The mature miRNA (nt 2–8 = GAGGUAG) and the DNA seed give the same groups.
        Assert.That(PerlRows("GENE1", "let-7/98", FindTargetScanSiteGroups("UGAGGUAGUAGGUUGUAUAGUU", utrs)), Is.EqualTo(perl));
        Assert.That(PerlRows("GENE1", "let-7/98", FindTargetScanSiteGroups("gaggtag", utrs)), Is.EqualTo(perl));
    }

    // groupThisPair never merges two existing groups: species 5's gapped 8mer (cols 35–56) overlaps the
    // sites of species 1/2 (group 1) and 3/4 (group 2) but only joins group 1, as the perl printed.
    [Test]
    public void ConstructedAlignment_GroupsAreNotMerged_EqualPerl()
    {
        string Row(int col, string site) => new string('G', col - 1) + site + new string('G', 71 - col - site.Length);
        var utrs = new Dictionary<string, string>
        {
            ["1"] = Row(30, "CUACCUCA"), ["2"] = Row(32, "CUACCUCA"), ["3"] = Row(50, "CUACCUCA"),
            ["4"] = Row(51, "CUACCUCA"), ["5"] = Row(35, "C--U--A--C--C--U--C--A"),
        };
        string[] perl =
        {
            "GENE2\tlet-7/98\t1\t30\t37\t30\t37\t1\t8mer-1a\tx\t8mer-1a\t1 2 5\t\t0",
            "GENE2\tlet-7/98\t2\t32\t39\t32\t39\t1\t8mer-1a\tx\t8mer-1a\t1 2 5\t\t0",
            "GENE2\tlet-7/98\t5\t35\t56\t35\t42\t1\t8mer-1a\tx\t8mer-1a\t1 2 5\t\t0",
            "GENE2\tlet-7/98\t3\t50\t57\t50\t57\t2\t8mer-1a\tx\t8mer-1a\t3 4\t\t0",
            "GENE2\tlet-7/98\t4\t51\t58\t51\t58\t2\t8mer-1a\tx\t8mer-1a\t3 4\t\t0",
        };
        Assert.That(PerlRows("GENE2", "let-7/98", FindTargetScanSiteGroups("GAGGUAG", utrs, new[] { "1", "2", "3", "4", "5" })), Is.EqualTo(perl));
    }

    [Test]
    public void SpeciesWithMiRna_RestrictsTheSearch()
    {
        var utrs = new Dictionary<string, string>
        {
            ["9606"] = new string('G', 20) + "CUACCUCA" + "GG",
            ["10090"] = new string('G', 20) + "CUACCUCA" + "GG",
        };
        var all = FindTargetScanSiteGroups("GAGGUAG", utrs);
        var humanOnly = FindTargetScanSiteGroups("GAGGUAG", utrs, new[] { "9606" });
        Assert.Multiple(() =>
        {
            Assert.That(all, Has.Count.EqualTo(1));
            Assert.That(all[0].Species, Is.EqualTo(new[] { "9606", "10090" }), "numeric order");
            Assert.That(humanOnly.Single().Species, Is.EqualTo(new[] { "9606" }));
            Assert.That(FindTargetScanSiteGroups("GAGGUAG", utrs, Array.Empty<string>()), Is.Empty);
        });
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        var utrs = new Dictionary<string, string> { ["9606"] = "ACGU" };
        Assert.Multiple(() =>
        {
            Assert.That(() => FindTargetScanSiteGroups(null!, utrs), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => FindTargetScanSiteGroups("GAGGUA", utrs), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => FindTargetScanSiteGroups("GAGNUAG", utrs), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => FindTargetScanSiteGroups("GAGGUAG", null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => FindTargetScanSiteGroups("GAGGUAG", new Dictionary<string, string> { [""] = "A" }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => AssignTargetScanUtrBin(new Dictionary<string, string> { ["10090"] = "ACGU" }), NUnit.Framework.Throws.ArgumentException, "no reference row");
            Assert.That(() => AssignTargetScanUtrBin(new Dictionary<string, string> { ["X"] = "ACGU" }, "X"), NUnit.Framework.Throws.ArgumentException, "reference not in tree");
            Assert.That(() => ScoreTargetScanConservation("GAGGUAG", utrs, utrBin: 11), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // targetscan_70_BL_bins.pl on a constructed 11-row alignment (lower-case row, a short row, N and
    // gaps, species 99999 not in the tree): "GENE3	2.55195311	3".
    [Test]
    public void ConstructedAlignment_UtrBin_EqualBinsScript()
    {
        var utrs = new Dictionary<string, string>
        {
            ["9606"] = "ACGUACGUAC-GUACGUACGU", ["9598"] = "NCG-ACNUA--GNAC-UANGU", ["9544"] = "NCGUNCGUNC-GNACGNACGN",
            ["10090"] = "NCGUA-GUAC-GUAC-UACGN", ["10116"] = "ncguacnuac-gnacguangu", ["9913"] = "NCGUACG-AC-GUANGUACGU",
            ["9615"] = "NCGUACGUNC-GUACGNACGU", ["9031"] = "NCGUACGUA--GUACGUANGU", ["8364"] = "NCGUACGUAC-GUACGUACGN",
            ["7897"] = "NCGUACGUAC--", ["99999"] = "NCGUACGUAC-GNACGUACGU",
        };
        var bin = AssignTargetScanUtrBin(utrs);
        Assert.Multiple(() =>
        {
            Assert.That(bin.MedianBranchLength.ToString("G15", Inv), Is.EqualTo("2.55195311"));
            Assert.That(bin.Bin, Is.EqualTo(3));
            Assert.That(AssignTargetScanUtrBin(new Dictionary<string, string> { ["9606"] = "---" }).Bin, Is.EqualTo(1), "no reference nucleotide");
            Assert.That(double.IsNaN(AssignTargetScanUtrBin(new Dictionary<string, string> { ["9606"] = "---" }).MedianBranchLength), Is.True);
        });
    }

    // A group's per-type species list plugs into context++ as the PCT conservation evidence.
    [Test]
    public void SiteGroup_FeedsContextPlusPlusPct_SameAsEndToEnd()
    {
        var (_, utrs) = KitGenes().Single(g => g.Gene == "FNDC3A");
        string[] species = KitFamilies()["miR-1/206"].Species;
        var result = ScoreTargetScanConservation("GGAAUGU", utrs, species);
        var site = result.Sites.First(s => s.Site.Species == "9606" && s.Site.SiteType == TargetSiteType.Seed8mer);
        var group = result.Groups.Single(g => g.GroupNumber == site.GroupNumber);
        double bl = ComputeBranchLengthScore(GetTargetScanBinTree(result.UtrBin.Bin), group.GetSpeciesWithSiteType(TargetSiteType.Seed8mer));
        Assert.Multiple(() =>
        {
            Assert.That(result.UtrBin.Bin, Is.EqualTo(2));
            Assert.That(group.GroupType, Is.EqualTo("7mer-m8+8mer-1a"));
            Assert.That(group.GetSpeciesWithSiteType(TargetSiteType.Seed8mer), Is.EqualTo(new[] { "9606", "9615", "10090", "10116" }));
            Assert.That(group.GetSpeciesWithSiteType(TargetSiteType.Seed7merM8), Is.EqualTo(new[] { "9031", "9606", "9615", "10090", "10116" }));
            Assert.That(bl, Is.EqualTo(site.BranchLength));
        });
    }

    private static int CountLeaves(Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.PhyloNode n) =>
        n.IsLeaf ? 1 : n.Children.Sum(CountLeaves);
}
