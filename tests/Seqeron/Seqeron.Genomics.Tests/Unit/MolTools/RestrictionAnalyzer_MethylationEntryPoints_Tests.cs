// Canonical Test File: RESTR-DIGEST-001 / RESTR-FIND-001 — the opt-in methylation mode on the entry points
// F68 left out (B08 round 3 G3 = F83: PartialDigest; G4 = F84: FindAllSites(EnzymeSet), raw-string FindSites,
// FindSites over an enzyme list).
// Methods Under Test: RestrictionAnalyzer.PartialDigest(..., methylases, names|enzymes),
//                     FindAllSites(seq, topology, set, methylases), FindSites(string, name, topology, methylases),
//                     FindSites(seq, topology, enzymes, methylases).
//
// Reference: EMBOSS `restrict -methylation -mfile Emethylsites.dat` (Ubuntu emboss 6.6.0+dfsg-12ubuntu2, executed;
// enzyme table generated from Biopython 1.88 `Bio.Restriction`; recipe in review-2026-09/B08.md F68 and F83/F84).
//
// G3 oracle (partial digest = the constructed set of F53, Rosalind PDPL): the surviving cut set is built
// independently of the library — Biopython 1.88 `Enzyme.search(Seq, linear=...)` gives each enzyme's cuts, a cut is
// kept when EMBOSS (run with and without -methylation) still reports the enzyme's site at `cut - fst5`; the
// partial-digest fragments are then every [cut_i, cut_j) interval over {0, cuts, L} (linear), or for k cuts every
// cut-to-cut arc k-1 .. 0 missed plus the intact circle (circular, k^2 + 1). The tests below build that set from the
// cut lists (constants from the Python run) and compare it with the library's fragments, member for member
// (start, length, missed cuts).
//
// TestSpec: TestSpecs/RESTR-FIND-001.md, TestSpecs/RESTR-DIGEST-001.md

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class RestrictionAnalyzer_MethylationEntryPoints_Tests
{
    private const string L1 = "AAAAGGATCCAAAATGATCAAAAAGATCAAAAACCAGGAAAACCTGGAAAAATTT";
    private const string L2 = "GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA";
    private const string L3 = "TTTTAACGGTACCGTGCTTTTGAATTCTTTTGCACGGTACCGTTTTTGGTACCTTT";
    // 500-bp random sequence with injected sites (seeded harness sequence #0), no methylase site spans the origin.
    private const string R0 = "CGAGCACGGTACCCGGGCTTCCGACCTTAGTACCGAGCTTAGATCCGTCATTGTTCGGTGCATGTAAGGAAACTAAATTGTGCATGGCGCACTTTCTATCTGACCCAGCAGGTACCGCACACACTACGGCCTCGAGGACGGCATTGAACTGGGCGGATGGTTCCTTCGAGATCTAGAGACACTTAAGCCCAGTCCCACCTGCATACGTTTAGCGCTAAGCTCAGCTAGCGAAATTCGTGTTGATTTAGGGTGAGCGATTTCGCAGCCCAAACGGTGTGGTGCCTTCCTGGGCTAACATGCACAGAGGGTTGCCTGTCGATGGTACGACTCGGTCTAATAAGTCGACCAACACTCTCGAGTCGGGAGTAAATTATGGGAGATCTTGCTGTTTGAACTTACCATTCAACCCCGGGTTGGGATGGGTACCAGAAAATAACTCGACTGTGCAGGGTATGTTCAAAATGCCGCAGACTACTGCAGATCCACTGTCGCGTCGTGAC";

    private static readonly string[] R0Enzymes = ["MboI", "BamHI", "BclI", "EcoRII", "KpnI", "EcoRI", "TaqI", "SmaI"];

    #region Constructed-set oracle (F53 method on the EMBOSS-surviving cuts)

    private static (int Start, int Length, int Missed)[] Constructed(int[] cuts, int length, bool circular)
    {
        var result = new List<(int, int, int)>();
        if (!circular)
        {
            var x = new List<int> { 0 };
            x.AddRange(cuts);
            x.Add(length);
            for (int i = 0; i < x.Count; i++)
                for (int j = i + 1; j < x.Count; j++)
                    result.Add((x[i], x[j] - x[i], j - i - 1));
        }
        else if (cuts.Length == 0)
        {
            result.Add((0, length, 0));
        }
        else
        {
            int k = cuts.Length;
            for (int i = 0; i < k; i++)
                for (int step = 1; step <= k; step++)
                {
                    int arc = ((cuts[(i + step) % k] - cuts[i]) % length + length) % length;
                    result.Add((cuts[i], arc == 0 ? length : arc, step - 1));
                }
            result.Add((0, length, k));
        }

        return [.. result.OrderBy(t => t)];
    }

    private static List<RestrictionEnzyme> Resolve(IEnumerable<string> names) =>
        [.. names.Select(n => RestrictionAnalyzer.GetEnzyme(n)!)];

    private static (int Start, int Length, int Missed)[] Actual(IEnumerable<PartialDigestFragment> fragments) =>
        [.. fragments.Select(f => (f.Fragment.StartPosition, f.Fragment.Length, f.MissedCuts)).OrderBy(t => t)];

    // (sequence, enzymes, methylation set, cuts surviving EMBOSS — same for linear and circular unless a second list is given)
    private static IEnumerable<TestCaseData> Surviving()
    {
        string[] l2 = ["EcoRI", "MboI", "SmaI", "TaqI"];
        // L2: Dam blocks MboI (10, 36); EcoRI 0/18, SmaI 27 (cut 30), TaqI 38/44 (cuts 39/45) stay.
        yield return new TestCaseData(L2, l2, "none", new[] { 1, 10, 19, 30, 36, 39, 45 }).SetName("L2_none");
        yield return new TestCaseData(L2, l2, "dam", new[] { 1, 19, 30, 39, 45 }).SetName("L2_dam");
        yield return new TestCaseData(L2, l2, "emboss", new[] { 1, 19, 30, 39, 45 }).SetName("L2_emboss");
        // L1: MboI/Sau3AI/BamHI/BclI (GATC windows) all blocked by Dam -> no cut at all; EcoRII needs Dcm.
        string[] l1 = ["MboI", "Sau3AI", "BamHI", "BclI"];
        yield return new TestCaseData(L1, l1, "none", new[] { 5, 15, 24 }).SetName("L1_none");
        yield return new TestCaseData(L1, l1, "dam", Array.Empty<int>()).SetName("L1_dam");
        yield return new TestCaseData(L1, l1.Append("EcoRII").ToArray(), "dam", new[] { 33, 42 }).SetName("L1_dam_EcoRII");
        yield return new TestCaseData(L1, l1.Append("EcoRII").ToArray(), "emboss", Array.Empty<int>()).SetName("L1_emboss_EcoRII");
        // L3: KpnI/Acc65I 7, 35, 47 -> only 47 survives; EcoRI 21 stays.
        string[] l3 = ["KpnI", "Acc65I", "EcoRI"];
        yield return new TestCaseData(L3, l3, "none", new[] { 8, 12, 22, 36, 40, 48, 52 }).SetName("L3_none");
        yield return new TestCaseData(L3, l3, "emboss", new[] { 22, 48, 52 }).SetName("L3_emboss");
        // R0 (500 bp, 8 enzymes): EMBOSS keeps 10 of 16 cuts with all five methylase entries, 12 with Dam only, 15 with Dcm only.
        yield return new TestCaseData(R0, R0Enzymes, "none",
            new[] { 12, 14, 41, 115, 132, 166, 169, 285, 316, 342, 355, 378, 410, 426, 438, 479 }).SetName("R0_none");
        yield return new TestCaseData(R0, R0Enzymes, "emboss",
            new[] { 12, 14, 115, 132, 166, 316, 342, 355, 410, 426 }).SetName("R0_emboss");
        yield return new TestCaseData(R0, R0Enzymes, "dam",
            new[] { 12, 14, 115, 132, 166, 285, 316, 342, 355, 410, 426, 438 }).SetName("R0_dam");
        yield return new TestCaseData(R0, R0Enzymes, "dcm",
            new[] { 12, 14, 41, 115, 132, 166, 169, 316, 342, 355, 378, 410, 426, 438, 479 }).SetName("R0_dcm");
    }

    [TestCaseSource(nameof(Surviving))]
    public void PartialDigest_LinearAndCircular_EqualTheConstructedSetOfTheEmbossSurvivingCuts(
        string seq, string[] enzymes, string methylation, int[] cuts)
    {
        var methylases = MethylaseSets.FromName(methylation);
        var dna = new DnaSequence(seq);

        Assert.Multiple(() =>
        {
            foreach (bool circular in new[] { false, true })
            {
                var topology = circular ? MoleculeTopology.Circular : MoleculeTopology.Linear;
                var byName = RestrictionAnalyzer.PartialDigest(dna, topology, null, methylases, enzymes);
                var byRecord = RestrictionAnalyzer.PartialDigest(
                    dna, topology, null, Resolve(enzymes), methylases);
                var expected = Constructed(cuts, seq.Length, circular);

                Assert.That(Actual(byName), Is.EqualTo(expected), $"{topology} by name");
                Assert.That(Actual(byRecord), Is.EqualTo(expected), $"{topology} by record");
                // count formula of F53: linear (k+1)(k+2)/2, circular k^2 + 1 (1 for k = 0)
                int k = cuts.Length;
                Assert.That(expected, Has.Length.EqualTo(circular ? (k == 0 ? 1 : k * k + 1) : (k + 1) * (k + 2) / 2));
            }
        });
    }

    [TestCaseSource(nameof(Surviving))]
    public void PartialDigest_CompleteDigestPrefix_IsTheMethylationAwareDigest(
        string seq, string[] enzymes, string methylation, int[] _)
    {
        var methylases = MethylaseSets.FromName(methylation);
        var dna = new DnaSequence(seq);
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var complete = RestrictionAnalyzer.Digest(dna, topology, methylases, enzymes).ToList();
            var prefix = RestrictionAnalyzer.PartialDigest(dna, topology, 0, methylases, enzymes)
                .Select(p => p.Fragment).ToList();
            Assert.That(prefix, Is.EqualTo(complete), $"{topology}");
        }
    }

    [Test]
    public void PartialDigest_BlockedSiteIsNotAMissedCut()
    {
        // L2 linear, Dam: 5 cuts remain, so the whole molecule misses 5 (not 7) cuts.
        var all = RestrictionAnalyzer.PartialDigest(
            new DnaSequence(L2), MoleculeTopology.Linear, null, MethylaseSets.Dam, "EcoRI", "MboI", "SmaI", "TaqI").ToList();
        Assert.That(all.Max(f => f.MissedCuts), Is.EqualTo(5));
        Assert.That(all, Has.Count.EqualTo(21));
        // no fragment boundary at the Dam-blocked MboI cuts 10 and 36
        Assert.That(all.Select(f => f.Fragment.StartPosition).Distinct(), Is.EquivalentTo(new[] { 0, 1, 19, 30, 39, 45 }));
    }

    [Test]
    public void PartialDigest_AllCutsBlocked_IsTheIntactMolecule()
    {
        var one = RestrictionAnalyzer.PartialDigest(
            new DnaSequence(L1), MoleculeTopology.Linear, null, MethylaseSets.Dam, "MboI", "Sau3AI", "BamHI", "BclI").Single();
        Assert.That((one.Fragment.StartPosition, one.Fragment.Length, one.MissedCuts), Is.EqualTo((0, 55, 0)));
        Assert.That(one.Fragment.Sequence, Is.EqualTo(L1));
        Assert.That(one.Fragment.LeftEnzyme, Is.Null);
    }

    [Test]
    public void PartialDigest_MaxMissedCuts_BoundsTheMethylatedCutList()
    {
        // k = 10 (R0, EMBOSS emboss list): fragments with <= 2 missed cuts = 11 + 10 + 9.
        var limited = RestrictionAnalyzer.PartialDigest(
            new DnaSequence(R0), MoleculeTopology.Linear, 2, MethylaseSets.Emboss, R0Enzymes).ToList();
        Assert.That(limited, Has.Count.EqualTo(30));
        Assert.That(limited.Max(f => f.MissedCuts), Is.EqualTo(2));
    }

    [Test]
    public void PartialDigest_FragmentLimit_CountsOnlyUnblockedCuts()
    {
        // 16 cuts unmethylated -> 153 fragments; the limit is on the surviving cut list, not the raw sites.
        Assert.That(RestrictionAnalyzer.PartialDigest(
            new DnaSequence(R0), MoleculeTopology.Linear, null, (IReadOnlyList<Methylase>?)null, R0Enzymes).Count(), Is.EqualTo(153));
        Assert.That(RestrictionAnalyzer.PartialDigest(
            new DnaSequence(R0), MoleculeTopology.Linear, null, MethylaseSets.Emboss, R0Enzymes).Count(), Is.EqualTo(66));
    }

    #endregion

    #region Default off (byte-identical)

    [Test]
    public void PartialDigest_NullOrEmptyMethylases_AreTheOverloadWithoutTheParameter()
    {
        Assert.Multiple(() =>
        {
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
            {
                foreach (int? max in new int?[] { null, 0, 2 })
                {
                    var dna = new DnaSequence(R0);
                    var plain = RestrictionAnalyzer.PartialDigest(dna, topology, max, R0Enzymes).ToList();
                    var enzymes = Resolve(R0Enzymes);
                    Assert.That(RestrictionAnalyzer.PartialDigest(dna, topology, max, (IReadOnlyList<Methylase>?)null, R0Enzymes),
                        Is.EqualTo(plain), $"{topology} {max} null");
                    Assert.That(RestrictionAnalyzer.PartialDigest(dna, topology, max, Array.Empty<Methylase>(), R0Enzymes),
                        Is.EqualTo(plain), $"{topology} {max} empty");
                    Assert.That(RestrictionAnalyzer.PartialDigest(dna, topology, max, enzymes, null),
                        Is.EqualTo(RestrictionAnalyzer.PartialDigest(dna, topology, max, enzymes)), $"{topology} {max} record");
                }
            }
        });
    }

    [Test]
    public void PartialDigest_MalformedMethylase_IsRejectedAtTheCall()
    {
        var bad = new[] { new Methylase("x", "GATC") };   // no 'm' marker
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.PartialDigest(
            new DnaSequence(L2), MoleculeTopology.Linear, null, bad, "EcoRI"));
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.PartialDigest(
            new DnaSequence(L2), MoleculeTopology.Linear, null, MethylaseSets.Dam, Array.Empty<string>()));
        // The overload must be named explicitly now that a raw-string sibling exists (B08 H27/F113):
        // a null DnaSequence throws, a null raw string yields no fragments.
        Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.PartialDigest(
            (DnaSequence)null!, MoleculeTopology.Linear, null, MethylaseSets.Dam, "EcoRI"));
        Assert.That(RestrictionAnalyzer.PartialDigest(
                (string?)null, MoleculeTopology.Linear, null, MethylaseSets.Dam, "EcoRI")
            .Select(f => f.Fragment.Length),
            Is.EqualTo(new[] { 0 }),
            "an empty molecule is one zero-length fragment, as new DnaSequence(\"\") also gives");
    }

    #endregion

    #region G4 (F84): FindAllSites(EnzymeSet), raw-string FindSites, FindSites over an enzyme list

    // G4 oracle: EMBOSS `restrict` over the 39 common enzymes (table generated from Biopython 1.88), one sequence per
    // call, with and without -methylation (-mfile: dam.dat = Dam only, dcm.dat = both Dcm lines, Emethylsites.dat = all
    // five), linear and -plasmid. Each row is the set of (enzyme, 0-based window start) EMBOSS reports, both strands
    // collapsed; the library's distinct (enzyme, position) set is the same for both topologies on these sequences (no
    // methylase site spans the origin). Dcm blocks nothing here (no common enzyme has a CCAGG/CCTGG window).
    private static IEnumerable<TestCaseData> FindAllSitesRows()
    {
        yield return new TestCaseData(L1, "none", "BamHI:4;DpnI:5,15,24;MboI:5,15,24;Sau3AI:5,15,24").SetName("L1_none");
        yield return new TestCaseData(L1, "dam", "").SetName("L1_dam");
        yield return new TestCaseData(L1, "dcm", "BamHI:4;DpnI:5,15,24;MboI:5,15,24;Sau3AI:5,15,24").SetName("L1_dcm");
        yield return new TestCaseData(L1, "emboss", "").SetName("L1_emboss");
        yield return new TestCaseData(L2, "none", "DpnI:10,36;EcoRI:0,18;HpaII:28;MboI:10,36;MspI:28;Sau3AI:10,36;SmaI:27;TaqI:38,44").SetName("L2_none");
        yield return new TestCaseData(L2, "dam", "EcoRI:0,18;HpaII:28;MspI:28;SmaI:27;TaqI:38,44").SetName("L2_dam");
        yield return new TestCaseData(L2, "dcm", "DpnI:10,36;EcoRI:0,18;HpaII:28;MboI:10,36;MspI:28;Sau3AI:10,36;SmaI:27;TaqI:38,44").SetName("L2_dcm");
        yield return new TestCaseData(L2, "emboss", "EcoRI:0,18;HpaII:28;MspI:28;SmaI:27;TaqI:38,44").SetName("L2_emboss");
        yield return new TestCaseData(L3, "none", "EcoRI:21;KpnI:7,35,47;RsaI:8,36,48").SetName("L3_none");
        yield return new TestCaseData(L3, "dam", "EcoRI:21;KpnI:7,35,47;RsaI:8,36,48").SetName("L3_dam");
        yield return new TestCaseData(L3, "dcm", "EcoRI:21;KpnI:7,35,47;RsaI:8,36,48").SetName("L3_dcm");
        yield return new TestCaseData(L3, "emboss", "EcoRI:21;KpnI:47;RsaI:48").SetName("L3_emboss");
        yield return new TestCaseData(R0, "none", "AluI:35,217,222;BglII:168,377;DpnI:41,169,378,479;HaeIII:127;HincII:340;HpaII:12,408;KpnI:7,110,421;MboI:41,169,378,479;MspI:12,408;NheI:223;PstI:474;RsaI:8,29,111,321,422;SalI:340;Sau3AI:41,169,378,479;SmaI:11,407;TaqI:131,165,315,341,354,437;XbaI:171;XhoI:130,353").SetName("R0_none");
        yield return new TestCaseData(R0, "dam", "AluI:35,217,222;HaeIII:127;HincII:340;HpaII:12,408;KpnI:7,110,421;MspI:12,408;NheI:223;PstI:474;RsaI:8,29,111,321,422;SalI:340;SmaI:11,407;TaqI:131,165,315,341,354,437;XbaI:171;XhoI:130,353").SetName("R0_dam");
        yield return new TestCaseData(R0, "dcm", "AluI:35,217,222;BglII:168,377;DpnI:41,169,378,479;HaeIII:127;HincII:340;HpaII:12,408;KpnI:7,110,421;MboI:41,169,378,479;MspI:12,408;NheI:223;PstI:474;RsaI:8,29,111,321,422;SalI:340;Sau3AI:41,169,378,479;SmaI:11,407;TaqI:131,165,315,341,354,437;XbaI:171;XhoI:130,353").SetName("R0_dcm");
        yield return new TestCaseData(R0, "emboss", "AluI:35,217,222;HaeIII:127;HincII:340;HpaII:12,408;KpnI:7,110,421;MspI:12,408;NheI:223;PstI:474;RsaI:8,29,111,321,422;SalI:340;SmaI:11,407;TaqI:131,165,315,341,354;XbaI:171;XhoI:130,353").SetName("R0_emboss");
    }

    private static string Format(IEnumerable<RestrictionSite> sites) =>
        string.Join(";", sites.Select(s => (s.Enzyme.Name, s.Position)).Distinct()
            .GroupBy(t => t.Name).OrderBy(g => g.Key, StringComparer.Ordinal)
            .Select(g => g.Key + ":" + string.Join(",", g.Select(t => t.Position).OrderBy(x => x))));

    [TestCaseSource(nameof(FindAllSitesRows))]
    public void FindAllSites_Common_EqualsEmbossSitesUnderTheSet(string seq, string set, string expected)
    {
        var methylases = MethylaseSets.FromName(set);
        var dna = new DnaSequence(seq);
        Assert.Multiple(() =>
        {
            Assert.That(Format(RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Linear, EnzymeSet.Common, methylases)),
                Is.EqualTo(expected), "linear");
            Assert.That(Format(RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Circular, EnzymeSet.Common, methylases)),
                Is.EqualTo(expected), "circular");
        });
    }

    [Test]
    public void FindAllSites_Methylation_DropsExactlyTheBlockedSites()
    {
        // L2 + Dam: EMBOSS removes the Dam windows (MboI 10, 36 ...) and keeps EcoRI 0,18.
        var dna = new DnaSequence(L2);
        var plain = RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Linear, EnzymeSet.Common).ToList();
        var dam = RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Linear, EnzymeSet.Common, MethylaseSets.Dam).ToList();
        Assert.That(dam, Is.SubsetOf(plain));
        Assert.That(dam.Count, Is.LessThan(plain.Count));
        Assert.That(dam.Where(s => s.Enzyme.Name == "EcoRI" && s.IsForwardStrand).Select(s => s.Position), Is.EqualTo(new[] { 0, 18 }));
        Assert.That(dam.Any(s => s.Enzyme.Name == "MboI"), Is.False);
    }

    [Test]
    public void FindAllSites_Rebase_IsThePerEnzymeFindSitesUnderTheSet()
    {
        // One mask for the molecule, applied to every enzyme of the table = the F68 (EMBOSS-confirmed) per-enzyme path.
        var dna = new DnaSequence(R0);
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            foreach (var set in new[] { MethylaseSets.Dam, MethylaseSets.Dcm, MethylaseSets.Emboss })
            {
                var all = RestrictionAnalyzer.FindAllSites(dna, topology, EnzymeSet.Rebase, set).ToList();
                var perEnzyme = RestrictionAnalyzer.FindSites(
                    dna, topology, RestrictionAnalyzer.RebaseEnzymes.Values, set).ToList();
                Assert.That(all, Is.EqualTo(perEnzyme), $"{topology}");
                Assert.That(all.Count, Is.LessThan(RestrictionAnalyzer.FindAllSites(dna, topology, EnzymeSet.Rebase).Count()));
            }
        }
    }

    [Test]
    public void FindAllSites_NullOrEmptyMethylases_AreTheOverloadWithoutTheParameter()
    {
        Assert.Multiple(() =>
        {
            foreach (var set in new[] { EnzymeSet.Common, EnzymeSet.Rebase })
            {
                foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
                {
                    var dna = new DnaSequence(R0);
                    var plain = RestrictionAnalyzer.FindAllSites(dna, topology, set).ToList();
                    Assert.That(RestrictionAnalyzer.FindAllSites(dna, topology, set, null), Is.EqualTo(plain), $"{set} {topology} null");
                    Assert.That(RestrictionAnalyzer.FindAllSites(dna, topology, set, Array.Empty<Methylase>()), Is.EqualTo(plain), $"{set} {topology} empty");
                }
            }

            var d = new DnaSequence(L2);
            Assert.That(RestrictionAnalyzer.FindAllSites(d, MoleculeTopology.Linear, EnzymeSet.Common, null),
                Is.EqualTo(RestrictionAnalyzer.FindAllSites(d)));
            Assert.That(RestrictionAnalyzer.FindAllSites(d, MoleculeTopology.Circular, EnzymeSet.Common, null),
                Is.EqualTo(RestrictionAnalyzer.FindAllSites(d, MoleculeTopology.Circular)));
        });
    }

    [Test]
    public void FindAllSites_Arguments_AreValidatedAtTheCall()
    {
        var dna = new DnaSequence(L2);
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Linear, (EnzymeSet)99, MethylaseSets.Dam));
        Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Linear, EnzymeSet.Common, new[] { new Methylase("x", "GATC") }));
        // Likewise here (B08 H27/F113).
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.FindAllSites((DnaSequence)null!, MoleculeTopology.Linear, EnzymeSet.Common, MethylaseSets.Dam));
        Assert.That(
            RestrictionAnalyzer.FindAllSites((string?)null, MoleculeTopology.Linear, EnzymeSet.Common, MethylaseSets.Dam),
            Is.Empty);
    }

    // Raw-string FindSites: ASCII upper-cased, methylase search on the same string.
    // EMBOSS on GGCCNNNNNGGCCAAGATCAAGAATTCAA: MboI 15 removed by Dam (and by emboss), EcoRI 21 kept. (EMBOSS never reports SfiI
    // here - a sequence N never matches in embPatRestrictScan - while Biopython/the library report it, F40; Dam does not touch it.)
    [Test]
    public void FindSites_RawString_CarriesTheMethylaseSet()
    {
        const string s = "GGCCNNNNNGGCCAAGATCAAGAATTCAA";
        Assert.Multiple(() =>
        {
            foreach (var m in new[] { MethylaseSets.Dam, MethylaseSets.DamDcm, MethylaseSets.Emboss })
            {
                foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
                {
                    Assert.That(RestrictionAnalyzer.FindSites(s, "MboI", topology, m), Is.Empty, $"{topology} MboI");
                    Assert.That(RestrictionAnalyzer.FindSites(s, "EcoRI", topology, m).Where(x => x.IsForwardStrand).Select(x => x.Position),
                        Is.EqualTo(new[] { 21 }), $"{topology} EcoRI");
                    Assert.That(RestrictionAnalyzer.FindSites(s, "SfiI", topology, m).Where(x => x.IsForwardStrand).Select(x => x.Position),
                        Is.EqualTo(new[] { 0 }));
                }
            }

            Assert.That(RestrictionAnalyzer.FindSites(s, "MboI", MoleculeTopology.Linear).Where(x => x.IsForwardStrand).Select(x => x.Position),
                Is.EqualTo(new[] { 15 }));
        });
    }

    [Test]
    public void FindSites_RawString_AgreesWithTheDnaSequenceFormAndFoldsCase()
    {
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            foreach (var name in new[] { "MboI", "EcoRI", "SmaI", "TaqI", "BamHI" })
            {
                var expected = RestrictionAnalyzer.FindSites(new DnaSequence(L2), name, topology, MethylaseSets.Dam).ToList();
                Assert.That(RestrictionAnalyzer.FindSites(L2, name, topology, MethylaseSets.Dam), Is.EqualTo(expected), $"{name} {topology}");
                Assert.That(RestrictionAnalyzer.FindSites(L2.ToLowerInvariant(), name, topology, MethylaseSets.Dam), Is.EqualTo(expected), $"{name} {topology} lower");
            }
        }

        Assert.That(RestrictionAnalyzer.FindSites(L2, "MboI", MoleculeTopology.Linear, MethylaseSets.Dam), Is.Empty);   // EMBOSS: MboI 10, 36 blocked
    }

    [Test]
    public void FindSites_RawString_DefaultOffAndEdges()
    {
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var plain = RestrictionAnalyzer.FindSites(L2, "MboI", topology).ToList();
            Assert.That(plain, Is.Not.Empty);
            Assert.That(RestrictionAnalyzer.FindSites(L2, "MboI", topology, null), Is.EqualTo(plain));
            Assert.That(RestrictionAnalyzer.FindSites(L2, "MboI", topology, Array.Empty<Methylase>()), Is.EqualTo(plain));
        }

        Assert.That(RestrictionAnalyzer.FindSites("", "MboI", MoleculeTopology.Linear, MethylaseSets.Dam), Is.Empty);
        Assert.That(RestrictionAnalyzer.FindSites((string)null!, "MboI", MoleculeTopology.Circular, MethylaseSets.Dam), Is.Empty);
        // an unknown enzyme and a malformed methylase are rejected eagerly, before the sequence is looked at
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites("", "NoSuchI", MoleculeTopology.Linear, MethylaseSets.Dam));
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(L2, "EcoRI", MoleculeTopology.Linear, new[] { new Methylase("x", "GATC") }));
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites("", "EcoRI", MoleculeTopology.Linear, new[] { new Methylase("x", "GATC") }));
        // ASCII-only fold: U+017F is not 'S' (the methylase search runs on the same folded string)
        Assert.That(RestrictionAnalyzer.FindSites("GGATC\u017FC", "BamHI", MoleculeTopology.Linear, MethylaseSets.Dam), Is.Empty);
    }

    [Test]
    public void FindSites_EnzymeList_CarriesTheMethylaseSet()
    {
        var enzymes = new[] { "EcoRI", "MboI", "SmaI", "TaqI" }.Select(n => RestrictionAnalyzer.GetEnzyme(n)!).ToList();
        var custom = new RestrictionEnzyme("DamProbe", "GATC", 0, 4, "test");
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var byName = RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, MethylaseSets.Dam, "EcoRI", "MboI", "SmaI", "TaqI").ToList();
            var byRecord = RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, enzymes, MethylaseSets.Dam).ToList();
            Assert.That(byRecord, Is.EqualTo(byName), $"{topology}");
            Assert.That(byRecord.Any(s => s.Enzyme.Name == "MboI"), Is.False);
            Assert.That(byRecord.Where(s => s.Enzyme.Name == "TaqI" && s.IsForwardStrand).Select(s => s.Position), Is.EqualTo(new[] { 38, 44 }));
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, new[] { custom }, MethylaseSets.Dam), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, new[] { custom }, null), Is.Not.Empty);
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, enzymes, null),
                Is.EqualTo(RestrictionAnalyzer.FindSites(new DnaSequence(L2), topology, enzymes)));
        }

        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(
            new DnaSequence(L2), MoleculeTopology.Linear, Array.Empty<RestrictionEnzyme>(), MethylaseSets.Dam));
    }

    #endregion
}
