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
        Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.PartialDigest(
            null!, MoleculeTopology.Linear, null, MethylaseSets.Dam, "EcoRI"));
    }

    #endregion
}
