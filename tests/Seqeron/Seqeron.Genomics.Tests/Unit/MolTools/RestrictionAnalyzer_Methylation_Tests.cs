// Canonical Test File: RESTR-FIND-001 / RESTR-DIGEST-001 — opt-in methylation-blocked recognition
// Methods Under Test: RestrictionAnalyzer.FindSites(DnaSequence, string|RestrictionEnzyme, MoleculeTopology,
//                     IReadOnlyList<Methylase>), Digest(..., methylases, ...), GetDigestSummary(...),
//                     CreateMap(...); MethylaseSets.
//
// Reference: EMBOSS `restrict -methylation -mfile Emethylsites.dat` (Ubuntu emboss 6.6.0+dfsg-12ubuntu2,
// installed locally; the C semantics read in EMBOSS nucleus/embpat.c: patRestrictReadMethyl,
// patRestrictMethylMod, embPatRestrictMatch, embPatRestrictScan). Every expected position below is that
// binary's own output ("Start" - 1, the 0-based start of the recognition window) on the very sequence
// shown, with an enzyme table built from Biopython 1.88 sites/cuts. The mask itself is locked with a
// one-window "NN" probe enzyme: it is reported exactly where the 2-base window holds no methylated base.
// The Emethylsites.dat data is https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/
// Emethylsites.dat (11 lines; fetched 2026-10-21, byte-identical to the Ubuntu package copy).
//
// EMBOSS semantics reproduced (and the two oddities it has, kept for parity and locked here):
//   * every m-marked base AND every N of the site is overwritten by N (EcoKI's 6-base spacer is masked),
//     on both strands; any enzyme window containing such a base is not a site;
//   * entries are visited last-file-entry-first (ajListPush), forward strand then reverse strand, each
//     seeing earlier Ns — for the palindromic GmATC only ONE A of a GATC is marked, so TaqI TCGA in
//     GATCGA is NOT blocked;
//   * this is a window rule, not an enzyme-sensitivity table: BamHI GGATCC is blocked by Dam.
// Deliberate difference (not EMBOSS-lockable): on a circular molecule a methylase site spanning the origin
// is found too (EMBOSS masks the linear string before doubling it).
//
// TestSpec: TestSpecs/RESTR-FIND-001.md, TestSpecs/RESTR-DIGEST-001.md

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class RestrictionAnalyzer_Methylation_Tests
{
    // EMBOSS probe enzyme: reported where a 2-base window has no methylated base.
    private static readonly RestrictionEnzyme Probe = new("PROBE", "NN", 1, 1, "test");

    private const string L1 = "AAAAGGATCCAAAATGATCAAAAAGATCAAAAACCAGGAAAACCTGGAAAAATTT";
    private const string L2 = "GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA";
    private const string L3 = "TTTTAACGGTACCGTGCTTTTGAATTCTTTTGCACGGTACCGTTTTTGGTACCTTT";

    private static int[] Positions(string s, string enzyme, IReadOnlyList<Methylase>? m,
        MoleculeTopology t = MoleculeTopology.Linear) =>
        RestrictionAnalyzer.FindSites(new DnaSequence(s), enzyme, t, m)
            .Where(x => x.IsForwardStrand).Select(x => x.Position).OrderBy(x => x).ToArray();

    #region Data (Emethylsites.dat)

    [Test]
    public void Emboss_IsTheFiveDataLinesOfEmethylsitesDat()
    {
        Assert.That(MethylaseSets.Emboss.Select(m => (m.Name, m.Site)), Is.EqualTo(new[]
        {
            ("dam", "GmATC"), ("dcm", "mCmCAGG"), ("dcm", "mCmCTGG"),
            ("MEcoKI", "mAmACNNNNNNGTGC"), ("MEcoKI", "GCmACNNNNNNGTT"),
        }));
        Assert.That(MethylaseSets.Dam, Has.Count.EqualTo(1));
        Assert.That(MethylaseSets.Dcm, Has.Count.EqualTo(2));
        Assert.That(MethylaseSets.DamDcm, Has.Count.EqualTo(3));
    }

    [TestCase(null, null)]
    [TestCase("", null)]
    [TestCase("none", null)]
    [TestCase("NONE", null)]
    [TestCase("dam", 1)]
    [TestCase("Dcm", 2)]
    [TestCase("dam+dcm", 3)]
    [TestCase("emboss", 5)]
    [TestCase(" Dam+DCM ", 3)]   // G11: ASCII fold + trim
    public void FromName_ResolvesSets(string? name, int? count)
    {
        var set = MethylaseSets.FromName(name);
        if (count is null) Assert.That(set, Is.Null);
        else Assert.That(set!, Has.Count.EqualTo(count));
    }

    [Test]
    public void FromName_UnknownName_Throws() =>
        Assert.Throws<ArgumentException>(() => MethylaseSets.FromName("cpg"));

    [TestCase("dɑm")]      // non-ASCII look-alikes are not folded onto the set names
    [TestCase("emboſs")]   // U+017F would fold to 'S' under ToUpperInvariant
    public void FromName_NonAsciiLookAlike_Throws(string name) =>
        Assert.Throws<ArgumentException>(() => MethylaseSets.FromName(name));

    [TestCase("GATC")]            // no m marker
    [TestCase("GATCm")]           // trailing marker
    [TestCase("GmmATC")]          // doubled marker
    [TestCase("GmAT!")]           // not IUPAC
    [TestCase("")]
    public void MalformedMethylaseSite_Throws(string site) =>
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(
            new DnaSequence("AAGATCAA"), "MboI", MoleculeTopology.Linear, [new Methylase("x", site)]).ToList());

    #endregion

    #region Default off

    [Test]
    public void NullOrEmptySet_IsIdenticalToTheOverloadWithoutIt()
    {
        IReadOnlyList<Methylase>? off = null; // typed: a bare null would also fit params string[]
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var dna = new DnaSequence(L1 + L2 + L3);
            foreach (var enzyme in new[] { "MboI", "BamHI", "BclI", "EcoRII", "EcoRI", "KpnI", "TaqI" })
            {
                var plain = RestrictionAnalyzer.FindSites(dna, enzyme, topology).ToList();
                Assert.That(RestrictionAnalyzer.FindSites(dna, enzyme, topology, null), Is.EqualTo(plain), enzyme);
                Assert.That(RestrictionAnalyzer.FindSites(dna, enzyme, topology, []), Is.EqualTo(plain), enzyme);
            }

            Assert.That(RestrictionAnalyzer.Digest(dna, topology, off, "EcoRI", "MboI"),
                Is.EqualTo(RestrictionAnalyzer.Digest(dna, topology, "EcoRI", "MboI")));
            Assert.That(RestrictionAnalyzer.CreateMap(dna, topology, off, "EcoRI", "MboI").Sites,
                Is.EqualTo(RestrictionAnalyzer.CreateMap(dna, topology, "EcoRI", "MboI").Sites));
            var s0 = RestrictionAnalyzer.GetDigestSummary(dna, topology, "EcoRI", "MboI");
            var s1 = RestrictionAnalyzer.GetDigestSummary(dna, topology, off, "EcoRI", "MboI");
            Assert.That(s1.FragmentSizes, Is.EqualTo(s0.FragmentSizes));
            Assert.That(s1.EnzymesUsed, Is.EqualTo(s0.EnzymesUsed));
        }
    }

    [Test]
    public void WithoutMethylation_L1_SitesAreTheEmbossNoMethylationOnes()
    {
        Assert.That(Positions(L1, "MboI", null), Is.EqualTo(new[] { 5, 15, 24 }));
        Assert.That(Positions(L1, "Sau3AI", null), Is.EqualTo(new[] { 5, 15, 24 }));
        Assert.That(Positions(L1, "BamHI", null), Is.EqualTo(new[] { 4 }));
        Assert.That(Positions(L1, "BclI", null), Is.EqualTo(new[] { 14 }));
        Assert.That(Positions(L1, "EcoRII", null), Is.EqualTo(new[] { 33, 42 }));
    }

    #endregion

    #region EMBOSS oracle — linear

    // `restrict -enzymes MboI,Sau3AI,BamHI,BclI,EcoRII -methylation` on L1: no hit at all.
    [TestCase("MboI")]
    [TestCase("Sau3AI")]
    [TestCase("BamHI")]   // window GGATCC contains the Dam-methylated A: blocked by EMBOSS (not by the bench literature)
    [TestCase("BclI")]
    [TestCase("EcoRII")]  // CCWGG: Dcm
    public void L1_EmbossBlocksEveryDamDcmWindow(string enzyme)
    {
        Assert.That(Positions(L1, enzyme, MethylaseSets.Emboss), Is.Empty);
        Assert.That(Positions(L1, enzyme, MethylaseSets.DamDcm), Is.Empty);
        Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L1), enzyme, MoleculeTopology.Linear,
            MethylaseSets.Emboss), Is.Empty, "reverse-strand records are blocked as well");
    }

    [Test]
    public void L1_DamAloneLeavesDcmSites_AndDcmAloneLeavesDamSites()
    {
        Assert.That(Positions(L1, "MboI", MethylaseSets.Dam), Is.Empty);
        Assert.That(Positions(L1, "EcoRII", MethylaseSets.Dam), Is.EqualTo(new[] { 33, 42 }));
        Assert.That(Positions(L1, "EcoRII", MethylaseSets.Dcm), Is.Empty);
        Assert.That(Positions(L1, "MboI", MethylaseSets.Dcm), Is.EqualTo(new[] { 5, 15, 24 }));
        Assert.That(Positions(L1, "BamHI", MethylaseSets.Dcm), Is.EqualTo(new[] { 4 }));
    }

    // L2, EMBOSS: EcoRI 0,18; SmaI 27; TaqI 38,44; MboI none (plain MboI: 10,36).
    [Test]
    public void L2_EmbossOutput()
    {
        var m = MethylaseSets.Emboss;
        Assert.That(Positions(L2, "MboI", null), Is.EqualTo(new[] { 10, 36 }));
        Assert.That(Positions(L2, "MboI", m), Is.Empty);
        Assert.That(Positions(L2, "EcoRI", m), Is.EqualTo(new[] { 0, 18 }));
        Assert.That(Positions(L2, "SmaI", m), Is.EqualTo(new[] { 27 }));
        // EMBOSS marks one A of GATC only: the TCGA window of GATCGA (starting at the T) is NOT blocked.
        Assert.That(Positions(L2, "TaqI", m), Is.EqualTo(new[] { 38, 44 }));
    }

    // L3, EMBOSS: EcoRI 21; Acc65I/KpnI 47 only — the GGTACC at 7 sits in the spacer of AAC(N6)GTGC and the one at
    // 35 in that of GCAC(N6)GTT (both EcoKI entries mask their N spacer).
    [TestCase("EcoRI", new[] { 21 }, new[] { 21 })]
    [TestCase("KpnI", new[] { 7, 35, 47 }, new[] { 47 })]
    [TestCase("Acc65I", new[] { 7, 35, 47 }, new[] { 47 })]
    public void L3_EcoKiSpacerIsMasked(string enzyme, int[] plain, int[] methylated)
    {
        Assert.That(Positions(L3, enzyme, null), Is.EqualTo(plain));
        Assert.That(Positions(L3, enzyme, MethylaseSets.Emboss), Is.EqualTo(methylated));
        Assert.That(Positions(L3, enzyme, MethylaseSets.DamDcm), Is.EqualTo(plain), "EcoKI is not in Dam+Dcm");
    }

    // The mask itself, from EMBOSS with the one-window "NN" probe (0-based window starts reported).
    [TestCase("AAGATCAAGATCAA", new[] { 0, 1, 4, 5, 6, 7, 10, 11, 12 })]            // one A per GATC
    [TestCase("GGCCAGGCCTGGCC", new[] { 0, 4, 5, 9, 10, 11, 12 })]                   // CmCAGG, CmCTGG
    [TestCase("CCAGGCCAGG", new[] { 2, 3, 7, 8 })]
    [TestCase("TGAACCGACCAGTGCTAACG", new[] { 0, 11, 12, 13, 14, 15, 16, 17, 18 })]   // AmAC + spacer
    [TestCase("GCACGGTACCGTTAACTT", new[] { 0, 10, 11, 12, 13, 14, 15, 16 })]         // GCmAC + spacer
    // CCTGG shares its first C with GTGC's last C; EMBOSS visits EcoKI before dcm, so both are masked.
    [TestCase("CAACCTCATTGTGCCTGGAGGG", new[] { 10, 11, 15, 16, 17, 18, 19, 20 })]
    public void Mask_EqualsTheEmbossProbe(string seq, int[] unmaskedWindowStarts)
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(seq), Probe, MoleculeTopology.Linear,
            MethylaseSets.Emboss).Where(x => x.IsForwardStrand).Select(x => x.Position).OrderBy(x => x);
        Assert.That(sites, Is.EqualTo(unmaskedWindowStarts));
    }

    [Test]
    public void RecognizedSequence_StillShowsTheRealBases()
    {
        var site = RestrictionAnalyzer.FindSites(new DnaSequence(L2), "EcoRI", MoleculeTopology.Linear,
            MethylaseSets.Emboss).First();
        Assert.That(site.RecognizedSequence, Is.EqualTo("GAATTC"));
    }

    [Test]
    public void EnzymeOverload_BehavesLikeTheNameOverload()
    {
        var mboi = RestrictionAnalyzer.GetEnzyme("MboI")!;
        Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L2), mboi, MoleculeTopology.Linear, MethylaseSets.Dam),
            Is.Empty);
        Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(L2), mboi, MoleculeTopology.Linear, null).Count(),
            Is.EqualTo(4));
    }

    #endregion

    #region EMBOSS oracle — circular

    // EMBOSS -plasmid, 46 bp: meth EcoRI 6,26 only. BglII 41 (AGATC|T across the origin, GATC wholly at the end)
    // and MboI 21,42 are blocked.
    [Test]
    public void Circular_WindowWrappingTheOrigin_ContainingMethylatedBase_IsBlocked()
    {
        const string c1 = "TAAAAAGAATTCAAAAAAAAAGATCAGAATTCAAAAGAAAAAGATC";
        Assert.That(Positions(c1, "BglII", null, MoleculeTopology.Circular), Is.EqualTo(new[] { 41 }));
        Assert.That(Positions(c1, "BglII", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.Empty);
        Assert.That(Positions(c1, "MboI", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.Empty);
        Assert.That(Positions(c1, "EcoRI", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.EqualTo(new[] { 6, 26 }));
    }

    // EMBOSS -plasmid, 39 bp: meth EcoRI 9 and 37 (GA|ATTC across the origin: unaffected), MboI 21 blocked.
    [Test]
    public void Circular_UnmethylatedWrappingSite_IsKept()
    {
        const string c2 = "ATTCAAAAAGAATTCAAAAAAGATCAAAAAAAAAAAAGA";
        Assert.That(Positions(c2, "EcoRI", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.EqualTo(new[] { 9, 37 }));
        Assert.That(Positions(c2, "MboI", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.Empty);
    }

    // Deliberate difference from EMBOSS: CC|AGG across the origin is a Dcm site on a circle. EMBOSS -plasmid
    // (masking the linear string) still reports EcoRII at 31; here it is blocked. EcoRI 7 agrees with EMBOSS.
    [Test]
    public void Circular_MethylaseSiteAcrossTheOrigin_IsFound()
    {
        const string c3 = "AGGAAAAGAATTCAAAAAAAAAGATCAAAAACC";
        Assert.That(Positions(c3, "EcoRII", null, MoleculeTopology.Circular), Is.EqualTo(new[] { 31 }));
        Assert.That(Positions(c3, "EcoRII", MethylaseSets.Dcm, MoleculeTopology.Circular), Is.Empty);
        Assert.That(Positions(c3, "EcoRI", MethylaseSets.Emboss, MoleculeTopology.Circular), Is.EqualTo(new[] { 7 }));
        Assert.That(Positions(c3, "MboI", MethylaseSets.Dcm, MoleculeTopology.Circular), Is.EqualTo(new[] { 22 }));
    }

    // A circle has no origin: rotating it must move every surviving site by the rotation, nothing else.
    [Test]
    public void Circular_RotationInvariant()
    {
        string seq = "AGGAAAAGAATTCAAAAAAAAAGATCAAAAACCTGGTTTTCAACGGTACCGTGCAAAAGATCGATAA";
        int n = seq.Length;
        foreach (var enzyme in new[] { "EcoRI", "EcoRII", "MboI", "KpnI", "TaqI", "BglII", "BclI" })
        {
            var baseline = Positions(seq, enzyme, MethylaseSets.Emboss, MoleculeTopology.Circular);
            for (int r = 1; r < n; r += 7)
            {
                string rot = seq[r..] + seq[..r];
                var shifted = Positions(rot, enzyme, MethylaseSets.Emboss, MoleculeTopology.Circular)
                    .Select(p => (p + r) % n).OrderBy(p => p);
                Assert.That(shifted, Is.EqualTo(baseline), $"{enzyme} rotated by {r}");
            }
        }
    }

    #endregion

    #region Digest / summary / map

    [Test]
    public void Digest_BlockedSitesDoNotCut()
    {
        var dna = new DnaSequence(L2);
        var plain = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, "EcoRI", "MboI").Select(f => f.Length).ToArray();
        var meth = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, MethylaseSets.Dam, "EcoRI", "MboI")
            .Select(f => f.Length).ToArray();
        var onlyEcoRi = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, "EcoRI").Select(f => f.Length).ToArray();

        Assert.That(plain.Length, Is.GreaterThan(meth.Length));
        Assert.That(meth, Is.EqualTo(onlyEcoRi));
        Assert.That(meth, Is.EqualTo(new[] { 1, 18, 32 }));
        Assert.That(meth.Sum(), Is.EqualTo(L2.Length));
    }

    [Test]
    public void Digest_Circular_AndEnzymeListOverload()
    {
        var dna = new DnaSequence(L2);
        var enzymes = new[] { RestrictionAnalyzer.GetEnzyme("EcoRI")!, RestrictionAnalyzer.GetEnzyme("MboI")! };
        var viaNames = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, MethylaseSets.Dam, "EcoRI", "MboI").ToList();
        var viaList = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, enzymes, MethylaseSets.Dam).ToList();
        Assert.That(viaList, Is.EqualTo(viaNames));
        Assert.That(viaNames.Select(f => f.Length), Is.EqualTo(new[] { 18, 33 }));
    }

    [Test]
    public void GetDigestSummary_AndCreateMap_CarryTheSet()
    {
        var dna = new DnaSequence(L2);
        var summary = RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Linear, MethylaseSets.Dam, "EcoRI", "MboI");
        Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 32, 18, 1 }));

        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, MethylaseSets.Dam, "EcoRI", "MboI");
        Assert.That(map.NonCutters, Does.Contain("MboI"));
        Assert.That(map.SitesByEnzyme.ContainsKey("MboI"), Is.False);
        Assert.That(map.SitesByEnzyme["EcoRI"], Is.EqualTo(new[] { 0, 0, 18, 18 }));

        var enzymes = new[] { RestrictionAnalyzer.GetEnzyme("EcoRI")!, RestrictionAnalyzer.GetEnzyme("MboI")! };
        Assert.That(RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, enzymes, MethylaseSets.Dam).Sites,
            Is.EqualTo(map.Sites));
        var viaList = RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Linear, enzymes, MethylaseSets.Dam);
        Assert.That(viaList.FragmentSizes, Is.EqualTo(summary.FragmentSizes));
        Assert.That(viaList.EnzymesUsed, Is.EqualTo(new[] { "EcoRI", "MboI" }));
    }

    #endregion
}
