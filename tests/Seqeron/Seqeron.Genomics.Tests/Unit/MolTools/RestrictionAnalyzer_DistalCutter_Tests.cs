// Canonical Test File: RESTR-FIND-001 / RESTR-DIGEST-001 — single-sided distal (Type IIS) cutters
// Methods Under Test: RestrictionAnalyzer.FindSites(DnaSequence, RestrictionEnzyme, MoleculeTopology),
//                     Digest(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>),
//                     PartialDigest(DnaSequence, MoleculeTopology, int?, IEnumerable<RestrictionEnzyme>)
//
// Reference: Biopython 1.88 Bio.Restriction (installed locally) — every expected value below is
// that library's own output: enzyme.search(Seq(s), linear=...) (1-based top-strand cuts),
// enzyme.catalyse(...) fragment lengths, and the end window between cut_on_watson and
// cut_on_crick = cut_on_watson - ovhg (AbstractCut._drop). Offset mapping:
// CutPositionForward = fst5, CutPositionReverse = size + fst3 (review 2026-09 B08 F54/F55).
//
// TestSpec: TestSpecs/RESTR-FIND-001.md, TestSpecs/RESTR-DIGEST-001.md

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class RestrictionAnalyzer_DistalCutter_Tests
{
    // Biopython 1.88 (site, fst5, fst3) → (site, fst5, size + fst3).
    private static readonly RestrictionEnzyme BsaI = new("BsaI", "GGTCTC", 7, 11, "test");   // GGTCTC(1/5): 7, 5
    private static readonly RestrictionEnzyme BsmBI = new("BsmBI", "CGTCTC", 7, 11, "test"); // CGTCTC(1/5): 7, 5
    private static readonly RestrictionEnzyme BbsI = new("BbsI", "GAAGAC", 8, 12, "test");   // GAAGAC(2/6): 8, 6
    private static readonly RestrictionEnzyme SapI = new("SapI", "GCTCTTC", 8, 11, "test");  // GCTCTTC(1/4): 8, 4
    private static readonly RestrictionEnzyme FokI = new("FokI", "GGATG", 14, 18, "test");   // GGATG(9/13): 14, 13
    private static readonly RestrictionEnzyme MlyI = new("MlyI", "GAGTC", 10, 10, "test");   // GAGTC(5/5): 10, 5
    private static readonly RestrictionEnzyme AcuI = new("AcuI", "CTGAAG", 22, 20, "test");  // CTGAAG(16/14): 22, 14
    private static readonly RestrictionEnzyme TspRI = new("TspRI", "CASTG", 7, -2, "test");  // CASTG: 7, -7

    private static List<DigestFragment> Digest(string s, MoleculeTopology t, RestrictionEnzyme e) =>
        RestrictionAnalyzer.Digest(new DnaSequence(s), t, new[] { e }).ToList();

    /// <summary>The top-strand cut of any site record (forward: CutPosition; reverse: see RestrictionSite).</summary>
    private static int TopCut(RestrictionSite s, int n) =>
        s.IsForwardStrand
            ? s.CutPosition
            : ((s.CutPosition + s.Enzyme.CutPositionForward - s.Enzyme.CutPositionReverse) % n + n) % n;

    #region Mapping

    [Test]
    public void DistalEnzyme_OffsetsMapToBiopythonFst5Fst3AndOvhg()
    {
        Assert.Multiple(() =>
        {
            // ovhg = fst5 - (size + fst3) = CutPositionForward - CutPositionReverse
            Assert.That(BsaI.CutPositionForward - BsaI.CutPositionReverse, Is.EqualTo(-4)); // BsaI.ovhg
            Assert.That(BsaI.OverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(SapI.CutPositionForward - SapI.CutPositionReverse, Is.EqualTo(-3)); // SapI.ovhg
            Assert.That(MlyI.OverhangType, Is.EqualTo(OverhangType.Blunt));                 // MlyI.ovhg 0
            Assert.That(AcuI.OverhangType, Is.EqualTo(OverhangType.ThreePrime));            // AcuI.ovhg 2
            // TspRI: fst5 7, fst3 -7 -> bottom cut 2 nt BEFORE the site; 9-nt 3' overhang NNCASTGNN
            Assert.That(TspRI.CutPositionForward - TspRI.CutPositionReverse, Is.EqualTo(9));
        });
    }

    #endregion

    #region FindSites

    [Test]
    public void FindSites_BsaI_Linear_ReverseStrandSiteCutsOnTheOtherSide()
    {
        const string s = "ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT";
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(s), BsaI, MoleculeTopology.Linear).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(2));
            var f = sites.Single(x => x.IsForwardStrand);
            var r = sites.Single(x => !x.IsForwardStrand);
            Assert.That((f.Position, f.CutPosition), Is.EqualTo((2, 9)), "GGTCTC at 2, cut 7 downstream");
            Assert.That((r.Position, r.CutPosition), Is.EqualTo((33, 32)),
                "GAGACC at 33: bottom-strand cut 1 nt 5' of the window");
            // Biopython BsaI.search(seq) == [10, 29] (1-based top-strand cuts)
            Assert.That(sites.Select(x => TopCut(x, s.Length) + 1).Order(), Is.EqualTo(new[] { 10, 29 }));
        });
    }

    [Test]
    public void FindSites_FokI_Circular_ReverseSiteNearOrigin_CutIsReducedIntoTheMolecule()
    {
        // CATCC (reverse FokI) at 0: top-strand cut 0+5-18 = -13, bottom 0+5-14 = -9 — both upstream
        // of the origin. Biopython FokI.search(seq, linear=False) == [30].
        const string s = "CATCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
        var site = RestrictionAnalyzer.FindSites(new DnaSequence(s), FokI, MoleculeTopology.Circular).Single();

        Assert.Multiple(() =>
        {
            Assert.That(site.IsForwardStrand, Is.False);
            Assert.That(site.CutPosition, Is.EqualTo(33), "-9 mod 42; was -9 before F54");
            Assert.That(TopCut(site, s.Length) + 1, Is.EqualTo(30));
        });
        Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(s), FokI, MoleculeTopology.Linear), Is.Empty,
            "Biopython FokI.search(seq) == [] on the linear molecule: the cut falls before base 0");
    }

    [TestCase("AAGGATGAAAAAAAAACCCCCTTT", 17)] // Biopython [17]: bottom cut 20 <= 23
    [TestCase("AAGGATGAAAAAAAAACCCC", null)]   // Biopython []: bottom cut 20 is the molecule end
    public void FindSites_FokI_Linear_SiteIsDroppedWhenTheCutFallsOffTheEnd(string s, int? biopython)
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(s), FokI, MoleculeTopology.Linear).ToList();
        Assert.That(sites.Select(x => TopCut(x, s.Length) + 1),
            Is.EqualTo(biopython is null ? Array.Empty<int>() : new[] { biopython.Value }));
    }

    [Test]
    public void FindSites_TspRI_BottomCutAfterBaseZero_IsKept_DocumentedDivergence()
    {
        // TspRI CASTG at 3: top cut 10, bottom cut 3-2 = 1 (the bond after base 0, which exists).
        // Biopython 1.88 TspRI.search(Seq("AAACAGTGAAAAAA")) == [] only because its stored
        // ovhg = 10 contradicts its own fst5/fst3 and 9-nt ovhgseq "NNCASTGNN" (F54).
        const string s = "AAACAGTGAAAAAA";
        var site = RestrictionAnalyzer.FindSites(new DnaSequence(s), TspRI, MoleculeTopology.Linear)
            .Single(x => x.IsForwardStrand);
        Assert.That(site.CutPosition, Is.EqualTo(10));

        var fragments = Digest(s, MoleculeTopology.Linear, TspRI);
        Assert.Multiple(() =>
        {
            Assert.That(fragments.Select(f => f.Length), Is.EqualTo(new[] { 10, 4 }));
            Assert.That(fragments[0].RightOverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(fragments[0].RightOverhang, Is.EqualTo("AACAGTGAA"), "the 9 bases [1, 10)");
        });
    }

    #endregion

    #region Digest

    // Expected: Biopython catalyse lengths; overhangs = window [min(w, w-ovhg), max) of the molecule,
    // 5' overhang → downstream left end reads the window, upstream right end its reverse complement.
    [TestCase("BsaI", "ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT", new[] { 9, 19, 13 }, new[] { "AATG", "CATT" }, new[] { "CATT", "AATG" })]
    [TestCase("BsmBI", "CGTCTCAGGAGAAAAAAAAAAACGCTTGAGACG", new[] { 7, 15, 11 }, new[] { "GGAG", "CGCT" }, new[] { "CTCC", "AGCG" })]
    [TestCase("BbsI", "GGAAGACAAGCTTTTTTTTTTAAACTGTCTTCC", new[] { 9, 11, 13 }, new[] { "GCTT", "TAAA" }, new[] { "AAGC", "TTTA" })]
    [TestCase("SapI", "AAGCTCTTCAGCATTTTTTTTGAAGAGCTT", new[] { 10, 7, 13 }, new[] { "GCA", "TTT" }, new[] { "TGC", "AAA" })]
    [TestCase("MlyI", "TTTTTGAGTCAAAAACCCCCGGGGGACTCAAA", new[] { 15, 4, 13 }, new[] { "", "" }, new[] { "", "" })]
    [TestCase("AcuI", "CTGAAGTTTTTTTTTTTTTTTTACGTAAAAAA", new[] { 22, 10 }, new[] { "AA" }, new[] { "TT" })]
    [TestCase("FokI", "AAGGATGAAAAAAAAACCCCCTTT", new[] { 16, 8 }, new[] { "CCCC" }, new[] { "GGGG" })]
    public void Digest_Linear_DistalCutter_MatchesBiopythonCatalyseAndEndWindows(
        string enzyme, string s, int[] lengths, string[] downstreamLeft, string[] upstreamRight)
    {
        var e = new[] { BsaI, BsmBI, BbsI, SapI, MlyI, AcuI, FokI }.Single(x => x.Name == enzyme);
        var fragments = Digest(s, MoleculeTopology.Linear, e);

        Assert.Multiple(() =>
        {
            Assert.That(fragments.Select(f => f.Length), Is.EqualTo(lengths));
            Assert.That(fragments.Skip(1).Select(f => f.LeftOverhang), Is.EqualTo(downstreamLeft));
            Assert.That(fragments.SkipLast(1).Select(f => f.RightOverhang), Is.EqualTo(upstreamRight));
            Assert.That(fragments.Skip(1).Select(f => f.LeftOverhangType), Is.All.EqualTo(e.OverhangType));
            Assert.That(fragments.Select(f => f.Sequence), Is.EqualTo(
                Enumerable.Range(0, lengths.Length).Select(i => s.Substring(lengths.Take(i).Sum(), lengths[i]))));
        });
    }

    [Test]
    public void Digest_BsaI_GoldenGateInsert_CarriesBothFourBaseEnds()
    {
        // Two inward-facing BsaI sites release the insert with its two 4-nt 5' overhangs taken from
        // the molecule (Biopython ovhgseq is only "NNNN"); fragment lengths = catalyse [9, 19, 13].
        var insert = Digest("ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT", MoleculeTopology.Linear, BsaI)[1];
        Assert.Multiple(() =>
        {
            Assert.That(insert.Sequence, Is.EqualTo("AATGCCCCCCCCCCCCCCG"));
            Assert.That((insert.LeftOverhangType, insert.LeftOverhang), Is.EqualTo((OverhangType.FivePrime, "AATG")));
            Assert.That((insert.RightOverhangType, insert.RightOverhang), Is.EqualTo((OverhangType.FivePrime, "AATG")));
            Assert.That(insert.LeftEnzymes, Is.EqualTo(new[] { "BsaI" }));
            Assert.That(insert.RightEnzymes, Is.EqualTo(new[] { "BsaI" }));
        });
    }

    [Test]
    public void Digest_BsaI_Circular_MatchesBiopython()
    {
        // Biopython BsaI.search(seq, linear=False) == [10, 29]; catalyse gives [22, 19] with the
        // origin-spanning fragment first (Seqeron: last, F27).
        var fragments = Digest("ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT", MoleculeTopology.Circular, BsaI);
        Assert.Multiple(() =>
        {
            Assert.That(fragments.Select(f => (f.StartPosition, f.Length)), Is.EqualTo(new[] { (9, 19), (28, 22) }));
            Assert.That(fragments[1].Sequence, Is.EqualTo("CATTTGAGACCATATGGTCTCA"));
        });
    }

    [Test]
    public void Digest_BsaI_Circular_SiteSpanningTheOrigin()
    {
        // GGTC|TC spans the origin of "TCAAATGCCCCCCCCATGGTC"; Biopython search(linear=False) == [4],
        // catalyse → one 21-nt fragment; the cut's window is AATG.
        var fragment = Digest("TCAAATGCCCCCCCCATGGTC", MoleculeTopology.Circular, BsaI).Single();
        Assert.Multiple(() =>
        {
            Assert.That((fragment.StartPosition, fragment.Length), Is.EqualTo((3, 21)));
            Assert.That(fragment.LeftOverhang, Is.EqualTo("AATG"));
            Assert.That(fragment.RightOverhang, Is.EqualTo("CATT"));
        });
    }

    [Test]
    public void Digest_FokI_Circular_ReverseSiteCuttingAcrossTheOrigin()
    {
        // Biopython FokI.search(seq, linear=False) == [30]; one cut linearises the 42-nt circle.
        var fragment = Digest("CATCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG", MoleculeTopology.Circular, FokI).Single();
        Assert.Multiple(() =>
        {
            Assert.That((fragment.StartPosition, fragment.Length), Is.EqualTo((29, 42)));
            Assert.That(fragment.LeftOverhang, Is.EqualTo("AAAA"));
            Assert.That(fragment.RightOverhang, Is.EqualTo("TTTT"));
        });
    }

    [Test]
    public void Digest_NonPalindromicInSiteCutter_ReverseSitesAreBoundariesToo()
    {
        // Before F55 the digest kept forward-strand records only, so the reverse-strand sites of a
        // non-palindromic enzyme were not cut. MlyI: GAGTC at 5 (cut 15) and GACTC at 22 (cut 19).
        var fragments = Digest("TTTTTGAGTCAAAAACCCCCGGGGGACTCAAA", MoleculeTopology.Linear, MlyI);
        Assert.That(fragments.Select(f => f.StartPosition), Is.EqualTo(new[] { 0, 15, 19 }));
    }

    [Test]
    public void Digest_EnzymeOverload_BuiltInEnzymeRecord_EqualsNameOverload()
    {
        var dna = new DnaSequence("GAATTCAAAGGATCCAAAGAATTCTTTAAGCTT");
        var ecoRI = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        var bamHI = RestrictionAnalyzer.GetEnzyme("BamHI")!;
        foreach (var t in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            Assert.That(RestrictionAnalyzer.Digest(dna, t, new[] { bamHI, ecoRI }).ToList(),
                Is.EqualTo(RestrictionAnalyzer.Digest(dna, t, "EcoRI", "BamHI").ToList()), t.ToString());
        }
    }

    [Test]
    public void Digest_EnzymeOverload_InvalidArguments_Throw()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(() => RestrictionAnalyzer.Digest(null!, MoleculeTopology.Linear, new[] { BsaI }),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, (IEnumerable<RestrictionEnzyme>)null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, Array.Empty<RestrictionEnzyme>()),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, new RestrictionEnzyme[] { null! }),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear,
                new[] { new RestrictionEnzyme("Bad", "GGXC", 1, 3, "test") }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Linear, -1, new[] { BsaI }),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Linear, null, Array.Empty<RestrictionEnzyme>()),
                NUnit.Framework.Throws.ArgumentException);
        });
    }

    #endregion

    #region PartialDigest

    [TestCase(MoleculeTopology.Linear)]
    [TestCase(MoleculeTopology.Circular)]
    public void PartialDigest_BsaI_ZeroMissedSubsetIsTheDigest_AndLengthsAreTheCutDifferences(MoleculeTopology t)
    {
        // Biopython cuts {9, 28} (0-based) on a 41-nt molecule.
        var dna = new DnaSequence("ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT");
        var partial = RestrictionAnalyzer.PartialDigest(dna, t, null, new[] { BsaI }).ToList();
        var complete = RestrictionAnalyzer.Digest(dna, t, new[] { BsaI }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(partial.Where(p => p.MissedCuts == 0).Select(p => p.Fragment), Is.EqualTo(complete));
            var expected = t == MoleculeTopology.Linear
                ? new[] { 9, 13, 19, 28, 32, 41 }   // all intervals over {0, 9, 28, 41}
                : new[] { 19, 22, 41, 41, 41 };     // k^2 + 1 = 5
            Assert.That(partial.Select(p => p.Fragment.Length).Order(), Is.EqualTo(expected));
        });
    }

    #endregion

    #region Two-sided (Type IIB) cutters — review 2026-09 B08 WP11b F56

    // Biopython 1.88 TwoCuts (fst5, size + fst3) and (scd5, size + scd3).
    private static readonly RestrictionEnzyme BaeI =
        new("BaeI", "ACNNNNGTAYC", -10, -15, "test") { SecondCut = (23, 18) };   // fst5 -10 fst3 -26 scd5 23 scd3 7
    private static readonly RestrictionEnzyme UcoMSI =
        new("UcoMSI", "GAGCTC", -7, -5, "test") { SecondCut = (11, 13) };        // palindromic; fst5 -7 fst3 -11 scd5 11 scd3 7
    private static readonly RestrictionEnzyme BcgI =
        new("BcgI", "CGANNNNNNTGC", -10, -12, "test") { SecondCut = (24, 22) };  // fst5 -10 fst3 -24 scd5 24 scd3 10

    private const string BaeIForward = "AAAAAAAAAAAAAAAAAAAAACAAAAGTACCTTTTTTTTTTTTTTTTTTTT"; // A20 + ACAAAAGTACC + T20

    [Test]
    public void TwoCut_FindSites_ForwardSite_ReportsBothCuts()
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(BaeIForward), BaeI, MoleculeTopology.Linear).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(x => (x.Position, x.IsForwardStrand, x.IsSecondCut, x.CutPosition)),
                Is.EqualTo(new[] { (20, true, false, 10), (20, true, true, 43) }));
            // Biopython BaeI.search(seq) == [11, 44]
            Assert.That(sites.Select(x => TopCutPair(x, BaeIForward.Length) + 1), Is.EqualTo(new[] { 11, 44 }));
            Assert.That(sites[1].CutOffsets, Is.EqualTo((23, 18)));
        });
    }

    [Test]
    public void TwoCut_Digest_Linear_BothSidesExcisesTheSite_EndsMatchBiopython()
    {
        var f = Digest(BaeIForward, MoleculeTopology.Linear, BaeI);

        Assert.Multiple(() =>
        {
            // BaeI.catalyse(seq) fragment lengths [10, 33, 8]
            Assert.That(f.Select(x => (x.StartPosition, x.Length)), Is.EqualTo(new[] { (0, 10), (10, 33), (43, 8) }));
            // ovhg 5 (3'): window between cut_on_watson and cut_on_watson - ovhg is AAAAA at 11, TTTTT at 44.
            Assert.That((f[0].RightOverhangType, f[0].RightOverhang), Is.EqualTo(((OverhangType?)OverhangType.ThreePrime, "AAAAA")));
            Assert.That((f[1].LeftOverhang, f[1].RightOverhang), Is.EqualTo(("TTTTT", "TTTTT")));
            Assert.That(f[2].LeftOverhang, Is.EqualTo("AAAAA"));
            Assert.That(f[1].Sequence, Does.Contain("ACAAAAGTACC"), "the excised middle fragment carries the site");
        });
    }

    [Test]
    public void TwoCut_Digest_ReverseStrandSite_CutsOnBothSides()
    {
        // Reverse complement of BaeIForward: BaeI.search == [14, 47], catalyse lengths [13, 33, 5].
        string s = DnaSequence.GetReverseComplementString(BaeIForward);
        var f = Digest(s, MoleculeTopology.Linear, BaeI);
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(s), BaeI).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.All(x => !x.IsForwardStrand), Is.True);
            Assert.That(sites.Select(x => TopCutPair(x, s.Length) + 1).Order(), Is.EqualTo(new[] { 14, 47 }));
            Assert.That(f.Select(x => x.Length), Is.EqualTo(new[] { 13, 33, 5 }));
            Assert.That((f[0].RightOverhang, f[1].RightOverhang), Is.EqualTo(("AAAAA", "TTTTT")));
        });
    }

    [Test]
    public void TwoCut_Linear_CutBeyondTheEnd_IsDroppedAlone()
    {
        // Site at 0: the upstream cut (-10) is off the molecule and dropped, the downstream one kept.
        // BaeI.search(Seq("ACAAAAGTACC" + "T"*20)) == [24]; catalyse lengths [23, 8].
        const string s = "ACAAAAGTACCTTTTTTTTTTTTTTTTTTTT";
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(s), BaeI).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Select(x => (x.IsSecondCut, x.CutPosition)), Is.EqualTo(new[] { (true, 23) }));
            Assert.That(Digest(s, MoleculeTopology.Linear, BaeI).Select(x => x.Length), Is.EqualTo(new[] { 23, 8 }));
        });
    }

    [Test]
    public void TwoCut_Circular_SiteSpanningTheOrigin()
    {
        // 41-nt circle, ACAAAA|GTACC across the origin. Biopython search(linear=False) == [26, 18]
        // (its one-pass fold leaves the list unsorted); cuts 17 and 25 (0-based), ends TTTTT both.
        const string c = "GTACCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACAAAA";
        var f = Digest(c, MoleculeTopology.Circular, BaeI);

        Assert.Multiple(() =>
        {
            Assert.That(f.Select(x => (x.StartPosition, x.Length)), Is.EqualTo(new[] { (17, 8), (25, 33) }));
            Assert.That(f.Select(x => x.LeftOverhang), Is.EqualTo(new[] { "AAAAA", "AAAAA" }));
            Assert.That(f.Select(x => x.RightOverhang), Is.EqualTo(new[] { "TTTTT", "TTTTT" }));
        });
    }

    [Test]
    public void TwoCut_Palindromic_OneCutPairPerWindow_NotFour()
    {
        // UcoMSI GAGCTC, palindromic: Biopython searches one strand, search == [14, 32]; fragments
        // 13, 18, 25 (its own catalyse agrees here; with two sites its unsorted list yields an empty
        // fragment — a Biopython defect not reproduced).
        const string s = "AAAAAAAAAAAAAAAAAAAAGAGCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        var f = Digest(s, MoleculeTopology.Linear, UcoMSI);

        Assert.Multiple(() =>
        {
            Assert.That(f.Select(x => x.Length), Is.EqualTo(new[] { 13, 18, 25 }));
            Assert.That(f[1].LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That((f[1].LeftOverhang, f[1].RightOverhang), Is.EqualTo(("AA", "TT")));
        });
    }

    [Test]
    public void TwoCut_PartialDigest_LinearIntervals()
    {
        var p = RestrictionAnalyzer.PartialDigest(new DnaSequence(BaeIForward), MoleculeTopology.Linear, null, new[] { BaeI }).ToList();
        // all intervals over {0, 10, 43, 51}
        Assert.That(p.Select(x => x.Fragment.Length).Order(), Is.EqualTo(new[] { 8, 10, 33, 41, 43, 51 }));
    }

    [Test]
    public void TwoCut_SecondPairWithDifferentStagger_Throws()
    {
        var bad = BaeI with { SecondCut = (23, 20) };
        Assert.Multiple(() =>
        {
            Assert.That(() => RestrictionAnalyzer.FindSites(new DnaSequence(BaeIForward), bad).ToList(),
                NUnit.Framework.Throws.ArgumentException.With.Message.Contains("stagger"));
            Assert.That(() => RestrictionAnalyzer.Digest(new DnaSequence(BaeIForward), MoleculeTopology.Linear, new[] { bad }),
                NUnit.Framework.Throws.ArgumentException);
        });
    }

    [Test]
    public void OneCut_Enzyme_HasNoSecondCut_AndRecordsAreNeverSecond()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.Enzymes.Values.All(e => e.SecondCut is null), Is.True);
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence("GGTCTCAAAAAAAAAA"), BsaI).All(x => !x.IsSecondCut), Is.True);
        });
    }

    /// <summary>Top-strand cut of a record, with that record's own cut pair.</summary>
    private static int TopCutPair(RestrictionSite s, int n)
    {
        var (fwd, rev) = s.CutOffsets;
        return s.IsForwardStrand ? s.CutPosition : ((s.CutPosition + fwd - rev) % n + n) % n;
    }

    #endregion

    #region Enzyme-level overhang and compatibility of distal ends — F57

    // Biopython 1.88 ovhgseq: the pattern window padded with N outside the site.
    [TestCase("GGTCTC", 7, 11, "NNNN")]        // BsaI
    [TestCase("GCTCTTC", 8, 11, "NNN")]        // SapI
    [TestCase("CTGAAG", 22, 20, "NN")]         // AcuI (3')
    [TestCase("CASTG", 7, -2, "NNCASTGNN")]    // TspRI
    [TestCase("GAATGC", 7, 5, "CN")]           // BsmI: half inside, half outside
    [TestCase("GAGTC", 10, 10, "")]            // MlyI blunt
    public void GetOverhangSequence_DistalEnzyme_MatchesBiopythonOvhgseq(string site, int f, int r, string expected)
    {
        Assert.That(RestrictionAnalyzer.GetOverhangSequence(new RestrictionEnzyme("x", site, f, r, "test")),
            Is.EqualTo(expected));
    }

    [Test]
    public void GetOverhangSequence_TwoCutEnzyme_IsTheFirstPairsWindow()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(BaeI), Is.EqualTo("NNNNN"));   // BaeI.ovhgseq
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(BcgI), Is.EqualTo("NN"));      // BcgI.ovhgseq
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(UcoMSI), Is.EqualTo("NN"));    // UcoMSI.ovhgseq
        });
    }

    // Every verdict is Biopython 1.88 `a % b` (equal to `b % a` for each of these pairs).
    [TestCase("BsaI", "EcoRI", true)]
    [TestCase("BsaI", "BamHI", true)]
    [TestCase("BsaI", "BsaI", true)]
    [TestCase("BsaI", "BsmBI", true)]
    [TestCase("BsaI", "PstI", false)]   // 3' TGCA
    [TestCase("BsaI", "SapI", false)]   // 3-nt
    [TestCase("BsaI", "EcoRV", false)]  // blunt
    [TestCase("MlyI", "EcoRV", true)]   // blunt distal with blunt
    [TestCase("FokI", "BamHI", true)]
    [TestCase("AcuI", "PstI", false)]   // 3' NN vs 3' TGCA: length
    [TestCase("AcuI", "BcgI", true)]    // 3' NN with a two-cut 3' NN
    [TestCase("BaeI", "BcgI", false)]   // 3' NNNNN vs NN
    [TestCase("TspRI", "TspRI", true)]
    public void AreCompatible_DistalAndTwoCutEnds_MatchBiopythonMod(string a, string b, bool expected)
    {
        var table = new Dictionary<string, RestrictionEnzyme>
        {
            ["BsaI"] = BsaI, ["BsmBI"] = BsmBI, ["SapI"] = SapI, ["MlyI"] = MlyI, ["FokI"] = FokI,
            ["AcuI"] = AcuI, ["TspRI"] = TspRI, ["BaeI"] = BaeI, ["BcgI"] = BcgI,
        };
        RestrictionEnzyme E(string n) => table.TryGetValue(n, out var e) ? e : RestrictionAnalyzer.GetEnzyme(n)!;

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.AreCompatible(E(a), E(b)), Is.EqualTo(expected), $"{a} % {b}");
            Assert.That(RestrictionAnalyzer.AreCompatible(E(b), E(a)), Is.EqualTo(expected), $"{b} % {a}");
        });
    }

    #endregion

    #region Custom enzymes in FindSites (multi), GetDigestSummary, CreateMap — WP11b part 3

    private const string BsaISeq = "ATGGTCTCAAATGCCCCCCCCCCCCCCGCATTTGAGACCAT"; // Biopython BsaI.search == [10, 29]

    [Test]
    public void FindSites_EnzymeCollection_IsTheConcatenationOfSingleSearches()
    {
        var dna = new DnaSequence(BsaISeq);
        var both = RestrictionAnalyzer.FindSites(dna, MoleculeTopology.Linear, new[] { BsaI, FokI }).ToList();

        Assert.That(both, Is.EqualTo(RestrictionAnalyzer.FindSites(dna, BsaI, MoleculeTopology.Linear)
            .Concat(RestrictionAnalyzer.FindSites(dna, FokI, MoleculeTopology.Linear)).ToList()));
    }

    [TestCase(MoleculeTopology.Linear, new[] { 19, 13, 9 })]
    [TestCase(MoleculeTopology.Circular, new[] { 22, 19 })]
    public void GetDigestSummary_EnzymeCollection_BsaI(MoleculeTopology t, int[] sizes)
    {
        // Cuts {9, 28} (0-based) on 41 nt: linear 9/19/13, circular 19/22.
        var summary = RestrictionAnalyzer.GetDigestSummary(new DnaSequence(BsaISeq), t, new[] { BsaI });

        Assert.Multiple(() =>
        {
            Assert.That(summary.FragmentSizes, Is.EqualTo(sizes));
            Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "BsaI" }));
            Assert.That(summary.Topology, Is.EqualTo(t));
        });
    }

    [Test]
    public void CreateMap_EnzymeCollection_CountsNonPalindromicReverseSites()
    {
        // BsaI: GGTCTC at 2 (forward) and GAGACC at 33 (reverse) — two physical sites;
        // FokI (GGATG/CATCC) absent.
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(BsaISeq), MoleculeTopology.Linear, new[] { BsaI, FokI });

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2));
            Assert.That(map.SitesByEnzyme["BsaI"], Is.EqualTo(new[] { 2, 33 }));
            Assert.That(map.UniqueCutters, Is.Empty);
            Assert.That(map.NonCutters, Is.EqualTo(new[] { "FokI" }));
        });
    }

    [Test]
    public void CreateMap_EnzymeCollection_TwoCutSiteCountsOnce()
    {
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(BaeIForward), MoleculeTopology.Linear, new[] { BaeI });

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(1));
            Assert.That(map.UniqueCutters, Is.EqualTo(new[] { "BaeI" }));
            Assert.That(map.Sites, Has.Count.EqualTo(2), "one record per cut");
        });
    }

    [Test]
    public void EnzymeCollectionOverloads_BuiltInRecords_EqualTheNameOverloads()
    {
        var dna = new DnaSequence("GAATTCAAAGGATCCAAAGAATTCAAACTGCAGAAAGCGGCCGCAAA");
        var names = new[] { "EcoRI", "BamHI", "PstI", "NotI", "SmaI" };
        var records = names.Select(n => RestrictionAnalyzer.GetEnzyme(n)!).ToArray();

        foreach (var t in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var byName = RestrictionAnalyzer.CreateMap(dna, t, names);
            var byRecord = RestrictionAnalyzer.CreateMap(dna, t, records);
            Assert.Multiple(() =>
            {
                Assert.That(byRecord.TotalSites, Is.EqualTo(byName.TotalSites), t.ToString());
                Assert.That(byRecord.UniqueCutters, Is.EqualTo(byName.UniqueCutters), t.ToString());
                Assert.That(byRecord.NonCutters, Is.EqualTo(byName.NonCutters), t.ToString());
                Assert.That(byRecord.Sites, Is.EqualTo(byName.Sites), t.ToString());
                Assert.That(RestrictionAnalyzer.GetDigestSummary(dna, t, records).FragmentSizes,
                    Is.EqualTo(RestrictionAnalyzer.GetDigestSummary(dna, t, names).FragmentSizes), t.ToString());
                Assert.That(RestrictionAnalyzer.FindSites(dna, t, records).ToList(),
                    Is.EqualTo(RestrictionAnalyzer.FindSites(dna, t, names).ToList()), t.ToString());
            });
        }
    }

    [Test]
    public void EnzymeCollectionOverloads_InvalidArguments_Throw()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(() => RestrictionAnalyzer.FindSites(dna, MoleculeTopology.Linear, Array.Empty<RestrictionEnzyme>()),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, Array.Empty<RestrictionEnzyme>()),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Linear, (IEnumerable<RestrictionEnzyme>)null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => RestrictionAnalyzer.CreateMap(null!, MoleculeTopology.Linear, new[] { BsaI }),
                NUnit.Framework.Throws.ArgumentNullException);
        });
    }

    #endregion
}
