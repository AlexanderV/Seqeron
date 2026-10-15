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
}
