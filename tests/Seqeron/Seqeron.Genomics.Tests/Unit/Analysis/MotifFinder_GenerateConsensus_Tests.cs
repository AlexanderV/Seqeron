// MOTIF-GENERATE-001 — IUPAC-Degenerate Consensus Generation
// Evidence: docs/Evidence/MOTIF-GENERATE-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-GENERATE-001.md
// Source: Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic
//         acid sequences: recommendations 1984. Nucleic Acids Research 13(9):3021. DOI 10.1093/nar/13.9.3021.
//         UCSC IUPAC ambiguity codes; Wikipedia "Nucleic acid notation" (NC-IUB 1984 table);
//         DECIPHER ConsensusSequence (threshold-consensus family, equal-abundance rule);
//         Cavener 1987 NAR 15(4):1353 via Biopython 1.88 Bio.motifs degenerate_consensus.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical test class for MOTIF-GENERATE-001: IUPAC-degenerate consensus generation. Verifies
/// <see cref="MotifFinder.GenerateConsensus(System.Collections.Generic.IEnumerable{string})"/>
/// against the NC-IUB 1984 set→symbol mapping and the documented 25% inclusion threshold.
/// Expected symbols are derived from the authoritative IUPAC table, not from the implementation.
/// </summary>
[TestFixture]
public class MotifFinder_GenerateConsensus_Tests
{
    #region GenerateConsensus — MUST (NC-IUB set→symbol mapping)

    // M1 — {A,G} both above threshold → R (purine). Source: NC-IUB 1984 / UCSC.
    // n=2, threshold=0.5; A=1>0.5, G=1>0.5 → set {A,G} → R.
    [Test]
    public void GenerateConsensus_ColumnAG_ReturnsR()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "G" });
        Assert.That(consensus, Is.EqualTo("R"),
            "NC-IUB: the set {A,G} (purine) is encoded by the IUPAC symbol R.");
    }

    // M2 — {C,T} → Y (pyrimidine). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCT_ReturnsY()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "T" });
        Assert.That(consensus, Is.EqualTo("Y"),
            "NC-IUB: the set {C,T} (pyrimidine) is encoded by the IUPAC symbol Y.");
    }

    // M3 — {C,G} → S (strong). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCG_ReturnsS()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "G" });
        Assert.That(consensus, Is.EqualTo("S"),
            "NC-IUB: the set {C,G} (strong) is encoded by the IUPAC symbol S.");
    }

    // M4 — {A,T} → W (weak). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAT_ReturnsW()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "T" });
        Assert.That(consensus, Is.EqualTo("W"),
            "NC-IUB: the set {A,T} (weak) is encoded by the IUPAC symbol W.");
    }

    // M5 — {G,T} → K (keto). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnGT_ReturnsK()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "G", "T" });
        Assert.That(consensus, Is.EqualTo("K"),
            "NC-IUB: the set {G,T} (keto) is encoded by the IUPAC symbol K.");
    }

    // M6 — {A,C} → M (amino). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAC_ReturnsM()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C" });
        Assert.That(consensus, Is.EqualTo("M"),
            "NC-IUB: the set {A,C} (amino) is encoded by the IUPAC symbol M.");
    }

    // M7 — {C,G,T} → B (not-A). n=3, threshold=0.75; each count=1>0.75. Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCGT_ReturnsB()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "G", "T" });
        Assert.That(consensus, Is.EqualTo("B"),
            "NC-IUB: the set {C,G,T} (not-A) is encoded by the IUPAC symbol B.");
    }

    // M8 — {A,G,T} → D (not-C). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAGT_ReturnsD()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "G", "T" });
        Assert.That(consensus, Is.EqualTo("D"),
            "NC-IUB: the set {A,G,T} (not-C) is encoded by the IUPAC symbol D.");
    }

    // M9 — {A,C,T} → H (not-G). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnACT_ReturnsH()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "T" });
        Assert.That(consensus, Is.EqualTo("H"),
            "NC-IUB: the set {A,C,T} (not-G) is encoded by the IUPAC symbol H.");
    }

    // M10 — {A,C,G} → V (not-T). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnACG_ReturnsV()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "G" });
        Assert.That(consensus, Is.EqualTo("V"),
            "NC-IUB: the set {A,C,G} (not-T) is encoded by the IUPAC symbol V.");
    }

    // M11 — Unanimous columns reproduce the input (INV-02): singleton set → standard base.
    [Test]
    public void GenerateConsensus_IdenticalSequences_ReturnsThatSequence()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ATGC", "ATGC", "ATGC" });
        Assert.That(consensus, Is.EqualTo("ATGC"),
            "Every column is unanimous; a singleton base set maps to its standard base (NC-IUB).");
    }

    // M12 — Multi-column mixed: col0 {A,G}→R, cols 1-3 unanimous → "RTGC". Source: NC-IUB 1984.
    [Test]
    public void GenerateConsensus_MixedFirstColumn_ReturnsRtgc()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ATGC", "GTGC" });
        Assert.That(consensus, Is.EqualTo("RTGC"),
            "Column 0 holds {A,G}→R; remaining columns are unanimous (T,G,C).");
    }

    // M13 — Strict 25% boundary (INV-05): a base at exactly the threshold is excluded.
    // n=4, threshold=1.0. col0/col1: A=4 → 'A'. col2: A=C=G=T=1, none >1.0 → all four tie at the
    // maximum → N (F13; Biopython degenerate_consensus of these rows is also "AANT"). col3: A=1,G=1,T=2 → only T(2)>1.0 → singleton {T} → 'T'.
    [Test]
    public void GenerateConsensus_ExactlyQuarterBoundary_ExcludesBaseAtThreshold()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "AAAA", "AAGT", "AACT", "AATT" });
        Assert.Multiple(() =>
        {
            Assert.That(consensus, Has.Length.EqualTo(4),
                "One IUPAC symbol per column (INV-01).");
            Assert.That(consensus[0], Is.EqualTo('A'), "Column 0 is all A.");
            Assert.That(consensus[1], Is.EqualTo('A'), "Column 1 is all A.");
            Assert.That(consensus[2], Is.EqualTo('N'),
                "Column 2 has four bases each at exactly 25% (count 1 = threshold 1.0); strict '>' " +
                "excludes all, and the four equally abundant bases are encoded as N.");
            Assert.That(consensus[3], Is.EqualTo('T'),
                "Column 3: only T (count 2 > threshold 1.0) passes; A and G at exactly 25% are dropped.");
        });
    }

    // M14 — Minority base below threshold is dropped before IUPAC encoding (DECIPHER threshold rule).
    // n=5, threshold=1.25. A=2>1.25, G=2>1.25, C=1≤1.25 dropped → set {A,G} → R (not a 3-base code).
    [Test]
    public void GenerateConsensus_MinorityBelowThreshold_DroppedYieldsR()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "A", "G", "G", "C" });
        Assert.That(consensus, Is.EqualTo("R"),
            "C (count 1 ≤ threshold 1.25) is below 25% and dropped; surviving {A,G} → R, not B/V.");
    }

    // M15 — No base passes the threshold → IUPAC code of the bases tied at the maximum count.
    // n=4, threshold=1.0; each base count=1, none >1.0; all four equally abundant → N.
    // Sources: DECIPHER ConsensusSequence ("degeneracy codes are always used in cases where multiple
    // characters are equally abundant"); Biopython 1.88 degenerate_consensus(["A","C","G","T"]) = "N".
    // Before F13 the method returned 'A' — a base with no more support than C, G or T.
    [Test]
    public void GenerateConsensus_FourEqualBases_ReturnsN()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "G", "T" });
        Assert.That(consensus, Is.EqualTo("N"),
            "Four equally abundant bases (none > 25%) are encoded as N, not an arbitrary 'A'.");
    }

    // M16 — Gap-diluted tie: n=4, threshold=1.0; A=1, C=1, two gaps → no base passes; A and C tie
    // at the maximum → {A,C} = M (was 'A' before F13). A single observed base → that base.
    [Test]
    public void GenerateConsensus_NoBasePasses_TiedBasesEncoded()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.GenerateConsensus(new[] { "A", "C", "-", "-" }), Is.EqualTo("M"),
                "A and C equally abundant, none > 25% → M.");
            Assert.That(MotifFinder.GenerateConsensus(new[] { "A", "-", "-", "-" }), Is.EqualTo("A"),
                "Only A observed → A.");
        });
    }

    // M17 — A column without any A/C/G/T (gaps / N) is unknown → N (was 'A' before F13).
    [Test]
    public void GenerateConsensus_ColumnWithoutAcgt_ReturnsN()
    {
        Assert.That(MotifFinder.GenerateConsensus(new[] { "A-", "AN", "A-" }), Is.EqualTo("AN"),
            "Column 1 has no A/C/G/T → N (any base, NC-IUB 1984), never a fabricated 'A'.");
    }

    #endregion

    #region GenerateConsensus — SHOULD (edge / invariants)

    // S1 — Case-insensitive: lowercase input gives the same result as upper-case (INV / normalisation).
    [Test]
    public void GenerateConsensus_LowercaseInput_SameAsUppercase()
    {
        var lower = MotifFinder.GenerateConsensus(new[] { "atgc", "gtgc" });
        Assert.That(lower, Is.EqualTo("RTGC"),
            "Input is upper-cased before counting; lowercase yields the same RTGC consensus.");
    }

    // S2 — Empty collection → "" (INV-06).
    [Test]
    public void GenerateConsensus_EmptyCollection_ReturnsEmpty()
    {
        var consensus = MotifFinder.GenerateConsensus(Array.Empty<string>());
        Assert.That(consensus, Is.Empty,
            "An empty collection has no columns to summarise → empty string.");
    }

    // S3 — Output length equals the first sequence's length (INV-01).
    [Test]
    public void GenerateConsensus_OutputLength_MatchesFirstSequence()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ACGTACG", "ACGTACG" });
        Assert.That(consensus, Has.Length.EqualTo(7),
            "One symbol per column; the column count equals the first sequence's length.");
    }

    #endregion

    #region GenerateConsensus — COULD (guards)

    // F13 — a null row threw NullReferenceException; unequal rows were silently truncated/padded.
    [Test]
    public void GenerateConsensus_NullElement_ThrowsArgumentException()
    {
        Assert.That(() => MotifFinder.GenerateConsensus(new[] { "ACGT", null!, "ACGT" }),
            NUnit.Framework.Throws.TypeOf<ArgumentException>());
    }

    // Biopython MultipleSeqAlignment / motifs.create reject rows of unequal length.
    [TestCase("ACG", "AC")]
    [TestCase("AC", "ACG")]
    public void GenerateConsensus_UnequalLengths_ThrowsArgumentException(string a, string b)
    {
        Assert.That(() => MotifFinder.GenerateConsensus(new[] { a, b }),
            NUnit.Framework.Throws.TypeOf<ArgumentException>());
    }

    // C1 — Null collection throws ArgumentNullException (documented guard).
    [Test]
    public void GenerateConsensus_Null_ThrowsArgumentNullException()
    {
        Assert.That(() => MotifFinder.GenerateConsensus(null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "A null sequence collection is rejected with ArgumentNullException.");
    }

    #endregion

    #region GenerateCavenerConsensus — Cavener 1987 / TRANSFAC / Biopython degenerate_consensus

    // Expected values computed with Biopython 1.88 Bio.motifs.create(rows).degenerate_consensus
    // (Cavener 1987, NAR 15(4):1353; rule text in Bio/motifs/matrix.py), not from this code.
    private static IEnumerable<TestCaseData> CavenerBiopythonCases()
    {
        // Biopython Tutorial "Sequence motif analysis" — m.degenerate_consensus = WACVC.
        yield return new TestCaseData(new[] { "TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC" }, "WACVC")
            .SetName("Cavener_BiopythonTutorial_WACVC");
        // Tutorial reverse complement — r.degenerate_consensus = GBGTW.
        yield return new TestCaseData(new[] { "TTGTA", "GCGTA", "GTGTA", "GGGTA", "GGGTT", "GCATT", "GCATT" }, "GBGTW")
            .SetName("Cavener_BiopythonTutorialReverseComplement_GBGTW");
        // Tutorial slice m[2:-1] — degenerate_consensus = CV.
        yield return new TestCaseData(new[] { "CA", "CG", "CA", "CC", "CC", "TG", "TG" }, "CV")
            .SetName("Cavener_BiopythonTutorialSlice_CV");
        // Rule branches (single / pair / triple / N).
        yield return new TestCaseData("AAAAAACCGT".Select(c => c.ToString()).ToArray(), "A").SetName("Cavener_Single_A6C2G1T1");
        yield return new TestCaseData(new[] { "A", "A", "A", "C" }, "A").SetName("Cavener_Single_A3C1");
        yield return new TestCaseData("AAAAACCCGT".Select(c => c.ToString()).ToArray(), "M").SetName("Cavener_Pair_A5C3_NotSingle");
        yield return new TestCaseData(new[] { "A", "A", "G", "G" }, "R").SetName("Cavener_Pair_A2G2");
        yield return new TestCaseData(new[] { "A", "C", "G" }, "V").SetName("Cavener_Triple_FourthAbsent");
        yield return new TestCaseData(new[] { "A", "A", "C", "G", "T" }, "N").SetName("Cavener_N_A2CGT");
        yield return new TestCaseData(new[] { "A", "C", "G", "T" }, "N").SetName("Cavener_N_FourEqual");
        yield return new TestCaseData(new[] { "AAAA", "AAGT", "AACT", "AATT" }, "AANT").SetName("Cavener_M13Rows_AANT");
        // Random alignments (seed 20260930), Biopython degenerate_consensus.
        yield return new TestCaseData(new[] { "TGTATTC", "TGAATAG", "AAATTTT", "AAATCGA", "TGAAAGG", "AAACCGG", "AAATAAT", "AAGATTA", "GTGAGTA" }, "WRAWNKN");
        yield return new TestCaseData(new[] { "ATTAGA", "ACTTGG" }, "AYTWGR");
        yield return new TestCaseData(new[] { "TAGCT", "CAAAT", "GAGGC", "TTTTT", "GGTAG", "TCTTC", "TGATC", "ATAAT" }, "NNDNY");
        yield return new TestCaseData(new[] { "GCGGGGGG", "CAGGGAGC" }, "SMGGGRGS");
        yield return new TestCaseData(new[] { "CCCCGCAGC", "CCCCCCCCC", "GCCCCCCAC", "CCCCCCCCC" }, "CCCCCCCVC");
        yield return new TestCaseData(new[] { "TCTTAA", "TTTTTT", "CTATTT", "CTTTTT", "TTATTA" }, "YTWTTW");
        yield return new TestCaseData(new[] { "TTCTATATT", "AATAATGTT", "TACAAGTAT", "TTCTATTAT" }, "TWCWATDWT");
        yield return new TestCaseData(new[] { "GGCTGAGGA", "TGGTCGGGG", "ACGGTTCAA", "CGAGTTAGG", "GCAGGCGAC", "GGCCGACAA" }, "NSVKKNSRR");
        yield return new TestCaseData(new[] { "CAGAACTGT", "CGTAGATTA", "ATGGCACAT", "CGCGAAATG", "CGACGCAGA", "ATGAGGTCC", "CCAAGACAG", "TGCGAGCGC", "TCTGCCTAC" }, "CNNRRMYNN");
        yield return new TestCaseData(new[] { "AAGGCAGTA", "GTGGGCGAA", "AAGATGAGA", "AGGAGGAGC", "GACAGAAGT", "AGCCAAACA", "ACGCTAGGA", "AGCAGAGGA", "GACCAGGGA" }, "RRSMNRRGA");
        yield return new TestCaseData(new[] { "CGGT", "GCTG", "GTGC" }, "SBKB");
        yield return new TestCaseData(new[] { "CCTTCCCC", "TCTTTCCC", "TCCTTTCC" }, "YCYTYYCC");
    }

    [TestCaseSource(nameof(CavenerBiopythonCases))]
    public void GenerateCavenerConsensus_EqualsBiopythonDegenerateConsensus(string[] rows, string expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.GenerateCavenerConsensus(rows), Is.EqualTo(expected));
            Assert.That(MotifFinder.GenerateCavenerConsensus(rows.Select(r => r.ToLowerInvariant())), Is.EqualTo(expected),
                "case-insensitive");
        });
    }

    [Test]
    public void GenerateCavenerConsensus_Empty_ReturnsEmpty()
    {
        Assert.That(MotifFinder.GenerateCavenerConsensus(Array.Empty<string>()), Is.Empty);
    }

    [Test]
    public void GenerateCavenerConsensus_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => MotifFinder.GenerateCavenerConsensus(null!),
                NUnit.Framework.Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => MotifFinder.GenerateCavenerConsensus(new[] { "ACGT", null! }),
                NUnit.Framework.Throws.TypeOf<ArgumentException>());
            Assert.That(() => MotifFinder.GenerateCavenerConsensus(new[] { "ACGT", "ACG" }),
                NUnit.Framework.Throws.TypeOf<ArgumentException>());
            Assert.That(() => MotifFinder.GenerateCavenerConsensus(new[] { "AC-T", "ACGT" }),
                NUnit.Framework.Throws.TypeOf<ArgumentException>(), "gaps are not A/C/G/T");
        });
    }

    #endregion
}
