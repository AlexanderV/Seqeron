// ONCO-CNA-003 — Homozygous (Deep) Deletion Detection
// Evidence: docs/Evidence/ONCO-CNA-003-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CNA-003.md
// Source: Cheng J et al. (2017). Pan-cancer homozygous deletions. Nat Commun 8:1221.
//         https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/  (homozygous = zero copies of both alleles, total CN 0)
//         cBioPortal discrete CNA: -2 Deep Deletion = homozygous; -1 shallow = heterozygous.
//         https://docs.cbioportal.org/file-formats/  https://docs.cbioportal.org/user-guide/faq/
//         CNVkit absolute_threshold integer-CN (defaults -1.1,-0.25,0.2,0.7); NCBI Gene arms.

using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberArmSegment;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_DetectHomozygousDeletions_Tests
{
    // Fixed arm length; arm fraction is irrelevant to deletion calling (only log2/CN matter).
    private const long Arm = 1_000_000;

    private static Segment Seg(string arm, double log2) => new(arm, 0, 1_000, Arm, log2);

    #region DetectHomozygousDeletions

    // M1 — log2 = -2.0 ⇒ CN 0 (DeepDeletion) ⇒ homozygous deletion.
    // Source: Cheng 2017 (total CN 0); cBioPortal -2; CNVkit (<= -1.1 ⇒ CN 0).
    [Test]
    public void DetectHomozygousDeletions_DeepDeletionSegment_IsReported()
    {
        var segments = new[] { Seg("9p", -2.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "log2 = -2.0 classifies to integer CN 0 (total CN 0 = both alleles lost), so it is a homozygous deletion.");
            Assert.That(result[0].Arm, Is.EqualTo("9p"),
                "The reported segment must be the input CN-0 segment.");
        });
    }

    // M2 — log2 = -0.5 ⇒ CN 1 (single-copy / heterozygous loss) ⇒ NOT homozygous.
    // Source: cBioPortal -1 = shallow/heterozygous (not -2); Cheng 2017 (one allele remains).
    [Test]
    public void DetectHomozygousDeletions_SingleCopyLoss_NotReported()
    {
        var segments = new[] { Seg("3p", -0.5) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "log2 = -0.5 is integer CN 1 (heterozygous single-copy loss), which is NOT a homozygous deletion.");
    }

    // M3 — log2 = 0.0 ⇒ CN 2 (diploid) ⇒ NOT reported.
    // Source: cBioPortal 0 = diploid.
    [Test]
    public void DetectHomozygousDeletions_NeutralSegment_NotReported()
    {
        var segments = new[] { Seg("10q", 0.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "A copy-number-neutral segment (CN 2) must not be reported as a homozygous deletion.");
    }

    // M4 — gain (log2 0.5 ⇒ CN 3) and amplification (log2 1.0 ⇒ CN >=4) ⇒ NOT reported.
    // Source: cBioPortal 1 = gain, 2 = amplification (neither is a deletion).
    [Test]
    public void DetectHomozygousDeletions_GainAndAmplification_NotReported()
    {
        var segments = new[] { Seg("8q", 0.5), Seg("17q", 1.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "Gain (CN 3) and amplification (CN >=4) are increases, not homozygous deletions.");
    }

    // M5 — mixed [CN0, CN1, CN0, CN2] ⇒ only the two CN-0 segments, in input order.
    // Source: INV-1 (CN 0 only); INV-3 (order-preserving filter).
    [Test]
    public void DetectHomozygousDeletions_MixedSet_ReturnsCn0InOrder()
    {
        var segments = new[]
        {
            Seg("9p", -2.0),  // CN 0
            Seg("3p", -0.5),  // CN 1
            Seg("10q", -1.5), // CN 0
            Seg("1q", 0.0),   // CN 2
        };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result.Select(s => s.Arm), Is.EqualTo(new[] { "9p", "10q" }),
                "Only the two CN-0 segments are reported, preserving input order.");
            Assert.That(result, Has.Count.EqualTo(2),
                "Exactly the two homozygous-deletion segments are reported.");
        });
    }

    // M6 — log2 exactly -1.1 (the deletion cutoff) ⇒ CN 0 (<= cutoff) ⇒ reported.
    // Source: CNVkit assigns CN by "less than or equal to each threshold in sequence".
    [Test]
    public void DetectHomozygousDeletions_Log2AtDeletionCutoff_IsReported()
    {
        var segments = new[] { Seg("9p", -1.1) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Has.Count.EqualTo(1),
            "log2 exactly at the deletion cutoff -1.1 is <= the cutoff, so CN 0 (homozygous).");
    }

    // M7 — log2 just above the cutoff (-1.0999) ⇒ CN 1 ⇒ NOT reported.
    // Source: CNVkit threshold boundary (> -1.1 ⇒ CN 1).
    [Test]
    public void DetectHomozygousDeletions_Log2JustAboveCutoff_NotReported()
    {
        var segments = new[] { Seg("9p", -1.0999) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "log2 = -1.0999 is above the -1.1 cutoff, so CN 1 (heterozygous), not homozygous.");
    }

    // M14 — custom thresholds: raising the deletion cutoff to -0.4 makes log2 = -0.5 ⇒ CN 0.
    // Source: CNVkit thresholds are parameters of the calling.
    [Test]
    public void DetectHomozygousDeletions_CustomThresholds_ShiftsCn0Boundary()
    {
        var segments = new[] { Seg("10q", -0.5) };
        var custom = new[] { -0.4, -0.25, 0.2, 0.7 }; // deletion cutoff now -0.4

        var withDefault = OncologyAnalyzer.DetectHomozygousDeletions(segments);
        var withCustom = OncologyAnalyzer.DetectHomozygousDeletions(segments, custom);

        Assert.Multiple(() =>
        {
            Assert.That(withDefault, Is.Empty,
                "Under default thresholds log2 = -0.5 is CN 1 (not homozygous).");
            Assert.That(withCustom, Has.Count.EqualTo(1),
                "With deletion cutoff -0.4, log2 = -0.5 (<= -0.4) is CN 0 (homozygous).");
        });
    }

    // S3 — ploidy parameter feeds CNVkit n = ploidy*2^log2; a triploid no-call/neutral is not homozygous.
    // Source: CNVkit n = ploidy*2^log2; NaN no-call returns neutral reference CN (rounded ploidy).
    [Test]
    public void DetectHomozygousDeletions_TriploidPloidy_RespectsBoundary()
    {
        // log2 = -2.0 is <= -1.1 (deletion cutoff is on log2, not ploidy-scaled) ⇒ CN 0 regardless of ploidy.
        var deep = new[] { Seg("9p", -2.0) };
        // NaN under triploid is a neutral no-call (CN = round(3) = 3), not homozygous.
        var noCall = new[] { Seg("9p", double.NaN) };

        var deepResult = OncologyAnalyzer.DetectHomozygousDeletions(deep, thresholds: null, ploidy: 3.0);
        var noCallResult = OncologyAnalyzer.DetectHomozygousDeletions(noCall, thresholds: null, ploidy: 3.0);

        Assert.Multiple(() =>
        {
            Assert.That(deepResult, Has.Count.EqualTo(1),
                "log2 <= -1.1 is CN 0 by the log2 threshold, independent of ploidy.");
            Assert.That(noCallResult, Is.Empty,
                "A NaN log2 is a neutral no-call (CN = rounded ploidy), not a homozygous deletion.");
        });
    }

    // C1 — empty input ⇒ empty result.
    [Test]
    public void DetectHomozygousDeletions_EmptyInput_ReturnsEmpty()
    {
        var result = OncologyAnalyzer.DetectHomozygousDeletions(Array.Empty<Segment>());

        Assert.That(result, Is.Empty, "An empty segment set yields no homozygous deletions.");
    }

    // C2 — null input ⇒ ArgumentNullException.
    [Test]
    public void DetectHomozygousDeletions_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(null!),
            "Null segments must throw ArgumentNullException.");
    }

    // C3 — invalid segment (End <= Start) ⇒ ArgumentException.
    [Test]
    public void DetectHomozygousDeletions_InvalidSegment_Throws()
    {
        var bad = new[] { new Segment("9p", 1_000, 1_000, Arm, -2.0) }; // End == Start

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(bad),
            "A segment with End <= Start must throw ArgumentException.");
    }

    // C3b — invalid segment (non-positive ArmLength) ⇒ ArgumentException.
    // Source: contract §3.3 — "A segment with non-positive ArmLength ... → ArgumentException" (ValidateArmSegment).
    [Test]
    public void DetectHomozygousDeletions_NonPositiveArmLength_Throws()
    {
        var bad = new[] { new Segment("9p", 0, 1_000, 0, -2.0) }; // ArmLength == 0

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(bad),
            "A segment with non-positive ArmLength must throw ArgumentException.");
    }

    // C6 — non-positive ploidy ⇒ ArgumentOutOfRangeException (via CallCopyNumber).
    // Source: contract §3.3 — "Non-positive ploidy → ArgumentOutOfRangeException".
    [Test]
    public void DetectHomozygousDeletions_NonPositivePloidy_Throws()
    {
        var segments = new[] { Seg("9p", -2.0) };

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(segments, thresholds: null, ploidy: 0.0),
            "Ploidy must be positive; ploidy = 0 must throw ArgumentOutOfRangeException.");
    }

    // C7 — thresholds not four strictly-ascending values ⇒ ArgumentException (via CallCopyNumber).
    // Source: contract §3.3 — "thresholds not four strictly ascending values → ArgumentException".
    [Test]
    public void DetectHomozygousDeletions_InvalidThresholds_Throws()
    {
        var segments = new[] { Seg("9p", -2.0) };
        var notAscending = new[] { -1.1, -1.1, 0.2, 0.7 }; // not strictly ascending
        var wrongCount = new[] { -1.1, -0.25, 0.2 };        // only three thresholds

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.DetectHomozygousDeletions(segments, notAscending),
                "Thresholds must be strictly ascending.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.DetectHomozygousDeletions(segments, wrongCount),
                "Exactly four thresholds are required.");
        });
    }

    // C5 — NaN log2 ⇒ neutral no-call ⇒ NOT reported.
    // Source: CNVkit NaN log2 is a no-call returning the neutral reference CN.
    [Test]
    public void DetectHomozygousDeletions_NaNLog2_NotReported()
    {
        var segments = new[] { Seg("9p", double.NaN) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "A NaN log2 is a neutral no-call (CN = ploidy), not a homozygous deletion.");
    }

    #endregion

    #region IsHomozygousDeletion

    // M8 — predicate: CN-0 segment true, CN-1 segment false.
    // Source: INV-1 (CN 0 = homozygous); INV-2 (CN 1 heterozygous, not).
    [Test]
    public void IsHomozygousDeletion_Cn0AndCn1_TrueThenFalse()
    {
        var cn0 = Seg("9p", -2.0);
        var cn1 = Seg("3p", -0.5);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(cn0), Is.True,
                "log2 = -2.0 (CN 0) is a homozygous deletion.");
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(cn1), Is.False,
                "log2 = -0.5 (CN 1) is a heterozygous loss, not homozygous.");
        });
    }

    // M8b — IsHomozygousDeletion at the inclusive boundary: log2 = -1.1 (<= cutoff) is true,
    // log2 just above (-1.0999) is false. Source: CNVkit "<= each threshold" boundary.
    [Test]
    public void IsHomozygousDeletion_DeletionCutoffBoundary_TrueThenFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(Seg("9p", -1.1)), Is.True,
                "log2 exactly -1.1 is <= the deletion cutoff, so CN 0 (homozygous).");
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(Seg("9p", -1.0999)), Is.False,
                "log2 = -1.0999 is above -1.1, so CN 1 (not homozygous).");
        });
    }

    // M8c — IsHomozygousDeletion validates the segment (End <= Start) ⇒ ArgumentException.
    // Source: contract §3.3 — IsHomozygousDeletion calls ValidateArmSegment.
    [Test]
    public void IsHomozygousDeletion_InvalidSegment_Throws()
    {
        var bad = new Segment("9p", 1_000, 1_000, Arm, -2.0); // End == Start

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.IsHomozygousDeletion(bad),
            "IsHomozygousDeletion must validate the segment (End <= Start ⇒ ArgumentException).");
    }

    // F1 (review 2026-09) — extreme positive log2 lies above the last CNVkit cutoff, so its copy number is
    // ceil(ploidy·2^log2) >= 4 (amplification), never CN 0. CNVkit absolute_threshold: log2 = 40 ⇒ CN
    // 2,199,023,255,552 (beyond Int32, so CallCopyNumber rejects it — ONCO-CNA-001 guard); log2 = +∞ ⇒ CNVkit
    // OverflowError. The CN-0 predicate must still answer false instead of throwing.
    // Source: CNVkit cnvlib/call.py absolute_threshold ("Above the last threshold value ... rounding up").
    [Test]
    public void IsHomozygousDeletion_Log2BeyondInt32CopyNumber_IsFalse_WhileCallCopyNumberSaturates()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(Seg("9p", 40.0)), Is.False,
                "log2 = 40 ⇒ CNVkit CN = ceil(2·2^40) = 2199023255552 (amplification), not CN 0.");
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(Seg("9p", double.PositiveInfinity)), Is.False,
                "log2 = +∞ is above every cutoff (unbounded amplification), not CN 0.");
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(Seg("9p", double.MaxValue)), Is.False,
                "log2 = double.MaxValue is above every cutoff, not CN 0.");
            Assert.That(OncologyAnalyzer.CallCopyNumber(40.0), Is.EqualTo(int.MaxValue),
                "ONCO-CNA-001 contract: an integer CN beyond Int32 saturates at Int32.MaxValue (Amplification).");
            Assert.That(OncologyAnalyzer.CallCopyNumber(double.PositiveInfinity), Is.EqualTo(int.MaxValue),
                "ONCO-CNA-001 contract: +∞ saturates at Int32.MaxValue (CNVkit raises OverflowError).");
        });
    }

    // F1 — the filter must not crash on a stream containing extreme amplifications; −∞ (zero copies) is CN 0.
    [Test]
    public void DetectHomozygousDeletions_StreamWithExtremeAmplifications_ReportsOnlyCn0()
    {
        var segments = new[]
        {
            Seg("17p", double.PositiveInfinity), // unbounded amplification
            Seg("9p", -2.0),                     // CN 0
            Seg("8q", 40.0),                     // CN 2.2e12
            Seg("10q", double.NegativeInfinity), // 2^-∞ = 0 copies ⇒ CN 0
        };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result.Select(s => s.Arm), Is.EqualTo(new[] { "9p", "10q" }),
            "Only the CN-0 segments are reported, in input order; extreme amplifications are excluded without throwing.");
    }

    // F2 (review 2026-09) — calling parameters are validated eagerly, even for an empty segment list
    // (consistent with DetectFocalAmplifications / ClassifyCopyNumbers).
    [Test]
    public void DetectHomozygousDeletions_EmptyInputWithInvalidParameters_Throws()
    {
        var empty = Array.Empty<Segment>();

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.DetectHomozygousDeletions(empty, new[] { -1.1, -0.25, 0.2 }),
                "Three thresholds are invalid even when there are no segments.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DetectHomozygousDeletions(empty, thresholds: null, ploidy: double.NaN),
                "A NaN ploidy is invalid even when there are no segments.");
        });
    }

    #endregion

    #region IdentifyDeletedTumorSuppressors

    // M9 — 17p ⇒ TP53. Source: NCBI Gene TP53 17p13.1.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm17p_MapsTp53()
    {
        var dels = new[] { Seg("17p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "TP53" }),
            "A homozygous deletion on 17p maps to TP53 (NCBI Gene 17p13.1).");
    }

    // M10 — 13q ⇒ RB1 and BRCA2 (both on 13q), in panel order RB1 before BRCA2.
    // Source: NCBI Gene RB1 13q14.2, BRCA2 13q13.1.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm13q_MapsRb1AndBrca2()
    {
        var dels = new[] { Seg("13q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "RB1", "BRCA2" }),
            "13q carries both RB1 (13q14.2) and BRCA2 (13q13.1); both reported in panel order.");
    }

    // M11 — 9p ⇒ CDKN2A. Source: NCBI Gene CDKN2A 9p21.3.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm9p_MapsCdkn2a()
    {
        var dels = new[] { Seg("9p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "CDKN2A" }),
            "A homozygous deletion on 9p maps to CDKN2A (NCBI Gene 9p21.3).");
    }

    // M12 — 10q ⇒ PTEN. Source: NCBI Gene PTEN 10q23.31.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm10q_MapsPten()
    {
        var dels = new[] { Seg("10q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "PTEN" }),
            "A homozygous deletion on 10q maps to PTEN (NCBI Gene 10q23.31).");
    }

    // M13 — 17q ⇒ BRCA1. Source: NCBI Gene BRCA1 17q21.31.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm17q_MapsBrca1()
    {
        var dels = new[] { Seg("17q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "BRCA1" }),
            "A homozygous deletion on 17q maps to BRCA1 (NCBI Gene 17q21.31).");
    }

    // S1 — non-panel arm (1p) ⇒ no gene.
    [Test]
    public void IdentifyDeletedTumorSuppressors_NonPanelArm_ReturnsEmpty()
    {
        var dels = new[] { Seg("1p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.Empty,
            "An arm with no panel tumour suppressor (1p) maps to no gene.");
    }

    // S2 — two deletions on 13q ⇒ RB1 and BRCA2 each once.
    [Test]
    public void IdentifyDeletedTumorSuppressors_DuplicateArm_GenesOnce()
    {
        var dels = new[] { Seg("13q", -2.0), Seg("13q", -1.5) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "RB1", "BRCA2" }),
            "Repeated deletions on the same arm report each gene once.");
    }

    // C4 — null ⇒ ArgumentNullException.
    [Test]
    public void IdentifyDeletedTumorSuppressors_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.IdentifyDeletedTumorSuppressors(null!),
            "Null deletions must throw ArgumentNullException.");
    }

    #endregion

    #region IdentifyDeletedTumorSuppressors — locus overlap (F49, GISTIC2 genes_at partial_hits = 1)

    // Locked from GNU Octave running GISTIC2 genes_at.m on the GRCh38 panel (GENCODE v22, GISTIC2 refgenes).
    [TestCase("17", 7_000_000, 50_000_000, new[] { "TP53", "BRCA1" })]
    [TestCase("13", 32_315_474, 32_315_474, new[] { "BRCA2" })]   // single base = BRCA2 start
    [TestCase("9", 1, 21_967_752, new string[0])]                  // ends one base before CDKN2A
    [TestCase("9", 1, 21_967_753, new[] { "CDKN2A" })]
    public void IdentifyDeletedTumorSuppressors_LocusOverlap_MatchesGistic2GenesAt(
        string chromosome, long start, long end, string[] expected)
    {
        var result = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(
            new[] { new OncologyAnalyzer.CopyNumberRegion(chromosome, start, end) });
        Assert.That(result, Is.EqualTo(expected),
            "Genes are reported iff their locus overlaps a deleted region (closed intervals, panel order).");
    }

    [Test]
    public void DefaultTumorSuppressorLoci_Gistic2Gencode22Coordinates()
    {
        Assert.That(OncologyAnalyzer.DefaultTumorSuppressorLoci, Is.EqualTo(new[]
        {
            new OncologyAnalyzer.GeneLocus("TP53", "17", 7_661_779, 7_687_550),
            new OncologyAnalyzer.GeneLocus("RB1", "13", 48_303_751, 48_481_986),
            new OncologyAnalyzer.GeneLocus("CDKN2A", "9", 21_967_753, 21_995_301),
            new OncologyAnalyzer.GeneLocus("PTEN", "10", 87_863_113, 87_971_930),
            new OncologyAnalyzer.GeneLocus("BRCA1", "17", 43_044_295, 43_125_483),
            new OncologyAnalyzer.GeneLocus("BRCA2", "13", 32_315_474, 32_400_266),
        }));
    }

    #endregion
}
