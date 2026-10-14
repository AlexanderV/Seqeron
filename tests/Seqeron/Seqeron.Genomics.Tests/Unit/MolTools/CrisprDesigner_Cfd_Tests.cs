namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the Doench et al. 2016 CFD (Cutting Frequency Determination) off-target score
/// (<see cref="CrisprDesigner.CalculateCfdScore"/>), unit CRISPR-OFF-001 (C7 residual).
///
/// Source of the model + matrices:
///   Doench, Fusi, Sullender, Hegde, et al. Nat Biotechnol 34:184-191 (2016), PMID 26780180.
///   The mismatch percent-activity matrix (240 entries = 12 mismatch types x 20 positions) and the
///   PAM-activity table (16 NGG-region dinucleotides) are shipped as the authoritative pickles
///   `mismatch_score.pkl` / `pam_scores.pkl` in the canonical `cfd-score-calculator.py`
///   (John Doench's lab; redistributed by CRISPOR maximilianh/crisporWebsite and bm2-lab/iGWOS).
///   The matrices were decoded to text and cross-checked element-by-element across BOTH repos
///   (240/240 + 16/16 identical) before being transcribed into the implementation.
///
/// Every expected value below was re-derived by an INDEPENDENT Python run of the decoded matrices
/// using the published `calc_cfd` formula (NOT read off this C# code). The first oracle is the
/// published iGWOS doctest value, tying the implementation to the reference's own documented output.
///
/// CFD algorithm (verbatim): score = product over the 20 protospacer positions of the per-position
/// mismatch penalty (1.0 where guide:off-target match) x the PAM-activity score; key for a mismatch
/// at index i is 'r'+guideBase(T->U)+':d'+complement(offTargetBase)+','+(i+1). Position 1 (index 0)
/// is the 5'/PAM-DISTAL end; position 20 is the 3'/PAM-PROXIMAL (seed) end. Perfect match + GG -> 1.0.
/// </summary>
[TestFixture]
public class CrisprDesigner_Cfd_Tests
{
    // The published/decoded matrix values are quoted to ~9-15 sig figs; C# double matches well
    // within this tolerance. Where an expected value is an exact matrix entry, the assertion
    // recomputes nothing (the constant is the EXPECTED, the C# result is the ACTUAL).
    private const double Tol = 1e-9;

    // A representative real-style guide used throughout (5' -> 3', PAM-distal first).
    private const string Guide = "GACGCATAAAGATGAGACGC";

    private static string Mutate(string seq, int position1Based, char newBase)
    {
        var chars = seq.ToCharArray();
        chars[position1Based - 1] = newBase;
        return new string(chars);
    }

    // -------------------------------------------------------------------------------------------
    // Perfect match -> 1.0
    // -------------------------------------------------------------------------------------------

    /// <summary>Perfect match against a canonical NGG PAM -> exactly 1.0 (Doench 2016 definition).</summary>
    [Test]
    public void Cfd_PerfectMatch_CanonicalPam_ReturnsExactlyOne()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "GG");
        Assert.That(score, Is.EqualTo(1.0));
    }

    /// <summary>Perfect match with a 3-nt PAM (NGG) -> 1.0 (only the last two PAM nt are scored).</summary>
    [Test]
    public void Cfd_PerfectMatch_ThreeNtPam_ReturnsExactlyOne()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "AGG");
        Assert.That(score, Is.EqualTo(1.0));
    }

    /// <summary>
    /// Published reference oracle (iGWOS `calcCfdScore` doctest):
    /// guide GGGGGGGGGGGGGGGGGGGG vs off GGGGGGGGGGGGGGGGGAAA + GG PAM -> 0.4635989007074176.
    /// This is the load-bearing cross-check tying this implementation to the reference's documented
    /// output (= rG:dT,18 x rG:dT,19 x rG:dT,20).
    /// </summary>
    [Test]
    public void Cfd_PublishedReferenceOracle_TripleMismatch_MatchesDoctest()
    {
        double score = CrisprDesigner.CalculateCfdScore("GGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAA", "GG");
        Assert.That(score, Is.EqualTo(0.4635989007074176).Within(Tol));
    }

    /// <summary>
    /// Published reference oracle (iGWOS `calcCfdScore` doctest, mixed-case input):
    /// guide GGGGGGGGGGGGGGGGGGGGGGG vs off "aaaaGaGaGGGGGGGGGGGGGGG" -> 0.5140384614450001.
    /// Exercises lowercase handling AND the reference's 23-mer-with-PAM contract (off[:-3] / off[-2:]).
    /// </summary>
    [Test]
    public void Cfd_PublishedReferenceOracle_MixedCase_MatchesDoctest()
    {
        // off (upper) = AAAAGAGAGGGGGGGGGGGGGGG ; PAM = GG ; protospacer = AAAAGAGAGGGGGGGGGGGG (20 nt).
        double score = CrisprDesigner.CalculateCfdScore("GGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGG", "gg");
        Assert.That(score, Is.EqualTo(0.5140384614450001).Within(Tol));
    }

    // -------------------------------------------------------------------------------------------
    // Single-mismatch -> exact matrix entry (each tied to the decoded source value)
    // -------------------------------------------------------------------------------------------

    /// <summary>
    /// Single mismatch at position 1 (PAM-distal). Guide base G, off base A -> key rG:dT,1 = 0.9.
    /// With a canonical GG PAM the whole CFD equals exactly that matrix entry.
    /// </summary>
    [Test]
    public void Cfd_SingleMismatch_Position1_EqualsMatrixEntry_rGdT1()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 1, 'A'), "GG");
        Assert.That(score, Is.EqualTo(0.9).Within(Tol)); // rG:dT,1
    }

    /// <summary>
    /// Single mismatch at position 5. Guide base C, off base A -> key rC:dT,5 = 0.571428571.
    /// </summary>
    [Test]
    public void Cfd_SingleMismatch_Position5_EqualsMatrixEntry_rCdT5()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 5, 'A'), "GG");
        Assert.That(score, Is.EqualTo(0.571428571).Within(Tol)); // rC:dT,5
    }

    /// <summary>
    /// Single mismatch at position 7. Guide base T (=rU), off base C -> complement G -> rU:dG,7 = 0.6875.
    /// Covers the rU (guide-T-as-U) key path AND the off-base complementation (off C -> dG, not dC).
    /// </summary>
    [Test]
    public void Cfd_SingleMismatch_GuideTAsU_EqualsMatrixEntry_rUdG7()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 7, 'C'), "GG");
        Assert.That(score, Is.EqualTo(0.6875).Within(Tol)); // rU:dG,7
    }

    /// <summary>
    /// Single mismatch at position 16, a NEAR-ZERO penalty: guide base G, off base T -> complement A
    /// -> rG:dA,16 = 0.0, so the whole CFD collapses to 0. Guards a high-impact (seed-ish) entry.
    /// </summary>
    [Test]
    public void Cfd_SingleMismatch_ZeroPenaltyEntry_CollapsesToZero_rGdA16()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 16, 'T'), "GG");
        Assert.That(score, Is.EqualTo(0.0).Within(Tol)); // rG:dA,16
    }

    /// <summary>
    /// Single mismatch at position 20 (PAM-proximal). Guide base C, off base A -> rC:dT,20 = 0.5.
    /// </summary>
    [Test]
    public void Cfd_SingleMismatch_Position20_EqualsMatrixEntry_rCdT20()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 20, 'A'), "GG");
        Assert.That(score, Is.EqualTo(0.5).Within(Tol)); // rC:dT,20
    }

    // -------------------------------------------------------------------------------------------
    // PAM-score application (canonical vs non-canonical)
    // -------------------------------------------------------------------------------------------

    /// <summary>Perfect match, GA PAM -> applies pam_scores["GA"] = 0.069444444 (non-canonical).</summary>
    [Test]
    public void Cfd_PerfectMatch_GaPam_AppliesPamScore()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "GA");
        Assert.That(score, Is.EqualTo(0.06944444400000001).Within(Tol));
    }

    /// <summary>Perfect match, AG PAM -> applies pam_scores["AG"] = 0.259259259.</summary>
    [Test]
    public void Cfd_PerfectMatch_AgPam_AppliesPamScore()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "AG");
        Assert.That(score, Is.EqualTo(0.25925925899999996).Within(Tol));
    }

    /// <summary>Perfect match, TG PAM -> applies pam_scores["TG"] = 0.038961039.</summary>
    [Test]
    public void Cfd_PerfectMatch_TgPam_AppliesPamScore()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "TG");
        Assert.That(score, Is.EqualTo(0.038961038999999996).Within(Tol));
    }

    /// <summary>Perfect match, AA PAM -> pam_scores["AA"] = 0.0, so the entire CFD is 0.</summary>
    [Test]
    public void Cfd_PerfectMatch_AaPam_ZeroPamCollapsesToZero()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Guide, "AA");
        Assert.That(score, Is.EqualTo(0.0));
    }

    // -------------------------------------------------------------------------------------------
    // Product-of-penalties (multi-mismatch) and combined mismatch x PAM
    // -------------------------------------------------------------------------------------------

    /// <summary>
    /// Two mismatches: position 1 (rG:dT,1 = 0.9) and position 20 (rC:dT,20 = 0.5). CFD is the PRODUCT
    /// of the per-position penalties: 0.9 x 0.5 = 0.45. Verifies multiplicative combination.
    /// </summary>
    [Test]
    public void Cfd_TwoMismatches_ProductOfPenalties_Equals045()
    {
        string off = Mutate(Mutate(Guide, 1, 'A'), 20, 'A');
        double score = CrisprDesigner.CalculateCfdScore(Guide, off, "GG");
        Assert.That(score, Is.EqualTo(0.9 * 0.5).Within(Tol)); // rG:dT,1 x rC:dT,20
    }

    /// <summary>
    /// Single mismatch (rG:dT,1 = 0.9) combined with a non-canonical GA PAM (0.069444444):
    /// CFD = 0.9 x 0.06944444400000001. Verifies mismatch penalty AND PAM score multiply together.
    /// </summary>
    [Test]
    public void Cfd_MismatchTimesNonCanonicalPam_MultipliesBoth()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 1, 'A'), "GA");
        Assert.That(score, Is.EqualTo(0.9 * 0.06944444400000001).Within(Tol));
    }

    // -------------------------------------------------------------------------------------------
    // ORIENTATION GUARD — fails if the position axis is reversed (the classic CFD bug)
    // -------------------------------------------------------------------------------------------

    /// <summary>
    /// Orientation guard. Using guide = C(x20), an off-target with a single A at position 1 yields
    /// rC:dT,1 = 1.0, while the SAME mismatch type at position 20 yields rC:dT,20 = 0.5. Because the
    /// two values differ AND the mismatch type is identical, this asserts that position 1 maps to
    /// matrix index 0 (PAM-distal) and position 20 to index 19 (PAM-proximal). If the position axis
    /// were reversed, the two expected values would swap (pos1 -> 0.5, pos20 -> 1.0) and BOTH
    /// assertions would fail. This is a genuine reversal-detector, not a code echo.
    /// </summary>
    [Test]
    public void Cfd_OrientationGuard_Position1VsPosition20_NotReversed()
    {
        const string guideC = "CCCCCCCCCCCCCCCCCCCC";

        double atPos1 = CrisprDesigner.CalculateCfdScore(guideC, Mutate(guideC, 1, 'A'), "GG");
        double atPos20 = CrisprDesigner.CalculateCfdScore(guideC, Mutate(guideC, 20, 'A'), "GG");

        Assert.That(atPos1, Is.EqualTo(1.0).Within(Tol), "rC:dT,1 (PAM-distal) must be 1.0");
        Assert.That(atPos20, Is.EqualTo(0.5).Within(Tol), "rC:dT,20 (PAM-proximal) must be 0.5");
        Assert.That(atPos1, Is.GreaterThan(atPos20), "position 1 must be more tolerated than position 20 for this mismatch");
    }

    // -------------------------------------------------------------------------------------------
    // Invariants / boundaries
    // -------------------------------------------------------------------------------------------

    /// <summary>CFD is always within [0, 1] for valid inputs (definition: product of values in [0,1]).</summary>
    [Test]
    public void Cfd_AnyValidInput_ScoreInUnitInterval()
    {
        string off = Mutate(Mutate(Mutate(Guide, 3, 'A'), 11, 'T'), 18, 'A');
        double score = CrisprDesigner.CalculateCfdScore(off /* swapped roles still valid 20-nt */, Guide, "GA");
        Assert.That(score, Is.GreaterThanOrEqualTo(0.0));
        Assert.That(score, Is.LessThanOrEqualTo(1.0));
    }

    /// <summary>Lowercase input is accepted (case-insensitive) and yields the same perfect-match 1.0.</summary>
    [Test]
    public void Cfd_LowercaseInput_TreatedAsUppercase()
    {
        double score = CrisprDesigner.CalculateCfdScore(Guide.ToLowerInvariant(), Guide.ToLowerInvariant(), "gg");
        Assert.That(score, Is.EqualTo(1.0));
    }

    /// <summary>Determinism: identical inputs produce identical outputs.</summary>
    [Test]
    public void Cfd_SameInput_Deterministic()
    {
        string off = Mutate(Guide, 8, 'C');
        double a = CrisprDesigner.CalculateCfdScore(Guide, off, "GG");
        double b = CrisprDesigner.CalculateCfdScore(Guide, off, "GG");
        Assert.That(a, Is.EqualTo(b));
    }

    // -------------------------------------------------------------------------------------------
    // Edge / error cases (documented contract)
    // -------------------------------------------------------------------------------------------

    [Test]
    public void Cfd_NullGuide_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateCfdScore(null!, Guide, "GG"));
    }

    [Test]
    public void Cfd_NullOffTarget_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateCfdScore(Guide, null!, "GG"));
    }

    [Test]
    public void Cfd_NullPam_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateCfdScore(Guide, Guide, null!));
    }

    [Test]
    public void Cfd_EmptyGuide_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateCfdScore("", Guide, "GG"));
    }

    [TestCase("ACGTACGTACGTACGTACG")]    // 19 nt
    [TestCase("ACGTACGTACGTACGTACGTA")]  // 21 nt
    public void Cfd_WrongLengthGuide_ThrowsArgumentException(string badGuide)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(badGuide, Guide, "GG"));
    }

    [TestCase("ACGTACGTACGTACGTACG")]    // 19 nt
    [TestCase("ACGTACGTACGTACGTACGTA")]  // 21 nt
    public void Cfd_WrongLengthOffTarget_ThrowsArgumentException(string badOff)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Guide, badOff, "GG"));
    }

    [TestCase("G")]      // 1 nt
    [TestCase("AGGG")]   // 4 nt
    public void Cfd_WrongLengthPam_ThrowsArgumentException(string badPam)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Guide, Guide, badPam));
    }

    [Test]
    public void Cfd_NonAcgtInGuide_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Mutate(Guide, 4, 'N'), Guide, "GG"));
    }

    [Test]
    public void Cfd_NonAcgtInOffTarget_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Guide, Mutate(Guide, 4, 'N'), "GG"));
    }

    [Test]
    public void Cfd_NonAcgtInPam_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Guide, Guide, "NG"));
    }

    // -------------------------------------------------------------------------------------------
    // Guide-level CFD aggregation (CRISPOR annotateOfftargets, crispor.py L1994-L1999)
    //
    //   default              guideCfdScore = calcMitGuideScore(sum(cfdScores))  = 100/(100+Σ)*100
    //   doCfdFix (fixCfd=1)  guideCfdScore = calcCfdGuideScore(sum(cfdScores))  = 100*100/Σ
    //
    // Under doCfdFix the reference ALSO rescales every per-hit CFD by x100 inside calcCfdScore
    // (crispor.py L2180-2182), so over a sum of RAW per-hit CFD values the two factors partly cancel:
    //   100*100 / (100 * Σ raw) = 100 / Σ raw.
    // This API always takes RAW CFD values in [0, 1], so ParkinsonNormalized must be 100/Σ raw —
    // writing 100*100/Σ raw here would be a factor-100 error.
    //
    // Expected values from an INDEPENDENT Python run of the reference formulas:
    //   raw {0.5, 0.25} (Σ=0.75): MitStyle 99.25558312655087 ; ParkinsonNormalized 133.33333333333334
    //   raw {1.0}       (Σ=1.0) : MitStyle 99.00990099009901 ; ParkinsonNormalized 100.0
    // -------------------------------------------------------------------------------------------

    [Test]
    public void CfdSpecificity_MitStyle_TwoHits_MatchesReference()
    {
        Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(new[] { 0.5, 0.25 }),
            Is.EqualTo(99.25558312655087).Within(Tol));
    }

    [Test]
    public void CfdSpecificity_ParkinsonNormalized_TwoHits_MatchesReference()
    {
        Assert.That(
            CrisprDesigner.CalculateCfdSpecificityScore(new[] { 0.5, 0.25 }, CfdAggregation.ParkinsonNormalized),
            Is.EqualTo(133.33333333333334).Within(Tol));
    }

    [Test]
    public void CfdSpecificity_SingleUnitHit_BothAggregationsMatchReference()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(new[] { 1.0 }),
                Is.EqualTo(99.00990099009901).Within(Tol));
            Assert.That(
                CrisprDesigner.CalculateCfdSpecificityScore(new[] { 1.0 }, CfdAggregation.ParkinsonNormalized),
                Is.EqualTo(100.0).Within(Tol));
        });
    }

    /// <summary>MitStyle is defined at Σ = 0 and yields the maximum, 100 (no off-targets).</summary>
    [Test]
    public void CfdSpecificity_MitStyle_EmptyAndZeroSum_Returns100()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(Array.Empty<double>()),
                Is.EqualTo(100.0).Within(1e-12));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(new[] { 0.0, 0.0 }),
                Is.EqualTo(100.0).Within(1e-12));
        });
    }

    /// <summary>
    /// ParkinsonNormalized divides by Σ, so Σ = 0 is undefined (the reference would raise
    /// ZeroDivisionError). It must throw rather than return a misleading 100 or infinity.
    /// </summary>
    [Test]
    public void CfdSpecificity_ParkinsonNormalized_ZeroSum_ThrowsArgumentOutOfRange()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.CalculateCfdSpecificityScore(
                Array.Empty<double>(), CfdAggregation.ParkinsonNormalized));
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.CalculateCfdSpecificityScore(
                new[] { 0.0 }, CfdAggregation.ParkinsonNormalized));
        });
    }

    [Test]
    public void CfdSpecificity_NullScores_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScore((IEnumerable<double>)null!));
    }

    // -------------------------------------------------------------------------------------------
    // Genome overload: each hit's CFD needs that hit's OWN PAM, read on the protospacer strand.
    // -------------------------------------------------------------------------------------------

    /// <summary>
    /// The off-target lies on the REVERSE strand and its PAM must be read on the protospacer strand.
    /// Genome = revComp("TACGCATAAAGATGAGACGC" + "AAG"), scanned with SpCas9-NAG, so the one hit is
    /// the reverse-strand protospacer "TACGCATAAAGATGAGACGC" with PAM "AAG".
    /// Per-hit CFD = mismatch rG:dA at position 1 (1.0) x PAM score for "AG" (0.25925925899999996)
    /// = 0.25925925899999996 → MitStyle aggregate 100/(100+0.2592592590)*100 = 99.74141115651946.
    /// If the PAM were reported reverse-complemented ("CTT" → key "TT", score 0.0) the per-hit CFD
    /// would collapse to 0 and the aggregate would be exactly 100, so this value pins the convention.
    /// </summary>
    [Test]
    public void CfdSpecificity_Genome_ReverseStrandHit_UsesProtospacerStrandPam()
    {
        var genome = new DnaSequence("CTTGCGTCTCATCTTTATGCGTA");

        var offTargets = CrisprDesigner
            .FindOffTargets(Guide, genome, 3, CrisprSystemType.SpCas9_NAG).ToList();
        double mitStyle = CrisprDesigner.CalculateCfdSpecificityScore(
            Guide, genome, 3, CrisprSystemType.SpCas9_NAG);
        double parkinson = CrisprDesigner.CalculateCfdSpecificityScore(
            Guide, genome, 3, CrisprSystemType.SpCas9_NAG, CfdAggregation.ParkinsonNormalized);

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].IsForwardStrand, Is.False);
            Assert.That(offTargets[0].Sequence, Is.EqualTo("TACGCATAAAGATGAGACGC"));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 0 }));

            Assert.That(mitStyle, Is.EqualTo(99.74141115651946).Within(Tol));
            Assert.That(mitStyle, Is.Not.EqualTo(100.0));
            Assert.That(parkinson, Is.EqualTo(385.71428610000004).Within(1e-6));
        });
    }

    /// <summary>No off-targets → the MitStyle aggregate is the maximum, 100.</summary>
    [Test]
    public void CfdSpecificity_Genome_NoOffTargets_Returns100()
    {
        var genome = new DnaSequence(Guide + "AGG"); // only the exact on-target
        Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(Guide, genome, 3),
            Is.EqualTo(100.0).Within(1e-12));
    }

    /// <summary>
    /// F10 (regression): the genome overload drives the PAM-aware scan directly instead of
    /// FindOffTargets, and used to skip FindOffTargets' own maxMismatches guard entirely. A negative
    /// cap then silently returned a PERFECT 100 (no hit can satisfy mismatches &lt;= -1) and a cap
    /// above 5 scanned unbounded. Both must now throw, exactly as FindOffTargets does.
    /// </summary>
    [TestCase(-1)]
    [TestCase(6)]
    [TestCase(99)]
    public void CfdSpecificity_Genome_MaxMismatchesOutOfRange_ThrowsArgumentOutOfRange(int maxMismatches)
    {
        var genome = new DnaSequence("CTTGCGTCTCATCTTTATGCGTA");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.CalculateCfdSpecificityScore(
                Guide, genome, maxMismatches, CrisprSystemType.SpCas9_NAG));
            // The MIT overload delegates to FindOffTargets and already enforced the same bound.
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.CalculateMitSpecificityScore(
                Guide, genome, maxMismatches, CrisprSystemType.SpCas9_NAG));
        });
    }

    /// <summary>
    /// CFD is defined only for a 20-nt protospacer, so SaCas9 (21 nt) is refused rather than fudged —
    /// the reference emits the CFD sentinel (cfdScore = -1) for SaCas9 hits too, and provides no
    /// length fallback for CFD (unlike calcHitScore's last-20 fudge for the MIT score).
    /// </summary>
    [Test]
    public void CfdSpecificity_Genome_SaCas9_Refused()
    {
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTCAAGAAT");

        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdSpecificityScore(
            "ACGTACGTACGTACGTACGTA", genome, 3, CrisprSystemType.SaCas9));
    }

    [Test]
    public void CfdSpecificity_Genome_InvalidArguments_Throw()
    {
        var genome = new DnaSequence("CTTGCGTCTCATCTTTATGCGTA");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScore(null!, genome));
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScore("", genome));
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScore(Guide, null!));
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScore("ACGT", genome));
        });
    }
}
