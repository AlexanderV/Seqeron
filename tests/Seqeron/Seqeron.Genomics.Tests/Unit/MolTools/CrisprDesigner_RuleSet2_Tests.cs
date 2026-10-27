using System.Globalization;
using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the Doench et al. 2016 "Rule Set 2" / Azimuth on-target efficacy score
/// (<see cref="CrisprDesigner.CalculateOnTargetRuleSet2(string)"/> and its gene-context overload).
///
/// Rule Set 2 is a trained scikit-learn GradientBoostingRegressor, not a published formula. The C#
/// implementation reads a sklearn-free reconstruction of Microsoft Research's Azimuth pickles
/// (V3_model_nopos / V3_model_full) and reproduces the exact featurization, including the Biopython
/// Tm_NN melting-temperature features and the CPython-2.7 dict column ordering used at training time.
///
/// ORACLE. Upstream's own regression fixture (azimuth/tests/1000guides.csv, columns `truth nopos` /
/// `truth pos`, embedded as TestData/Azimuth/*_oracle.csv column `upstream`) is the independent,
/// third-party oracle: the C# score must reproduce it on ALL 947 guides of BOTH models, to 1e-8
/// (upstream prints 9 decimals; the measured maximum is 5e-10).
///
/// The same CSVs carry `ref_score`, the prediction of the project's Python reference
/// (scripts/azimuth/extract_azimuth_model.py, which loads the shipped pickles without scikit-learn).
/// It agrees with `upstream` on all 947 rows, and the extractor refuses to write these files otherwise.
///
/// HISTORY (audit round 4, H1). An earlier revision of this suite called the fixture "INTERNALLY
/// INCONSISTENT with the shipped pickles for ~38% of rows" and asserted against `upstream` only on the
/// `agrees == 1` subset (585 nopos / 637 full). That was wrong: the 362 (nopos) and 310 (full)
/// excluded rows were exactly the rows where the NGGX one-hot was written in the wrong alphabet order
/// (A,T,C,G instead of upstream's alphabetical A,C,G,T). The subset filter hid the defect.
/// </summary>
[TestFixture]
public class CrisprDesigner_RuleSet2_Tests
{
    private const int ExpectedRows = 947;

    // C# reproduces the reference exactly except for last-ULP differences in Math.Log within Tm; the
    // stored ref_score and upstream truth are rounded to 9 decimals. 1e-8 covers both with margin.
    private const double RefTol = 1e-8;
    private const double UpstreamTol = 1e-8;

    private sealed record NoPos(string Guide, double Ref, double Upstream);
    private sealed record Full(string Guide, int AaCut, double PctPep, double Ref, double Upstream);

    private static IEnumerable<string> ReadResourceLines(string fileName)
    {
        var name = $"Seqeron.Genomics.Tests.TestData.Azimuth.{fileName}";
        var asm = Assembly.GetExecutingAssembly();
        using var stream = asm.GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        string? line;
        bool header = true;
        while ((line = reader.ReadLine()) is not null)
        {
            if (header) { header = false; continue; }
            if (line.Length > 0) yield return line;
        }
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static List<NoPos> LoadNoPos() => ReadResourceLines("nopos_oracle.csv")
        .Select(l => l.Split(','))
        .Select(c => new NoPos(c[0], D(c[1]), D(c[2])))
        .ToList();

    private static List<Full> LoadFull() => ReadResourceLines("full_oracle.csv")
        .Select(l => l.Split(','))
        .Select(c => new Full(c[0], int.Parse(c[1], CultureInfo.InvariantCulture), D(c[2]), D(c[3]), D(c[4])))
        .ToList();

    // -------------------------------------------------------------------------------------------
    // Primary oracle: C# reproduces the verified Python reference for every guide.
    // -------------------------------------------------------------------------------------------

    [Test]
    public void RuleSet2_NoPos_ReproducesVerifiedReference_AllGuides()
    {
        var rows = LoadNoPos();
        rows.Should().HaveCount(ExpectedRows);

        double maxErr = 0;
        string worst = "";
        foreach (var r in rows)
        {
            double got = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide);
            double err = Math.Abs(got - r.Ref);
            if (err > maxErr) { maxErr = err; worst = r.Guide; }
        }
        maxErr.Should().BeLessThan(RefTol, $"worst guide was {worst}");
    }

    [Test]
    public void RuleSet2_Full_ReproducesVerifiedReference_AllGuides()
    {
        var rows = LoadFull();
        rows.Should().HaveCount(ExpectedRows);

        double maxErr = 0;
        string worst = "";
        foreach (var r in rows)
        {
            double got = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, r.AaCut, r.PctPep);
            double err = Math.Abs(got - r.Ref);
            if (err > maxErr) { maxErr = err; worst = r.Guide; }
        }
        maxErr.Should().BeLessThan(RefTol, $"worst guide was {worst}");
    }

    // -------------------------------------------------------------------------------------------
    // Independent oracle: upstream's own fixture (azimuth/tests/1000guides.csv), every row, both models.
    // -------------------------------------------------------------------------------------------

    [Test]
    public void RuleSet2_NoPos_MatchesUpstreamFixture_AllGuides()
    {
        var rows = LoadNoPos();
        rows.Should().HaveCount(ExpectedRows);

        foreach (var r in rows)
        {
            double got = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide);
            got.Should().BeApproximately(r.Upstream, UpstreamTol, $"guide {r.Guide}");
        }
    }

    [Test]
    public void RuleSet2_Full_MatchesUpstreamFixture_AllGuides()
    {
        var rows = LoadFull();
        rows.Should().HaveCount(ExpectedRows);

        foreach (var r in rows)
        {
            double got = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, r.AaCut, r.PctPep);
            got.Should().BeApproximately(r.Upstream, UpstreamTol, $"guide {r.Guide}");
        }
    }

    /// <summary>
    /// The NGGX block (PAM N at offset 24 x the +1 base at offset 27) is the one block in alphabetical
    /// A,C,G,T order; the others keep A,T,C,G. The all-guides tests above only constrain that ordering
    /// for the cells the fixture actually hits, so pin that it hits all sixteen, in both models.
    /// </summary>
    [Test]
    public void RuleSet2_Fixture_ExercisesAllSixteenNggxCells_InBothModels()
    {
        LoadNoPos().Select(r => $"{r.Guide[24]}{r.Guide[27]}").Distinct().Should().HaveCount(16);
        LoadFull().Select(r => $"{r.Guide[24]}{r.Guide[27]}").Distinct().Should().HaveCount(16);
    }

    // -------------------------------------------------------------------------------------------
    // Behavioural properties.
    // -------------------------------------------------------------------------------------------

    [Test]
    public void RuleSet2_Score_IsInUnitInterval()
    {
        foreach (var r in LoadNoPos())
            CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide).Should().BeInRange(0.0, 1.0);
    }

    [Test]
    public void RuleSet2_IsCaseInsensitive()
    {
        const string g = "AAAAAAAAAAAAAAAAAAATGCAGCGGGAG";
        CrisprDesigner.CalculateOnTargetRuleSet2(g.ToLowerInvariant())
            .Should().Be(CrisprDesigner.CalculateOnTargetRuleSet2(g));
    }

    [Test]
    public void RuleSet2_IsDeterministic()
    {
        const string g = "AAAAAAAAAAAAAAAAAAATGCAGCGGGAG";
        CrisprDesigner.CalculateOnTargetRuleSet2(g)
            .Should().Be(CrisprDesigner.CalculateOnTargetRuleSet2(g));
    }

    [Test]
    public void RuleSet2_FullModel_AddsGeneContextSignal_OverNoPos()
    {
        // The full model incorporates gene position, so it generally differs from the sequence-only score.
        var r = LoadFull().First();
        double nopos = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide);
        double full = CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, r.AaCut, r.PctPep);
        full.Should().NotBe(nopos);
    }

    // -------------------------------------------------------------------------------------------
    // Input validation.
    // -------------------------------------------------------------------------------------------

    /// <summary>
    /// A null reference is a missing argument (<see cref="ArgumentNullException"/>); an empty string
    /// is a wrong VALUE and is reported like every other malformed context, i.e. as a plain
    /// <see cref="ArgumentException"/> carrying the length message. The reference separates the two
    /// the same way (<c>model_comparison.py</c> L523 <c>assert len(seq[0]) > 0</c> vs L524
    /// <c>assert isinstance(seq[0], basestring)</c>). This test previously locked
    /// <see cref="ArgumentNullException"/> for <c>""</c>, which was the wart, not the contract.
    /// </summary>
    [Test]
    public void RuleSet2_NullOrEmpty_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(null!));
        Assert.Throws<ArgumentNullException>(
            () => CrisprDesigner.CalculateOnTargetRuleSet2(null!, 10, 25.0));

        var empty = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(""));
        empty.Should().NotBeOfType<ArgumentNullException>();
        empty!.Message.Should().Contain("got 0");
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetRuleSet2("", 10, 25.0))
            .Should().NotBeOfType<ArgumentNullException>();
    }

    [TestCase("ACGTACGT")]                          // too short
    [TestCase("AAAAAAAAAAAAAAAAAAATGCAGCGGGAGTT")]   // too long (31)
    public void RuleSet2_WrongLength_Throws(string seq)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(seq));
    }

    [Test]
    public void RuleSet2_NonAcgt_Throws()
    {
        // 30-mer with an 'N'.
        Assert.Throws<ArgumentException>(
            () => CrisprDesigner.CalculateOnTargetRuleSet2("AAAAAAAAAAAAAAAAAANTGCAGCGGGAG"));
    }

    /// <summary>G11: the case fold is ASCII-only; 'ſ' (U+017F) is reported as itself, not as the 'S' it used to fold to.</summary>
    [Test]
    public void RuleSet2_LongS_IsReportedAsItself_AndLowerCaseIsStillFolded()
    {
        const string valid = "ACGTACGTACGTACGTACGTACGTAGGTAC";
        var e = Assert.Throws<ArgumentException>(
            () => CrisprDesigner.CalculateOnTargetRuleSet2(valid.Remove(3, 1).Insert(3, "ſ")));
        e!.Message.Should().Contain("'ſ'");
        CrisprDesigner.CalculateOnTargetRuleSet2(valid.ToLowerInvariant())
            .Should().Be(CrisprDesigner.CalculateOnTargetRuleSet2(valid));
    }

    // Round 4 H22. CRISPOR scores every SpCas9-class PAM with Azimuth's `pam_audit=False` ("do not check for NGG PAM",
    // crisporEffScores.py L1174-1178), so the public scorer no longer rejects a non-GG 30-mer. Oracle: the project's
    // tree walker over the ORIGINAL V3_model_*.pickle with the sorted-NGGX featurizer (the same one that reproduces
    // Azimuth's 947-row fixture to 5e-10), executed on each context; differential on 300 random 30-mers (281 with a
    // non-GG PAM) for both models: max |C# - walker| = 0.0.
    [TestCase("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 181, 50.0, 0.31571342050275053, 0.34306783305466304)]
    [TestCase("ACGTGACGCATAAAGATGAGACGCTAGTGA", 181, 50.0, 0.6717138839537462, 0.6591691589095564)]   // NAG
    [TestCase("ACGTGACGCATAAAGATGAGACGCCGATGA", 181, 50.0, 0.7128489249311819, 0.6493995783608925)]   // NGA
    [TestCase("TTAGGCTAGCTAGGCTAGCATCGAAAGCAT", 12, 75.0, 0.5356722721642655, 0.4858958755113498)]    // NAG
    [TestCase("CCGTTTGCAAGCTAGCTAGGCTACGGAATA", 0, 10.0, 0.6009388007571157, 0.613497372682248)]      // NGA
    public void RuleSet2_NonNggPam_IsScored_AsCrisporDoesWithoutPamAudit(
        string context, int aaCut, double percentPeptide, double expectedNoPos, double expectedFull)
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateOnTargetRuleSet2(context), Is.EqualTo(expectedNoPos).Within(1e-9));
            Assert.That(CrisprDesigner.CalculateOnTargetRuleSet2(context, aaCut, percentPeptide),
                Is.EqualTo(expectedFull).Within(1e-9));
        });
    }

    /// <summary>
    /// H22 end to end: SpCas9-NAG / SpCas9-NGA guides now carry a Context30Mer and an OnTargetScore, so
    /// <c>GuideRnaRanking.OnTargetRuleSet2</c> orders them instead of degenerating to position order.
    /// </summary>
    [Test]
    public void DesignedGuides_OfTheNagAndNgaSystems_CarryARuleSet2Score()
    {
        const string guide = "GACGCATAAAGATGAGACGC";
        foreach (var (system, pam) in new[] { (CrisprSystemType.SpCas9_NAG, "TAG"), (CrisprSystemType.SpCas9_NGA, "CGA") })
        {
            var text = "ACGT" + guide + pam + "TGA" + "TTCATCAGTA";
            var designed = CrisprDesigner.DesignGuideRnas(new DnaSequence(text), 0, text.Length - 1, system,
                GuideRnaParameters.Default with { MinScore = 0 }).Where(c => c.Sequence == guide).ToList();

            Assert.That(designed, Has.Count.EqualTo(1), system.ToString());
            Assert.That(designed[0].Context30Mer, Is.EqualTo("ACGT" + guide + pam + "TGA"));
            Assert.That(designed[0].OnTargetScore,
                Is.EqualTo(CrisprDesigner.CalculateOnTargetRuleSet2("ACGT" + guide + pam + "TGA")).Within(1e-12));
        }
    }

    // Round 4 H21: the amino-acid cut position is a residue count; a negative one is refused (documented), and the
    // -1 percent-peptide sentinel makes the argument irrelevant, so a negative cut position is accepted there.
    [Test]
    public void RuleSet2_AminoAcidCutPosition_NegativeIsRefused_UnlessThePercentPeptideSentinelIgnoresIt()
    {
        const string context = "ACGTGACGCATAAAGATGAGACGCAGGTGA";
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(context, -1, 50.0));
            Assert.That(CrisprDesigner.CalculateOnTargetRuleSet2(context, -1, -1.0),
                Is.EqualTo(CrisprDesigner.CalculateOnTargetRuleSet2(context)));
            Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(null!));
            Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetRuleSet2(""));
        });
    }

    // Round 4 H23: a malformed model blob fails as InvalidOperationException - including a header whose
    // nodeCount * 24 overflows a 32-bit int (89 478 486 * 24 = 2 147 483 664 wrapped negative, the old truncation
    // guard passed and the node slice threw ArgumentOutOfRangeException) and a node or tree start pointing outside.
    private static byte[] Blob(int treeCount, int nodeCount, int featureCount, bool genePosition,
        IReadOnlyList<int> treeStarts, IReadOnlyList<(double thr, int left, int right, int feature)> nodes)
    {
        using var ms = new MemoryStream();
        using var w = new BinaryWriter(ms);
        w.Write(0x32535241u); w.Write((ushort)1); w.Write((ushort)(genePosition ? 1 : 0));
        w.Write(treeCount); w.Write(nodeCount); w.Write(featureCount); w.Write(0.5); w.Write(0.1);
        foreach (int t in treeStarts) w.Write(t);
        foreach (var n in nodes) { w.Write(n.thr); w.Write(n.left); w.Write(n.right); w.Write(n.feature); w.Write(0); }
        return ms.ToArray();
    }

    [Test]
    public void ModelBlob_Malformed_FailsAsInvalidOperationException()
    {
        var leaf = (0.25, -1, -1, -1);
        var ok = Blob(1, 1, 627, false, new[] { 0 }, new[] { leaf });
        Assert.DoesNotThrow(() => AzimuthRuleSet2.ValidateModelBlob(ok, "ok", expectGenePosition: false));

        Assert.Multiple(() =>
        {
            // nodeCount * 24 overflows int: header alone, 36 bytes
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(
                Blob(1, 89_478_486, 627, false, new[] { 0 }, Array.Empty<(double, int, int, int)>()), "overflow", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(
                Blob(1, 1, 627, false, new[] { 5 }, new[] { leaf }), "treeStart", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(
                Blob(1, 1, 627, false, new[] { -1 }, new[] { leaf }), "treeStartNegative", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(
                Blob(1, 3, 627, false, new[] { 0 }, new[] { (0.5, 1, 2, 9999), leaf, leaf }), "feature", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(
                Blob(1, 3, 627, false, new[] { 0 }, new[] { (0.5, 1, 7, 3), leaf, leaf }), "child", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(ok.AsSpan(0, 20), "truncatedHeader", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(ok.AsSpan(0, ok.Length - 1), "truncatedBody", false));
            Assert.Throws<InvalidOperationException>(() => AzimuthRuleSet2.ValidateModelBlob(ok, "geneFlag", expectGenePosition: true));
        });
    }

    // -------------------------------------------------------------------------------------------
    // Tm featurization: AzimuthRuleSet2 delegates its four melting-temperature features to the
    // canonical ThermoConstants.CalculateNearestNeighborTm instead of carrying its own DNA_NN3 table.
    // Azimuth's featurization.py calls the obsolete Biopython MeltingTemp.Tm_staluc (Tm_feature,
    // L381-L384), which Biopython defines as Tm_NN(s, dnac1=dnac/2.0, dnac2=dnac/2.0, Na=saltc) with
    // dnac=50 / saltc=50 (Bio/SeqUtils/MeltingTemp.py L1087/L1127 at tag biopython-176; the function
    // is gone in 1.88) -- i.e. Tm_NN's own defaults nn_table=DNA_NN3, dnac1 = dnac2 = 25 nM,
    // Na = 50 mM, saltcorr = 5, which are exactly the defaults of the canonical method. 1.76's
    // Tm_staluc and 1.88's Tm_NN were executed side by side on every segment below and agree
    // bitwise. The values below were produced by REAL Biopython 1.88
    // (python3 -c "from Bio.SeqUtils import MeltingTemp as mt; mt.Tm_NN(s, nn_table=mt.DNA_NN3,
    // dnac1=25, dnac2=25, Na=50, saltcorr=5)") and are asserted BIT-EXACTLY, so a drift in either
    // implementation breaks this test and would change the Rule Set 2 score.
    // -------------------------------------------------------------------------------------------

    [TestCase("ACGTACGTACGTACGTACGTACGTAGGACG", 62.55983015288416)]  // the whole 30-mer feature
    [TestCase("TAGG", -67.90079099878676)]                           // seq[19:24)-style short window
    [TestCase("ACGTA", -34.426168537504026)]
    [TestCase("AAAAA", -50.31930198408023)]
    [TestCase("GCGCGCGC", 36.45901389113385)]
    public void RuleSet2_TmFeature_IsBiopythonTmNnBitExact(string segment, double biopythonTm)
    {
        ThermoConstants.CalculateNearestNeighborTm(segment).Should().Be(biopythonTm);
    }

    /// <summary>
    /// The Rule Set 2 score depends on the four Tm features, so this guards the delegation end to
    /// end: perturbing the canonical Tm would move these scores away from the oracle.
    /// </summary>
    [Test]
    public void RuleSet2_TmFeature_FeedsTheScore_OracleStillReproduced()
    {
        var rows = LoadNoPos().Take(25).ToList();
        rows.Should().HaveCount(25);
        foreach (var r in rows)
            CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide).Should().BeApproximately(r.Ref, RefTol);
    }

    // -------------------------------------------------------------------------------------------
    // Azimuth's "no gene position" sentinel (percent_peptide == -1).
    //
    // The reference does NOT feed -1 into the gene-context model. azimuth/model_comparison.py:
    //   L539  if np.any(percent_peptide == -1) or (percent_peptide is None and aa_cut is None):
    //             model_name = 'V3_model_nopos.pickle'
    //   L563  if np.all(percent_peptide != -1) and (percent_peptide is not None and aa_cut is not None):
    //             gene_position = DataFrame(... percent_peptide, aa_cut ...)
    //         else: gene_position = DataFrame(... -1, -1 ...)
    // i.e. the sentinel selects V3_model_nopos AND discards the caller's gene position, so the
    // answer is exactly the sequence-only score. The library used to answer with V3_model_full and
    // percent_peptide = -1 / "Percent Peptide <50%" = 1, which is a score the reference never
    // returns for any input (e.g. 0.608099450578 instead of 0.500929648 for the guide below).
    // -------------------------------------------------------------------------------------------

    [Test]
    public void RuleSet2_PercentPeptideSentinel_SelectsTheNoPosModel_ForEveryOracleGuide()
    {
        var rows = LoadNoPos();
        rows.Should().HaveCount(ExpectedRows);
        foreach (var r in rows)
        {
            // Any amino-acid cut position must be irrelevant once the sentinel is passed, because
            // the nopos model has no gene-position columns at all.
            CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, 181, -1.0)
                .Should().Be(CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide));
            CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, 0, -1.0)
                .Should().Be(CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide));
            CrisprDesigner.CalculateOnTargetRuleSet2(r.Guide, 181, -1.0)
                .Should().BeApproximately(r.Ref, RefTol);
        }
    }

    /// <summary>
    /// The one guide where the two models are furthest apart in the oracle: the sentinel must give
    /// the verified nopos reference score (nopos_oracle.csv row 1, ref_score 0.500929648 = upstream's `truth nopos`), not the
    /// full-model value the implementation used to return for percent_peptide = -1.
    /// </summary>
    [Test]
    public void RuleSet2_PercentPeptideSentinel_ReturnsTheVerifiedNoPosScore_NotTheFullModelOne()
    {
        const string guide = "AAAAAAAAAAAAAAAAAAATGCAGCGGGAG"; // full_oracle/nopos_oracle row 1
        double sentinel = CrisprDesigner.CalculateOnTargetRuleSet2(guide, 181, -1.0);

        sentinel.Should().BeApproximately(0.500929648, RefTol);       // verified nopos reference
        sentinel.Should().NotBeApproximately(0.608099, 1e-4);          // pre-fix full-model answer
        CrisprDesigner.CalculateOnTargetRuleSet2(guide, 181, 50.0)
            .Should().BeApproximately(0.570371166, RefTol);               // full model (upstream `truth pos`)
    }

    /// <summary>
    /// percent_peptide outside the documented [0, 100] percentage domain is a caller error. The
    /// reference only asserts it is a real number (model_comparison.py L534) and has no range check,
    /// so this is a Seqeron-side contract (the XML doc already promised [0, 100]); -1 stays accepted
    /// because the reference uses it as its own sentinel (L539/L563).
    /// </summary>
    [TestCase(-2.0)]
    [TestCase(-1.0001)]
    [TestCase(-0.5)]
    [TestCase(100.0001)]
    [TestCase(101.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void RuleSet2_PercentPeptideOutOfRange_Throws(double percentPeptide)
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CrisprDesigner.CalculateOnTargetRuleSet2(
                "AAAAAAAAAAAAAAAAAAATGCAGCGGGAG", 181, percentPeptide));
    }

    /// <summary>
    /// "Amino Acid Cut position" is a residue count. The reference asserts only that it is real
    /// (model_comparison.py, <c>np.isreal(aa_cut)</c>); azimuth/tests/1000guides.csv (column "AA cut", 947 guides)
    /// spans the integers 0..499 (0 occurs twice, so the bound is &gt;= 0, not &gt;= 1) and its percent
    /// peptide spans 0..99. A negative count with a real percent_peptide is a caller error; with the
    /// percent_peptide sentinel -1 the argument is ignored by the reference, so it stays unchecked there.
    /// </summary>
    [TestCase(-1)]
    [TestCase(-2)]
    [TestCase(int.MinValue)]
    public void RuleSet2_AminoAcidCutPositionNegative_Throws(int aminoAcidCutPosition)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => CrisprDesigner.CalculateOnTargetRuleSet2(
                "AAAAAAAAAAAAAAAAAAATGCAGCGGGAG", aminoAcidCutPosition, 50.0));
        Assert.That(ex!.ParamName, Is.EqualTo("aminoAcidCutPosition"));
        // sentinel percent_peptide: the reference replaces both columns by -1, so aa_cut is not looked at
        Assert.That(
            CrisprDesigner.CalculateOnTargetRuleSet2("AAAAAAAAAAAAAAAAAAATGCAGCGGGAG", aminoAcidCutPosition, -1.0),
            Is.EqualTo(CrisprDesigner.CalculateOnTargetRuleSet2("AAAAAAAAAAAAAAAAAAATGCAGCGGGAG")));
    }

    [TestCase(0)]      // occurs in the reference fixture
    [TestCase(181)]
    [TestCase(499)]    // the fixture's maximum
    [TestCase(100000)] // no upper bound in the reference
    public void RuleSet2_AminoAcidCutPositionInRange_IsAccepted(int aminoAcidCutPosition)
    {
        CrisprDesigner.CalculateOnTargetRuleSet2("AAAAAAAAAAAAAAAAAAATGCAGCGGGAG", aminoAcidCutPosition, 50.0)
            .Should().BeInRange(0.0, 1.0);
    }

    [TestCase(0.0)]
    [TestCase(0.18)]   // the fractional value used in Azimuth's own README example
    [TestCase(50.0)]
    [TestCase(100.0)]
    [TestCase(-1.0)]   // the sentinel
    public void RuleSet2_PercentPeptideInRange_IsAccepted(double percentPeptide)
    {
        CrisprDesigner.CalculateOnTargetRuleSet2("AAAAAAAAAAAAAAAAAAATGCAGCGGGAG", 181, percentPeptide)
            .Should().BeInRange(0.0, 1.0);
    }

    // -------------------------------------------------------------------------------------------
    // Model-blob shape. The featurizer writes fixed column offsets (gene position at 1 / 607 / 629
    // of 630 for the full model; 627 columns for nopos), so a blob whose header declares a different
    // shape must be rejected by the LOADER as an InvalidOperationException -- the same failure mode
    // as the magic-number check -- instead of reaching the featurizer and throwing
    // IndexOutOfRangeException from FillGenePosition.
    // -------------------------------------------------------------------------------------------

    private static (Type Type, MethodInfo Load) AzimuthInternals()
    {
        var type = typeof(CrisprDesigner).Assembly
            .GetType("Seqeron.Genomics.MolTools.AzimuthRuleSet2", throwOnError: true)!;
        var load = type.GetMethod("Load", BindingFlags.NonPublic | BindingFlags.Static)!;
        return (type, load);
    }

    [TestCase("azimuth_rs2_nopos.bin", false, 627, false)]
    [TestCase("azimuth_rs2_full.bin", true, 630, true)]
    public void RuleSet2_ModelBlob_DeclaresTheShapeTheFeaturizerAssumes(
        string fileName, bool expectGenePosition, int expectedFeatureCount, bool expectedFlag)
    {
        var (type, load) = AzimuthInternals();
        object model = load.Invoke(null, [fileName, expectGenePosition])!;
        var modelType = model.GetType();
        modelType.GetField("FeatureCount")!.GetValue(model).Should().Be(expectedFeatureCount);
        modelType.GetField("HasGenePosition")!.GetValue(model).Should().Be(expectedFlag);
        type.Should().NotBeNull();
    }

    [TestCase("azimuth_rs2_nopos.bin", true)]   // nopos blob offered as the gene-context model
    [TestCase("azimuth_rs2_full.bin", false)]   // full blob offered as the sequence-only model
    public void RuleSet2_ModelBlob_WrongDeclaredShape_ThrowsInvalidOperation(
        string fileName, bool expectGenePosition)
    {
        var (_, load) = AzimuthInternals();
        var thrown = Assert.Throws<TargetInvocationException>(
            () => load.Invoke(null, [fileName, expectGenePosition]));
        thrown!.InnerException.Should().BeOfType<InvalidOperationException>();
    }
}
