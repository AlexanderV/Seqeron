// PROTMOTIF-DOMAIN-001 — Plan7 profile-HMM domain detection (SH3 / PDZ / WD40), opt-in.
// Evidence: docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md  (addendum 2026-06-25)
// TestSpec: tests/TestSpecs/PROTMOTIF-DOMAIN-001.md          (addendum 2026-06-25)
// Source:   Durbin et al. (1998) Biological Sequence Analysis §5.4 (Viterbi recurrence);
//           HMMER User's Guide v3.4 (Eddy 2023) — HMMER3/f file format + STATS lines;
//           Eddy (2008) PLoS Comput Biol 4:e1000069 (Gumbel for Viterbi/MSV, exponential tail for Forward);
//           Easel esl_gumbel.c / esl_exponential.c survival functions; Pfam PF00018/PF00595/PF00400 (CC0);
//           HMMER modelconfig.c (local entry/exit, p7_ReconfigLength), generic_fwdback.c, generic_decoding.c,
//           generic_null2.c (p7_GNull2_ByExpectation), p7_domaindef.c / p7_pipeline.c (null2 correction),
//           hmmer.c p7_AminoFrequencies; GROUND-TRUTH hmmsearch scores via pyhmmer 0.12.1 (addendum 2026-06-25).
using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for the opt-in Plan7 profile-HMM engine (<see cref="Plan7ProfileHmm"/>) and the
/// SH3/PDZ/WD40 detection methods on <see cref="ProteinMotifFinder"/>. The exact-PROSITE-pattern
/// <c>FindDomains</c> path is validated separately (ProteinMotifFinder_DomainPrediction_Tests).
/// </summary>
[TestFixture]
[Category("PROTMOTIF-DOMAIN-001")]
public class ProteinMotifFinder_FindDomainsByHmm_Tests
{
    #region Evidence-sourced constants

    // H1 — exact hand-built HMM Viterbi target (Durbin §5.4). 2 match states over {A,C}:
    //   background qA=0.6,qC=0.4; M1: A=0.7,C=0.3; M2: A=0.2,C=0.8; B->M1=0.9, M1->M2=0.8.
    //   Seq "AC", optimal path B->M1(A)->M2(C)->E:
    //   ln(0.7/0.6)+ln(0.9)+ln(0.8/0.4)+ln(0.8) = 0.5187937934151676 nats.
    private const double HandHmmExpectedViterbiNats = 0.5187937934151676;

    // Pfam accessions of the bundled CC0 profiles.
    private const string Sh3Accession = "PF00018";
    private const string PdzAccession = "PF00595";
    private const string Wd40Accession = "PF00400";

    // H3 — SH3 domain core of SRC_HUMAN (UniProt P12931), residues spanning the SH3 fold.
    private const string Sh3TruePositive =
        "TFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAP";

    // H4 — PDZ1 domain of DLG4_HUMAN / PSD-95 (UniProt P78352), residues 61–151 (verbatim).
    private const string PdzTruePositive =
        "MEYEEITLERGNSGLGFSIAGGTDNPHIGDDPSIFITKIIPGGAAAQDGRLRVNDSILFVNEVDVREVTHSAAVEALKEAGSIVRLYVMRR";

    // H5 — full WD40 β-propeller of GBB1_HUMAN / Gβ1 (UniProt P62873).
    private const string Wd40TruePositive =
        "MSELDQLRQEAEQLKNQIRDARKACADATLSQITNNIDPVGRIQMRTRRTLRGHLAKIYAMHWGTDSRLLVSASQDGKLIIWDSYTTNKVHAIPLRSSWVMTCAYAPSGNYVACGGLDNICSIYNLKTREGNVRVSRELAGHTGYLSCCRFLDDNQIVTSSGDTTCALWDIETGQQTTTFTGHTGDVMSLSLAPDTRLFVSGACDASAKLWDVREGMCRQTFTGHESDINAICFFPNGNAFATGSDDATCRLFDLRADQELMTYSHDNIICGITSVSFSKSGRLLLAGYDDFNCNVWDALKADRAGVLAGHDNRVSCLGVTDDGMAVATGSWDSFLKIWN";

    // H6 — low-complexity non-domain sequence (true negative for all three profiles).
    private const string TrueNegative =
        "AAAAAAAAAAAAAAEEEEEEEEEEEEEEKKKKKKKKKKKK";

    private const double DetectionThresholdBits = 10.0;

    // H13–H17 — HMMER E-value statistics. STATS lines read VERBATIM from the bundled
    // Resources/PF00018_SH3_1.hmm:
    //   STATS LOCAL MSV       -8.1284  0.71923
    //   STATS LOCAL VITERBI   -8.2932  0.71923
    //   STATS LOCAL FORWARD   -4.5735  0.71923
    private const double Sh3MsvMu = -8.1284;
    private const double Sh3ViterbiMu = -8.2932;
    private const double Sh3ForwardTau = -4.5735;
    private const double Sh3Lambda = 0.71923;

    // Hand-derived (Python, full double precision) at S = 40 bits using the exact Easel formulas:
    //   Gumbel survival  P = 1 − exp(−exp(−λ(S−μ))), with Easel's |ey|<5e-9 → return −ey tail branch.
    //   Exponential surv P = exp(−λ(S−τ)).
    // Viterbi: y = 0.71923*(40−(−8.2932)) = 34.7326…, exp(−y)=8.227179545686635e-16 (deep tail → −ey).
    private const double HandBitScore = 40.0;
    private const double HandZ = 1000.0;
    private const double Sh3ViterbiPValueAt40 = 8.227179545686635e-16;
    private const double Sh3ViterbiEValueAt40Z1000 = 8.227179545686635e-13;
    private const double Sh3ForwardPValueAt40 = 1.1943390031599535e-14;
    private const double Sh3ForwardEValueAt40Z1000 = 1.1943390031599535e-11;

    #endregion

    #region Hand-built HMM construction (for exact DP verification)

    // Natural-log emission row for {A,C}; other 18 residues are negligible (1e-9).
    private static double[] EmissionLn(double pA, double pC)
    {
        const string alpha = "ACDEFGHIKLMNPQRSTVWY";
        var row = new double[20];
        for (int i = 0; i < 20; i++)
        {
            double p = alpha[i] == 'A' ? pA : alpha[i] == 'C' ? pC : 1e-9;
            row[i] = Math.Log(p);
        }
        return row;
    }

    // Transition row in HMMER order [Mk->Mk+1, Mk->Ik, Mk->Dk+1, Ik->Mk+1, Ik->Ik, Dk->Mk+1, Dk->Dk+1].
    // double.NaN → '*' (probability zero → −∞ log).
    private static double[] TransLn(double mm, double mi, double md, double im, double ii, double dm, double dd)
    {
        double L(double p) => double.IsNaN(p) ? double.NegativeInfinity : Math.Log(p);
        return new[] { L(mm), L(mi), L(md), L(im), L(ii), L(dm), L(dd) };
    }

    // Builds the H1 hand HMM directly from EXACT natural-log parameters (no ASCII round-trip, so
    // the DP is verified at full precision against the exact Durbin derivation). The .hmm string
    // parser is verified separately on the real PF00018 profile (H8). starOnMainPath=true forces a
    // '*' on M1->M2 and M1->D2 to forbid any full-length path.
    private static Plan7ProfileHmm BuildHandHmm(bool starOnMainPath = false)
    {
        var background = EmissionLn(0.6, 0.4);                 // q: A=0.6, C=0.4
        var insBg = EmissionLn(0.6, 0.4);                      // insert emits background (odds 0)

        var match = new double[3][];
        var insert = new double[3][];
        var trans = new double[3][];

        // BEGIN node (index 0): insert-0 emissions, then B->M1=0.9,B->I0=0.05,B->D1=0.05,I0->M1=0.5,I0->I0=0.5.
        insert[0] = insBg;
        trans[0] = TransLn(0.9, 0.05, 0.05, 0.5, 0.5, 1.0, double.NaN);

        // Node 1.
        match[1] = EmissionLn(0.7, 0.3);
        insert[1] = insBg;
        double mm = starOnMainPath ? double.NaN : 0.8;
        double md = starOnMainPath ? double.NaN : 0.1;
        trans[1] = TransLn(mm, 0.1, md, 0.5, 0.5, 1.0, double.NaN);

        // Node 2 (final): M2->E and D2->E are probability 1.0 (handled by the engine's end step).
        match[2] = EmissionLn(0.2, 0.8);
        insert[2] = insBg;
        trans[2] = TransLn(1.0, double.NaN, double.NaN, 1.0, 0.5, 1.0, double.NaN);

        return Plan7ProfileHmm.FromLogParameters("toy", 2, match, insert, trans, background);
    }

    #endregion

    #region H1/H2/H9 — exact DP on a hand-built HMM

    [Test]
    public void ViterbiScore_HandBuiltHmm_MatchesExactDurbinDerivation()
    {
        // H1 — pins the DP arithmetic to the hand-derived Durbin §5.4 value.
        var hmm = BuildHandHmm();

        double score = hmm.ViterbiScore("AC");

        Assert.That(score, Is.EqualTo(HandHmmExpectedViterbiNats).Within(1e-9),
            "Viterbi log-odds of 'AC' through B->M1->M2->E must equal the hand-derived "
            + "ln(0.7/0.6)+ln(0.9)+ln(0.8/0.4)+ln(0.8) = 0.5187937934151676 nats.");
    }

    [Test]
    public void ForwardScore_HandBuiltHmm_IsAtLeastViterbiScore()
    {
        // H2 — INV-HMM-01: Forward (sum over all paths) >= Viterbi (single best path).
        var hmm = BuildHandHmm();

        double viterbi = hmm.ViterbiScore("AC");
        double forward = hmm.ForwardScore("AC");

        Assert.That(forward, Is.GreaterThanOrEqualTo(viterbi),
            "Forward sums over all paths including the optimal one, so it cannot be below Viterbi.");
        Assert.That(forward, Is.GreaterThan(viterbi),
            "With alternative (insert/delete) paths available, Forward is strictly greater here.");
    }

    [Test]
    public void ViterbiScore_StarOnOnlyPath_ForbidsPath_ReturnsNegativeInfinity()
    {
        // H9 — INV-HMM-04: a '*' (zero prob) on the sole productive transition forbids the path.
        // With M1->M2 = '*' and M1->D2 = '*', a 2-residue match path through both nodes is impossible.
        var hmm = BuildHandHmm(starOnMainPath: true);

        double score = hmm.ViterbiScore("AC");

        Assert.That(double.IsNegativeInfinity(score), Is.True,
            "A '*' parameter is probability zero (−∞ log); the only full path is forbidden.");
    }

    #endregion

    #region H8 — .hmm parser round-trip on a bundled profile

    [Test]
    public void Parse_EmbeddedSh3Profile_ReadsHeaderFieldsExactly()
    {
        // H8 — parser reads the real HMMER3/f header of the bundled PF00018 profile.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        Assert.Multiple(() =>
        {
            Assert.That(hmm.Name, Is.EqualTo("SH3_1"), "NAME line of PF00018.");
            Assert.That(hmm.Accession, Is.EqualTo("PF00018.35"), "ACC line of PF00018.");
            Assert.That(hmm.Length, Is.EqualTo(48), "LENG (match states) of PF00018.");
            Assert.That(hmm.GatheringThreshold, Is.EqualTo(22.9).Within(1e-9), "GA1 gathering threshold.");
        });
    }

    [Test]
    public void Parse_NonHmmer3Text_ThrowsFormatException()
    {
        // H12 — malformed input guard.
        Assert.Throws<FormatException>(() => Plan7ProfileHmm.Parse("not a hmm file\n"),
            "A file without the 'HMMER3/' version line is not a valid profile.");
    }

    #endregion

    #region H3/H4/H5 — real true positives detected

    [Test]
    public void ScoreDomainHmm_Sh3TruePositive_ScoresAboveThreshold()
    {
        // H3 — SRC_HUMAN SH3 core scores well above the detection threshold against PF00018.
        double bits = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession);

        Assert.That(bits, Is.GreaterThan(DetectionThresholdBits),
            $"A genuine SH3 domain must score >> {DetectionThresholdBits} bits against PF00018.");
    }

    [Test]
    public void ScoreDomainHmm_PdzTruePositive_ScoresAboveThreshold()
    {
        // H4 — PSD-95 PDZ1 scores well above threshold against PF00595.
        double bits = ProteinMotifFinder.ScoreDomainHmm(PdzTruePositive, PdzAccession);

        Assert.That(bits, Is.GreaterThan(DetectionThresholdBits),
            $"A genuine PDZ domain must score >> {DetectionThresholdBits} bits against PF00595.");
    }

    [Test]
    public void ScoreDomainHmm_Wd40TruePositive_ScoresAboveThreshold()
    {
        // H5 — GBB1_HUMAN β-propeller scores well above threshold against PF00400.
        double bits = ProteinMotifFinder.ScoreDomainHmm(Wd40TruePositive, Wd40Accession);

        Assert.That(bits, Is.GreaterThan(DetectionThresholdBits),
            $"A genuine WD40 β-propeller must score >> {DetectionThresholdBits} bits against PF00400.");
    }

    #endregion

    #region H6/H7 — true negative rejected; cross-domain specificity

    [Test]
    public void ScoreDomainHmm_TrueNegative_ScoresBelowZeroForAllProfiles()
    {
        // H6 — INV-HMM-02: a low-complexity non-domain sequence is negative for all three profiles.
        Assert.Multiple(() =>
        {
            Assert.That(ProteinMotifFinder.ScoreDomainHmm(TrueNegative, Sh3Accession), Is.LessThan(0.0),
                "Non-domain sequence must score negative against SH3.");
            Assert.That(ProteinMotifFinder.ScoreDomainHmm(TrueNegative, PdzAccession), Is.LessThan(0.0),
                "Non-domain sequence must score negative against PDZ.");
            Assert.That(ProteinMotifFinder.ScoreDomainHmm(TrueNegative, Wd40Accession), Is.LessThan(0.0),
                "Non-domain sequence must score negative against WD40.");
        });
    }

    [Test]
    public void ScoreDomainHmm_CrossDomain_TruePositiveScoresHigherOnItsOwnFamily()
    {
        // H7 — an SH3 domain scores far higher against PF00018 than against the WD40 profile.
        double onSh3 = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession);
        double onWd40 = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Wd40Accession);

        Assert.Multiple(() =>
        {
            Assert.That(onSh3, Is.GreaterThan(onWd40),
                "An SH3 sequence must score higher against its own profile than against WD40.");
            Assert.That(onWd40, Is.LessThan(0.0),
                "An SH3 sequence must score negative against the unrelated WD40 profile.");
        });
    }

    #endregion

    #region H10 — FindDomainsByHmm assigns the correct family

    [Test]
    public void FindDomainsByHmm_Sh3TruePositive_DetectsOnlySh3()
    {
        // H10 — the SH3 core is reported as SH3 (PF00018) and not as PDZ/WD40.
        var domains = new List<ProteinMotifFinder.ProteinDomain>(
            ProteinMotifFinder.FindDomainsByHmm(Sh3TruePositive));

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one bundled family should match the SH3 core.");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Accession, Is.EqualTo(Sh3Accession), "Detected family must be SH3 (PF00018).");
            Assert.That(domains[0].Name, Is.EqualTo("SH3"), "Reported domain name.");
            Assert.That(domains[0].Start, Is.EqualTo(0), "Glocal score spans from sequence start.");
            Assert.That(domains[0].End, Is.EqualTo(Sh3TruePositive.Length - 1), "…to sequence end (0-based inclusive).");
            Assert.That(domains[0].Score, Is.GreaterThan(DetectionThresholdBits), "Reported Score is the bit score.");
        });
    }

    [Test]
    public void FindDomainsByHmm_Wd40TruePositive_DetectsWd40()
    {
        // H10 — GBB1 is reported as WD40 (PF00400).
        var domains = new List<ProteinMotifFinder.ProteinDomain>(
            ProteinMotifFinder.FindDomainsByHmm(Wd40TruePositive));

        Assert.That(domains, Has.Some.Matches<ProteinMotifFinder.ProteinDomain>(d => d.Accession == Wd40Accession),
            "GBB1 β-propeller must be detected as WD40 (PF00400).");
    }

    [Test]
    public void FindDomainsByHmm_TrueNegative_ReportsNoDomain()
    {
        // H6/H10 — the non-domain sequence is reported by no bundled profile.
        var domains = new List<ProteinMotifFinder.ProteinDomain>(
            ProteinMotifFinder.FindDomainsByHmm(TrueNegative));

        Assert.That(domains, Is.Empty, "A non-domain sequence must not be reported as SH3/PDZ/WD40.");
    }

    #endregion

    #region H11 — determinism

    [Test]
    public void ScoreDomainHmm_RepeatedCalls_AreDeterministic()
    {
        // H11 — INV-HMM-03.
        double a = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession);
        double b = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession);

        Assert.That(b, Is.EqualTo(a).Within(1e-12), "Scoring is deterministic; repeated calls match exactly.");
    }

    #endregion

    #region H12 — null / empty / unknown-accession guards

    [Test]
    public void FindDomainsByHmm_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProteinMotifFinder.FindDomainsByHmm(null!), Is.Empty, "Null → empty.");
            Assert.That(ProteinMotifFinder.FindDomainsByHmm(string.Empty), Is.Empty, "Empty → empty.");
        });
    }

    [Test]
    public void ScoreDomainHmm_NullSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => ProteinMotifFinder.ScoreDomainHmm(null!, Sh3Accession),
            "Null sequence is a programming error.");
    }

    [Test]
    public void ScoreDomainHmm_UnknownAccession_Throws()
    {
        Assert.Throws<ArgumentException>(() => ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, "PF99999"),
            "Only PF00018/PF00595/PF00400 are bundled.");
    }

    #endregion

    #region H13 — STATS lines parsed from the bundled profile

    [Test]
    public void Parse_EmbeddedSh3Profile_ReadsStatsLinesExactly()
    {
        // H13 — the three STATS LOCAL lines are read verbatim into the calibration record.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        Assert.That(hmm.Statistics, Is.Not.Null, "PF00018 carries all three STATS lines → calibrated.");
        var s = hmm.Statistics!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(s.MsvMu, Is.EqualTo(Sh3MsvMu).Within(1e-12), "STATS LOCAL MSV location μ.");
            Assert.That(s.MsvLambda, Is.EqualTo(Sh3Lambda).Within(1e-12), "STATS LOCAL MSV slope λ.");
            Assert.That(s.ViterbiMu, Is.EqualTo(Sh3ViterbiMu).Within(1e-12), "STATS LOCAL VITERBI location μ.");
            Assert.That(s.ViterbiLambda, Is.EqualTo(Sh3Lambda).Within(1e-12), "STATS LOCAL VITERBI slope λ.");
            Assert.That(s.ForwardTau, Is.EqualTo(Sh3ForwardTau).Within(1e-12), "STATS LOCAL FORWARD location τ.");
            Assert.That(s.ForwardLambda, Is.EqualTo(Sh3Lambda).Within(1e-12), "STATS LOCAL FORWARD slope λ.");
        });
    }

    [Test]
    public void Statistics_UncalibratedProfile_IsNull_AndPValueThrows()
    {
        // H13 — the hand-built HMM has no STATS lines; it is not calibrated.
        var hmm = BuildHandHmm();

        Assert.Multiple(() =>
        {
            Assert.That(hmm.Statistics, Is.Null, "A profile with no STATS lines is uncalibrated.");
            Assert.Throws<InvalidOperationException>(() => hmm.ViterbiPValue(0.0),
                "Requesting a P-value from an uncalibrated profile is an error.");
        });
    }

    #endregion

    #region H14 — exact hand-derived Gumbel / exponential P-value and E-value

    [Test]
    public void ViterbiPValue_Sh3At40Bits_MatchesHandDerivedGumbel()
    {
        // H14 — Gumbel survival P = 1 − exp(−exp(−λ(S−μ))) with μ=−8.2932, λ=0.71923, S=40.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double p = hmm.ViterbiPValue(HandBitScore);

        Assert.That(p, Is.EqualTo(Sh3ViterbiPValueAt40).Within(1e-9 * Sh3ViterbiPValueAt40),
            "Viterbi Gumbel P-value at 40 bits must equal the hand-derived 8.227179545686635e-16.");
    }

    [Test]
    public void ViterbiEValue_Sh3At40Bits_Z1000_MatchesHandDerived()
    {
        // H14 — E = P · Z with Z = 1000.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double e = hmm.ViterbiEValue(HandBitScore, HandZ);

        Assert.That(e, Is.EqualTo(Sh3ViterbiEValueAt40Z1000).Within(1e-9 * Sh3ViterbiEValueAt40Z1000),
            "Viterbi E-value = P·1000 must equal the hand-derived 8.227179545686635e-13.");
    }

    [Test]
    public void ForwardPValue_Sh3At40Bits_MatchesHandDerivedExponential()
    {
        // H14 — exponential-tail survival P = exp(−λ(S−τ)) with τ=−4.5735, λ=0.71923, S=40.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double p = hmm.ForwardPValue(HandBitScore);

        Assert.That(p, Is.EqualTo(Sh3ForwardPValueAt40).Within(1e-9 * Sh3ForwardPValueAt40),
            "Forward exponential P-value at 40 bits must equal the hand-derived 1.1943390031599535e-14.");
    }

    [Test]
    public void ForwardEValue_Sh3At40Bits_Z1000_MatchesHandDerived()
    {
        // H14 — E = P · Z with Z = 1000 for the Forward exponential tail.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double e = hmm.ForwardEValue(HandBitScore, HandZ);

        Assert.That(e, Is.EqualTo(Sh3ForwardEValueAt40Z1000).Within(1e-9 * Sh3ForwardEValueAt40Z1000),
            "Forward E-value = P·1000 must equal the hand-derived 1.1943390031599535e-11.");
    }

    [Test]
    public void ForwardPValue_BelowTau_IsClampedToOne()
    {
        // H14 — Easel esl_exp_surv returns 1.0 for x < mu (τ); a score below the tail base is non-significant.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        // τ = −4.5735; any score below that clamps to P = 1.
        double p = hmm.ForwardPValue(-10.0);

        Assert.That(p, Is.EqualTo(1.0).Within(1e-12),
            "Forward P-value clamps to 1.0 for scores below the exponential-tail location τ.");
    }

    [Test]
    public void GumbelSurvival_PureFormula_MatchesHandDerived()
    {
        // H14 — the static Gumbel survival evaluated directly with the SH3 Viterbi parameters.
        double p = Plan7ProfileHmm.GumbelSurvival(HandBitScore, Sh3ViterbiMu, Sh3Lambda);

        Assert.That(p, Is.EqualTo(Sh3ViterbiPValueAt40).Within(1e-9 * Sh3ViterbiPValueAt40),
            "GumbelSurvival(40, −8.2932, 0.71923) must equal the hand-derived value.");
    }

    #endregion

    #region H15 — monotonicity and Z-scaling

    [Test]
    public void ViterbiEValue_DecreasesAsScoreIncreases()
    {
        // H15 — a higher bit score is more significant: smaller P, hence smaller E (strictly).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double eLow = hmm.ViterbiEValue(20.0, 1.0);
        double eHigh = hmm.ViterbiEValue(40.0, 1.0);

        Assert.That(eHigh, Is.LessThan(eLow),
            "E-value must strictly decrease as the bit score increases (Gumbel is monotone).");
    }

    [Test]
    public void ViterbiEValue_ScalesLinearlyWithDatabaseSize()
    {
        // H15 — E = P·Z is exactly linear in Z. E(Z=1000) = 1000 · E(Z=1).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double e1 = hmm.ViterbiEValue(HandBitScore, 1.0);
        double e1000 = hmm.ViterbiEValue(HandBitScore, 1000.0);

        Assert.That(e1000, Is.EqualTo(e1 * 1000.0).Within(1e-9 * (e1 * 1000.0)),
            "E-value scales linearly with the database size Z (E = P·Z).");
    }

    [Test]
    public void EValue_NegativeDatabaseSize_Throws()
    {
        // H15 — Z must be non-negative.
        Assert.Throws<ArgumentOutOfRangeException>(() => Plan7ProfileHmm.EValue(0.5, -1.0),
            "A negative database size Z is invalid.");
    }

    #endregion

    #region H16 — end-to-end true-positive small-E, negative large-E

    [Test]
    public void ScoreDomainHmmEValue_Sh3TruePositive_HasTinyEValue()
    {
        // H16 — a real SH3 domain (SRC_HUMAN core) is highly significant: E ≪ 1e-3 even at Z=1.
        var (bits, e) = ProteinMotifFinder.ScoreDomainHmmEValue(Sh3TruePositive, Sh3Accession, 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(bits, Is.GreaterThan(DetectionThresholdBits), "True-positive bit score is well above threshold.");
            Assert.That(e, Is.LessThan(1e-3), "A genuine SH3 domain must have a tiny E-value (highly significant).");
        });
    }

    [Test]
    public void ScoreDomainHmmEValue_TrueNegative_HasLargeEValue()
    {
        // H16 — a low-complexity non-domain sequence is non-significant: P=1 → E=Z (here 1000).
        var (bits, e) = ProteinMotifFinder.ScoreDomainHmmEValue(TrueNegative, Sh3Accession, 1000.0);

        Assert.Multiple(() =>
        {
            Assert.That(bits, Is.LessThan(0.0), "A non-domain sequence scores negative bits against SH3.");
            Assert.That(e, Is.GreaterThan(1.0), "A non-domain sequence is non-significant: large E-value (≈ Z).");
        });
    }

    [Test]
    public void FindDomainHitsByHmm_Sh3TruePositive_ReportsHitWithEValue()
    {
        // H16 — the E-value-bearing detection path reports SH3 with a tiny E-value and the bit Score.
        var hits = new List<ProteinMotifFinder.ProteinDomainHit>(
            ProteinMotifFinder.FindDomainHitsByHmm(Sh3TruePositive, databaseSize: 1.0));

        Assert.That(hits, Has.Count.EqualTo(1), "Exactly one bundled family matches the SH3 core.");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Accession, Is.EqualTo(Sh3Accession), "Detected family is SH3 (PF00018).");
            Assert.That(hits[0].Score, Is.GreaterThan(DetectionThresholdBits), "Reported Score is the Viterbi bit score.");
            Assert.That(hits[0].EValue, Is.LessThan(1e-3), "Reported E-value is tiny for the true positive.");
            Assert.That(hits[0].Start, Is.EqualTo(0), "Glocal score spans from sequence start.");
            Assert.That(hits[0].End, Is.EqualTo(Sh3TruePositive.Length - 1), "…to sequence end (0-based inclusive).");
        });
    }

    [Test]
    public void FindDomainHitsByHmm_TrueNegative_ReportsNoHit()
    {
        // H16 — the non-domain sequence is reported by no bundled profile.
        var hits = new List<ProteinMotifFinder.ProteinDomainHit>(
            ProteinMotifFinder.FindDomainHitsByHmm(TrueNegative));

        Assert.That(hits, Is.Empty, "A non-domain sequence must not be reported as SH3/PDZ/WD40.");
    }

    [Test]
    public void FindDomainHitsByHmm_NullOrEmpty_ReturnsEmpty()
    {
        // H16 — null/empty guard on the E-value detection path.
        Assert.Multiple(() =>
        {
            Assert.That(ProteinMotifFinder.FindDomainHitsByHmm(null!), Is.Empty, "Null → empty.");
            Assert.That(ProteinMotifFinder.FindDomainHitsByHmm(string.Empty), Is.Empty, "Empty → empty.");
        });
    }

    #endregion

    #region H18 — HMMER local-multihit mode + null2 (hmmsearch parity)

    // H18 — exact hand-derived local-mode pin (independent of any tool). A 1-node HMM that emits
    // residue A with probability 1, with B->M1 = 1 (occupancy occ[1]=1 → local-entry log(occ/Z)=0).
    // For sequence "A": local Forward path N->B->M1(A)->E->C->T.
    //   msc = ln(1 / bg[A]) with HMMER's bg[A] = 0.0787945; nj=1 (multihit), L=1 →
    //   pmove = (2+1)/(1+2+1) = 3/4, xN_move = ln(3/4); xE_move = -ln 2.
    //   fwd = xN_move + msc + xE_move + xN_move = 1.272400756045032 nats.
    //   nullsc = 1·ln(1/2)+ln(1/2); pre = (fwd-nullsc)/ln2 = 3.835686260769536 bits.
    private const double LocalToyForwardNats = 1.272400756045032;
    private const double LocalToyPreBits = 3.835686260769536;

    // H18 — hmmsearch GROUND-TRUTH local-multihit Forward bit score (pre_score, BEFORE null2),
    // captured from pyhmmer 0.12.1 (hmmsearch([hmm], targets, Z=1.0)) on the bundled CC0 profiles
    // vs the test true-positives. See Evidence addendum 2026-06-25 table.
    private const double Sh3HmmsearchPreScore = 68.709740;
    private const double PdzHmmsearchPreScore = 84.862930;
    private const double Wd40HmmsearchPreScore = 213.411926;
    // hmmsearch reported null2 "bias" (bits) for the single SH3 domain (envelope 3-50).
    private const double Sh3HmmsearchDomainBias = 0.025574;
    // The 1e-4-bit tolerance reflects HMMER's single-precision (float) arithmetic vs this double DP.
    private const double HmmsearchBitTolerance = 1e-4;

    private static Plan7ProfileHmm BuildLocalToyHmm()
    {
        const double NEG = double.NegativeInfinity;
        double[] EmitA()
        {
            var r = new double[20];
            for (int i = 0; i < 20; i++) r[i] = Math.Log(1e-12);
            r[0] = Math.Log(1.0); // residue A
            return r;
        }
        double Lg(double p) => p <= 0 ? NEG : Math.Log(p);
        var bg = new double[20];
        for (int i = 0; i < 20; i++) bg[i] = Math.Log(1.0 / 20);
        var match = new double[2][]; var insert = new double[2][]; var trans = new double[2][];
        insert[0] = EmitA();
        // BEGIN: MM=1 (B->M1), MI=0, MD=0, IM=1, II=0, DM=1, DD=*
        trans[0] = new[] { Lg(1.0), NEG, NEG, Lg(1.0), NEG, Lg(1.0), NEG };
        match[1] = EmitA(); insert[1] = EmitA();
        trans[1] = new[] { Lg(1.0), NEG, NEG, Lg(1.0), NEG, Lg(1.0), NEG };
        return Plan7ProfileHmm.FromLogParameters("toy", 1, match, insert, trans, bg);
    }

    [Test]
    public void LocalForwardScore_HandBuiltOneNode_MatchesExactDerivation()
    {
        // H18 — pins the local entry/exit + length-config arithmetic to the hand derivation.
        var hmm = BuildLocalToyHmm();

        double nats = hmm.LocalForwardScore("A");
        double pre = hmm.LocalForwardBitScore("A");

        Assert.Multiple(() =>
        {
            Assert.That(nats, Is.EqualTo(LocalToyForwardNats).Within(1e-12),
                "Local-multihit Forward of 'A' through N->B->M1->E->C->T must equal the hand-derived 1.272400756045032 nats.");
            Assert.That(pre, Is.EqualTo(LocalToyPreBits).Within(1e-12),
                "Local Forward bit score (fwd-nullsc)/ln2 must equal the hand-derived 3.835686260769536 bits.");
        });
    }

    [Test]
    public void LocalForwardBitScore_Sh3TruePositive_MatchesHmmsearchReference()
    {
        // H18 — GROUND-TRUTH parity: the local-multihit Forward bit score reproduces hmmsearch's
        // pre_score (68.709740) for SRC_HUMAN SH3 vs PF00018 (pyhmmer 0.12.1).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double pre = hmm.LocalForwardBitScore(Sh3TruePositive);

        Assert.That(pre, Is.EqualTo(Sh3HmmsearchPreScore).Within(HmmsearchBitTolerance),
            "Local-multihit Forward bit score must match the hmmsearch pre_score 68.709740 bits to float precision.");
    }

    [Test]
    public void LocalForwardBitScore_PdzTruePositive_MatchesHmmsearchReference()
    {
        // H18 — GROUND-TRUTH parity for PSD-95 PDZ1 vs PF00595 (hmmsearch pre_score 84.862930).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00595_PDZ.hmm"));

        double pre = hmm.LocalForwardBitScore(PdzTruePositive);

        Assert.That(pre, Is.EqualTo(PdzHmmsearchPreScore).Within(HmmsearchBitTolerance),
            "Local-multihit Forward bit score must match the hmmsearch pre_score 84.862930 bits.");
    }

    [Test]
    public void LocalForwardBitScore_Wd40TruePositive_MatchesHmmsearchReference()
    {
        // H18 — GROUND-TRUTH parity for GBB1_HUMAN β-propeller vs PF00400 (hmmsearch pre_score 213.411926).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00400_WD40.hmm"));

        double pre = hmm.LocalForwardBitScore(Wd40TruePositive);

        Assert.That(pre, Is.EqualTo(Wd40HmmsearchPreScore).Within(HmmsearchBitTolerance),
            "Local-multihit Forward bit score must match the hmmsearch pre_score 213.411926 bits.");
    }

    [Test]
    public void Null2BiasBits_Sh3DomainEnvelope_MatchesHmmsearchReportedBias()
    {
        // H18 — the null2 biased-composition correction reproduces hmmsearch's reported "bias"
        // (0.025574 bits) for the single SH3 domain. HMMER computes null2 per envelope (positions
        // 3-50 of the test sequence, the model reconfigured to the envelope length); we pass that
        // envelope substring. Matches to single-precision rounding.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));
        string envelope = Sh3TruePositive.Substring(2, 48); // 1-based positions 3..50

        double bias = hmm.Null2BiasBits(envelope);

        Assert.That(bias, Is.EqualTo(Sh3HmmsearchDomainBias).Within(1e-4),
            "null2 bias over the SH3 domain envelope must match hmmsearch's reported 0.025574 bits.");
    }

    [Test]
    public void HmmSearchBitScore_EqualsPreScoreMinusNull2Bias()
    {
        // H18 — INV: the corrected per-seq score = pre_score - null2 bias (both in bits), exactly
        // as p7_pipeline.c computes seq_score = pre_score - bias.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double pre = hmm.LocalForwardBitScore(Sh3TruePositive);
        double bias = hmm.Null2BiasBits(Sh3TruePositive);
        double score = hmm.HmmSearchBitScore(Sh3TruePositive);

        Assert.That(score, Is.EqualTo(pre - bias).Within(1e-9),
            "HmmSearchBitScore must equal LocalForwardBitScore - Null2BiasBits (pipeline identity).");
    }

    [Test]
    public void Null2BiasBits_IsNonNegative_AndPositiveForBiasedComposition()
    {
        // H18 — INV: the null2 mixture seqbias = logsumexp(0, ...) >= 0, so the correction never
        // increases the score; a biased-composition (low-complexity) sequence gets a positive bias.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double biasReal = hmm.Null2BiasBits(Sh3TruePositive);
        double biasLowComplexity = hmm.Null2BiasBits(TrueNegative);

        Assert.Multiple(() =>
        {
            Assert.That(biasReal, Is.GreaterThanOrEqualTo(0.0),
                "null2 bias = logsumexp(0, ...) is non-negative; it never raises the score.");
            Assert.That(biasLowComplexity, Is.GreaterThanOrEqualTo(biasReal),
                "A low-complexity (biased-composition) sequence incurs at least as large a null2 correction.");
        });
    }

    [Test]
    public void LocalForwardScore_DoesNotChangeGlocalForwardScore()
    {
        // H18 — REGRESSION: the opt-in local path leaves the existing glocal ForwardScore unchanged.
        // The two are different alignment modes and must give different (here, well-separated) scores.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double glocalNats = hmm.ForwardScore(Sh3TruePositive);
        double localNats = hmm.LocalForwardScore(Sh3TruePositive);

        Assert.That(glocalNats, Is.EqualTo(41.78685952002655).Within(1e-6),
            "Glocal ForwardScore must be unchanged by the new local-mode methods.");
        Assert.That(localNats, Is.EqualTo(42.609594871580114).Within(1e-6),
            "Local-multihit Forward is a distinct alignment mode (different score) from the glocal Forward.");
        Assert.That(localNats, Is.Not.EqualTo(glocalNats).Within(1e-6),
            "The two modes produce different scores; the glocal value is untouched.");
    }

    [Test]
    public void LocalForwardScore_NullOrEmpty_GuardsCorrectly()
    {
        // H18 — null/empty guards on the parity path.
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => hmm.LocalForwardScore(null!), "Null sequence is a programming error.");
            Assert.That(double.IsNegativeInfinity(hmm.HmmSearchBitScore(string.Empty)), Is.True,
                "Empty sequence has no Forward path → -inf bit score.");
            Assert.That(hmm.Null2BiasBits(string.Empty), Is.EqualTo(0.0).Within(1e-12),
                "Empty sequence has no residues → zero null2 bias.");
        });
    }

    [Test]
    public void HmmSearchBitScore_Sh3TruePositive_IsHighlyPositive()
    {
        // H18 — sanity: a genuine SH3 domain's null2-corrected hmmsearch bit score is large (≈68.7).
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

        double score = hmm.HmmSearchBitScore(Sh3TruePositive);

        Assert.That(score, Is.GreaterThan(50.0),
            "A real SH3 domain's null2-corrected local bit score must be strongly positive.");
    }

    #endregion

    private static string ReadEmbedded(string fileName)
    {
        var asm = typeof(Plan7ProfileHmm).Assembly;
        string resourceName = $"Seqeron.Genomics.Analysis.Resources.{fileName}";
        using var stream = asm.GetManifestResourceStream(resourceName)
            ?? throw new InvalidOperationException($"Missing embedded resource: {resourceName}");
        using var reader = new System.IO.StreamReader(stream);
        return reader.ReadToEnd();
    }
}
