// PROTMOTIF-DOMAIN-001 — Plan7 profile-HMM domain detection (SH3 / PDZ / WD40), opt-in.
// Evidence: docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md  (addendum 2026-06-25)
// TestSpec: tests/TestSpecs/PROTMOTIF-DOMAIN-001.md          (addendum 2026-06-25)
// Source:   Durbin et al. (1998) Biological Sequence Analysis §5.4 (Viterbi recurrence);
//           HMMER User's Guide v3.4 (Eddy 2023) — HMMER3/f file format; Pfam PF00018/PF00595/PF00400 (CC0).
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
