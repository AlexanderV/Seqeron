// SEQ-PI-001 — Isoelectric point: EMBOSS iep / Biopython window conventions and iep qualifiers
// (review 2026-09, FIN-B03 WP5, B03 audit rows A9/A10, entries F43/F44).
// Evidence: docs/Evidence/SEQ-PI-001-Evidence.md  TestSpec: tests/TestSpecs/SEQ-PI-001.md
// Sources (reference implementations executed 2026-10-14):
//   • EMBOSS 6.6.0 `iep` binary (Debian emboss 6.6.0): `iep -sequence q.fa -outfile stdout -auto
//     [-notermini] [-lysinemodified N] [-disulphides N]` — "Isoelectric Point = x.xxxx" / "None".
//     Semantics from nucleus/embiep.c (embIepCompC count adjustment, embIepPhConverge [1,14]
//     bisection to 1e-4 returning the lower bound; "None" when the charge has the same strict sign
//     at pH 1 and 14), https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/nucleus/embiep.c
//   • Biopython 1.88 Bio.SeqUtils.IsoelectricPoint.IsoelectricPoint(seq).pi() (defaults pH=7.775,
//     min_=4.05, max_=12, bracket 1e-4, no range extension) — printed with repr().

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_IsoelectricPointOptions_Tests
{
    // iep prints the pI with 4 decimals (%-6.4lf).
    private const double IepPrint = 0.00005;
    // Biopython's bisection path is reproduced step by step; values agree to the last digits.
    private const double BiopythonExact = 1e-12;

    private const string LacIEcoli =
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVATSSLALHAPSQIVAAIKSRADQLG" +
        "ASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIA" +
        "LLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDD" +
        "TEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLAPNTQTASPRALADSLMQLARQVSRLESGQ";

    private static readonly string Poly40R = new('R', 40);

    private static SequenceStatistics.IsoelectricPointOptions Emboss(
        bool termini = true, int lysinesModified = 0, int disulphides = 0) => new()
    {
        Scale = SequenceStatistics.PkaScale.Emboss,
        Window = SequenceStatistics.IsoelectricPointWindow.Emboss,
        Termini = termini,
        LysinesModified = lysinesModified,
        Disulphides = disulphides,
    };

    private static readonly SequenceStatistics.IsoelectricPointOptions Biopython = new()
    {
        Scale = SequenceStatistics.PkaScale.Bjellqvist,
        Window = SequenceStatistics.IsoelectricPointWindow.Biopython,
    };

    #region EMBOSS window ([1,14], "None")

    // iep 6.6.0 default qualifiers.
    [TestCase("A", 5.5500)]
    [TestCase("GGG", 5.5500)]
    [TestCase("DDDD", 3.2279)]
    [TestCase("DDDDDDDD", 2.9549)]
    [TestCase("KKKK", 11.2772)]
    [TestCase("KKKKKRRRR", 12.9956)]
    [TestCase("RRRRRRRR", 13.3450)]
    [TestCase("PETER", 4.2577)]
    [TestCase("INGAR", 10.0000)]
    [TestCase("CCKK", 8.5598)]
    [TestCase("ACCKCKKA", 8.8115)]
    [TestCase("KDDK", 5.9023)]
    [TestCase("MAEGEITTFT", 3.6135)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 6.9681)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.5678)]
    public void EmbossWindow_DefaultQualifiers_MatchesIep(string seq, double iep)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Emboss(), out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), $"iep {seq}");
        Assert.That(Math.Round(pi, 2), Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint(seq)).Within(1e-9),
            "inside [1,14] the EMBOSS window agrees with the library pI at 2 dp");
    }

    [Test]
    public void EmbossWindow_LacIEcoli_MatchesIepDocExample()
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(LacIEcoli, Emboss(), out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(6.8385).Within(IepPrint), "EMBOSS iep documentation: LACI_ECOLI 6.8385");
    }

    // Swiss-Prot IFNA2_HUMAN (P01563) from the EMBOSS test database (test/swiss/seq.dat) — the iep
    // documentation usage example 2 (`-disulphides 2 -lysinemodified 2` = 5.7240), plus the other
    // qualifier combinations executed with iep 6.6.0; Biopython 1.88 pi() = 6.32307643890381.
    private const string Ifna2Human =
        "MALTFALLVALLVLSCKSSCSVGCDLPQTHSLGSRRTLMLLAQMRKISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTK" +
        "DSSAAWDETLLDKFYTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE";

    [TestCase(true, 0, 0, 6.6891)]
    [TestCase(true, 2, 2, 5.7240)]   // EMBOSS iep documentation example
    [TestCase(false, 0, 0, 6.7573)]
    [TestCase(false, 2, 2, 5.7289)]
    public void EmbossWindow_Ifna2Human_MatchesIep(bool termini, int disulphides, int lysines, double iep)
    {
        var options = Emboss(termini: termini, lysinesModified: lysines, disulphides: disulphides);
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(Ifna2Human, options, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), "iep 6.6.0 IFNA2_HUMAN");
    }

    [Test]
    public void BiopythonWindow_Ifna2Human_MatchesBiopythonPi()
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(Ifna2Human, Biopython, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(6.32307643890381).Within(BiopythonExact));
    }

    // 40 R: charge > 0 at both pH 1 and pH 14 → iep "Isoelectric Point = None"; the library window
    // reports the [0,14] edge instead.
    [Test]
    public void EmbossWindow_Poly40Arginine_NoPi()
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(Poly40R, Emboss(), out double pi), Is.False);
        Assert.That(pi, Is.NaN);
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(Poly40R, null, out double lib), Is.True);
        Assert.That(lib, Is.EqualTo(14.0).Within(1e-9), "library convention: window edge");
    }

    // iep -notermini.
    [TestCase("A", 7.5000)]      // no charged group: zero charge everywhere → first midpoint 7.5
    [TestCase("GGG", 7.5000)]
    [TestCase("PETER", 4.0999)]
    [TestCase("CCKK", 9.6500)]
    [TestCase("ACCKCKKA", 9.6500)]
    [TestCase("KDDK", 7.3499)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 7.4939)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 10.5115)]
    public void EmbossWindow_NoTermini_MatchesIep(string seq, double iep)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Emboss(termini: false), out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), $"iep -notermini {seq}");
    }

    // iep -notermini → "None" (only acidic or only basic side chains).
    [TestCase("DDDD")]
    [TestCase("DDDDDDDD")]
    [TestCase("KKKK")]
    [TestCase("KKKKKRRRR")]
    [TestCase("RRRRRRRR")]
    [TestCase("INGAR")]
    public void EmbossWindow_NoTermini_NoPi(string seq)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Emboss(termini: false), out double pi), Is.False,
            $"iep -notermini {seq}: Isoelectric Point = None");
        Assert.That(pi, Is.NaN);
    }

    // iep -lysinemodified 2 (fewer than 2 K → all K removed, embIepCompC).
    [TestCase("KKKK", 10.8004)]
    [TestCase("KKKKKRRRR", 12.9883)]
    [TestCase("CCKK", 5.5096)]
    [TestCase("ACCKCKKA", 7.6820)]
    [TestCase("KDDK", 3.4918)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 5.4304)]  // 1 K < 2
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 7.8001)]
    [TestCase("PETER", 4.2577)]                 // no K: unchanged
    public void EmbossWindow_LysinesModified2_MatchesIep(string seq, double iep)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Emboss(lysinesModified: 2), out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), $"iep -lysinemodified 2 {seq}");
    }

    // iep -disulphides 1 (fewer than 2 C → all C removed, embIepCompC).
    [TestCase("CCKK", 10.8004)]
    [TestCase("ACCKCKKA", 10.5059)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 7.0017)]  // 1 C < 2
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 11.1011)]
    [TestCase("KDDK", 5.9023)]                  // no C: unchanged
    public void EmbossWindow_Disulphides1_MatchesIep(string seq, double iep)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Emboss(disulphides: 1), out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), $"iep -disulphides 1 {seq}");
    }

    // iep -lysinemodified 1 -disulphides 1 -notermini.
    [TestCase("ACCKCKKA", 10.8043)]
    [TestCase("KDDK", 3.9000)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 5.4144)]
    [TestCase("A", 7.5000)]
    public void EmbossWindow_Combined_MatchesIep(string seq, double iep)
    {
        var options = Emboss(termini: false, lysinesModified: 1, disulphides: 1);
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, options, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(iep).Within(IepPrint), $"iep -lysinemodified 1 -disulphides 1 -notermini {seq}");
    }

    [TestCase("CCKK")]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC")]
    public void EmbossWindow_Combined_NoPi(string seq)
    {
        var options = Emboss(termini: false, lysinesModified: 1, disulphides: 1);
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, options, out _), Is.False,
            $"iep -lysinemodified 1 -disulphides 1 -notermini {seq}: None");
    }

    #endregion

    #region Biopython window ([4.05,12] from 7.775)

    // Biopython 1.88 IsoelectricPoint(seq).pi() (== ProteinAnalysis(seq).isoelectric_point()).
    [TestCase("A", 5.57001667022705)]
    [TestCase("AG", 5.57001667022705)]
    [TestCase("GGG", 5.525000190734862)]
    [TestCase("D", 4.299381065368652)]
    [TestCase("E", 4.599263572692871)]
    [TestCase("DE", 4.0851549148559565)]
    [TestCase("DDDD", 4.0500284194946286)]          // true root 3.52: clamped at the window edge
    [TestCase("DDDDDDDD", 4.0500284194946286)]
    [TestCase("EEEE", 4.0500284194946286)]
    [TestCase("MAEGEITTFT", 4.0500284194946286)]
    [TestCase("KKKK", 10.477740287780762)]
    [TestCase("KKKKKRRRR", 11.999967765808105)]     // basic edge
    [TestCase("RRRRRRRR", 11.999967765808105)]      // true root 12.85
    [TestCase("PETER", 4.532079887390137)]
    [TestCase("peter", 4.532079887390137)]
    [TestCase("INGAR", 9.750021171569824)]
    [TestCase("AKD", 6.13152904510498)]
    [TestCase("CCKK", 8.9037446975708)]
    [TestCase("ACCKCKKA", 9.109205436706542)]
    [TestCase("KDDK", 5.959761619567872)]
    [TestCase("ACDEFGHIKLMNPQRSTVWY", 6.784552192687988)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", 9.390158653259277)]
    public void BiopythonWindow_Bjellqvist_MatchesBiopythonPi(string seq, double biopython)
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, Biopython, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(biopython).Within(BiopythonExact), $"Biopython 1.88 pi() {seq}");
    }

    [Test]
    public void BiopythonWindow_Poly40Arginine_ReturnsUpperEdge()
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(Poly40R, Biopython, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(11.999967765808105).Within(BiopythonExact), "Biopython 1.88 pi() of R×40");
    }

    // Library window keeps the true root where Biopython clamps.
    [TestCase("DDDD", 3.52)]
    [TestCase("RRRRRRRR", 12.85)]
    public void DefaultWindow_Bjellqvist_KeepsTrueRootOutsideBiopythonWindow(string seq, double root)
    {
        var options = new SequenceStatistics.IsoelectricPointOptions { Scale = SequenceStatistics.PkaScale.Bjellqvist };
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, options, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(root).Within(1e-9));
    }

    #endregion

    #region Defaults, library window with qualifiers, guards

    [TestCase("PETER")]
    [TestCase("KKKK")]
    [TestCase("ACDEFGHIKLMNPQRSTVWY")]
    [TestCase("BZ")]
    [TestCase(LacIEcoli)]
    public void DefaultOptions_EqualExistingOverloads(string seq)
    {
        foreach (var scale in new[] { SequenceStatistics.PkaScale.Emboss, SequenceStatistics.PkaScale.Bjellqvist })
        {
            var options = new SequenceStatistics.IsoelectricPointOptions { Scale = scale };
            Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, options, out double pi), Is.True);
            Assert.That(pi, Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint(seq, scale)));
            foreach (double pH in new[] { 1.0, 4.5, 7.0, 9.3, 14.0 })
                Assert.That(SequenceStatistics.CalculateNetCharge(seq, pH, options),
                    Is.EqualTo(SequenceStatistics.CalculateNetCharge(seq, pH, scale)));
        }
    }

    // Library window ([0,14], 2 dp) with iep qualifiers = iep value rounded to 2 dp.
    [TestCase("KKKK", 0, 2, 10.80)]       // iep -lysinemodified 2 10.8004
    [TestCase("CCKK", 1, 0, 10.80)]       // iep -disulphides 1 10.8004
    [TestCase("ACCKCKKA", 1, 0, 10.51)]   // iep 10.5059
    public void DefaultWindow_Qualifiers_MatchIepRounded(string seq, int disulphides, int lysines, double expected)
    {
        var options = new SequenceStatistics.IsoelectricPointOptions { Disulphides = disulphides, LysinesModified = lysines };
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(seq, options, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void DefaultWindow_NoTermini_MatchesIepRounded()
    {
        var options = new SequenceStatistics.IsoelectricPointOptions { Termini = false };
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint("PETER", options, out double pi), Is.True);
        Assert.That(pi, Is.EqualTo(4.10).Within(1e-9), "iep -notermini PETER 4.0999");
    }

    // -notermini: side-chain charge only; zero for a peptide without ionizable side chains.
    [Test]
    public void NetCharge_NoTermini_ZeroWithoutIonizableSideChains_AndSideChainsOnlyOtherwise()
    {
        var noTermini = new SequenceStatistics.IsoelectricPointOptions { Termini = false };
        Assert.That(SequenceStatistics.CalculateNetCharge("AGLV", 3.0, noTermini), Is.EqualTo(0.0));
        double k = 1.0 / (1.0 + Math.Pow(10, 7.0 - 10.8)); // EMBOSS Epk.dat K 10.8
        Assert.That(SequenceStatistics.CalculateNetCharge("AKA", 7.0, noTermini), Is.EqualTo(k).Within(1e-12));
    }

    // Bjellqvist terminal-residue pKs are preserved: a modified N-terminal K is neutralised in place
    // (charge of "GAD"), not deleted (which would expose the N-terminal A, pK 7.59).
    [Test]
    public void NetCharge_Bjellqvist_ModifiedTerminalLysine_KeepsTerminalPk()
    {
        var options = new SequenceStatistics.IsoelectricPointOptions
        {
            Scale = SequenceStatistics.PkaScale.Bjellqvist,
            LysinesModified = 1,
        };
        double charge = SequenceStatistics.CalculateNetCharge("KAD", 7.0, options);
        Assert.That(charge, Is.EqualTo(SequenceStatistics.CalculateNetCharge("GAD", 7.0, SequenceStatistics.PkaScale.Bjellqvist)));
        Assert.That(charge, Is.Not.EqualTo(SequenceStatistics.CalculateNetCharge("AD", 7.0, SequenceStatistics.PkaScale.Bjellqvist)));
    }

    [Test]
    public void EmptyInput_WindowConventions()
    {
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint("", null, out double lib), Is.True);
        Assert.That(lib, Is.EqualTo(7.0), "library sentinel");
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint("", Emboss(), out _), Is.False);
        Assert.That(SequenceStatistics.TryCalculateIsoelectricPoint(null!, Biopython, out _), Is.False,
            "Biopython raises IndexError on an empty sequence");
        Assert.That(SequenceStatistics.CalculateNetCharge("", 7.0, Emboss(termini: false)), Is.EqualTo(0.0));
    }

    [Test]
    public void InvalidOptions_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.TryCalculateIsoelectricPoint("K", new() { LysinesModified = -1 }, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.TryCalculateIsoelectricPoint("C", new() { Disulphides = -1 }, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.TryCalculateIsoelectricPoint("K", new() { Window = (SequenceStatistics.IsoelectricPointWindow)42 }, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.CalculateNetCharge("K", 7.0, new SequenceStatistics.IsoelectricPointOptions { Scale = (SequenceStatistics.PkaScale)42 }));
    }

    #endregion
}
