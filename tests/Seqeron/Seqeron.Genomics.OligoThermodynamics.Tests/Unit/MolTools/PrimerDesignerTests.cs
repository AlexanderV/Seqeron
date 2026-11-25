// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesignerTests
{



    [Test]
    public void DefaultParameters_HasReasonableValues()
    {
        var param = ThermodynamicPrimerDesigner.DefaultParameters;

        Assert.That(param.MinLength, Is.EqualTo(18));
        Assert.That(param.MaxLength, Is.EqualTo(25));
        Assert.That(param.OptimalLength, Is.EqualTo(20));
        Assert.That(param.MinGcContent, Is.EqualTo(40));
        Assert.That(param.MaxGcContent, Is.EqualTo(60));
        Assert.That(param.MinTm, Is.EqualTo(57));
        Assert.That(param.MaxTm, Is.EqualTo(63));
        Assert.That(param.MaxHomopolymer, Is.EqualTo(4));
    }



    /// <summary>
    /// Kills survived ||→&amp;&amp; mutation on EvaluatePrimer line 145:
    /// <c>if (tm &lt; param.MinTm || tm &gt; param.MaxTm)</c>
    /// Tests case where ONLY tm &lt; MinTm is true (Tm is NOT above MaxTm).
    /// With ||: issue flagged. With &amp;&amp;: no issue (both must be true).
    /// </summary>
    [Test]
    public void EvaluatePrimer_TmOnlyBelowMin_FlagsTmIssue()
    {
        // 20bp, 50% GC → Primer3-default Tm = 56.817°C (primer3-py 2.3.1 calc_tm)
        // Primer3 MinTm=57: tm(56.8) < 57 is TRUE, tm(56.8) > 63 is FALSE
        string primer = "GAACTCGTGAACTCGTGAAC"; // 20bp, 10 GC
        var param = new PrimerParameters(
            MinLength: 18, MaxLength: 25, OptimalLength: 20,
            MinGcContent: 0, MaxGcContent: 100,
            MinTm: 57, MaxTm: 63,
            OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
            Avoid3PrimeGC: false, Check3PrimeStability: false);

        var candidate = ThermodynamicPrimerDesigner.EvaluatePrimer(primer, 0, true, param);

        Assert.That(candidate.Issues, Has.Some.Contains("Tm"),
            "Tm below MinTm should flag issue even when Tm is not above MaxTm");
    }

    /// <summary>
    /// Kills survived ||→&amp;&amp; mutation (reverse clause): tests case where
    /// ONLY tm &gt; MaxTm is true (Tm is NOT below MinTm).
    /// </summary>
    [Test]
    public void EvaluatePrimer_TmOnlyAboveMax_FlagsTmIssue()
    {
        // 20bp, 90% GC → Primer3-default Tm = 79.82°C (primer3-py 2.3.1 calc_tm)
        // Primer3 MaxTm=63: tm(79.8) > 63 is TRUE, tm(79.8) < 57 is FALSE
        string primer = "GCGCGCGCGCGCGCGCGCAT"; // 20bp, 18 GC
        var param = new PrimerParameters(
            MinLength: 18, MaxLength: 25, OptimalLength: 20,
            MinGcContent: 0, MaxGcContent: 100,
            MinTm: 57, MaxTm: 63,
            OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
            Avoid3PrimeGC: false, Check3PrimeStability: false);

        var candidate = ThermodynamicPrimerDesigner.EvaluatePrimer(primer, 0, true, param);

        Assert.That(candidate.Issues, Has.Some.Contains("Tm"),
            "Tm above MaxTm should flag issue even when Tm is not below MinTm");
    }

    /// <summary>
    /// Kills survived &lt;→&lt;= boundary mutation on MinTm comparison.
    /// When Tm == MinTm, the original &lt; returns false (no issue).
    /// With &lt;= mutation, it returns true (issue added).
    /// </summary>
    [Test]
    public void EvaluatePrimer_TmExactlyAtMinTm_NoTmIssue()
    {
        string primer = "ATGCGATCGATCGATCGATC"; // 20bp
        double exactTm = ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(primer);

        var param = new PrimerParameters(
            MinLength: 18, MaxLength: 25, OptimalLength: 20,
            MinGcContent: 0, MaxGcContent: 100,
            MinTm: exactTm, MaxTm: 100,
            OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
            Avoid3PrimeGC: false, Check3PrimeStability: false);

        var candidate = ThermodynamicPrimerDesigner.EvaluatePrimer(primer, 0, true, param);

        Assert.That(candidate.Issues, Has.None.Contains("Tm"),
            "Tm exactly at MinTm boundary should not flag issue");
    }

    /// <summary>
    /// Kills survived &gt;→&gt;= boundary mutation on MaxTm comparison.
    /// When Tm == MaxTm, the original &gt; returns false (no issue).
    /// With &gt;= mutation, it returns true (issue added).
    /// </summary>
    [Test]
    public void EvaluatePrimer_TmExactlyAtMaxTm_NoTmIssue()
    {
        string primer = "ATGCGATCGATCGATCGATC"; // 20bp
        double exactTm = ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(primer);

        var param = new PrimerParameters(
            MinLength: 18, MaxLength: 25, OptimalLength: 20,
            MinGcContent: 0, MaxGcContent: 100,
            MinTm: 0, MaxTm: exactTm,
            OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
            Avoid3PrimeGC: false, Check3PrimeStability: false);

        var candidate = ThermodynamicPrimerDesigner.EvaluatePrimer(primer, 0, true, param);

        Assert.That(candidate.Issues, Has.None.Contains("Tm"),
            "Tm exactly at MaxTm boundary should not flag issue");
    }
}
