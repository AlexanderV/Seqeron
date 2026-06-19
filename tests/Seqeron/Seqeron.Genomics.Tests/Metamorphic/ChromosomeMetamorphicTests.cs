using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Chromosome area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-TELO-001 — telomere analysis (Chromosome).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 48.
///
/// API under test (ChromosomeAnalyzer.AnalyzeTelomeres):
///   Scans the 3' end for tandem TTAGGG repeats (≥70 % per-hexamer identity) counting inward,
///   and the 5' start for the reverse-complement CCCTAA; reports the terminal repeat-tract
///   LENGTHS and purities at each end.
///
/// Relations (derived from that definition, NOT from output):
///   • MON (more repeats ⇒ longer telomere): a 3' tract of k TTAGGG units (with a
///          non-telomeric interior) is measured as length 6k, strictly increasing in k; the
///          5' CCCTAA tract behaves symmetrically.
///   • INV (interior flank doesn't affect the core): inserting non-telomeric sequence BETWEEN
///          the two terminal tracts leaves both telomere lengths and purities unchanged — the
///          ends (the biological "core" of the measurement) are untouched.
///   • SHIFT/strand-duality: this API reports terminal LENGTHS, not interior coordinates, so
///          the checklist's positional-shift relation is realised as the strand duality that
///          underlies telomere calling — reverse-complementing the chromosome maps the 3'
///          TTAGGG tract onto a 5' CCCTAA tract, swapping the 5' and 3' measurements exactly.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ChromosomeMetamorphicTests
{
    #region Helpers

    private const string TeloRepeat = "TTAGGG";              // 3' tract
    private const string TeloRepeatRc = "CCCTAA";            // 5' tract (reverse complement)

    private static string Repeat(string unit, int count) => string.Concat(Enumerable.Repeat(unit, count));

    /// <summary>Non-telomeric filler whose every hexamer falls well below the 70 % match threshold for both tracts.</summary>
    private static string Filler(int length)
    {
        var sb = new System.Text.StringBuilder(length);
        while (sb.Length < length) sb.Append("GC");
        return sb.ToString(0, length);
    }

    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    #endregion

    #region MON — more terminal repeats give a longer telomere

    [Test]
    [Description("MON: a 3' tract of k TTAGGG units is measured as telomere length 6k, strictly increasing in k; the 5' CCCTAA tract behaves symmetrically.")]
    public void AnalyzeTelomeres_MoreRepeats_IncreasesLength()
    {
        int prev3 = int.MinValue, prev5 = int.MinValue;

        foreach (int k in new[] { 1, 2, 5, 10, 20, 50 })
        {
            var three = ChromosomeAnalyzer.AnalyzeTelomeres("chr", Filler(60) + Repeat(TeloRepeat, k));
            three.TelomereLength3Prime.Should().Be(6 * k,
                because: $"the 3' end carries {k} contiguous TTAGGG units, each contributing 6 bp to the measured telomere");
            three.TelomereLength3Prime.Should().BeGreaterThan(prev3,
                because: "adding TTAGGG units to the 3' end strictly lengthens the measured telomere");
            prev3 = three.TelomereLength3Prime;

            var five = ChromosomeAnalyzer.AnalyzeTelomeres("chr", Repeat(TeloRepeatRc, k) + Filler(60));
            five.TelomereLength5Prime.Should().Be(6 * k,
                because: $"the 5' end carries {k} contiguous CCCTAA units (the reverse complement of TTAGGG)");
            five.TelomereLength5Prime.Should().BeGreaterThan(prev5,
                because: "adding CCCTAA units to the 5' end strictly lengthens the measured telomere");
            prev5 = five.TelomereLength5Prime;
        }
    }

    #endregion

    #region INV — interior non-telomeric sequence does not change the telomeres

    [Test]
    [Description("INV: inserting non-telomeric sequence between the two terminal tracts leaves both telomere lengths and purities unchanged.")]
    public void AnalyzeTelomeres_InteriorFlank_PreservesTelomeres()
    {
        const int a = 30, b = 40;   // 5' and 3' tract repeat counts
        string fivePrime = Repeat(TeloRepeatRc, a);
        string threePrime = Repeat(TeloRepeat, b);

        var baseline = ChromosomeAnalyzer.AnalyzeTelomeres("chr", fivePrime + Filler(100) + threePrime);

        foreach (int extra in new[] { 200, 1000, 5000 })
        {
            var grown = ChromosomeAnalyzer.AnalyzeTelomeres("chr", fivePrime + Filler(100 + extra) + threePrime);

            grown.TelomereLength5Prime.Should().Be(baseline.TelomereLength5Prime,
                because: "the 5' terminal tract is untouched by interior insertions, so its measured length is unchanged");
            grown.TelomereLength3Prime.Should().Be(baseline.TelomereLength3Prime,
                because: "the 3' terminal tract is untouched by interior insertions, so its measured length is unchanged");
            grown.RepeatPurity5Prime.Should().BeApproximately(baseline.RepeatPurity5Prime, 1e-12,
                because: "the 5' tract content is unchanged, so its repeat purity is unchanged");
            grown.RepeatPurity3Prime.Should().BeApproximately(baseline.RepeatPurity3Prime, 1e-12,
                because: "the 3' tract content is unchanged, so its repeat purity is unchanged");
        }
    }

    #endregion

    #region SHIFT/strand-duality — reverse complement swaps the 5' and 3' telomeres

    [Test]
    [Description("SHIFT/strand-duality: reverse-complementing the chromosome maps the 3' TTAGGG tract onto a 5' CCCTAA tract, swapping the 5' and 3' telomere measurements.")]
    public void AnalyzeTelomeres_ReverseComplement_SwapsEnds()
    {
        const int a = 25, b = 60;   // asymmetric tracts so the swap is observable
        string seq = Repeat(TeloRepeatRc, a) + Filler(120) + Repeat(TeloRepeat, b);

        var forward = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
        var reversed = ChromosomeAnalyzer.AnalyzeTelomeres("chr", RevComp(seq));

        reversed.TelomereLength5Prime.Should().Be(forward.TelomereLength3Prime,
            because: "reverse-complementing turns the 3' TTAGGG tract into a 5' CCCTAA tract of equal length");
        reversed.TelomereLength3Prime.Should().Be(forward.TelomereLength5Prime,
            because: "reverse-complementing turns the 5' CCCTAA tract into a 3' TTAGGG tract of equal length");
        reversed.RepeatPurity5Prime.Should().BeApproximately(forward.RepeatPurity3Prime, 1e-12,
            because: "the swapped tracts carry the same residues, so their purities swap too");
        reversed.RepeatPurity3Prime.Should().BeApproximately(forward.RepeatPurity5Prime, 1e-12,
            because: "the swapped tracts carry the same residues, so their purities swap too");
    }

    #endregion
}
