using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Epigenetics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-CPG-001 — CpG ratio / sites / island detection (Epigenetics).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 85.
///
/// API under test (EpigeneticsAnalyzer.CalculateCpGObservedExpected / FindCpGSites /
///                 FindCpGIslands):
///   CpG O/E ratio = (#CpG dinucleotides) / (C·G / length). CpG sites are the positions of
///   "CG"; islands (Gardiner-Garden & Frommer) are long windows with high GC% and O/E ratio.
///
/// Relations (derived from the dinucleotide counting, NOT from output):
///   • MON  (more CG dinucleotides ⇒ higher ratio): at fixed C, G and length, arranging the
///          bases into more adjacent CG dinucleotides raises the O/E ratio.
///   • SHIFT (prepend flank shifts positions): a non-CG prepended flank shifts every CpG-site
///          position by the flank length.
///   • INV  (non-CG flank ⇒ same island detection): appending an AT-rich, CG-free flank leaves
///          the detected island (its span and GC%) unchanged.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class EpigeneticsMetamorphicTests
{
    #region EPIGEN-CPG-001 MON — more CG dinucleotides raise the O/E ratio

    [Test]
    [Description("MON: holding C, G and length fixed, rearranging the bases into more adjacent CG dinucleotides strictly raises the CpG observed/expected ratio.")]
    public void CpGRatio_MoreCgDinucleotides_HigherRatio()
    {
        // All four sequences have exactly 3 C, 3 G, length 6 — only the # of CG dinucleotides differs.
        var byCgCount = new[]
        {
            "GGGCCC", // 0 CpG
            "CCCGGG", // 1 CpG
            "CCGCGG", // 2 CpG
            "CGCGCG", // 3 CpG
        };

        double previous = double.MinValue;
        foreach (var seq in byCgCount)
        {
            double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            ratio.Should().BeGreaterThan(previous,
                because: $"'{seq}' packs more CG dinucleotides than the previous arrangement at the same composition");
            previous = ratio;
        }
    }

    #endregion

    #region EPIGEN-CPG-001 SHIFT — prepending a non-CG flank shifts CpG-site positions

    [Test]
    [Description("SHIFT: a prepended flank with no CG dinucleotide (and no CG at the junction) shifts every CpG-site position by exactly the flank length.")]
    public void FindCpGSites_PrependFlank_ShiftsPositions()
    {
        const string core = "ATCGATCGAT"; // CG dinucleotides at indices 2 and 6
        var basePositions = EpigeneticsAnalyzer.FindCpGSites(core).ToList();
        basePositions.Should().NotBeEmpty(because: "the core contains CG dinucleotides");

        foreach (int k in new[] { 2, 4, 6 })
        {
            string flank = string.Concat(Enumerable.Repeat("AT", k / 2)); // no CG, junction 'T'+'A' is not CG
            EpigeneticsAnalyzer.FindCpGSites(flank + core).Should().Equal(basePositions.Select(p => p + k),
                because: $"a CG-free {k}-base flank shifts every CpG-site position by {k}");
        }
    }

    #endregion

    #region EPIGEN-CPG-001 INV — a non-CG flank leaves island detection unchanged

    [Test]
    [Description("INV: appending an AT-rich, CG-free flank still detects the CpG island over the CG-dense core — the island remains called (start 0, covering the whole core, GC%/ratio above threshold). The exact 3' boundary may extend by the sliding window, but detection of the core island is preserved.")]
    public void FindCpGIslands_NonCgFlank_CoreIslandStillDetected()
    {
        // A 300-bp perfect CpG island.
        const int coreLength = 300;
        string island = string.Concat(Enumerable.Repeat("CG", coreLength / 2));
        EpigeneticsAnalyzer.FindCpGIslands(island).Should().ContainSingle(i => i.Start == 0);

        foreach (var flank in new[]
                 {
                     string.Concat(Enumerable.Repeat("AT", 50)),   // AT-rich, CG-free, low GC → not itself an island
                     string.Concat(Enumerable.Repeat("TA", 100)),
                 })
        {
            var detected = EpigeneticsAnalyzer.FindCpGIslands(island + flank).Single(i => i.Start == 0);

            detected.End.Should().BeGreaterThanOrEqualTo(coreLength,
                because: "the CG-dense core is still entirely covered by the detected island");
            detected.GcContent.Should().BeGreaterThanOrEqualTo(0.5, because: "the detected island still passes the GC% criterion");
            detected.CpGRatio.Should().BeGreaterThanOrEqualTo(0.6, because: "the detected island still passes the CpG O/E criterion");
        }
    }

    #endregion
}
