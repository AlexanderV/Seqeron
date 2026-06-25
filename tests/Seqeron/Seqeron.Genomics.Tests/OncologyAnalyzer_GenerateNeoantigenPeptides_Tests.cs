// ONCO-NEO-001 — Neoantigen Candidate Peptide Window Generation
// Evidence: docs/Evidence/ONCO-NEO-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-NEO-001.md
// Source: Hundal J et al. (2020). pVACtools. Cancer Immunol Res 8(3):409-420.
//         https://doi.org/10.1158/2326-6066.CIR-19-0401
//         Li Y et al. (2020). ProGeo-neo. BMC Med Genomics 13:52.
//         https://doi.org/10.1186/s12920-020-0683-4
//         Wells DK et al. (2020). TESLA. Cell 183(3):818-834.
//         https://doi.org/10.1016/j.cell.2020.09.015
//
// Expected windows below are derived BY HAND from the windowing definition (every length-k window
// of the mutant protein that spans the mutated residue; k in 8..14 for class I (NetMHCpan-4.1 window);
// wild-type agretope at
// the same coordinates) — NOT copied from the implementation output.
// Reference protein (1-based): M1 K2 T3 A4 Y5 I6 A7 K8 Q9 R10 S11 T12 V13 W14 L15 N16 D17 E18 F19 G20 H21.

using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_GenerateNeoantigenPeptides_Tests
{
    private const string WildType = "MKTAYIAKQRSTVWLNDEFGH"; // length 21

    #region GenerateNeoantigenPeptides

    // M1 — Y5C, single length k=9. p=5 is < k-1 from the N-terminus, so starts s=1..5 → 5 windows. Li (2020).
    [Test]
    public void GenerateNeoantigenPeptides_SingleLengthNearNTerminus_FiveWindows()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 9, maxLength: 9);

        int[] expectedStarts = { 1, 2, 3, 4, 5 };
        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(5),
                "p=5 with k=9: starts s in [max(1,5-8)=1, min(5,21-9+1=13)=5] → 5 windows (Li 2020)");
            Assert.That(peptides.Select(p => p.StartPosition), Is.EqualTo(expectedStarts),
                "windows start at positions 1..5, each spanning position 5");
            Assert.That(peptides.All(p => p.Length == 9), Is.True, "all windows are 9-mers");
        });
    }

    // M2 — Y5C, default range k=8..14 (NetMHCpan-4.1 class I window). Each of the 7 lengths yields 5 windows
    // (p=5 near N-terminus; for every k≤14 the right clamp min(5,21-k+1) is still 5) → 35 total. Reynisson (2020).
    [Test]
    public void GenerateNeoantigenPeptides_DefaultClassIRange_ThirtyFivePeptides()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5);

        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(35),
                "k=8..14, 5 windows each (p=5) → 35 candidate peptides (Reynisson 2020 NetMHCpan-4.1 8-14mer; Li 2020 windowing)");
            for (int k = 8; k <= 14; k++)
            {
                int len = k;
                Assert.That(peptides.Count(p => p.Length == len), Is.EqualTo(5), $"5 {len}-mers");
            }
            // Ordered by length ascending then start ascending (INV-06).
            Assert.That(peptides.First().Length, Is.EqualTo(8), "ordered length-ascending: first is an 8-mer");
            Assert.That(peptides.Last().Length, Is.EqualTo(14), "ordered length-ascending: last is a 14-mer");
        });
    }

    // M3 — Y5C, first 8-mer window: mutant MKTACIAK vs WT MKTAYIAK, offset 4 (C replaces Y). Wells (2020) agretope.
    [Test]
    public void GenerateNeoantigenPeptides_FirstEightMer_ExactMutantAndWildTypePair()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 8, maxLength: 8);

        OncologyAnalyzer.NeoantigenPeptide first = peptides[0];
        Assert.Multiple(() =>
        {
            Assert.That(first.StartPosition, Is.EqualTo(1), "first 8-mer starts at position 1");
            Assert.That(first.MutantPeptide, Is.EqualTo("MKTACIAK"),
                "mutant 8-mer at positions 1-8 with Y5→C (derived from WT MKTAYIAK)");
            Assert.That(first.WildTypePeptide, Is.EqualTo("MKTAYIAK"),
                "wild-type agretope at the same coordinates (Wells 2020)");
            Assert.That(first.MutationOffset, Is.EqualTo(4), "mutated residue (position 5) is at 0-based offset 4");
            Assert.That(first.MutantPeptide[first.MutationOffset], Is.EqualTo('C'), "mutant residue at the offset");
            Assert.That(first.WildTypePeptide[first.MutationOffset], Is.EqualTo('Y'), "original residue at the offset");
        });
    }

    // M3b — Y5C, last 8-mer window: start 5 → mutant CIAKQRST vs WT YIAKQRST, offset 0. Li (2020).
    [Test]
    public void GenerateNeoantigenPeptides_LastEightMer_MutationAtWindowStart()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 8, maxLength: 8);

        OncologyAnalyzer.NeoantigenPeptide last = peptides[^1];
        Assert.Multiple(() =>
        {
            Assert.That(last.StartPosition, Is.EqualTo(5), "last 8-mer starts at position 5 (mutation at window start)");
            Assert.That(last.MutantPeptide, Is.EqualTo("CIAKQRST"), "positions 5-12 with Y5→C");
            Assert.That(last.WildTypePeptide, Is.EqualTo("YIAKQRST"), "wild-type agretope");
            Assert.That(last.MutationOffset, Is.EqualTo(0), "mutation is the first residue of this window");
        });
    }

    // M4 — every default peptide spans the mutation and differs from WT only at the offset. Li/Wells (2020).
    [Test]
    public void GenerateNeoantigenPeptides_AllPeptides_SpanMutationAndDifferOnlyAtOffset()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5);

        Assert.Multiple(() =>
        {
            foreach (OncologyAnalyzer.NeoantigenPeptide p in peptides)
            {
                // INV-02: window spans the mutation at protein position 5.
                Assert.That(p.StartPosition + p.MutationOffset, Is.EqualTo(5),
                    $"window starting at {p.StartPosition} must place the mutation at protein position 5");
                Assert.That(p.MutationOffset, Is.InRange(0, p.Length - 1),
                    "mutation offset lies within the peptide");
                // INV-03/04: mutant and WT equal length, differ only at the offset.
                Assert.That(p.MutantPeptide.Length, Is.EqualTo(p.WildTypePeptide.Length),
                    "mutant and wild-type peptides have equal length");
                int diffCount = 0;
                for (int i = 0; i < p.Length; i++)
                {
                    if (p.MutantPeptide[i] != p.WildTypePeptide[i])
                    {
                        diffCount++;
                    }
                }
                Assert.That(diffCount, Is.EqualTo(1),
                    "mutant and wild-type differ at exactly one residue (the substitution)");
                Assert.That(p.MutantPeptide[p.MutationOffset], Is.EqualTo('C'), "mutant residue C at the offset");
                Assert.That(p.WildTypePeptide[p.MutationOffset], Is.EqualTo('Y'), "wild-type residue Y at the offset");
            }
        });
    }

    // M5 — fully interior mutation V13A at p=13 with k=9: starts s in [5,13] → exactly 9 windows. INV-05.
    [Test]
    public void GenerateNeoantigenPeptides_InteriorMutation_ExactlyKWindows()
    {
        // p=13 is >= k-1 (=8) from both ends (13>=9 and 21-13=8>=8), so a full k=9 windows exist.
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'A', 13, minLength: 9, maxLength: 9);

        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(9),
                "interior mutation (>= k-1 from both ends) yields exactly k=9 windows (INV-05)");
            Assert.That(peptides.First().StartPosition, Is.EqualTo(5),
                "first start = p-k+1 = 13-9+1 = 5");
            Assert.That(peptides.Last().StartPosition, Is.EqualTo(13),
                "last start = p = 13");
            Assert.That(peptides.All(p => p.StartPosition + p.MutationOffset == 13), Is.True,
                "all windows span position 13");
        });
    }

    // M6 — N-terminal mutation M1V, k=9: only one window fits (start 1). ProGeo-neo "if possible".
    [Test]
    public void GenerateNeoantigenPeptides_NTerminalMutation_SingleWindow()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'V', 1, minLength: 9, maxLength: 9);

        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(1),
                "mutation at position 1: only the window starting at 1 spans it (ProGeo-neo 'if possible')");
            Assert.That(peptides[0].StartPosition, Is.EqualTo(1), "single window starts at position 1");
            Assert.That(peptides[0].MutantPeptide, Is.EqualTo("VKTAYIAKQ"), "M1→V over positions 1-9");
            Assert.That(peptides[0].WildTypePeptide, Is.EqualTo("MKTAYIAKQ"), "wild-type agretope");
            Assert.That(peptides[0].MutationOffset, Is.EqualTo(0), "mutation at the first residue");
        });
    }

    // M7 — C-terminal mutation H21R, k=8: only one window fits, ending at position 21 (start 14). ProGeo-neo.
    [Test]
    public void GenerateNeoantigenPeptides_CTerminalMutation_SingleWindow()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'R', 21, minLength: 8, maxLength: 8);

        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(1),
                "mutation at position 21 (last): only the window starting at 14 fits and spans it");
            Assert.That(peptides[0].StartPosition, Is.EqualTo(14), "single 8-mer window starts at position 14");
            Assert.That(peptides[0].MutantPeptide, Is.EqualTo("WLNDEFGR"), "positions 14-21 with H21→R");
            Assert.That(peptides[0].WildTypePeptide, Is.EqualTo("WLNDEFGH"), "wild-type agretope");
            Assert.That(peptides[0].MutationOffset, Is.EqualTo(7), "mutation at the last residue of the 8-mer");
        });
    }

    #endregion

    #region Validation and edge cases

    // S1 — null protein → ArgumentNullException.
    [Test]
    public void GenerateNeoantigenPeptides_NullProtein_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(null!, 'C', 5),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "a null protein sequence is invalid input");
    }

    // S2 — empty protein → ArgumentException.
    [Test]
    public void GenerateNeoantigenPeptides_EmptyProtein_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(string.Empty, 'C', 1),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(), "an empty protein sequence cannot yield peptides");
    }

    // S3 — mutant residue equals the wild-type residue → not a substitution → ArgumentException.
    [Test]
    public void GenerateNeoantigenPeptides_NonSubstitution_Throws()
    {
        // Wild-type residue at position 5 is 'Y'; passing 'Y' is not a missense mutation.
        Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'Y', 5),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(), "a missense mutation requires a different amino acid (Hundal 2020)");
    }

    // S4 — mutation position out of range → ArgumentOutOfRangeException.
    [Test]
    public void GenerateNeoantigenPeptides_PositionOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 0),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "position 0 is below the 1-based range");
            Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 22),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "position 22 exceeds the protein length 21");
        });
    }

    // S5 — invalid length range → ArgumentException.
    [Test]
    public void GenerateNeoantigenPeptides_InvalidLengthRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 0, maxLength: 8),
                NUnit.Framework.Throws.TypeOf<ArgumentException>(), "minimum peptide length below 1 is invalid");
            Assert.That(() => OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 11, maxLength: 8),
                NUnit.Framework.Throws.TypeOf<ArgumentException>(), "maxLength < minLength is invalid");
        });
    }

    // C1 — protein shorter than every requested length → empty result (no window fits).
    [Test]
    public void GenerateNeoantigenPeptides_ProteinShorterThanAllLengths_Empty()
    {
        // Protein length 5; default request k=8..14 — none fit.
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides("MKTAY", 'C', 5);

        Assert.That(peptides, Is.Empty, "no class I window (8-14mer) fits a 5-residue protein");
    }

    // C2 — single length subset (min==max==10) returns only 10-mers.
    [Test]
    public void GenerateNeoantigenPeptides_SingleLengthSubset_OnlyThatLength()
    {
        IReadOnlyList<OncologyAnalyzer.NeoantigenPeptide> peptides =
            OncologyAnalyzer.GenerateNeoantigenPeptides(WildType, 'C', 5, minLength: 10, maxLength: 10);

        Assert.Multiple(() =>
        {
            Assert.That(peptides.Count, Is.EqualTo(5), "p=5 with k=10: starts 1..5 → 5 windows");
            Assert.That(peptides.All(p => p.Length == 10), Is.True, "only 10-mers requested and returned");
        });
    }

    #endregion
}
