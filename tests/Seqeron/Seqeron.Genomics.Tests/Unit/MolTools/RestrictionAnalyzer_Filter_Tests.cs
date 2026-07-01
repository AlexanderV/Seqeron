// RESTR-FILTER-001 — Enzyme Filtering
// Evidence: docs/Evidence/RESTR-FILTER-001-Evidence.md
// TestSpec: tests/TestSpecs/RESTR-FILTER-001.md
// Source: Wikipedia "Sticky and blunt ends" (2026); Wikipedia "Restriction enzyme" (2026);
//         Wikipedia "List of restriction enzyme cutting sites" (2026); NEB/REBASE KpnI R0142, EcoRI-HF R3101.

using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Canonical test class for RESTR-FILTER-001: Restriction Enzyme Filtering.
/// Tests GetBluntCutters, GetStickyCutters, and GetEnzymesByCutLength (range + single-length)
/// of RestrictionAnalyzer.
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Filter_Tests
{
    #region GetBluntCutters

    // M1 — Blunt set contains center-cutters (SmaI/EcoRV/AluI/HaeIII) and excludes
    //      overhang producers (EcoRI 5', KpnI 3'). Evidence: Wikipedia "Sticky and blunt ends"
    //      (blunt = both strands terminate in a base pair); "Restriction enzyme" (SmaI blunt).
    [Test]
    public void GetBluntCutters_KnownEnzymes_IncludesBluntExcludesSticky()
    {
        // Act
        var names = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(names, Does.Contain("SmaI"), "SmaI (CCC^GGG) is a center cut → blunt");
            Assert.That(names, Does.Contain("EcoRV"), "EcoRV (GAT^ATC) is a center cut → blunt");
            Assert.That(names, Does.Contain("AluI"), "AluI (AG^CT) is a center cut → blunt");
            Assert.That(names, Does.Contain("HaeIII"), "HaeIII (GG^CC) is a center cut → blunt");
            Assert.That(names, Does.Not.Contain("EcoRI"), "EcoRI leaves a 5' overhang → sticky, not blunt");
            Assert.That(names, Does.Not.Contain("KpnI"), "KpnI leaves a 3' overhang → sticky, not blunt");
        });
    }

    // M1 — Every enzyme returned by GetBluntCutters must satisfy the blunt criterion cf == cr.
    [Test]
    public void GetBluntCutters_AllResults_HaveEqualForwardAndReverseCut()
    {
        // Act
        var blunt = RestrictionAnalyzer.GetBluntCutters().ToList();

        // Assert: blunt iff the two strands are cut at the same offset (center cut).
        Assert.That(blunt.All(e => e.CutPositionForward == e.CutPositionReverse), Is.True,
            "A blunt cutter must cut both strands at the same position (Wikipedia, Restriction enzyme)");
    }

    // M1 — EXACT blunt-cutter set. Every name below is the documented center-cut (blunt) enzyme,
    //      each cut site cross-checked against REBASE / NEB / Wikipedia this session:
    //        AluI  AG^CT   [blunt]   DpnI  Gm6A^TC [blunt]  EcoRV GAT^ATC [blunt]
    //        HaeIII GG^CC  [blunt]   HincII GTY^RAC [blunt]  RsaI  GT^AC   [blunt]
    //        ScaI  AGT^ACT [blunt]   SmaI  CCC^GGG [blunt]   StuI  AGG^CCT [blunt]
    //        SwaI  ATTT^AAAT [blunt]
    //      A locked exact-set assertion: a wrong cut-position for ANY one enzyme (mis-classifying
    //      it as sticky, or pulling a sticky enzyme in) fails here — not just the 4-of-10 spot check.
    //      Sources: NEB R0176 DpnI (Gm6A^TC blunt); NEB R0103 HincII (GTY^RAC blunt);
    //      NEB R0167 RsaI (GT^AC blunt); ScaI/StuI/SwaI center cut (REBASE/Thermo/research);
    //      Wikipedia "Restriction enzyme" (SmaI/EcoRV blunt examples).
    [Test]
    public void GetBluntCutters_ReturnsExactlyTheDocumentedBluntSet()
    {
        // Arrange: the complete, externally-sourced set of blunt (center-cut) enzymes in the library.
        var expectedBlunt = new[]
        {
            "AluI", "DpnI", "EcoRV", "HaeIII", "HincII", "RsaI", "ScaI", "SmaI", "StuI", "SwaI"
        };

        // Act
        var blunt = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToHashSet();

        // Assert: exact set equality — no missing blunt cutter, no sticky cutter leaking in.
        Assert.That(blunt, Is.EquivalentTo(expectedBlunt),
            "GetBluntCutters must return exactly the documented center-cut (blunt) enzymes");
    }

    #endregion

    #region GetStickyCutters

    // M2 — Sticky set contains overhang producers (EcoRI 5', KpnI/PstI 3') and excludes
    //      blunt cutters (SmaI/EcoRV). Evidence: Wikipedia "Restriction enzyme"; NEB/REBASE KpnI.
    [Test]
    public void GetStickyCutters_KnownEnzymes_IncludesStickyExcludesBlunt()
    {
        // Act
        var names = RestrictionAnalyzer.GetStickyCutters().Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(names, Does.Contain("EcoRI"), "EcoRI (G^AATTC) leaves a 5' overhang → sticky");
            Assert.That(names, Does.Contain("KpnI"), "KpnI (GGTAC^C) leaves a 3' overhang → sticky");
            Assert.That(names, Does.Contain("PstI"), "PstI (CTGCA^G) leaves a 3' overhang → sticky");
            Assert.That(names, Does.Not.Contain("SmaI"), "SmaI is blunt, not sticky");
            Assert.That(names, Does.Not.Contain("EcoRV"), "EcoRV is blunt, not sticky");
        });
    }

    #endregion

    #region Blunt/Sticky Partition

    // M3 — Blunt and sticky sets are disjoint and partition the full library (INV-1).
    //      Evidence: Wikipedia "Sticky and blunt ends" — every end is blunt or an overhang.
    [Test]
    public void BluntAndStickyCutters_PartitionTheLibrary()
    {
        // Arrange
        var library = RestrictionAnalyzer.Enzymes.Values.Select(e => e.Name).ToHashSet();
        var blunt = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToHashSet();
        var sticky = RestrictionAnalyzer.GetStickyCutters().Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blunt.Overlaps(sticky), Is.False, "Blunt and sticky sets must be disjoint");
            Assert.That(blunt.Count + sticky.Count, Is.EqualTo(library.Count),
                "Blunt + sticky counts must sum to the full library size (total partition)");
            var union = new System.Collections.Generic.HashSet<string>(blunt);
            union.UnionWith(sticky);
            Assert.That(union, Is.EquivalentTo(library), "Blunt ∪ Sticky must equal the full library");
        });
    }

    #endregion

    #region GetEnzymesByCutLength (range)

    // M4 — Range [6,6] returns exactly the 6-bp recognition enzymes.
    //      Evidence: Wikipedia "List of restriction enzyme cutting sites" (EcoRI/BamHI/PstI = 6 bp).
    [Test]
    public void GetEnzymesByCutLength_Range6To6_ReturnsOnly6Cutters()
    {
        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(6, 6).ToList();
        var names = result.Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.All(e => e.RecognitionLength == 6), Is.True, "All returned sites must be 6 bp");
            Assert.That(names, Does.Contain("EcoRI"), "EcoRI (GAATTC) is a 6-cutter");
            Assert.That(names, Does.Contain("BamHI"), "BamHI (GGATCC) is a 6-cutter");
            Assert.That(names, Does.Contain("PstI"), "PstI (CTGCAG) is a 6-cutter");
            Assert.That(names, Does.Not.Contain("AluI"), "AluI (AGCT) is a 4-cutter, excluded");
            Assert.That(names, Does.Not.Contain("NotI"), "NotI (GCGGCCGC) is an 8-cutter, excluded");
        });
    }

    // M5 — Range [4,8] returns the full library EXCEPT the interrupted-palindrome SfiI.
    //      The 4–8 nt range applies to undivided Type II sites (Wikipedia, "Restriction enzyme");
    //      SfiI is a divided palindrome GGCCNNNN^NGGCC (13-nt recognition string, PMC PMC548270),
    //      so it lies outside [4,8] and is correctly excluded.
    [Test]
    public void GetEnzymesByCutLength_Range4To8_ReturnsAllUndividedSites_ExcludesSfiI()
    {
        // Arrange: every other library enzyme has a 4–8 nt undivided recognition string.
        int total = RestrictionAnalyzer.Enzymes.Count;

        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(4, 8).ToList();
        var names = result.Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(total - 1),
                "[4,8] returns every enzyme except the single 13-nt interrupted palindrome");
            Assert.That(names, Does.Not.Contain("SfiI"),
                "SfiI (GGCCNNNNNGGCC, length 13) is a divided palindrome outside the 4–8 nt range");
            Assert.That(result.All(e => e.RecognitionLength is >= 4 and <= 8), Is.True,
                "All returned recognition strings must be 4–8 nt");
        });
    }

    // M6 — Range above the 8-nt maximum returns nothing.
    //      Evidence: Wikipedia "Restriction enzyme" — recognition sites are at most 8 nt.
    [Test]
    public void GetEnzymesByCutLength_Range9To10_ReturnsEmpty()
    {
        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(9, 10).ToList();

        // Assert
        Assert.That(result, Is.Empty, "No recognition site exceeds 8 nt, so [9,10] is empty");
    }

    // M7 — Range [4,4] returns only the 4-bp recognition enzymes.
    //      Evidence: Wikipedia "List of restriction enzyme cutting sites" (AluI/HaeIII/TaqI = 4 bp).
    [Test]
    public void GetEnzymesByCutLength_Range4To4_ReturnsOnly4Cutters()
    {
        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(4, 4).ToList();
        var names = result.Select(e => e.Name).ToHashSet();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.All(e => e.RecognitionLength == 4), Is.True, "All returned sites must be 4 bp");
            Assert.That(names, Does.Contain("AluI"), "AluI (AGCT) is a 4-cutter");
            Assert.That(names, Does.Contain("HaeIII"), "HaeIII (GGCC) is a 4-cutter");
            Assert.That(names, Does.Contain("TaqI"), "TaqI (TCGA) is a 4-cutter");
            Assert.That(names, Does.Not.Contain("EcoRI"), "EcoRI is a 6-cutter, excluded");
        });
    }

    // S1 — Inverted range (min > max) returns empty (INV-5).
    [Test]
    public void GetEnzymesByCutLength_MinGreaterThanMax_ReturnsEmpty()
    {
        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(8, 4).ToList();

        // Assert
        Assert.That(result, Is.Empty, "An inverted range describes an empty interval → no enzymes");
    }

    // S2 — Range [L,L] equals the single-length overload for L in {4,6,8} (INV-4).
    [Test]
    public void GetEnzymesByCutLength_EqualBounds_MatchesSingleLengthOverload()
    {
        Assert.Multiple(() =>
        {
            foreach (int length in new[] { 4, 6, 8 })
            {
                var range = RestrictionAnalyzer.GetEnzymesByCutLength(length, length).Select(e => e.Name).ToHashSet();
                var single = RestrictionAnalyzer.GetEnzymesByCutLength(length).Select(e => e.Name).ToHashSet();
                Assert.That(range, Is.EquivalentTo(single),
                    $"Range [{length},{length}] must equal the single-length overload for {length}");
            }
        });
    }

    // C1 — Non-positive bounds return empty (no site has length ≤ 0).
    [Test]
    public void GetEnzymesByCutLength_NonPositiveBounds_ReturnsEmpty()
    {
        // Act
        var result = RestrictionAnalyzer.GetEnzymesByCutLength(-1, 0).ToList();

        // Assert
        Assert.That(result, Is.Empty, "No recognition site has length ≤ 0");
    }

    #endregion
}
