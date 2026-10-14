// RESTR-FILTER-001 — Enzyme Filtering
// Evidence: docs/Evidence/RESTR-FILTER-001-Evidence.md
// TestSpec: tests/TestSpecs/RESTR-FILTER-001.md
// Source: Wikipedia "Sticky and blunt ends" (2026); Wikipedia "Restriction enzyme" (2026);
//         Wikipedia "List of restriction enzyme cutting sites" (2026); NEB/REBASE KpnI R0142, EcoRI-HF R3101.
// Oracle (RESTR-FILTER-001, this session): Biopython 1.88 Bio.Restriction — per enzyme `.size`,
//         `.ovhg`, `.ovhgseq`, `.is_blunt()`, `.is_5overhang()`, `.is_3overhang()`,
//         `.equischizomers(batch)` / `.neoschizomers(batch)` / `.isoschizomers(batch)`, and the
//         `%` compatible-ends operator. Audited: 39/39 enzymes on site string, size, cut
//         positions, blunt/sticky, 5'/3' polarity and overhang sequence; 1521/1521 ordered pairs
//         on compatibility; 75 466/75 466 pairs of the 389 REBASE enzymes with a fully-defined
//         overhang. Every value asserted below was read off that oracle.

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

    #region GetOverhangSequence (REBASE / Biopython ovhgseq)

    // F28 — The enzyme-level overhang sequence is REBASE/Biopython `ovhgseq`: the top-strand bases
    //       between the two strand cuts, 5'→3'. Values read from Biopython 1.88 this session.
    [Test]
    [TestCase("EcoRI", "AATT")]   // G^AATT_C   ovhg -4
    [TestCase("BamHI", "GATC")]   // G^GATC_C   ovhg -4
    [TestCase("BglII", "GATC")]   // A^GATC_T   ovhg -4
    [TestCase("Sau3AI", "GATC")]  // N^GATC_N   ovhg -4
    [TestCase("MboI", "GATC")]    // N^GATC_N   ovhg -4
    [TestCase("SalI", "TCGA")]    // G^TCGA_C   ovhg -4
    [TestCase("XhoI", "TCGA")]    // C^TCGA_G   ovhg -4
    [TestCase("NheI", "CTAG")]
    [TestCase("XbaI", "CTAG")]
    [TestCase("SpeI", "CTAG")]
    [TestCase("AvrII", "CTAG")]
    [TestCase("HindIII", "AGCT")] // A^AGCT_T   ovhg -4
    [TestCase("NotI", "GGCC")]    // GC^GGCC_GC ovhg -4
    [TestCase("AscI", "CGCG")]
    [TestCase("NdeI", "TA")]      // CA^TA_TG   ovhg -2
    [TestCase("ClaI", "CG")]      // AT^CG_AT   ovhg -2
    [TestCase("TaqI", "CG")]      // T^CG_A     ovhg -2
    [TestCase("MspI", "CG")]
    [TestCase("HpaII", "CG")]
    [TestCase("AgeI", "CCGG")]
    [TestCase("NcoI", "CATG")]
    [TestCase("PstI", "TGCA")]    // C_TGCA^G   ovhg +4
    [TestCase("SacI", "AGCT")]    // G_AGCT^C   ovhg +4
    [TestCase("KpnI", "GTAC")]    // G_GTAC^C   ovhg +4
    [TestCase("SphI", "CATG")]    // G_CATG^C   ovhg +4
    [TestCase("ApaI", "GGCC")]    // G_GGCC^C   ovhg +4
    [TestCase("FseI", "CCGG")]    // GG_CCGG^CC ovhg +4
    [TestCase("PacI", "AT")]      // TTA_AT^TAA ovhg +2
    [TestCase("SfiI", "NNN")]     // GGCCN_NNN^NGGCC ovhg +3 — the reference keeps the N codes
    [TestCase("EcoRV", "")]       // GAT^_ATC   ovhg 0 (blunt)
    [TestCase("SmaI", "")]
    [TestCase("HincII", "")]      // GTY^_RAC   ovhg 0 (blunt) although the site is degenerate
    [TestCase("SwaI", "")]        // ATTT^_AAAT ovhg 0
    [TestCase("DpnI", "")]        // GA^_TC     ovhg 0
    public void GetOverhangSequence_MatchesBiopythonOvhgseq(string enzyme, string expected)
    {
        Assert.That(RestrictionAnalyzer.GetOverhangSequence(enzyme), Is.EqualTo(expected),
            $"{enzyme} overhang must equal Biopython Bio.Restriction.{enzyme}.ovhgseq");
    }

    // F28 — Every built-in entry: |overhang| == |ovhg| and the overhang is empty exactly for the
    //       blunt cutters. Locks the whole table, not a sample.
    [Test]
    public void GetOverhangSequence_WholeTable_LengthEqualsCutStaggerAndEmptyOnlyWhenBlunt()
    {
        Assert.Multiple(() =>
        {
            foreach (var enzyme in RestrictionAnalyzer.Enzymes.Values)
            {
                string overhang = RestrictionAnalyzer.GetOverhangSequence(enzyme.Name);
                int stagger = Math.Abs(enzyme.CutPositionForward - enzyme.CutPositionReverse);

                Assert.That(overhang, Has.Length.EqualTo(stagger),
                    $"{enzyme.Name}: overhang length must equal |ovhg|");
                Assert.That(overhang.Length == 0, Is.EqualTo(enzyme.IsBluntEnd),
                    $"{enzyme.Name}: an empty overhang must mean blunt and vice versa");
            }
        });
    }

    // F28 — An unknown name throws, like every other name-taking member of the class.
    [Test]
    public void GetOverhangSequence_UnknownEnzyme_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetOverhangSequence("Nope"));
        Assert.That(ex!.Message, Does.Contain("Unknown enzyme: Nope"));
    }

    // F28 — A cutter that cuts OUTSIDE its own recognition sequence (Type IIS-style distal cut)
    //       has no enzyme-level overhang: it depends on the flanking molecule, so the record
    //       overload reports null rather than slicing out of range.
    [Test]
    public void GetOverhangSequence_DistalCutter_IsNotDerivableFromThePattern()
    {
        // BsaI is GGTCTC(1/5): the cut is 1 nt past the site on the top strand and 5 nt past it on
        // the bottom strand, i.e. offsets 7 and 11 into a 6-nt recognition string.
        var distal = new RestrictionEnzyme("BsaI-like", "GGTCTC", 7, 11, "test");

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(distal), Is.Null,
                "the overhang window lies outside the recognition sequence");
            Assert.That(RestrictionAnalyzer.AreCompatible(distal, distal), Is.False,
                "an undetermined overhang is compatible with nothing (Biopython NotDefined)");
        });
    }

    #endregion

    #region AreCompatible — ligation rule

    // F29 — Confirmed compatible families. Overhang sequences and the verdicts are Biopython's
    //       (`BamHI % BglII` etc. all True this session).
    [Test]
    [TestCase("BamHI", "BglII", "GATC")]
    [TestCase("BamHI", "Sau3AI", "GATC")]
    [TestCase("BglII", "MboI", "GATC")]
    [TestCase("MboI", "Sau3AI", "GATC")]
    [TestCase("SalI", "XhoI", "TCGA")]
    [TestCase("NheI", "XbaI", "CTAG")]
    [TestCase("XbaI", "SpeI", "CTAG")]
    [TestCase("SpeI", "AvrII", "CTAG")]
    [TestCase("ClaI", "TaqI", "CG")]
    [TestCase("MspI", "HpaII", "CG")]
    [TestCase("TaqI", "MspI", "CG")]
    public void AreCompatible_SameFivePrimeOverhang_IsCompatible(string e1, string e2, string end)
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.AreCompatible(e1, e2), Is.True);
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(e1), Is.EqualTo(end));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(e2), Is.EqualTo(end));
        });
    }

    // F29 — EcoRI's AATT is unique in this table, so EcoRI is compatible with itself and nothing
    //       else. (REBASE MfeI, CAATTG, is EcoRI's AATT partner but is not a built-in here.)
    [Test]
    public void AreCompatible_UniqueOverhang_OnlyCompatibleWithItself()
    {
        var partners = RestrictionAnalyzer.FindCompatibleEnzymes("EcoRI").Select(t => t.Enzyme).ToList();

        Assert.That(partners, Is.EquivalentTo(new[] { "EcoRI" }),
            "no other built-in leaves a 5' AATT overhang");
    }

    // F29 — THE POLARITY RULE: identical overhang letters with opposite polarity do NOT ligate —
    //       both single-stranded extensions would lie on the same strand of the joint. All four
    //       such pairs that exist in the table are asserted, and Biopython agrees on all four.
    [Test]
    [TestCase("HindIII", "SacI", "AGCT")]  // 5' AGCT vs 3' AGCT
    [TestCase("NcoI", "SphI", "CATG")]     // 5' CATG vs 3' CATG
    [TestCase("NotI", "ApaI", "GGCC")]     // 5' GGCC vs 3' GGCC
    [TestCase("AgeI", "FseI", "CCGG")]     // 5' CCGG vs 3' CCGG
    public void AreCompatible_SameLettersOppositePolarity_IsNotCompatible(string fivePrime, string threePrime, string letters)
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(fivePrime), Is.EqualTo(letters));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(threePrime), Is.EqualTo(letters));
            Assert.That(RestrictionAnalyzer.GetEnzyme(fivePrime)!.OverhangType,
                Is.EqualTo(OverhangType.FivePrime));
            Assert.That(RestrictionAnalyzer.GetEnzyme(threePrime)!.OverhangType,
                Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(RestrictionAnalyzer.AreCompatible(fivePrime, threePrime), Is.False,
                $"a 5' {letters} extension cannot anneal to a 3' {letters} extension");
            Assert.That(RestrictionAnalyzer.AreCompatible(threePrime, fivePrime), Is.False,
                "and the rule is symmetric");
        });
    }

    // F29 — Blunt ligates to blunt; blunt never ligates to a sticky end.
    [Test]
    public void AreCompatible_BluntAndSticky_FollowTheEndChemistry()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI"), Is.True, "blunt + blunt");
            Assert.That(RestrictionAnalyzer.AreCompatible("HincII", "SwaI"), Is.True, "blunt + blunt");
            Assert.That(RestrictionAnalyzer.AreCompatible("EcoRV", "EcoRI"), Is.False, "blunt + 5'");
            Assert.That(RestrictionAnalyzer.AreCompatible("SmaI", "PstI"), Is.False, "blunt + 3'");
        });
    }

    // F29 — Different overhang length is never compatible, same polarity or not. SfiI's 3-nt
    //       degenerate 3' overhang is NOT reported compatible with the 4-nt 3' overhangs, which is
    //       Biopython's verdict too (`SfiI % PstI` and `SfiI % ApaI` are both False — length).
    [Test]
    [TestCase("SfiI", "PstI")]
    [TestCase("SfiI", "ApaI")]
    [TestCase("SfiI", "KpnI")]
    [TestCase("NdeI", "EcoRI")]
    [TestCase("ClaI", "BamHI")]
    [TestCase("PacI", "PstI")]
    public void AreCompatible_DifferentOverhangLength_IsNotCompatible(string e1, string e2)
    {
        Assert.That(RestrictionAnalyzer.AreCompatible(e1, e2), Is.False);
    }

    // F29 — A degenerate overhang is matched code-wise, not as a literal string, so SfiI's "NNN"
    //       is compatible with itself (Biopython: `SfiI % SfiI` is True) and with any other 3-nt
    //       3' overhang. The second case is not reachable from the built-in table — no other
    //       built-in has a 3-nt 3' overhang — so it is exercised through the record overload with
    //       REBASE values: Psp03I/FmuI leave a 3' "GNC" and PssI a 3' "GNC"; BthCI leaves "CNG".
    [Test]
    public void AreCompatible_DegenerateOverhang_MatchesCodeWiseNotAsAString()
    {
        // SfiI GGCCNNNNNGGCC 8/5 → 3' "NNN"; a hypothetical 3-nt 3' cutter leaving "GTC".
        var gtc = new RestrictionEnzyme("GTC-3p", "AAGTCAA", 5, 2, "test");

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence("SfiI"), Is.EqualTo("NNN"));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(gtc), Is.EqualTo("GTC"));
            Assert.That(RestrictionAnalyzer.AreCompatible("SfiI", "SfiI"), Is.True,
                "Biopython SfiI % SfiI is True");
            Assert.That(RestrictionAnalyzer.AreCompatible(RestrictionAnalyzer.GetEnzyme("SfiI")!, gtc), Is.True,
                "NNN intersects GTC at every position, so some pair of SfiI/GTC-3p sites ligates");
            Assert.That(gtc.OverhangType, Is.EqualTo(OverhangType.ThreePrime));
        });
    }

    // F30 — The reverse-complement branch. Two same-polarity overhangs also ligate when one is the
    //       reverse complement of the other, because each cut offers BOTH W and revcomp(W) (the
    //       downstream fragment's end and the upstream fragment's end). Unreachable from the
    //       built-in table — every built-in window is its own reverse complement — but real in
    //       REBASE: `BsmI` ("CN", 3') and `BsrI` ("GN", 3') are compatible only through it
    //       (BsmI can leave a 3' GG and BsrI a 3' CC), one of 12 such REBASE pairs found this
    //       session. Biopython's plain `ovhgseq` equality misses all 12.
    [Test]
    public void AreCompatible_ReverseComplementaryOverhangs_AreCompatible()
    {
        var ca = new RestrictionEnzyme("CA-5p", "AAACAAA", 3, 5, "test");   // 5' "CA"
        var tg = new RestrictionEnzyme("TG-5p", "AAATGAA", 3, 5, "test");   // 5' "TG" = revcomp("CA")
        var ga = new RestrictionEnzyme("GA-5p", "AAAGAAA", 3, 5, "test");   // 5' "GA"
        var cn = new RestrictionEnzyme("CN-3p", "AAACNAA", 5, 3, "test");   // 3' "CN"
        var gn = new RestrictionEnzyme("GN-3p", "AAAGNAA", 5, 3, "test");   // 3' "GN"

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(ca), Is.EqualTo("CA"));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(tg), Is.EqualTo("TG"));

            Assert.That(RestrictionAnalyzer.AreCompatible(ca, tg), Is.True,
                "5' CA and 5' TG are reverse complements, so the two ends anneal");
            Assert.That(RestrictionAnalyzer.AreCompatible(tg, ca), Is.True, "and symmetrically");
            Assert.That(RestrictionAnalyzer.AreCompatible(ca, ga), Is.False,
                "5' CA matches neither GA nor revcomp(GA) = TC");

            // The REBASE BsmI/BsrI shape: 3' "CN" with 3' "GN" — revcomp("GN") = "NC", which
            // intersects "CN" position by position.
            Assert.That(RestrictionAnalyzer.AreCompatible(cn, gn), Is.True,
                "BsmI(3' CN)/BsrI(3' GN) are compatible through the reverse-complement branch");

            // Polarity still dominates: the same letters with opposite polarity stay incompatible.
            Assert.That(RestrictionAnalyzer.AreCompatible(ca, new RestrictionEnzyme("TG-3p", "AAATGAA", 5, 3, "test")),
                Is.False, "5' CA cannot anneal to a 3' TG");
        });
    }

    // F29 — Null enzyme definitions are rejected rather than silently answering false.
    [Test]
    public void AreCompatible_NullEnzymeDefinition_Throws()
    {
        var ecoRI = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.AreCompatible(null!, ecoRI));
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.AreCompatible(ecoRI, null!));
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.GetOverhangSequence((RestrictionEnzyme)null!));
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.GetEnzyme(null!));
        });
    }

    #endregion

    #region FindCompatibleEnzymes

    // F29 — The whole-table enumeration: 64 unordered pairs on this table, 45 of them blunt×blunt
    //       (10 blunt cutters → C(10,2) = 45) and 19 sticky, grouped by the shared overhang
    //       GATC (4 enzymes → 6 pairs), CTAG (4 → 6), CG (4 → 6) and TCGA (2 → 1).
    //       Every verdict agrees with Biopython's `%` operator (1521/1521 ordered pairs).
    [Test]
    public void FindCompatibleEnzymes_WholeTable_MatchesTheReferenceGrouping()
    {
        var pairs = RestrictionAnalyzer.FindCompatibleEnzymes().ToList();
        var byEnd = pairs.GroupBy(p => p.CompatibleEnd).ToDictionary(g => g.Key, g => g.Count());

        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(64));
            Assert.That(byEnd["blunt"], Is.EqualTo(45), "C(10,2) blunt pairs");
            Assert.That(byEnd["GATC"], Is.EqualTo(6), "BamHI/BglII/MboI/Sau3AI");
            Assert.That(byEnd["CTAG"], Is.EqualTo(6), "NheI/XbaI/SpeI/AvrII");
            Assert.That(byEnd["CG"], Is.EqualTo(6), "ClaI/TaqI/MspI/HpaII");
            Assert.That(byEnd["TCGA"], Is.EqualTo(1), "SalI/XhoI");
            Assert.That(byEnd.Keys, Has.Count.EqualTo(5), "no other shared end exists in the table");

            // No pair may be self-paired, and every reported pair must re-verify.
            Assert.That(pairs.All(p => p.Enzyme1 != p.Enzyme2), Is.True);
            Assert.That(pairs.All(p => RestrictionAnalyzer.AreCompatible(p.Enzyme1, p.Enzyme2)), Is.True);
        });
    }

    // F32 — The per-enzyme query (Biopython's `compatible_end()` notion, without the reference's
    //       own 5'-overhang-only filter): it includes the enzyme itself, so an enzyme with a unique
    //       overhang is no longer invisible. Exact sets from Biopython's `%` over the 39 built-ins.
    [Test]
    [TestCase("BamHI", "BamHI,BglII,MboI,Sau3AI")]
    [TestCase("Sau3AI", "BamHI,BglII,MboI,Sau3AI")]
    [TestCase("SalI", "SalI,XhoI")]
    [TestCase("TaqI", "ClaI,HpaII,MspI,TaqI")]
    [TestCase("SfiI", "SfiI")]
    [TestCase("PstI", "PstI")]
    [TestCase("KpnI", "KpnI")]
    [TestCase("EcoRV", "AluI,DpnI,EcoRV,HaeIII,HincII,RsaI,ScaI,SmaI,StuI,SwaI")]
    public void FindCompatibleEnzymes_PerEnzyme_MatchesTheReferenceSet(string enzyme, string expected)
    {
        var partners = RestrictionAnalyzer.FindCompatibleEnzymes(enzyme).Select(t => t.Enzyme).ToList();

        Assert.That(partners, Is.EquivalentTo(expected.Split(',')),
            $"{enzyme}'s compatible set must match Biopython's `%` over the built-in table");
    }

    // F32 — The per-enzyme query is consistent with the pairwise predicate over the whole table,
    //       and the predicate is symmetric for all 39 × 39 ordered pairs.
    [Test]
    public void FindCompatibleEnzymes_PerEnzyme_AgreesWithAreCompatible_AndIsSymmetric()
    {
        var names = RestrictionAnalyzer.Enzymes.Values.Select(e => e.Name).ToList();

        Assert.Multiple(() =>
        {
            foreach (string a in names)
            {
                var reported = RestrictionAnalyzer.FindCompatibleEnzymes(a).Select(t => t.Enzyme).ToHashSet();
                var expected = names.Where(b => RestrictionAnalyzer.AreCompatible(a, b)).ToHashSet();
                Assert.That(reported, Is.EquivalentTo(expected), $"{a}: per-enzyme query vs predicate");

                foreach (string b in names)
                {
                    Assert.That(RestrictionAnalyzer.AreCompatible(a, b),
                        Is.EqualTo(RestrictionAnalyzer.AreCompatible(b, a)), $"symmetry for {a}/{b}");
                }
            }
        });
    }

    [Test]
    public void FindCompatibleEnzymes_PerEnzyme_UnknownName_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => RestrictionAnalyzer.FindCompatibleEnzymes("Nope").ToList());
    }

    #endregion

    #region Isoschizomers / equischizomers / neoschizomers

    // F33 — Biopython definitions: equischizomer = same site AND same cut positions;
    //       neoschizomer = same site, different cut; isoschizomer = the union, excluding the
    //       enzyme itself. Expected sets are Biopython's own, computed with `batch=` restricted to
    //       these 39 enzymes this session.
    [Test]
    [TestCase("MspI", "HpaII", "HpaII", "")]
    [TestCase("HpaII", "MspI", "MspI", "")]
    [TestCase("MboI", "DpnI,Sau3AI", "Sau3AI", "DpnI")]
    [TestCase("Sau3AI", "DpnI,MboI", "MboI", "DpnI")]
    [TestCase("DpnI", "MboI,Sau3AI", "", "MboI,Sau3AI")]
    [TestCase("EcoRI", "", "", "")]
    [TestCase("HincII", "", "", "")]
    [TestCase("SfiI", "", "", "")]
    public void FindIsoschizomers_MatchesBiopythonWithinTheBuiltInTable(
        string enzyme, string iso, string equi, string neo)
    {
        static string[] Expected(string csv) => csv.Length == 0 ? [] : csv.Split(',');
        static string[] Actual(IEnumerable<RestrictionEnzyme> result) =>
            result.Select(e => e.Name).OrderBy(n => n, StringComparer.Ordinal).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(Actual(RestrictionAnalyzer.FindIsoschizomers(enzyme)), Is.EquivalentTo(Expected(iso)));
            Assert.That(Actual(RestrictionAnalyzer.FindEquischizomers(enzyme)), Is.EquivalentTo(Expected(equi)));
            Assert.That(Actual(RestrictionAnalyzer.FindNeoschizomers(enzyme)), Is.EquivalentTo(Expected(neo)));
        });
    }

    // F33 — Structural invariants over the whole table: iso = equi ⊎ neo, the queried enzyme is
    //       never its own isoschizomer, an equischizomer always leaves the same end (so it is
    //       compatible), and a neoschizomer by definition does not cut the same way.
    [Test]
    public void FindIsoschizomers_WholeTable_PartitionAndEndConsistency()
    {
        Assert.Multiple(() =>
        {
            foreach (var enzyme in RestrictionAnalyzer.Enzymes.Values)
            {
                var iso = RestrictionAnalyzer.FindIsoschizomers(enzyme.Name).Select(e => e.Name).ToHashSet();
                var equi = RestrictionAnalyzer.FindEquischizomers(enzyme.Name).Select(e => e.Name).ToHashSet();
                var neo = RestrictionAnalyzer.FindNeoschizomers(enzyme.Name).Select(e => e.Name).ToHashSet();

                Assert.That(equi.Overlaps(neo), Is.False, $"{enzyme.Name}: equi and neo are disjoint");
                Assert.That(equi.Union(neo), Is.EquivalentTo(iso), $"{enzyme.Name}: iso = equi ∪ neo");
                Assert.That(iso, Does.Not.Contain(enzyme.Name), $"{enzyme.Name} is not its own isoschizomer");

                foreach (string same in equi)
                {
                    Assert.That(RestrictionAnalyzer.GetOverhangSequence(same),
                        Is.EqualTo(RestrictionAnalyzer.GetOverhangSequence(enzyme.Name)),
                        $"{enzyme.Name}/{same} are equischizomers, so they leave the same end");
                    Assert.That(RestrictionAnalyzer.AreCompatible(enzyme.Name, same), Is.True);
                }
            }
        });
    }

    [Test]
    public void FindIsoschizomers_UnknownName_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindIsoschizomers("Nope").ToList());
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindEquischizomers("Nope").ToList());
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindNeoschizomers("Nope").ToList());
        });
    }

    #endregion

    #region Table audit against the reference

    // F34 — The whole built-in table re-audited against Biopython 1.88 in one locked test: site
    //       string, size, both cut offsets (Biopython `fst5` and `size + fst3`), blunt/sticky,
    //       5'/3' polarity and overhang sequence. 39 rows, every value read off the oracle this
    //       session; 0 wrong entries were found.
    [Test]
    public void EnzymeTable_EveryEntry_MatchesTheReferenceRow()
    {
        // name, site, cutForward (fst5), cutReverse (size + fst3), overhangType, ovhgseq
        (string Name, string Site, int Cf, int Cr, OverhangType Type, string Ovhg)[] reference =
        [
            ("EcoRI", "GAATTC", 1, 5, OverhangType.FivePrime, "AATT"),
            ("BamHI", "GGATCC", 1, 5, OverhangType.FivePrime, "GATC"),
            ("HindIII", "AAGCTT", 1, 5, OverhangType.FivePrime, "AGCT"),
            ("XhoI", "CTCGAG", 1, 5, OverhangType.FivePrime, "TCGA"),
            ("SalI", "GTCGAC", 1, 5, OverhangType.FivePrime, "TCGA"),
            ("PstI", "CTGCAG", 5, 1, OverhangType.ThreePrime, "TGCA"),
            ("SphI", "GCATGC", 5, 1, OverhangType.ThreePrime, "CATG"),
            ("KpnI", "GGTACC", 5, 1, OverhangType.ThreePrime, "GTAC"),
            ("SacI", "GAGCTC", 5, 1, OverhangType.ThreePrime, "AGCT"),
            ("XbaI", "TCTAGA", 1, 5, OverhangType.FivePrime, "CTAG"),
            ("NcoI", "CCATGG", 1, 5, OverhangType.FivePrime, "CATG"),
            ("NdeI", "CATATG", 2, 4, OverhangType.FivePrime, "TA"),
            ("NheI", "GCTAGC", 1, 5, OverhangType.FivePrime, "CTAG"),
            ("SpeI", "ACTAGT", 1, 5, OverhangType.FivePrime, "CTAG"),
            ("AvrII", "CCTAGG", 1, 5, OverhangType.FivePrime, "CTAG"),
            ("ClaI", "ATCGAT", 2, 4, OverhangType.FivePrime, "CG"),
            ("AgeI", "ACCGGT", 1, 5, OverhangType.FivePrime, "CCGG"),
            ("NotI", "GCGGCCGC", 2, 6, OverhangType.FivePrime, "GGCC"),
            ("SfiI", "GGCCNNNNNGGCC", 8, 5, OverhangType.ThreePrime, "NNN"),
            ("PacI", "TTAATTAA", 5, 3, OverhangType.ThreePrime, "AT"),
            ("AscI", "GGCGCGCC", 2, 6, OverhangType.FivePrime, "CGCG"),
            ("FseI", "GGCCGGCC", 6, 2, OverhangType.ThreePrime, "CCGG"),
            ("SwaI", "ATTTAAAT", 4, 4, OverhangType.Blunt, ""),
            ("MspI", "CCGG", 1, 3, OverhangType.FivePrime, "CG"),
            ("HpaII", "CCGG", 1, 3, OverhangType.FivePrime, "CG"),
            ("TaqI", "TCGA", 1, 3, OverhangType.FivePrime, "CG"),
            ("AluI", "AGCT", 2, 2, OverhangType.Blunt, ""),
            ("RsaI", "GTAC", 2, 2, OverhangType.Blunt, ""),
            ("HaeIII", "GGCC", 2, 2, OverhangType.Blunt, ""),
            ("DpnI", "GATC", 2, 2, OverhangType.Blunt, ""),
            ("Sau3AI", "GATC", 0, 4, OverhangType.FivePrime, "GATC"),
            ("MboI", "GATC", 0, 4, OverhangType.FivePrime, "GATC"),
            ("EcoRV", "GATATC", 3, 3, OverhangType.Blunt, ""),
            ("SmaI", "CCCGGG", 3, 3, OverhangType.Blunt, ""),
            ("HincII", "GTYRAC", 3, 3, OverhangType.Blunt, ""),
            ("ScaI", "AGTACT", 3, 3, OverhangType.Blunt, ""),
            ("StuI", "AGGCCT", 3, 3, OverhangType.Blunt, ""),
            ("BglII", "AGATCT", 1, 5, OverhangType.FivePrime, "GATC"),
            ("ApaI", "GGGCCC", 5, 1, OverhangType.ThreePrime, "GGCC"),
        ];

        Assert.That(reference, Has.Length.EqualTo(RestrictionAnalyzer.Enzymes.Count),
            "the reference row set must cover the whole built-in table");

        Assert.Multiple(() =>
        {
            foreach (var row in reference)
            {
                var enzyme = RestrictionAnalyzer.GetEnzyme(row.Name);
                Assert.That(enzyme, Is.Not.Null, $"{row.Name} must be in the table");
                Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(row.Site), $"{row.Name} site");
                Assert.That(enzyme.RecognitionLength, Is.EqualTo(row.Site.Length), $"{row.Name} size");
                Assert.That(enzyme.CutPositionForward, Is.EqualTo(row.Cf), $"{row.Name} fst5");
                Assert.That(enzyme.CutPositionReverse, Is.EqualTo(row.Cr), $"{row.Name} size+fst3");
                Assert.That(enzyme.OverhangType, Is.EqualTo(row.Type), $"{row.Name} polarity");
                Assert.That(enzyme.IsBluntEnd, Is.EqualTo(row.Type == OverhangType.Blunt), $"{row.Name} blunt");
                Assert.That(RestrictionAnalyzer.GetOverhangSequence(row.Name), Is.EqualTo(row.Ovhg),
                    $"{row.Name} ovhgseq");
            }
        });
    }

    #endregion
}
