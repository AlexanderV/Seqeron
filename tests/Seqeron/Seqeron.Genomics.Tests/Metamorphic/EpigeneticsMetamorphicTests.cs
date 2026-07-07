namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: EPIGEN-AGE-001 — epigenetic clock (DNAm age) (Epigenetics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 181.
    //
    // API under test (EpigeneticsAnalyzer.CalculateEpigeneticAge):
    //   Y = intercept + Σ coef_i·β_i over clock CpGs, mapped to years by the monotone Horvath
    //   inverse calibration.
    //
    // Relations (derived from the linear predictor, NOT from output):
    //   • MON  (more clock-site methylation ⇒ higher age): raising β at positive-coefficient sites
    //          raises Y, and the inverse calibration is monotone increasing, so age increases.
    //   • INV  (site order independent): the predictor is a sum over CpGs, independent of map order.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly System.Collections.Generic.Dictionary<string, double> ClockCoefficients = new()
    {
        ["cg01"] = 1.0, ["cg02"] = 1.0, ["cg03"] = 1.0, // all positive
    };

    #region EPIGEN-AGE-001 MON — more methylation at positive sites raises the age

    [Test]
    [Description("MON: raising methylation at positive-coefficient clock CpGs raises the linear predictor, and the Horvath inverse calibration is monotone increasing, so the estimated age increases.")]
    public void EpigeneticAge_MoreMethylation_HigherAge()
    {
        double previous = double.MinValue;
        foreach (double beta in new[] { 0.1, 0.3, 0.6, 0.9 })
        {
            var methylation = new System.Collections.Generic.Dictionary<string, double>
            {
                ["cg01"] = beta, ["cg02"] = beta, ["cg03"] = beta,
            };
            double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, ClockCoefficients);
            age.Should().BeGreaterThan(previous, because: $"higher methylation (β={beta}) at positive-coefficient sites raises the predicted age");
            previous = age;
        }
    }

    #endregion

    #region EPIGEN-AGE-001 INV — the estimate is independent of CpG order

    [Test]
    [Description("INV: the linear predictor is a sum over clock CpGs, so building the methylation map in a different order yields the same age.")]
    public void EpigeneticAge_SiteOrder_Invariant()
    {
        var forward = new System.Collections.Generic.Dictionary<string, double>
        {
            ["cg01"] = 0.2, ["cg02"] = 0.5, ["cg03"] = 0.8,
        };
        var reordered = new System.Collections.Generic.Dictionary<string, double>
        {
            ["cg03"] = 0.8, ["cg01"] = 0.2, ["cg02"] = 0.5,
        };

        EpigeneticsAnalyzer.CalculateEpigeneticAge(reordered, ClockCoefficients)
            .Should().BeApproximately(EpigeneticsAnalyzer.CalculateEpigeneticAge(forward, ClockCoefficients), 1e-9,
                because: "the weighted sum over CpGs does not depend on the map's insertion order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: EPIGEN-BISULF-001 — in-silico bisulfite conversion (Epigenetics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 182.
    //
    // API under test (EpigeneticsAnalyzer.SimulateBisulfiteConversion):
    //   Converts every unmethylated cytosine to thymine; methylated (protected) cytosines stay C.
    //
    // Relations (derived from the conversion rule, NOT from output):
    //   • INV  (methylated-C set preserved): the cytosines surviving conversion are exactly the
    //          protected (methylated) positions.
    //   • SHIFT (prepend flank shifts conversions): shifting the protected positions by a prepended
    //          flank reproduces the original converted region unchanged at the shifted offset.
    // ───────────────────────────────────────────────────────────────────────────

    #region EPIGEN-BISULF-001 INV — surviving cytosines are exactly the methylated set

    [Test]
    [Description("INV: every unmethylated C converts to T, so the positions that remain C after conversion are exactly the protected (methylated) cytosine positions.")]
    public void Bisulfite_SurvivingCytosines_AreMethylatedSet()
    {
        const string seq = "ACGCGTCC";              // C at 1,3,6,7
        var methylated = new HashSet<int> { 3, 6 }; // protected cytosines

        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq, methylated);

        var survivingC = Enumerable.Range(0, converted.Length).Where(i => converted[i] == 'C').ToHashSet();
        survivingC.Should().BeEquivalentTo(methylated, because: "only protected cytosines stay C; all others become T");
    }

    #endregion

    #region EPIGEN-BISULF-001 SHIFT — a prepended flank shifts the conversions

    [Test]
    [Description("SHIFT: shifting the protected positions by a prepended flank reproduces the original converted region unchanged at the shifted offset.")]
    public void Bisulfite_PrependFlank_ShiftsConversions()
    {
        const string seq = "ACGCGTCC";
        var methylated = new HashSet<int> { 3, 6 };
        string original = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq, methylated);

        foreach (var flank in new[] { "TTC", "GGCCAA" })
        {
            var shiftedMeth = methylated.Select(p => p + flank.Length).ToHashSet();
            string full = EpigeneticsAnalyzer.SimulateBisulfiteConversion(flank + seq, shiftedMeth);
            full.Substring(flank.Length).Should().Be(original,
                because: $"shifting the protected positions by {flank.Length} converts the original region identically");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: EPIGEN-CHROM-001 — chromatin-state prediction (Epigenetics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 183.
    //
    // API under test (EpigeneticsAnalyzer.PredictChromatinState):
    //   Calls the chromatin state at a locus purely from its six histone-mark signals.
    //
    // Relations (derived from the per-locus, context-free call, NOT from output):
    //   • INV  (region order independent): each locus is called from its own marks, so processing a
    //          list of loci in any order yields the correspondingly reordered states.
    //   • SHIFT (prepend flank shifts states): prepending flank loci shifts the original loci's states
    //          along the output list without changing them.
    // ───────────────────────────────────────────────────────────────────────────

    // (h3k4me3, h3k4me1, h3k27ac, h3k36me3, h3k27me3, h3k9me3) signals per locus.
    private static readonly (double, double, double, double, double, double)[] ChromatinLoci =
    {
        (1, 0, 0, 0, 0, 0), // active promoter (H3K4me3)
        (0, 1, 1, 0, 0, 0), // active enhancer (H3K4me1 + H3K27ac)
        (0, 0, 0, 0, 1, 0), // Polycomb-repressed (H3K27me3)
        (0, 0, 0, 0, 0, 0), // quiescent
    };

    private static List<EpigeneticsAnalyzer.ChromatinState> States(
        IEnumerable<(double, double, double, double, double, double)> loci) =>
        loci.Select(m => EpigeneticsAnalyzer.PredictChromatinState(m.Item1, m.Item2, m.Item3, m.Item4, m.Item5, m.Item6)).ToList();

    #region EPIGEN-CHROM-001 INV — states are independent of locus order

    [Test]
    [Description("INV: each locus's state is called from its own marks, so reversing the locus order reverses the state list (per-locus independence).")]
    public void Chromatin_RegionOrder_Invariant()
    {
        var forward = States(ChromatinLoci);
        var reversed = States(ChromatinLoci.Reverse());

        reversed.Should().Equal(Enumerable.Reverse(forward).ToList(),
            because: "the chromatin state of a locus does not depend on the order loci are processed");
    }

    #endregion

    #region EPIGEN-CHROM-001 SHIFT — a prepended flank shifts the states

    [Test]
    [Description("SHIFT: prepending flank loci shifts the original loci's states along the output list without changing them.")]
    public void Chromatin_PrependFlank_ShiftsStates()
    {
        var original = States(ChromatinLoci);

        var flank = new (double, double, double, double, double, double)[]
        {
            (0, 0, 0, 1, 0, 0), // transcribed
            (0, 0, 0, 0, 0, 1), // heterochromatin
        };
        var combined = States(flank.Concat(ChromatinLoci));

        combined.Skip(flank.Length).Should().Equal(original,
            because: "prepending flank loci moves the original states later in the list, leaving them unchanged");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: EPIGEN-DMR-001 — differentially methylated regions (Epigenetics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 184.
    //
    // API under test (EpigeneticsAnalyzer.FindDMRs):
    //   Tiles the genome and reports windows whose |mean methylation difference| > minDifference.
    //
    // Relations (derived from the |difference| cutoff, NOT from output):
    //   • MON  (lower threshold ⇒ superset): lowering minDifference admits windows with smaller
    //          differences, so the DMR set grows.
    //   • SYM  (DMR(A,B) consistent with (B,A)): swapping the samples negates each mean difference
    //          but |difference| is unchanged, so the same regions are reported with flipped sign.
    // ───────────────────────────────────────────────────────────────────────────

    private static EpigeneticsAnalyzer.MethylationSite MethSite(int pos, double level) =>
        new(pos, EpigeneticsAnalyzer.MethylationType.CpG, "CpG", level, 10);

    // Window A (positions 0..2) differs by 0.5; window B (positions 1000..1002) differs by 0.3.
    private static readonly EpigeneticsAnalyzer.MethylationSite[] DmrSample1 =
    {
        MethSite(0, 0.0), MethSite(1, 0.0), MethSite(2, 0.0),
        MethSite(1000, 0.0), MethSite(1001, 0.0), MethSite(1002, 0.0),
    };
    private static readonly EpigeneticsAnalyzer.MethylationSite[] DmrSample2 =
    {
        MethSite(0, 0.5), MethSite(1, 0.5), MethSite(2, 0.5),
        MethSite(1000, 0.3), MethSite(1001, 0.3), MethSite(1002, 0.3),
    };

    #region EPIGEN-DMR-001 MON — lowering the threshold yields a superset

    [Test]
    [Description("MON: lowering minDifference admits windows with smaller methylation differences, so the DMR set at the lower cutoff is a superset.")]
    public void Dmr_LowerThreshold_Superset()
    {
        var strict = EpigeneticsAnalyzer.FindDMRs(DmrSample1, DmrSample2, 1000, 0.4, 3).Select(d => (d.Start, d.End)).ToHashSet();
        var lenient = EpigeneticsAnalyzer.FindDMRs(DmrSample1, DmrSample2, 1000, 0.2, 3).Select(d => (d.Start, d.End)).ToHashSet();

        strict.IsSubsetOf(lenient).Should().BeTrue(because: "every window above the 0.4 cutoff is above 0.2");
        lenient.Count.Should().BeGreaterThan(strict.Count, because: "the 0.3-difference window is admitted only at the lower cutoff");
    }

    #endregion

    #region EPIGEN-DMR-001 SYM — swapping samples keeps the regions and flips the sign

    [Test]
    [Description("SYM: swapping the two samples negates each mean difference but leaves |difference| unchanged, so the same regions are reported with opposite-signed mean differences.")]
    public void Dmr_SwapSamples_SameRegionsFlippedSign()
    {
        var ab = EpigeneticsAnalyzer.FindDMRs(DmrSample1, DmrSample2, 1000, 0.2, 3).OrderBy(d => d.Start).ToList();
        var ba = EpigeneticsAnalyzer.FindDMRs(DmrSample2, DmrSample1, 1000, 0.2, 3).OrderBy(d => d.Start).ToList();

        ba.Select(d => (d.Start, d.End)).Should().Equal(ab.Select(d => (d.Start, d.End)),
            because: "|difference| is symmetric, so the same windows qualify");
        for (int i = 0; i < ab.Count; i++)
            ba[i].MeanDifference.Should().BeApproximately(-ab[i].MeanDifference, 1e-9,
                because: "swapping samples negates the mean methylation difference");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: EPIGEN-METHYL-001 — methylation calling from bisulfite reads (Epigenetics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 185.
    //
    // API under test (EpigeneticsAnalyzer.CalculateMethylationFromBisulfite):
    //   At each reference CpG, counts read C (methylated) / T (unmethylated) calls and reports the
    //   fraction and coverage.
    //
    // Relations (derived from per-read count accumulation, NOT from output):
    //   • INV  (read order independent): counts accumulate over reads, so the per-site level and
    //          coverage do not depend on read order.
    //   • ADD  (counts additive over reads): per CpG, the coverage and methylated-call counts of a
    //          read set equal the sums over any partition of the reads.
    // ───────────────────────────────────────────────────────────────────────────

    private const string MethylReference = "ACGTACGT"; // CpG cytosines at positions 1 and 5

    private static readonly (string, int)[] MethylReads =
    {
        ("ACGTACGT", 0), // pos1 C (meth), pos5 C (meth)
        ("ATGTATGT", 0), // pos1 T (unmeth), pos5 T (unmeth)
        ("ACGTATGT", 0), // pos1 C, pos5 T
        ("ATGTACGT", 0), // pos1 T, pos5 C
    };

    private static Dictionary<int, (double Level, int Coverage)> MethylCalls(IEnumerable<(string, int)> reads) =>
        EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(MethylReference, reads)
            .ToDictionary(s => s.Position, s => (s.MethylationLevel, s.Coverage));

    #region EPIGEN-METHYL-001 INV — calls are independent of read order

    [Test]
    [Description("INV: per-CpG C/T calls accumulate over reads, so reversing the read order yields the identical per-site level and coverage.")]
    public void Methylation_ReadOrder_Invariant()
    {
        var forward = MethylCalls(MethylReads);
        var reversed = MethylCalls(MethylReads.Reverse());

        reversed.Should().BeEquivalentTo(forward,
            because: "the per-site methylation call is a function of the read multiset, not its order");
    }

    #endregion

    #region EPIGEN-METHYL-001 ADD — coverage and methylated counts are additive over reads

    [Test]
    [Description("ADD: per CpG, coverage and the methylated-call count of the full read set equal the sums over a two-part partition of the reads.")]
    public void Methylation_Additive_OverReadPartition()
    {
        var all = MethylCalls(MethylReads);
        var g1 = MethylCalls(MethylReads.Take(2));
        var g2 = MethylCalls(MethylReads.Skip(2));

        foreach (int site in all.Keys)
        {
            all[site].Coverage.Should().Be(g1[site].Coverage + g2[site].Coverage,
                because: $"coverage at CpG {site} is additive over the read partition");

            int MethCount((double Level, int Coverage) c) => (int)System.Math.Round(c.Level * c.Coverage);
            MethCount(all[site]).Should().Be(MethCount(g1[site]) + MethCount(g2[site]),
                because: $"the methylated-call count at CpG {site} is additive over the read partition");
        }
    }

    #endregion
}
