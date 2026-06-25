namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the MolTools area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("MolTools")]
public class MolToolsCombinatorialTests
{
    /// <summary>Deterministic well-mixed ACGT sequence (LCG) so PAMs/sites occur at many positions.</summary>
    private static string DiverseDna(int n, uint seed = 0x12345678u)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CRISPR-PAM-001 — PAM-site discovery (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 18.
    // Dimensions: pamType(3: SpCas9/Cas12a/SaCas9) × strand(3: +/−/both) × seqLen(3).
    //             Full grid 3×3×3 = 27 cells.
    //
    // Model (Jinek et al. 2012; Zetsche et al. 2015): a PAM site is a position
    // whose IUPAC PAM motif matches one strand AND whose protospacer (guide-length)
    // window lies fully in the sequence. The PAM's placement relative to the target
    // differs by system: SpCas9 (NGG) and SaCas9 (NNGRRT) put the PAM 3′ of the
    // target; Cas12a (TTTV) puts it 5′. FindPamSites scans BOTH strands, reporting
    // forward-strand coordinates with an IsForwardStrand flag.
    //
    // The combinatorial point: pamType and strand interact — the PAM pattern and
    // its target side change with the system, and the strand selector partitions
    // the hits. Every reported site's PAM must satisfy the system's IUPAC motif
    // (on the strand it was found), the target window must have guide length, and
    // the forward/reverse partition must reproduce an independent strand-wise count.
    // ═══════════════════════════════════════════════════════════════════════

    public enum Strand { Forward, Reverse, Both }

    private static bool PamMatches(string seq, int i, string pam)
    {
        if (i + pam.Length > seq.Length) return false;
        for (int j = 0; j < pam.Length; j++)
            if (!IupacHelper.MatchesIupac(seq[i + j], pam[j]))
                return false;
        return true;
    }

    /// <summary>Independent count of PAM hits on one strand string, mirroring the production target-window bound.</summary>
    private static int BruteForcePamCount(string strandSeq, CrisprSystem system)
    {
        int count = 0, pamLen = system.PamSequence.Length, g = system.GuideLength;
        for (int i = 0; i + pamLen <= strandSeq.Length; i++)
        {
            if (!PamMatches(strandSeq, i, system.PamSequence)) continue;
            int targetStart = system.PamAfterTarget ? i - g : i + pamLen;
            int targetEnd = system.PamAfterTarget ? i - 1 : targetStart + g - 1;
            if (targetStart >= 0 && targetEnd < strandSeq.Length) count++;
        }
        return count;
    }

    [Test, Combinatorial]
    public void CrisprPam_SitesValid_StrandPartitionMatchesCount(
        [Values(CrisprSystemType.SpCas9, CrisprSystemType.Cas12a, CrisprSystemType.SaCas9)] CrisprSystemType pamType,
        [Values(Strand.Forward, Strand.Reverse, Strand.Both)] Strand strand,
        [Values(30, 80, 200)] int seqLen)
    {
        string seq = DiverseDna(seqLen);
        var dna = new DnaSequence(seq);
        var system = CrisprDesigner.GetSystem(pamType);

        var all = CrisprDesigner.FindPamSites(dna, pamType).ToList();
        var selected = strand switch
        {
            Strand.Forward => all.Where(s => s.IsForwardStrand).ToList(),
            Strand.Reverse => all.Where(s => !s.IsForwardStrand).ToList(),
            _ => all,
        };

        // Strand partition is exhaustive and disjoint.
        all.Count(s => s.IsForwardStrand).Should().Be(BruteForcePamCount(seq, system));
        all.Count(s => !s.IsForwardStrand).Should().Be(BruteForcePamCount(RevComp(seq), system));
        (all.Count(s => s.IsForwardStrand) + all.Count(s => !s.IsForwardStrand)).Should().Be(all.Count);

        foreach (var site in selected)
        {
            site.Position.Should().BeInRange(0, seqLen - system.PamSequence.Length);
            site.TargetSequence.Length.Should().Be(system.GuideLength, "the protospacer is guide-length");
            site.System.PamSequence.Should().Be(system.PamSequence);

            // The PAM, read on the strand it was found on, satisfies the IUPAC motif.
            string pamOnStrand = site.IsForwardStrand ? site.PamSequence : RevComp(site.PamSequence);
            for (int j = 0; j < system.PamSequence.Length; j++)
                IupacHelper.MatchesIupac(pamOnStrand[j], system.PamSequence[j]).Should().BeTrue(
                    $"PAM char {j} of \"{pamOnStrand}\" must match motif {system.PamSequence}");

            if (site.IsForwardStrand)
                site.PamSequence.Should().Be(seq.Substring(site.Position, system.PamSequence.Length));
        }
    }

    /// <summary>
    /// Interaction witness: the PAM motif's target side depends on the system —
    /// SpCas9 (NGG) places the protospacer 5′ of the PAM, Cas12a (TTTV) places it
    /// 3′. Verified on a forward-strand hit of each system.
    /// </summary>
    [Test]
    public void CrisprPam_TargetSide_DependsOnSystem()
    {
        var seq = new DnaSequence(DiverseDna(120));

        var cas9 = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9).First(s => s.IsForwardStrand);
        cas9.TargetStart.Should().Be(cas9.Position - cas9.System.GuideLength, "Cas9 PAM is 3′ of the target");

        var cas12a = CrisprDesigner.FindPamSites(seq, CrisprSystemType.Cas12a).First(s => s.IsForwardStrand);
        cas12a.TargetStart.Should().Be(cas12a.Position + cas12a.System.PamSequence.Length, "Cas12a PAM is 5′ of the target");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CRISPR-GUIDE-001 — Guide-RNA design & specificity (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 19.
    // Dimensions: pamType(3) × guideLen(3) × maxOff(3) × scoringMethod(2).
    //             Full grid 3×3×3×2 = 54 cells.
    //
    // Note on the grid: guideLen is system-determined (SpCas9 20, SaCas9 21,
    // Cas12a 23), so it is not an independent axis — it co-varies with pamType.
    // The on-target scoringMethod (Doench 2014 / Rule Set 2) is defined only for
    // the 20-nt SpCas9 model, so it cannot be crossed with non-Cas9 systems.
    // The row is therefore exercised by two coherent grids: the design/off-target
    // pipeline over pamType × maxOff × seqLen, and on-target scoring over
    // scoringMethod × context.
    //
    // Model: a designed guide is the protospacer of a PAM site; its off-target
    // burden is the set of OTHER PAM-flanked windows within maxMismatches of it
    // (Hsu et al. 2013). Off-targets exclude the exact (0-mismatch) on-target and
    // grow monotonically with maxMismatches.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void CrisprGuide_OffTargetBurden_ValidAndMonotone(
        [Values(CrisprSystemType.SpCas9, CrisprSystemType.SaCas9, CrisprSystemType.Cas12a)] CrisprSystemType pamType,
        [Values(1, 2, 3)] int maxOff,
        [Values(60, 120, 240)] int seqLen)
    {
        var system = CrisprDesigner.GetSystem(pamType);
        var genome = new DnaSequence(DiverseDna(seqLen));
        string guide = genome.Sequence.Substring(0, system.GuideLength);

        var offs = CrisprDesigner.FindOffTargets(guide, genome, maxOff, pamType).ToList();

        foreach (var ot in offs)
        {
            ot.Mismatches.Should().BeInRange(1, maxOff, "off-targets exclude the exact match and respect the cap");
            ot.Sequence.Length.Should().Be(system.GuideLength);
            CountMm(guide, ot.Sequence).Should().Be(ot.Mismatches);
            ot.MismatchPositions.Count.Should().Be(ot.Mismatches);
        }

        // Monotone growth in the mismatch budget (maxOff interaction).
        var atCap = offs.Select(o => (o.Position, o.IsForwardStrand)).ToHashSet();
        if (maxOff < 3)
        {
            var bigger = CrisprDesigner.FindOffTargets(guide, genome, maxOff + 1, pamType)
                .Select(o => (o.Position, o.IsForwardStrand)).ToHashSet();
            atCap.Should().BeSubsetOf(bigger, "allowing more mismatches cannot remove an off-target");
        }

        CrisprDesigner.CalculateSpecificityScore(guide, genome, pamType).Should().BeInRange(0.0, 100.0);
    }

    private static int CountMm(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    public enum OnTargetMethod { Doench2014, RuleSet2 }

    [Test, Combinatorial]
    public void CrisprGuide_OnTargetScore_RangeAndDeterminism(
        [Values(OnTargetMethod.Doench2014, OnTargetMethod.RuleSet2)] OnTargetMethod method,
        [Values("GC", "AT", "MIX")] string contextKind)
    {
        // 30-mer = 4 upstream + 20 protospacer + NGG PAM (offsets 25-26 = GG) + 3 downstream.
        string proto = contextKind switch
        {
            "GC" => "GCGCGCGCGCGCGCGCGCGC",
            "AT" => "ATATATATATATATATATAT",
            _ => "ACGTACGTACGTACGTACGT",
        };
        string context = "ACGT" + proto + "AGG" + "ACG"; // 4+20+3+3 = 30, NGG PAM

        (double lo, double hi) = method == OnTargetMethod.Doench2014 ? (0.0, 100.0) : (0.0, 1.0);
        double Score(string c) => method == OnTargetMethod.Doench2014
            ? CrisprDesigner.CalculateOnTargetDoench2014(c)
            : CrisprDesigner.CalculateOnTargetRuleSet2(c);

        double s = Score(context);
        s.Should().BeInRange(lo, hi, $"{method} score is bounded");
        Score(context).Should().Be(s, "scoring is deterministic");
    }

    /// <summary>
    /// Interaction witness: the two on-target models are sequence-sensitive — a
    /// GC-rich and an AT-rich protospacer do not receive identical scores.
    /// </summary>
    [Test]
    public void CrisprGuide_OnTargetScore_IsSequenceSensitive()
    {
        string Ctx(string p) => "ACGT" + p + "AGGACG";
        string gc = Ctx("GCGCGCGCGCGCGCGCGCGC");
        string at = Ctx("ATATATATATATATATATAT");

        CrisprDesigner.CalculateOnTargetDoench2014(gc).Should().NotBe(CrisprDesigner.CalculateOnTargetDoench2014(at));
        CrisprDesigner.CalculateOnTargetRuleSet2(gc).Should().NotBe(CrisprDesigner.CalculateOnTargetRuleSet2(at));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CRISPR-OFF-001 — Off-target scoring (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 20.
    // Dimensions: maxMismatch(3) × seedLen(3) × scoringMethod(2: MIT/CFD).
    //             Full grid 3×3×2 = 18 cells.
    //
    // Model (Hsu et al. 2013 MIT; Doench et al. 2016 CFD): an off-target hit score
    // rates how likely a mismatched site is still cut. An exact match is the
    // maximum (MIT 100, CFD 1.0) and every score is bounded ([0,100] / [0,1]).
    // CFD is a strict product of per-position mismatch factors ≤ 1, so it is
    // MONOTONE non-increasing as mismatches are added; MIT's published
    // position-weight vector is empirically irregular (some positions carry ~0
    // weight, and the peak is not the most PAM-proximal base), so MIT is NOT
    // monotone in mismatch position — only the aggregate specificity is monotone
    // in off-target burden.
    //
    // The combinatorial point: maxMismatch (count), seedLen (the PAM-proximal
    // window the mismatches fall in) and the scoring model interact. Within a seed
    // window, CFD must be monotone as mismatches accumulate, and a fully-mismatched
    // site must score below an exact match under either model.
    // ═══════════════════════════════════════════════════════════════════════

    public enum OffTargetMethod { Mit, Cfd }

    private const string Guide20 = "GACCTGCAGTACGTTGCAAC"; // arbitrary 20-nt guide (pos 0 = PAM-distal, 19 = PAM-proximal)

    private static double OffScore(OffTargetMethod method, string offTarget20) => method switch
    {
        OffTargetMethod.Mit => CrisprDesigner.CalculateMitHitScore(Guide20, offTarget20),
        OffTargetMethod.Cfd => CrisprDesigner.CalculateCfdScore(Guide20, offTarget20, "AGG"),
        _ => throw new ArgumentOutOfRangeException(nameof(method)),
    };

    private static string Mutate(string s, IEnumerable<int> positions)
    {
        char[] c = s.ToCharArray();
        foreach (int p in positions) c[p] = c[p] == 'A' ? 'C' : 'A';
        return new string(c);
    }

    [Test, Combinatorial]
    public void CrisprOff_ExactIsMax_CfdMonotoneInSeed(
        [Values(1, 2, 3)] int maxMismatch,
        [Values(3, 6, 12)] int seedLen,
        [Values(OffTargetMethod.Mit, OffTargetMethod.Cfd)] OffTargetMethod method)
    {
        double max = method == OffTargetMethod.Mit ? 100.0 : 1.0;

        // Exact match is the maximum, and scoring is deterministic.
        double exact = OffScore(method, Guide20);
        exact.Should().BeApproximately(max, 1e-9);
        OffScore(method, Guide20).Should().Be(exact);

        // Accumulate up to maxMismatch mismatches inside the PAM-proximal seed window.
        int seedStart = 20 - seedLen;
        var seedPositions = Enumerable.Range(seedStart, seedLen).Take(maxMismatch).ToList();

        double prev = max;
        foreach (int k in Enumerable.Range(1, seedPositions.Count))
        {
            double s = OffScore(method, Mutate(Guide20, seedPositions.Take(k)));
            s.Should().BeInRange(0.0, max);
            if (method == OffTargetMethod.Cfd)
                s.Should().BeLessThanOrEqualTo(prev + 1e-12, "CFD is a product of per-position factors ≤ 1, hence monotone");
            prev = s;
        }

        // A fully-mismatched site is below an exact match under either model.
        OffScore(method, Mutate(Guide20, Enumerable.Range(0, 20))).Should().BeLessThan(max);
    }

    /// <summary>
    /// Interaction witness: the aggregate MIT specificity score (Hsu 2013 guide
    /// score) decreases monotonically as off-target burden grows.
    /// </summary>
    [Test]
    public void CrisprOff_MitAggregateSpecificity_DecreasesWithBurden()
    {
        double none = CrisprDesigner.CalculateMitSpecificityScore(Array.Empty<double>());
        double few = CrisprDesigner.CalculateMitSpecificityScore(new[] { 10.0, 5.0 });
        double many = CrisprDesigner.CalculateMitSpecificityScore(new[] { 10.0, 5.0, 40.0, 30.0 });

        none.Should().Be(100.0, "no off-targets ⇒ perfect specificity");
        few.Should().BeLessThan(none);
        many.Should().BeLessThan(few);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-TM-001 — Primer melting-temperature calculation (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 21.
    // Dimensions: method(2) × saltConc(3) × primerLen(3). Full grid 2×3×3 = 18.
    //
    // Note on the grid (cf. CRISPR-GUIDE-001 above): the checklist names the method
    // axis "basic/SantaLucia", but the implemented Tm models are (a) a salt-free
    // base Tm whose FORMULA is itself length-selected — Wallace's rule 2·(A+T)+4·(G+C)
    // for short oligos (<14 valid nt) and the Marmur-Doty GC% formula
    // 64.9 + 41·(#GC − 16.4)/N for ≥14 nt — and (b) that base Tm plus a
    // Schildkraut-Lifson salt correction 16.6·log10([Na⁺]). (A SantaLucia
    // nearest-neighbour model exists only as the ΔG-based 3′-stability metric, not a
    // Tm.) So method = {Basic, SaltCorrected}, and primerLen straddles the
    // Wallace↔Marmur-Doty switch (10 nt → Wallace; 14, 24 nt → Marmur-Doty).
    //
    // The combinatorial point: method and saltConc INTERACT. Under SaltCorrected the
    // salt axis shifts Tm by +16.6·log10([Na⁺]/1000) and is monotone increasing in
    // [Na⁺]; under Basic the salt axis is INERT (identical Tm at every saltConc).
    // primerLen interacts with method by selecting which base formula applies.
    // — Wallace 1979 NAR 6:3543; Marmur & Doty 1962 JMB 5:109; Schildkraut & Lifson
    //   1965 Biopolymers 3:195.
    // ═══════════════════════════════════════════════════════════════════════

    public enum TmMethod { Basic, SaltCorrected }

    /// <summary>Deterministic 50%-GC primer of length n ("ACGT…").</summary>
    private static string PrimerOfLen(int n) => string.Concat(Enumerable.Range(0, n).Select(i => "ACGT"[i % 4]));

    /// <summary>Independent re-derivation of the documented base Tm (no salt).</summary>
    private static double ExpectedBaseTm(string seq)
    {
        int at = seq.Count(c => c is 'A' or 'T'), gc = seq.Count(c => c is 'G' or 'C');
        int n = at + gc;
        if (n == 0) return 0;
        return n < 14 ? 2 * at + 4 * gc : Math.Max(0, 64.9 + 41.0 * (gc - 16.4) / n);
    }

    [Test, Combinatorial]
    public void PrimerTm_MatchesDocumentedFormula_SaltActiveOnlyForSaltCorrected(
        [Values(TmMethod.Basic, TmMethod.SaltCorrected)] TmMethod method,
        [Values(10.0, 50.0, 200.0)] double saltMm,
        [Values(10, 14, 24)] int primerLen)
    {
        string seq = PrimerOfLen(primerLen);
        double baseTm = ExpectedBaseTm(seq);

        double actual = method == TmMethod.Basic
            ? PrimerDesigner.CalculateMeltingTemperature(seq)
            : PrimerDesigner.CalculateMeltingTemperatureWithSalt(seq, saltMm);

        double expected = method == TmMethod.Basic
            ? baseTm                                                       // salt axis inert
            : Math.Round(baseTm + 16.6 * Math.Log10(saltMm / 1000.0), 1);  // Schildkraut-Lifson

        actual.Should().BeApproximately(expected, 1e-9,
            $"{method} Tm of a {primerLen}-mer at {saltMm} mM follows the documented formula");

        // Only the base Tm is clamped at 0; the additive salt correction may legitimately
        // drive a low-Tm short primer below 0 at very low [Na⁺] (e.g. 10-mer at 10 mM → −3.2 °C).
        if (method == TmMethod.Basic)
            actual.Should().BeGreaterThanOrEqualTo(0.0, "the base melting temperature is clamped non-negative");
    }

    /// <summary>
    /// Interaction witness: saltConc drives Tm only under SaltCorrected. The base
    /// model is constant across [Na⁺]; the salt-corrected model rises strictly with it.
    /// </summary>
    [Test]
    public void PrimerTm_SaltAxis_InertForBasic_MonotoneForSaltCorrected()
    {
        string seq = PrimerOfLen(24);

        double b10 = PrimerDesigner.CalculateMeltingTemperature(seq);
        double b200 = PrimerDesigner.CalculateMeltingTemperature(seq);
        b10.Should().Be(b200, "the base Tm has no salt term");

        double s10 = PrimerDesigner.CalculateMeltingTemperatureWithSalt(seq, 10.0);
        double s50 = PrimerDesigner.CalculateMeltingTemperatureWithSalt(seq, 50.0);
        double s200 = PrimerDesigner.CalculateMeltingTemperatureWithSalt(seq, 200.0);
        s10.Should().BeLessThan(s50);
        s50.Should().BeLessThan(s200, "Tm rises with [Na⁺] via +16.6·log10([Na⁺])");

        // At 1 M Na⁺ the correction vanishes, so salt-corrected ≡ base Tm.
        PrimerDesigner.CalculateMeltingTemperatureWithSalt(seq, 1000.0)
            .Should().BeApproximately(b10, 0.05);
    }

    /// <summary>
    /// Interaction witness: the primerLen axis selects the formula — a 12-mer uses
    /// Wallace (all-GC ⇒ 4·12 = 48 °C); a 20-mer with 10 GC uses Marmur-Doty
    /// (64.9 + 41·(10−16.4)/20 = 51.78 °C).
    /// </summary>
    [Test]
    public void PrimerTm_FormulaSwitch_WorkedExamples()
    {
        PrimerDesigner.CalculateMeltingTemperature("GCGCGCGCGCGC")  // 12 nt, Wallace
            .Should().BeApproximately(48.0, 1e-9);

        PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACACGTACGTAC")  // 20 nt, 10 GC, Marmur-Doty
            .Should().BeApproximately(51.78, 0.01);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-DESIGN-001 — Primer-candidate acceptance (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 22.
    // Dimensions: minLen(3) × maxLen(3) × gcRange(3) × tmRange(3). Grid 3⁴ = 81.
    //
    // Model (Primer3 / Untergasser 2012 acceptance semantics): a candidate primer
    // is accepted iff it lies inside EVERY configured window simultaneously —
    // length ∈ [MinLength, MaxLength], GC% ∈ [MinGcContent, MaxGcContent] and
    // Tm ∈ [MinTm, MaxTm] (plus structural filters held constant here). Acceptance
    // is therefore the logical CONJUNCTION of the per-axis membership tests, and a
    // violation of any single window must surface its own diagnostic.
    //
    // The combinatorial point: the four windows interact multiplicatively. The probe
    // primer (20 nt, 50% GC, Tm 51.78 °C) is placed so each axis straddles its three
    // windows — minLen {15,20,22} and maxLen {18,20,25} bracket length 20 on both
    // sides, the GC windows straddle 50%, the Tm windows straddle 51.78 °C — so every
    // axis genuinely flips acceptance and the AND is exercised across the grid.
    // ═══════════════════════════════════════════════════════════════════════

    // 20-mer, 50% GC, Tm 51.78 °C, hairpin-free, max homopolymer 2 (verified independently).
    private const string CleanPrimer20 = "GACGCTGTCTGAGACTAGAA";

    private static readonly (double Lo, double Hi)[] GcWindows = { (30, 45), (45, 55), (55, 70) };
    private static readonly (double Lo, double Hi)[] TmWindows = { (40, 50), (50, 60), (52, 62) };

    /// <summary>Permissive structural filters so only the length/GC/Tm windows can gate acceptance.</summary>
    private static PrimerParameters WithWindows(int minLen, int maxLen, double gcLo, double gcHi, double tmLo, double tmHi) =>
        PrimerDesigner.DefaultParameters with
        {
            MinLength = minLen,
            MaxLength = maxLen,
            MinGcContent = gcLo,
            MaxGcContent = gcHi,
            MinTm = tmLo,
            MaxTm = tmHi,
            MaxHomopolymer = 1000,
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false,
            Avoid3PrimeGC = false,
        };

    [Test, Combinatorial]
    public void PrimerDesign_AcceptanceIsConjunctionOfLengthGcTmWindows(
        [Values(15, 20, 22)] int minLen,
        [Values(18, 20, 25)] int maxLen,
        [Values(0, 1, 2)] int gcWin,
        [Values(0, 1, 2)] int tmWin)
    {
        var (gcLo, gcHi) = GcWindows[gcWin];
        var (tmLo, tmHi) = TmWindows[tmWin];

        var cand = PrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true,
            WithWindows(minLen, maxLen, gcLo, gcHi, tmLo, tmHi));

        double gc = PrimerDesigner.CalculateGcContent(CleanPrimer20);
        double tm = PrimerDesigner.CalculateMeltingTemperature(CleanPrimer20);
        int len = CleanPrimer20.Length;

        bool lenOk = len >= minLen && len <= maxLen;
        bool gcOk = gc >= gcLo && gc <= gcHi;
        bool tmOk = tm >= tmLo && tm <= tmHi;

        cand.IsValid.Should().Be(lenOk && gcOk && tmOk,
            $"accept ⟺ len∈[{minLen},{maxLen}] ∧ GC∈[{gcLo},{gcHi}] ∧ Tm∈[{tmLo},{tmHi}]");

        // Every violated window contributes exactly its own diagnostic; satisfied windows stay silent.
        cand.Issues.Any(i => i.StartsWith("Length")).Should().Be(!lenOk);
        cand.Issues.Any(i => i.StartsWith("GC content")).Should().Be(!gcOk);
        cand.Issues.Any(i => i.StartsWith("Tm")).Should().Be(!tmOk);
    }

    /// <summary>
    /// Interaction witness: from a config that accepts the primer, narrowing ANY
    /// single window past the primer's measured value flips acceptance to reject —
    /// confirming each axis is an independent necessary condition.
    /// </summary>
    [Test]
    public void PrimerDesign_EachWindow_IndependentlyGatesAcceptance()
    {
        var wide = WithWindows(1, 100, 0, 100, 0, 200);
        PrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide).IsValid
            .Should().BeTrue("a clean primer is accepted under an all-permissive config");

        double gc = PrimerDesigner.CalculateGcContent(CleanPrimer20);
        double tm = PrimerDesigner.CalculateMeltingTemperature(CleanPrimer20);

        PrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MinGcContent = gc + 5 })
            .IsValid.Should().BeFalse("a GC floor above the primer's GC excludes it");
        PrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MaxTm = tm - 1 })
            .IsValid.Should().BeFalse("a Tm ceiling below the primer's Tm excludes it");
        PrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MaxLength = CleanPrimer20.Length - 1 })
            .IsValid.Should().BeFalse("a length ceiling below the primer's length excludes it");
    }

    /// <summary>
    /// Worked example through the end-to-end design pipeline: any pair
    /// <see cref="PrimerDesigner.DesignPrimers"/> returns must honour every
    /// configured window on BOTH primers (acceptance is enforced during selection).
    /// </summary>
    [Test]
    public void PrimerDesign_DesignedPair_HonoursEveryWindow()
    {
        var template = new DnaSequence(DiverseDna(600));
        // Tm window chosen to match the Marmur-Doty scale for 18–25-mers.
        var param = PrimerDesigner.DefaultParameters with
        {
            MinLength = 18, MaxLength = 25, MinGcContent = 35, MaxGcContent = 65, MinTm = 45, MaxTm = 65,
        };

        var result = PrimerDesigner.DesignPrimers(template, 280, 320, param);

        foreach (var p in new[] { result.Forward, result.Reverse }.Where(p => p is not null))
        {
            p!.Length.Should().BeInRange(param.MinLength, param.MaxLength);
            p.GcContent.Should().BeInRange(param.MinGcContent, param.MaxGcContent);
            p.MeltingTemperature.Should().BeInRange(param.MinTm, param.MaxTm);
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-STRUCT-001 — Primer secondary-structure analysis (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 23.
    // Spec: tests/TestSpecs/PRIMER-STRUCT-001.md (hairpin / primer-dimer / 3′-end ΔG
    //       stability / homopolymer / dinucleotide battery).
    //
    // Axis mapping (deviation from the checklist, documented per spec): the implemented
    // structure routines — HasHairpinPotential, HasPrimerDimer, Calculate3PrimeStability
    // — are SALT- and TEMPERATURE-INDEPENDENT (the SantaLucia 1998 ΔG°37 is the fixed
    // reference state at 1 M NaCl/37 °C; salt enters the model only through Tm, which is
    // covered by PRIMER-TM-001). So the checklist's saltConc/tempC axes do not exist in
    // this code; the parameters that actually gate structure detection are the structural
    // stringency knobs. The grid is therefore built on the cleanest crisp interaction —
    // the primer-dimer decision boundary — over 3′-complementarity × minComplementarity ×
    // primerLen, with hairpin and 3′-stability covered as theory-anchored witnesses.
    //
    // Model (Wikipedia Primer-dimer; Primer3): two primers dimerize when their 3′ ends
    // are complementary. HasPrimerDimer compares the 3′ window (≤8 nt) of primer1 against
    // the 3′ window of primer2 and reports a dimer iff the complementary-base count meets
    // minComplementarity. So detection = (3′-complementary count ≥ minComplementarity),
    // and — crucially — depends ONLY on the 3′ window, not on total primer length.
    // ═══════════════════════════════════════════════════════════════════════

    // Each pair is engineered so the 8-base 3′ comparison window holds EXACTLY K complementary
    // bases: primer1's 3′ window is fixed (DimerP1Window); E2 is the first 8 bases that
    // reverse-complement(primer2) must present, complementary to primer1 in its 3′-most K positions.
    private const string DimerP1Window = "GACTGACT";
    private static readonly (int K, string E2)[] DimerWindows = { (2, "AAACAAGA"), (4, "AAACCTGA"), (6, "AAGACTGA") };

    private static (string P1, string P2) MakeDimerPair(string e2, int length)
    {
        string p1 = new string('T', length - 8) + DimerP1Window;          // 5′ filler outside the 3′ window
        string p2 = RevComp(e2 + new string('A', length - 8));            // revComp(p2) starts with E2
        return (p1, p2);
    }

    [Test, Combinatorial]
    public void PrimerStruct_DimerDetection_ThresholdOnThreePrimeComplementarity(
        [Values(0, 1, 2)] int windowIdx,
        [Values(3, 4, 5)] int minComplementarity,
        [Values(8, 12, 20)] int primerLen)
    {
        var (k, e2) = DimerWindows[windowIdx];
        var (p1, p2) = MakeDimerPair(e2, primerLen);

        PrimerDesigner.HasPrimerDimer(p1, p2, minComplementarity)
            .Should().Be(k >= minComplementarity,
                $"a {k}-base 3′-complementary window dimerizes iff {k} ≥ minComplementarity({minComplementarity})");
    }

    /// <summary>
    /// Interaction witness: each engineered pair has EXACTLY K complementary 3′ bases
    /// (detected at threshold K, rejected at K+1), and detection is invariant to the
    /// 5′ length — confirming primer-dimer is governed solely by the 3′ window.
    /// </summary>
    [Test]
    public void PrimerStruct_Dimer_HasExactCount_AndIgnoresPrimerLength()
    {
        foreach (var (k, e2) in DimerWindows)
        {
            foreach (int len in new[] { 8, 12, 20 })
            {
                var (p1, p2) = MakeDimerPair(e2, len);
                PrimerDesigner.HasPrimerDimer(p1, p2, k).Should().BeTrue($"the 3′ window has {k} complementary bases");
                PrimerDesigner.HasPrimerDimer(p1, p2, k + 1).Should().BeFalse($"the 3′ window has only {k} complementary bases");
            }
        }
    }

    /// <summary>
    /// Hairpin witness (stem-loop theory): a perfect stem-loop (4-bp stem, 3-nt loop) is
    /// detected; a non-self-complementary oligo is not; and the structural length guard
    /// (length ≥ 2·minStem + minLoop) suppresses an otherwise-real hairpin when the
    /// stringency exceeds what the length can hold. — spec invariant #3.
    /// </summary>
    [Test]
    public void PrimerStruct_Hairpin_DetectsStemLoop_AndRespectsLengthGuard()
    {
        const string hairpin = "GGGGAAACCCC"; // 4-bp stem (GGGG/CCCC) + 3-nt loop (AAA), length 11
        PrimerDesigner.HasHairpinPotential(hairpin).Should().BeTrue("a 4-bp stem with a 3-nt loop folds");
        PrimerDesigner.HasHairpinPotential("AAAAAAAAAAA").Should().BeFalse("a poly-A oligo is not self-complementary");

        // Requiring a 5-bp stem needs length ≥ 2·5+3 = 13 > 11, so the guard returns false.
        PrimerDesigner.HasHairpinPotential(hairpin, minStemLength: 5).Should().BeFalse("length 11 < 2·5+3");
    }

    /// <summary>
    /// 3′-stability witness (SantaLucia 1998 + Primer3): the most stable 5-mer GCGCG and
    /// the least stable TATAT match their published ΔG°37, and GC-rich 3′ ends are strictly
    /// more stable (more negative) than AT-rich ones. — spec M16/M17, invariant #4.
    /// </summary>
    [Test]
    public void PrimerStruct_ThreePrimeStability_MatchesSantaLucia()
    {
        PrimerDesigner.Calculate3PrimeStability("GCGCG").Should().BeApproximately(-6.86, 1e-9);
        PrimerDesigner.Calculate3PrimeStability("TATAT").Should().BeApproximately(-0.86, 1e-9);
        PrimerDesigner.Calculate3PrimeStability("GCGCG")
            .Should().BeLessThan(PrimerDesigner.Calculate3PrimeStability("TATAT"), "GC-rich 3′ ends are more stable");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-DESIGN-001 — Hybridization-probe design (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 24.
    // Dimensions: minLen(3) × maxLen(3) × tmRange(3) × gcRange(3). Grid 3⁴ = 81.
    //
    // Model (microarray/qPCR probe selection): ProbeDesigner.DesignProbes scans every
    // window in [MinLength, MaxLength], keeps those whose GC fraction lies within the
    // configured band (with a ±0.1 early-rejection slack), scores each (penalising GC,
    // Tm, homopolymer, self-complementarity, structure deviations) and returns the top
    // probes. Two of the four windows are HARD constraints on any returned probe —
    // length ∈ [MinLength, MaxLength] (loop bound) and GC ∈ [MinGc−0.1, MaxGc+0.1]
    // (early rejection) — whereas the Tm window is a SOFT scoring term (−0.3 penalty),
    // so out-of-Tm probes can still be returned. Tm itself is the salt-adjusted formula
    // 81.5 + 16.6·log10[Na⁺] + 41·GC − 600/N (≥14 nt) — Howley 1979 / SantaLucia.
    //
    // The combinatorial point: the four windows interact during the scan. The grid
    // asserts the HARD invariants hold for EVERY returned probe under EVERY (len×len×
    // gc×tm) configuration — i.e. the search never emits a probe outside the length or
    // GC bands regardless of how the other windows are set. The soft Tm axis and the GC
    // ±0.1 slack are pinned by dedicated witnesses below.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string ProbeTarget = DiverseDna(400);
    private static readonly (double Lo, double Hi)[] ProbeGcWindows = { (0.30, 0.50), (0.40, 0.60), (0.50, 0.70) };
    private static readonly (double Lo, double Hi)[] ProbeTmWindows = { (50, 65), (55, 70), (60, 80) };

    private static ProbeDesigner.ProbeParameters MkProbeParam(
        int minLen, int maxLen, double gcLo, double gcHi, double tmLo, double tmHi) =>
        new(MinLength: minLen, MaxLength: maxLen, MinTm: tmLo, MaxTm: tmHi,
            MinGc: gcLo, MaxGc: gcHi, MaxHomopolymer: 100,
            AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

    [Test, Combinatorial]
    public void ProbeDesign_ReturnedProbes_HonourLengthAndGcWindows(
        [Values(18, 20, 22)] int minLen,
        [Values(24, 28, 32)] int maxLen,
        [Values(0, 1, 2)] int gcWin,
        [Values(0, 1, 2)] int tmWin)
    {
        var (gcLo, gcHi) = ProbeGcWindows[gcWin];
        var (tmLo, tmHi) = ProbeTmWindows[tmWin];

        var probes = ProbeDesigner.DesignProbes(
            ProbeTarget, MkProbeParam(minLen, maxLen, gcLo, gcHi, tmLo, tmHi), maxProbes: 20).ToList();

        foreach (var p in probes)
        {
            p.Sequence.Length.Should().BeInRange(minLen, maxLen, "probe length is hard-bounded by the window");
            p.GcContent.Should().BeInRange(gcLo - 0.1 - 1e-9, gcHi + 0.1 + 1e-9,
                "GC honours the window within the ±0.1 early-rejection slack");
        }
    }

    /// <summary>
    /// Interaction witness: the Tm window is the soft scoring axis. With every other
    /// term held neutral, moving the Tm band from "includes all probe Tms" to "excludes
    /// them" lowers each common probe's score by exactly the 0.3 Tm penalty.
    /// </summary>
    [Test]
    public void ProbeDesign_TmWindow_ContributesExactScoringPenalty()
    {
        var tmIn = MkProbeParam(20, 20, 0.0, 1.0, 0, 200);    // GC + Tm both always in range
        var tmOut = tmIn with { MinTm = 200, MaxTm = 300 };    // only the Tm term now fires

        var inByStart = ProbeDesigner.DesignProbes(ProbeTarget, tmIn, maxProbes: 50).ToDictionary(p => p.Start);
        var outByStart = ProbeDesigner.DesignProbes(ProbeTarget, tmOut, maxProbes: 50).ToDictionary(p => p.Start);

        var common = inByStart.Keys.Intersect(outByStart.Keys).ToList();
        common.Should().NotBeEmpty("a uniform Tm penalty preserves ranking, so the same probes are returned");
        foreach (int start in common)
            (inByStart[start].Score - outByStart[start].Score).Should().BeApproximately(0.3, 1e-6,
                "an out-of-Tm probe loses exactly the 0.3 Tm penalty");
    }

    /// <summary>
    /// Interaction witness: the GC band carries a ±0.1 early-rejection slack. A 65%-GC
    /// 20-mer is accepted when 0.65 ≤ MaxGc+0.1 (band 0.40–0.60) but rejected once it
    /// exceeds the slack (band 0.40–0.50 ⇒ ceiling 0.60).
    /// </summary>
    [Test]
    public void ProbeDesign_GcWindow_AppliesTenPercentSlack()
    {
        const string gc65 = "GCGCGCGCGCGCGAATATAT"; // 20 nt, 13 GC = 0.65

        ProbeDesigner.DesignProbes(gc65, MkProbeParam(20, 20, 0.40, 0.60, 0, 200), maxProbes: 5)
            .Should().ContainSingle("0.65 is within MaxGc(0.60)+0.1 slack");
        ProbeDesigner.DesignProbes(gc65, MkProbeParam(20, 20, 0.40, 0.50, 0, 200), maxProbes: 5)
            .Should().BeEmpty("0.65 exceeds MaxGc(0.50)+0.1 slack");
    }

    /// <summary>
    /// Worked example through the qPCR-defaults pipeline: every returned probe respects
    /// the strict length window and the result set is capped at maxProbes.
    /// </summary>
    [Test]
    public void ProbeDesign_QpcrDefaults_RespectLengthWindowAndCap()
    {
        var param = ProbeDesigner.Defaults.qPCR; // 20–30 nt
        var probes = ProbeDesigner.DesignProbes(ProbeTarget, param, maxProbes: 10).ToList();

        probes.Count.Should().BeLessThanOrEqualTo(10);
        foreach (var p in probes)
            p.Sequence.Length.Should().BeInRange(param.MinLength, param.MaxLength);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-VALID-001 — Probe validation (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 25.
    // Spec: tests/TestSpecs/PROBE-VALID-001.md — canonical methods ValidateProbe &
    //       CheckSpecificity.
    //
    // Axis mapping (documented deviation): the canonical validators have NO gc/tm
    // windows — those belong to design (PROBE-DESIGN-001). Probe VALIDATION judges a
    // probe against a reference for (a) cross-hybridization / off-target specificity and
    // (b) self-complementarity / secondary structure. So the genuine validation axes are
    // off-target multiplicity, the probe's self-complementarity, and the self-comp
    // threshold (selfCompMax) — the grid uses offTargetCount × selfCompProbe ×
    // selfCompThreshold.
    //
    // Model (microarray cross-hybridization; off-target editing): specificity is a strict
    // function of exact hit count — 0 hits ⇒ 0.0 (probe doesn't bind), 1 ⇒ 1.0 (unique),
    // N ⇒ 1/N (invariants #4-6). A probe is reported invalid when it accumulates issues
    // (>1 off-target site, self-comp above threshold, secondary structure) UNLESS the
    // lenient clause holds (≤1 hit AND self-comp ≤ 0.4). The combinatorial point: the
    // IsValid decision composes the off-target and self-complementarity checks, and the
    // self-comp threshold interacts with the probe's self-comp to gate that issue.
    // ═══════════════════════════════════════════════════════════════════════

    // selfComp = fraction of positions Watson-Crick-paired with the mirror position; all
    // three probes are secondary-structure-free by construction (verified independently).
    private static readonly (string Seq, double SelfComp)[] ValidationProbes =
    {
        ("AAAAAAAAAAAAAAAAAAAA", 0.0),
        ("TGGCGCGGGGTAACGCGCGC", 0.5),
        ("ACGTACGTACGTACGTACGT", 1.0),
    };

    /// <summary>Reference holding exactly <paramref name="k"/> exact copies of the probe, C-padded, G-spaced.</summary>
    private static string BuildOffTargetReference(string probe, int k)
    {
        string core = k == 0 ? "" : string.Join(new string('G', 30), Enumerable.Repeat(probe, k));
        return new string('C', 40) + core + new string('C', 40);
    }

    [Test, Combinatorial]
    public void ProbeValid_SpecificityAndIssues_FollowValidationRules(
        [Values(0, 1, 3)] int offTargetCount,
        [Values(0, 1, 2)] int probeIdx,
        [Values(0.25, 0.40, 0.60)] double selfCompThreshold)
    {
        var (probe, expSelfComp) = ValidationProbes[probeIdx];
        string reference = BuildOffTargetReference(probe, offTargetCount);

        var v = ProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0,
            selfComplementarityThreshold: selfCompThreshold);

        // Specificity invariants (#4/#5/#6) — depend solely on off-target multiplicity.
        double expSpec = offTargetCount == 0 ? 0.0 : offTargetCount == 1 ? 1.0 : 1.0 / offTargetCount;
        v.OffTargetHits.Should().Be(offTargetCount);
        v.SpecificityScore.Should().BeApproximately(expSpec, 1e-9);

        // Self-complementarity measured exactly and bounded to [0,1] (#2).
        v.SelfComplementarity.Should().BeApproximately(expSelfComp, 1e-9);
        v.SelfComplementarity.Should().BeInRange(0.0, 1.0);
        v.HasSecondaryStructure.Should().BeFalse("the three probes are structure-free by construction");

        // IsValid composes the off-target and self-comp checks (with the lenient clause).
        bool offIssue = offTargetCount > 1;
        bool selfIssue = expSelfComp > selfCompThreshold;
        int issueCount = (offIssue ? 1 : 0) + (selfIssue ? 1 : 0);
        bool expectedValid = issueCount == 0 || (offTargetCount <= 1 && expSelfComp <= 0.4);

        v.IsValid.Should().Be(expectedValid);
        v.Issues.Any(i => i.Contains("off-target")).Should().Be(offIssue);
        v.Issues.Any(i => i.StartsWith("Self-complementarity")).Should().Be(selfIssue);
    }

    /// <summary>
    /// Interaction witness: at exact matching the suffix-tree CheckSpecificity score
    /// agrees with ValidateProbe's specificity for every probe × hit-count — the two
    /// specificity routes are consistent.
    /// </summary>
    [Test]
    public void ProbeValid_CheckSpecificity_AgreesWithValidateProbe()
    {
        foreach (var (probe, _) in ValidationProbes)
            foreach (int k in new[] { 0, 1, 3 })
            {
                string reference = BuildOffTargetReference(probe, k);
                double viaTree = ProbeDesigner.CheckSpecificity(probe, global::SuffixTree.SuffixTree.Build(reference));
                double viaValidate = ProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0).SpecificityScore;
                viaTree.Should().BeApproximately(viaValidate, 1e-9, $"both score {k} exact hits identically");
            }
    }

    /// <summary>
    /// Interaction witness: the mismatch budget governs off-target sensitivity — a site
    /// carrying one substitution is missed under strict matching but found once one
    /// mismatch is permitted. — spec S4.
    /// </summary>
    [Test]
    public void ProbeValid_MaxMismatches_GovernsOffTargetSensitivity()
    {
        string probe = ValidationProbes[1].Seq;
        char[] mutated = probe.ToCharArray();
        mutated[10] = mutated[10] == 'A' ? 'C' : 'A';                 // exactly one substitution
        string reference = new string('C', 30) + new string(mutated) + new string('C', 30);

        ProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0).OffTargetHits
            .Should().Be(0, "strict matching misses a 1-mismatch site");
        ProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 1).OffTargetHits
            .Should().Be(1, "allowing one mismatch finds it");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RESTR-DIGEST-001 — Restriction digest simulation (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 27.
    // Spec: tests/TestSpecs/RESTR-DIGEST-001.md (canonical RestrictionAnalyzer.Digest).
    // Dimensions: enzyme(4) × topology(2: linear/circular) × fragments(3). Grid 4×2×3 = 24.
    //
    // Model (restriction digest; gel electrophoresis): cutting a molecule at k distinct
    // forward-strand sites yields k+1 fragments when LINEAR (two free ends) but exactly
    // k fragments when CIRCULAR (a plasmid; the origin-spanning piece joins the first and
    // last cut), with the special case that an uncut circle is one full-length fragment
    // (Addgene Plasmids 101). In every case fragment lengths SUM to the sequence length.
    //
    // The combinatorial point: topology and cut count interact — the same enzyme on the
    // same sequence yields a different fragment count purely from topology, while the
    // length-conservation invariant holds across all enzyme × topology × cut-count cells.
    // The "fragments" axis is realised as the number of engineered cut sites (0/1/2).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Sequence with exactly <paramref name="n"/> copies of a (palindromic) site in C-filler.</summary>
    private static string BuildDigestSequence(string recognition, int n)
    {
        string filler = new string('C', 10);
        var sb = new System.Text.StringBuilder(filler);
        for (int i = 0; i < n; i++) sb.Append(recognition).Append(filler);
        return sb.ToString();
    }

    [Test, Combinatorial]
    public void RestrDigest_FragmentCount_FollowsTopology_AndConservesLength(
        [Values("EcoRI", "BamHI", "HindIII", "EcoRV")] string enzyme,
        [Values(MoleculeTopology.Linear, MoleculeTopology.Circular)] MoleculeTopology topology,
        [Values(0, 1, 2)] int cutCount)
    {
        string recognition = RestrictionAnalyzer.GetEnzyme(enzyme)!.RecognitionSequence;
        var dna = new DnaSequence(BuildDigestSequence(recognition, cutCount));

        // Construction guard: exactly cutCount distinct forward-strand cut sites.
        RestrictionAnalyzer.FindSites(dna, enzyme).Count(s => s.IsForwardStrand)
            .Should().Be(cutCount, "the C-filler must not add or hide sites");

        var fragments = RestrictionAnalyzer.Digest(dna, topology, enzyme).ToList();

        int expected = topology == MoleculeTopology.Linear ? cutCount + 1 : (cutCount == 0 ? 1 : cutCount);
        fragments.Count.Should().Be(expected,
            $"{topology} digest with {cutCount} cuts yields {expected} fragments");

        fragments.Sum(f => f.Length).Should().Be(dna.Length, "fragment lengths conserve total length");
        fragments.Should().OnlyContain(f => f.Length > 0, "every fragment is non-empty");

        if (topology == MoleculeTopology.Linear)
        {
            fragments[0].LeftEnzyme.Should().BeNull("the 5′ terminus has no upstream cut");
            fragments[^1].RightEnzyme.Should().BeNull("the 3′ terminus has no downstream cut");
        }
    }

    /// <summary>
    /// Interaction witness: for k ≥ 1 cuts a linear molecule yields exactly one more
    /// fragment than the circular molecule of the same sequence (the linear ends fuse
    /// into the origin-spanning fragment when circular).
    /// </summary>
    [Test]
    public void RestrDigest_Linear_YieldsExactlyOneMoreFragmentThanCircular()
    {
        foreach (string enzyme in new[] { "EcoRI", "BamHI", "HindIII", "EcoRV" })
            foreach (int k in new[] { 1, 2 })
            {
                string recognition = RestrictionAnalyzer.GetEnzyme(enzyme)!.RecognitionSequence;
                var dna = new DnaSequence(BuildDigestSequence(recognition, k));

                int linear = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, enzyme).Count();
                int circular = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, enzyme).Count();
                (linear - circular).Should().Be(1, $"{enzyme}, {k} cuts: linear = circular + 1");
            }
    }

    /// <summary>
    /// Worked example: a single EcoRI site splits a linear molecule into two fragments
    /// whose sequences are the exact substrings around the cut, with the cut enzyme named
    /// on the inner ends; the digest summary's sizes are sorted descending and conserve length.
    /// </summary>
    [Test]
    public void RestrDigest_SingleEcoRI_WorkedExample()
    {
        var dna = new DnaSequence(BuildDigestSequence("GAATTC", 1)); // CCCCCCCCCC GAATTC CCCCCCCCCC
        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, "EcoRI").ToList();

        fragments.Should().HaveCount(2);
        string.Concat(fragments.Select(f => f.Sequence)).Should().Be(dna.Sequence, "fragments reassemble the template");
        fragments[0].RightEnzyme.Should().Be("EcoRI");
        fragments[1].LeftEnzyme.Should().Be("EcoRI");

        var summary = RestrictionAnalyzer.GetDigestSummary(dna, "EcoRI");
        summary.FragmentSizes.Should().BeInDescendingOrder("gel-ordering convention");
        summary.FragmentSizes.Sum().Should().Be(dna.Length);
        summary.LargestFragment.Should().BeGreaterThanOrEqualTo(summary.SmallestFragment);
        summary.EnzymesUsed.Should().Contain("EcoRI");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RESTR-FILTER-001 — Restriction-enzyme filtering (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 224.
    // Spec: tests/TestSpecs/RESTR-FILTER-001.md (canonical GetEnzymesByCutLength / GetBluntCutters /
    //       GetStickyCutters). ADVANCED §10.
    // Dimensions: enzyme(3) × criteria(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Type II restriction enzymes; Wikipedia sticky/blunt ends): recognition sites are 4–8 bp;
    // an enzyme cuts blunt (both strands at the same offset) or sticky (a 5'/3' overhang). Blunt and
    // sticky cutters partition the catalogue.
    //
    // Axis mapping (documented): criteria → the filter {ByCutLength range, Blunt, Sticky}; enzyme → the
    // recognition-length band {4, 6, 8} used by the range filter. The combinatorial point: each filter
    // returns exactly the enzymes satisfying its predicate, and blunt ⊎ sticky = the full set.
    // ═══════════════════════════════════════════════════════════════════════

    public enum RestrFilter { ByCutLength, Blunt, Sticky }

    [Test, Combinatorial]
    public void RestrFilter_ReturnsMatchingEnzymes_AcrossLengthBandAndCriteria(
        [Values(4, 6, 8)] int lengthBand,
        [Values(RestrFilter.ByCutLength, RestrFilter.Blunt, RestrFilter.Sticky)] RestrFilter criteria)
    {
        switch (criteria)
        {
            case RestrFilter.ByCutLength:
                RestrictionAnalyzer.GetEnzymesByCutLength(lengthBand, lengthBand)
                    .Should().OnlyContain(e => e.RecognitionSequence.Length == lengthBand,
                        "the range filter returns only enzymes whose recognition length is in [min,max]");
                // The single-length overload agrees with the degenerate range.
                RestrictionAnalyzer.GetEnzymesByCutLength(lengthBand).Select(e => e.Name)
                    .Should().BeEquivalentTo(RestrictionAnalyzer.GetEnzymesByCutLength(lengthBand, lengthBand).Select(e => e.Name));
                break;
            case RestrFilter.Blunt:
                RestrictionAnalyzer.GetBluntCutters().Should().OnlyContain(e => e.IsBluntEnd, "blunt cutters cut both strands at one offset");
                break;
            default:
                RestrictionAnalyzer.GetStickyCutters().Should().OnlyContain(e => !e.IsBluntEnd, "sticky cutters leave an overhang");
                break;
        }
    }

    /// <summary>
    /// Interaction witness — blunt and sticky cutters partition the catalogue (disjoint, and together
    /// the whole set), and the range filter is monotone in its bounds.
    /// </summary>
    [Test]
    public void RestrFilter_BluntStickyPartition_AndRangeMonotone()
    {
        var blunt = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToHashSet();
        var sticky = RestrictionAnalyzer.GetStickyCutters().Select(e => e.Name).ToHashSet();
        blunt.Should().NotIntersectWith(sticky, "an enzyme is either blunt or sticky, not both");

        var all = RestrictionAnalyzer.GetEnzymesByCutLength(1, 100).Select(e => e.Name).ToHashSet();
        blunt.Union(sticky).Should().BeEquivalentTo(all, "blunt ⊎ sticky = the full catalogue");

        int narrow = RestrictionAnalyzer.GetEnzymesByCutLength(6, 6).Count();
        int wide = RestrictionAnalyzer.GetEnzymesByCutLength(4, 8).Count();
        wide.Should().BeGreaterThanOrEqualTo(narrow, "a wider length range admits no fewer enzymes");
    }
}
