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
}
