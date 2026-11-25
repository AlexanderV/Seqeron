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

            // PamSequence is always read on the strand it was found on, so it satisfies the IUPAC motif directly.
            string pamOnStrand = site.PamSequence;
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
    // 64.9 + 41·(#GC − 16.4)/N for ≥14 nt (both assume 50 mM Na⁺) — and (b) the
    // OligoCalc salt-adjusted Tm (Kibbe 2007): Wallace + 16.6·log10([Na⁺]/0.050 M) for
    // <14 nt, 100.5 + 41·GC/N − 820/N + 16.6·log10([Na⁺] M) for ≥14 nt. (A SantaLucia
    // nearest-neighbour model exists only as the ΔG-based 3′-stability metric, not a
    // Tm.) So method = {Basic, SaltCorrected}, and primerLen straddles the
    // Wallace↔Marmur-Doty switch (10 nt → Wallace; 14, 24 nt → Marmur-Doty).
    //
    // The combinatorial point: method and saltConc INTERACT. Under SaltCorrected the
    // salt axis shifts Tm by +16.6 °C per decade of [Na⁺] and is monotone increasing in
    // [Na⁺]; under Basic the salt axis is INERT (identical Tm at every saltConc).
    // primerLen interacts with method by selecting which base formula applies.
    // — Wallace 1979 NAR 6:3543; Marmur & Doty 1962 JMB 5:109; Schildkraut & Lifson
    //   1965 Biopolymers 3:195.
    // ═══════════════════════════════════════════════════════════════════════

    public enum TmMethod { Basic, SaltCorrected }

    /// <summary>Deterministic 50%-GC primer of length n ("ACGT…").</summary>
    private static string PrimerOfLen(int n) => string.Concat(Enumerable.Range(0, n).Select(i => "ACGT"[i % 4]));

    /// <summary>Independent re-derivation of the OligoCalc salt-adjusted Tm.</summary>
    private static double ExpectedSaltTm(string seq, double saltMm)
    {
        int at = seq.Count(c => c is 'A' or 'T'), gc = seq.Count(c => c is 'G' or 'C');
        int n = at + gc;
        double m = saltMm / 1000.0;
        return n < 14
            ? 2 * at + 4 * gc + 16.6 * Math.Log10(m / 0.050)
            : 100.5 + 41.0 * gc / n - 820.0 / n + 16.6 * Math.Log10(m);
    }

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
            : Math.Round(ExpectedSaltTm(seq, saltMm), 1);                  // OligoCalc salt adjusted

        actual.Should().BeApproximately(expected, 1e-9,
            $"{method} Tm of a {primerLen}-mer at {saltMm} mM follows the documented formula");

        // Only the base Tm is clamped at 0 (the salt-adjusted formula is reported as computed).
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

        // Short oligo: the Wallace rule is defined at 50 mM, so salt-adjusted(50 mM) ≡ base Tm.
        string shortSeq = PrimerOfLen(10);
        PrimerDesigner.CalculateMeltingTemperatureWithSalt(shortSeq, 50.0)
            .Should().BeApproximately(PrimerDesigner.CalculateMeltingTemperature(shortSeq), 1e-9);
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
    // primer (20 nt, 50% GC, Primer3-default Tm 56.43 °C) is placed so each axis straddles its three
    // windows — minLen {15,20,22} and maxLen {18,20,25} bracket length 20 on both
    // sides, the GC windows straddle 50%, the Tm windows straddle 56.43 °C — so every
    // axis genuinely flips acceptance and the AND is exercised across the grid.
    // ═══════════════════════════════════════════════════════════════════════

    // 20-mer, 50% GC, Primer3-default Tm 56.4298 °C (primer3-py 2.3.1 calc_tm), hairpin-free,
    // max homopolymer 2 (verified independently).
    private const string CleanPrimer20 = "GACGCTGTCTGAGACTAGAA";

    private static readonly (double Lo, double Hi)[] GcWindows = { (30, 45), (45, 55), (55, 70) };
    private static readonly (double Lo, double Hi)[] TmWindows = { (45, 55), (55, 65), (57, 67) };

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
    // Model (Primer3 alignment-mode PRIMER_PAIR_COMPL_END, libprimer3.cc characterize_pair +
    // dpal.c DPAL_GLOBAL_END): two primers dimerize when their 3′ ends are complementary, i.e. the
    // 3′-terminal K bases of primer2 are the reverse complement of the 3′-terminal K bases of
    // primer1. HasPrimerDimer reports a dimer iff that 3′-anchored complementarity score meets
    // minComplementarity, independently of the 5′ length.
    // ═══════════════════════════════════════════════════════════════════════

    // primer1 = A-filler + DimerP1Window; primer2 = A-filler + revcomp(last K bases of the window).
    // The A fillers cannot pair with each other, so the compl_end score is exactly K (verified
    // against Primer3's dpal.c compiled from source: 200/400/600 for K = 2/4/6 at lengths 8/12/20).
    private const string DimerP1Window = "GACTGACT";
    private static readonly int[] DimerWindows = { 2, 4, 6 };

    private static (string P1, string P2) MakeDimerPair(int k, int length)
    {
        string p1 = new string('A', length - 8) + DimerP1Window;
        string p2 = new string('A', length - k) + RevComp(DimerP1Window[^k..]);
        return (p1, p2);
    }

    [Test, Combinatorial]
    public void PrimerStruct_DimerDetection_ThresholdOnThreePrimeComplementarity(
        [Values(0, 1, 2)] int windowIdx,
        [Values(3, 4, 5)] int minComplementarity,
        [Values(8, 12, 20)] int primerLen)
    {
        int k = DimerWindows[windowIdx];
        var (p1, p2) = MakeDimerPair(k, primerLen);

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
        foreach (int k in DimerWindows)
        {
            foreach (int len in new[] { 8, 12, 20 })
            {
                var (p1, p2) = MakeDimerPair(k, len);
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
    // Model (microarray/qPCR probe selection): ThermodynamicProbeDesigner.DesignProbes scans every
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
    // IsValid decision composes the off-target and self-complementarity checks. Since audit
    // round 2 (B07 A6) the fallback self-complementarity issue is Primer3's alignment-mode
    // self_any/self_end > 12.00; the legacy selfCompThreshold axis must have no effect.
    // ═══════════════════════════════════════════════════════════════════════

    // selfComp = fraction of positions Watson-Crick-paired with the mirror position.
    // SelfAny = Primer3 alignment-mode self_any (dpal.c + align(), compiled): 0.00, 10.00, 20.00 (self_end equal).
    // Hairpin = the fallback stem-loop screen PrimerDesigner.HasHairpinPotential (exact ≥ 4-bp stem, loop ≥ 3 nt;
    // audit round 5, A5-2): poly-A has no stem; GCGC (pos 2) · GCGC (pos 16) and ACGT (pos 0) · ACGT (pos 8) are exact
    // 4-bp stems with ≥ 3-nt loops. (The former private scan tested a loop of exactly 3 nt only and missed both.)
    private static readonly (string Seq, double SelfComp, double SelfAny, bool Hairpin)[] ValidationProbes =
    {
        ("AAAAAAAAAAAAAAAAAAAA", 0.0, 0.0, false),
        ("TGGCGCGGGGTAACGCGCGC", 0.5, 10.0, true),
        ("ACGTACGTACGTACGTACGT", 1.0, 20.0, true),
    };

    /// <summary>Reference holding exactly <paramref name="k"/> exact copies of the probe, C-padded, G-spaced.</summary>
    private static string BuildOffTargetReference(string probe, int k)
    {
        string core = k == 0 ? "" : string.Join(new string('G', 30), Enumerable.Repeat(probe, k));
        return new string('C', 40) + core + new string('C', 40);
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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-NNTM-001 — Nearest-neighbour salt / mismatch / dangling-end Tm (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 240.
    // Spec: tests/TestSpecs/PRIMER-NNTM-001.md (PrimerDesigner.CalculateMeltingTemperatureNNMismatch). ADVANCED §10.
    //
    // Sources: SantaLucia (1998) unified NN; Owczarzy et al. (2004) monovalent / (2008) divalent salt
    //   corrections; Allawi & SantaLucia (1997) internal mismatches.
    //
    // Model: Tm = ΔH/(ΔS + R·ln(C_T/x)) − 273.15, then a salt correction. A duplex is scored through
    //   the mismatch-capable NN path against either the perfect WC complement or a single central
    //   mismatch. The four axes carry signed, documented effects.
    //
    // Dimensions: len(3) × [Na⁺](2) × [Mg²⁺](2) × mismatch(2). Grid 3×2×2×2 = 24 (exhaustive).
    //
    // The combinatorial point: length, monovalent salt, divalent salt and a mismatch each move Tm in a
    // known direction, and the effects are independent. Every cell isolates its own axis against the
    // one-step-lower setting: longer duplex ⇒ higher Tm; higher [Na⁺] ⇒ higher Tm; added [Mg²⁺] ⇒
    // higher Tm; a central mismatch ⇒ lower Tm.
    // ═══════════════════════════════════════════════════════════════════════

    private const string NnTemplate = "ATGCGTACATGCGTACATGCGTAC"; // 24 nt, ~50% GC
    private static char NnComp(char b) => b switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => b };
    private static string PerfectBottom(string top) => new(top.Select(NnComp).ToArray());

    private static string MismatchBottom(string top)
    {
        var arr = top.Select(NnComp).ToArray();
        int mid = arr.Length / 2;
        arr[mid] = arr[mid] == 'A' ? 'C' : 'A'; // break complementarity at the centre → internal mismatch
        return new string(arr);
    }

    private static double NnTm(int len, double na, double mg, bool mismatch)
    {
        string top = NnTemplate[..len];
        string bottom = mismatch ? MismatchBottom(top) : PerfectBottom(top);
        return PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            top, bottom, sodiumMolar: na, magnesiumMolar: mg,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent);
    }

    [Test, Combinatorial]
    public void NnTm_AxisEffects_AcrossLengthSaltMagnesiumAndMismatch(
        [Values(12, 18, 24)] int len,
        [Values(0.05, 0.2)] double sodium,
        [Values(0.0, 0.003)] double magnesium,
        [Values(false, true)] bool mismatch)
    {
        double tm = NnTm(len, sodium, magnesium, mismatch);
        double.IsNaN(tm).Should().BeFalse("a computable NN duplex yields a finite Tm");
        tm.Should().BeInRange(-40.0, 160.0, "Tm is in a physical range");

        if (sodium > 0.05)
            tm.Should().BeGreaterThan(NnTm(len, 0.05, magnesium, mismatch), "higher [Na⁺] raises Tm (Owczarzy 2004)");
        if (magnesium > 0.0)
            tm.Should().BeGreaterThan(NnTm(len, sodium, 0.0, mismatch), "added [Mg²⁺] raises Tm (Owczarzy 2008)");
        if (mismatch)
            tm.Should().BeLessThan(NnTm(len, sodium, magnesium, false), "a central mismatch destabilises the duplex (Allawi & SantaLucia)");
        if (len > 12)
            tm.Should().BeGreaterThan(NnTm(12, sodium, magnesium, mismatch), "a longer duplex stacks more ⇒ higher Tm");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-HAIRPIN-001 — Unimolecular hairpin folder / secondary-structure Tm (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 241.
    // Spec: tests/TestSpecs/PRIMER-HAIRPIN-001.md
    //       (PrimerDesigner.FindMostStableHairpin / CalculateHairpinThermodynamicsNtthal). ADVANCED §10.
    //
    // Sources: SantaLucia & Hicks (2004) NN hairpin model (stem stacks + loop ΔG); primer3 ntthal.
    //
    // Model: the most stable hairpin's ΔG37 is the NN sum over the stem stacks plus a loop ΔG penalty.
    //   A planted G^stem · T^loop · C^stem hairpin folds a stem of exactly `stem` base pairs around a
    //   `loop`-nt loop. More stem stacks lower ΔG37; a longer loop raises (destabilises) ΔG37.
    //
    // Dimensions: stemLen(3) × loopLen(3) × specialLoop(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // Axis mapping (documented): specialLoop is realised as the FOLDER axis — the plain NN folder
    // (FindMostStableHairpin, no special-loop bonus) vs the primer3 ntthal model
    // (CalculateHairpinThermodynamicsNtthal, which bundles the special tri/tetraloop ΔG bonuses).
    // Both fold the same planted stem; the bonus itself is pinned by a dedicated witness.
    // ═══════════════════════════════════════════════════════════════════════

    public enum HairpinEngine { NnFolder, Ntthal }

    private static string Hairpin(int stem, int loop) => new string('G', stem) + new string('T', loop) + new string('C', stem);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-DIMER-001 — ntthal self / hetero-dimer Tm (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 242.
    // Spec: tests/TestSpecs/PRIMER-DIMER-001.md (PrimerDesigner.FindMostStableDimer). ADVANCED §10.
    //
    // Sources: SantaLucia (1998) NN; primer3 ntthal thermodynamic alignment.
    //
    // Model: the most stable duplex between two strands is the contiguous WC alignment whose ΔH/ΔS is
    //   the NN sum over its consecutive base-pair steps. Fully complementary strands form a dimer whose
    //   base-pair count equals the length of the shorter complementary span.
    //
    // Dimensions: lenA(3) × lenB(3) × self/hetero(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // Axis mapping (documented): in SELF mode the single strand is a G·C palindrome of length lenA
    // (lenB inert); in HETERO mode strand A is poly-G(lenA) and strand B is poly-C(lenB), so the
    // optimal contiguous G·C dimer is exactly min(lenA, lenB). Every cell forms a stable dimer (ΔG37 < 0).
    // ═══════════════════════════════════════════════════════════════════════

    public enum DimerKind { Self, Hetero }

    private static string DimerPalindrome(int len) => string.Concat(Enumerable.Range(0, len).Select(i => "GC"[i % 2])); // self-complementary

    [Test, Combinatorial]
    public void Dimer_FormsStableDuplex_AcrossLengthsAndSelfHetero(
        [Values(6, 8, 10)] int lenA,
        [Values(6, 8, 10)] int lenB,
        [Values(DimerKind.Self, DimerKind.Hetero)] DimerKind kind)
    {
        if (kind == DimerKind.Self)
        {
            string s = DimerPalindrome(lenA); // lenB inert (documented)
            var d = PrimerDesigner.FindMostStableDimer(s, s);
            d.Should().NotBeNull("a self-complementary strand forms a self-dimer");
            d!.Value.BasePairs.Should().Be(lenA, "a G·C palindrome self-pairs over its full length");
            d.Value.DeltaG37.Should().BeLessThan(0.0, "a fully paired G·C self-dimer is stable");
        }
        else
        {
            string a = new('G', lenA);
            string b = new('C', lenB);                        // poly-G · poly-C → contiguous G·C duplex
            var d = PrimerDesigner.FindMostStableDimer(a, b);
            d.Should().NotBeNull("complementary strands form a hetero-dimer");
            d!.Value.BasePairs.Should().Be(Math.Min(lenA, lenB), "the dimer spans the shorter complementary region");
            d.Value.DeltaG37.Should().BeLessThan(0.0, "a contiguous WC hetero-dimer is stable");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-LNATM-001 — LNA-adjusted NN Tm + MGB probe design (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 243.
    // Spec: tests/TestSpecs/PROBE-LNATM-001.md
    //       (PrimerDesigner.CalculateMeltingTemperatureNNLna / ProbeDesigner.EvaluateMgbProbeDesign). ADVANCED §10.
    //
    // Sources: McTigue et al. (2004) LNA NN increments; Kutyavin et al. (2000) MGB probe design rules.
    //
    // Model: an LNA-modified duplex adds per-position stabilising increments to the DNA NN sum; with no
    //   LNA positions the Tm reduces EXACTLY to the plain DNA NN Tm. MGB probes are designed shorter
    //   (12–20mer, 3'-attached).
    //
    // Dimensions: len(3) × lnaCount(3) × mgbFlag(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // The combinatorial point: probe length and LNA count both raise Tm; the zero-LNA case is the plain
    // DNA Tm identity; and (when mgbFlag is set) the MGB length-window rule reflects the probe length.
    // ═══════════════════════════════════════════════════════════════════════

    private const string LnaSeq = "GCATGCATGCATGCATGCATGC";          // 22 nt
    private static readonly int[] LnaPositions = { 4, 5, 8, 9 };       // interior G/C positions (stabilising LNA)

    [Test, Combinatorial]
    public void LnaTm_LengthAndLnaRaiseTm_MgbWindowReflectsLength(
        [Values(14, 18, 22)] int len,
        [Values(0, 2, 4)] int lnaCount,
        [Values(false, true)] bool mgbFlag)
    {
        string probe = LnaSeq[..len];
        int[] positions = LnaPositions.Take(lnaCount).ToArray();

        double tm = PrimerDesigner.CalculateMeltingTemperatureNNLna(probe, positions);
        double tm0 = PrimerDesigner.CalculateMeltingTemperatureNNLna(probe, System.Array.Empty<int>());

        if (lnaCount == 0)
            tm.Should().BeApproximately(ThermoConstants.CalculateNearestNeighborTm(probe, parameterSet: NnParameterSet.AllawiSantaLucia1997,
                dnac1: 500, dnac2: 0, selfComplementary: true, sodium: 50, saltCorrection: NnSaltCorrection.Owczarzy2004, gasConstant: 1.9872), 1e-9,
                "zero LNA reduces to the SantaLucia (1998) unified DNA NN Tm (the base set of the LNA models); the GCAT-repeat probes are self-complementary");
        else
            tm.Should().BeGreaterThan(tm0, "stabilising LNA increments raise Tm above the unmodified probe");

        if (len > 14)
            PrimerDesigner.CalculateMeltingTemperatureNNLna(LnaSeq[..14], positions)
                .Should().BeLessThan(tm, "a longer probe (same LNA set) has a higher Tm");

        if (mgbFlag)
        {
            var mgb = ProbeDesigner.EvaluateMgbProbeDesign(probe);
            mgb.MgbAttachmentEnd.Should().Be("3'", "MGB is attached at the 3' end (Kutyavin 2000)");
            mgb.LengthInMgbRange.Should().Be(len is >= 12 and <= 20, "the MGB 12–20mer window reflects the probe length");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-EVALUE-001 — Karlin–Altschul off-target E-value / bit score (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 244.
    // Spec: tests/TestSpecs/PROBE-EVALUE-001.md (ProbeDesigner.ComputeKarlinAltschul). ADVANCED §10.
    //
    // Sources: Karlin & Altschul (1990) PNAS 87:2264; Altschul et al. (1990) BLAST.
    //
    // Model: bit score S' = (λS − ln K)/ln 2 and expectation E = m·n·2^(−S') = K·m·n·e^(−λS).
    //
    // Dimensions: rawScore(3) × dbSize(3) × K(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // The combinatorial point: every cell reproduces the closed-form bit/E identity from the returned λ,
    // and the three axes move E in known directions — a higher raw score lowers E, a larger database
    // raises E, and a larger K lowers the bit score (hence raises E).
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void KarlinAltschul_BitAndEValue_AcrossRawScoreDbSizeAndK(
        [Values(30.0, 45.0, 60.0)] double rawScore,
        [Values(1_000_000L, 100_000_000L, 3_000_000_000L)] long dbSize,
        [Values(0.711, 0.3)] double k)
    {
        const int queryLength = 25;
        var s = ProbeDesigner.ComputeKarlinAltschul(rawScore, queryLength, dbSize, k: k);

        double expectedBit = (s.Lambda * s.RawScore - Math.Log(k)) / Math.Log(2.0);
        s.BitScore.Should().BeApproximately(expectedBit, 1e-9, "bit = (λS − ln K)/ln 2");
        double expectedE = (double)queryLength * dbSize * Math.Pow(2.0, -s.BitScore);
        s.EValue.Should().BeApproximately(expectedE, Math.Abs(expectedE) * 1e-9 + 1e-300, "E = m·n·2^(−bit)");
        s.EValue.Should().BeGreaterThanOrEqualTo(0.0, "an E-value is non-negative");

        // Axis effects (same query length / K, isolating one axis).
        if (rawScore > 30.0)
            s.EValue.Should().BeLessThan(ProbeDesigner.ComputeKarlinAltschul(30.0, queryLength, dbSize, k: k).EValue,
                "a higher raw score is more significant ⇒ lower E");
        if (dbSize > 1_000_000L)
            s.EValue.Should().BeGreaterThan(ProbeDesigner.ComputeKarlinAltschul(rawScore, queryLength, 1_000_000L, k: k).EValue,
                "a larger search space ⇒ higher E");
        if (k < 0.711)
            s.BitScore.Should().BeGreaterThan(ProbeDesigner.ComputeKarlinAltschul(rawScore, queryLength, dbSize, k: 0.711).BitScore,
                "a smaller K (larger −ln K) raises the bit score");
    }
}
