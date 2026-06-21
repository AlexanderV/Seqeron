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
}
