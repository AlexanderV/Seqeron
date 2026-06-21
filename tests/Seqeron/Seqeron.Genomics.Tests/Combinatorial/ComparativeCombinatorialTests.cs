namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Comparative Genomics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Comparative")]
public class ComparativeCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (ANIb) (Comparative)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 131.
    // Spec: tests/TestSpecs/COMPGEN-ANI-001.md (ComparativeGenomics.CalculateANI).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Goris et al. (2007) ANIb — the query is cut into consecutive non-overlapping fragments of
    // fragmentLength nt; each is ungapped-aligned to the reference for its best match; fragments with
    // identity > minIdentity (0.30) and alignable fraction ≥ minAlignableFraction (0.70) contribute; ANI is
    // the mean of qualifying fragment identities.
    //
    // Checklist axes kmerSize(3) × genomeLen(3) × divergence(3) map onto the real knobs:
    //   • kmerSize   → fragmentLength ∈ {5, 10, 20}.
    //   • genomeLen  → sequence length ∈ {40, 80, 120}.
    //   • divergence → the fraction of query bases mutated from the reference ∈ {0.0, 0.1, 0.3}.
    // Grid = 3³ = 27 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: ANI is a JOINT function of fragment size, genome length and divergence. Each
    // cell re-derives the ANIb procedure independently from the constructed sequences and checks production
    // matches, plus the structural bounds (ANI ∈ [0,1], identical → 1.0).
    // ═══════════════════════════════════════════════════════════════════════

    private const string Bases = "ACGT";

    /// <summary>
    /// For every (fragment length, genome length, divergence) the production ANI equals the ANIb procedure
    /// re-derived independently from the constructed query/reference, and lies in [0, 1].
    /// </summary>
    [Test, Combinatorial]
    public void CalculateANI_FragmentLengthGenomeDivergenceGrid_MatchesAnibProcedure(
        [Values(5, 10, 20)] int fragmentLength,
        [Values(40, 80, 120)] int genomeLen,
        [Values(0.0, 0.1, 0.3)] double divergence)
    {
        string reference = BuildVariedDna(genomeLen);
        string query = MutateAtDensity(reference, divergence);

        double expected = GroundTruthAni(query, reference, fragmentLength, 0.30, 0.70);

        double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength, 0.30, 0.70);

        ani.Should().BeApproximately(expected, 1e-9, "ANI = mean qualifying-fragment identity (Goris ANIb)");
        ani.Should().BeInRange(0.0, 1.0, "[INV-1] ANI ∈ [0, 1]");
        if (divergence == 0.0)
            ani.Should().BeApproximately(1.0, 1e-9, "[INV-2] identical sequences → ANI = 1.0");
    }

    /// <summary>
    /// Interaction witness (divergence axis monotonicity): at fixed fragment and genome length, higher
    /// query–reference divergence gives a lower (or equal) ANI. Source: Goris et al. (2007).
    /// </summary>
    [Test]
    public void CalculateANI_DivergenceAxis_IsMonotoneDecreasing()
    {
        string reference = BuildVariedDna(80);

        double identical = ComparativeGenomics.CalculateANI(reference, reference, 10);
        double low = ComparativeGenomics.CalculateANI(MutateAtDensity(reference, 0.1), reference, 10);
        double high = ComparativeGenomics.CalculateANI(MutateAtDensity(reference, 0.3), reference, 10);

        identical.Should().Be(1.0);
        identical.Should().BeGreaterThanOrEqualTo(low);
        low.Should().BeGreaterThanOrEqualTo(high);
    }

    /// <summary>
    /// Interaction witness (fragmentation, INV-4): the query is cut into floor(N / fragmentLength)
    /// consecutive non-overlapping fragments; a trailing partial fragment is ignored. With an identical
    /// reference every fragment is perfect so ANI = 1.0 regardless of the trailing remainder. Source:
    /// Goris et al. (2007) "consecutive fragments".
    /// </summary>
    [Test]
    public void CalculateANI_TrailingPartialFragment_IsIgnored()
    {
        string reference = BuildVariedDna(45); // 45 = 4×10 + 5 (trailing partial)
        ComparativeGenomics.CalculateANI(reference, reference, 10).Should().Be(1.0, "4 full fragments, partial ignored");
    }

    /// <summary>Witness (INV-5): empty input yields ANI 0.</summary>
    [Test]
    public void CalculateANI_EmptyInput_IsZero()
    {
        ComparativeGenomics.CalculateANI("", BuildVariedDna(40), 10).Should().Be(0.0);
        ComparativeGenomics.CalculateANI(BuildVariedDna(40), "", 10).Should().Be(0.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Helpers — engineered constructs + independent ANIb ground truth
    // ───────────────────────────────────────────────────────────────────────

    /// <summary>Builds a low-repeat DNA sequence so that the best ungapped placement is at the aligned offset.</summary>
    private static string BuildVariedDna(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = Bases[(i * 7 + (i / 4) * 3 + 1) % 4];
        return new string(chars);
    }

    /// <summary>Mutates every (1/divergence)-th base to a different base; divergence 0 returns the input unchanged.</summary>
    private static string MutateAtDensity(string sequence, double divergence)
    {
        if (divergence <= 0.0) return sequence;
        int stride = (int)Math.Round(1.0 / divergence);
        var chars = sequence.ToCharArray();
        for (int i = 0; i < chars.Length; i += stride)
            chars[i] = Bases[(Array.IndexOf(Bases.ToCharArray(), chars[i]) + 1) % 4]; // flip to next base
        return new string(chars);
    }

    /// <summary>Independent ANIb ground truth (Goris et al. 2007): mean identity of qualifying fragments.</summary>
    private static double GroundTruthAni(string query, string reference, int fragLen, double minId, double minFrac)
    {
        if (string.IsNullOrEmpty(query) || string.IsNullOrEmpty(reference)) return 0.0;
        string q = query.ToUpperInvariant();
        string r = reference.ToUpperInvariant();

        var qualifying = new List<double>();
        for (int start = 0; start + fragLen <= q.Length; start += fragLen)
        {
            (double id, double frac) = BestUngappedMatch(q.Substring(start, fragLen), r);
            if (id > minId && frac >= minFrac)
                qualifying.Add(id);
        }
        return qualifying.Count > 0 ? qualifying.Average() : 0.0;
    }

    private static (double identity, double alignableFraction) BestUngappedMatch(string fragment, string reference)
    {
        int fragLen = fragment.Length;
        if (reference.Length < fragLen) return (0.0, 0.0);

        int best = 0;
        for (int offset = 0; offset + fragLen <= reference.Length; offset++)
        {
            int matches = 0;
            for (int k = 0; k < fragLen; k++)
                if (fragment[k] == reference[offset + k]) matches++;
            if (matches > best) best = matches;
        }
        return (Math.Min((double)best / fragLen, 1.0), 1.0);
    }
}
