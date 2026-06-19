using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Comparative-genomics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (ANIb, Goris et al. 2007).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 131.
///
/// API under test (ComparativeGenomics.CalculateANI):
///   The query genome is cut into consecutive fragments of fragmentLength; each fragment is
///   ungapped-best-placed in the reference; ANI is the mean identity of the qualifying fragments,
///   reported as a fraction in [0, 1] (1.0 ≡ 100 %).
///
/// Relations (derived from the ANIb definition, NOT from output):
///   • INV  (self-identity is maximal): every fragment of a genome is a perfect substring of
///          itself, so ANI(A,A)=1.0 regardless of fragmentation.
///   • MON  (more mutations ⇒ lower ANI): each substituted base lowers the matching count of its
///          fragment, so adding substitutions monotonically lowers ANI.
///   • SYM  (reciprocal symmetry under equal-length, full-length alignment): ANIb is asymmetric in
///          general because only the query is fragmented (Goris 2007 / pyani). For equal-length
///          genomes aligned as a single full-length fragment, the only placement is the diagonal
///          (offset 0), so both directions count the same matched positions and ANI(A,B)=ANI(B,A).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ComparativeMetamorphicTests
{
    // A fixed pseudo-genome over {A,C,G,T}; its own .Length is used as the full-length fragment size.
    private const string Genome = "ACGTTGCAACGTGGATCCGTACGATCGATTACAGGCATTAGCATCGTA";

    // Substitutes the first <paramref name="count"/> bases with a guaranteed-different base,
    // producing a copy that differs from the original at exactly <paramref name="count"/> positions.
    private static string Substitute(string seq, int count)
    {
        char[] arr = seq.ToCharArray();
        for (int i = 0; i < count; i++)
            arr[i] = arr[i] == 'A' ? 'C' : 'A';
        return new string(arr);
    }

    #region COMPGEN-ANI-001 INV — identical genomes give the maximal ANI of 1.0

    [Test]
    [Description("INV: every fragment of a genome is a perfect substring of itself, so ANI(A,A)=1.0 (the [0,1] form of 100 %), here exercised over multiple fragments.")]
    public void Ani_IdenticalGenomes_IsOne()
    {
        // fragmentLength 12 over the length-48 genome → 4 fragments, each a perfect self-match.
        ComparativeGenomics.CalculateANI(Genome, Genome, fragmentLength: 12).Should().Be(1.0,
            because: "each consecutive fragment occurs verbatim in the genome, so every per-fragment identity is 1.0 and so is their mean");
    }

    #endregion

    #region COMPGEN-ANI-001 MON — more substitutions lower ANI

    [Test]
    [Description("MON: each substituted base reduces the matched-position count, so introducing progressively more substitutions monotonically lowers ANI.")]
    public void Ani_MoreSubstitutions_LowerAni()
    {
        double previous = double.MaxValue;
        foreach (int mutations in new[] { 0, 1, 3, 6 })
        {
            // Full-length single fragment: ANI = (L − mutations)/L exactly, so the trend is strict.
            double ani = ComparativeGenomics.CalculateANI(Genome, Substitute(Genome, mutations), fragmentLength: Genome.Length);
            ani.Should().BeLessThan(previous, because: $"{mutations} substitutions leave fewer matched bases than {previous} did");
            previous = ani;
        }
    }

    #endregion

    #region COMPGEN-ANI-001 SYM — reciprocal symmetry for equal-length full-length alignment

    [Test]
    [Description("SYM: for equal-length genomes aligned as a single full-length fragment the only placement is the diagonal, so the two reciprocal directions count the same matched positions and ANI(A,B)=ANI(B,A).")]
    public void Ani_EqualLengthFullFragment_Symmetric()
    {
        foreach (int mutations in new[] { 2, 4, 8 })
        {
            string a = Genome;
            string b = Substitute(Genome, mutations);

            ComparativeGenomics.CalculateANI(a, b, fragmentLength: Genome.Length)
                .Should().Be(ComparativeGenomics.CalculateANI(b, a, fragmentLength: Genome.Length),
                    because: "matched-position count is symmetric and, at equal length, only the offset-0 placement exists");
        }
    }

    #endregion
}
