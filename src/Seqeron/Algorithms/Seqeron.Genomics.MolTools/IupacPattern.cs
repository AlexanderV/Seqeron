namespace Seqeron.Genomics.MolTools;

/// <summary>
/// The one per-position rule by which a site pattern (a restriction recognition site, a PAM) is matched against a
/// sequence in this assembly, composed from the canonical <see cref="IupacHelper"/> predicates (no second ambiguity
/// table here). Used by the restriction scans, the methylation mask and the PAM scan, which each carried a private copy.
/// </summary>
/// <remarks>
/// <para>
/// The rule is the one of both references. Biopython 1.88 <c>Bio.Restriction</c> compiles a site's <c>N</c> to the regex
/// <c>.</c> and every other code to the class of its A/C/G/T expansion (<c>SfiI.compsite</c> is <c>GGCC.....GGCC</c>,
/// <c>HincII.compsite</c> is <c>GT[CT][AG]AC</c>), and the sequence keeps its own IUPAC letters; CRISPOR's
/// <c>patMatch</c> (crispor.py L2347-2392) says <c>if patChar=="N": continue</c> and otherwise needs a concrete base of the
/// code's set. So a pattern <c>N</c> accepts <b>any</b> sequence code - an <c>N</c> included (<c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c>)
/// - while a pattern A/C/G/T or degenerate code accepts only a concrete base it expands to: an ambiguity code in the sequence
/// still blocks it (<c>EcoRI.search(Seq("AGAANTCA")) == []</c>, <c>HincII.search(Seq("AGTYRACA")) == []</c>; a sequence <c>N</c>
/// at a PAM's <c>G</c> or <c>R</c> is no PAM). A sequence character outside the 15 IUPAC nucleotide codes (and any lower-case
/// letter: callers fold to upper case first) matches nothing, not even a pattern <c>N</c> - Biopython would let a site <c>N</c>
/// match any letter, <c>X</c> included, a deliberate narrowing.
/// </para>
/// <para>
/// <see cref="IupacHelper.MatchesIupac"/> alone would let a sequence <c>N</c> block a pattern's <c>N</c> position
/// (it answers for A/C/G/T only), which is why the pattern <c>N</c> is decided here.
/// </para>
/// </remarks>
internal static class IupacPattern
{
    /// <summary>Whether <paramref name="sequenceBase"/> satisfies the pattern code <paramref name="patternCode"/>.</summary>
    public static bool MatchesCode(char sequenceBase, char patternCode) =>
        patternCode == 'N'
            ? IupacHelper.IsNucleotideCode(sequenceBase)
            : IupacHelper.MatchesIupac(sequenceBase, patternCode);

    /// <summary>
    /// Whether <paramref name="pattern"/> matches <paramref name="seq"/> at <paramref name="position"/> without wrapping;
    /// <see langword="false"/> when the window does not fit.
    /// </summary>
    public static bool MatchesAt(string seq, int position, string pattern)
    {
        if (position < 0 || position + pattern.Length > seq.Length)
            return false;

        for (int i = 0; i < pattern.Length; i++)
        {
            if (!MatchesCode(seq[position + i], pattern[i]))
                return false;
        }

        return true;
    }

    /// <summary>
    /// Whether <paramref name="pattern"/> matches the window of <c>pattern.Length</c> bases starting at
    /// <paramref name="position"/> of a circular sequence (indices wrap modulo <c>seq.Length</c>).
    /// </summary>
    public static bool MatchesAtCircular(string seq, int position, string pattern)
    {
        int n = seq.Length;
        for (int k = 0; k < pattern.Length; k++)
        {
            if (!MatchesCode(seq[(position + k) % n], pattern[k]))
                return false;
        }

        return true;
    }
}
