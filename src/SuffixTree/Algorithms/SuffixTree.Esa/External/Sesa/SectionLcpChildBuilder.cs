namespace SuffixTree.Esa.External;

/// <summary>
/// PH4 (base mode) and PH5 (layout A) as <see cref="LcpChildBuilder"/> runs them, with the outputs as sections of one container file instead of three files of their own (work package 1b.7e): the suffix array is read from a
/// byte range of a file (<see cref="ExtLcpOptions.SaOffset"/>), the CHILD section is a range of the container (known up front: <c>4 (n + 1)</c> bytes), and the LCP and its overflow list are placed by a callback once
/// the width is known (after B2, before the first LCP entry is written), so that the container can grow behind the sections that were known after PH0. Nothing is copied. The stage creates and removes its own scratch
/// subdirectories <c>ph4</c> and <c>ph5</c>; the container file is the owner's (it is neither completed nor deleted here, also not on a failure).
/// </summary>
internal static class SectionLcpChildBuilder
{
    /// <summary>
    /// Builds the LCP, its overflow list and the child table of <paramref name="text"/> from the suffix array at byte <paramref name="saOffset"/> of <paramref name="saPath"/>. The child table is written to
    /// <c>[childOffset, childOffset + 4 (n + 1))</c> of <paramref name="file"/>, which must exist at that length when the stage starts; <paramref name="place"/> returns where the LCP and the overflow list go.
    /// </summary>
    /// <exception cref="InvalidDataException">The suffix array is not a permutation of the text positions (wrong range, <c>SA[0] != n</c>, an entry out of range or repeated).</exception>
    internal static LcpChildResult Build(string saPath, long saOffset, MappedTextReader text, MappedOutputFile file, long childOffset, Func<LcpPlacementRequest, LcpPlacement> place, LcpChildOptions options)
        => throw new NotImplementedException();
}
