namespace SuffixTree.Esa.External;

/// <summary>
/// Where the sections go, as far as they are known after PH0 (ADR-004 section 4.3): every size except the LCP and its overflow list is a function of n, the width b and sigma, so TEXT, ALPHABET, COUNTS, SA, CHILD, PSI and
/// BUILD_INFO are placed up front, in this order after the header page, each at the alignment of ADR-003 section 3; LCP follows them, because its width is known only after PH4-B2. <see cref="ReservedLength"/> is the
/// length at which the file is created: everything up to the end of the LCP at its narrowest width (n + 1 bytes), which no text can undercut.
/// </summary>
internal sealed record SesaPlan
{
    internal long N { get; init; }

    internal int Bits { get; init; }

    internal int Sigma { get; init; }

    internal bool Psi { get; init; }

    internal long TextOffset { get; init; }

    internal long TextLength { get; init; }

    internal long AlphabetOffset { get; init; }

    internal long AlphabetLength { get; init; }

    internal long CountsOffset { get; init; }

    internal long CountsLength { get; init; }

    internal long SaOffset { get; init; }

    internal long SaLength { get; init; }

    internal long ChildOffset { get; init; }

    internal long ChildLength { get; init; }

    /// <summary>0 when the build has no PSI section.</summary>
    internal long PsiOffset { get; init; }

    internal long PsiLength { get; init; }

    internal long BuildInfoOffset { get; init; }

    /// <summary>Bytes reserved for the BUILD_INFO text (its length is known only at the end).</summary>
    internal long BuildInfoCapacity { get; init; }

    internal long LcpOffset { get; init; }

    /// <summary>The length at which the container file is created.</summary>
    internal long ReservedLength { get; init; }
}

/// <summary>The part of the layout that the width of the LCP decides.</summary>
internal sealed record SesaTail
{
    internal LcpKind Kind { get; init; }

    internal int LcpWidth { get; init; }

    internal long LcpLength { get; init; }

    internal long OverflowCount { get; init; }

    /// <summary>0 when there is no overflow list.</summary>
    internal long OverflowOffset { get; init; }

    internal long OverflowLength { get; init; }

    /// <summary>The length of the file: the end of the last section.</summary>
    internal long FinalLength { get; init; }
}

/// <summary>The arithmetic of <see cref="SesaPlan"/> and <see cref="SesaTail"/>.</summary>
internal static class SesaLayout
{
    /// <summary>Bytes reserved for BUILD_INFO.</summary>
    internal const int BuildInfoCapacity = 4096;

    /// <summary>The plan for a text described by <paramref name="info"/>, with or without the PSI section.</summary>
    internal static SesaPlan Plan(PackedTextInfo info, bool psi) => throw new NotImplementedException();

    /// <summary>The LCP and overflow sections once the width and the number of escapes are known.</summary>
    internal static SesaTail PlaceTail(SesaPlan plan, LcpKind kind, long overflowCount) => throw new NotImplementedException();
}
