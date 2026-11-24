using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The checked bridges from the 64-bit views and the 64-bit search engine of stage S2 to the places that still read 32-bit numbers: where such a place needs an <see cref="int"/> from a <see cref="long"/>, the conversion is one of these two and
/// never an unchecked cast. A length or a count that the 32-bit place cannot hold is a <see cref="NotSupportedException"/> that names the member (and the 64-bit one to use), a value that was read from a file and does not fit is an
/// <see cref="InvalidDataException"/> (the file is damaged for an engine whose text fits the integers). Work package S2.4 removed the bridges of the search side (the descent, <c>EsaStatistics</c>, <c>LeftmostOccurrence</c>, the constructor of
/// <see cref="MappedSearch"/>: they are 64-bit now); what stays is the 32-bit surface that the public members keep until S2.6 (<c>EsaIndex.Length</c>, <c>CountOccurrences</c>, <c>Traverse(ISuffixTreeVisitor)</c>, the facade) and the navigation
/// side that S2.5 replaces (<c>EsaNavigator</c>, <c>EsaNavigationPack</c>, the 32-bit members of <c>SaView</c>, <c>MappedChild</c> and <c>LcpPyramid</c>); the list is in the log of PROGRESS-S2.md.
/// </summary>
internal static class IntBridge
{
    /// <summary>A length or a count (<paramref name="length"/>) as an <see cref="int"/>.</summary>
    /// <exception cref="NotSupportedException">The length is above <see cref="int.MaxValue"/> (or negative: no length is).</exception>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int Length(long length, string member)
    {
        if ((ulong)length > int.MaxValue) ThrowLength(length, member);
        return (int)length;
    }

    /// <summary>A rank, a position or a value that came from a container as an <see cref="int"/>.</summary>
    /// <exception cref="InvalidDataException">The value does not fit: the file is damaged.</exception>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int Value(long value, string what)
    {
        if ((ulong)value > int.MaxValue) ThrowValue(value, what);
        return (int)value;
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void ThrowLength(long length, string member) =>
        throw new NotSupportedException($"{member} is {length}, more than the 32-bit engine holds ({int.MaxValue}): use the 64-bit member (stage S2 replaces the 32-bit engine)");

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void ThrowValue(long value, string what) =>
        throw new InvalidDataException($"{what} is {value}, which is not a 32-bit number: the file is damaged");
}
