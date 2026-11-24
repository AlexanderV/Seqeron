using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The checked bridges from the 64-bit views of work package S2.3 to the 32-bit engine that stays until S2.4 and S2.5 (<see cref="MappedSearch"/>, <c>EsaIndex.*</c>, <see cref="LcpIntervalScan"/>, <c>EsaNavigator</c>, <c>EsaNavigationPack</c>): where the
/// engine needs an <see cref="int"/> from a <see cref="long"/>, the conversion is one of these two and never an unchecked cast. A length or a count that the 32-bit engine cannot hold is a <see cref="NotSupportedException"/> that names the
/// member (the facade says the same for the contract: use the long member), a value that was read from a file and does not fit is an <see cref="InvalidDataException"/> (the file is damaged for an engine whose text fits the integers). Every use
/// of this class is a place that S2.4 / S2.5 removes; the list is in the log of PROGRESS-S2.md.
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
