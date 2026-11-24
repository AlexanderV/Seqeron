using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The sections of <c>a^n</c> as readers that cost a few instructions per read (work package S2.4): <c>SA[r] = n - r</c>, the packed child cell <c>Down[r] = r + 1</c> (0 at the last rank), <c>LCP[r] = r - 1</c> (0 for the first two ranks) and the
/// text of one letter. The long descent of <c>a^k</c> for k near n takes k steps of six reads each; the general readers (<see cref="PosGen{TWidth}"/>) compute the closed form through a virtual call and need minutes for 2^33 steps, these a
/// fraction of that. They are checked against the closed form (<see cref="WideEngineLargeTests"/>) and are always checked readers (a rank outside 0..n is an <see cref="InvalidDataException"/>).
/// </summary>
internal static class PowerReaders
{
    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void Outside(long index, long n) => throw new InvalidDataException($"a section of {n + 1} elements is read at element {index}: the file is damaged");

    internal readonly struct Sa(long n) : IPosReader
    {
        public static int ElementBytes => 5;

        public long Count => n + 1;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ulong Get(long index)
        {
            if ((ulong)index > (ulong)n) Outside(index, n);
            return (ulong)(n - index);
        }
    }

    internal readonly struct Child(long n) : IPosReader
    {
        public static int ElementBytes => 5;

        public long Count => n + 1;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ulong Get(long index)
        {
            if ((ulong)index > (ulong)n) Outside(index, n);
            return index < n ? (ulong)(index + 1) : 0;
        }
    }

    internal readonly struct Lcp(long n) : IPosReader
    {
        public static int ElementBytes => 8;

        public long Count => n + 1;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ulong Get(long index)
        {
            if ((ulong)index > (ulong)n) Outside(index, n);
            return index <= 1 ? 0 : (ulong)(index - 1);
        }
    }

    internal readonly struct Text(long n) : IWideText
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Symbol(long position)
        {
            if ((ulong)position >= (ulong)n) throw new ArgumentOutOfRangeException(nameof(position), position, $"the text has {n} symbols");
            return 0;
        }
    }
}
