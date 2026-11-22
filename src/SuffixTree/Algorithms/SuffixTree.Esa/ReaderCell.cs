#pragma warning disable S2325 // tests-first stub: the members are instance members of the implementation that follows
using System.Runtime.InteropServices;

namespace SuffixTree.Esa;

/// <summary>
/// The thread-exclusive slot of the readers of the file-backed trees (work package 1b.9g). Tests-first stub: the members exist so that the tests of the registry compile and fail; the implementation follows in the next commit.
/// </summary>
[StructLayout(LayoutKind.Explicit, Size = 8)]
internal sealed class ReaderCell
{
    internal static ReaderCell Current => throw new NotImplementedException();

    internal static long NextId() => throw new NotImplementedException();

    internal static int CellCount => throw new NotImplementedException();

    internal int Depth => throw new NotImplementedException();

    internal void Publish(long id) => throw new NotImplementedException();

    internal void Retract() => throw new NotImplementedException();

    internal static bool AnyInside(long id) => throw new NotImplementedException();
}
