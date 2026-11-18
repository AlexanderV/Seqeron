using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>Distribution tuple of B1 (base mode): the position <c>SA[i]</c>, its predecessor <c>SA[i-1]</c> (the Phi value) and the rank <c>i</c>; 12 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal readonly struct PhiTuple(uint position, uint phi, uint rank)
{
    internal readonly uint Position = position;
    internal readonly uint Phi = phi;
    internal readonly uint Rank = rank;
}

/// <summary>Distribution tuple of B2: a rank and the LCP value at that rank (<c>PLCP[SA[i]]</c>); 8 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal readonly struct RankTuple(uint rank, uint value)
{
    internal readonly uint Rank = rank;
    internal readonly uint Value = value;
}
