namespace SuffixTree.Esa.Tests.Large;

/// <summary>STUB of work package S2.1 (tests first): the closed form of <c>(a^k b)*</c> is written after its tests.</summary>
internal sealed class PeriodicAkB : ClosedFormText
{
    internal PeriodicAkB(long n, long k) : base(n) => throw new NotImplementedException("S2.1 stub");

    /// <summary><c>(ab)^m</c>: the case k = 1 with n = 2m.</summary>
    internal static PeriodicAkB AbRepeat(long m) => new(2 * m, 1);

    internal override string Name => throw new NotImplementedException();
    internal override long Period => throw new NotImplementedException();
    internal override int Sigma => throw new NotImplementedException();
    internal override long SymbolCount(int code) => throw new NotImplementedException();
    internal override int InfiniteCodeAt(long position) => throw new NotImplementedException();
    internal override long SuffixAt(long rank) => throw new NotImplementedException();
    internal override long RankOf(long position) => throw new NotImplementedException();
    internal override long Lcp(long rank) => throw new NotImplementedException();
    internal override Int128 LcpSum => throw new NotImplementedException();
    internal override long InternalCount => throw new NotImplementedException();
    internal override long LrsDepth => throw new NotImplementedException();
    internal override long LrsRank => throw new NotImplementedException();
    internal override RankRange RangeOf(ReadOnlySpan<char> pattern) => throw new NotImplementedException();
    internal override RankRange RangeOfRun(long length) => throw new NotImplementedException();
    internal override RankRange RangeOfFactor(long start, long length) => throw new NotImplementedException();
    internal override long Down(long i) => throw new NotImplementedException();
    internal override long Up(long j) => throw new NotImplementedException();
    internal override long Next(long k) => throw new NotImplementedException();
}
