namespace Seqeron.Genomics.Alignment;

/// <summary>The sequence helpers of the MUMmer programs (<c>tigrinc.cc</c>): the IUPAC-aware lower-case complement and the reverse complement of a 1-based sequence array (<c>Reverse_Complement (S, 1, Len)</c>), shared by <c>postnuc</c> and the <c>show-*</c> ports.</summary>
internal static class MummerSequence
{
    /// <summary><c>Complement</c> of <c>tigrinc.cc</c>: A/C/G/T and the IUPAC codes (R-Y, S, W, M-K, B-V, D-H); every other character is lower-cased and kept.</summary>
    public static byte Complement(byte ch) => char.ToLowerInvariant((char)ch) switch
    {
        'a' => (byte)'t', 'c' => (byte)'g', 'g' => (byte)'c', 't' => (byte)'a', 'r' => (byte)'y', 'y' => (byte)'r', 's' => (byte)'s', 'w' => (byte)'w',
        'm' => (byte)'k', 'k' => (byte)'m', 'b' => (byte)'v', 'd' => (byte)'h', 'h' => (byte)'d', 'v' => (byte)'b', var c => (byte)c,
    };

    /// <summary>The reverse complement of a 1-based sequence array (element 0 and the element after the last base are padding); the input is not changed.</summary>
    public static byte[] ReverseComplement(byte[] oneBased)
    {
        long len = oneBased.Length - 2;
        var rev = (byte[])oneBased.Clone();
        long lo = 1, hi = len;
        while (lo <= hi)
        {
            byte ch = rev[hi];
            rev[hi] = Complement(rev[lo]);
            rev[lo] = Complement(ch);
            lo++;
            hi--;
        }

        return rev;
    }
}
