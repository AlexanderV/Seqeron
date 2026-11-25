namespace Seqeron.Genomics.MolTools;

/// <summary>Self-complementarity (reverse-complement palindrome) test of an oligo, shared by the nearest-neighbour Tm paths.</summary>
public static class OligoSymmetry
{
    /// <summary>
    /// True when the oligo is its own reverse complement: even length and every mirrored pair is a Watson–Crick pair
    /// (A/T, C/G). Non-ACGT characters never pair with an A/T/C/G (two non-ACGT characters pass); the test is case-insensitive.
    /// The symmetry correction of the DNA duplex initiation (SantaLucia 1998, ΔS° −1.4 cal/(K·mol)) applies when this holds.
    /// </summary>
    public static bool IsSymmetric(string oligo)
    {
        ArgumentNullException.ThrowIfNull(oligo);
        int len = oligo.Length;
        if (len % 2 != 0)
            return false;
        for (int i = 0; i < len / 2; i++)
        {
            char x = char.ToUpperInvariant(oligo[i]);
            char y = char.ToUpperInvariant(oligo[len - 1 - i]);
            if (!MirrorPairOk(x, y))
                return false;
        }

        return true;
    }

    /// <summary>The string reversed (3′→5′ reading of the second strand).</summary>
    public static string Reverse(string s)
    {
        ArgumentNullException.ThrowIfNull(s);
        var c = s.ToCharArray();
        Array.Reverse(c);
        return new string(c);
    }

    private static bool MirrorPairOk(char x, char y) => x switch
    {
        'A' => y == 'T',
        'T' => y == 'A',
        'C' => y == 'G',
        'G' => y == 'C',
        _ => y is not ('A' or 'T' or 'C' or 'G'),
    };
}
