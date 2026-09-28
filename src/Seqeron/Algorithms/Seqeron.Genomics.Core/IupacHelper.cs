namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Helper methods for IUPAC nucleotide code matching.
    /// Implements the IUPAC-IUB 1970/NC-IUB 1984 nucleotide ambiguity code standard.
    /// Only the 15 standard codes (A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V) are accepted.
    /// </summary>
    public static class IupacHelper
    {
        /// <summary>
        /// Determines whether <paramref name="code"/> is one of the 15 IUPAC DNA nucleotide codes
        /// A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V (upper case only) — exactly the set accepted by
        /// <see cref="MatchesIupac"/> and Biopython <c>IUPACData.ambiguous_dna_letters</c>.
        /// </summary>
        public static bool IsNucleotideCode(char code) => code is
            'A' or 'C' or 'G' or 'T' or 'N' or 'R' or 'Y' or 'S' or 'W' or 'K' or 'M' or 'B' or 'D' or 'H' or 'V';

        /// <summary>
        /// Determines if a nucleotide matches an IUPAC ambiguity code.
        /// </summary>
        /// <param name="nucleotide">The nucleotide to check (A, C, G, T).</param>
        /// <param name="iupacCode">The IUPAC code (A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V).</param>
        /// <returns>True if the nucleotide matches the IUPAC code.</returns>
        public static bool MatchesIupac(char nucleotide, char iupacCode) => iupacCode switch
        {
            'A' => nucleotide == 'A',
            'C' => nucleotide == 'C',
            'G' => nucleotide == 'G',
            'T' => nucleotide == 'T',
            'N' => nucleotide is 'A' or 'C' or 'G' or 'T',
            'R' => nucleotide is 'A' or 'G',          // puRine
            'Y' => nucleotide is 'C' or 'T',          // pYrimidine
            'S' => nucleotide is 'G' or 'C',          // Strong
            'W' => nucleotide is 'A' or 'T',          // Weak
            'K' => nucleotide is 'G' or 'T',          // Keto
            'M' => nucleotide is 'A' or 'C',          // aMino
            'B' => nucleotide is 'C' or 'G' or 'T',   // not A
            'D' => nucleotide is 'A' or 'G' or 'T',   // not C
            'H' => nucleotide is 'A' or 'C' or 'T',   // not G
            'V' => nucleotide is 'A' or 'C' or 'G',   // not T
            _ => throw new ArgumentOutOfRangeException(
                nameof(iupacCode), iupacCode,
                "Not a valid IUPAC nucleotide code. Valid codes: A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V.")
        };
    }
}
