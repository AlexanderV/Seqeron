namespace SuffixTree.Genomics
{
    /// <summary>
    /// Helper methods for IUPAC nucleotide code matching.
    /// </summary>
    public static class IupacHelper
    {
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
            _ => nucleotide == iupacCode
        };
    }
}
