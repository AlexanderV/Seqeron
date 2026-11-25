using System.Collections.Frozen;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Authoritative Turner 2004 nearest-neighbor free-energy parameters (kcal/mol at 37°C),
    /// shared by every algorithm that folds RNA. Centralising the tables here guarantees that
    /// the RNA secondary-structure folder (<c>Seqeron.Genomics.Analysis.RnaSecondaryStructure</c>)
    /// and the pre-miRNA hairpin energy model (<c>Seqeron.Genomics.Annotation.MiRnaAnalyzer</c>)
    /// always score with the identical parameter set and can never silently drift apart.
    /// </summary>
    /// <remarks>
    /// Source: Nearest Neighbor Database (NNDB), Turner 2004 set —
    /// https://rna.urmc.rochester.edu/NNDB/turner04/
    /// </remarks>
    public static class Turner2004Parameters
    {
        /// <summary>
        /// Watson-Crick and GU-wobble nearest-neighbor stacking energies (kcal/mol at 37°C).
        /// Key format: <c>5'-XY-3' / 3'-X'Y'-5'</c> where X pairs with X' and Y with Y'.
        /// Source: NNDB — wc-parameters.html and gu-parameters.html.
        /// </summary>
        public static readonly IReadOnlyDictionary<string, double> StackingEnergies = new Dictionary<string, double>
        {
            // Watson-Crick stacking (10 unique + 6 symmetric = 16)
            { "AA/UU", -0.93 }, { "UU/AA", -0.93 },
            { "AU/UA", -1.10 },
            { "UA/AU", -1.33 },
            { "CU/GA", -2.08 }, { "AG/UC", -2.08 },
            { "CA/GU", -2.11 }, { "UG/AC", -2.11 },
            { "GU/CA", -2.24 }, { "AC/UG", -2.24 },
            { "GA/CU", -2.35 }, { "UC/AG", -2.35 },
            { "CG/GC", -2.36 },
            { "GG/CC", -3.26 }, { "CC/GG", -3.26 },
            { "GC/CG", -3.42 },
            // GU wobble stacking (11 unique + 9 symmetric = 20)
            { "AG/UU", -0.55 }, { "UU/GA", -0.55 },
            { "UG/AU", -1.00 }, { "UA/GU", -1.00 },
            { "GA/UU", -1.27 }, { "UU/AG", -1.27 },
            { "AU/UG", -1.36 }, { "GU/UA", -1.36 },
            { "CG/GU", -1.41 }, { "UG/GC", -1.41 },
            { "GG/CU", -1.53 }, { "UC/GG", -1.53 },
            { "CU/GG", -2.11 }, { "GG/UC", -2.11 },
            { "GU/CG", -2.51 }, { "GC/UG", -2.51 },
            { "GG/UU", -0.50 }, { "UU/GG", -0.50 },
            { "UG/GU", +0.30 },
            { "GU/UG", +1.29 },
        }.ToFrozenDictionary();

        /// <summary>
        /// Hairpin loop initiation energies by loop size (kcal/mol at 37°C).
        /// Sizes 3-9 are experimentally determined; 10-30 are extrapolated via
        /// ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9). Source: NNDB — loop.txt.
        /// </summary>
        public static readonly IReadOnlyDictionary<int, double> HairpinLoopInitiation = new Dictionary<int, double>
        {
            {  3, 5.4 }, {  4, 5.6 }, {  5, 5.7 }, {  6, 5.4 }, {  7, 6.0 },
            {  8, 5.5 }, {  9, 6.4 }, { 10, 6.5 }, { 11, 6.6 }, { 12, 6.7 },
            { 13, 6.8 }, { 14, 6.9 }, { 15, 6.9 }, { 16, 7.0 }, { 17, 7.1 },
            { 18, 7.1 }, { 19, 7.2 }, { 20, 7.2 }, { 21, 7.3 }, { 22, 7.3 },
            { 23, 7.4 }, { 24, 7.4 }, { 25, 7.5 }, { 26, 7.5 }, { 27, 7.5 },
            { 28, 7.6 }, { 29, 7.6 }, { 30, 7.7 }
        }.ToFrozenDictionary();

        /// <summary>
        /// Terminal mismatch stacking energies (kcal/mol at 37°C), the NNDB Turner 2004 "Terminal
        /// Mismatches" table (RNAstructure <c>tstack</c>). Key = closing5' + mismatch 3' of closing5' +
        /// mismatch 5' of closing3' + closing3' (helix end 5'-c5 x-3' / 3'-c3 y-5'). Used for exterior and
        /// multibranch helix ends (identical, entry for entry, to ViennaRNA 2.7.2 <c>mismatchExt37</c> /
        /// <c>mismatchM37</c>) and, with the NNDB hairpin bonuses added, for hairpin first mismatches.
        /// Interior loops use their own tables (<c>Turner2004InteriorLoopTables</c>).
        /// Source: NNDB rna_2004_terminal_mismatches.html (turner_2004_terminal_mismatch_dg.txt).
        /// </summary>
        public static readonly IReadOnlyDictionary<string, double> TerminalMismatchEnergies = new Dictionary<string, double>
        {
            // Closing pair A-U (5'AX / 3'UY)
            { "AAAU", -0.8 }, { "AACU", -1.0 }, { "AAGU", -0.8 }, { "AAUU", -1.0 },
            { "ACAU", -0.6 }, { "ACCU", -0.7 }, { "ACGU", -0.6 }, { "ACUU", -0.7 },
            { "AGAU", -0.8 }, { "AGCU", -1.0 }, { "AGGU", -0.8 }, { "AGUU", -1.0 },
            { "AUAU", -0.6 }, { "AUCU", -0.8 }, { "AUGU", -0.6 }, { "AUUU", -0.8 },
            // Closing pair C-G (5'CX / 3'GY)
            { "CAAG", -1.5 }, { "CACG", -1.5 }, { "CAGG", -1.4 }, { "CAUG", -1.5 },
            { "CCAG", -1.0 }, { "CCCG", -1.1 }, { "CCGG", -1.0 }, { "CCUG", -0.8 },
            { "CGAG", -1.4 }, { "CGCG", -1.5 }, { "CGGG", -1.6 }, { "CGUG", -1.5 },
            { "CUAG", -1.0 }, { "CUCG", -1.4 }, { "CUGG", -1.0 }, { "CUUG", -1.2 },
            // Closing pair G-C (5'GX / 3'CY)
            { "GAAC", -1.1 }, { "GACC", -1.5 }, { "GAGC", -1.3 }, { "GAUC", -1.5 },
            { "GCAC", -1.1 }, { "GCCC", -0.7 }, { "GCGC", -1.1 }, { "GCUC", -0.5 },
            { "GGAC", -1.6 }, { "GGCC", -1.5 }, { "GGGC", -1.4 }, { "GGUC", -1.5 },
            { "GUAC", -1.1 }, { "GUCC", -1.0 }, { "GUGC", -1.1 }, { "GUUC", -0.7 },
            // Closing pair G-U (5'GX / 3'UY)
            { "GAAU", -0.3 }, { "GACU", -1.0 }, { "GAGU", -0.8 }, { "GAUU", -1.0 },
            { "GCAU", -0.6 }, { "GCCU", -0.7 }, { "GCGU", -0.6 }, { "GCUU", -0.7 },
            { "GGAU", -0.6 }, { "GGCU", -1.0 }, { "GGGU", -0.8 }, { "GGUU", -1.0 },
            { "GUAU", -0.6 }, { "GUCU", -0.8 }, { "GUGU", -0.6 }, { "GUUU", -0.6 },
            // Closing pair U-A (5'UX / 3'AY)
            { "UAAA", -1.0 }, { "UACA", -0.8 }, { "UAGA", -1.1 }, { "UAUA", -0.8 },
            { "UCAA", -0.7 }, { "UCCA", -0.6 }, { "UCGA", -0.7 }, { "UCUA", -0.5 },
            { "UGAA", -1.1 }, { "UGCA", -0.8 }, { "UGGA", -1.2 }, { "UGUA", -0.8 },
            { "UUAA", -0.7 }, { "UUCA", -0.6 }, { "UUGA", -0.7 }, { "UUUA", -0.5 },
            // Closing pair U-G (5'UX / 3'GY)
            { "UAAG", -1.0 }, { "UACG", -0.8 }, { "UAGG", -1.1 }, { "UAUG", -0.8 },
            { "UCAG", -0.7 }, { "UCCG", -0.6 }, { "UCGG", -0.7 }, { "UCUG", -0.5 },
            { "UGAG", -0.5 }, { "UGCG", -0.8 }, { "UGGG", -0.8 }, { "UGUG", -0.8 },
            { "UUAG", -0.7 }, { "UUCG", -0.6 }, { "UUGG", -0.7 }, { "UUUG", -0.5 },
        }.ToFrozenDictionary();

        /// <summary>NNDB special hairpin loops (tri/tetra/hexaloops): TOTAL loop free energy at 37 °C (kcal/mol) keyed closing5' + loop + closing3'; replaces initiation + mismatch (NNDB hairpin-special-parameters).</summary>
        public static readonly IReadOnlyDictionary<string, double> SpecialHairpinLoops = new Dictionary<string, double>
        {
            // Triloops (closing + 3nt loop + closing = 5-mer key)
            { "CAACG", 6.8 },
            { "GUUAC", 6.9 },
            // Tetraloops (closing + 4nt loop + closing = 6-mer key)
            { "CCUCGG", 2.5 },
            { "CUCCGG", 2.7 },
            { "CUACGG", 2.8 },
            { "CUGCGG", 2.8 },
            { "CCAAGG", 3.3 },
            { "CCCAGG", 3.4 },
            { "CCGAGG", 3.5 },
            { "CUUAGG", 3.5 },
            { "CCGCGG", 3.6 },
            { "CUAAGG", 3.6 },
            { "CCUAGG", 3.7 },
            { "CCACGG", 3.7 },
            { "CUCAGG", 3.7 },
            { "CUUCGG", 3.7 },
            { "CUUUGG", 3.7 },
            { "CAACGG", 5.5 },
            // Hexaloops (closing + 6nt loop + closing = 8-mer key)
            { "ACAGUGUU", 1.8 },
            { "ACAGUACU", 2.8 },
            { "ACAGUGCU", 2.9 },
            { "ACAGUGAU", 3.6 },
        }.ToFrozenDictionary();

        /// <summary>Internal-loop initiation (kcal/mol, 37 °C), n = 4..30 unpaired nucleotides (NNDB turner_2004_internal_loop_dg.txt, = ViennaRNA 2.7.2 internal_loop37).</summary>
        public static readonly IReadOnlyDictionary<int, double> InternalLoopInitiation = new Dictionary<int, double>
        {
            {  4, 1.1 }, {  5, 2.0 }, {  6, 2.0 }, {  7, 2.1 }, {  8, 2.3 },
            {  9, 2.4 }, { 10, 2.5 }, { 11, 2.6 }, { 12, 2.7 }, { 13, 2.8 },
            { 14, 2.9 }, { 15, 2.9 }, { 16, 3.0 }, { 17, 3.1 }, { 18, 3.1 },
            { 19, 3.2 }, { 20, 3.3 }, { 21, 3.3 }, { 22, 3.4 }, { 23, 3.4 },
            { 24, 3.5 }, { 25, 3.5 }, { 26, 3.5 }, { 27, 3.6 }, { 28, 3.6 },
            { 29, 3.7 }, { 30, 3.7 }
        }.ToFrozenDictionary();

        /// <summary>Bulge-loop initiation (kcal/mol, 37 °C), n = 1..30 (NNDB turner_2004 loop.txt BULGE column).</summary>
        public static readonly IReadOnlyDictionary<int, double> BulgeLoopInitiation = new Dictionary<int, double>
        {
            {  1, 3.8 }, {  2, 2.8 }, {  3, 3.2 }, {  4, 3.6 }, {  5, 4.0 },
            {  6, 4.4 }, {  7, 4.6 }, {  8, 4.7 }, {  9, 4.8 }, { 10, 4.9 },
            { 11, 5.0 }, { 12, 5.1 }, { 13, 5.2 }, { 14, 5.3 }, { 15, 5.4 },
            { 16, 5.4 }, { 17, 5.5 }, { 18, 5.5 }, { 19, 5.6 }, { 20, 5.7 },
            { 21, 5.7 }, { 22, 5.8 }, { 23, 5.8 }, { 24, 5.8 }, { 25, 5.9 },
            { 26, 5.9 }, { 27, 6.0 }, { 28, 6.0 }, { 29, 6.0 }, { 30, 6.1 }
        }.ToFrozenDictionary();

        /// <summary>3' dangling-end energies (kcal/mol, 37 °C), key closing5' + dangling base + closing3' (the dangling base 3' of closing5'); NNDB rna_2004_dangling_ends = ViennaRNA 2.7.2 dangle3_37.</summary>
        public static readonly IReadOnlyDictionary<string, double> DanglingEnd3Energies = new Dictionary<string, double>
        {
            // Closing A-U
            { "AAU", -0.8 }, { "ACU", -0.5 }, { "AGU", -0.8 }, { "AUU", -0.6 },
            // Closing C-G
            { "CAG", -1.7 }, { "CCG", -0.8 }, { "CGG", -1.7 }, { "CUG", -1.2 },
            // Closing G-C
            { "GAC", -1.1 }, { "GCC", -0.4 }, { "GGC", -1.3 }, { "GUC", -0.6 },
            // Closing G-U
            { "GAU", -0.8 }, { "GCU", -0.5 }, { "GGU", -0.8 }, { "GUU", -0.6 },
            // Closing U-A
            { "UAA", -0.7 }, { "UCA", -0.1 }, { "UGA", -0.7 }, { "UUA", -0.1 },
            // Closing U-G
            { "UAG", -0.7 }, { "UCG", -0.1 }, { "UGG", -0.7 }, { "UUG", -0.1 },
        }.ToFrozenDictionary();

        /// <summary>5' dangling-end energies (kcal/mol, 37 °C), key closing5' + dangling base + closing3' (the dangling base 5' of closing3'); NNDB rna_2004_dangling_ends = ViennaRNA 2.7.2 dangle5_37.</summary>
        public static readonly IReadOnlyDictionary<string, double> DanglingEnd5Energies = new Dictionary<string, double>
        {
            // Closing A-U
            { "AAU", -0.3 }, { "ACU", -0.1 }, { "AGU", -0.2 }, { "AUU", -0.2 },
            // Closing C-G
            { "CAG", -0.2 }, { "CCG", -0.3 }, { "CGG",  0.0 }, { "CUG",  0.0 },
            // Closing G-C
            { "GAC", -0.5 }, { "GCC", -0.3 }, { "GGC", -0.2 }, { "GUC", -0.1 },
            // Closing G-U
            { "GAU", -0.3 }, { "GCU", -0.1 }, { "GGU", -0.2 }, { "GUU", -0.2 },
            // Closing U-A
            { "UAA", -0.3 }, { "UCA", -0.3 }, { "UGA", -0.4 }, { "UUA", -0.2 },
            // Closing U-G
            { "UAG", -0.3 }, { "UCG", -0.3 }, { "UGG", -0.4 }, { "UUG", -0.2 },
        }.ToFrozenDictionary();


        /// <summary>
        /// Terminal AU/GU penalty (kcal/mol at 37°C), applied at each end of a helix that
        /// terminates with an AU/UA or GU/UG base pair.
        /// Source: NNDB — "Per AU end" (wc-parameters.html), "Per GU end" (gu-parameters.html).
        /// </summary>
        public const double TerminalAuGuPenalty = 0.45;

        /// <summary>
        /// Returns <c>true</c> if the base pair (<paramref name="b1"/>, <paramref name="b2"/>)
        /// is AU/UA or GU/UG and therefore incurs the <see cref="TerminalAuGuPenalty"/> at a
        /// helix terminus. Bases are case-insensitive and DNA T is read as U (ViennaRNA
        /// <c>vrna_nucleotide_encode</c>), so A·T and G·T count as A·U and G·U.
        /// </summary>
        public static bool IsTerminalPenaltyPair(char b1, char b2)
        {
            char x = Normalize(b1), y = Normalize(b2);
            return (x == 'A' && y == 'U') || (x == 'U' && y == 'A') ||
                   (x == 'G' && y == 'U') || (x == 'U' && y == 'G');
        }

        private static char Normalize(char b)
        {
            char u = char.ToUpperInvariant(b);
            return u == 'T' ? 'U' : u;
        }
    }
}
