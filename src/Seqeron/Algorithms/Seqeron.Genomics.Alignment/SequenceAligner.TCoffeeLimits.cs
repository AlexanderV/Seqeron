namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The documented resource limit of the quadratic pair dynamic programs of the T-Coffee port (audit A5-2): biphasic_pair_wise, glocal_pair_wise, gotoh_pair_wise_lgp, gotoh_pair_wise / _sw (traceback), viterbi_pair_wise and
// the banded Gotoh of the FASTA methods (fast_pair, cfasta_pair_wise, ...) refuse an input whose matrices exceed TCoffeePairDpMaxBytes BEFORE anything is allocated (64-bit arithmetic: no OutOfMemoryException, no 32-bit overflow).
// The program has no limit of its own: measured with the installed t_coffee 13.41 under a 4 GB address-space limit, two random 20 000-base DNA sequences make `-method biphasic_pair` and `-method fast_pair` stop with
// "ERROR: FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)" (exit status 1) after 2 - 3 s, and two 40 000-base sequences make `-method slow_pair -dp_mode gotoh_pair_wise` print the same message and die with a segmentation
// fault (exit status 139) and `-dp_mode viterbi_pair_wise` stop with the message (exit status 1).
// =====================================================================================================================

public static partial class SequenceAligner
{
    /// <summary>
    /// The matrices of ONE pair dynamic program of the T-Coffee port (<c>biphasic_pair</c> 48 bytes per cell: twelve int matrices; <c>glocal_pair_wise</c> 40; <c>gotoh_pair_wise_lgp</c> 32; the traceback of <c>gotoh_pair_wise</c> /
    /// <c>gotoh_pair_wise_sw</c> 4; <c>viterbi_pair_wise</c> 44 (five Viterbi floats, five trace-back ints and the match probability per cell); the FASTA band 20 bytes per (row x flagged diagonal)): 2^31 bytes = 2 GiB, e.g. biphasic_pair up to
    /// about 6 700 x 6 700 residues, the Gotoh traceback up to about 23 000 x 23 000. <c>t_coffee</c> 13.41 has no limit of its own and stops with "FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)" (a segmentation fault for
    /// <c>gotoh_pair_wise</c> at 40 000 x 40 000).
    /// </summary>
    public const long TCoffeePairDpMaxBytes = 1L << 31;

    /// <summary>
    /// The residue contract of the stand-alone T-Coffee pair stages (<c>slow_pair</c>, <c>biphasic_pair</c>, <c>lalign_id_pair</c>): letters only (as <c>MultipleAlignTCoffeeProgram</c> requires: "Sequences must be strings of letters"). The program's FASTA reader
    /// drops a <c>*</c>, a <c>-</c> or a blank (measured with 13.41: the rows are those of the cleaned sequence) and the program dies with "ERROR: COREDUMP" for a non-ASCII letter (U+00E9), so the stage entry points refuse such a character with a
    /// clear <see cref="ArgumentException"/> instead of an <see cref="IndexOutOfRangeException"/> in the scoring matrix.
    /// </summary>
    private static void RequireTcResidues(string sequence, string paramName, string stage)
    {
        for (int i = 0; i < sequence.Length; i++)
            if (!char.IsAsciiLetter(sequence[i]))
                throw new ArgumentException(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"T-Coffee {stage}: the character U+{(int)sequence[i]:X4} at index {i} is not a letter. t_coffee removes '*', '-' and blanks when it reads its input and stops (\"ERROR: COREDUMP\") for a non-ASCII letter; clean the sequence first."), paramName);
    }

    private static void RequireTcPairDp(string method, string keeps, long bytesPerCell, long rows, long columns, string unit, string? alternative = null) =>
        RequireTcPairDpCells(method, keeps, bytesPerCell, rows * columns, string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{rows - 1:N0} x {columns - 1:N0} {unit}"), alternative);

    /// <summary><paramref name="cells"/> cells of <paramref name="bytesPerCell"/> bytes (both are 64-bit: the product is compared by division, the byte count of the message is a 128-bit number).</summary>
    private static void RequireTcPairDpCells(string method, string keeps, long bytesPerCell, long cells, string size, string? alternative = null)
    {
        if (cells <= TCoffeePairDpMaxBytes / bytesPerCell)
            return;
        Int128 bytes = (Int128)cells * bytesPerCell;
        throw ResourceLimit($"T-Coffee {method}", $"{keeps} ({bytesPerCell} bytes per cell)", string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{size} = {bytes:N0} bytes"), TCoffeePairDpMaxBytes,
            nameof(SequenceAligner) + "." + nameof(TCoffeePairDpMaxBytes), "t_coffee itself allocates the same matrices and stops with \"FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)\"",
            alternative ?? "use the slow_pair or lalign_id_pair primary library (linear space) for sequences of this length");
    }
}
