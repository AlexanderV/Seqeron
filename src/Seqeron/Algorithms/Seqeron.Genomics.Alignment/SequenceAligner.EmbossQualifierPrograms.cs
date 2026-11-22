using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// A5-7 (B09 F172): needle / water / stretcher / matcher with the sequence-input qualifiers (-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein of every input).
// The programs trim both sequences (ajSeqTrim) and work on the trimmed residues; the report adds the cut residues to the offsets (ajSeqGetOffset) and counts down for a reversed sequence.
public static partial class SequenceAligner
{
    private static string MatrixFileName(SubstitutionMatrix matrix) => matrix.Name.StartsWith('E') ? matrix.Name : "E" + matrix.Name;

    private static EmbossPairAlignment ProgramPair(AlignmentResult r, EmbossSequence first, EmbossSequenceView v1, EmbossSequence second, EmbossSequenceView v2)
    {
        var pair = EmbossPairAlignment.From(r, first.Name, v1, second.Name, v2);
        return pair with { First = pair.First with { Description = first.Description, Usa = first.Name }, Second = pair.Second with { Description = second.Description, Usa = second.Name } };
    }

    /// <summary>
    /// EMBOSS <c>needle</c> on two sequences read with their qualifiers: the ranges and orientation of <see cref="EmbossSequence.Qualifiers"/> are applied (<see cref="EmbossInputSequence.Read"/>), the alignment is
    /// <see cref="NeedleAlign"/> of the trimmed sequences and the result reports the offsets / reversal of the program. The default matrix follows the type of the first sequence (forced type included), and the second
    /// sequence must be of the matching type, as the ACD of needle requires.
    /// </summary>
    /// <param name="first"><c>-asequence</c> with its <c>-sbegin1</c> ... qualifiers.</param>
    /// <param name="second"><c>-bsequence</c> with its <c>-sbegin2</c> ... qualifiers.</param>
    /// <param name="matrix">The scoring matrix; null = EDNAFULL / EBLOSUM62 by the first sequence.</param>
    /// <param name="gapOpen"><c>-gapopen</c>.</param>
    /// <param name="gapExtend"><c>-gapextend</c>.</param>
    /// <param name="endWeight"><c>-endweight</c>.</param>
    /// <param name="endOpen"><c>-endopen</c>.</param>
    /// <param name="endExtend"><c>-endextend</c>.</param>
    /// <param name="brief"><c>-brief</c>: false adds the <c>Longest_Identity</c> ... sub-header lines.</param>
    public static EmbossProgramResult NeedleProgram(EmbossSequence first, EmbossSequence second, SubstitutionMatrix? matrix = null, double gapOpen = 10, double gapExtend = 0.5,
        bool endWeight = false, double endOpen = 10, double endExtend = 0.5, bool brief = true)
    {
        var (v1, v2, protein) = EmbossPairInput.Read(first, second);
        matrix ??= protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
        var r = NeedleAlign(v1.Residues, v2.Residues, matrix, gapOpen, gapExtend, endWeight, endOpen, endExtend);
        var pair = ProgramPair(r, first, v1, second, v2);
        if (!brief)
            pair = pair with { ExtraSubHeader = SimilaritySubHeader(r.AlignedSequence1, r.AlignedSequence2, v1.Residues.Length, v2.Residues.Length, matrix) };
        return new EmbossProgramResult([pair], "", MatrixFileName(matrix));
    }

    /// <summary>EMBOSS <c>water</c> on two sequences read with their qualifiers (see <see cref="NeedleProgram"/>); no alignment above score 0 gives no alignment, as the program prints none.</summary>
    public static EmbossProgramResult WaterProgram(EmbossSequence first, EmbossSequence second, SubstitutionMatrix? matrix = null, double gapOpen = 10, double gapExtend = 0.5)
    {
        var (v1, v2, protein) = EmbossPairInput.Read(first, second);
        matrix ??= protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
        var r = WaterAlign(v1.Residues, v2.Residues, matrix, gapOpen, gapExtend);
        if (r.AlignedSequence1.Length == 0)
            return new EmbossProgramResult([], "", MatrixFileName(matrix));
        return new EmbossProgramResult([ProgramPair(r, first, v1, second, v2)], "", MatrixFileName(matrix));
    }

    /// <summary>EMBOSS <c>stretcher</c> on two sequences read with their qualifiers (see <see cref="NeedleProgram"/>); the default penalties are 16 / 4 for nucleotide and 12 / 2 for protein, chosen by the first sequence.</summary>
    public static EmbossProgramResult StretcherProgram(EmbossSequence first, EmbossSequence second, SubstitutionMatrix? matrix = null, int? gapOpen = null, int? gapExtend = null)
    {
        var (v1, v2, protein) = EmbossPairInput.Read(first, second);
        matrix ??= protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
        var r = StretcherAlign(v1.Residues, v2.Residues, matrix, gapOpen ?? (protein ? 12 : 16), gapExtend ?? (protein ? 2 : 4));
        return new EmbossProgramResult([ProgramPair(r, first, v1, second, v2)], "", MatrixFileName(matrix));
    }

    /// <summary>
    /// EMBOSS <c>matcher</c> on two sequences read with their qualifiers (see <see cref="NeedleProgram"/>; the penalties are integers: 16 / 4 for nucleotide, 14 / 4 for protein): the best <paramref name="alternatives"/> local alignments (Huang-Miller SIM, the program's tie order).
    /// </summary>
    /// <exception cref="InvalidOperationException">Fewer than <paramref name="alternatives"/> alignments exist: the program aborts.</exception>
    /// <exception cref="ArgumentException"><paramref name="gapOpen"/> / <paramref name="gapExtend"/> is not an integer (matcher: <c>Invalid integer value</c>), or the score range of the 32-bit program is exceeded.</exception>
    /// <remarks>
    /// Any <c>-gapopen</c> is legal, also one below <c>-gapextend</c> (the first gap residue then costs less than a further one; SIM's <c>q = gapopen − gapextend</c> is negative: measured on matcher 6.6.0 with
    /// <c>-gapopen 0 .. 3 -gapextend 4</c>, A6-3). Negative values are reset to the minimum (0 / 0 / 1) with the program's warning, which is returned in <see cref="EmbossProgramResult.StandardError"/>; a value the program reads as UB (its 32-bit scores
    /// overflow: <c>-gapopen 2147483647</c>) is refused by the integer range check of the aligner.
    /// </remarks>
    public static EmbossProgramResult MatcherProgram(EmbossSequence first, EmbossSequence second, SubstitutionMatrix? matrix = null, double? gapOpen = null, double gapExtend = 4, int alternatives = 1)
    {
        var (v1, v2, protein) = EmbossPairInput.Read(first, second);
        matrix ??= protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
        // matcher.acd: -gapopen / -gapextend / -alternatives are integers with the minima 0 / 0 / 1 (below: warning and reset), no maximum; a gap of length k costs gapopen + (k - 1) * gapextend, also when gapopen is below gapextend
        var warnings = new System.Text.StringBuilder();
        alternatives = EmbossInteger(alternatives, 1, nameof(alternatives), "-alternatives", warnings);              // the order of matcher.acd: alternatives, gapopen, gapextend
        int open = EmbossInteger(gapOpen ?? (protein ? 14 : 16), 0, nameof(gapOpen), "-gapopen", warnings);      // 16 for nucleotide, 14 for protein
        int extend = EmbossInteger(gapExtend, 0, nameof(gapExtend), "-gapextend", warnings);
        var scoring = new MatrixScoring(matrix, -((double)open - extend), -(double)extend);
        var list = LocalAlignAffineTopN(v1.Residues, v2.Residues, scoring, alternatives, 0, TopNTieOrder.Sim);
        if (list.Count < alternatives)
            throw new InvalidOperationException($"EMBOSS matcher aborts: only {list.Count} alignment(s) exist for -alternatives {alternatives}.");
        return new EmbossProgramResult([.. list.Select(r => ProgramPair(r, first, v1, second, v2))], "", MatrixFileName(matrix), warnings.ToString());
    }

    /// <summary>
    /// EMBOSS <c>merger</c> on two sequences read with their qualifiers (see <see cref="NeedleProgram"/>). Measured on 6.6.0: merger does not trim its sequences (the whole, possibly reverse-complemented, sequences are merged and the
    /// whole merged sequence is written; <c>-sbegin</c> / <c>-send</c> have no effect on it), the case options are overridden by the program's own lower-casing, and the conflict positions count the sequences as used.
    /// </summary>
    /// <exception cref="ArgumentException">A sequence is empty after reading, or the forced types contradict the program's input types.</exception>
    public static MergerResult MergerProgram(EmbossSequence first, EmbossSequence second, SubstitutionMatrix? matrix = null, double gapOpen = 50, double gapExtend = 5, CancellationToken cancellationToken = default)
    {
        var (v1, v2, protein) = EmbossPairInput.Read(first, second);
        return Merger(v1.Full, v2.Full, matrix ?? (protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull), gapOpen, gapExtend, cancellationToken);
    }
}
