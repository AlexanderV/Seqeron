using Seqeron.Genomics.Infrastructure;
using System.Text;

namespace Seqeron.Genomics.Analysis;

public static partial class MotifFinder
{
    #region EMBOSS consambig

    // EMBOSS 6.6.0 Ebases.iub / Eresidues.iub binary codes (ajBaseAlphaToBin / ajResidueAlphaToBin): nucleotides A=1 C=2 G=4 T=8 (U = T, X = N = 15),
    // amino acids one bit each in the order of the data file (A=1 C=2 D=4 E=8 F=16 G=32 H=64 I=128 K=256 L=512 M=1024 N=2048 P=4096 Q=8192 R=16384 S=32768
    // T=65536 V=131072 W=262144 Y=524288 U=1048576 O=2097152), B = D|N, Z = E|Q, J = I|L, X = the twenty standard residues.
    private const string AminoBitOrder = "ACDEFGHIKLMNPQRSTVWYUO";
    private const int AminoBitX = (1 << 20) - 1;

    private static int NucleotideBits(char upper)
    {
        if (upper == 'X')
            return 15;
        if (upper == 'U')
            upper = 'T';
        int bits = 0;
        const string bases = "ACGT";
        for (int i = 0; i < bases.Length; i++)
            if (Seqeron.Genomics.Core.IupacHelper.IsNucleotideCode(upper) && Seqeron.Genomics.Core.IupacHelper.MatchesIupac(bases[i], upper))
                bits |= 1 << i;
        return bits;
    }

    private static int AminoBits(char upper) => upper switch
    {
        'B' => AminoBits('D') | AminoBits('N'),
        'Z' => AminoBits('E') | AminoBits('Q'),
        'J' => AminoBits('I') | AminoBits('L'),
        'X' => AminoBitX,
        _ => AminoBitOrder.IndexOf(upper) is var i and >= 0 ? 1 << i : 0,
    };

    /// <summary>
    /// Ambiguity-code consensus of an alignment — EMBOSS <c>consambig</c> (EMBOSS 6.6.0 <c>emboss/consambig.c</c>, <c>ajAlignConsAmbig</c> / <c>ajAlignConsAmbigNuc</c> / <c>ajAlignConsAmbigProt</c> and the
    /// code tables <c>ajBaseAlphaToBin</c>, <c>ajBaseBinToAlpha</c>, <c>ajResidueAlphaToBin</c>, <c>ajResidueBinToAlpha</c>). Per column the binary codes of all residues are OR-ed
    /// (nucleotides: A = 1, C = 2, G = 4, T = 8, so R = A|G, N = any; amino acids: one bit each, B = D|N, Z = E|Q, J = I|L, X = all twenty) and the union is written as
    /// the one letter that stands for exactly that set; a column with a gap (<c>-</c>, also <c>.</c> and <c>~</c>) is written in lower case. A nucleotide union is always an IUPAC code
    /// (a column without a base, such as an all-gap column, gives <c>x</c>, written <c>n</c> as <c>consambig</c> sets the output type to nucleotide); an amino-acid union that is no letter of the
    /// table (for example A|D) gives <c>X</c>, and a column with no scored residue (all gaps, or only <c>*</c>) gives the NUL character of the program's <c>ajResidueBinToAlpha(0)</c>.
    /// Unlike <see cref="GenerateEmbossConsensus(IEnumerable{string}, ConsensusResidueType, float?, int, float?, IReadOnlyList{float}?)"/> (the plurality
    /// <c>cons</c>) there is no vote: a single deviating residue widens the code. Input normalisation is the one of <c>cons</c> (upper case, <c>.</c> <c>~</c> → gap, <c>?</c> → X, X → N in nucleotide rows).
    /// </summary>
    /// <param name="alignedSequences">At least two aligned sequences.</param>
    /// <param name="residueType">Nucleotide, Protein, or Auto (decided from the first sequence as the program does).</param>
    /// <param name="padRaggedRows">Pad rows shorter than the longest with trailing gaps (<c>ajSeqsetFill</c>, what the program does for every aligned set) instead of rejecting them.</param>
    /// <returns>The consensus, one symbol per column.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">Fewer than two sequences (the program dies: "Insufficient sequences"), a null element, unequal lengths without padding, or an invalid character.</exception>
    public static string GenerateEmbossAmbiguousConsensus(
        IEnumerable<string> alignedSequences, ConsensusResidueType residueType = ConsensusResidueType.Auto, bool padRaggedRows = true)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);
        if (residueType is not (ConsensusResidueType.Nucleotide or ConsensusResidueType.Protein or ConsensusResidueType.Auto))
            throw new ArgumentOutOfRangeException(nameof(residueType));

        List<string> rows;
        if (padRaggedRows)
        {
            rows = ReadEmbossRows(MaterializeRows(alignedSequences, nameof(alignedSequences)));
            int len = rows.Count == 0 ? 0 : rows.Max(r => r.Length);
            for (int i = 0; i < rows.Count; i++)
                rows[i] = rows[i].PadRight(len, '-');
        }
        else
        {
            rows = ReadEmbossRows(MaterializeRows(alignedSequences, nameof(alignedSequences)));
            if (rows.Any(r => r.Length != rows[0].Length))
                throw new ArgumentException("All sequences must have the same length.", nameof(alignedSequences));
        }

        if (rows.Count < 2)
            throw new ArgumentException($"Insufficient sequences ({rows.Count}) to create a consensus; EMBOSS consambig requires at least 2.", nameof(alignedSequences));

        // measured on consambig 6.6.0: -snucleotide forces a nucleotide set, but -sprotein does not make a nucleotide-looking set protein: the set type is the type of the first sequence
        bool nucleotide = residueType == ConsensusResidueType.Nucleotide || EmbossSequenceReader.IsNucleotide(rows[0]);
        var chars = NormalizedAlignmentChars(rows, residueType, residueType == ConsensusResidueType.Nucleotide, nameof(alignedSequences));
        int mlen = rows[0].Length;
        var consensus = new StringBuilder(mlen);
        for (int k = 0; k < mlen; k++)
        {
            int union = 0;
            bool gap = false;
            for (int s = 0; s < chars.Length; s++)
            {
                char c = chars[s][k];
                gap |= c == '-';
                union |= nucleotide ? NucleotideBits(c) : AminoBits(c);
            }

            char res = nucleotide ? NucleotideLetter(union) : AminoLetter(union);
            if (gap)
                res = char.ToLowerInvariant(res);
            if (nucleotide && res is 'X' or 'x')
                res = res == 'X' ? 'N' : 'n'; // ajSeqSetNuc on the output sequence
            consensus.Append(res);
        }

        return consensus.ToString();
    }

    private static char NucleotideLetter(int union) => "XACMGRSVTWYHKDBN"[union]; // iubbases of ajbase.c, indexed by the OR-ed code 0..15

    private static char AminoLetter(int union)
    {
        if (union == 0)
            return '\0'; // ajResidueBinToAlpha(0): the first table entry with code 0 is the empty one
        if (union == AminoBitX)
            return 'X';
        foreach (char letter in "ABCDEFGHIJKLMNOPQRSTUVWYZ")
            if (AminoBits(letter) == union)
                return letter;
        return 'X';
    }

    #endregion
}
