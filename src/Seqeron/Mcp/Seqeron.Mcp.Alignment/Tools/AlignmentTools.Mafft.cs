using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-19 (stage 1): MAFFT 7.505, the alignment of two sequences (port of the program, verified against the installed binary); delegates to SequenceAligner.MultipleAlignMafft.
public partial class AlignmentTools
{
    [McpServerTool(Name = "multiple_align_mafft", Title = "MSA — MAFFT 7.505, two sequences (port of the program)", ReadOnly = true)]
    [Description("MAFFT 7.505's alignment of TWO sequences (Katoh & Standley 2013, MBE 30:772; mafft --anysymbol in.fa), ported step by step from the program's source (stage 1 of the MAFFT port; three or more sequences need the 6-mer distances, the guide tree and the group-group dynamic program and are refused with the pointer B09.md A2-19b): the alphabet is guessed from the letters (nucleotide when more than 75 % of the letters are A, C, G, T, U or N; sequenceType 'nuc' / 'amino' forces it), with anySymbol (default, as T-Coffee's mafft_pair calls it) letters outside the alphabet - digits, '*', '#', B, Z, J, O, U in proteins - are aligned as X (n for nucleotides) and restored in the output, '-' characters of the input are removed; scoring = the BLOSUM matrix (bl 30, 45, 50, 62 = default, 80, 90, 100 or 0) rescaled to 600 / mean diagonal with the offset (ep) subtracted, for nucleotides the 200-PAM Kimura matrix (kappa 2) with IUPAC codes by averaging, the gap opening penalty op (default 1.53, scaled by 0.6, or 1.8 for nucleotides, to an integer), no gap extension penalty and free end gaps (the script passes -O); dynamic program = G__align11 (a gap-open-only affine recurrence over 'best earlier cell' trackers, tie rules >= for the trackers, > for the cells) with the end of the alignment chosen in the last row / column. retree / fft / reorder / thread are accepted (no effect for two sequences). The default FFT-NS-2 strategy gives the same rows as this dynamic program on every committed pair, up to 1 500 residues. VERIFICATION: rows equal the installed Debian mafft 7.505-1 on 884 / 884 committed pairs (protein, DNA / RNA, IUPAC and unusual symbols, op / ep / bl / nuc / amino, the plain command without --anysymbol where nucleotides come back in lower case and amino acids in upper case). T-Coffee's mafft_pair uses this port.")]
    public static MultipleAlignmentResultDto MultipleAlignMafft(
        [Description("Exactly two sequences (letters; with anySymbol any character except blanks and '-' is a residue).")] string[] sequences,
        [Description("'auto' (default: guessed from the letters), 'nuc' (mafft --nuc) or 'amino' (mafft --amino).")] string sequenceType = "auto",
        [Description("mafft --op: gap opening penalty (default 1.53).")] double gapOpen = 1.53,
        [Description("mafft --ep: offset (default 0.0; a negative value adds to every substitution score, a positive one subtracts).")] double offset = 0.0,
        [Description("mafft --bl: BLOSUM matrix of the amino-acid scoring: 30, 45, 50, 62 (default), 80, 90, 100 or 0.")] int blosum = 62,
        [Description("true (default): mafft --anysymbol, the original letters are restored in the output. false: the plain command (letters and '.' only; nucleotides lower case, amino acids upper case in the rows).")] bool anySymbol = true,
        [Description("mafft --retree (accepted, no effect for two sequences).")] int retree = 2,
        [Description("mafft --fft (true, default) / --nofft (false): accepted, no effect on the committed pairs.")] bool fft = true,
        [Description("mafft --thread (accepted, no effect).")] int threads = 1,
        [Description("mafft --reorder (accepted, no effect for two sequences).")] bool reorder = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequences);
        var options = new global::Seqeron.Genomics.Alignment.MafftOptions
        {
            SequenceType = (sequenceType ?? "auto").Trim().ToLowerInvariant() switch
            {
                "auto" or "" => global::Seqeron.Genomics.Alignment.MafftSequenceType.Auto,
                "nuc" => global::Seqeron.Genomics.Alignment.MafftSequenceType.Nucleotide,
                "amino" => global::Seqeron.Genomics.Alignment.MafftSequenceType.Amino,
                _ => throw new ArgumentException($"sequenceType must be 'auto', 'nuc' or 'amino', got '{sequenceType}'.", nameof(sequenceType)),
            },
            GapOpen = gapOpen,
            Offset = offset,
            Blosum = blosum,
            AnySymbol = anySymbol,
            Retree = retree,
            Fft = fft,
            Threads = threads,
            Reorder = reorder,
            CancellationToken = cancellationToken,
        };
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignMafft(sequences, options);
        return new MultipleAlignmentResultDto(r.AlignedSequences, r.Consensus, r.TotalScore, r.ExactTotalScore);
    }
}
