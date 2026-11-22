using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-19 / A2-19b: MAFFT 7.505, progressive strategies (port of the program, verified against the installed binary); delegates to SequenceAligner.RunMafft.
public partial class AlignmentTools
{
    [McpServerTool(Name = "multiple_align_mafft", Title = "MSA — MAFFT 7.505, progressive strategies FFT-NS-1 / -2, NW-NS-1 / -2 (port of the program)", ReadOnly = true)]
    [Description("MAFFT 7.505 (Katoh & Standley 2013, MBE 30:772; mafft [--anysymbol] [options] in.fa), ported step by step from the program's source, for ANY number of sequences: the progressive strategies FFT-NS-2 (default), FFT-NS-1 (retree 1), NW-NS-2 (fft false) and NW-NS-1 (fft false, retree 1; retree 3 builds the guide tree three times). The pipeline is the program's (script mafft -> disttbfast): the alphabet is guessed from the letters (nucleotide when more than 75 % of the letters are A, C, G, T, U or N; sequenceType 'nuc' / 'amino' forces it); with anySymbol letters outside the alphabet - digits, '*', '#', B, Z, J, O, U in proteins - are aligned as X (n for nucleotides) and restored in the output (the plain command stops on them), '-' characters of the input are removed; the distance of two sequences is the number of shared 6-mers (10-mers with tupleSize 10, nucleotides only; amino acids in the six classes of the program) scaled by the lengths, the guide tree is the UPGMA of MAFFT (mixed linkage: treeMixing 0.1 of average and 0.9 of minimum linkage; 1 = --averagelinkage, 0 = --minimumlinkage) with the weights of the sequences from the branch lengths (unweight switches them off), the groups are aligned by the fast Fourier transform of Falign (the polarity and volume of the amino acids, the letter 'a' of nucleotides: 20 candidate shifts, windows with a high score as anchors, a chain of anchors, the segments between them aligned by the dynamic program) when fft is true, and by the group-group dynamic program A__align (gap penalties scaled by the gap frequencies of the columns of the two groups, history of the profiles of the children for groups of more than 20 sequences; legacyGapPenalty = --leavegappyregion switches the scaling off) or the dynamic program G__align11 for two single sequences; for retree 2 the distances are then recomputed from the first alignment, a second tree is built and the sequences are aligned again. Scoring = the BLOSUM matrix (bl 30, 45, 50, 62 = default, 80, 90, 100 or 0), the JTT or the transmembrane PAM matrix (scoringModel 'jtt' / 'tm' with pam steps, default 200) rescaled to 600 / mean diagonal with the offset (ep) subtracted, for nucleotides the Kimura matrix (kappa 2, 200 PAM; pam sets the distance; frequencyModel uses the letter frequencies of the data; wildcardN scores n like a wildcard) with IUPAC codes by averaging; gap opening penalty op (default 1.53, scaled by 0.6, or 1.8 for nucleotides, to an integer; opdist = --opdist for the distance step), gap extension exp (default 0), free end gaps unless terminalGapPenalty. reorder returns the sequences in the order of the guide tree (--reorder). thread is accepted: the program's alignment does not depend on it; its progress text with threads is thread-order dependent and is written as of one thread. The result carries the program's standard error (the banner lines of disttbfast and the strategy text of the script) unless quiet. Not available (refused): --retree 0, the memory-saving mode for groups longer than 30 000 columns; the iterative (maxiterate), pair-library (localpair, globalpair, genafpair), parttree, add / merge / seed strategies are registered stages of the port. VERIFICATION: rows and standard error equal the installed Debian mafft 7.505-1 on every committed run of scripts/alignment/generate_mafft_progressive_oracle.py (6 530 sets in 28 strategies, at least 210 per strategy: protein and DNA, 3 - 45 sequences, tiny and unrelated sequences, unusual symbols, every option above), and on 884 / 884 pairs (generate_mafft_pair_oracle.py). T-Coffee's mafft_pair uses this port.")]
    public static MafftRunResultDto MultipleAlignMafft(
        [Description("Two or more sequences (letters; with anySymbol any character except blanks and '-' is a residue). Like the program: an empty list gives an empty result and ONE sequence is written back (gaps dropped; without anySymbol nucleotides in lower and amino acids in upper case).")] string[] sequences,
        [Description("'auto' (default: guessed from the letters), 'nuc' (mafft --nuc) or 'amino' (mafft --amino).")] string sequenceType = "auto",
        [Description("mafft --op: gap opening penalty (default 1.53).")] double gapOpen = 1.53,
        [Description("mafft --ep: offset (default 0.0; a negative value adds to every substitution score, a positive one subtracts).")] double offset = 0.0,
        [Description("mafft --bl: BLOSUM matrix of the amino-acid scoring: 30, 45, 50, 62 (default), 80, 90, 100 or 0.")] int blosum = 62,
        [Description("true (default): mafft --anysymbol, the original letters are restored in the output. false: the plain command (letters and '.' only; nucleotides lower case, amino acids upper case in the rows).")] bool anySymbol = true,
        [Description("mafft --retree: how many times the guide tree is built (1 = FFT-NS-1, 2 = FFT-NS-2 default, 3; 0 is refused). Two sequences always use 1.")] int retree = 2,
        [Description("mafft --fft (true, default) / --nofft (false, NW-NS).")] bool fft = true,
        [Description("mafft --thread (default 0; accepted, no effect on the alignment).")] int threads = 0,
        [Description("mafft --reorder: the rows in the order of the guide tree.")] bool reorder = false,
        [Description("mafft --exp: gap extension penalty (default 0).")] double gapExtension = 0.0,
        [Description("mafft --opdist: gap opening penalty of the distance step (default = gapOpen).")] double? gapOpenDistance = null,
        [Description("'blosum' (default; mafft --bl), 'jtt' (--jtt / --kimura) or 'tm' (--tm): the scoring matrix of amino acids.")] string scoringModel = "blosum",
        [Description("PAM steps of --jtt / --tm / --kimura (default 200).")] int pam = 200,
        [Description("mafft --kappa: transition / transversion ratio of the nucleotide matrix (default: the program's 2).")] int? kappa = null,
        [Description("mafft --fmodel: letter frequencies of the data in the scoring matrix.")] bool frequencyModel = false,
        [Description("mafft --nwildcard: n scores like a wildcard in nucleotide sequences.")] bool wildcardN = false,
        [Description("mafft --fft given explicitly: the FFT is used for every pair of groups.")] bool forceFft = false,
        [Description("mafft --6merpair (6, default) / --10merpair (10, nucleotides only).")] int tupleSize = 6,
        [Description("mafft --mixedlinkage x (0 = --minimumlinkage, 1 = --averagelinkage; default 0.1).")] double? treeMixing = null,
        [Description("mafft --leavegappyregion / --legacygappenalty.")] bool legacyGapPenalty = false,
        [Description("mafft --termgappenalty: penalise gaps at the ends.")] bool terminalGapPenalty = false,
        [Description("mafft --unweight: no weighting of the sequences by the guide tree.")] bool unweight = false,
        [Description("mafft --quiet: no standard error text.")] bool quiet = false,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
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
            GapExtension = gapExtension,
            GapOpenDistance = gapOpenDistance,
            ScoringModel = (scoringModel ?? "blosum").Trim().ToLowerInvariant() switch
            {
                "blosum" or "" => global::Seqeron.Genomics.Alignment.MafftScoringModel.Blosum,
                "jtt" or "kimura" => global::Seqeron.Genomics.Alignment.MafftScoringModel.Jtt,
                "tm" => global::Seqeron.Genomics.Alignment.MafftScoringModel.Tm,
                _ => throw new ArgumentException($"scoringModel must be 'blosum', 'jtt' or 'tm', got '{scoringModel}'.", nameof(scoringModel)),
            },
            Pam = pam,
            Kappa = kappa,
            FrequencyModel = frequencyModel,
            WildcardN = wildcardN,
            ForceFft = forceFft,
            TupleSize = tupleSize,
            TreeMixing = treeMixing,
            LegacyGapPenalty = legacyGapPenalty,
            TerminalGapPenalty = terminalGapPenalty,
            Unweight = unweight,
            Quiet = quiet,
            CancellationToken = cancellationToken,
        };
        var run = global::Seqeron.Genomics.Alignment.SequenceAligner.RunMafft(sequences, options);
        var r = run.Alignment;
        var dto = new MafftRunResultDto(r.AlignedSequences, r.Consensus, r.TotalScore, r.ExactTotalScore, run.StandardError);
        if (string.IsNullOrWhiteSpace(format))
            return dto;
        double score = dto.ExactTotalScore != 0 ? dto.ExactTotalScore : dto.TotalScore;
        return dto with { Formatted = FormatRowsForOutput(format, dto.AlignedSequences, names, score) };
    }
}
