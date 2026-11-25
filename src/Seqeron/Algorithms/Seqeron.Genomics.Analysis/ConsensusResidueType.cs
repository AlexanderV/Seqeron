namespace Seqeron.Genomics.Analysis;

/// <summary>Residue type of an alignment passed to the consensus and statistics functions.</summary>
public enum ConsensusResidueType
{
    /// <summary>Nucleotides: EMBOSS <c>EDNAFULL</c> matrix, no-consensus symbol <c>N</c>.</summary>
    Nucleotide,

    /// <summary>Amino acids: EMBOSS <c>EBLOSUM62</c> matrix, no-consensus symbol <c>X</c>.</summary>
    Protein,

    /// <summary>
    /// Decided from the alignment as the EMBOSS 6.6.0 <c>cons</c> program does without
    /// <c>-snucleotide</c>/<c>-sprotein</c>: the sequence-set type is the type of the first
    /// sequence (<c>ajSeqsetFromList</c>/<c>ajSeqType</c>), which is nucleotide when every
    /// character of that sequence is in <c>ACGTU</c> + <c>BDHKMNRSVWXY?</c> + <c>.~-</c>
    /// (case-insensitive, <c>ajSeqTypeGapnucS</c>; an empty first sequence counts as nucleotide),
    /// otherwise protein. <c>ajSeqsetIsNuc</c> then selects the no-consensus symbol and
    /// <c>$(acdprotein)</c> the matrix. As in the EMBOSS reader, each sequence is also typed on its
    /// own: <c>?</c> becomes <c>X</c> (<c>gapany</c>), and in a nucleotide-looking sequence
    /// <c>ajSeqSetNuc</c> then turns <c>X</c> into <c>N</c>; other sequences keep <c>X</c>.
    /// </summary>
    Auto,
}
