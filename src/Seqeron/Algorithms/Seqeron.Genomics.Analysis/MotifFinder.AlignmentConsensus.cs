using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>Residue type of an alignment passed to <see cref="MotifFinder.GenerateEmbossConsensus(IEnumerable{string}, ConsensusResidueType, float?, int, float?, IReadOnlyList{float}?)"/>.</summary>
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

/// <summary>
/// Alignment-consensus algorithms of reference tools: EMBOSS <c>cons</c> (scoring-matrix
/// plurality consensus), Biopython <c>SummaryInfo.dumb_consensus</c> (majority threshold) and the
/// configurable-threshold form of <see cref="GenerateConsensus(IEnumerable{string})"/>.
/// </summary>
public static partial class MotifFinder
{
    #region EMBOSS cons

    // The EMBOSS 6.6.0 data files EDNAFULL (= NCBI NUC.4.4) and EBLOSUM62 are the canonical tables
    // SubstitutionMatrix.EdnaFull / Blosum62 (Infrastructure; cell-for-cell equal, B09 ALIGN-STATS-001 F14).

    /// <summary>
    /// Creates a consensus sequence from a multiple alignment with the EMBOSS <c>cons</c>
    /// algorithm (EMBOSS 6.6.0 <c>nucleus/embcons.c</c> <c>embConsCalc</c>, Tim Carver 2001; driver
    /// <c>emboss/cons.c</c>, defaults from <c>emboss/acd/cons.acd</c>). Per column:
    /// <list type="number">
    /// <item>each residue i gets the score Σ<sub>j≠i</sub> M(r<sub>i</sub>, r<sub>j</sub>)·w<sub>j</sub>
    /// over the non-gap residues of the column (M = <c>EDNAFULL</c> for nucleotides,
    /// <c>EBLOSUM62</c> for proteins); the highest-scoring residue is the candidate (first maximum,
    /// but a gap incumbent is displaced by a later residue of equal score);</item>
    /// <item>its <i>positive matches</i> = Σ w<sub>j</sub> over the residues j (itself included)
    /// with M(candidate, r<sub>j</sub>) &gt; 0;</item>
    /// <item>if positive matches ≥ <paramref name="plurality"/> the candidate is emitted, otherwise the
    /// no-consensus symbol (<c>N</c> nucleotide, <c>X</c> protein); if positive matches ≤
    /// <paramref name="setcase"/> the symbol is lower-cased (this applies to <c>N</c>/<c>X</c> too);</item>
    /// <item>if <paramref name="identity"/> &gt; 0 and fewer than <paramref name="identity"/>
    /// sequences carry the residue with the most positive matches (ties → more identical weight),
    /// the column is the upper-case no-consensus symbol.</item>
    /// </list>
    /// All score arithmetic is single-precision as in EMBOSS, so outputs are identical to
    /// <c>cons -plurality P -identity I -setcase S</c>.
    /// </summary>
    /// <remarks>
    /// Input normalisation follows the EMBOSS sequence reader of <c>cons -snucleotide</c> /
    /// <c>-sprotein</c>: letters are upper-cased, <c>.</c> and <c>~</c> become the gap <c>-</c>,
    /// <c>X</c> is read as <c>N</c> for nucleotides, and <c>?</c> becomes <c>X</c> (for nucleotides an
    /// unscored code; if it is emitted it is written as <c>N</c>, as <c>cons</c> writes the
    /// consensus through <c>ajSeqSetNuc</c>). Characters absent from the matrix (gaps, <c>*</c> in
    /// DNA, <c>J</c>/<c>O</c>/<c>U</c> in protein) contribute no score and no positive matches, exactly
    /// as <c>embConsCalc</c> treats code 0; they can still be emitted when <paramref name="plurality"/>
    /// ≤ 0. Rows of unequal length are rejected here; the <c>cons</c> program pads them with
    /// trailing gaps first (<c>ajSeqsetFill</c>) — use the overload with <c>padRaggedRows</c>.
    /// Cross-checked character-for-character against the EMBOSS 6.6.0 <c>cons</c> binary
    /// (docs/Evidence/MOTIF-CONS-001-Evidence.md).
    /// </remarks>
    /// <param name="alignedSequences">At least two aligned sequences of equal length.</param>
    /// <param name="residueType">
    /// Selects the matrix (<c>EDNAFULL</c>/<c>EBLOSUM62</c>) and the no-consensus symbol;
    /// <see cref="ConsensusResidueType.Auto"/> decides both from the first sequence as <c>cons</c> does.
    /// </param>
    /// <param name="plurality">
    /// Minimum positive-match weight for a consensus residue (<c>-plurality</c>); default half the
    /// total sequence weight.
    /// </param>
    /// <param name="identity">Required number of identical residues at a position (<c>-identity</c>, default 0 = off).</param>
    /// <param name="setcase">
    /// Positive-match weight at or below which the output is lower case (<c>-setcase</c>); default
    /// half the total sequence weight.
    /// </param>
    /// <param name="weights">
    /// Optional per-sequence weights (MSF <c>Weight:</c>; default 1.0 each); finite and ≥ 0.
    /// </param>
    /// <returns>The consensus, one symbol per alignment column.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// Fewer than two sequences, a null element, unequal lengths, a character that is not a letter,
    /// <c>*</c>, <c>?</c> or gap (<c>-</c> <c>.</c> <c>~</c>), or a weight count different from the
    /// sequence count.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Negative <paramref name="identity"/>, a NaN threshold, or a negative / non-finite weight.
    /// </exception>
    public static string GenerateEmbossConsensus(
        IEnumerable<string> alignedSequences,
        ConsensusResidueType residueType = ConsensusResidueType.Nucleotide,
        float? plurality = null,
        int identity = 0,
        float? setcase = null,
        IReadOnlyList<float>? weights = null)
        => GenerateEmbossConsensus(alignedSequences, false, residueType, plurality, identity, setcase, weights);

    /// <summary>
    /// <see cref="GenerateEmbossConsensus(IEnumerable{string}, ConsensusResidueType, float?, int, float?, IReadOnlyList{float}?)"/>
    /// with the EMBOSS sequence-set padding option: when <paramref name="padRaggedRows"/> is true,
    /// rows shorter than the longest are padded at the end with gaps (<c>-</c>), as
    /// <c>ajSeqsetFill</c> (EMBOSS 6.6.0 <c>ajax/core/ajseq.c</c>) does for every <c>aligned: "Y"</c>
    /// sequence-set input of <c>cons</c> (<c>ajax/acd/ajacd.c</c>) before <c>embConsCalc</c>;
    /// when false, rows of unequal length are rejected (the original overload's contract).
    /// </summary>
    /// <param name="alignedSequences">At least two aligned sequences.</param>
    /// <param name="padRaggedRows">Pad shorter rows with trailing gaps instead of rejecting them.</param>
    /// <param name="residueType">Nucleotide, Protein, or Auto (decided from the first sequence like <c>cons</c>).</param>
    /// <param name="plurality">See the primary overload.</param>
    /// <param name="identity">See the primary overload.</param>
    /// <param name="setcase">See the primary overload.</param>
    /// <param name="weights">See the primary overload.</param>
    /// <returns>The consensus, one symbol per (padded) alignment column.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">As the primary overload (unequal lengths only when not padding).</exception>
    /// <exception cref="ArgumentOutOfRangeException">As the primary overload.</exception>
    public static string GenerateEmbossConsensus(
        IEnumerable<string> alignedSequences,
        bool padRaggedRows,
        ConsensusResidueType residueType = ConsensusResidueType.Nucleotide,
        float? plurality = null,
        int identity = 0,
        float? setcase = null,
        IReadOnlyList<float>? weights = null)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);
        if (residueType is not (ConsensusResidueType.Nucleotide or ConsensusResidueType.Protein or ConsensusResidueType.Auto))
            throw new ArgumentOutOfRangeException(nameof(residueType));
        ArgumentOutOfRangeException.ThrowIfNegative(identity);
        if (plurality is float p && float.IsNaN(p))
            throw new ArgumentOutOfRangeException(nameof(plurality), "Plurality cannot be NaN.");
        if (setcase is float sc && float.IsNaN(sc))
            throw new ArgumentOutOfRangeException(nameof(setcase), "Setcase cannot be NaN.");

        List<string> rows;
        if (padRaggedRows)
        {
            rows = ReadEmbossRows(MaterializeRows(alignedSequences, nameof(alignedSequences)));
            int len = rows.Count == 0 ? 0 : rows.Max(r => r.Length);
            for (int i = 0; i < rows.Count; i++)
                rows[i] = rows[i].PadRight(len, '-'); // ajSeqsetFill: append '-' × (Len − own length)
        }
        else
        {
            rows = ReadEmbossRows(MaterializeRows(alignedSequences, nameof(alignedSequences)));
            if (rows.Any(r => r.Length != rows[0].Length))
                throw new ArgumentException("All sequences must have the same length.", nameof(alignedSequences));
        }

        if (rows.Count < 2)
            throw new ArgumentException(
                $"Insufficient sequences ({rows.Count}) to create a consensus; EMBOSS cons requires at least 2.",
                nameof(alignedSequences));
        return EmbossConsensusCore(rows, residueType, plurality, identity, setcase, weights, gaps: false, keepCase: false, writeAsNucleotide: true, null, nameof(alignedSequences));
    }

    /// <summary>
    /// The EMBOSS <c>embConsCalc</c> consensus (<c>nucleus/embcons.c</c>) of already materialised rows, shared by <c>cons</c> and <c>showalign</c>.
    /// </summary>
    /// <param name="rows">The materialised, equal-length sequence rows (upper / lower case as given; gap characters allowed).</param>
    /// <param name="residueType">Nucleotide or protein scoring of the rows (<see cref="ConsensusResidueType"/>).</param>
    /// <param name="plurality"><c>-plurality</c>: minimum positive-match weight (null = half the total weight).</param>
    /// <param name="identity"><c>-identity</c>: required number of identical residues at a position (0 = off).</param>
    /// <param name="setcase"><c>-setcase</c>: weight at or below which the consensus letter is lower case (null = half the total weight).</param>
    /// <param name="weights">One weight per row (null = all 1).</param>
    /// <param name="gaps"><c>embConsCalc</c>'s <c>gaps</c>: gap characters count as a residue (code 0) in the identical / matching weights, so a column can come out as a gap (showalign <c>-gaps</c>, default true there; <c>cons</c> passes false).</param>
    /// <param name="keepCase">Keep the case of the input letters (showalign does not upper-case the set; <c>cons</c> does): the consensus letter is the letter of the winning row.</param>
    /// <param name="writeAsNucleotide"><c>cons</c> writes the consensus through <c>ajSeqSetNuc</c>, which exchanges <c>x</c> for <c>n</c>; showalign shows the raw string.</param>
    /// <param name="matrixOverride">Scoring matrix (showalign <c>-matrix</c>); null = EDNAFULL for a nucleotide set, EBLOSUM62 otherwise.</param>
    /// <param name="paramName">Parameter name reported by the argument exceptions of the calling program (<c>nameof</c> of its own parameter).</param>
    internal static string EmbossConsensusCore(
        List<string> rows, ConsensusResidueType residueType, float? plurality, int identity, float? setcase, IReadOnlyList<float>? weights,
        bool gaps, bool keepCase, bool writeAsNucleotide, SubstitutionMatrix? matrixOverride, string paramName)
    {
        int nseqs = rows.Count;

        var w = new float[nseqs];
        if (weights is null)
        {
            Array.Fill(w, 1.0f);
        }
        else
        {
            if (weights.Count != nseqs)
                throw new ArgumentException("One weight per sequence is required.", nameof(weights));
            for (int i = 0; i < nseqs; i++)
            {
                float wi = weights[i];
                if (!float.IsFinite(wi) || wi < 0f)
                    throw new ArgumentOutOfRangeException(nameof(weights), "Weights must be finite and non-negative.");
                w[i] = wi;
            }
        }

        // ajSeqsetGetTotweight: float accumulation; cons.acd default "@( $(sequence.totweight) / 2)".
        float totalWeight = 0f;
        for (int i = 0; i < nseqs; i++)
            totalWeight += w[i];
        float fplural = plurality ?? totalWeight / 2f;
        float fsetcase = setcase ?? totalWeight / 2f;

        bool nucleotide = residueType == ConsensusResidueType.Auto
            ? EmbossSequenceReader.IsNucleotide(rows[0])
            : residueType == ConsensusResidueType.Nucleotide;
        var table = matrixOverride ?? (nucleotide ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62);
        string labels = table.Alphabet;
        var matrix = new int[labels.Length, labels.Length];
        for (int a = 0; a < labels.Length; a++)
            for (int b = 0; b < labels.Length; b++)
                matrix[a, b] = table.GetScore(labels[a], labels[b]);
        char nocon = nucleotide ? 'N' : 'X';

        // ajSeqcvtNewStr: label i → code i + 1 (both cases); every other character → 0.
        var codeOf = new int[128];
        for (int i = 0; i < labels.Length; i++)
        {
            codeOf[char.ToUpperInvariant(labels[i])] = i + 1;
            codeOf[char.ToLowerInvariant(labels[i])] = i + 1;
        }

        int matsize = labels.Length + 1;

        var chars = NormalizedAlignmentChars(rows, residueType, nucleotide, paramName, keepCase);
        int mlen = rows[0].Length;
        var codes = new int[nseqs][];
        for (int s = 0; s < nseqs; s++)
        {
            codes[s] = new int[mlen];
            for (int k = 0; k < mlen; k++)
                codes[s][k] = codeOf[chars[s][k]];
        }

        var identical = new float[matsize];
        var matching = new float[matsize];
        var score = new float[nseqs];
        var consensus = new StringBuilder(mlen);

        for (int k = 0; k < mlen; k++)
        {
            if ((k & 1023) == 0)
                AlignmentCancellation.ThrowIfRequested();
            char res = nocon;
            Array.Clear(identical);
            Array.Clear(matching);
            Array.Clear(score);

            // Column scores (gaps = false in cons: code 0 contributes nothing; embConsCalc with gaps counts the gap code too).
            for (int i = 0; i < nseqs; i++)
            {
                int m1 = codes[i][k];
                if (m1 != 0 || gaps)
                    identical[m1] += w[i];

                for (int j = i + 1; j < nseqs; j++)
                {
                    int m2 = codes[j][k];
                    if (m1 != 0 && m2 != 0)
                    {
                        float mv = matrix[m1 - 1, m2 - 1];
                        float contri = mv * w[j] + score[i];
                        float contrj = mv * w[i] + score[j];
                        score[i] = contri;
                        score[j] = contrj;
                    }
                }
            }

            int highindex = -1;
            float max = -(float)int.MaxValue;
            for (int i = 0; i < nseqs; i++)
            {
                if (score[i] > max || (score[i] == max && highindex >= 0 && chars[highindex][k] == '-'))
                {
                    highindex = i;
                    max = score[i];
                }
            }

            // Positive matches of each residue in the column (the residue itself included).
            for (int i = 0; i < nseqs; i++)
            {
                int m1 = codes[i][k];
                if (matching[m1] != 0f)
                    continue;
                for (int j = 0; j < nseqs; j++)
                {
                    int m2 = codes[j][k];
                    if (m1 != 0 && m2 != 0 && matrix[m1 - 1, m2 - 1] > 0)
                        matching[m1] += w[j];
                    if (gaps && m1 == 0 && m2 == 0)
                        matching[m1] += w[j];
                }
            }

            int matchingMaxIndex = 0;
            for (int i = 0; i < nseqs; i++)
            {
                int m1 = codes[i][k];
                if (matching[m1] > matching[matchingMaxIndex])
                    matchingMaxIndex = m1;
                else if (matching[m1] == matching[matchingMaxIndex] && identical[m1] > identical[matchingMaxIndex])
                    matchingMaxIndex = m1;
            }

            // Plurality check on the highest-scoring residue.
            int hm = codes[highindex][k];
            if (matching[hm] >= fplural)
                res = chars[highindex][k];
            if (matching[hm] <= fsetcase)
                res = char.ToLowerInvariant(res);

            if (identity > 0)
            {
                int same = 0;
                for (int i = 0; i < nseqs; i++)
                {
                    if (codes[i][k] == matchingMaxIndex)
                        same++;
                }
                if (same < identity)
                    res = nocon;
            }

            // cons.c writes the consensus as a nucleotide sequence when ajSeqsetIsNuc, and
            // ajSeqSetNuc exchanges x/X for n/N: an unscored X (from '?' under -snucleotide, or
            // from a protein-looking row under Auto) emitted with plurality <= 0 is written as N.
            if (writeAsNucleotide && nucleotide && res is 'X' or 'x')
                res = res == 'X' ? 'N' : 'n';

            consensus.Append(res);
        }

        return consensus.ToString();
    }

    /// <summary>The rows of an alignment as the EMBOSS file reader keeps them (<see cref="EmbossSequenceReader.CleanGapped(string)"/>: letters and <c>* . ~ ? # + -</c>; digits, blanks and every other character vanish, which can make rows unequal: <c>ajSeqsetFill</c> pads them).</summary>
    internal static List<string> ReadEmbossRows(List<string> rows)
    {
        for (int i = 0; i < rows.Count; i++)
            rows[i] = EmbossSequenceReader.CleanGapped(rows[i]);
        return rows;
    }

    /// <summary>
    /// The alignment as the EMBOSS sequence reader hands it to <c>cons</c> / <c>consambig</c>: upper case, gaps <c>.</c> <c>~</c> → <c>-</c>, <c>?</c> → <c>X</c> (<c>N</c> in a nucleotide-typed row), <c>X</c> → <c>N</c> in a nucleotide row.
    /// Auto: every row is typed on its own (ajSeqType); an explicit type applies to all rows.
    /// </summary>
    internal static char[][] NormalizedAlignmentChars(List<string> rows, ConsensusResidueType residueType, bool nucleotide, string paramName, bool keepCase = false)
    {
        int mlen = rows[0].Length;
        var chars = new char[rows.Count][];
        for (int s = 0; s < rows.Count; s++)
        {
            AlignmentCancellation.ThrowIfRequested();
            chars[s] = new char[mlen];
            // Auto: every sequence is typed on its own when read (ajSeqType -> ajSeqSetNuc turns
            // X into N in a nucleotide-looking sequence only); explicit types apply to all rows.
            bool rowNucleotide = residueType == ConsensusResidueType.Auto
                ? EmbossSequenceReader.IsNucleotide(rows[s])
                : nucleotide;
            for (int k = 0; k < mlen; k++)
                chars[s][k] = NormalizeEmbossResidue(rows[s][k], rowNucleotide, residueType == ConsensusResidueType.Auto, s, k, paramName, keepCase);
        }

        return chars;
    }

    /// <summary>
    /// EMBOSS sequence-reader normalisation of one alignment character (ajSeqsetFmtUpper, gap
    /// conversion <c>.~</c> → <c>-</c>, <c>ajseqtype.c</c> ambiguity conversion).
    /// </summary>
    /// <param name="c">Input character.</param>
    /// <param name="nucleotide">The sequence is read as nucleotide (ajSeqSetNuc: X → N).</param>
    /// <param name="typedOnRead">
    /// The type was found by ajSeqType after the <c>gapany</c> conversion (Auto), so <c>?</c> → X →
    /// N; with an explicit <c>-snucleotide</c> the sequence is set nucleotide first
    /// (ajSeqTypeCheckIn) and the later <c>?</c> → X stays an unscored X.
    /// </param>
    /// <param name="row">Row index (error message).</param>
    /// <param name="column">Column index (error message).</param>
    /// <param name="paramName">Parameter name (error message).</param>
    /// <param name="keepCase">Keep a lower-case letter lower case (showalign does not upper-case the alignment).</param>
    private static char NormalizeEmbossResidue(char c, bool nucleotide, bool typedOnRead, int row, int column, string paramName, bool keepCase = false)
    {
        if (c is '-' or '.' or '~')
            return '-';
        if (c == '*')
            return '*';
        if (c == '?')                     // (# and + are read as ? by EmbossSequenceReader.CleanGapped)
            return nucleotide && typedOnRead ? 'N' : 'X';
        bool lower = c is >= 'a' and <= 'z';
        if (lower)
            c = (char)(c - 'a' + 'A');
        if (c is >= 'A' and <= 'Z')
        {
            char letter = nucleotide && c == 'X' ? 'N' : c;
            return lower && keepCase ? char.ToLowerInvariant(letter) : letter;
        }
        throw new ArgumentException(
            $"Invalid alignment character '{c}' at position {column} in sequence {row}.",
            paramName);
    }

    #endregion

    #region Biopython dumb_consensus

    /// <summary>
    /// Majority-threshold consensus with the semantics of Biopython
    /// <c>Bio.Align.AlignInfo.SummaryInfo.dumb_consensus(threshold, ambiguous, require_multiple)</c>
    /// (Biopython ≤ 1.85; deprecated in 1.82 and removed in later releases). Per column the
    /// residues other than the gaps <c>-</c> and <c>.</c> are counted (case-sensitive, any
    /// alphabet — DNA, RNA or protein); if exactly one residue has the maximum count and that count
    /// divided by the number of non-gap residues is ≥ <paramref name="threshold"/>, it is emitted;
    /// otherwise (ties, below threshold, all-gap column, or a single non-gap residue when
    /// <paramref name="requireMultiple"/> is set) <paramref name="ambiguous"/> is emitted.
    /// </summary>
    /// <param name="alignedSequences">Aligned sequences of equal length.</param>
    /// <param name="threshold">Required fraction of the non-gap residues (Biopython default 0.7).</param>
    /// <param name="ambiguous">Symbol for columns without consensus (Biopython default <c>X</c>).</param>
    /// <param name="requireMultiple">
    /// When true, a column with exactly one non-gap residue gives <paramref name="ambiguous"/>.
    /// </param>
    /// <returns>The consensus; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A null element or sequences of unequal length.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> is NaN.</exception>
    public static string GenerateDumbConsensus(
        IEnumerable<string> alignedSequences,
        double threshold = 0.7,
        char ambiguous = 'X',
        bool requireMultiple = false)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);
        if (double.IsNaN(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), "Threshold cannot be NaN.");

        List<string> rows = MaterializeAligned(alignedSequences, nameof(alignedSequences));
        if (rows.Count == 0) return "";

        int length = rows[0].Length;
        var consensus = new StringBuilder(length);
        var counts = new Dictionary<char, int>();

        for (int n = 0; n < length; n++)
        {
            counts.Clear();
            int numAtoms = 0;
            foreach (string row in rows)
            {
                char c = row[n];
                if (c != '-' && c != '.')
                {
                    counts[c] = counts.GetValueOrDefault(c) + 1;
                    numAtoms++;
                }
            }

            char best = '\0';
            int maxSize = 0;
            int maxAtoms = 0;
            foreach (var (atom, count) in counts)
            {
                if (count > maxSize)
                {
                    best = atom;
                    maxSize = count;
                    maxAtoms = 1;
                }
                else if (count == maxSize)
                {
                    maxAtoms++;
                }
            }

            if (requireMultiple && numAtoms == 1)
                consensus.Append(ambiguous);
            else if (maxAtoms == 1 && (double)maxSize / numAtoms >= threshold)
                consensus.Append(best);
            else
                consensus.Append(ambiguous);
        }

        return consensus.ToString();
    }

    #endregion

    #region Configurable IUPAC threshold

    /// <summary>
    /// <see cref="GenerateConsensus(IEnumerable{string})"/> with a configurable per-base inclusion
    /// threshold: at each column the bases whose count is strictly greater than
    /// <paramref name="inclusionThreshold"/> × n (n = number of sequences) form the NC-IUB symbol;
    /// when no base passes, the bases tied at the maximum count are encoded (a column without
    /// A/C/G/T gives <c>N</c>). <paramref name="inclusionThreshold"/> = 0.25 is bit-identical to the
    /// parameterless overload; 0 includes every base present; values ≥ the maximum column
    /// frequency reduce every column to its (tied) most frequent bases.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <param name="inclusionThreshold">Per-base frequency cut, in [0, 1].</param>
    /// <returns>Consensus over the 15 IUPAC symbols; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A null element, or sequences of unequal length.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="inclusionThreshold"/> is NaN or outside [0, 1].</exception>
    public static string GenerateConsensus(IEnumerable<string> sequences, double inclusionThreshold)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (double.IsNaN(inclusionThreshold) || inclusionThreshold < 0 || inclusionThreshold > 1)
            throw new ArgumentOutOfRangeException(nameof(inclusionThreshold), "Threshold must be in [0, 1].");

        return GenerateConsensusCore(sequences, inclusionThreshold);
    }

    #endregion
}
