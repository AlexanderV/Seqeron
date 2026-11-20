using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// EMBOSS 6.6 needleall (audit A2-28c, B09) and the needle / needleall -nobrief similarity sub-header (nucleus/embaln.c embAlignCalcSimilarity).
public static partial class SequenceAligner
{
    /// <summary>
    /// EMBOSS <c>needleall</c>: the global alignment (<see cref="NeedleAlign"/>: the program's own Needleman-Wunsch with affine gaps, end-gap options and tie order) of every sequence of the set
    /// <paramref name="setA"/> against every sequence of <paramref name="setB"/>, in the order of the program: for each B sequence, every A sequence. A pair whose score is not above <paramref name="minScore"/> is not reported;
    /// the program writes <c>Alignment score (S) is less than minimum score(M) for sequences A vs B</c> to its error file instead (returned in <see cref="EmbossProgramResult.ErrorText"/>).
    /// </summary>
    /// <param name="setA">The <c>-asequence</c> set (first member of every pair).</param>
    /// <param name="setB">The <c>-bsequence</c> sequences (second member of every pair).</param>
    /// <param name="matrix">Residue-pair scores (<c>-datafile</c>); default EDNAFULL when every sequence of <paramref name="setA"/> is nucleotide, EBLOSUM62 otherwise (the program decides on <c>-asequence</c> alone).</param>
    /// <param name="gapOpen">Cost of a gap of length one (<c>-gapopen</c>, default 10).</param>
    /// <param name="gapExtend">Cost of each further gap position (<c>-gapextend</c>, default 0.5).</param>
    /// <param name="endWeight"><c>-endweight</c>: apply the end-gap penalties.</param>
    /// <param name="endOpen"><c>-endopen</c> (default 10).</param>
    /// <param name="endExtend"><c>-endextend</c> (default 0.5).</param>
    /// <param name="minScore"><c>-minscore</c>: minimum score to report an alignment (default 0: the program's acd says 1.0 in its help text but the default value is 0).</param>
    /// <param name="brief"><c>-brief</c> (default true); false adds the <c>Longest_Identity</c> / <c>Longest_Similarity</c> / <c>Shortest_Identity</c> / <c>Shortest_Similarity</c> lines to the sub-header.</param>
    public static EmbossProgramResult NeedleAll(IReadOnlyList<EmbossSequence> setA, IReadOnlyList<EmbossSequence> setB, SubstitutionMatrix? matrix = null, double gapOpen = 10, double gapExtend = 0.5,
        bool endWeight = false, double endOpen = 10, double endExtend = 0.5, double minScore = 0, bool brief = true)
    {
        ArgumentNullException.ThrowIfNull(setA);
        ArgumentNullException.ThrowIfNull(setB);
        IReadOnlyList<EmbossSequenceView>? viewsA = null, viewsB = null;
        bool? proteinByFirst = null;
        if (EmbossSetInput.Any(setA, setB))
        {
            // A5-7: the qualifiers (-sbegin -send -sreverse ...) of the inputs: asequence has the ACD type any, bsequence stopprotein / nucleotide by the first sequence of asequence (acdprotein)
            (setA, viewsA) = EmbossSetInput.ReadAll(setA, EmbossSequenceType.Any, EmbossInputKind.Set);
            proteinByFirst = viewsA.Count > 0 && !viewsA[0].IsNucleotide;
            (setB, viewsB) = EmbossSetInput.ReadAll(setB, proteinByFirst == true ? EmbossSequenceType.StopProtein : EmbossSequenceType.Nucleotide);
        }
        else
        {
            setA = EmbossSequenceReader.Clean(setA);
            setB = EmbossSequenceReader.Clean(setB);
        }

        RequireEmbossRange(minScore, -10, 100, nameof(minScore), "-minscore", "needleall");
        RequireEmbossRange(gapOpen, 1, 100, nameof(gapOpen), "-gapopen", "needleall");
        RequireEmbossRange(gapExtend, 0, 10, nameof(gapExtend), "-gapextend", "needleall");
        RequireEmbossRange(endOpen, 1, 100, nameof(endOpen), "-endopen", "needleall");
        RequireEmbossRange(endExtend, 0, 10, nameof(endExtend), "-endextend", "needleall");
        matrix ??= EmbossSetInput.DefaultMatrix(proteinByFirst, setA);
        float minscore = (float)minScore;
        var alignments = new List<EmbossPairAlignment>();
        var errors = new StringBuilder();
        for (int bi = 0; bi < setB.Count; bi++)
        {
            var b = setB[bi];
            for (int ai = 0; ai < setA.Count; ai++)
            {
                var a = setA[ai];
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair (the pair itself polls per column)
                var r = NeedleAlign(a.Residues, b.Residues, matrix, gapOpen, gapExtend, endWeight, endOpen, endExtend);
                float score = (float)r.ExactScore;
                if (score > minscore)
                {
                    var pair = viewsA is null ? EmbossPairAlignment.From(r, a.Name, a.Residues.Length, b.Name, b.Residues.Length) : EmbossPairAlignment.From(r, a.Name, viewsA[ai], b.Name, viewsB![bi]);
                    pair = pair with
                    {
                        First = pair.First with { Description = a.Description, Usa = a.Name },
                        Second = pair.Second with { Description = b.Description, Usa = b.Name },
                        ExtraSubHeader = brief ? "" : SimilaritySubHeader(r.AlignedSequence1, r.AlignedSequence2, a.Residues.Length, b.Residues.Length, matrix),
                    };
                    alignments.Add(pair);
                }
                else
                {
                    errors.Append(CultureInfo.InvariantCulture, $"Alignment score ({CFormat.Fixed(score, 1)}) is less than minimum score({CFormat.Fixed(minscore, 1)}) for sequences {a.Name} vs {b.Name}\n");
                }
            }
        }

        return new EmbossProgramResult(alignments, errors.ToString(), MatrixFileName(matrix));
    }

    /// <summary>
    /// The sub-header lines of <c>needle -nobrief</c> / <c>needleall -nobrief</c> (<c>embAlignCalcSimilarity</c>): percent identity and similarity of the aligned columns (<c>Longest_*</c>: relative to the columns
    /// with no gap) and relative to the longer sequence (<c>Shortest_*</c>), in single precision.
    /// </summary>
    /// <param name="row1">First aligned row.</param>
    /// <param name="row2">Second aligned row.</param>
    /// <param name="length1">Length of the first sequence.</param>
    /// <param name="length2">Length of the second sequence.</param>
    /// <param name="matrix">The scoring matrix (a column is similar when its score is positive).</param>
    public static string SimilaritySubHeader(string row1, string row2, int length1, int length2, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(row1);
        ArgumentNullException.ThrowIfNull(row2);
        ArgumentNullException.ThrowIfNull(matrix);
        if (row1.Length != row2.Length)
            throw new ArgumentException("The aligned rows must have the same length.", nameof(row2));
        var table = new EmbossScoreTable(matrix);
        string fm = row1.Replace('-', '.').ToUpperInvariant();
        string fn = row2.Replace('-', '.').ToUpperInvariant();
        float id = 0f;
        float sim = 0f;
        int gaps = 0;
        for (int i = 0; i < fm.Length; ++i)
        {
            if (fm[i] == '.' || fn[i] == '.')
            {
                ++gaps;
                continue;
            }

            float match = table.Score(fm[i], fn[i]);
            if (fm[i] == fn[i])
            {
                ++id;
                ++sim;
                continue;
            }

            if (match > 0.0)
                ++sim;
        }

        int max = length1 > length2 ? length1 : length2;
        float idx = (id / max) * 100f;
        float simx = (sim / max) * 100f;
        float scale = 100f / (fm.Length - gaps);
        id *= scale;
        sim *= scale;
        return string.Create(CultureInfo.InvariantCulture,
            $"Longest_Identity = {Percent(id)}%\nLongest_Similarity = {Percent(sim)}%\nShortest_Identity = {Percent(idx)}%\nShortest_Similarity = {Percent(simx)}%");

        // %5.2f; a pair with no gap-free column divides 0 by 0 (x86 produces the negative quiet NaN, printed "-nan")
        static string Percent(float v) => (float.IsNaN(v) ? "-nan" : CFormat.Fixed(v, 2)).PadLeft(5);
    }
}
