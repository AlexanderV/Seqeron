using System.Collections;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The alignments <see cref="PairwiseAligner.Align"/> found (Biopython's <c>PairwiseAlignments</c>): all alignments of the optimal score, each as a <see cref="BioCoordAlignment"/> with the <c>score</c> attribute.</summary>
public sealed class PairwiseAlignments : IReadOnlyList<BioCoordAlignment>
{
    private readonly IReadOnlyList<BioCoordAlignment> _alignments;

    internal PairwiseAlignments(double score, IReadOnlyList<BioCoordAlignment> alignments)
    {
        Score = score;
        _alignments = alignments;
    }

    /// <summary>The optimal score (<c>alignments.score</c>).</summary>
    public double Score { get; }

    /// <inheritdoc/>
    public int Count => _alignments.Count;

    /// <inheritdoc/>
    public BioCoordAlignment this[int index] => _alignments[index];

    /// <inheritdoc/>
    public IEnumerator<BioCoordAlignment> GetEnumerator() => _alignments.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}

public sealed partial class PairwiseAligner
{
    private const int MaxAlignments = 1000;

    private static void RequireStrand(char strand)
    {
        if (strand is not ('+' or '-'))
            throw new ArgumentException("strand must be '+' or '-'.", nameof(strand));
    }


    // the alphabet-wide scalar matrix of match / mismatch / wildcard scoring (the aligners score letters through SubstitutionMatrix when the scores are not integers or a wildcard letter exists)
    private SubstitutionMatrix ScalarMatrix(string a, string b)
    {
        var letters = new SortedSet<char>(a.Concat(b));
        if (Wildcard is { } w)
            letters.Add(char.ToUpperInvariant(w));
        string alphabet = new([.. letters]);
        var scores = new decimal[alphabet.Length, alphabet.Length];
        for (int i = 0; i < alphabet.Length; i++)
        {
            for (int j = 0; j < alphabet.Length; j++)
            {
                bool wild = Wildcard is { } wc && (alphabet[i] == char.ToUpperInvariant(wc) || alphabet[j] == char.ToUpperInvariant(wc));
                double score = i == j ? MatchScore : MismatchScore;
                scores[i, j] = wild ? 0m : (decimal)score;
            }
        }

        return new SubstitutionMatrix("scalar", alphabet, scores);
    }

    private GapScoreProfile Profile(char strand)
    {
        GapScore Slot(int s) => GapScore.FromBiopython(_open[s], _extend[s]);
        // strand '-': the left and right end gaps change places (Bio.Align swaps them for the reverse strand)
        int left = strand == '+' ? 1 : 2, right = 3 - left;
        return new GapScoreProfile(Interior1: Slot(3), Start1: Slot(3 + left), End1: Slot(3 + right), Interior2: Slot(0), Start2: Slot(left), End2: Slot(right));
    }

    private static AlignmentResult DpSingle(string a, string b, ScoringMatrix scoring, GapScoreProfile profile, bool local, bool affine)
    {
        if (local)
            return affine ? SequenceAligner.LocalAlignAffine(a, b, scoring, profile) : SequenceAligner.LocalAlign(a, b, scoring, profile);
        return affine ? SequenceAligner.GlobalAlignAffine(a, b, scoring, profile) : SequenceAligner.GlobalAlign(a, b, scoring, profile);
    }

    private static AlignmentResult DpSingle(string a, string b, MatrixScoring scoring, GapScoreProfile profile, bool local, bool affine)
    {
        if (local)
            return affine ? SequenceAligner.LocalAlignAffine(a, b, scoring, profile) : SequenceAligner.LocalAlign(a, b, scoring, profile);
        return affine ? SequenceAligner.GlobalAlignAffine(a, b, scoring, profile) : SequenceAligner.GlobalAlign(a, b, scoring, profile);
    }

    private static IReadOnlyList<AlignmentResult> DpAll(string a, string b, ScoringMatrix scoring, GapScoreProfile profile, bool local, bool affine)
    {
        if (local)
            return affine ? SequenceAligner.LocalAlignAffineCoOptimal(a, b, scoring, profile, MaxAlignments) : SequenceAligner.LocalAlignCoOptimal(a, b, scoring, profile, MaxAlignments);
        return affine ? SequenceAligner.GlobalAlignAffineCoOptimal(a, b, scoring, profile, MaxAlignments) : SequenceAligner.GlobalAlignCoOptimal(a, b, scoring, profile, MaxAlignments);
    }

    private static IReadOnlyList<AlignmentResult> DpAll(string a, string b, MatrixScoring scoring, GapScoreProfile profile, bool local, bool affine)
    {
        if (local)
            return affine ? SequenceAligner.LocalAlignAffineCoOptimal(a, b, scoring, profile, MaxAlignments) : SequenceAligner.LocalAlignCoOptimal(a, b, scoring, profile, MaxAlignments);
        return affine ? SequenceAligner.GlobalAlignAffineCoOptimal(a, b, scoring, profile, MaxAlignments) : SequenceAligner.GlobalAlignCoOptimal(a, b, scoring, profile, MaxAlignments);
    }

    private (double Score, AlignmentResult? Single, IReadOnlyList<AlignmentResult> All) RunDp(string a, string b, char strand, bool all)
    {
        string ua = a.ToUpperInvariant(), ub = b.ToUpperInvariant();
        var profile = Profile(strand);
        bool local = Mode == PairwiseAlignerMode.Local;
        bool integers = SubstitutionMatrix is null && Wildcard is null && MatchScore == Math.Round(MatchScore) && MismatchScore == Math.Round(MismatchScore)
            && Enumerable.Range(0, 6).All(i => _open[i] == Math.Round(_open[i]) && _extend[i] == Math.Round(_extend[i]));
        bool affine = IsAffine;
        AlignmentResult single;
        IReadOnlyList<AlignmentResult> list = [];
        if (integers)
        {
            var scoring = new ScoringMatrix((int)MatchScore, (int)MismatchScore, (int)(_open[3] - _extend[3]), (int)_extend[3]);
            single = DpSingle(ua, ub, scoring, profile, local, affine);
            if (all)
                list = DpAll(ua, ub, scoring, profile, local, affine);
        }
        else
        {
            var matrix = SubstitutionMatrix ?? ScalarMatrix(ua, ub);
            var scoring = new MatrixScoring(matrix, _open[3] - _extend[3], _extend[3]);
            single = DpSingle(ua, ub, scoring, profile, local, affine);
            if (all)
                list = DpAll(ua, ub, scoring, profile, local, affine);
        }

        return (single.ExactScore, single, list.Count > 0 ? list : [single]);
    }

    private (double Score, List<(char A, char B)> Columns) RunFogsaaScored(string a, string b, char strand)
    {
        string ua = a.ToUpperInvariant(), ub = b.ToUpperInvariant();
        if (SubstitutionMatrix is { } matrix)
        {
            double match = double.NegativeInfinity, mismatch = double.PositiveInfinity;
            decimal first = 0;
            bool any = false;
            foreach (char x in matrix.Alphabet)
            {
                foreach (char y in matrix.Alphabet)
                {
                    decimal v = matrix.GetValue(x, y);
                    if (!any)
                    {
                        first = v;
                        match = mismatch = (double)first;
                        any = true;
                    }

                    if ((double)v > match)
                        match = (double)v;
                    else if ((double)v < mismatch)
                        mismatch = (double)v;
                }
            }

            return RunFogsaa(ua, ub, strand, (x, y) => (double)matrix.GetValue(x, y), match, mismatch);
        }

        char? wildcard = Wildcard is { } w ? char.ToUpperInvariant(w) : null;
        double m = MatchScore, mm = MismatchScore;
        double Compare(char x, char y)
        {
            if (wildcard is { } c && (x == c || y == c))
                return 0;
            return x == y ? m : mm;
        }

        return RunFogsaa(ua, ub, strand, Compare, m, mm);
    }

    /// <summary>
    /// <c>aligner.score(sequenceA, sequenceB, strand)</c>: the optimal alignment score. Letters are compared case-insensitively (the sequences are upper-cased); strand <c>-</c> aligns A with the reverse complement of B.
    /// </summary>
    public double Score(string sequenceA, string sequenceB, char strand = '+')
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        RequireStrand(strand);
        string b = strand == '-' ? AlignmentCoordinateFormats.ReverseComplement(sequenceB) : sequenceB;
        if (Mode == PairwiseAlignerMode.Fogsaa)
            return RunFogsaaScored(sequenceA, b, strand).Score;
        if (UsesWaterman)
            return WsbScore(sequenceA.ToUpperInvariant(), b.ToUpperInvariant(), strand);
        if (UsesSmithWaterman)
            return SmithWatermanScore(sequenceA.ToUpperInvariant(), b.ToUpperInvariant());
        if (UsesGotohLocal)
            return GotohLocalScore(sequenceA.ToUpperInvariant(), b.ToUpperInvariant());
        return RunDp(sequenceA, b, strand, all: false).Score;
    }

    /// <summary>
    /// <c>aligner.align(sequenceA, sequenceB, strand)</c>: every alignment of the optimal score (one alignment in <c>fogsaa</c> mode) as coordinates of the target (A) and the query (B; on strand <c>-</c> the query coordinates decrease). The enumeration is the
    /// library's co-optimal enumeration (the same SET of alignments as Biopython's; the order is not Biopython's: compare with <see cref="AlignmentOperations.Compare"/>). A local alignment without a positive score gives no alignment and the score 0.
    /// </summary>
    public PairwiseAlignments Align(string sequenceA, string sequenceB, char strand = '+')
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        RequireStrand(strand);
        string b = strand == '-' ? AlignmentCoordinateFormats.ReverseComplement(sequenceB) : sequenceB;
        double score;
        List<AlignmentResult> results;
        if (Mode == PairwiseAlignerMode.Fogsaa)
        {
            var (s, columns) = RunFogsaaScored(sequenceA, b, strand);
            score = s;
            string rowA = new([.. columns.Select(c => c.A)]), rowB = new([.. columns.Select(c => c.B)]);
            results = [new AlignmentResult(rowA, rowB, (int)Math.Round(s, MidpointRounding.AwayFromZero), AlignmentType.Global, 0, 0, sequenceA.Length - 1, b.Length - 1) { ExactScore = s }];
        }
        else if (UsesWaterman || UsesGotohLocal)
        {
            string wa = sequenceA.ToUpperInvariant(), wb = b.ToUpperInvariant();
            (double, List<AlignmentResult>) run;
            if (UsesWaterman)
                run = WsbAlign(wa, wb, strand);
            else if (UsesSmithWaterman)
                run = SmithWatermanAlign(wa, wb);
            else
                run = GotohLocalAlign(wa, wb);
            var (s, found) = run;
            score = s;
            results = [.. found.Where(r => r.AlignedSequence1.Length > 0)];
            if (Mode == PairwiseAlignerMode.Local && results.Count == 0)
                return new PairwiseAlignments(score, []);
        }
        else
        {
            var (s, _, all) = RunDp(sequenceA, b, strand, all: true);
            score = s;
            results = [.. all.Where(r => r.AlignedSequence1.Length > 0)];
            if (Mode == PairwiseAlignerMode.Local && results.Count == 0)
                return new PairwiseAlignments(0, []);
        }

        var alignments = new List<BioCoordAlignment>();
        string ua = sequenceA.ToUpperInvariant(), ub = b.ToUpperInvariant();
        foreach (var r in results)
        {
            var bridged = r.ToBioCoordAlignment("target", ua, "query", ub);
            var target = new BioCoordSequence { Id = null, Length = ua.Length, Parts = bridged.Sequences[0].Parts };
            var original = strand == '-' ? AlignmentCoordinateFormats.ReverseComplement(ub) : ub;
            var query = new BioCoordSequence { Id = null, Length = original.Length, Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, original)], original.Length) };
            var coordinates = bridged.Coordinates.Select(row => (long[])row.Clone()).ToArray();
            if (strand == '-')
                for (int k = 0; k < coordinates[1].Length; k++)
                    coordinates[1][k] = original.Length - coordinates[1][k];
            alignments.Add(new BioCoordAlignment { Sequences = [target, query], Coordinates = coordinates, Attributes = bridged.Attributes });
        }

        return new PairwiseAlignments(score, alignments);
    }
}
