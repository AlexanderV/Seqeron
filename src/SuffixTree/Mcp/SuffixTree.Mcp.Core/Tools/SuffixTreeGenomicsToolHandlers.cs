namespace SuffixTree.Mcp.Core.Tools;

internal static class SuffixTreeGenomicsToolHandlers
{
    public static FindLongestRepeatResult FindLongestRepeat(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var dna = new global::Seqeron.Genomics.Core.DnaSequence(sequence);
        var result = GenomicAnalyzer.FindLongestRepeat(dna);

        return new FindLongestRepeatResult(result.Sequence, result.Positions.ToArray(), result.Length);
    }

    public static FindLongestCommonRegionResult FindLongestCommonRegion(string sequence1, string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var dna1 = new global::Seqeron.Genomics.Core.DnaSequence(sequence1);
        var dna2 = new global::Seqeron.Genomics.Core.DnaSequence(sequence2);
        var result = GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);

        return new FindLongestCommonRegionResult(result.Sequence, result.PositionInFirst, result.PositionInSecond, result.Length);
    }

    public static CalculateSimilarityResult CalculateSimilarity(string sequence1, string sequence2, int kmerSize)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var dna1 = new global::Seqeron.Genomics.Core.DnaSequence(sequence1);
        var dna2 = new global::Seqeron.Genomics.Core.DnaSequence(sequence2);
        var similarity = GenomicAnalyzer.CalculateSimilarity(dna1, dna2, kmerSize);

        return new CalculateSimilarityResult(similarity);
    }

    public static HammingDistanceResult HammingDistance(string sequence1, string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));
        if (sequence1.Length != sequence2.Length)
            throw new ArgumentException("Sequences must have equal length for Hamming distance");

        var distance = ApproximateMatcher.HammingDistance(sequence1, sequence2);
        return new HammingDistanceResult(distance);
    }

    public static EditDistanceResult EditDistance(string sequence1, string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var distance = ApproximateMatcher.EditDistance(sequence1, sequence2);
        return new EditDistanceResult(distance);
    }

    public static CountApproximateOccurrencesResult CountApproximateOccurrences(string sequence, string pattern, int maxMismatches)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        if (maxMismatches < 0)
            throw new ArgumentException("MaxMismatches cannot be negative", nameof(maxMismatches));

        var count = ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches);
        return new CountApproximateOccurrencesResult(count);
    }
}
