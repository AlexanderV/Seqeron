using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// RSAT <c>oligo-analysis -return mseq,proba</c> with degenerate words (<c>-oneN</c> / <c>-onedeg</c>): matching-sequence
/// statistics of words with one N or one IUPAC code. Ported from rsa-tools/rsat-code master 10043f2
/// <c>perl-scripts/oligo-analysis</c> v1.169 (<c>CountOligos</c> mseq block, <c>Degenerate</c>, <c>CalcExpected</c>,
/// <c>CalcProba</c>) and <c>perl-scripts/lib/RSA.disco.lib</c> <c>NbPossibleOligos</c>.
/// </summary>
public static partial class MotifFinder
{
    /// <summary>
    /// <see cref="FindSharedMotifs(IEnumerable{DnaSequence}, int, int, OligoBackgroundModel, OligoStrandMode, double)"/> for
    /// degenerate words (RSAT <c>-oneN</c> / <c>-onedeg</c> with <c>-return mseq,proba</c>): every word of length k is replaced
    /// by the words with one position replaced by each IUPAC code containing its residue there (N, or R Y W S M K H B V D N),
    /// and each degenerate word D (or pair D|D' with both strands) with mseq ≥ <paramref name="minSequences"/> is scored.
    /// </summary>
    /// <remarks>
    /// <para>mseq(D) is the number of sequences containing at least one word matching D (with both strands: matching D or D'),
    /// i.e. the union of the matching sequences of the words of D — RSAT's definition of mseq ("number of matching
    /// sequences"; <c>ms_freq</c> = "proportion of matching sequences (sequences with at least one occurrence)").
    /// exp_freq(D) = Σ exp_freq(w) over the words w matching D (per strand, then the pair sum; RSAT <c>CalcExpected</c> with
    /// the code probability = the sum of its residues' probabilities), optionally corrected by the pseudo-frequency
    /// ψ / NPO; P₁, exp_ms, ms_P, ms_E and ms_sig as for plain words, with NPO = k·|codes|·4^(k−1) (one strand) or
    /// NPO − (NPO − P)/2 with P = 4^(k/2) for even k, 0 for odd k (both strands; RSAT <c>NbPossibleOligos</c>).</para>
    /// <para>Deliberate difference from the RSAT code: RSAT 1.169 <c>Degenerate</c> sums the per-word counts
    /// (<c>$deg_mseq{$deg} += $patterns{$pattern_seq}->{mseq}</c>), which counts a sequence once per matching word and
    /// can exceed the number of sequences (3 sequences ACGTACGGATCC / ATGCATGAAC / ACGATGTT, k = 3, <c>-onedeg</c>:
    /// mseq = 4 &gt; 3, and RSAT stops with "Successes (4) cannot be higher than trials (3)"); as for occurrences the RSAT
    /// path returns nothing at all (it reads the never-filled global <c>%IUPAC</c>). The union, which is the documented
    /// quantity, is used here.</para>
    /// Results are in order of first occurrence (sequence, then first word, position, code).
    /// </remarks>
    /// <param name="sequences">Input DNA sequences (no null elements).</param>
    /// <param name="k">Word length (≥ 1).</param>
    /// <param name="minSequences">Matching-sequence quorum (≥ 1).</param>
    /// <param name="background">Background model.</param>
    /// <param name="strands">Single strand or reverse-complement pairs.</param>
    /// <param name="pseudoFrequency">RSAT <c>-pseudo</c> ψ in [0, 1].</param>
    /// <param name="degeneracy">Degenerate-word mode; <see cref="OligoDegeneracy.None"/> = the plain-word overload exactly.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or <paramref name="background"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1, <paramref name="minSequences"/> &lt; 1, ψ outside [0, 1], an unknown degeneracy, or an invalid Markov order.</exception>
    /// <exception cref="ArgumentException">A null sequence, or the background assigns probability 0 to an observed degenerate word.</exception>
    public static SharedMotifAnalysisResult FindSharedMotifs(
        IEnumerable<DnaSequence> sequences,
        int k,
        int minSequences,
        OligoBackgroundModel background,
        OligoStrandMode strands,
        double pseudoFrequency,
        OligoDegeneracy degeneracy)
    {
        if (!Enum.IsDefined(degeneracy))
            throw new ArgumentOutOfRangeException(nameof(degeneracy), degeneracy, "Unknown degeneracy mode.");
        if (degeneracy == OligoDegeneracy.None)
            return FindSharedMotifs(sequences, k, minSequences, background, strands, pseudoFrequency);

        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(background);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(minSequences, 1);
        ValidatePseudoFrequency(pseudoFrequency, nameof(pseudoFrequency));
        background.ValidateFor(k, nameof(background));

        bool both = strands == OligoStrandMode.Both;
        List<string> seqs = SequenceStrings(sequences);

        char[] codes = degeneracy == OligoDegeneracy.OneN ? OneNCodes : OneDegenerateCodes;
        long possiblePositions = 0;
        var rcCache = new Dictionary<string, string>(StringComparer.Ordinal);
        var matching = new Dictionary<string, (string Rc, List<int> Indices)>(StringComparer.Ordinal);
        var buffer = new char[k];
        for (int i = 0; i < seqs.Count; i++)
        {
            string s = seqs[i];
            if (s.Length < k) continue;
            possiblePositions += s.Length - k + 1;
            foreach (string word in s.AsSpan().CountKmersSpan(k).Keys)
            {
                for (int l = 0; l < k; l++)
                {
                    foreach (char code in codes)
                    {
                        if (!IupacHelper.MatchesIupac(word[l], code)) continue;
                        word.CopyTo(0, buffer, 0, k);
                        buffer[l] = code;
                        AddMatchingSequence(matching, new string(buffer), i, both, rcCache);
                    }
                }
            }
        }

        double npo = RsatPossibleOligos(k, both, codes.Length);
        double logNpo = Math.Log(npo);
        if (!double.IsFinite(npo))
        {
            // k·|codes|·4^(k−1) beyond the double range; the palindrome term is negligible against it.
            logNpo = Math.Log(k) + Math.Log(codes.Length) + (k - 1) * Ln4;
            if (both) logNpo -= Math.Log(2.0);
        }
        var motifs = new List<SignificantSharedMotif>();
        var shared = matching.Where(e => e.Value.Indices.Count >= minSequences).ToList();
        if (shared.Count > 0)
        {
            Func<string, double> logWord = background.CreateLogProbability(seqs, k, both);
            int sequenceCount = seqs.Count;
            double positionsPerSequence = (double)possiblePositions / sequenceCount;
            foreach (var (word, (rc, indices)) in shared)
            {
                var pattern = new GroupedOligo(word, rc, 0, 0, 0, new List<OligoOccurrence>());
                double logP = PatternLogFrequency(pattern, logWord, codes, pseudoFrequency, logNpo);
                if (double.IsNegativeInfinity(logP))
                    throw new ArgumentException($"The background model assigns probability 0 to the observed word {word}.");
                motifs.Add(SharedMotifStatistics(word, both ? rc : null, indices, sequenceCount, positionsPerSequence, logP, logNpo));
            }
        }

        return new SharedMotifAnalysisResult(
            Motifs: motifs.AsReadOnly(),
            OligoLength: k,
            Strands: strands,
            SequenceCount: seqs.Count,
            PossiblePositions: possiblePositions,
            PossibleOligos: npo)
        {
            Degeneracy = degeneracy,
        };
    }
}
