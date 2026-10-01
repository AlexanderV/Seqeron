namespace Seqeron.Genomics.Analysis;

/// <summary>
/// RSAT <c>oligo-analysis -seqtype dna|prot|other</c> for word statistics on raw sequence strings (DNA with undefined
/// residues, protein, or any text). Ported from rsa-tools/rsat-code master 10043f2 <c>perl-scripts/oligo-analysis</c> v1.169:
/// <c>ReadArguments</c> (<c>-seqtype</c>), the alphabet / accepted-residue block, <c>CountOligos</c> (residue counts,
/// discarded patterns and residues), <c>sub alphabet</c>, <c>CalcAlphabet</c>, <c>CalcExpected</c> (equiprobable branch,
/// <c>-pseudo</c>), <c>perl-scripts/lib/RSA.disco.lib</c> <c>NbPossibleOligos</c>, <c>perl-scripts/lib/RSA.seq.lib</c>
/// <c>FoldSequence</c> / <c>OverlapCoeff</c>.
/// </summary>
public static partial class MotifFinder
{
    /// <summary>
    /// RSAT <c>oligo-analysis</c> word statistics of raw sequence strings with the sequence type of
    /// <see cref="OligoAnalysisOptions.SequenceType"/> (RSAT <c>-seqtype</c>); all other options and every output column are
    /// as in <see cref="AnalyzeOligos(IEnumerable{DnaSequence}, int, OligoAnalysisOptions?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>Input handling (RSAT <c>FoldSequence</c> and <c>CountOligos</c>): white space (space, tab, CR, LF, FF, VT) is
    /// removed; ASCII letters are case-folded — to upper case for DNA and protein, to lower case for other (RSAT reports
    /// lower-case words); other characters are kept as they are (Perl <c>lc</c> on byte strings folds only A–Z).</para>
    /// <list type="bullet">
    /// <item><see cref="OligoSequenceType.Dna"/>: alphabet A C G T. A window containing any other residue (N, IUPAC code, gap,
    /// …) is discarded: it is neither an occurrence nor an overlap, and nb_pos = n counts only the remaining windows
    /// (RSAT <c>$nb_possible_pos -= $discarded_occurrences</c>); residues outside the alphabet are not counted for the input
    /// Bernoulli model. Without such residues the result equals the <see cref="DnaSequence"/> overload bit for bit.
    /// Every DNA option applies.</item>
    /// <item><see cref="OligoSequenceType.Protein"/>: alphabet of the 20 amino acids A C D E F G H I K L M N P Q R S T V W Y
    /// (RSAT <c>@protein_alphabet</c>; X, B, Z, U, O, * and any other character discard their windows), |A| = 20.</item>
    /// <item><see cref="OligoSequenceType.Other"/>: every residue is valid; the alphabet is the set of distinct residues of the
    /// sequences of length ≥ k (RSAT <c>sub alphabet</c>: <c>sort keys %residue_occ</c>).</item>
    /// </list>
    /// <para>For protein and other sequences (RSAT sets <c>$sum_rc = 0</c> with <c>-seqtype prot|other</c>): single strand only;
    /// NPO = |A|^k (no reverse-complement or degenerate correction, RSAT <c>NbPossibleOligos</c> applies those to DNA only),
    /// used by <c>-pseudo</c>; equiprobable exp_freq = 1 / |A|^k; the input Bernoulli residue probabilities are
    /// count(r) / Σ over the alphabet residues of the sequences of length ≥ k (RSAT <c>CalcAlphabet</c>); input Markov
    /// chains and the lexicon use the counted words exactly as for DNA. The overlap coefficient uses the input Bernoulli
    /// probabilities (input Bernoulli, lexicon) or 1 / |A| (equiprobable). Deliberate difference: for Markov backgrounds
    /// RSAT <c>OverlapCoeff</c> falls back to "equiprobable nucleotides" (a, c, g, t = ¼, every other residue 0, so a
    /// protein word such as LLL has overlap coefficient 1); the equiprobable residues of the alphabet, 1 / |A|, are used here.</para>
    /// <para>Not available for protein / other (ArgumentException): both strands (no reverse complement), degenerate IUPAC
    /// words (nucleotide codes; RSAT's degenerate NPO is DNA-only), calibration tables (RSAT <c>calibrate-oligos</c> has no
    /// sequence type and <c>ReadCalibration</c> applies DNA reverse-complement inference), and the DNA-specific background
    /// models <see cref="OligoBackgroundModel.Bernoulli"/> and <see cref="OligoBackgroundModel.MarkovFromOligoFrequencies"/>.</para>
    /// <para>Positions refer to the sequences after white-space removal.</para>
    /// </remarks>
    /// <param name="sequences">Input sequences (no null elements).</param>
    /// <param name="k">Word length (≥ 1).</param>
    /// <param name="options">Options; <see cref="OligoAnalysisOptions.SequenceType"/> selects the alphabet (default DNA).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1 or an invalid option value.</exception>
    /// <exception cref="ArgumentException">A null sequence, or an option not available for the sequence type.</exception>
    public static OligoAnalysisReport AnalyzeOligoStrings(IEnumerable<string> sequences, int k, OligoAnalysisOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        options ??= new OligoAnalysisOptions();
        ValidateOligoOptions(k, options);

        var type = options.SequenceType;
        var seqs = new List<string>();
        foreach (string? raw in sequences)
        {
            if (raw is null)
                throw new ArgumentException($"Sequence at index {seqs.Count} is null.", nameof(sequences));
            seqs.Add(FoldResidueString(raw, type));
        }

        var alphabet = OligoResidueAlphabet.For(type, seqs, k);
        if (type == OligoSequenceType.Dna && !seqs.Any(s => s.Length >= k && alphabet.InvalidPrefixCounts(s)[^1] > 0))
            alphabet = null; // pure ACGT: the DnaSequence code path
        return AnalyzeOligosCore(seqs, k, options, alphabet);
    }

    /// <summary>Protein / other sequence types: the options RSAT supports for them (see <see cref="AnalyzeOligoStrings"/>).</summary>
    private static void ValidateResidueAlphabetOptions(OligoAnalysisOptions options)
    {
        string type = options.SequenceType.ToString().ToLowerInvariant();
        if (options.Strands == OligoStrandMode.Both)
            throw new ArgumentException(
                $"Both strands need reverse complements; RSAT -seqtype {type} counts a single strand.", nameof(options));
        if (options.Degeneracy != OligoDegeneracy.None)
            throw new ArgumentException(
                $"Degenerate words use IUPAC nucleotide codes and are available for DNA only (sequence type {type}).", nameof(options));
        if (options.Calibration is not null)
            throw new ArgumentException(
                $"Calibration tables are DNA word tables (RSAT calibrate-oligos); not available for sequence type {type}.", nameof(options));
        if (!options.Background.IsAlphabetGeneric)
            throw new ArgumentException(
                $"Sequence type {type} supports the equiprobable, input Bernoulli, input Markov and lexicon backgrounds only.", nameof(options));
    }

    // RSAT FoldSequence (s/\s+//g) + lc: white space removed, ASCII letters folded (upper case for DNA / protein).
    private static string FoldResidueString(string raw, OligoSequenceType type)
    {
        var sb = new System.Text.StringBuilder(raw.Length);
        bool lower = type == OligoSequenceType.Other;
        foreach (char c in raw)
        {
            if (c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v') continue;
            if (lower && c is >= 'A' and <= 'Z') sb.Append((char)(c + 32));
            else if (!lower && c is >= 'a' and <= 'z') sb.Append((char)(c - 32));
            else sb.Append(c);
        }
        return sb.ToString();
    }
}

/// <summary>Sequence type of RSAT <c>oligo-analysis</c> (<c>-seqtype</c>).</summary>
public enum OligoSequenceType
{
    /// <summary>DNA (default): A, C, G, T; windows with other residues are discarded.</summary>
    Dna,

    /// <summary>Protein: the 20 amino acids A C D E F G H I K L M N P Q R S T V W Y; single strand.</summary>
    Protein,

    /// <summary>Any text: every residue is valid, the alphabet is the set of observed residues; single strand.</summary>
    Other,
}

/// <summary>Residue alphabet of an RSAT <c>-seqtype</c> analysis.</summary>
internal sealed class OligoResidueAlphabet
{
    private const string ProteinResidues = "ACDEFGHIKLMNPQRSTVWY";
    private readonly HashSet<char>? _valid; // null = every residue valid (other)

    private OligoResidueAlphabet(OligoSequenceType type, IReadOnlyList<char> residues, HashSet<char>? valid)
    {
        Type = type;
        Residues = residues;
        _valid = valid;
    }

    public OligoSequenceType Type { get; }

    public bool IsDna => Type == OligoSequenceType.Dna;

    /// <summary>Alphabet residues in ordinal order.</summary>
    public IReadOnlyList<char> Residues { get; }

    /// <summary>RSAT <c>alphabet_size</c>.</summary>
    public int Size => Residues.Count;

    /// <summary>The alphabet of <paramref name="type"/>; for other, the distinct residues of the sequences of length ≥ k.</summary>
    public static OligoResidueAlphabet For(OligoSequenceType type, IReadOnlyList<string> sequences, int k)
    {
        switch (type)
        {
            case OligoSequenceType.Dna:
                return new OligoResidueAlphabet(type, MotifFinder.AcgtBases, new HashSet<char>(MotifFinder.AcgtBases));
            case OligoSequenceType.Protein:
                return new OligoResidueAlphabet(type, ProteinResidues.ToCharArray(), new HashSet<char>(ProteinResidues));
            default:
                var observed = new SortedSet<char>(Comparer<char>.Create((a, b) => a.CompareTo(b)));
                foreach (string s in sequences)
                {
                    if (s.Length < k) continue;
                    foreach (char c in s) observed.Add(c);
                }
                return new OligoResidueAlphabet(type, observed.ToArray(), null);
        }
    }

    public bool IsValid(char c) => _valid is null || _valid.Contains(c);

    /// <summary>p[i] = number of residues outside the alphabet in s[0..i); a window [pos, pos + k) is valid iff p[pos + k] = p[pos].</summary>
    public int[] InvalidPrefixCounts(string s)
    {
        var p = new int[s.Length + 1];
        for (int i = 0; i < s.Length; i++)
            p[i + 1] = p[i] + (IsValid(s[i]) ? 0 : 1);
        return p;
    }

    /// <summary>RSAT <c>NbPossibleOligos</c> for protein / other: |A|^k.</summary>
    public double PossibleOligos(int k) => Math.Pow(Size, k);

    /// <summary>ln |A|^k, finite for every k.</summary>
    public double LogPossibleOligos(int k)
    {
        double npo = PossibleOligos(k);
        return double.IsFinite(npo) ? Math.Log(npo) : k * Math.Log(Size);
    }

    /// <summary>RSAT <c>CalcAlphabet</c>: count(r) / Σ count over the alphabet residues of the sequences of length ≥ k.</summary>
    public Dictionary<char, double> InputResidueProbabilities(IReadOnlyList<string> sequences, int k)
    {
        var counts = new Dictionary<char, long>();
        long total = 0;
        foreach (string s in sequences)
        {
            if (s.Length < k) continue;
            foreach (char c in s)
            {
                if (!IsValid(c)) continue;
                counts[c] = counts.GetValueOrDefault(c) + 1;
                total++;
            }
        }
        return counts.ToDictionary(e => e.Key, e => (double)e.Value / total);
    }

    /// <summary>
    /// Residue probabilities of the overlap coefficient (RSAT <c>%residue_proba</c> in <c>OverlapCoeff</c>): input Bernoulli
    /// (input Bernoulli and lexicon backgrounds), otherwise the equiprobable 1 / |A| (equiprobable; Markov — see
    /// <see cref="MotifFinder.AnalyzeOligoStrings"/> for RSAT's DNA-only fallback).
    /// </summary>
    public Dictionary<char, double> OverlapResidueProbabilities(OligoBackgroundModel background, IReadOnlyList<string> sequences, int k)
    {
        if (!background.IsEquiprobable && !background.IsMarkov)
            return InputResidueProbabilities(sequences, k);
        double q = 1.0 / Size;
        return Residues.ToDictionary(c => c, _ => q);
    }
}
