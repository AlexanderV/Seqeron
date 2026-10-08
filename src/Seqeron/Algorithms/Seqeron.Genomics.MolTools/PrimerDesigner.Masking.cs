using System.Buffers.Binary;

namespace Seqeron.Genomics.MolTools;

// Primer3 template masking (PRIMER_MASK_TEMPLATE, PRIMER_MASK_FAILURE_RATE, PRIMER_MASK_5P/3P_DIRECTION,
// PRIMER_MASK_KMERLIST_PATH/PREFIX, PRIMER_WT_MASK_FAILURE_RATE) — primer3-py 2.3.1 masker.c / masker.h
// (Kõressaar et al. 2018, "Primer3_masker: integrating masking of template sequence with primer design software",
// Bioinformatics 34:1937) as libprimer3.c drives it: create_default_formula_parameters (the two GenomeTester4 k-mer
// lists <prefix>_11.list / <prefix>_16.list, coefficients 0.1772 / 0.239, intercept −4.336), calculate_scores,
// mask_oligo_region, read_and_mask_sequence (both strands separately, soft masking), calc_and_check_oligo_features
// (failure_rate, is_lowercase_masked) and p_obj_fn (weights.failure_rate).
public static partial class PrimerDesigner
{
    /// <summary>Primer3 default PRIMER_MASK_FAILURE_RATE (0.1, <c>pr_set_default_global_args_1</c>).</summary>
    public const double Primer3MaskFailureRate = 0.1;

    /// <summary>Primer3 default PRIMER_MASK_5P_DIRECTION (1 nucleotide).</summary>
    public const int Primer3MaskFivePrimeDirection = 1;

    /// <summary>Primer3 default PRIMER_MASK_3P_DIRECTION (0 nucleotides).</summary>
    public const int Primer3MaskThreePrimeDirection = 0;

    // masker.h: DEFAULT_WORD_LEN_1/_2, DEFAULT_COEF_1/_2, DEFAULT_INTERCEPT, MAX_BUFFER_SIZE.
    internal const int MaskWordLength1 = 11;
    internal const int MaskWordLength2 = 16;
    internal const double MaskCoefficient1 = 0.1772;
    internal const double MaskCoefficient2 = 0.239;
    internal const double MaskIntercept = -4.336;
    private const int MaskBufferSize = 5000;

    /// <summary>
    /// Primer3's predicted PCR failure rate of a primer (<c>calc_and_check_oligo_features</c> → masker.c
    /// <c>calculate_scores</c>; the <c>failure_rate</c> of PRIMER_WT_MASK_FAILURE_RATE): with s = 0.1772·ln n₁₁ +
    /// 0.239·ln n₁₆, where n₁₁ / n₁₆ are the genome counts of the primer's 3′-terminal 11-mer / 16-mer (the k-mer's own
    /// count if listed, otherwise its reverse complement's, otherwise 1 — Primer3's lists omit k-mers seen once),
    /// failure rate = e^(s − 4.336) / (1 + e^(s − 4.336)), and 0 when s = 0 (both k-mers absent) or the primer is shorter
    /// than 16 nt (masker window).
    /// </summary>
    /// <param name="primer">Primer sequence 5′→3′ (A/C/G/T, case-insensitive; a right primer as synthesized).</param>
    /// <param name="kmerLists">The 11-mer and 16-mer count lists.</param>
    /// <returns>The failure rate in [0, 1).</returns>
    /// <exception cref="ArgumentNullException">A null argument.</exception>
    /// <exception cref="ArgumentException">A base other than A/C/G/T among the last 16.</exception>
    public static double CalculateMaskFailureRatePrimer3(string primer, PrimerMaskingKmerLists kmerLists)
    {
        ArgumentNullException.ThrowIfNull(primer);
        ArgumentNullException.ThrowIfNull(kmerLists);
        if (primer.Length < MaskWordLength2)
            return 0.0;
        ulong word = 0;
        for (int i = primer.Length - MaskWordLength2; i < primer.Length; i++)
            word = (word << 2) | MaskNucleotideValue(primer[i], nameof(primer));
        // calc_and_check_oligo_features: op.rev = op.fwd ("not used in this calculation"); failure_rate = score_fwd.
        return kmerLists.Scores(word, word, MaskWordLength2).Forward;
    }

    /// <summary>
    /// Primer3's soft-masked template (masker.c <c>read_and_mask_sequence</c> with libprimer3's settings: both strands
    /// masked separately, soft masking). Every 16-nt window is scored as in <see cref="CalculateMaskFailureRatePrimer3"/>
    /// — on the forward strand from its last 11 / 16 bases, on the reverse strand from the reverse complement (the 11-mer
    /// at the window start); when a strand's score exceeds <paramref name="failureRate"/> (and it is non-zero) the
    /// forward copy lower-cases the window's last base and <paramref name="fivePrimeDirection"/> − 1 bases before it plus
    /// the <paramref name="threePrimeDirection"/> following bases, the reverse copy the window's first base, the
    /// <paramref name="fivePrimeDirection"/> − 1 bases after it and the <paramref name="threePrimeDirection"/> bases before
    /// it (the 3′ end of a primer on that strand and its 5′ / 3′ neighbourhood). Reproduces masker.c's 5000-character ring
    /// buffer exactly.
    /// </summary>
    /// <param name="template">Template (A/C/G/T, upper case).</param>
    /// <param name="kmerLists">The 11-mer and 16-mer count lists.</param>
    /// <param name="failureRate">PRIMER_MASK_FAILURE_RATE (default 0.1; 0 masks nothing).</param>
    /// <param name="fivePrimeDirection">PRIMER_MASK_5P_DIRECTION (default 1, ≥ 0).</param>
    /// <param name="threePrimeDirection">PRIMER_MASK_3P_DIRECTION (default 0, 0–4984).</param>
    /// <returns>The forward-strand and reverse-strand masked copies (masked bases in lower case).</returns>
    public static (string Forward, string Reverse) MaskTemplatePrimer3(
        string template, PrimerMaskingKmerLists kmerLists, double failureRate = Primer3MaskFailureRate,
        int fivePrimeDirection = Primer3MaskFivePrimeDirection, int threePrimeDirection = Primer3MaskThreePrimeDirection)
    {
        ArgumentNullException.ThrowIfNull(template);
        ArgumentNullException.ThrowIfNull(kmerLists);
        ValidateMaskSettings(failureRate, fivePrimeDirection, threePrimeDirection, nameof(failureRate));
        foreach (char c in template)
            MaskNucleotideValue(c, nameof(template));
        return new TemplateMasker(kmerLists, failureRate, fivePrimeDirection, threePrimeDirection).Mask(template);
    }

    internal static void ValidateMaskSettings(double failureRate, int fivePrimeDirection, int threePrimeDirection, string paramName)
    {
        if (!double.IsFinite(failureRate))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_MASK_FAILURE_RATE must be finite.");
        // Primer3 stores PRIMER_MASK_3P_DIRECTION in an unsigned counter and sizes the masking window with it: negative
        // values (and a window beyond masker.c's 5000-character buffer) are undefined there, so they are rejected.
        if (fivePrimeDirection < 0 || threePrimeDirection < 0 || threePrimeDirection > MaskBufferSize - MaskWordLength2)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_MASK_5P_DIRECTION must be ≥ 0 and PRIMER_MASK_3P_DIRECTION in [0, 4984].");
    }

    // masker.c get_nucl_value for A/C/G/T (U = T).
    private static ulong MaskNucleotideValue(char c, string paramName) => c switch
    {
        'A' or 'a' => 0UL,
        'C' or 'c' => 1UL,
        'G' or 'g' => 2UL,
        'T' or 't' or 'U' or 'u' => 3UL,
        _ => throw new ArgumentException($"Masking needs A/C/G/T, found '{c}'.", paramName),
    };

    // Per-design masking state: the soft-masked forward / reverse copies of the included region.
    internal sealed class MaskContext
    {
        private readonly string _forward, _reverse;
        private readonly int _offset;

        internal MaskContext(string seq, int incStart, int incEnd, PrimerPairOptions opt)
        {
            Lists = opt.MaskKmerLists!;
            _offset = incStart;
            (_forward, _reverse) = new TemplateMasker(Lists, opt.MaskFailureRate, opt.MaskFivePrimeDirection,
                opt.MaskThreePrimeDirection).Mask(seq.Substring(incStart, incEnd - incStart));
        }

        internal PrimerMaskingKmerLists Lists { get; }

        // is_lowercase_masked at the primer's 3' end: left primer → trimmed_masked_seq[last base],
        // right primer → trimmed_masked_seq_r[leftmost base].
        internal bool ThreePrimeEndMasked(int position, int length, bool isForward) => isForward
            ? char.IsLower(_forward[position + length - 1 - _offset])
            : char.IsLower(_reverse[position - _offset]);
    }

    // Literal port of masker.c read_and_mask_sequence / add_char_to_buffer / empty_buffer / mask_oligo_region for
    // mdir = both_separately, do_soft_masking = 1, abs_cutoff = 0 and an A/C/G/T input (no FASTA header, whitespace or
    // masked characters).
    private sealed class TemplateMasker(PrimerMaskingKmerLists lists, double failureRate, int m5p, int m3p)
    {
        private readonly char[] _buffer = new char[MaskBufferSize];
        private readonly bool[] _maskFwd = new bool[MaskBufferSize], _maskRev = new bool[MaskBufferSize];
        private int _ri, _wi, _ei;
        private uint _mi;

        private static int Forward(int i) => i == MaskBufferSize - 1 ? 0 : i + 1;
        private static int Back(int i) => i == 0 ? MaskBufferSize - 1 : i - 1;

        public (string Forward, string Reverse) Mask(string seq)
        {
            const int w = MaskWordLength2;
            ulong mask = (1UL << (2 * w)) - 1;
            var fwdOut = new System.Text.StringBuilder(seq.Length);
            var revOut = new System.Text.StringBuilder(seq.Length);
            // initialize_masking_buffer(word_length + m3p).
            _ei = MaskBufferSize - (w + m3p) + 1;
            ulong wordFwd = 0, wordRev = 0;
            int currentLength = 0;
            bool initRound = true;
            foreach (char c in seq)
            {
                if (!initRound && _wi == _ri)
                    Empty(fwdOut, revOut, flushAll: false);
                initRound = false;
                AddChar(c);
                ulong nv = MaskNucleotideValue(c, nameof(seq));
                wordFwd = (wordFwd << 2) | nv;
                wordRev = (wordRev >> 2) | ((~nv & 3) << ((w - 1) * 2));
                currentLength++;
                if (currentLength > w)
                {
                    wordFwd &= mask;
                    wordRev &= mask;
                    currentLength = w;
                }
                if (currentLength == w)
                    MaskOligoRegion(wordFwd, wordRev);
            }
            Empty(fwdOut, revOut, flushAll: true);
            return (fwdOut.ToString(), revOut.ToString());
        }

        // add_char_to_buffer for a nucleotide (non_nucleotide_positions stay 0).
        private void AddChar(char c)
        {
            _buffer[_wi] = c;
            _maskFwd[_wi] = false;
            _maskRev[_wi] = false;
            if (_mi > 0)
            {
                _maskFwd[_wi] = true;
                _mi--;
            }
            _ei = Forward(_ei);
            _wi = Forward(_wi);
        }

        private void Empty(System.Text.StringBuilder fwd, System.Text.StringBuilder rev, bool flushAll)
        {
            int end = flushAll ? _wi : _ei;
            while (_ri != end)
            {
                char b = _buffer[_ri];
                if (b >= 'a') // COND0: already lower case
                {
                    fwd.Append(b);
                    rev.Append(b);
                }
                else
                {
                    fwd.Append(_maskFwd[_ri] ? (char)(b + 32) : b);
                    rev.Append(_maskRev[_ri] ? (char)(b + 32) : b);
                }
                _ri = Forward(_ri);
            }
        }

        private void MaskOligoRegion(ulong wordFwd, ulong wordRev)
        {
            var (scoreFwd, scoreRev) = lists.Scores(wordFwd, wordRev, MaskWordLength2);
            if (failureRate != 0 && scoreFwd > failureRate)
            {
                int masked = 0, i = Back(_wi);
                while (masked < m5p)
                {
                    if (!_maskFwd[i]) _maskFwd[i] = true;
                    masked++;
                    i = Back(i);
                }
                _mi = (uint)m3p;
            }
            if (failureRate != 0 && scoreRev > failureRate)
            {
                int masked = 0, i = Back(_ei);
                while (masked < m5p + m3p)
                {
                    if (!_maskRev[i]) _maskRev[i] = true;
                    masked++;
                    i = Forward(i);
                }
            }
        }
    }
}

/// <summary>
/// The genome k-mer count lists of Primer3's template masker (PRIMER_MASK_KMERLIST_PATH / PRIMER_MASK_KMERLIST_PREFIX: the
/// GenomeTester4 <c>glistmaker</c> lists <c>&lt;prefix&gt;_11.list</c> and <c>&lt;prefix&gt;_16.list</c>), used by
/// <see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>, <see cref="PrimerDesigner.MaskTemplatePrimer3"/> and
/// <see cref="PrimerPairOptions.MaskTemplate"/>. A k-mer is looked up as given, then as its reverse complement; one found in
/// neither (or with count 0) counts as 1, as in masker.c <c>get_frequency_of_canonical_oligo</c>.
/// </summary>
public sealed class PrimerMaskingKmerLists
{
    private readonly Dictionary<ulong, uint> _list11, _list16;

    /// <summary>Builds the lists from caller-supplied k-mer → genome count maps.</summary>
    /// <param name="kmers11">11-mer counts (keys: 11 A/C/G/T, case-insensitive; counts ≥ 0).</param>
    /// <param name="kmers16">16-mer counts (keys: 16 A/C/G/T; counts ≥ 0).</param>
    /// <exception cref="ArgumentException">An empty list (Primer3 "List file contains no kmers"), a key of the wrong
    /// length or with another character, a negative count, or a key given twice in different case.</exception>
    public PrimerMaskingKmerLists(IReadOnlyDictionary<string, int> kmers11, IReadOnlyDictionary<string, int> kmers16)
    {
        ArgumentNullException.ThrowIfNull(kmers11);
        ArgumentNullException.ThrowIfNull(kmers16);
        _list11 = Build(kmers11, PrimerDesigner.MaskWordLength1, nameof(kmers11));
        _list16 = Build(kmers16, PrimerDesigner.MaskWordLength2, nameof(kmers16));
    }

    private PrimerMaskingKmerLists(Dictionary<ulong, uint> list11, Dictionary<ulong, uint> list16)
    {
        _list11 = list11;
        _list16 = list16;
    }

    /// <summary>Number of 11-mers in the list.</summary>
    public int Count11 => _list11.Count;

    /// <summary>Number of 16-mers in the list.</summary>
    public int Count16 => _list16.Count;

    /// <summary>
    /// Reads the two GenomeTester4 binary lists Primer3 loads for PRIMER_MASK_KMERLIST_PATH / _PREFIX: the file names are
    /// <c>{kmerListsPath}{prefix}_11.list</c> and <c>{kmerListsPath}{prefix}_16.list</c> (plain concatenation, as masker.c
    /// <c>create_formula_parameters_from_list_file_prefix</c> builds them — end the path with a separator). Format
    /// (masker.c <c>create_formula_parameters_from_list_file_name</c>): little-endian magic "GT4C" at byte 0, the k-mer
    /// length (uint32) at byte 12, the number of k-mers (uint32) at byte 16, the header size (uint64) at byte 32, then
    /// (uint64 2-bit-packed k-mer, A=0 C=1 G=2 T=3, uint32 count) records of 12 bytes.
    /// </summary>
    /// <param name="kmerListsPath">PRIMER_MASK_KMERLIST_PATH (directory with a trailing separator).</param>
    /// <param name="prefix">PRIMER_MASK_KMERLIST_PREFIX (Primer3 default "homo_sapiens").</param>
    /// <exception cref="ArgumentException">A missing file, a file that is not a list ("Given file is not a list file"), a list
    /// with no k-mers ("List file contains no kmers") or with a k-mer length other than 11 / 16.</exception>
    public static PrimerMaskingKmerLists FromGenomeTester4Files(string kmerListsPath, string prefix = "homo_sapiens")
    {
        ArgumentNullException.ThrowIfNull(kmerListsPath);
        ArgumentNullException.ThrowIfNull(prefix);
        return new PrimerMaskingKmerLists(
            ReadList($"{kmerListsPath}{prefix}_{PrimerDesigner.MaskWordLength1}.list", PrimerDesigner.MaskWordLength1),
            ReadList($"{kmerListsPath}{prefix}_{PrimerDesigner.MaskWordLength2}.list", PrimerDesigner.MaskWordLength2));
    }

    private static Dictionary<ulong, uint> ReadList(string path, int k)
    {
        if (!File.Exists(path))
            throw new ArgumentException($"Cannot find list file {path} (Primer3 masker).", nameof(path));
        byte[] data = File.ReadAllBytes(path);
        const uint magic = ('G' << 24) | ('T' << 16) | ('4' << 8) | 'C';
        if (data.Length < 40 || BinaryPrimitives.ReadUInt32LittleEndian(data) != magic)
            throw new ArgumentException($"Given file is not a list file: {path}", nameof(path));
        uint length = BinaryPrimitives.ReadUInt32LittleEndian(data.AsSpan(12));
        uint words = BinaryPrimitives.ReadUInt32LittleEndian(data.AsSpan(16));
        ulong header = BinaryPrimitives.ReadUInt64LittleEndian(data.AsSpan(32));
        if (words == 0)
            throw new ArgumentException($"List file contains no kmers: {path}", nameof(path));
        if (length != k)
            throw new ArgumentException($"{path}: k-mer length {length}, expected {k}.", nameof(path));
        if (header > (ulong)data.Length || (ulong)data.Length - header < 12UL * words)
            throw new ArgumentException($"{path}: truncated list file.", nameof(path));
        var list = new Dictionary<ulong, uint>((int)words);
        for (long r = 0; r < words; r++)
        {
            int at = checked((int)header + (int)(12 * r));
            ulong word = BinaryPrimitives.ReadUInt64LittleEndian(data.AsSpan(at));
            uint count = BinaryPrimitives.ReadUInt32LittleEndian(data.AsSpan(at + 8));
            // binary_search finds one record per word; a sorted list holds each word once (first kept otherwise).
            list.TryAdd(word, count);
        }
        return list;
    }

    private static Dictionary<ulong, uint> Build(IReadOnlyDictionary<string, int> kmers, int k, string paramName)
    {
        if (kmers.Count == 0)
            throw new ArgumentException("List contains no kmers (Primer3 masker).", paramName);
        var list = new Dictionary<ulong, uint>(kmers.Count);
        foreach (var (kmer, count) in kmers)
        {
            if (kmer is null || kmer.Length != k)
                throw new ArgumentException($"Every k-mer must have length {k}.", paramName);
            if (count < 0)
                throw new ArgumentException("K-mer counts must be ≥ 0.", paramName);
            ulong word = 0;
            foreach (char c in kmer)
                word = (word << 2) | c switch
                {
                    'A' or 'a' => 0UL,
                    'C' or 'c' => 1UL,
                    'G' or 'g' => 2UL,
                    'T' or 't' => 3UL,
                    _ => throw new ArgumentException($"K-mer '{kmer}' contains '{c}' (A/C/G/T only).", paramName),
                };
            if (!list.TryAdd(word, (uint)count))
                throw new ArgumentException($"K-mer '{kmer}' is given twice.", paramName);
        }
        return list;
    }

    // masker.c get_frequency_of_canonical_oligo: the word's count, else its reverse complement's, else 1.
    private static uint Frequency(Dictionary<ulong, uint> list, ulong word, int k)
    {
        if (list.TryGetValue(word, out uint f) && f != 0)
            return f;
        return list.TryGetValue(ReverseComplement(word, k), out uint r) && r != 0 ? r : 1;
    }

    private static ulong ReverseComplement(ulong word, int k)
    {
        word = ~word;
        ulong rc = 0;
        for (int i = 0; i < k; i++)
        {
            rc = (rc << 2) | (word & 3);
            word >>= 2;
        }
        return rc;
    }

    // masker.c calculate_scores (both_separately, abs_cutoff 0, the default formula: mm0 coefficients only) for a window of
    // windowLength nucleotides: lists of the window length count the forward word for both strands; shorter lists count
    // the trailing k-mer of the forward word and of the reverse-complement word separately; each strand score becomes the
    // logistic failure rate unless it is exactly 0.
    internal (double Forward, double Reverse) Scores(ulong wordFwd, ulong wordRev, int windowLength)
    {
        double fwd = 0.0, rev = 0.0;
        Add(_list11, PrimerDesigner.MaskWordLength1, PrimerDesigner.MaskCoefficient1);
        Add(_list16, PrimerDesigner.MaskWordLength2, PrimerDesigner.MaskCoefficient2);
        const double b = PrimerDesigner.MaskIntercept;
        if (fwd != 0) fwd = Math.Exp(fwd + b) / (1 + Math.Exp(fwd + b));
        if (rev != 0) rev = Math.Exp(rev + b) / (1 + Math.Exp(rev + b));
        return (fwd, rev);

        void Add(Dictionary<ulong, uint> list, int k, double coefficient)
        {
            ulong mask = (1UL << (2 * k)) - 1;
            if (windowLength == k)
            {
                // score += mm0·log(count) + mm0_2·log(count)² (mm0_2 = 0); added to both strands.
                double lc = Math.Log(Frequency(list, wordFwd & mask, k));
                double score = 0.0 + (coefficient * lc + 0.0 * lc * lc);
                fwd += score;
                rev += score;
            }
            else
            {
                double lf = Math.Log(Frequency(list, wordFwd & mask, k));
                fwd += coefficient * lf + 0.0 * lf * lf;
                double lr = Math.Log(Frequency(list, wordRev & mask, k));
                rev += coefficient * lr + 0.0 * lr * lr;
            }
        }
    }
}
