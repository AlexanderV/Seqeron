using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Core;

/// <summary>
/// Base interface for all biological sequences.
/// </summary>
public interface ISequence
{
    /// <summary>
    /// Gets the sequence data as a string.
    /// </summary>
    string Sequence { get; }

    /// <summary>
    /// Gets the length of the sequence.
    /// </summary>
    int Length { get; }

    /// <summary>
    /// Gets the sequence type.
    /// </summary>
    SequenceType Type { get; }

    /// <summary>
    /// Gets the valid alphabet for this sequence type.
    /// </summary>
    IReadOnlySet<char> Alphabet { get; }

    /// <summary>
    /// Gets the character at the specified position.
    /// </summary>
    char this[int index] { get; }

    /// <summary>
    /// Gets a subsequence.
    /// </summary>
    ISequence Subsequence(int start, int length);

    /// <summary>
    /// Validates the sequence against the alphabet.
    /// </summary>
    bool IsValid();

    /// <summary>
    /// Gets the complement of the sequence (for nucleotides).
    /// </summary>
    ISequence? GetComplement();

    /// <summary>
    /// Gets the reverse of the sequence.
    /// </summary>
    ISequence GetReverse();

    /// <summary>
    /// Gets the reverse complement (for nucleotides).
    /// </summary>
    ISequence? GetReverseComplement();
}

/// <summary>
/// Sequence type enumeration.
/// </summary>
public enum SequenceType
{
    Dna,
    Rna,
    Protein,
    IupacDna,
    IupacRna,
    Quality
}

/// <summary>
/// Base class for biological sequences with common functionality.
/// </summary>
public abstract class SequenceBase : ISequence
{
    protected readonly string _sequence;

    protected SequenceBase(string sequence)
    {
        _sequence = sequence?.ToUpperInvariant() ?? throw new ArgumentNullException(nameof(sequence));
    }

    public string Sequence => _sequence;
    public int Length => _sequence.Length;
    public abstract SequenceType Type { get; }
    public abstract IReadOnlySet<char> Alphabet { get; }

    public char this[int index] => _sequence[index];

    public abstract ISequence Subsequence(int start, int length);

    public virtual bool IsValid()
    {
        foreach (char c in _sequence)
        {
            if (!Alphabet.Contains(c))
                return false;
        }
        return true;
    }

    public abstract ISequence? GetComplement();
    public abstract ISequence GetReverse();
    public abstract ISequence? GetReverseComplement();

    public override string ToString() => _sequence;
    public override int GetHashCode() => _sequence.GetHashCode();
    public override bool Equals(object? obj) => obj is SequenceBase other && _sequence == other._sequence;
}

/// <summary>
/// IUPAC DNA sequence with ambiguity codes support.
/// Supports: A, C, G, T, N (any), R (A/G), Y (C/T), W (A/T), S (G/C), K (G/T), M (A/C),
/// B (C/G/T), D (A/G/T), H (A/C/T), V (A/C/G)
/// </summary>
public class IupacDnaSequence : SequenceBase
{
    private static readonly HashSet<char> _alphabet = new()
    {
        'A', 'C', 'G', 'T', 'U',  // Standard bases
        'N',                      // Any base
        'R', 'Y',                 // Purine (A/G), Pyrimidine (C/T)
        'W', 'S',                 // Weak (A/T), Strong (G/C)
        'K', 'M',                 // Keto (G/T), Amino (A/C)
        'B', 'D', 'H', 'V',       // 3-base ambiguity
        '-', '.'                  // Gaps
    };

    private static readonly Dictionary<char, char> _complements = new()
    {
        ['A'] = 'T',
        ['T'] = 'A',
        ['G'] = 'C',
        ['C'] = 'G',
        ['U'] = 'A',
        ['N'] = 'N',
        ['R'] = 'Y',
        ['Y'] = 'R',  // Purine <-> Pyrimidine
        ['W'] = 'W',
        ['S'] = 'S',  // Weak and Strong are self-complementary
        ['K'] = 'M',
        ['M'] = 'K',  // Keto <-> Amino
        ['B'] = 'V',
        ['V'] = 'B',  // B(CGT) <-> V(ACG)
        ['D'] = 'H',
        ['H'] = 'D',  // D(AGT) <-> H(ACT)
        ['-'] = '-',
        ['.'] = '.'
    };

    private static readonly Dictionary<char, char[]> _expansions = new()
    {
        ['A'] = new[] { 'A' },
        ['C'] = new[] { 'C' },
        ['G'] = new[] { 'G' },
        ['T'] = new[] { 'T' },
        ['U'] = new[] { 'T' },
        ['N'] = new[] { 'A', 'C', 'G', 'T' },
        ['R'] = new[] { 'A', 'G' },
        ['Y'] = new[] { 'C', 'T' },
        ['W'] = new[] { 'A', 'T' },
        ['S'] = new[] { 'G', 'C' },
        ['K'] = new[] { 'G', 'T' },
        ['M'] = new[] { 'A', 'C' },
        ['B'] = new[] { 'C', 'G', 'T' },
        ['D'] = new[] { 'A', 'G', 'T' },
        ['H'] = new[] { 'A', 'C', 'T' },
        ['V'] = new[] { 'A', 'C', 'G' }
    };

    public IupacDnaSequence(string sequence) : base(sequence) { }

    public override SequenceType Type => SequenceType.IupacDna;
    public override IReadOnlySet<char> Alphabet => _alphabet;

    public override ISequence Subsequence(int start, int length)
        => new IupacDnaSequence(_sequence.Substring(start, length));

    public override ISequence? GetComplement()
    {
        var complement = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            if (_complements.TryGetValue(_sequence[i], out char comp))
                complement[i] = comp;
            else
                complement[i] = 'N';
        }
        return new IupacDnaSequence(new string(complement));
    }

    public override ISequence GetReverse()
    {
        var reversed = new char[Length];
        for (int i = 0; i < Length; i++)
            reversed[i] = _sequence[Length - 1 - i];
        return new IupacDnaSequence(new string(reversed));
    }

    public override ISequence? GetReverseComplement()
    {
        var result = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            char c = _sequence[Length - 1 - i];
            result[i] = _complements.TryGetValue(c, out char comp) ? comp : 'N';
        }
        return new IupacDnaSequence(new string(result));
    }

    /// <summary>
    /// Expands IUPAC code to possible bases.
    /// </summary>
    public static char[] ExpandCode(char iupacCode)
    {
        return _expansions.TryGetValue(char.ToUpperInvariant(iupacCode), out var bases)
            ? bases
            : new[] { 'N' };
    }

    /// <summary>
    /// Gets IUPAC code from a set of bases.
    /// </summary>
    public static char GetIupacCode(IEnumerable<char> bases)
    {
        var baseSet = new HashSet<char>(bases.Select(char.ToUpperInvariant));

        if (baseSet.SetEquals(new[] { 'A' })) return 'A';
        if (baseSet.SetEquals(new[] { 'C' })) return 'C';
        if (baseSet.SetEquals(new[] { 'G' })) return 'G';
        if (baseSet.SetEquals(new[] { 'T' })) return 'T';
        if (baseSet.SetEquals(new[] { 'A', 'G' })) return 'R';
        if (baseSet.SetEquals(new[] { 'C', 'T' })) return 'Y';
        if (baseSet.SetEquals(new[] { 'A', 'T' })) return 'W';
        if (baseSet.SetEquals(new[] { 'G', 'C' })) return 'S';
        if (baseSet.SetEquals(new[] { 'G', 'T' })) return 'K';
        if (baseSet.SetEquals(new[] { 'A', 'C' })) return 'M';
        if (baseSet.SetEquals(new[] { 'C', 'G', 'T' })) return 'B';
        if (baseSet.SetEquals(new[] { 'A', 'G', 'T' })) return 'D';
        if (baseSet.SetEquals(new[] { 'A', 'C', 'T' })) return 'H';
        if (baseSet.SetEquals(new[] { 'A', 'C', 'G' })) return 'V';

        return 'N';
    }

    /// <summary>
    /// Checks if pattern matches sequence at position (with IUPAC wildcards).
    /// </summary>
    public bool MatchesAt(string pattern, int position)
    {
        if (position < 0 || position + pattern.Length > Length)
            return false;

        for (int i = 0; i < pattern.Length; i++)
        {
            char p = char.ToUpperInvariant(pattern[i]);
            char s = _sequence[position + i];

            if (!CodesMatch(p, s))
                return false;
        }

        return true;
    }

    /// <summary>
    /// Checks if two IUPAC codes can represent the same base.
    /// </summary>
    public static bool CodesMatch(char code1, char code2)
    {
        var bases1 = _expansions.TryGetValue(code1, out var b1) ? b1 : new[] { code1 };
        var bases2 = _expansions.TryGetValue(code2, out var b2) ? b2 : new[] { code2 };

        return bases1.Intersect(bases2).Any();
    }

    /// <summary>
    /// Finds all positions where pattern matches (with IUPAC support).
    /// </summary>
    public IEnumerable<int> FindPattern(string pattern)
    {
        for (int i = 0; i <= Length - pattern.Length; i++)
        {
            if (MatchesAt(pattern, i))
                yield return i;
        }
    }

    /// <summary>
    /// Calculates ambiguity level (1.0 = no ambiguity, 0.0 = all N).
    /// </summary>
    public double GetAmbiguityLevel()
    {
        if (Length == 0) return 1.0;

        int unambiguous = _sequence.Count(c => c == 'A' || c == 'C' || c == 'G' || c == 'T');
        return unambiguous / (double)Length;
    }

    /// <summary>
    /// Generates all possible concrete sequences from IUPAC sequence.
    /// Warning: exponential complexity for many ambiguous positions.
    /// </summary>
    public IEnumerable<string> ExpandAll(int maxResults = 1000)
    {
        var results = new List<string> { "" };

        foreach (char c in _sequence)
        {
            var expansions = ExpandCode(c);
            var newResults = new List<string>();

            foreach (var prefix in results)
            {
                foreach (var b in expansions)
                {
                    newResults.Add(prefix + b);
                    if (newResults.Count >= maxResults)
                    {
                        foreach (var r in newResults)
                            yield return r;
                        yield break;
                    }
                }
            }

            results = newResults;
        }

        foreach (var r in results)
            yield return r;
    }
}

/// <summary>
/// Sequence with per-base quality scores (FASTQ).
/// </summary>
public class QualitySequence : SequenceBase
{
    private static readonly HashSet<char> _alphabet = new()
    {
        'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'
    };

    private readonly byte[] _qualities;

    public QualitySequence(string sequence, byte[] qualities) : base(sequence)
    {
        if (qualities.Length != sequence.Length)
            throw new ArgumentException("Quality array must match sequence length");

        _qualities = qualities;
    }

    public QualitySequence(string sequence, string qualityString, int phredOffset = 33)
        : base(sequence)
    {
        _qualities = new byte[sequence.Length];
        for (int i = 0; i < qualityString.Length && i < sequence.Length; i++)
        {
            _qualities[i] = (byte)Math.Max(0, qualityString[i] - phredOffset);
        }
    }

    public override SequenceType Type => SequenceType.Quality;
    public override IReadOnlySet<char> Alphabet => _alphabet;

    /// <summary>
    /// Gets the quality scores.
    /// </summary>
    public IReadOnlyList<byte> Qualities => _qualities;

    /// <summary>
    /// Gets quality at position.
    /// </summary>
    public byte GetQuality(int index) => _qualities[index];

    /// <summary>
    /// Gets mean quality score.
    /// </summary>
    public double MeanQuality => _qualities.Average(q => (double)q);

    /// <summary>
    /// Gets the quality string (Phred+33 encoding).
    /// </summary>
    public string GetQualityString(int phredOffset = 33)
    {
        var chars = new char[_qualities.Length];
        for (int i = 0; i < _qualities.Length; i++)
        {
            chars[i] = (char)(_qualities[i] + phredOffset);
        }
        return new string(chars);
    }

    public override ISequence Subsequence(int start, int length)
    {
        var subQual = new byte[length];
        Array.Copy(_qualities, start, subQual, 0, length);
        return new QualitySequence(_sequence.Substring(start, length), subQual);
    }

    public override ISequence? GetComplement()
    {
        var complement = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            char comp = SequenceExtensions.GetComplementBase(_sequence[i]);
            complement[i] = comp == _sequence[i] ? 'N' : comp; // Unknown bases become N
        }
        return new QualitySequence(new string(complement), _qualities);
    }

    public override ISequence GetReverse()
    {
        var reversed = new char[Length];
        var revQual = new byte[Length];
        for (int i = 0; i < Length; i++)
        {
            reversed[i] = _sequence[Length - 1 - i];
            revQual[i] = _qualities[Length - 1 - i];
        }
        return new QualitySequence(new string(reversed), revQual);
    }

    public override ISequence? GetReverseComplement()
    {
        var result = new char[Length];
        var revQual = new byte[Length];
        for (int i = 0; i < Length; i++)
        {
            char c = _sequence[Length - 1 - i];
            char comp = SequenceExtensions.GetComplementBase(c);
            result[i] = comp == c ? 'N' : comp; // Unknown bases become N
            revQual[i] = _qualities[Length - 1 - i];
        }
        return new QualitySequence(new string(result), revQual);
    }

    /// <summary>
    /// Trims low quality bases from ends.
    /// </summary>
    public QualitySequence TrimByQuality(byte minQuality = 20)
    {
        int start = 0;
        int end = Length - 1;

        while (start < Length && _qualities[start] < minQuality)
            start++;

        while (end > start && _qualities[end] < minQuality)
            end--;

        if (start > end)
            return new QualitySequence("", Array.Empty<byte>());

        int len = end - start + 1;
        var trimmedQual = new byte[len];
        Array.Copy(_qualities, start, trimmedQual, 0, len);

        return new QualitySequence(_sequence.Substring(start, len), trimmedQual);
    }

    /// <summary>
    /// Masks low quality bases with N.
    /// </summary>
    public QualitySequence MaskLowQuality(byte minQuality = 20)
    {
        var masked = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            masked[i] = _qualities[i] >= minQuality ? _sequence[i] : 'N';
        }
        return new QualitySequence(new string(masked), _qualities);
    }

    /// <summary>
    /// Gets error probability from Phred score.
    /// </summary>
    public static double PhredToErrorProbability(byte phred)
        => Math.Pow(10, -phred / 10.0);

    /// <summary>
    /// Gets Phred score from error probability.
    /// </summary>
    public static byte ErrorProbabilityToPhred(double errorProb)
        => (byte)Math.Max(0, Math.Min(93, -10 * Math.Log10(errorProb)));

    /// <summary>
    /// Calculates expected number of errors.
    /// </summary>
    public double ExpectedErrors()
        => _qualities.Sum(q => PhredToErrorProbability(q));
}
