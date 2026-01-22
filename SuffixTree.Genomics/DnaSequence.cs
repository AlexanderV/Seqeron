using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics
{
    /// <summary>
    /// Represents a DNA sequence with validation and common operations.
    /// Valid nucleotides: A (Adenine), C (Cytosine), G (Guanine), T (Thymine).
    /// </summary>
    public sealed class DnaSequence
    {
        private readonly string _sequence;
        private SuffixTree? _suffixTree;

        /// <summary>
        /// Creates a new DNA sequence from a string.
        /// </summary>
        /// <param name="sequence">DNA sequence string (case-insensitive).</param>
        /// <exception cref="ArgumentException">Thrown if sequence contains invalid characters.</exception>
        public DnaSequence(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
            {
                _sequence = string.Empty;
                return;
            }

            var normalized = sequence.ToUpperInvariant();
            ValidateSequence(normalized);
            _sequence = normalized;
        }

        /// <summary>
        /// Gets the DNA sequence string.
        /// </summary>
        public string Sequence => _sequence;

        /// <summary>
        /// Gets the length of the sequence.
        /// </summary>
        public int Length => _sequence.Length;

        /// <summary>
        /// Gets or builds the suffix tree for this sequence.
        /// </summary>
        public SuffixTree SuffixTree => _suffixTree ??= SuffixTree.Build(_sequence);

        /// <summary>
        /// Gets the complement of this DNA sequence.
        /// A ↔ T, C ↔ G
        /// </summary>
        public DnaSequence Complement()
        {
            var sb = new StringBuilder(_sequence.Length);
            foreach (char c in _sequence)
            {
                sb.Append(c switch
                {
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => c
                });
            }
            return new DnaSequence(sb.ToString());
        }

        /// <summary>
        /// Gets the reverse complement of this DNA sequence.
        /// This is important for double-stranded DNA analysis.
        /// </summary>
        public DnaSequence ReverseComplement()
        {
            var sb = new StringBuilder(_sequence.Length);
            for (int i = _sequence.Length - 1; i >= 0; i--)
            {
                sb.Append(_sequence[i] switch
                {
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => _sequence[i]
                });
            }
            return new DnaSequence(sb.ToString());
        }

        /// <summary>
        /// Calculates GC content (percentage of G and C nucleotides).
        /// Higher GC content = higher melting temperature.
        /// </summary>
        public double GcContent()
        {
            if (_sequence.Length == 0) return 0;

            int gcCount = _sequence.Count(c => c == 'G' || c == 'C');
            return (double)gcCount / _sequence.Length * 100;
        }

        /// <summary>
        /// Transcribes DNA to RNA (T → U).
        /// </summary>
        public string Transcribe()
        {
            return _sequence.Replace('T', 'U');
        }

        /// <summary>
        /// Gets the nucleotide at the specified position.
        /// </summary>
        public char this[int index] => _sequence[index];

        /// <summary>
        /// Gets a subsequence (substring) of the DNA.
        /// </summary>
        public DnaSequence Subsequence(int start, int length)
        {
            return new DnaSequence(_sequence.Substring(start, length));
        }

        public override string ToString() => _sequence;

        public override bool Equals(object? obj) =>
            obj is DnaSequence other && _sequence == other._sequence;

        public override int GetHashCode() => _sequence.GetHashCode();

        private static void ValidateSequence(string sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                char c = sequence[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                {
                    throw new ArgumentException(
                        $"Invalid nucleotide '{c}' at position {i}. Valid nucleotides: A, C, G, T.",
                        nameof(sequence));
                }
            }
        }

        /// <summary>
        /// Tries to create a DNA sequence, returning false if invalid.
        /// </summary>
        public static bool TryCreate(string sequence, out DnaSequence? result)
        {
            try
            {
                result = new DnaSequence(sequence);
                return true;
            }
            catch (ArgumentException)
            {
                result = null;
                return false;
            }
        }

        /// <summary>
        /// Gets the reverse complement of a DNA sequence string.
        /// Static helper method for use when a full DnaSequence object is not needed.
        /// </summary>
        /// <param name="sequence">DNA sequence string.</param>
        /// <returns>Reverse complement string.</returns>
        public static string GetReverseComplementString(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return sequence;

            var result = new char[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
            {
                result[sequence.Length - 1 - i] = char.ToUpperInvariant(sequence[i]) switch
                {
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    'U' => 'A', // Support RNA too
                    _ => sequence[i]
                };
            }
            return new string(result);
        }
    }
}
