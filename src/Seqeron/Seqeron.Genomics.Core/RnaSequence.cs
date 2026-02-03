using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Represents an RNA sequence with validation and common operations.
    /// Valid nucleotides: A (Adenine), C (Cytosine), G (Guanine), U (Uracil).
    /// </summary>
    public sealed class RnaSequence
    {
        private readonly string _sequence;
        private SuffixTree.SuffixTree? _suffixTree;

        /// <summary>
        /// Creates a new RNA sequence from a string.
        /// </summary>
        /// <param name="sequence">RNA sequence string (case-insensitive).</param>
        /// <exception cref="ArgumentException">Thrown if sequence contains invalid characters.</exception>
        public RnaSequence(string sequence)
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
        /// Gets the RNA sequence string.
        /// </summary>
        public string Sequence => _sequence;

        /// <summary>
        /// Gets the length of the sequence.
        /// </summary>
        public int Length => _sequence.Length;

        /// <summary>
        /// Gets or builds the suffix tree for this sequence.
        /// </summary>
        public SuffixTree.SuffixTree SuffixTree => _suffixTree ??= global::SuffixTree.SuffixTree.Build(_sequence);

        /// <summary>
        /// Gets the complement of this RNA sequence.
        /// A ↔ U, C ↔ G
        /// </summary>
        public RnaSequence Complement()
        {
            var sb = new StringBuilder(_sequence.Length);
            foreach (char c in _sequence)
            {
                sb.Append(c switch
                {
                    'A' => 'U',
                    'U' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => c
                });
            }
            return new RnaSequence(sb.ToString());
        }

        /// <summary>
        /// Gets the reverse complement of this RNA sequence.
        /// </summary>
        public RnaSequence ReverseComplement()
        {
            var sb = new StringBuilder(_sequence.Length);
            for (int i = _sequence.Length - 1; i >= 0; i--)
            {
                sb.Append(_sequence[i] switch
                {
                    'A' => 'U',
                    'U' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => _sequence[i]
                });
            }
            return new RnaSequence(sb.ToString());
        }

        /// <summary>
        /// Calculates GC content (percentage of G and C nucleotides).
        /// </summary>
        public double GcContent() => _sequence.CalculateGcContentFast();

        /// <summary>
        /// Reverse transcribes RNA to DNA (U → T).
        /// </summary>
        public DnaSequence ReverseTranscribe()
        {
            return new DnaSequence(_sequence.Replace('U', 'T'));
        }

        /// <summary>
        /// Finds codons (triplets) in the RNA sequence.
        /// </summary>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <returns>Enumerable of codon strings.</returns>
        public IEnumerable<string> GetCodons(int frame = 0)
        {
            if (frame < 0 || frame > 2)
                throw new ArgumentOutOfRangeException(nameof(frame), "Frame must be 0, 1, or 2.");

            for (int i = frame; i + 3 <= _sequence.Length; i += 3)
            {
                yield return _sequence.Substring(i, 3);
            }
        }

        /// <summary>
        /// Gets the nucleotide at the specified position.
        /// </summary>
        public char this[int index] => _sequence[index];

        /// <summary>
        /// Gets a subsequence (substring) of the RNA.
        /// </summary>
        public RnaSequence Subsequence(int start, int length)
        {
            return new RnaSequence(_sequence.Substring(start, length));
        }

        /// <summary>
        /// Calculates AU content (percentage of A and U nucleotides).
        /// </summary>
        public double AuContent()
        {
            if (_sequence.Length == 0) return 0;

            int auCount = _sequence.Count(c => c == 'A' || c == 'U');
            return (double)auCount / _sequence.Length * 100;
        }

        /// <summary>
        /// Creates an RNA sequence from a DNA sequence (transcription).
        /// </summary>
        public static RnaSequence FromDna(DnaSequence dna)
        {
            return new RnaSequence(dna.Sequence.Replace('T', 'U'));
        }

        public override string ToString() => _sequence;

        public override bool Equals(object? obj) =>
            obj is RnaSequence other && _sequence == other._sequence;

        public override int GetHashCode() => _sequence.GetHashCode();

        private static void ValidateSequence(string sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                char c = sequence[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'U')
                {
                    throw new ArgumentException(
                        $"Invalid nucleotide '{c}' at position {i}. Valid nucleotides: A, C, G, U.",
                        nameof(sequence));
                }
            }
        }

        /// <summary>
        /// Tries to create an RNA sequence, returning false if invalid.
        /// </summary>
        public static bool TryCreate(string sequence, out RnaSequence? result)
        {
            try
            {
                result = new RnaSequence(sequence);
                return true;
            }
            catch (ArgumentException)
            {
                result = null;
                return false;
            }
        }
    }
}
