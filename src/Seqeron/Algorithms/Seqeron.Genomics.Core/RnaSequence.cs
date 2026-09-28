using System.Diagnostics.CodeAnalysis;

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
        /// A ↔ U, C ↔ G. Delegates to the canonical <see cref="SequenceExtensions.GetRnaComplementBase(char)"/>.
        /// </summary>
        public RnaSequence Complement()
        {
            var result = new char[_sequence.Length];
            for (int i = 0; i < _sequence.Length; i++)
            {
                result[i] = SequenceExtensions.GetRnaComplementBase(_sequence[i]);
            }
            return new RnaSequence(new string(result));
        }

        /// <summary>
        /// Gets the reverse complement of this RNA sequence.
        /// Delegates per base to the canonical <see cref="SequenceExtensions.GetRnaComplementBase(char)"/>.
        /// </summary>
        public RnaSequence ReverseComplement()
        {
            var result = new char[_sequence.Length];
            for (int i = 0; i < _sequence.Length; i++)
            {
                result[i] = SequenceExtensions.GetRnaComplementBase(_sequence[_sequence.Length - 1 - i]);
            }
            return new RnaSequence(new string(result));
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
            return GetCodonsCore(frame);
        }

        private IEnumerable<string> GetCodonsCore(int frame)
        {
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
            // Canonical predicate: SequenceExtensions.IndexOfInvalidRna (SEQ-VALID-001).
            int i = sequence.AsSpan().IndexOfInvalidRna();
            if (i >= 0)
            {
                throw new ArgumentException(
                    $"Invalid nucleotide '{sequence[i]}' at position {i}. Valid nucleotides: A, C, G, U.",
                    nameof(sequence));
            }
        }

        /// <summary>
        /// Tries to create an RNA sequence, returning false if invalid.
        /// </summary>
        public static bool TryCreate(string sequence, [NotNullWhen(true)] out RnaSequence? result)
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
