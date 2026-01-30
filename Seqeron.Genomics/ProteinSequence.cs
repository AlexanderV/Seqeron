using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics
{
    /// <summary>
    /// Represents a protein (amino acid) sequence with validation and common operations.
    /// Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (20 standard amino acids).
    /// Also supports: * (stop codon), X (unknown amino acid).
    /// </summary>
    public sealed class ProteinSequence
    {
        /// <summary>
        /// Standard 20 amino acids (single-letter codes).
        /// </summary>
        public static readonly IReadOnlySet<char> StandardAminoAcids = new HashSet<char>
        {
            'A', // Alanine
            'C', // Cysteine
            'D', // Aspartic acid
            'E', // Glutamic acid
            'F', // Phenylalanine
            'G', // Glycine
            'H', // Histidine
            'I', // Isoleucine
            'K', // Lysine
            'L', // Leucine
            'M', // Methionine
            'N', // Asparagine
            'P', // Proline
            'Q', // Glutamine
            'R', // Arginine
            'S', // Serine
            'T', // Threonine
            'V', // Valine
            'W', // Tryptophan
            'Y'  // Tyrosine
        };

        /// <summary>
        /// All valid characters (standard amino acids + stop + unknown).
        /// </summary>
        public static readonly IReadOnlySet<char> ValidCharacters = new HashSet<char>(
            StandardAminoAcids.Concat(new[] { '*', 'X' })
        );

        /// <summary>
        /// Amino acid properties.
        /// </summary>
        public static readonly IReadOnlyDictionary<char, AminoAcidProperties> Properties = new Dictionary<char, AminoAcidProperties>
        {
            ['A'] = new("Alanine", "Ala", 89.09, AminoAcidType.Nonpolar),
            ['C'] = new("Cysteine", "Cys", 121.16, AminoAcidType.Polar),
            ['D'] = new("Aspartic acid", "Asp", 133.10, AminoAcidType.Acidic),
            ['E'] = new("Glutamic acid", "Glu", 147.13, AminoAcidType.Acidic),
            ['F'] = new("Phenylalanine", "Phe", 165.19, AminoAcidType.Nonpolar),
            ['G'] = new("Glycine", "Gly", 75.07, AminoAcidType.Nonpolar),
            ['H'] = new("Histidine", "His", 155.16, AminoAcidType.Basic),
            ['I'] = new("Isoleucine", "Ile", 131.17, AminoAcidType.Nonpolar),
            ['K'] = new("Lysine", "Lys", 146.19, AminoAcidType.Basic),
            ['L'] = new("Leucine", "Leu", 131.17, AminoAcidType.Nonpolar),
            ['M'] = new("Methionine", "Met", 149.21, AminoAcidType.Nonpolar),
            ['N'] = new("Asparagine", "Asn", 132.12, AminoAcidType.Polar),
            ['P'] = new("Proline", "Pro", 115.13, AminoAcidType.Nonpolar),
            ['Q'] = new("Glutamine", "Gln", 146.15, AminoAcidType.Polar),
            ['R'] = new("Arginine", "Arg", 174.20, AminoAcidType.Basic),
            ['S'] = new("Serine", "Ser", 105.09, AminoAcidType.Polar),
            ['T'] = new("Threonine", "Thr", 119.12, AminoAcidType.Polar),
            ['V'] = new("Valine", "Val", 117.15, AminoAcidType.Nonpolar),
            ['W'] = new("Tryptophan", "Trp", 204.23, AminoAcidType.Nonpolar),
            ['Y'] = new("Tyrosine", "Tyr", 181.19, AminoAcidType.Polar)
        };

        private readonly string _sequence;
        private SuffixTree.SuffixTree? _suffixTree;

        /// <summary>
        /// Creates a new protein sequence from a string.
        /// </summary>
        /// <param name="sequence">Protein sequence string (case-insensitive).</param>
        /// <exception cref="ArgumentException">Thrown if sequence contains invalid characters.</exception>
        public ProteinSequence(string sequence)
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
        /// Gets the protein sequence string.
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
        /// Gets the amino acid at the specified position.
        /// </summary>
        public char this[int index] => _sequence[index];

        /// <summary>
        /// Gets a subsequence of the protein.
        /// </summary>
        public ProteinSequence Subsequence(int start, int length)
        {
            return new ProteinSequence(_sequence.Substring(start, length));
        }

        /// <summary>
        /// Calculates the molecular weight of the protein in Daltons.
        /// Uses average isotopic masses.
        /// </summary>
        public double MolecularWeight()
        {
            if (_sequence.Length == 0) return 0;

            // Sum of amino acid weights minus water (18.015 Da) for each peptide bond
            double weight = 0;
            int peptideBonds = 0;

            foreach (char aa in _sequence)
            {
                if (Properties.TryGetValue(aa, out var props))
                {
                    weight += props.MolecularWeight;
                    peptideBonds++;
                }
            }

            // Subtract water for peptide bond formation (n-1 peptide bonds for n amino acids)
            if (peptideBonds > 1)
            {
                weight -= (peptideBonds - 1) * 18.015;
            }

            return Math.Round(weight, 2);
        }

        /// <summary>
        /// Calculates the theoretical isoelectric point (pI) using pKa values.
        /// Simplified calculation using Henderson-Hasselbalch approximation.
        /// </summary>
        public double IsoelectricPoint()
        {
            if (_sequence.Length == 0) return 0;

            // Count charged amino acids
            int nTerm = 1; // N-terminus (pKa ~9.69)
            int cTerm = 1; // C-terminus (pKa ~2.34)
            int asp = _sequence.Count(c => c == 'D');  // pKa 3.9
            int glu = _sequence.Count(c => c == 'E');  // pKa 4.1
            int cys = _sequence.Count(c => c == 'C');  // pKa 8.3
            int tyr = _sequence.Count(c => c == 'Y');  // pKa 10.1
            int his = _sequence.Count(c => c == 'H');  // pKa 6.0
            int lys = _sequence.Count(c => c == 'K');  // pKa 10.5
            int arg = _sequence.Count(c => c == 'R');  // pKa 12.5

            // Binary search for pI
            double pHLow = 0;
            double pHHigh = 14;
            double pI = 7;

            while (pHHigh - pHLow > 0.01)
            {
                pI = (pHLow + pHHigh) / 2;
                double charge = CalculateCharge(pI, nTerm, cTerm, asp, glu, cys, tyr, his, lys, arg);

                if (charge > 0)
                    pHLow = pI;
                else
                    pHHigh = pI;
            }

            return Math.Round(pI, 2);
        }

        private static double CalculateCharge(double pH, int nTerm, int cTerm,
            int asp, int glu, int cys, int tyr, int his, int lys, int arg)
        {
            // Positive charges
            double positive =
                nTerm * (1 / (1 + Math.Pow(10, pH - 9.69))) +
                his * (1 / (1 + Math.Pow(10, pH - 6.0))) +
                lys * (1 / (1 + Math.Pow(10, pH - 10.5))) +
                arg * (1 / (1 + Math.Pow(10, pH - 12.5)));

            // Negative charges
            double negative =
                cTerm * (1 / (1 + Math.Pow(10, 2.34 - pH))) +
                asp * (1 / (1 + Math.Pow(10, 3.9 - pH))) +
                glu * (1 / (1 + Math.Pow(10, 4.1 - pH))) +
                cys * (1 / (1 + Math.Pow(10, 8.3 - pH))) +
                tyr * (1 / (1 + Math.Pow(10, 10.1 - pH)));

            return positive - negative;
        }

        /// <summary>
        /// Counts the occurrences of each amino acid.
        /// </summary>
        public IReadOnlyDictionary<char, int> AminoAcidComposition()
        {
            var composition = new Dictionary<char, int>();
            foreach (char aa in _sequence)
            {
                if (!composition.TryAdd(aa, 1))
                    composition[aa]++;
            }
            return composition;
        }

        /// <summary>
        /// Calculates hydropathicity index (GRAVY - Grand Average of Hydropathy).
        /// Positive = hydrophobic, negative = hydrophilic.
        /// Uses Kyte-Doolittle scale.
        /// </summary>
        public double Gravy()
        {
            if (_sequence.Length == 0) return 0;

            var hydropathy = new Dictionary<char, double>
            {
                ['A'] = 1.8,
                ['C'] = 2.5,
                ['D'] = -3.5,
                ['E'] = -3.5,
                ['F'] = 2.8,
                ['G'] = -0.4,
                ['H'] = -3.2,
                ['I'] = 4.5,
                ['K'] = -3.9,
                ['L'] = 3.8,
                ['M'] = 1.9,
                ['N'] = -3.5,
                ['P'] = -1.6,
                ['Q'] = -3.5,
                ['R'] = -4.5,
                ['S'] = -0.8,
                ['T'] = -0.7,
                ['V'] = 4.2,
                ['W'] = -0.9,
                ['Y'] = -1.3
            };

            double sum = 0;
            int count = 0;

            foreach (char aa in _sequence)
            {
                if (hydropathy.TryGetValue(aa, out double value))
                {
                    sum += value;
                    count++;
                }
            }

            return count > 0 ? Math.Round(sum / count, 3) : 0;
        }

        /// <summary>
        /// Calculates the percentage of a specific amino acid type.
        /// </summary>
        public double TypePercentage(AminoAcidType type)
        {
            if (_sequence.Length == 0) return 0;

            int count = _sequence.Count(aa =>
                Properties.TryGetValue(aa, out var props) && props.Type == type);

            return Math.Round((double)count / _sequence.Length * 100, 2);
        }

        /// <summary>
        /// Converts to three-letter code representation.
        /// </summary>
        public string ToThreeLetterCode()
        {
            var sb = new StringBuilder();
            foreach (char aa in _sequence)
            {
                if (Properties.TryGetValue(aa, out var props))
                {
                    if (sb.Length > 0) sb.Append('-');
                    sb.Append(props.ThreeLetterCode);
                }
                else if (aa == '*')
                {
                    if (sb.Length > 0) sb.Append('-');
                    sb.Append("Ter");
                }
                else if (aa == 'X')
                {
                    if (sb.Length > 0) sb.Append('-');
                    sb.Append("Xaa");
                }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Finds all occurrences of a motif pattern in the protein sequence.
        /// </summary>
        public IEnumerable<int> FindMotif(string pattern)
        {
            if (string.IsNullOrEmpty(pattern))
                yield break;

            var normalizedPattern = pattern.ToUpperInvariant();

            for (int i = 0; i <= _sequence.Length - normalizedPattern.Length; i++)
            {
                if (_sequence.Substring(i, normalizedPattern.Length) == normalizedPattern)
                    yield return i;
            }
        }

        public override string ToString() => _sequence;

        public override bool Equals(object? obj) =>
            obj is ProteinSequence other && _sequence == other._sequence;

        public override int GetHashCode() => _sequence.GetHashCode();

        private static void ValidateSequence(string sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                char c = sequence[i];
                if (!ValidCharacters.Contains(c))
                {
                    throw new ArgumentException(
                        $"Invalid amino acid '{c}' at position {i}. Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (and *, X).",
                        nameof(sequence));
                }
            }
        }

        /// <summary>
        /// Tries to create a protein sequence, returning false if invalid.
        /// </summary>
        public static bool TryCreate(string sequence, out ProteinSequence? result)
        {
            try
            {
                result = new ProteinSequence(sequence);
                return true;
            }
            catch (ArgumentException)
            {
                result = null;
                return false;
            }
        }
    }

    /// <summary>
    /// Type of amino acid based on side chain properties.
    /// </summary>
    public enum AminoAcidType
    {
        Nonpolar,  // A, F, G, I, L, M, P, V, W
        Polar,     // C, N, Q, S, T, Y
        Acidic,    // D, E
        Basic      // H, K, R
    }

    /// <summary>
    /// Properties of an amino acid.
    /// </summary>
    public readonly record struct AminoAcidProperties(
        string Name,
        string ThreeLetterCode,
        double MolecularWeight,
        AminoAcidType Type);
}
