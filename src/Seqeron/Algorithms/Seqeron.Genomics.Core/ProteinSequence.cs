using System.Diagnostics.CodeAnalysis;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Represents a protein (amino acid) sequence with validation and common operations.
    /// Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (20 standard amino acids).
    /// Also supports: * (stop codon), X (unknown amino acid) and the IUPAC-IUBMB ambiguity codes
    /// B (Asx = D/N), Z (Glx = E/Q) and J (Xle = I/L), which <see cref="GeneticCode.Translate(string)"/>
    /// emits for ambiguous codons (Biopython <c>Bio.Seq._translate_str</c> / <c>IUPACData.extended_protein_letters</c>).
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
        /// IUPAC-IUBMB ambiguous amino-acid codes with their three-letter symbols
        /// (Biopython <c>IUPACData.protein_letters_1to3_extended</c>): B = Asx (D or N),
        /// Z = Glx (E or Q), J = Xle (I or L). Produced by translating ambiguous codons
        /// such as RAY, SAR and MTH (see <see cref="GeneticCode.Translate(string)"/>).
        /// </summary>
        public static readonly IReadOnlyDictionary<char, string> AmbiguousAminoAcids = new Dictionary<char, string>
        {
            ['B'] = "Asx",
            ['Z'] = "Glx",
            ['J'] = "Xle",
        };

        /// <summary>
        /// All valid characters (standard amino acids + stop + unknown + ambiguous B/Z/J).
        /// </summary>
        public static readonly IReadOnlySet<char> ValidCharacters = new HashSet<char>(
            StandardAminoAcids.Concat(new[] { '*', 'X' }).Concat(AmbiguousAminoAcids.Keys)
        );

        /// <summary>
        /// Amino acid properties. <see cref="AminoAcidProperties.MolecularWeight"/> is the average mass of
        /// the free amino acid from the canonical table <see cref="ProteinPhysicochemistry.AverageAminoAcidMasses"/>
        /// (Biopython <c>IUPACData.protein_weights</c>).
        /// </summary>
        public static readonly IReadOnlyDictionary<char, AminoAcidProperties> Properties = new Dictionary<char, AminoAcidProperties>
        {
            ['A'] = new("Alanine", "Ala", Mass('A'), AminoAcidType.Nonpolar),
            ['C'] = new("Cysteine", "Cys", Mass('C'), AminoAcidType.Polar),
            ['D'] = new("Aspartic acid", "Asp", Mass('D'), AminoAcidType.Acidic),
            ['E'] = new("Glutamic acid", "Glu", Mass('E'), AminoAcidType.Acidic),
            ['F'] = new("Phenylalanine", "Phe", Mass('F'), AminoAcidType.Nonpolar),
            ['G'] = new("Glycine", "Gly", Mass('G'), AminoAcidType.Nonpolar),
            ['H'] = new("Histidine", "His", Mass('H'), AminoAcidType.Basic),
            ['I'] = new("Isoleucine", "Ile", Mass('I'), AminoAcidType.Nonpolar),
            ['K'] = new("Lysine", "Lys", Mass('K'), AminoAcidType.Basic),
            ['L'] = new("Leucine", "Leu", Mass('L'), AminoAcidType.Nonpolar),
            ['M'] = new("Methionine", "Met", Mass('M'), AminoAcidType.Nonpolar),
            ['N'] = new("Asparagine", "Asn", Mass('N'), AminoAcidType.Polar),
            ['P'] = new("Proline", "Pro", Mass('P'), AminoAcidType.Nonpolar),
            ['Q'] = new("Glutamine", "Gln", Mass('Q'), AminoAcidType.Polar),
            ['R'] = new("Arginine", "Arg", Mass('R'), AminoAcidType.Basic),
            ['S'] = new("Serine", "Ser", Mass('S'), AminoAcidType.Polar),
            ['T'] = new("Threonine", "Thr", Mass('T'), AminoAcidType.Polar),
            ['V'] = new("Valine", "Val", Mass('V'), AminoAcidType.Nonpolar),
            ['W'] = new("Tryptophan", "Trp", Mass('W'), AminoAcidType.Nonpolar),
            ['Y'] = new("Tyrosine", "Tyr", Mass('Y'), AminoAcidType.Polar)
        };

        private static double Mass(char aminoAcid) => ProteinPhysicochemistry.AverageAminoAcidMasses[aminoAcid];

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
        /// Calculates the average-isotopic molecular weight of the protein in daltons: sum of the
        /// free amino-acid masses minus one water per peptide bond (Biopython
        /// <c>molecular_weight(seq, "protein")</c>; Expasy Compute pI/Mw), unrounded. Delegates to the
        /// canonical <see cref="ProteinPhysicochemistry.MolecularWeight"/> (the implementation shared with
        /// <c>SequenceStatistics.CalculateMolecularWeight</c>). Ambiguous residues (B, Z, J, X) and '*'
        /// carry no mass and are skipped. Returns 0 for an empty sequence.
        /// </summary>
        public double MolecularWeight() => ProteinPhysicochemistry.MolecularWeight(_sequence);

        /// <summary>
        /// Calculates the theoretical isoelectric point (pI) on the EMBOSS <c>iep</c> pK scale
        /// (Epk.dat: N-terminus 7.5, C-terminus 3.6, C 8.5, D 3.9, E 4.1, H 6.5, K 10.8, R 12.5, Y 10.1),
        /// located by bisection to 1e-9 pH and rounded to two decimals. Delegates to the canonical
        /// <see cref="ProteinPhysicochemistry.IsoelectricPoint"/> (shared with
        /// <c>SequenceStatistics.CalculateIsoelectricPoint</c>). Returns 0 for an empty sequence.
        /// </summary>
        public double IsoelectricPoint() =>
            _sequence.Length == 0 ? 0 : ProteinPhysicochemistry.IsoelectricPoint(_sequence);

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
        /// Calculates the grand average of hydropathy (GRAVY): the mean Kyte–Doolittle (1982) value
        /// over the standard residues, unrounded (Biopython <c>ProteinAnalysis.gravy()</c>).
        /// Positive = hydrophobic, negative = hydrophilic. Delegates to the canonical
        /// <see cref="ProteinPhysicochemistry.Gravy"/> (shared with
        /// <c>SequenceStatistics.CalculateHydrophobicity</c>); B, Z, J, X and '*' are skipped.
        /// </summary>
        public double Gravy() => ProteinPhysicochemistry.Gravy(_sequence);

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
                else if (AmbiguousAminoAcids.TryGetValue(aa, out var ambiguous))
                {
                    if (sb.Length > 0) sb.Append('-');
                    sb.Append(ambiguous);
                }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Finds all (overlapping) occurrences of an exact motif in the protein sequence, in ascending
        /// position order (case-insensitive). Uses the sequence's <see cref="SuffixTree"/>, the same
        /// exact-matching engine as <c>MotifFinder.FindExactMotif</c> (PAT-EXACT-001).
        /// </summary>
        public IEnumerable<int> FindMotif(string pattern)
        {
            if (string.IsNullOrEmpty(pattern) || pattern.Length > _sequence.Length)
                return Array.Empty<int>();

            return SuffixTree.FindAllOccurrences(pattern.ToUpperInvariant()).OrderBy(p => p);
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
                        $"Invalid amino acid '{c}' at position {i}. Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (and *, X, B, Z, J).",
                        nameof(sequence));
                }
            }
        }

        /// <summary>
        /// Tries to create a protein sequence, returning false if invalid.
        /// </summary>
        public static bool TryCreate(string sequence, [NotNullWhen(true)] out ProteinSequence? result)
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
