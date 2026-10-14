using System.Diagnostics.CodeAnalysis;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Represents a protein (amino acid) sequence with validation and common operations.
    /// Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (20 standard amino acids).
    /// Also supports: * (stop codon), X (unknown amino acid), the IUPAC-IUBMB ambiguity codes
    /// B (Asx = D/N), Z (Glx = E/Q) and J (Xle = I/L), which <see cref="GeneticCode.Translate(string)"/>
    /// emits for ambiguous codons (Biopython <c>Bio.Seq._translate_str</c> / <c>IUPACData.extended_protein_letters</c>),
    /// and the 21st/22nd proteinogenic amino acids U (Sec, selenocysteine) and O (Pyl, pyrrolysine)
    /// (Biopython <c>IUPACData.extended_protein_letters</c> = "ACDEFGHIKLMNPQRSTVWYBXZJUO").
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
        /// Non-standard (genetically encoded by stop-codon recoding) amino acids with their
        /// three-letter symbols (Biopython <c>IUPACData.protein_letters_1to3_extended</c>):
        /// U = Sec (selenocysteine, UGA recoding), O = Pyl (pyrrolysine, UAG recoding).
        /// Their average masses are in <see cref="ProteinPhysicochemistry.AverageAminoAcidMasses"/>
        /// (Biopython <c>IUPACData.protein_weights</c>). <see cref="Translator"/> never emits them
        /// (standard-table translation yields '*' for UGA/UAG).
        /// </summary>
        public static readonly IReadOnlyDictionary<char, string> NonStandardAminoAcids = new Dictionary<char, string>
        {
            ['U'] = "Sec",
            ['O'] = "Pyl",
        };

        /// <summary>
        /// All valid characters: Biopython <c>IUPACData.extended_protein_letters</c>
        /// ("ACDEFGHIKLMNPQRSTVWYBXZJUO": 20 standard + ambiguous B/Z/J + unknown X + Sec U / Pyl O)
        /// plus the stop symbol '*' — 27 characters.
        /// </summary>
        public static readonly IReadOnlySet<char> ValidCharacters = new HashSet<char>(
            StandardAminoAcids.Concat(new[] { '*', 'X' }).Concat(AmbiguousAminoAcids.Keys).Concat(NonStandardAminoAcids.Keys)
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

            // ASCII-only case folding (TRANS-PROT-001): ToUpperInvariant folds U+017F 'ſ' → 'S',
            // which would accept "Mſ" as "MS"; scikit-bio Protein rejects non-ASCII letters.
            var normalized = SequenceExtensions.ToUpperAscii(sequence);
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
        /// <c>SequenceStatistics.CalculateMolecularWeight</c>). U (Sec, 168.0532 Da) and O (Pyl,
        /// 255.3134 Da) are weighed as in Biopython; ambiguous residues (B, Z, J, X) and '*'
        /// carry no mass and are skipped. Returns 0 for an empty sequence.
        /// </summary>
        public double MolecularWeight() => ProteinPhysicochemistry.MolecularWeight(_sequence);

        /// <summary>
        /// Calculates the theoretical isoelectric point (pI) on the EMBOSS <c>iep</c> pK scale
        /// (Epk.dat: N-terminus 7.5, C-terminus 3.6, C 8.5, D 3.9, E 4.1, H 6.5, K 10.8, R 12.5, Y 10.1),
        /// located by bisection to 1e-9 pH and rounded to two decimals. Delegates to the canonical
        /// <see cref="ProteinPhysicochemistry.IsoelectricPoint"/> (shared with
        /// <c>SequenceStatistics.CalculateIsoelectricPoint</c>). U (Sec) and O (Pyl) have no pK in
        /// Epk.dat and are treated as neutral, exactly as EMBOSS 6.6 <c>iep</c> and Biopython
        /// <c>IsoelectricPoint</c> do. Returns 0 for an empty sequence.
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
        /// <c>SequenceStatistics.CalculateHydrophobicity</c>); B, Z, J, X, U, O and '*' have no
        /// Kyte–Doolittle value and are skipped (Biopython raises <c>KeyError</c> on them).
        /// </summary>
        public double Gravy() => ProteinPhysicochemistry.Gravy(_sequence);

        /// <summary>
        /// Residue-class membership of EMBOSS <c>pepstats</c> (EMBOSS 6.6.0 <c>emboss/pepstats.c</c>,
        /// <c>pepstats_issize</c>), read from the property columns of the EMBOSS data file
        /// <c>Eamino.dat</c> (Tiny, Small, Aliphatic, Aromatic, Non-polar, Polar, Charge≠0, +ve, −ve).
        /// <c>Eamino.dat</c> copies U (Sec) from C, O (Pyl) from K and J (Xle) from I/L, so those letters
        /// are members although pepstats' printed labels (e.g. "Basic (H+K+R)") omit them. X and '*' belong
        /// to no class. B and Z are listed, but pepstats never counts them as such: they are first
        /// re-assigned to D/N and E/Q (see <see cref="TypeCount"/>).
        /// </summary>
        private static readonly IReadOnlyDictionary<AminoAcidType, string> PepstatsClassMembers =
            new Dictionary<AminoAcidType, string>
            {
                [AminoAcidType.Tiny] = "ACGSTU",
                [AminoAcidType.Small] = "ABCDGNPSTUV",
                [AminoAcidType.Aliphatic] = "AIJLV",
                [AminoAcidType.Aromatic] = "FHWY",
                [AminoAcidType.Nonpolar] = "ACFGIJLMPUVWY",
                [AminoAcidType.Polar] = "DEHKNOQRSTZ",
                [AminoAcidType.Charged] = "BDEHKORZ",
                [AminoAcidType.Basic] = "HKOR",
                [AminoAcidType.Acidic] = "BDEZ",
            };

        /// <summary>
        /// Number of residues in residue class <paramref name="type"/>, exactly as the "Number" column of the
        /// EMBOSS <c>pepstats</c> "Property / Residues / Number / Mole%" table.
        /// <para>
        /// Residue letters are tallied as in EMBOSS <c>embIepCompC</c> (<c>nucleus/embiep.c</c>): the B count
        /// is split by Dayhoff frequency, <c>j = (int)(0.5 + nB·5.5/9.8)</c> to D and <c>nB − j</c> to N; the Z
        /// count likewise, <c>j = (int)(0.5 + nZ·6.0/9.9)</c> to E and <c>nZ − j</c> to Q. The class sum then runs
        /// over the <c>Eamino.dat</c> members (U counts as C, O as K, J as I/L; X and '*' in no class).
        /// </para>
        /// </summary>
        /// <param name="type">Residue class (any <see cref="AminoAcidType"/> member is a pepstats class).</param>
        public int TypeCount(AminoAcidType type)
        {
            if (!PepstatsClassMembers.TryGetValue(type, out var members))
                throw new ArgumentOutOfRangeException(nameof(type), type, "Unknown residue class.");

            var counts = new int[26];
            foreach (char aa in _sequence)
                if (aa is >= 'A' and <= 'Z')
                    counts[aa - 'A']++;

            int nB = counts['B' - 'A'];
            if (nB > 0)
            {
                int toD = (int)(0.5 + nB * 5.5 / 9.8);
                counts['D' - 'A'] += toD;
                counts['N' - 'A'] += nB - toD;
                counts['B' - 'A'] = 0;
            }

            int nZ = counts['Z' - 'A'];
            if (nZ > 0)
            {
                int toE = (int)(0.5 + nZ * 6.0 / 9.9);
                counts['E' - 'A'] += toE;
                counts['Q' - 'A'] += nZ - toE;
                counts['Z' - 'A'] = 0;
            }

            int sum = 0;
            foreach (char m in members)
                sum += counts[m - 'A'];
            return sum;
        }

        /// <summary>
        /// Mole percent of residue class <paramref name="type"/> per EMBOSS <c>pepstats</c>
        /// ("Property / Residues / Number / Mole%" table): <c>100 · <see cref="TypeCount"/> / len</c>,
        /// unrounded (pepstats prints it with 3 decimals).
        /// <para>
        /// The denominator <c>len</c> is the pepstats "Residues" value: the sequence length, where a single
        /// terminal '*' is not counted (the EMBOSS sequence reader drops it), while X, B, Z, J, U, O and any
        /// other '*' are counted (pepstats: "MKRHDE*" → 6 residues; "AK**" → 3; "AAXXHK" → 6).
        /// </para>
        /// <para>
        /// Behaviour change (review 2026-09, B02 F47): this method previously used an unsourced textbook
        /// 4-class polarity grouping (C and Y "polar"; charged residues not "polar"; U/O in no class), divided
        /// by the full length including a terminal '*', and rounded to 2 decimals. It is superseded by the
        /// pepstats definition for all classes: <see cref="AminoAcidType.Acidic"/> (D+E) and
        /// <see cref="AminoAcidType.Basic"/> (H+K+R, now +O) keep their membership;
        /// <see cref="AminoAcidType.Nonpolar"/> is pepstats "Non-polar" (adds C, Y, U, J);
        /// <see cref="AminoAcidType.Polar"/> is pepstats "Polar" (adds D, E, H, K, R, O; drops C, Y).
        /// The per-residue <see cref="AminoAcidProperties.Type"/> label in <see cref="Properties"/> is unchanged
        /// and is not used here.
        /// </para>
        /// </summary>
        /// <param name="type">Residue class (Tiny, Small, Aliphatic, Aromatic, Nonpolar, Polar, Charged, Basic, Acidic).</param>
        /// <returns>Mole% in [0, 100]; 0 for an empty sequence.</returns>
        public double TypePercentage(AminoAcidType type)
        {
            int count = TypeCount(type);
            int len = _sequence.Length;
            if (len > 0 && _sequence[len - 1] == '*') len--;
            if (len == 0) return 0;
            return 100.0 * count / len;
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
                else if (NonStandardAminoAcids.TryGetValue(aa, out var nonStandard))
                {
                    if (sb.Length > 0) sb.Append('-');
                    sb.Append(nonStandard);
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

            return SuffixTree.FindAllOccurrences(SequenceExtensions.ToUpperAscii(pattern)).OrderBy(p => p);
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
                        $"Invalid amino acid '{c}' at position {i}. Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y (and *, X, B, Z, J, U, O).",
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
    /// Amino-acid residue class. As an argument to <see cref="ProteinSequence.TypePercentage"/> /
    /// <see cref="ProteinSequence.TypeCount"/> every member is one of the nine (overlapping) EMBOSS
    /// <c>pepstats</c> residue classes, with membership from <c>Eamino.dat</c>.
    /// The first four members are also used as the single per-residue label
    /// <see cref="AminoAcidProperties.Type"/> in <see cref="ProteinSequence.Properties"/>
    /// (textbook grouping: Nonpolar A F G I L M P V W; Polar C N Q S T Y; Acidic D E; Basic H K R).
    /// </summary>
    public enum AminoAcidType
    {
        /// <summary>pepstats "Non-polar" (A+C+F+G+I+L+M+P+V+W+Y, plus U, J).</summary>
        Nonpolar,
        /// <summary>pepstats "Polar" (D+E+H+K+N+Q+R+S+T+Z, plus O).</summary>
        Polar,
        /// <summary>pepstats "Acidic" (B+D+E+Z).</summary>
        Acidic,
        /// <summary>pepstats "Basic" (H+K+R, plus O).</summary>
        Basic,
        /// <summary>pepstats "Tiny" (A+C+G+S+T, plus U).</summary>
        Tiny,
        /// <summary>pepstats "Small" (A+B+C+D+G+N+P+S+T+V, plus U).</summary>
        Small,
        /// <summary>pepstats "Aliphatic" (A+I+L+V, plus J).</summary>
        Aliphatic,
        /// <summary>pepstats "Aromatic" (F+H+W+Y).</summary>
        Aromatic,
        /// <summary>pepstats "Charged" (B+D+E+H+K+R+Z, plus O): Eamino.dat charge ≠ 0.</summary>
        Charged
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
