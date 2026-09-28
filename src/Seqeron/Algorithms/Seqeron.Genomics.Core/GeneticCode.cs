using System.Collections.Concurrent;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Represents a genetic code table for translating codons to amino acids.
    /// All NCBI translation tables (gc.prt, Version 4.6: tables 1–6, 9–16, 21–33) are
    /// supported and built verbatim from the NCBI <c>ncbieaa</c> / <c>sncbieaa</c> strings.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Source: NCBI genetic code table <c>gc.prt</c> Version 4.6 (Elzanowski &amp; Ostell,
    /// "The Genetic Codes", NCBI Taxonomy; distributed with the NCBI C++ Toolkit,
    /// <c>src/objects/seqfeat/gc.prt</c>). Biopython <c>Bio.Data.CodonTable</c> is a
    /// faithful copy of the same file.
    /// </para>
    /// <para>
    /// Stop codons are the codons marked '*' in either the <c>ncbieaa</c> (AAs) or the
    /// <c>sncbieaa</c> (Starts) line; start codons are those marked 'M' in <c>sncbieaa</c>.
    /// In tables 27, 28 and 31 some codons are context-dependent ("dual-coding": stop or
    /// amino acid); as in Biopython, <see cref="Translate"/> returns their amino acid while
    /// <see cref="IsStopCodon"/> reports them as stops.
    /// </para>
    /// <para>
    /// IUPAC-ambiguous codons (R, Y, S, W, K, M, B, D, H, V, N) are resolved exactly as in
    /// Biopython's ambiguous codon tables: the codon is expanded to all concrete codons;
    /// if all are stops the result is '*'; if stops and amino acids are mixed the result is
    /// 'X'; a single amino acid is returned as-is; otherwise the most specific IUPAC
    /// ambiguous amino-acid code covering the set is returned (B = D/N, Z = E/Q,
    /// J = I/L), falling back to 'X'. Example: GCN → A, TAR → *, RAY → B, TAN → X.
    /// </para>
    /// </remarks>
    public sealed class GeneticCode
    {
        // Codon order of the NCBI ncbieaa/sncbieaa strings: Base1..Base3 each iterate T, C, A, G.
        private const string NcbiBaseOrder = "UCAG";

        // IUPAC nucleotide codes (RNA alphabet; T is normalised to U before lookup) and the
        // concrete bases each one stands for (IUPAC-IUB 1984, Cornish-Bowden 1985 NAR 13:3021),
        // derived from the canonical IupacHelper code set and matcher (no private copy of the
        // ambiguity table), with T spelled U.
        private static readonly IReadOnlyDictionary<char, string> IupacExpansion = BuildIupacExpansion();

        private static Dictionary<char, string> BuildIupacExpansion()
        {
            const string concreteDnaBases = "ACGT";
            var expansion = new Dictionary<char, string>();
            for (char code = 'A'; code <= 'Z'; code++)
            {
                if (!IupacHelper.IsNucleotideCode(code))
                    continue;

                var bases = concreteDnaBases
                    .Where(b => IupacHelper.MatchesIupac(b, code))
                    .Select(ToRnaBase)
                    .ToArray();
                expansion[ToRnaBase(code)] = new string(bases);
            }

            return expansion;

            static char ToRnaBase(char c) => c == 'T' ? 'U' : c;
        }

        // IUPAC ambiguous amino-acid codes with two meanings (IUPAC-IUB JCBN 1984; J added for I/L),
        // tried before the catch-all 'X' (Biopython IUPACData.extended_protein_values).
        private static readonly (char Code, string Meaning)[] AmbiguousAminoAcids =
        {
            ('B', "DN"), ('J', "IL"), ('Z', "EQ"),
        };

        private readonly IReadOnlyDictionary<string, char> _codonTable;
        private readonly IReadOnlySet<string> _startCodons;
        private readonly IReadOnlySet<string> _stopCodons;
        private readonly ConcurrentDictionary<string, char> _ambiguousCache = new();

        /// <summary>
        /// Gets the name of this genetic code (first NCBI name of the table).
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Gets the NCBI table number.
        /// </summary>
        public int TableNumber { get; }

        /// <summary>
        /// Gets the codon to amino acid mapping (64 unambiguous RNA codons).
        /// </summary>
        public IReadOnlyDictionary<string, char> CodonTable => _codonTable;

        /// <summary>
        /// Gets the (unambiguous, RNA) start codons.
        /// </summary>
        public IReadOnlySet<string> StartCodons => _startCodons;

        /// <summary>
        /// Gets the (unambiguous, RNA) stop codons.
        /// </summary>
        public IReadOnlySet<string> StopCodons => _stopCodons;

        private GeneticCode(string name, int tableNumber,
            Dictionary<string, char> codonTable,
            HashSet<string> startCodons,
            HashSet<string> stopCodons)
        {
            Name = name;
            TableNumber = tableNumber;
            _codonTable = codonTable;
            _startCodons = startCodons;
            _stopCodons = stopCodons;
        }

        /// <summary>
        /// Translates a codon to an amino acid.
        /// </summary>
        /// <param name="codon">Three-letter codon (RNA or DNA, IUPAC ambiguity codes allowed).</param>
        /// <returns>
        /// Single-letter amino acid code, '*' for a stop codon, an IUPAC ambiguous amino-acid
        /// code (B, Z, J) or 'X' for an ambiguous codon (see class remarks).
        /// </returns>
        /// <exception cref="ArgumentException">The codon is not 3 characters or contains a non-IUPAC symbol.</exception>
        public char Translate(string codon)
        {
            if (string.IsNullOrEmpty(codon) || codon.Length != 3)
                throw new ArgumentException("Codon must be exactly 3 characters.", nameof(codon));

            var normalized = Normalize(codon);

            if (_codonTable.TryGetValue(normalized, out char aa))
                return aa;

            if (!IsIupacCodon(normalized))
                throw new ArgumentException($"Unknown codon: {codon}", nameof(codon));

            return _ambiguousCache.GetOrAdd(normalized, TranslateAmbiguous);
        }

        // Biopython Bio.Data.CodonTable.AmbiguousForwardTable + Seq._translate_str semantics.
        private char TranslateAmbiguous(string codon)
        {
            var aminoAcids = new HashSet<char>();
            bool anyStop = false;

            foreach (var concrete in Expand(codon))
            {
                char aa = _codonTable[concrete];
                if (aa == '*')
                    anyStop = true;
                else
                    aminoAcids.Add(aa);
            }

            if (anyStop)
                return aminoAcids.Count == 0 ? '*' : 'X'; // all stops → stop; possible stop → X

            if (aminoAcids.Count == 1)
                return aminoAcids.First();

            foreach (var (code, meaning) in AmbiguousAminoAcids)
            {
                if (aminoAcids.All(a => meaning.Contains(a)))
                    return code;
            }

            return 'X';
        }

        private static string Normalize(string codon) => codon.ToUpperInvariant().Replace('T', 'U');

        private static bool IsIupacCodon(string codon)
        {
            foreach (char c in codon)
            {
                if (!IupacExpansion.ContainsKey(c))
                    return false;
            }
            return true;
        }

        private static IEnumerable<string> Expand(string codon)
        {
            foreach (char b1 in IupacExpansion[codon[0]])
                foreach (char b2 in IupacExpansion[codon[1]])
                    foreach (char b3 in IupacExpansion[codon[2]])
                        yield return new string(new[] { b1, b2, b3 });
        }

        // True when the codon is valid IUPAC and every concrete codon it stands for is in the set
        // (Biopython CodonTable.list_ambiguous_codons semantics for ambiguous start/stop lists).
        private static bool AllExpansionsIn(string codon, IReadOnlySet<string> set)
        {
            if (string.IsNullOrEmpty(codon) || codon.Length != 3)
                return false;

            var normalized = Normalize(codon);
            if (set.Contains(normalized))
                return true;
            if (!IsIupacCodon(normalized))
                return false;

            return Expand(normalized).All(set.Contains);
        }

        /// <summary>
        /// Checks if a codon is a start codon. An IUPAC-ambiguous codon is a start codon only
        /// if every concrete codon it stands for is a start codon (e.g. YTG in table 1).
        /// </summary>
        public bool IsStartCodon(string codon) => AllExpansionsIn(codon, _startCodons);

        /// <summary>
        /// Checks if a codon is a stop codon. An IUPAC-ambiguous codon is a stop codon only
        /// if every concrete codon it stands for is a stop codon (e.g. TAR, TRA in table 1).
        /// </summary>
        public bool IsStopCodon(string codon) => AllExpansionsIn(codon, _stopCodons);

        /// <summary>
        /// Gets all (unambiguous, RNA) codons that encode a specific amino acid, in NCBI codon order.
        /// For '*' this returns the table's stop codons (including dual-coding stops of tables 27/28/31).
        /// </summary>
        public IEnumerable<string> GetCodonsForAminoAcid(char aminoAcid)
        {
            var upperAa = char.ToUpperInvariant(aminoAcid);
            foreach (var kvp in _codonTable)
            {
                if (kvp.Value == upperAa || (upperAa == '*' && _stopCodons.Contains(kvp.Key)))
                    yield return kvp.Key;
            }
        }

        #region NCBI Genetic Codes

        private sealed record NcbiTable(int Id, string Name, string Aas, string Starts);

        // NCBI gc.prt Version 4.6 — id, first name, ncbieaa (AAs), sncbieaa (Starts); verbatim.
        private static readonly NcbiTable[] NcbiTables =
        {
            new(1, "Standard",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**--*----M---------------M----------------------------"),
            new(2, "Vertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                "----------**--------------------MMMM----------**---M------------"),
            new(3, "Yeast Mitochondrial",
                "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**----------------------MM---------------M------------"),
            new(4, "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--MM------**-------M------------MMMM---------------M------------"),
            new(5, "Invertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                "---M------**--------------------MMMM---------------M------------"),
            new(6, "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
                "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------"),
            new(9, "Echinoderm Mitochondrial; Flatworm Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "----------**-----------------------M---------------M------------"),
            new(10, "Euplotid Nuclear",
                "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**-----------------------M----------------------------"),
            new(11, "Bacterial, Archaeal and Plant Plastid",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**--*----M------------MMMM---------------M------------"),
            new(12, "Alternative Yeast Nuclear",
                "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*----M---------------M----------------------------"),
            new(13, "Ascidian Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                "---M------**----------------------MM---------------M------------"),
            new(14, "Alternative Flatworm Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "-----------*-----------------------M----------------------------"),
            new(15, "Blepharisma Macronuclear",
                "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------*---*--------------------M----------------------------"),
            new(16, "Chlorophycean Mitochondrial",
                "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------*---*--------------------M----------------------------"),
            new(21, "Trematode Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "----------**-----------------------M---------------M------------"),
            new(22, "Scenedesmus obliquus Mitochondrial",
                "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "------*---*---*--------------------M----------------------------"),
            new(23, "Thraustochytrium Mitochondrial",
                "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--*-------**--*-----------------M--M---------------M------------"),
            new(24, "Rhabdopleuridae Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                "---M------**-------M---------------M---------------M------------"),
            new(25, "Candidate Division SR1 and Gracilibacteria",
                "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**-----------------------M---------------M------------"),
            new(26, "Pachysolen tannophilus Nuclear",
                "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*----M---------------M----------------------------"),
            new(27, "Karyorelict Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------"),
            new(28, "Condylostoma Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*--------------------M----------------------------"),
            new(29, "Mesodinium Nuclear",
                "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------"),
            new(30, "Peritrich Nuclear",
                "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------"),
            new(31, "Blastocrithidia Nuclear",
                "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**-----------------------M----------------------------"),
            new(32, "Balanophoraceae Plastid",
                "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------*---*----M------------MMMM---------------M------------"),
            new(33, "Cephalodiscidae Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                "---M-------*-------M---------------M---------------M------------"),
        };

        private static readonly IReadOnlyDictionary<int, GeneticCode> ByTableNumber =
            NcbiTables.ToDictionary(t => t.Id, Build);

        /// <summary>
        /// Gets the NCBI translation table numbers supported by <see cref="GetByTableNumber"/>
        /// (all tables of NCBI gc.prt Version 4.6).
        /// </summary>
        public static IReadOnlyList<int> SupportedTableNumbers { get; } =
            NcbiTables.Select(t => t.Id).ToArray();

        /// <summary>
        /// Standard genetic code (NCBI Table 1).
        /// Used by most organisms.
        /// </summary>
        public static GeneticCode Standard { get; } = ByTableNumber[1];

        /// <summary>
        /// Vertebrate mitochondrial genetic code (NCBI Table 2).
        /// </summary>
        public static GeneticCode VertebrateMitochondrial { get; } = ByTableNumber[2];

        /// <summary>
        /// Yeast mitochondrial genetic code (NCBI Table 3).
        /// </summary>
        public static GeneticCode YeastMitochondrial { get; } = ByTableNumber[3];

        /// <summary>
        /// Bacterial, archaeal and plant plastid genetic code (NCBI Table 11).
        /// </summary>
        public static GeneticCode BacterialPlastid { get; } = ByTableNumber[11];

        /// <summary>
        /// Gets a genetic code by NCBI table number (see <see cref="SupportedTableNumbers"/>).
        /// </summary>
        /// <exception cref="ArgumentException">The number is not an NCBI translation table.</exception>
        public static GeneticCode GetByTableNumber(int tableNumber) =>
            ByTableNumber.TryGetValue(tableNumber, out var code)
                ? code
                : throw new ArgumentException($"Unknown genetic code table: {tableNumber}", nameof(tableNumber));

        private static GeneticCode Build(NcbiTable t)
        {
            var table = new Dictionary<string, char>(64);
            var starts = new HashSet<string>();
            var stops = new HashSet<string>();

            int i = 0;
            foreach (char b1 in NcbiBaseOrder)
                foreach (char b2 in NcbiBaseOrder)
                    foreach (char b3 in NcbiBaseOrder)
                    {
                        var codon = new string(new[] { b1, b2, b3 });
                        table[codon] = t.Aas[i];
                        if (t.Aas[i] == '*' || t.Starts[i] == '*')
                            stops.Add(codon);
                        if (t.Starts[i] == 'M')
                            starts.Add(codon);
                        i++;
                    }

            return new GeneticCode(t.Name, t.Id, table, starts, stops);
        }

        #endregion
    }
}
