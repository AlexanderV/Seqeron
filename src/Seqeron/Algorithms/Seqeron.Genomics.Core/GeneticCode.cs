using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Represents a genetic code table for translating codons to amino acids.
    /// Supports multiple genetic codes (Standard, Mitochondrial, etc.).
    /// </summary>
    public sealed class GeneticCode
    {
        private readonly IReadOnlyDictionary<string, char> _codonTable;
        private readonly IReadOnlySet<string> _startCodons;
        private readonly IReadOnlySet<string> _stopCodons;

        /// <summary>
        /// Gets the name of this genetic code.
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Gets the NCBI table number.
        /// </summary>
        public int TableNumber { get; }

        /// <summary>
        /// Gets the codon to amino acid mapping.
        /// </summary>
        public IReadOnlyDictionary<string, char> CodonTable => _codonTable;

        /// <summary>
        /// Gets the start codons.
        /// </summary>
        public IReadOnlySet<string> StartCodons => _startCodons;

        /// <summary>
        /// Gets the stop codons.
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
        /// <param name="codon">Three-letter codon (RNA or DNA).</param>
        /// <returns>Single-letter amino acid code, or '*' for stop codon.</returns>
        public char Translate(string codon)
        {
            if (string.IsNullOrEmpty(codon) || codon.Length != 3)
                throw new ArgumentException("Codon must be exactly 3 characters.", nameof(codon));

            var normalized = codon.ToUpperInvariant().Replace('T', 'U');

            if (_codonTable.TryGetValue(normalized, out char aa))
                return aa;

            throw new ArgumentException($"Unknown codon: {codon}", nameof(codon));
        }

        /// <summary>
        /// Checks if a codon is a start codon.
        /// </summary>
        public bool IsStartCodon(string codon)
        {
            if (string.IsNullOrEmpty(codon) || codon.Length != 3)
                return false;

            var normalized = codon.ToUpperInvariant().Replace('T', 'U');
            return _startCodons.Contains(normalized);
        }

        /// <summary>
        /// Checks if a codon is a stop codon.
        /// </summary>
        public bool IsStopCodon(string codon)
        {
            if (string.IsNullOrEmpty(codon) || codon.Length != 3)
                return false;

            var normalized = codon.ToUpperInvariant().Replace('T', 'U');
            return _stopCodons.Contains(normalized);
        }

        /// <summary>
        /// Gets all codons that encode a specific amino acid.
        /// </summary>
        public IEnumerable<string> GetCodonsForAminoAcid(char aminoAcid)
        {
            var upperAa = char.ToUpperInvariant(aminoAcid);
            foreach (var kvp in _codonTable)
            {
                if (kvp.Value == upperAa)
                    yield return kvp.Key;
            }
        }

        #region Standard Genetic Codes

        /// <summary>
        /// Standard genetic code (NCBI Table 1).
        /// Used by most organisms.
        /// </summary>
        public static GeneticCode Standard { get; } = CreateStandard();

        /// <summary>
        /// Vertebrate mitochondrial genetic code (NCBI Table 2).
        /// </summary>
        public static GeneticCode VertebrateMitochondrial { get; } = CreateVertebrateMitochondrial();

        /// <summary>
        /// Yeast mitochondrial genetic code (NCBI Table 3).
        /// </summary>
        public static GeneticCode YeastMitochondrial { get; } = CreateYeastMitochondrial();

        /// <summary>
        /// Bacterial and plant plastid genetic code (NCBI Table 11).
        /// </summary>
        public static GeneticCode BacterialPlastid { get; } = CreateBacterialPlastid();

        /// <summary>
        /// Gets a genetic code by NCBI table number.
        /// </summary>
        public static GeneticCode GetByTableNumber(int tableNumber) => tableNumber switch
        {
            1 => Standard,
            2 => VertebrateMitochondrial,
            3 => YeastMitochondrial,
            11 => BacterialPlastid,
            _ => throw new ArgumentException($"Unknown genetic code table: {tableNumber}", nameof(tableNumber))
        };

        private static GeneticCode CreateStandard()
        {
            var table = CreateStandardCodonTable();
            var startCodons = new HashSet<string> { "AUG" };
            var stopCodons = new HashSet<string> { "UAA", "UAG", "UGA" };

            return new GeneticCode("Standard", 1, table, startCodons, stopCodons);
        }

        private static GeneticCode CreateVertebrateMitochondrial()
        {
            // Start with standard table and apply differences
            var table = CreateStandardCodonTable();

            // Vertebrate mitochondrial differences:
            table["AUA"] = 'M'; // Isoleucine → Methionine
            table["UGA"] = 'W'; // Stop → Tryptophan
            table["AGA"] = '*'; // Arginine → Stop
            table["AGG"] = '*'; // Arginine → Stop

            var startCodons = new HashSet<string> { "AUG", "AUA", "AUU", "AUC" };
            var stopCodons = new HashSet<string> { "UAA", "UAG", "AGA", "AGG" };

            return new GeneticCode("Vertebrate Mitochondrial", 2, table, startCodons, stopCodons);
        }

        private static GeneticCode CreateYeastMitochondrial()
        {
            // Start with standard table and apply differences
            var table = CreateStandardCodonTable();

            // Yeast mitochondrial differences:
            table["CUU"] = 'T'; // Leucine → Threonine
            table["CUC"] = 'T'; // Leucine → Threonine
            table["CUA"] = 'T'; // Leucine → Threonine
            table["CUG"] = 'T'; // Leucine → Threonine
            table["AUA"] = 'M'; // Isoleucine → Methionine
            table["UGA"] = 'W'; // Stop → Tryptophan

            var startCodons = new HashSet<string> { "AUG", "AUA" };
            var stopCodons = new HashSet<string> { "UAA", "UAG" };

            return new GeneticCode("Yeast Mitochondrial", 3, table, startCodons, stopCodons);
        }

        private static GeneticCode CreateBacterialPlastid()
        {
            // Same as standard genetic code, just different start codons
            var table = CreateStandardCodonTable();

            // Bacteria can use alternative start codons
            var startCodons = new HashSet<string> { "AUG", "GUG", "UUG" };
            var stopCodons = new HashSet<string> { "UAA", "UAG", "UGA" };

            return new GeneticCode("Bacterial, Archaeal and Plant Plastid", 11, table, startCodons, stopCodons);
        }

        /// <summary>
        /// Creates the standard codon table as a base for other genetic codes.
        /// </summary>
        private static Dictionary<string, char> CreateStandardCodonTable()
        {
            return new Dictionary<string, char>
            {
                // Phenylalanine (F)
                ["UUU"] = 'F',
                ["UUC"] = 'F',
                // Leucine (L)
                ["UUA"] = 'L',
                ["UUG"] = 'L',
                ["CUU"] = 'L',
                ["CUC"] = 'L',
                ["CUA"] = 'L',
                ["CUG"] = 'L',
                // Isoleucine (I)
                ["AUU"] = 'I',
                ["AUC"] = 'I',
                ["AUA"] = 'I',
                // Methionine (M) - Start
                ["AUG"] = 'M',
                // Valine (V)
                ["GUU"] = 'V',
                ["GUC"] = 'V',
                ["GUA"] = 'V',
                ["GUG"] = 'V',
                // Serine (S)
                ["UCU"] = 'S',
                ["UCC"] = 'S',
                ["UCA"] = 'S',
                ["UCG"] = 'S',
                ["AGU"] = 'S',
                ["AGC"] = 'S',
                // Proline (P)
                ["CCU"] = 'P',
                ["CCC"] = 'P',
                ["CCA"] = 'P',
                ["CCG"] = 'P',
                // Threonine (T)
                ["ACU"] = 'T',
                ["ACC"] = 'T',
                ["ACA"] = 'T',
                ["ACG"] = 'T',
                // Alanine (A)
                ["GCU"] = 'A',
                ["GCC"] = 'A',
                ["GCA"] = 'A',
                ["GCG"] = 'A',
                // Tyrosine (Y)
                ["UAU"] = 'Y',
                ["UAC"] = 'Y',
                // Stop codons (*)
                ["UAA"] = '*',
                ["UAG"] = '*',
                ["UGA"] = '*',
                // Histidine (H)
                ["CAU"] = 'H',
                ["CAC"] = 'H',
                // Glutamine (Q)
                ["CAA"] = 'Q',
                ["CAG"] = 'Q',
                // Asparagine (N)
                ["AAU"] = 'N',
                ["AAC"] = 'N',
                // Lysine (K)
                ["AAA"] = 'K',
                ["AAG"] = 'K',
                // Aspartic acid (D)
                ["GAU"] = 'D',
                ["GAC"] = 'D',
                // Glutamic acid (E)
                ["GAA"] = 'E',
                ["GAG"] = 'E',
                // Cysteine (C)
                ["UGU"] = 'C',
                ["UGC"] = 'C',
                // Tryptophan (W)
                ["UGG"] = 'W',
                // Arginine (R)
                ["CGU"] = 'R',
                ["CGC"] = 'R',
                ["CGA"] = 'R',
                ["CGG"] = 'R',
                ["AGA"] = 'R',
                ["AGG"] = 'R',
                // Glycine (G)
                ["GGU"] = 'G',
                ["GGC"] = 'G',
                ["GGA"] = 'G',
                ["GGG"] = 'G'
            };
        }

        #endregion
    }
}
