using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace Seqeron.Genomics.Infrastructure
{
    /// <summary>
    /// Read-only per-residue amino-acid scales for sliding-window profiles (ExPASy ProtScale /
    /// Biopython <c>ProteinAnalysis.protein_scale</c>), keyed by upper-case one-letter code of the 20
    /// standard residues. Values are Biopython 1.88 <c>Bio.SeqUtils.ProtParamData</c> verbatim
    /// (review 2026-09 B03 F52). Kyte–Doolittle stays canonical in
    /// <see cref="ProteinPhysicochemistry.KyteDoolittle"/>; <see cref="KyteDoolittle"/> only forwards to it.
    /// </summary>
    public static class ProteinScales
    {
        /// <summary>Kyte &amp; Doolittle (1982) hydropathy — the canonical
        /// <see cref="ProteinPhysicochemistry.KyteDoolittle"/> table (Biopython <c>kd</c>).</summary>
        public static IReadOnlyDictionary<char, double> KyteDoolittle => ProteinPhysicochemistry.KyteDoolittle;

        /// <summary>
        /// Hopp &amp; Woods (1981) hydrophilicity, Proc Natl Acad Sci USA 78:3824-3828 (Biopython <c>hw</c>;
        /// positive = hydrophilic).
        /// </summary>
        public static IReadOnlyDictionary<char, double> HoppWoods { get; } = Freeze(new Dictionary<char, double>
        {
            ['A'] = -0.5, ['C'] = -1.0, ['D'] = 3.0, ['E'] = 3.0, ['F'] = -2.5, ['G'] = 0.0, ['H'] = -0.5,
            ['I'] = -1.8, ['K'] = 3.0, ['L'] = -1.8, ['M'] = -1.3, ['N'] = 0.2, ['P'] = 0.0, ['Q'] = 0.2,
            ['R'] = 3.0, ['S'] = 0.3, ['T'] = -0.4, ['V'] = -1.5, ['W'] = -3.4, ['Y'] = -2.3,
        });

        /// <summary>
        /// Emini surface fractional probability (surface accessibility; Vergoten &amp; Theophanides,
        /// Biomolecular Structure and Dynamics p. 138, 1997) — Biopython <c>em</c>. Not the Eisenberg
        /// consensus hydrophobicity scale, whatever the two-letter key may suggest.
        /// </summary>
        public static IReadOnlyDictionary<char, double> Emini { get; } = Freeze(new Dictionary<char, double>
        {
            ['A'] = 0.815, ['C'] = 0.394, ['D'] = 1.283, ['E'] = 1.445, ['F'] = 0.695, ['G'] = 0.714, ['H'] = 1.18,
            ['I'] = 0.603, ['K'] = 1.545, ['L'] = 0.603, ['M'] = 0.714, ['N'] = 1.296, ['P'] = 1.236, ['Q'] = 1.348,
            ['R'] = 1.475, ['S'] = 1.115, ['T'] = 1.184, ['V'] = 0.606, ['W'] = 0.808, ['Y'] = 1.089,
        });

        /// <summary>
        /// Janin (1979) interior-to-surface transfer energy scale (Biopython <c>ja</c>).
        /// </summary>
        public static IReadOnlyDictionary<char, double> Janin { get; } = Freeze(new Dictionary<char, double>
        {
            ['A'] = 0.28, ['C'] = 0.97, ['D'] = -0.52, ['E'] = -1.01, ['F'] = 0.46, ['G'] = 0.43, ['H'] = -0.31,
            ['I'] = 0.6, ['K'] = -1.62, ['L'] = 0.6, ['M'] = 0.43, ['N'] = -0.55, ['P'] = -0.42, ['Q'] = -0.69,
            ['R'] = -1.14, ['S'] = -0.19, ['T'] = -0.32, ['V'] = 0.6, ['W'] = 0.29, ['Y'] = -0.15,
        });

        /// <summary>
        /// Normalized average flexibility parameters (B-values), Vihinen, Torkkila &amp; Riikonen (1994)
        /// Proteins 19(2):141-149 (Biopython <c>Flex</c>).
        /// </summary>
        public static IReadOnlyDictionary<char, double> Flexibility { get; } = Freeze(new Dictionary<char, double>
        {
            ['A'] = 0.984, ['C'] = 0.906, ['D'] = 1.068, ['E'] = 1.094, ['F'] = 0.915, ['G'] = 1.031, ['H'] = 0.95,
            ['I'] = 0.927, ['K'] = 1.102, ['L'] = 0.935, ['M'] = 0.952, ['N'] = 1.048, ['P'] = 1.049, ['Q'] = 1.037,
            ['R'] = 1.008, ['S'] = 1.046, ['T'] = 0.997, ['V'] = 0.931, ['W'] = 0.904, ['Y'] = 0.929,
        });

        private static IReadOnlyDictionary<char, double> Freeze(Dictionary<char, double> table) =>
            new ReadOnlyDictionary<char, double>(table);
    }
}
