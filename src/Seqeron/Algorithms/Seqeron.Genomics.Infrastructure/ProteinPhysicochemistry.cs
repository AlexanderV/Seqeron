namespace Seqeron.Genomics.Infrastructure
{
    /// <summary>
    /// pK set used by <see cref="ProteinPhysicochemistry.IsoelectricPoint"/> and
    /// <see cref="ProteinPhysicochemistry.NetCharge"/>.
    /// </summary>
    public enum ProteinPkaScale
    {
        /// <summary>
        /// EMBOSS <c>iep</c> default <c>Epk.dat</c> (EMBOSS 6.6.0): N-terminus 7.5, C-terminus 3.6,
        /// C 8.5, D 3.9, E 4.1, H 6.5, K 10.8, R 12.5, Y 10.1; ambiguity codes B/Z are split into
        /// D/N and E/Q by Dayhoff frequencies exactly as <c>embIepCompC</c> does.
        /// </summary>
        Emboss,

        /// <summary>
        /// Bjellqvist et al. 1993/1994 (ExPASy Compute pI/Mw; Biopython
        /// <c>Bio.SeqUtils.IsoelectricPoint</c>): side chains C 9.0, D 4.05, E 4.45, H 5.98,
        /// K 10.0, R 12.0, Y 10.0; N-terminus 7.5 unless the N-terminal residue is
        /// A 7.59 / M 7.0 / S 6.93 / P 8.36 / T 6.82 / V 7.44 / E 7.7; C-terminus 3.55 unless the
        /// C-terminal residue is D 4.55 / E 4.75.
        /// </summary>
        Bjellqvist
    }

    /// <summary>
    /// Layer-neutral home of the protein physico-chemical calculations (average molecular weight,
    /// isoelectric point / net charge, Kyte–Doolittle GRAVY) and their residue tables, so that
    /// <c>Seqeron.Genomics.Core.ProteinSequence</c> (which cannot reference the Analysis project) and
    /// <c>SequenceStatistics</c> share one implementation (review 2026-09, DUP_MAP §7; B03 requests
    /// R2/R5/R8). The formulas, tables and conventions are those of the validated
    /// <c>SequenceStatistics.CalculateMolecularWeight</c> / <c>CalculateIsoelectricPoint</c> /
    /// <c>CalculateNetCharge</c> / <c>CalculateHydrophobicity</c> (units SEQ-MW-001, SEQ-PI-001,
    /// SEQ-HYDRO-001), reproduced exactly so the Analysis methods can delegate here.
    /// </summary>
    public static class ProteinPhysicochemistry
    {
        #region Molecular weight

        /// <summary>
        /// Average isotopic mass of one water molecule (Da): Biopython <c>Bio.SeqUtils.molecular_weight</c>
        /// uses 18.0153 (Expasy FindMod "Other mass values": H2O = 18.01524).
        /// </summary>
        public const double AverageWaterMass = 18.0153;

        // Average molecular masses of the 20 standard free amino acids plus selenocysteine (U) and
        // pyrrolysine (O) (Da). Source: Biopython Bio/Data/IUPACData.py `protein_weights`
        // ("Mass data taken from PubChem").
        private static readonly Dictionary<char, double> AminoAcidWeights = new()
        {
            { 'A', 89.0932 },  { 'C', 121.1582 }, { 'D', 133.1027 }, { 'E', 147.1293 },
            { 'F', 165.1891 }, { 'G', 75.0666 },  { 'H', 155.1546 }, { 'I', 131.1729 },
            { 'K', 146.1876 }, { 'L', 131.1729 }, { 'M', 149.2113 }, { 'N', 132.1179 },
            { 'O', 255.3134 }, { 'P', 115.1305 }, { 'Q', 146.1445 }, { 'R', 174.201 },
            { 'S', 105.0926 }, { 'T', 119.1192 }, { 'U', 168.0532 }, { 'V', 117.1463 },
            { 'W', 204.2252 }, { 'Y', 181.1885 }
        };

        /// <summary>
        /// Average molecular masses (Da) of the free amino acids (Biopython <c>IUPACData.protein_weights</c>):
        /// the 20 standard residues plus U (Sec) and O (Pyl), keyed by upper-case one-letter code.
        /// </summary>
        public static IReadOnlyDictionary<char, double> AverageAminoAcidMasses => AminoAcidWeights;

        /// <summary>
        /// Average-isotopic molecular weight of a linear protein (Da): sum of the free amino-acid
        /// masses minus (n − 1) waters (Biopython <c>molecular_weight(seq, "protein")</c>; Expasy
        /// Compute pI/Mw). Case-insensitive; symbols without a mass (B, Z, J, X, '*', gaps, …) are
        /// skipped (Biopython raises <c>ValueError</c>). Returns 0 for null/empty input or when no
        /// residue has a mass.
        /// </summary>
        public static double MolecularWeight(string? proteinSequence)
        {
            if (string.IsNullOrEmpty(proteinSequence))
                return 0;

            double weight = 0;
            int residues = 0;

            foreach (char aa in proteinSequence.ToUpperInvariant())
            {
                if (AminoAcidWeights.TryGetValue(aa, out double aaWeight))
                {
                    weight += aaWeight;
                    residues++;
                }
            }

            if (residues == 0)
                return 0;

            // One water is lost per peptide bond; n residues form (n − 1) bonds.
            return weight - (residues - 1) * AverageWaterMass;
        }

        #endregion

        #region Isoelectric point

        // EMBOSS 6.6.0 emboss/data/Epk.dat (the file read by `iep`; embIepPkReadFile defaults
        // amino = 7.50). The Epk.dat listing on the EMBOSS iep web page (Amino 8.6) is stale — the
        // page's own worked outputs (LACI_ECOLI pI 6.8385) are only reproduced with Amino 7.5.
        private const double EmbossNTerminusPka = 7.5;
        private const double EmbossCTerminusPka = 3.6;

        // Ionizable side chains: pKa and sign (+1 basic, -1 acidic).
        private static readonly Dictionary<char, (double pKa, int charge)> EmbossSideChains = new()
        {
            { 'C', (8.5, -1) }, { 'D', (3.9, -1) }, { 'E', (4.1, -1) }, { 'Y', (10.1, -1) },
            { 'H', (6.5, 1) },  { 'K', (10.8, 1) }, { 'R', (12.5, 1) }
        };

        // Bjellqvist et al. 1993 (Electrophoresis 14:1023) / 1994 (Electrophoresis 15:529), as
        // implemented by Biopython Bio/SeqUtils/IsoelectricPoint.py and ExPASy Compute pI/Mw.
        private static readonly Dictionary<char, (double pKa, int charge)> BjellqvistSideChains = new()
        {
            { 'C', (9.0, -1) }, { 'D', (4.05, -1) }, { 'E', (4.45, -1) }, { 'Y', (10.0, -1) },
            { 'H', (5.98, 1) }, { 'K', (10.0, 1) },  { 'R', (12.0, 1) }
        };

        private const double BjellqvistNTerminusPka = 7.5;
        private const double BjellqvistCTerminusPka = 3.55;

        private static readonly Dictionary<char, double> BjellqvistNTerminalResiduePka = new()
        {
            { 'A', 7.59 }, { 'M', 7.0 }, { 'S', 6.93 }, { 'P', 8.36 }, { 'T', 6.82 }, { 'V', 7.44 }, { 'E', 7.7 }
        };

        private static readonly Dictionary<char, double> BjellqvistCTerminalResiduePka = new()
        {
            { 'D', 4.55 }, { 'E', 4.75 }
        };

        // Bisection window and convergence: the root is located to 1e-9 pH and then rounded, so the
        // returned value is the correctly rounded pI.
        private const double MinPh = 0.0;
        private const double MaxPh = 14.0;
        private const double PiBisectionPrecision = 1e-9;

        // pI returned for empty/null input (pI is undefined for a zero-length protein).
        private const double NeutralPhDefault = 7.0;

        private const int PiDecimalPlaces = 2;

        /// <summary>
        /// Theoretical isoelectric point: the pH at which the Henderson–Hasselbalch net charge
        /// (<see cref="NetCharge"/>) is zero, located by bisection over [0, 14] to 1e-9 pH and rounded
        /// to two decimals. EMBOSS scale = EMBOSS 6.6.0 <c>iep</c> defaults; Bjellqvist scale = Biopython
        /// <c>IsoelectricPoint.pi()</c> / ExPASy Compute pI (Biopython clamps its search to [4.05, 12]).
        /// Residues without an ionizable side chain are ignored. Null/empty returns 7.0.
        /// </summary>
        public static double IsoelectricPoint(string? proteinSequence, ProteinPkaScale scale = ProteinPkaScale.Emboss)
        {
            if (string.IsNullOrEmpty(proteinSequence))
                return NeutralPhDefault;

            var model = BuildChargeModel(proteinSequence, scale);

            double pHLow = MinPh;
            double pHHigh = MaxPh;

            while (pHHigh - pHLow > PiBisectionPrecision)
            {
                double pH = (pHLow + pHHigh) / 2.0;

                // Net charge is monotonically non-increasing in pH: positive ⇒ pI is higher.
                if (model.NetCharge(pH) > 0)
                    pHLow = pH;
                else
                    pHHigh = pH;
            }

            return Math.Round((pHLow + pHHigh) / 2.0, PiDecimalPlaces);
        }

        /// <summary>
        /// Net charge at a given pH (Henderson–Hasselbalch): basic groups (N-terminus, K, R, H)
        /// contribute +1/(1+10^(pH−pKa)), acidic groups (C-terminus, D, E, C, Y) −1/(1+10^(pKa−pH));
        /// each terminus is counted once (EMBOSS <c>embIepGetCharge</c>; Biopython <c>charge_at_pH</c>).
        /// Returns 0 for null/empty input.
        /// </summary>
        public static double NetCharge(string? proteinSequence, double pH, ProteinPkaScale scale = ProteinPkaScale.Emboss)
        {
            if (string.IsNullOrEmpty(proteinSequence))
                return 0.0;

            return BuildChargeModel(proteinSequence, scale).NetCharge(pH);
        }

        private readonly record struct ChargeModel(
            double NTerminusPka,
            double CTerminusPka,
            Dictionary<char, (double pKa, int charge)> SideChains,
            Dictionary<char, int> Counts)
        {
            public double NetCharge(double pH)
            {
                double charge = 1.0 / (1.0 + Math.Pow(10, pH - NTerminusPka));
                charge -= 1.0 / (1.0 + Math.Pow(10, CTerminusPka - pH));

                foreach (var (aa, count) in Counts)
                {
                    var (pKa, sign) = SideChains[aa];
                    if (sign > 0)
                        charge += count / (1.0 + Math.Pow(10, pH - pKa));   // basic group
                    else
                        charge -= count / (1.0 + Math.Pow(10, pKa - pH));   // acidic group
                }

                return charge;
            }
        }

        private static ChargeModel BuildChargeModel(string proteinSequence, ProteinPkaScale scale)
        {
            string upper = proteinSequence.ToUpperInvariant();
            var sideChains = scale == ProteinPkaScale.Bjellqvist ? BjellqvistSideChains : EmbossSideChains;

            var counts = new Dictionary<char, int>();
            int countB = 0, countZ = 0;
            foreach (char aa in upper)
            {
                if (sideChains.ContainsKey(aa))
                    counts[aa] = counts.GetValueOrDefault(aa) + 1;
                else if (aa == 'B')
                    countB++;
                else if (aa == 'Z')
                    countZ++;
            }

            if (scale == ProteinPkaScale.Bjellqvist)
            {
                // Terminal-residue-specific pKs (Bjellqvist 1994; Biopython _update_pKs_tables).
                double nPka = BjellqvistNTerminalResiduePka.GetValueOrDefault(upper[0], BjellqvistNTerminusPka);
                double cPka = BjellqvistCTerminalResiduePka.GetValueOrDefault(upper[^1], BjellqvistCTerminusPka);
                return new ChargeModel(nPka, cPka, sideChains, counts);
            }

            // EMBOSS embIepCompC: B = D or N, Z = E or Q, split by Dayhoff frequencies
            // (D 5.5 / N 4.3; E 6.0 / Q 3.9), rounding half up via (int)(0.5 + x).
            if (countB > 0)
            {
                int asp = (int)(0.5 + countB * 5.5 / 9.8);
                if (asp > 0) counts['D'] = counts.GetValueOrDefault('D') + asp;
            }
            if (countZ > 0)
            {
                int glu = (int)(0.5 + countZ * 6.0 / 9.9);
                if (glu > 0) counts['E'] = counts.GetValueOrDefault('E') + glu;
            }

            return new ChargeModel(EmbossNTerminusPka, EmbossCTerminusPka, sideChains, counts);
        }

        #endregion

        #region Hydropathy

        // Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132; identical to Biopython
        // Bio.SeqUtils.ProtParamData.kd.
        private static readonly Dictionary<char, double> KyteDoolittleScale = new()
        {
            { 'A', 1.8 },  { 'R', -4.5 }, { 'N', -3.5 }, { 'D', -3.5 },
            { 'C', 2.5 },  { 'E', -3.5 }, { 'Q', -3.5 }, { 'G', -0.4 },
            { 'H', -3.2 }, { 'I', 4.5 },  { 'L', 3.8 },  { 'K', -3.9 },
            { 'M', 1.9 },  { 'F', 2.8 },  { 'P', -1.6 }, { 'S', -0.8 },
            { 'T', -0.7 }, { 'W', -0.9 }, { 'Y', -1.3 }, { 'V', 4.2 }
        };

        /// <summary>Kyte–Doolittle (1982) hydropathy index of the 20 standard residues.</summary>
        public static IReadOnlyDictionary<char, double> KyteDoolittle => KyteDoolittleScale;

        /// <summary>
        /// Grand average of hydropathy (GRAVY; Kyte &amp; Doolittle 1982; Biopython
        /// <c>ProteinAnalysis.gravy()</c>): mean Kyte–Doolittle value over the recognised residues,
        /// unrounded. Case-insensitive; residues outside the 20 standard ones are skipped (Biopython
        /// raises <c>KeyError</c>). Returns 0 for null/empty input or no recognised residue.
        /// </summary>
        public static double Gravy(string? proteinSequence)
        {
            if (string.IsNullOrEmpty(proteinSequence))
                return 0;

            double sum = 0;
            int count = 0;

            foreach (char aa in proteinSequence.ToUpperInvariant())
            {
                if (KyteDoolittleScale.TryGetValue(aa, out double value))
                {
                    sum += value;
                    count++;
                }
            }

            return count > 0 ? sum / count : 0;
        }

        #endregion
    }
}
