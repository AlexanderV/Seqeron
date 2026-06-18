using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the sequence-statistics algorithm group (<see cref="SequenceStatistics"/>).
///
/// This file is the single home for the Statistics block of checklist 01 (rows #121–#130, and the
/// profile/codon rows #227, #232, #234). Each test unit lives in its own <c>#region</c>. Oracles are
/// derived INDEPENDENTLY from the cited theory/doc, never routed through the production result, so a
/// self-consistent-but-wrong production formula is still caught.
///
/// Test Units: SEQ-COMPOSITION-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Statistics")]
public class SequenceStatisticsProperties
{
    #region SEQ-COMPOSITION-001 — Nucleotide Composition (counts, GC/AT fractions, skews)

    // -------------------------------------------------------------------------
    // Theory (nucleotide composition):
    //   • Every character is tallied into exactly one bucket (A/T/G/C/U/N/Other), so the seven
    //     counts sum to the sequence length.                                          (P counts = length)
    //   • Over the canonical bases (total = A+T+G+C+U): GcContent = (G+C)/total and
    //     AtContent = (A+T+U)/total, which partition the canonical bases ⇒ sum to 1.   (P Σ fractions = 1)
    //   • GcContent, AtContent ∈ [0,1]; GcSkew = (G−C)/(G+C), AtSkew = (A−T)/(A+T) ∈ [−1,1]. (R)
    //
    // The bucketing and the fraction formulae are recomputed independently here.
    // -------------------------------------------------------------------------

    private const double CompTolerance = 1e-9;

    /// <summary>A string over a mixed alphabet: canonical bases, U/N, lower-case, and non-base symbols.</summary>
    private static Arbitrary<string> MixedSequenceArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'U', 'N', 'a', 'c', 'g', 't', 'u', 'n', 'X', '-', '5')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// P (checklist "counts sum = length"): the seven nucleotide buckets (A, T, G, C, U, N, Other) account
    /// for every character exactly once, so they sum to the sequence length and each is non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_CountsPartitionTheSequence()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
            int sum = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU + comp.CountN + comp.CountOther;
            bool nonNeg = comp.CountA >= 0 && comp.CountT >= 0 && comp.CountG >= 0 && comp.CountC >= 0
                          && comp.CountU >= 0 && comp.CountN >= 0 && comp.CountOther >= 0;
            return (sum == comp.Length && comp.Length == seq.Length && nonNeg)
                .Label($"counts sum {sum} ≠ length {comp.Length} (seq len {seq.Length})");
        });
    }

    /// <summary>
    /// R (checklist "each fraction ∈ [0,1]") + P (checklist "Σ fractions = 1.0"): GcContent and AtContent lie
    /// in [0,1] and sum to 1 whenever the sequence has a canonical base (else both 0); the GC/AT skews lie in
    /// [−1,1]. GcContent matches the independent (G+C)/(A+T+G+C+U) oracle.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_FractionsInUnitRange_AndSumToOne()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
            int canonical = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU;

            bool inRange = comp.GcContent is >= 0.0 and <= 1.0 && comp.AtContent is >= 0.0 and <= 1.0
                           && comp.GcSkew is >= -1.0 and <= 1.0 && comp.AtSkew is >= -1.0 and <= 1.0;

            bool sumOk;
            bool gcOracleOk;
            if (canonical > 0)
            {
                sumOk = Math.Abs(comp.GcContent + comp.AtContent - 1.0) < CompTolerance;
                double oracleGc = (double)(comp.CountG + comp.CountC) / canonical;
                gcOracleOk = Math.Abs(comp.GcContent - oracleGc) < CompTolerance;
            }
            else
            {
                sumOk = comp.GcContent == 0.0 && comp.AtContent == 0.0;
                gcOracleOk = true;
            }

            return (inRange && sumOk && gcOracleOk)
                .Label($"GC={comp.GcContent}, AT={comp.AtContent}, canonical={canonical}");
        });
    }

    /// <summary>D (determinism): nucleotide composition is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
            SequenceStatistics.CalculateNucleotideComposition(seq)
                .Equals(SequenceStatistics.CalculateNucleotideComposition(seq))
                .Label("CalculateNucleotideComposition is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: a GC-only sequence is 100% GC / 0% AT; a mixed sequence partitions correctly; the empty
    /// sequence yields all-zero counts and fractions; non-base symbols fall into Other.
    /// </summary>
    [Test]
    [Category("Property")]
    public void NucleotideComposition_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var gc = SequenceStatistics.CalculateNucleotideComposition("GCGCGC");
            Assert.That(gc.GcContent, Is.EqualTo(1.0).Within(CompTolerance), "All G/C ⇒ GC 100%.");
            Assert.That(gc.AtContent, Is.EqualTo(0.0).Within(CompTolerance), "No A/T/U ⇒ AT 0%.");

            var mixed = SequenceStatistics.CalculateNucleotideComposition("AATTGGCCNNXX");
            Assert.That(mixed.CountA + mixed.CountT + mixed.CountG + mixed.CountC + mixed.CountU + mixed.CountN + mixed.CountOther,
                Is.EqualTo(12), "All 12 characters are bucketed.");
            Assert.That(mixed.CountN, Is.EqualTo(2), "Two N bases.");
            Assert.That(mixed.CountOther, Is.EqualTo(2), "Two non-base symbols (X) ⇒ Other.");
            Assert.That(mixed.GcContent + mixed.AtContent, Is.EqualTo(1.0).Within(CompTolerance), "GC% + AT% = 1.");

            var empty = SequenceStatistics.CalculateNucleotideComposition("");
            Assert.That(empty.Length, Is.EqualTo(0));
            Assert.That(empty.GcContent, Is.EqualTo(0.0));
        });
    }

    #endregion

    #region SEQ-DINUC-001 — Dinucleotide Frequencies (Karlin genomic signature)

    // -------------------------------------------------------------------------
    // Theory (Karlin genomic signature, PMC126251):
    //   • f_XY = count(XY) / total, total = number of ATGCU-only adjacent dinucleotides.
    //   • Each frequency ≥ 0 (in (0,1]); the frequencies sum to 1 when any valid dinucleotide exists.
    //   • For a pure ATGCU sequence of length N every adjacent pair is valid, so total = N−1
    //     (the "Σ dinucleotide counts = length−1" invariant).
    //
    // The adjacent-pair counts and the normalization are recomputed independently.
    // -------------------------------------------------------------------------

    private static Arbitrary<string> PureDnaRnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'U')
            .ArrayOf()
            .Where(a => a.Length >= 2)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// R (checklist "each frequency ≥ 0") + P (sum to 1): every dinucleotide frequency is in (0,1] and the
    /// frequencies sum to 1 whenever the sequence contains at least one ATGCU dinucleotide; each key is a
    /// two-character ATGCU dinucleotide. (Karlin genomic signature)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DinucleotideFrequencies_NonNegative_AndSumToOne()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

            bool valuesOk = freq.Values.All(f => f is > 0.0 and <= 1.0 + CompTolerance);
            bool keysOk = freq.Keys.All(k => k.Length == 2 && k.All(c => "ATGCU".Contains(c)));
            bool sumOk = freq.Count == 0 || Math.Abs(freq.Values.Sum() - 1.0) < CompTolerance;
            return (valuesOk && keysOk && sumOk)
                .Label($"dinuc freqs sum {freq.Values.Sum()} over {freq.Count} keys");
        });
    }

    /// <summary>
    /// P (checklist "Σ dinucleotide counts = length−1"): for a pure ATGCU sequence of length N every adjacent
    /// pair is a valid dinucleotide, so the implied counts (f_XY × total) sum to N−1 and each frequency equals
    /// the independent adjacent-pair-count / (N−1) oracle. (Karlin genomic signature)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DinucleotideFrequencies_PureSequence_MatchesAdjacentPairOracle()
    {
        return Prop.ForAll(PureDnaRnaArbitrary(), seq =>
        {
            var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

            var counts = new Dictionary<string, int>();
            string upper = seq.ToUpperInvariant();
            for (int i = 0; i < upper.Length - 1; i++)
            {
                string d = upper.Substring(i, 2);
                counts[d] = counts.GetValueOrDefault(d) + 1;
            }

            int total = upper.Length - 1; // every adjacent pair is valid in a pure ATGCU sequence
            bool keysMatch = freq.Count == counts.Count;
            bool valuesMatch = counts.All(kv =>
                freq.TryGetValue(kv.Key, out double f) && Math.Abs(f - (double)kv.Value / total) < CompTolerance);
            bool countsSumToLengthMinusOne = counts.Values.Sum() == total;

            return (keysMatch && valuesMatch && countsSumToLengthMinusOne)
                .Label($"freq keys {freq.Count} vs counts {counts.Count}; total {total} (len {seq.Length})");
        });
    }

    /// <summary>D (determinism): dinucleotide frequencies are identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property DinucleotideFrequencies_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var a = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
            var b = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
            return (a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out double v) && v == kv.Value))
                .Label("CalculateDinucleotideFrequencies is not deterministic for identical input");
        });
    }

    /// <summary>
    /// Anchors: "AAAA" ⇒ a single AA dinucleotide at frequency 1.0 (3 of 3 positions); a sequence shorter
    /// than 2 yields no dinucleotides.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DinucleotideFrequencies_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var aaaa = SequenceStatistics.CalculateDinucleotideFrequencies("AAAA");
            Assert.That(aaaa["AA"], Is.EqualTo(1.0).Within(CompTolerance), "All three dinucleotides are AA ⇒ freq 1.");
            Assert.That(aaaa.Values.Sum(), Is.EqualTo(1.0).Within(CompTolerance));

            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies("A"), Is.Empty, "Length < 2 ⇒ no dinucleotides.");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(""), Is.Empty);
        });
    }

    #endregion

    #region SEQ-HYDRO-001 — Hydrophobicity (Kyte-Doolittle GRAVY)

    // -------------------------------------------------------------------------
    // Theory (Kyte & Doolittle 1982; Biopython ProtParam kd):
    //   • GRAVY = mean Kyte-Doolittle hydropathy over the recognized residues.
    //   • Each KD value ∈ [−4.5 (R), 4.5 (I)], so the mean ∈ [−4.5, 4.5] and is finite.   (R)
    //   • Appending the most hydrophobic residue (I=4.5) cannot lower GRAVY; appending the
    //     least (R=−4.5) cannot raise it (more hydrophobic residues → higher mean).        (M)
    //
    // The KD scale and the mean are reconstructed independently here.
    // -------------------------------------------------------------------------

    private static readonly Dictionary<char, double> KdScaleOracle = new()
    {
        ['A'] = 1.8, ['R'] = -4.5, ['N'] = -3.5, ['D'] = -3.5, ['C'] = 2.5, ['E'] = -3.5, ['Q'] = -3.5,
        ['G'] = -0.4, ['H'] = -3.2, ['I'] = 4.5, ['L'] = 3.8, ['K'] = -3.9, ['M'] = 1.9, ['F'] = 2.8,
        ['P'] = -1.6, ['S'] = -0.8, ['T'] = -0.7, ['W'] = -0.9, ['Y'] = -1.3, ['V'] = 4.2,
    };

    /// <summary>Protein strings over the 20 standard residues plus non-standard symbols (B/Z/X, skipped).</summary>
    private static Arbitrary<string> ProteinArbitrary() =>
        Gen.Elements("ARNDCEQGHILKMFPSTWYVBZX".ToCharArray())
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    private static double OracleGravy(string seq)
    {
        double sum = 0;
        int count = 0;
        foreach (char aa in seq.ToUpperInvariant())
        {
            if (KdScaleOracle.TryGetValue(aa, out double v))
            {
                sum += v;
                count++;
            }
        }

        return count > 0 ? sum / count : 0.0;
    }

    /// <summary>
    /// R (checklist "score finite within scale range"): GRAVY equals the independent mean Kyte-Doolittle
    /// hydropathy over the recognized residues, lies in [−4.5, 4.5], and is finite (0 when no residue is
    /// recognized). (Kyte & Doolittle 1982)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hydrophobicity_EqualsMeanKyteDoolittle_InScaleRange()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            double gravy = SequenceStatistics.CalculateHydrophobicity(seq);
            double oracle = OracleGravy(seq);
            return (Math.Abs(gravy - oracle) < CompTolerance && double.IsFinite(gravy)
                    && gravy is >= -4.5 - CompTolerance and <= 4.5 + CompTolerance)
                .Label($"GRAVY {gravy} ≠ oracle {oracle}");
        });
    }

    /// <summary>
    /// M (checklist "more hydrophobic residues → higher mean"): appending the most hydrophobic residue (I,
    /// 4.5) never decreases GRAVY, and appending the least hydrophobic (R, −4.5) never increases it — a mean
    /// moves toward any appended value beyond its current level. (Kyte & Doolittle 1982)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hydrophobicity_AppendingExtremeResidue_IsMonotone()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            double baseGravy = SequenceStatistics.CalculateHydrophobicity(seq);
            double withIle = SequenceStatistics.CalculateHydrophobicity(seq + "I");
            double withArg = SequenceStatistics.CalculateHydrophobicity(seq + "R");
            return (withIle >= baseGravy - CompTolerance && withArg <= baseGravy + CompTolerance)
                .Label($"base {baseGravy}, +I {withIle}, +R {withArg}");
        });
    }

    /// <summary>D (determinism): GRAVY is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property Hydrophobicity_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
            (SequenceStatistics.CalculateHydrophobicity(seq) == SequenceStatistics.CalculateHydrophobicity(seq))
                .Label("CalculateHydrophobicity is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: all-Ile is the scale maximum 4.5, all-Arg the minimum −4.5, "AAAA" is 1.8, and an empty or
    /// fully non-standard sequence is 0. (Kyte & Doolittle 1982)
    /// </summary>
    [Test]
    [Category("Property")]
    public void Hydrophobicity_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateHydrophobicity("IIII"), Is.EqualTo(4.5).Within(CompTolerance), "All Ile ⇒ max 4.5.");
            Assert.That(SequenceStatistics.CalculateHydrophobicity("RRRR"), Is.EqualTo(-4.5).Within(CompTolerance), "All Arg ⇒ min −4.5.");
            Assert.That(SequenceStatistics.CalculateHydrophobicity("AAAA"), Is.EqualTo(1.8).Within(CompTolerance), "Ala = 1.8.");
            Assert.That(SequenceStatistics.CalculateHydrophobicity(""), Is.EqualTo(0.0), "Empty ⇒ 0.");
            Assert.That(SequenceStatistics.CalculateHydrophobicity("BZX"), Is.EqualTo(0.0), "No recognized residue ⇒ 0.");
        });
    }

    #endregion

    #region SEQ-MW-001 — Protein Molecular Weight (Expasy / Biopython average mass)

    // -------------------------------------------------------------------------
    // Theory (Expasy Compute pI/Mw; Biopython protein_weights):
    //   • MW = Σ(free residue masses) − (n − 1)·water, removing one water per peptide bond.
    //   • MW > 0 for any sequence with a recognized residue; 0 for empty / all-unknown.   (R)
    //   • Appending a residue raises MW by (residue mass − water) > 0.                     (M longer → higher)
    //   • Concatenation: MW(A+B) = MW(A) + MW(B) − water (one extra joining bond).         (P additive)
    //
    // The residue masses and the water-loss correction are reconstructed independently.
    // -------------------------------------------------------------------------

    private const double WaterMassOracle = 18.0153;

    private static readonly Dictionary<char, double> AaWeightOracle = new()
    {
        ['A'] = 89.0932, ['C'] = 121.1582, ['D'] = 133.1027, ['E'] = 147.1293, ['F'] = 165.1891,
        ['G'] = 75.0666, ['H'] = 155.1546, ['I'] = 131.1729, ['K'] = 146.1876, ['L'] = 131.1729,
        ['M'] = 149.2113, ['N'] = 132.1179, ['P'] = 115.1305, ['Q'] = 146.1445, ['R'] = 174.201,
        ['S'] = 105.0926, ['T'] = 119.1192, ['V'] = 117.1463, ['W'] = 204.2252, ['Y'] = 181.1885,
    };

    private static (double mw, int residues) OracleMolecularWeight(string seq)
    {
        double sum = 0;
        int n = 0;
        foreach (char aa in seq.ToUpperInvariant())
        {
            if (AaWeightOracle.TryGetValue(aa, out double w))
            {
                sum += w;
                n++;
            }
        }

        return n == 0 ? (0.0, 0) : (sum - (n - 1) * WaterMassOracle, n);
    }

    /// <summary>
    /// R (checklist "MW &gt; 0 for non-empty") + P (additive over residues): the molecular weight equals the
    /// independent Σ(residue masses) − (n−1)·water and is strictly positive when at least one residue is
    /// recognized (0 otherwise). (Expasy / Biopython)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeight_EqualsResidueSumMinusPeptideWater()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            double mw = SequenceStatistics.CalculateMolecularWeight(seq);
            (double oracle, int residues) = OracleMolecularWeight(seq);
            bool formulaOk = Math.Abs(mw - oracle) <= 1e-6 * Math.Max(1.0, oracle);
            bool positivity = residues > 0 ? mw > 0.0 : mw == 0.0;
            return (formulaOk && positivity).Label($"MW {mw} ≠ oracle {oracle} (residues {residues})");
        });
    }

    /// <summary>
    /// M (checklist "longer sequence → higher MW"): appending a recognized residue strictly increases the
    /// molecular weight (each residue adds its mass minus one water, ≥ Gly − water &gt; 0). (Expasy / Biopython)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeight_AppendingResidue_StrictlyIncreases()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            double baseMw = SequenceStatistics.CalculateMolecularWeight(seq);
            double longer = SequenceStatistics.CalculateMolecularWeight(seq + "G");
            return (longer > baseMw).Label($"appending G did not increase MW: {baseMw} → {longer}");
        });
    }

    /// <summary>
    /// P (checklist "additive over residues"): concatenating two peptides each with a recognized residue adds
    /// their weights minus one water for the new joining bond — MW(A+B) = MW(A) + MW(B) − water. (Expasy)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeight_Concatenation_RemovesOneJoiningWater()
    {
        var arb = (from a in ProteinArbitrary().Generator
                   from b in ProteinArbitrary().Generator
                   select (a, b)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            (_, int nA) = OracleMolecularWeight(t.a);
            (_, int nB) = OracleMolecularWeight(t.b);
            if (nA == 0 || nB == 0)
            {
                return true.ToProperty(); // the joining-water relation needs a residue on each side
            }

            double mwA = SequenceStatistics.CalculateMolecularWeight(t.a);
            double mwB = SequenceStatistics.CalculateMolecularWeight(t.b);
            double mwAb = SequenceStatistics.CalculateMolecularWeight(t.a + t.b);
            return (Math.Abs(mwAb - (mwA + mwB - WaterMassOracle)) < 1e-6 * Math.Max(1.0, mwAb))
                .Label($"MW(A+B)={mwAb} ≠ MW(A)+MW(B)−water = {mwA + mwB - WaterMassOracle}");
        });
    }

    /// <summary>D (determinism): protein molecular weight is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeight_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
            (SequenceStatistics.CalculateMolecularWeight(seq) == SequenceStatistics.CalculateMolecularWeight(seq))
                .Label("CalculateMolecularWeight is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: glycine "G" is a free amino acid (75.0666 Da, no peptide bond); "GG" loses one water; empty
    /// and all-unknown sequences are 0. (Expasy / Biopython)
    /// </summary>
    [Test]
    [Category("Property")]
    public void MolecularWeight_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateMolecularWeight("G"), Is.EqualTo(75.0666).Within(1e-6), "Free Gly.");
            Assert.That(SequenceStatistics.CalculateMolecularWeight("GG"), Is.EqualTo(2 * 75.0666 - 18.0153).Within(1e-6),
                "Gly-Gly loses one water for the peptide bond.");
            Assert.That(SequenceStatistics.CalculateMolecularWeight(""), Is.EqualTo(0.0), "Empty ⇒ 0.");
            Assert.That(SequenceStatistics.CalculateMolecularWeight("BZX"), Is.EqualTo(0.0), "No recognized residue ⇒ 0.");
        });
    }

    #endregion

    #region SEQ-PI-001 — Isoelectric Point (EMBOSS pKa / Henderson-Hasselbalch)

    // -------------------------------------------------------------------------
    // Theory (EMBOSS iep; Henderson-Hasselbalch; Osorio 2015):
    //   • pI is the pH where the net charge is zero, found by bisection over [0,14].   (R pI ∈ [0,14])
    //   • Net charge is monotone decreasing in pH; at the returned pI the independent
    //     net-charge oracle is ≈ 0 (within the bisection bracket × max slope).          (P charge ≈ 0)
    //
    // The pKa scale, termini and the net-charge formula are reconstructed independently.
    // -------------------------------------------------------------------------

    private static readonly Dictionary<char, (double pKa, int sign)> IonizableOracle = new()
    {
        ['D'] = (3.9, -1), ['E'] = (4.1, -1), ['C'] = (8.5, -1), ['Y'] = (10.1, -1),
        ['H'] = (6.5, 1), ['K'] = (10.8, 1), ['R'] = (12.5, 1),
    };

    private const double NTerminusPkaOracle = 8.6;
    private const double CTerminusPkaOracle = 3.6;

    private static (double charge, int groups) OracleNetCharge(string seq, double pH)
    {
        double charge = 1.0 / (1.0 + Math.Pow(10, pH - NTerminusPkaOracle))
                        - 1.0 / (1.0 + Math.Pow(10, CTerminusPkaOracle - pH));
        int groups = 2; // both termini
        foreach (char aa in seq.ToUpperInvariant())
        {
            if (IonizableOracle.TryGetValue(aa, out var g))
            {
                groups++;
                if (g.sign > 0)
                {
                    charge += 1.0 / (1.0 + Math.Pow(10, pH - g.pKa));
                }
                else
                {
                    charge -= 1.0 / (1.0 + Math.Pow(10, g.pKa - pH));
                }
            }
        }

        return (charge, groups);
    }

    private static Arbitrary<string> BoundedProteinArbitrary(int maxLen) =>
        (from n in Gen.Choose(0, maxLen)
         from chars in Gen.Elements("ARNDCEQGHILKMFPSTWYV".ToCharArray()).ArrayOf(n)
         select new string(chars)).ToArbitrary();

    /// <summary>
    /// R (checklist "pI ∈ [0,14]"): the isoelectric point of any protein lies within the standard pH window
    /// [0, 14]. (EMBOSS iep bisection bounds)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPoint_InPhWindow()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);
            return (pi is >= 0.0 and <= 14.0).Label($"pI {pi} outside [0,14]");
        });
    }

    /// <summary>
    /// P (checklist "net charge at pI ≈ 0"): the independent Henderson-Hasselbalch net charge evaluated at the
    /// returned pI is ≈ 0 — within the bisection bracket (0.01 pH) times the maximum charge slope
    /// (ln10/4 per ionizable group). (EMBOSS iep / Osorio 2015)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPoint_NetChargeIsApproximatelyZero()
    {
        return Prop.ForAll(BoundedProteinArbitrary(25), seq =>
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);
            if (string.IsNullOrEmpty(seq))
            {
                return true.ToProperty(); // empty input returns the neutral 7.0 sentinel, not a charge zero
            }

            (double charge, int groups) = OracleNetCharge(seq, pi);
            // |charge(pI)| ≤ maxSlope·bracket; maxSlope = (ln10/4)·groups, bracket ≤ 0.01, +rounding margin.
            double tolerance = Math.Log(10) / 4.0 * groups * 0.02 + 0.05;
            return (Math.Abs(charge) < tolerance)
                .Label($"net charge {charge} at pI {pi} exceeds tolerance {tolerance} (groups {groups})");
        });
    }

    /// <summary>D (determinism): the isoelectric point is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPoint_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
            (SequenceStatistics.CalculateIsoelectricPoint(seq) == SequenceStatistics.CalculateIsoelectricPoint(seq))
                .Label("CalculateIsoelectricPoint is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: an acidic protein (poly-Asp) has a low pI (&lt; 7), a basic one (poly-Lys) a high pI (&gt; 7),
    /// and the empty sentinel is the neutral 7.0. (EMBOSS iep)
    /// </summary>
    [Test]
    [Category("Property")]
    public void IsoelectricPoint_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("DDDDDD"), Is.LessThan(7.0), "Acidic ⇒ low pI.");
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint("KKKKKK"), Is.GreaterThan(7.0), "Basic ⇒ high pI.");
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint(""), Is.EqualTo(7.0), "Empty ⇒ neutral sentinel 7.0.");
        });
    }

    #endregion

    #region SEQ-SECSTRUCT-001 — Secondary Structure Propensity (Chou-Fasman)

    // -------------------------------------------------------------------------
    // Theory (Chou & Fasman 1978):
    //   • Per-window mean conformational propensities Pa (helix), Pb (sheet), Pt (turn).
    //   • All published propensities are positive ⇒ every emitted propensity ≥ 0.        (R ≥ 0)
    //   • A pure-residue sequence of length N with window w yields N−w+1 windows, each of
    //     which assigns a dominant class H/E/C (argmax of the three propensities).         (P every residue → H/E/C)
    //
    // The Chou-Fasman table and the per-window mean are reconstructed independently.
    // -------------------------------------------------------------------------

    private static readonly Dictionary<char, (double Helix, double Sheet, double Turn)> ChouFasmanOracle = new()
    {
        ['A'] = (1.42, 0.83, 0.66), ['R'] = (0.98, 0.93, 0.95), ['N'] = (0.67, 0.89, 1.56), ['D'] = (1.01, 0.54, 1.46),
        ['C'] = (0.70, 1.19, 1.19), ['E'] = (1.51, 0.37, 0.74), ['Q'] = (1.11, 1.10, 0.98), ['G'] = (0.57, 0.75, 1.56),
        ['H'] = (1.00, 0.87, 0.95), ['I'] = (1.08, 1.60, 0.47), ['L'] = (1.21, 1.30, 0.59), ['K'] = (1.14, 0.74, 1.01),
        ['M'] = (1.45, 1.05, 0.60), ['F'] = (1.13, 1.38, 0.60), ['P'] = (0.57, 0.55, 1.52), ['S'] = (0.77, 0.75, 1.43),
        ['T'] = (0.83, 1.19, 0.96), ['W'] = (1.08, 1.37, 0.96), ['Y'] = (0.69, 1.47, 1.14), ['V'] = (1.06, 1.70, 0.50),
    };

    private static Arbitrary<(string seq, int window)> SecStructProblemArbitrary() =>
        (from n in Gen.Choose(1, 20)
         from chars in Gen.Elements("ARNDCEQGHILKMFPSTWYV".ToCharArray()).ArrayOf(n)
         from w in Gen.Choose(1, n)
         select (new string(chars), w)).ToArbitrary();

    /// <summary>
    /// R (checklist "each propensity ≥ 0") + formula: each window's (Helix, Sheet, Turn) equals the
    /// independent mean Chou-Fasman propensity over its residues, and every value is finite and ≥ 0.
    /// (Chou & Fasman 1978)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SecondaryStructure_PerWindowMeans_AreNonNegative_AndMatchOracle()
    {
        return Prop.ForAll(SecStructProblemArbitrary(), t =>
        {
            var windows = SequenceStatistics.PredictSecondaryStructure(t.seq, t.window).ToList();
            bool ok = true;
            for (int i = 0; ok && i < windows.Count; i++)
            {
                double h = 0, s = 0, turn = 0;
                for (int j = 0; j < t.window; j++)
                {
                    var p = ChouFasmanOracle[t.seq[i + j]];
                    h += p.Helix; s += p.Sheet; turn += p.Turn;
                }

                h /= t.window; s /= t.window; turn /= t.window;
                var w = windows[i];
                ok &= Math.Abs(w.Helix - h) < CompTolerance && Math.Abs(w.Sheet - s) < CompTolerance && Math.Abs(w.Turn - turn) < CompTolerance
                      && w.Helix >= 0 && w.Sheet >= 0 && w.Turn >= 0
                      && double.IsFinite(w.Helix) && double.IsFinite(w.Sheet) && double.IsFinite(w.Turn);
            }

            return ok.Label($"window means mismatch or negative ({windows.Count} windows, w={t.window})");
        });
    }

    /// <summary>
    /// P (checklist "every residue assigned H/E/C"): a pure-residue sequence of length N with window w yields
    /// exactly N−w+1 windows, and each window admits a dominant class (one of helix/sheet/turn is the maximum).
    /// (Chou & Fasman 1978 sliding window)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SecondaryStructure_WindowCount_AndClassAssignment()
    {
        return Prop.ForAll(SecStructProblemArbitrary(), t =>
        {
            var windows = SequenceStatistics.PredictSecondaryStructure(t.seq, t.window).ToList();
            bool countOk = windows.Count == t.seq.Length - t.window + 1;
            bool classOk = windows.All(w =>
            {
                double max = Math.Max(w.Helix, Math.Max(w.Sheet, w.Turn));
                return max == w.Helix || max == w.Sheet || max == w.Turn; // a class is always assignable
            });
            return (countOk && classOk).Label($"windows {windows.Count} ≠ N−w+1 = {t.seq.Length - t.window + 1}");
        });
    }

    /// <summary>
    /// Metamorphic: a homopolymer yields a constant profile — every window mean equals that residue's
    /// propensity (mean of identical values). (Chou & Fasman 1978)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SecondaryStructure_Homopolymer_IsConstantProfile()
    {
        var arb = (from aa in Gen.Elements("ARNDCEQGHILKMFPSTWYV".ToCharArray())
                   from reps in Gen.Choose(1, 12)
                   from w in Gen.Choose(1, reps)
                   select (seq: new string(aa, reps), aa, w)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var windows = SequenceStatistics.PredictSecondaryStructure(t.seq, t.w).ToList();
            var p = ChouFasmanOracle[t.aa];
            return windows.All(win =>
                Math.Abs(win.Helix - p.Helix) < CompTolerance
                && Math.Abs(win.Sheet - p.Sheet) < CompTolerance
                && Math.Abs(win.Turn - p.Turn) < CompTolerance)
                .Label($"homopolymer {t.aa} profile not constant at its propensity");
        });
    }

    /// <summary>D (determinism): the secondary-structure profile is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property SecondaryStructure_IsDeterministic()
    {
        return Prop.ForAll(SecStructProblemArbitrary(), t =>
            SequenceStatistics.PredictSecondaryStructure(t.seq, t.window)
                .SequenceEqual(SequenceStatistics.PredictSecondaryStructure(t.seq, t.window))
                .Label("PredictSecondaryStructure is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: a window over "AAAAA" (w=5) is alanine's propensity (1.42, 0.83, 0.66); a window larger than
    /// the sequence yields no profile; the empty sequence yields no profile. (Chou & Fasman 1978)
    /// </summary>
    [Test]
    [Category("Property")]
    public void SecondaryStructure_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var ala = SequenceStatistics.PredictSecondaryStructure("AAAAA", 5).ToList();
            Assert.That(ala, Has.Count.EqualTo(1));
            Assert.That(ala[0].Helix, Is.EqualTo(1.42).Within(CompTolerance));
            Assert.That(ala[0].Sheet, Is.EqualTo(0.83).Within(CompTolerance));
            Assert.That(ala[0].Turn, Is.EqualTo(0.66).Within(CompTolerance));

            Assert.That(SequenceStatistics.PredictSecondaryStructure("AA", 5), Is.Empty, "Window > length ⇒ no profile.");
            Assert.That(SequenceStatistics.PredictSecondaryStructure("", 3), Is.Empty, "Empty ⇒ no profile.");
        });
    }

    #endregion

    #region SEQ-STATS-001 / SEQ-SUMMARY-001 — Sequence Summary Statistics

    // -------------------------------------------------------------------------
    // Theory (summary statistics over a nucleotide sequence):
    //   • Composition counts ≥ 0; over {A,T,G,C,U,N} they sum to the sequence length.   (#127 R, P)
    //   • GcContent (fraction, ≡ GC%/100) ∈ [0,1].                                       (#127/#128 R GC%∈[0,100])
    //   • Reported Length equals the input length; summary fields agree with the         (#128 P length)
    //     standalone composition/GC computation.                                          (consistency)
    //
    // Verified against the standalone CalculateNucleotideComposition (a sibling already
    // pinned to an independent oracle in SEQ-COMPOSITION-001).
    // -------------------------------------------------------------------------

    private static Arbitrary<string> AtgcunArbitrary() =>
        Gen.Elements('A', 'T', 'G', 'C', 'U', 'N', 'a', 't', 'g', 'c', 'u', 'n')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// SEQ-STATS-001 — R (counts ≥ 0) + P (Σ counts = length over A/T/G/C/U/N): every composition count is
    /// non-negative, and for a sequence drawn from {A,T,G,C,U,N} the six counts sum to the sequence length;
    /// the summary GC fraction is in [0,1] (≡ GC% ∈ [0,100]).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SequenceSummary_CompositionCountsPartition_AndGcInRange()
    {
        return Prop.ForAll(AtgcunArbitrary(), seq =>
        {
            var summary = SequenceStatistics.SummarizeNucleotideSequence(seq);
            bool nonNeg = summary.Composition.Values.All(v => v >= 0);
            bool sumOk = summary.Composition.Values.Sum() == seq.Length; // all chars are A/T/G/C/U/N
            bool gcOk = summary.GcContent is >= 0.0 and <= 1.0;
            return (nonNeg && sumOk && gcOk)
                .Label($"counts sum {summary.Composition.Values.Sum()} ≠ length {seq.Length}; GC {summary.GcContent}");
        });
    }

    /// <summary>
    /// SEQ-SUMMARY-001 — P (reported length = sequence length) + R (GC% ∈ [0,100]) + consistency: the summary
    /// reports the input length, a GC fraction in [0,1], and its GcContent/Composition agree with the
    /// standalone <c>CalculateNucleotideComposition</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SequenceSummary_ReportsLength_GcInRange_AndAgreesWithComposition()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var summary = SequenceStatistics.SummarizeNucleotideSequence(seq);
            var comp = SequenceStatistics.CalculateNucleotideComposition(seq);

            bool lengthOk = summary.Length == seq.Length;
            bool gcOk = summary.GcContent is >= 0.0 and <= 1.0 && Math.Abs(summary.GcContent - comp.GcContent) < CompTolerance;
            bool compOk = summary.Composition['A'] == comp.CountA && summary.Composition['T'] == comp.CountT
                          && summary.Composition['G'] == comp.CountG && summary.Composition['C'] == comp.CountC
                          && summary.Composition['U'] == comp.CountU && summary.Composition['N'] == comp.CountN;
            return (lengthOk && gcOk && compOk)
                .Label($"summary length {summary.Length} vs {seq.Length}; GC {summary.GcContent} vs {comp.GcContent}");
        });
    }

    /// <summary>D (determinism): the sequence summary is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property SequenceSummary_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var a = SequenceStatistics.SummarizeNucleotideSequence(seq);
            var b = SequenceStatistics.SummarizeNucleotideSequence(seq);
            return (a.Length == b.Length && a.GcContent == b.GcContent && a.Entropy == b.Entropy
                    && a.Complexity == b.Complexity && a.MeltingTemperature == b.MeltingTemperature
                    && a.Composition.Count == b.Composition.Count
                    && a.Composition.All(kv => b.Composition.TryGetValue(kv.Key, out int v) && v == kv.Value))
                .Label("SummarizeNucleotideSequence is not deterministic for identical input");
        });
    }

    /// <summary>
    /// Anchors: "ATGCGC" reports length 6, GC fraction 4/6; null and empty both report length 0 and GC 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void SequenceSummary_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var s = SequenceStatistics.SummarizeNucleotideSequence("ATGCGC");
            Assert.That(s.Length, Is.EqualTo(6));
            Assert.That(s.GcContent, Is.EqualTo(4.0 / 6.0).Within(CompTolerance), "G/C = 4 of 6.");
            Assert.That(s.Composition.Values.Sum(), Is.EqualTo(6));

            var nul = SequenceStatistics.SummarizeNucleotideSequence(null);
            Assert.That(nul.Length, Is.EqualTo(0));
            Assert.That(nul.GcContent, Is.EqualTo(0.0));
            Assert.That(SequenceStatistics.SummarizeNucleotideSequence("").Length, Is.EqualTo(0));
        });
    }

    #endregion
}
