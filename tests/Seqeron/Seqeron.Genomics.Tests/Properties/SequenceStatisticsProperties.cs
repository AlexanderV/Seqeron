using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the sequence-statistics algorithm group (<see cref="SequenceStatistics"/>).
///
/// This file is the single home for the Statistics block of checklist 01 (rows #121–#130, and the
/// profile/codon rows #227, #232, #234). Each test unit lives in its own <c>#region</c>. Oracles are
/// derived INDEPENDENTLY from the cited theory/doc, never routed through the production result, so a
/// self-consistent-but-wrong production formula is still caught.
///
/// Test Units: SEQ-COMPOSITION-001, SEQ-CODON-FREQ-001, SEQ-ENTROPY-PROFILE-001, SEQ-GC-PROFILE-001
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

    private const double NTerminusPkaOracle = 7.5; // EMBOSS 6.6.0 Epk.dat "Amino" (the value iep actually uses)
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

    #region SEQ-THERMO-001 — DNA Duplex Thermodynamics (nearest-neighbor)

    // -------------------------------------------------------------------------
    // Theory (Allawi & SantaLucia 1997; SantaLucia 1998):
    //   • ΔH°, ΔS°, ΔG° from the unified nearest-neighbor model; all finite for len ≥ 2.   (R finite)
    //   • ΔG°₃₇ = ΔH° − T·ΔS°/1000 with T = 310.15 K (internal consistency).
    //   • A GC-rich duplex is more stable (lower, more-negative ΔG) than an AT-rich one
    //     of the same length.                                                              (M more GC → lower ΔG)
    //
    // The ΔG = ΔH − T·ΔS/1000 relation is reconstructed independently (within rounding).
    // -------------------------------------------------------------------------

    private const double ThermoRefTempKelvin = 310.15;

    private static Arbitrary<string> DuplexDnaArbitrary() =>
        (from n in Gen.Choose(2, 30)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         select new string(chars)).ToArbitrary();

    /// <summary>
    /// R (checklist "ΔG, ΔH, ΔS finite"): for any DNA duplex (length ≥ 2) the nearest-neighbor ΔH°, ΔS°, ΔG°
    /// and Tm are all finite, and ΔG° equals the independent ΔH° − T·ΔS°/1000 within rounding.
    /// (Allawi & SantaLucia 1997)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Thermodynamics_AreFinite_AndDeltaGIsConsistent()
    {
        return Prop.ForAll(DuplexDnaArbitrary(), seq =>
        {
            var t = SequenceStatistics.CalculateThermodynamics(seq);
            bool finite = double.IsFinite(t.DeltaH) && double.IsFinite(t.DeltaS)
                          && double.IsFinite(t.DeltaG) && double.IsFinite(t.MeltingTemperature);
            double oracleG = t.DeltaH - (ThermoRefTempKelvin * t.DeltaS / 1000.0);
            bool consistent = Math.Abs(t.DeltaG - oracleG) < 0.02; // both fields rounded to 2 dp
            return (finite && consistent).Label($"ΔH={t.DeltaH}, ΔS={t.DeltaS}, ΔG={t.DeltaG} vs {oracleG}");
        });
    }

    /// <summary>
    /// M (checklist "more GC → more stable, lower ΔG"): at any matched length a GC-rich duplex (GCGC…) has a
    /// strictly lower (more negative) ΔG° than an AT-rich duplex (ATAT…) — GC pairs stabilize the duplex.
    /// (Allawi & SantaLucia 1997)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Thermodynamics_GcRich_IsMoreStableThanAtRich()
    {
        return Prop.ForAll(Gen.Choose(2, 30).ToArbitrary(), n =>
        {
            string gc = string.Concat(Enumerable.Range(0, n).Select(i => i % 2 == 0 ? 'G' : 'C'));
            string at = string.Concat(Enumerable.Range(0, n).Select(i => i % 2 == 0 ? 'A' : 'T'));
            double gcG = SequenceStatistics.CalculateThermodynamics(gc).DeltaG;
            double atG = SequenceStatistics.CalculateThermodynamics(at).DeltaG;
            return (gcG < atG).Label($"ΔG(GC)={gcG} not < ΔG(AT)={atG} at length {n}");
        });
    }

    /// <summary>D (determinism): duplex thermodynamics are identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property Thermodynamics_IsDeterministic()
    {
        return Prop.ForAll(DuplexDnaArbitrary(), seq =>
            SequenceStatistics.CalculateThermodynamics(seq).Equals(SequenceStatistics.CalculateThermodynamics(seq))
                .Label("CalculateThermodynamics is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: a GC duplex has negative (stable) ΔG; a sequence shorter than 2 yields all-zero properties.
    /// (Allawi & SantaLucia 1997)
    /// </summary>
    [Test]
    [Category("Property")]
    public void Thermodynamics_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var gc = SequenceStatistics.CalculateThermodynamics("GCGCGCGC");
            Assert.That(gc.DeltaG, Is.LessThan(0.0), "A GC-rich duplex is thermodynamically stable (ΔG < 0).");
            Assert.That(double.IsFinite(gc.DeltaH) && double.IsFinite(gc.DeltaS), Is.True);

            var tooShort = SequenceStatistics.CalculateThermodynamics("A");
            Assert.That(tooShort.DeltaH, Is.EqualTo(0.0));
            Assert.That(tooShort.DeltaG, Is.EqualTo(0.0));
        });
    }

    #endregion

    #region SEQ-TM-001 — Melting Temperature (Wallace / Marmur-Doty)

    // -------------------------------------------------------------------------
    // Theory (Wallace rule; Marmur-Doty):
    //   • Short oligos (length < 14): Tm = 2·(A+T) + 4·(G+C).
    //   • Longer: Marmur-Doty Tm increasing in GC fraction. Both give Tm ≥ 0.            (R Tm ≥ 0)
    //   • More GC ⇒ higher Tm (G/C contribute 4 vs 2 in Wallace; positive GC coefficient
    //     in Marmur-Doty).                                                                (M more GC → higher Tm)
    //
    // The Wallace rule is reconstructed independently for the short-oligo regime.
    // -------------------------------------------------------------------------

    private static Arbitrary<string> DnaAnyLengthArbitrary() =>
        (from n in Gen.Choose(1, 40)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         select new string(chars)).ToArbitrary();

    /// <summary>
    /// R (checklist "Tm ≥ 0"): the melting temperature is non-negative for any DNA sequence, and for short
    /// oligos (length &lt; 14) it equals the independent Wallace rule 2·(A+T) + 4·(G+C). (Wallace rule)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_NonNegative_AndWallaceForShortOligos()
    {
        return Prop.ForAll(DnaAnyLengthArbitrary(), seq =>
        {
            double tm = SequenceStatistics.CalculateMeltingTemperature(seq);
            bool nonNeg = tm >= 0.0;

            bool wallaceOk = true;
            if (seq.Length < 14)
            {
                int at = seq.Count(c => c is 'A' or 'T');
                int gc = seq.Count(c => c is 'G' or 'C');
                wallaceOk = Math.Abs(tm - (2 * at + 4 * gc)) < CompTolerance;
            }

            return (nonNeg && wallaceOk).Label($"Tm {tm} (len {seq.Length})");
        });
    }

    /// <summary>
    /// M (checklist "more GC → higher Tm"): at any matched length a GC-only sequence has a strictly higher
    /// melting temperature than an AT-only sequence. (Wallace 4 vs 2; Marmur-Doty positive GC coefficient)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_GcRich_IsHigherThanAtRich()
    {
        return Prop.ForAll(Gen.Choose(1, 40).ToArbitrary(), n =>
        {
            double gcTm = SequenceStatistics.CalculateMeltingTemperature(new string('G', n));
            double atTm = SequenceStatistics.CalculateMeltingTemperature(new string('A', n));
            return (gcTm > atTm).Label($"Tm(G×{n})={gcTm} not > Tm(A×{n})={atTm}");
        });
    }

    /// <summary>D (determinism): melting temperature is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IsDeterministic()
    {
        return Prop.ForAll(DnaAnyLengthArbitrary(), seq =>
            (SequenceStatistics.CalculateMeltingTemperature(seq) == SequenceStatistics.CalculateMeltingTemperature(seq))
                .Label("CalculateMeltingTemperature is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: Wallace rule "AATT" ⇒ 8 °C, "GGCC" ⇒ 16 °C; empty ⇒ 0. (Wallace rule)
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperature_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateMeltingTemperature("AATT"), Is.EqualTo(8.0).Within(CompTolerance), "2·4 = 8.");
            Assert.That(SequenceStatistics.CalculateMeltingTemperature("GGCC"), Is.EqualTo(16.0).Within(CompTolerance), "4·4 = 16.");
            Assert.That(SequenceStatistics.CalculateMeltingTemperature(""), Is.EqualTo(0.0), "Empty ⇒ 0.");
        });
    }

    #endregion

    #region Profile generators

    private static Arbitrary<string> ProfileDnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    /// <summary>DNA of length a multiple of 3 (≥ 3), all ACGT — full codons in frame 0.</summary>
    private static Arbitrary<string> CodonDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 3)
            .Select(a => new string(a, 0, a.Length - a.Length % 3)).Where(s => s.Length >= 3)
            .ToArbitrary();

    #endregion

    #region SEQ-CODON-FREQ-001: R: each freq ≥ 0; P: Σ frequencies = 1 (Σ counts = len/3); D: deterministic

    // CalculateCodonFrequencies returns count/total per codon, so the frequencies sum to 1 over the
    // len/3 in-frame codons (Kazusa CUTG; ambiguous codons excluded — none here as input is ACGT).

    /// <summary>
    /// INV-1 (R + P): every codon frequency is in [0,1] and the frequencies sum to 1 (i.e. the codon
    /// counts sum to the number of in-frame codons, len/3).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonFrequencies_NonNegative_SumToOne()
    {
        return Prop.ForAll(CodonDnaArbitrary(), seq =>
        {
            var freq = SequenceStatistics.CalculateCodonFrequencies(seq);
            bool nonNeg = freq.Values.All(v => v is >= 0.0 and <= 1.0);
            bool sumsToOne = Math.Abs(freq.Values.Sum() - 1.0) < 1e-9;
            return (nonNeg && sumsToOne).Label($"freqs non-negative={nonNeg}, sum={freq.Values.Sum()}");
        });
    }

    /// <summary>INV-2 (D): codon frequencies are deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property CodonFrequencies_AreDeterministic()
    {
        return Prop.ForAll(CodonDnaArbitrary(), seq =>
        {
            var a = SequenceStatistics.CalculateCodonFrequencies(seq);
            var b = SequenceStatistics.CalculateCodonFrequencies(seq);
            return (a.Count == b.Count && a.All(kv => Math.Abs(kv.Value - b[kv.Key]) < 1e-12))
                .Label("CalculateCodonFrequencies must be deterministic");
        });
    }

    /// <summary>
    /// INV-3 (review 2026-09 B03 F15): RNA U is read as T (EMBOSS cusp, CodonW), so the lower-case RNA
    /// spelling of a CDS gives exactly the DNA table in every frame.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonFrequencies_RnaSpelling_EqualsDnaSpelling()
    {
        return Prop.ForAll(CodonDnaArbitrary(), seq =>
        {
            string rna = seq.Replace('T', 'U').ToLowerInvariant();
            bool ok = Enumerable.Range(0, 3).All(f =>
            {
                var dna = SequenceStatistics.CalculateCodonFrequencies(seq, f);
                var r = SequenceStatistics.CalculateCodonFrequencies(rna, f);
                return dna.Count == r.Count && dna.All(kv => r.TryGetValue(kv.Key, out double v) && v == kv.Value);
            });
            return ok.Label("RNA spelling must equal DNA spelling");
        });
    }

    #endregion

    #region SEQ-ENTROPY-PROFILE-001: R: each entropy ≥ 0; P: profile length = len−w+1; D: deterministic

    private const int ProfileWindow = 10;

    /// <summary>
    /// INV-1 (R + P): the entropy profile has len−w+1 entries (step 1), each a non-negative Shannon
    /// entropy.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EntropyProfile_LengthAndNonNegative()
    {
        return Prop.ForAll(ProfileDnaArbitrary(20), seq =>
        {
            var profile = SequenceStatistics.CalculateEntropyProfile(seq, ProfileWindow, 1).ToList();
            return (profile.Count == seq.Length - ProfileWindow + 1 && profile.All(e => e >= -1e-12))
                .Label($"entropy profile length {profile.Count}, expected {seq.Length - ProfileWindow + 1}");
        });
    }

    /// <summary>INV-2 (D): the entropy profile is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property EntropyProfile_IsDeterministic()
    {
        return Prop.ForAll(ProfileDnaArbitrary(20), seq =>
            SequenceStatistics.CalculateEntropyProfile(seq, ProfileWindow, 1)
                .SequenceEqual(SequenceStatistics.CalculateEntropyProfile(seq, ProfileWindow, 1))
                .Label("CalculateEntropyProfile must be deterministic"));
    }

    #endregion

    #region SEQ-GC-PROFILE-001: R: each GC% ∈ [0,100]; P: profile length = len−w+1; D: deterministic

    /// <summary>
    /// INV-1 (R + P): the GC profile has len−w+1 entries (step 1), each a GC percentage in [0,100].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcProfile_LengthAndRange()
    {
        return Prop.ForAll(ProfileDnaArbitrary(20), seq =>
        {
            var profile = SequenceStatistics.CalculateGcContentProfile(seq, ProfileWindow, 1).ToList();
            return (profile.Count == seq.Length - ProfileWindow + 1 && profile.All(g => g is >= 0.0 and <= 100.0))
                .Label($"GC profile length {profile.Count}, expected {seq.Length - ProfileWindow + 1}");
        });
    }

    /// <summary>INV-2 (D): the GC profile is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property GcProfile_IsDeterministic()
    {
        return Prop.ForAll(ProfileDnaArbitrary(20), seq =>
            SequenceStatistics.CalculateGcContentProfile(seq, ProfileWindow, 1)
                .SequenceEqual(SequenceStatistics.CalculateGcContentProfile(seq, ProfileWindow, 1))
                .Label("CalculateGcContentProfile must be deterministic"));
    }

    #endregion

    #region Review 2026-09 B03 — new behaviour (heavy tier)

    // -------------------------------------------------------------------------
    // Invariants of the behaviour introduced by review batch B03 (docs/Validation/review-2026-09/B03.md).
    // Each follows from the model definition, not from observed output:
    //   F3–F6  pI / net charge: every Henderson–Hasselbalch term is strictly decreasing in pH, so the net
    //          charge is monotone non-increasing on both pK scales; the pI is the unique root, located by
    //          bisection to 1e-9 and rounded to 0.01, so the charge changes sign within ±0.005 of it;
    //          |charge| ≤ ionizable groups (each group carries at most one charge).
    //   F8     edge-weighted hydropathy: a normalised weighted mean (Σw·kd/Σw) — homopolymer windows
    //          return the residue's kd, every value lies in [min kd, max kd] = [−4.5, 4.5], count n−W+1.
    //   F2     ds / circular nucleotide MW (Biopython molecular_weight): ds = ss(seq) + ss(complement);
    //          circular = linear − one water (18.0153 Da) per strand.
    //   F9/F10 NN thermodynamics: Biopython _check (U→T, non-ACGT dropped); the symmetry correction adds
    //          ΔS −1.4 and leaves ΔH unchanged.
    //   F12/F20 basic Tm: switch on the A/C/G/T/U count, U read as T, other symbols ignored.
    //   F13    Karlin ρ*: computed on seq + revcomp(seq), hence revcomp-invariant and ρ*_XY = ρ*_revcomp(XY).
    //   F19    GC profile ambiguity modes: Remove = default when the window holds no S/W; values in range.
    //   F21    linguistic complexity = canonical SequenceComplexity (Σ V_k / Σ V_max,k) ∈ [0, 1].
    // -------------------------------------------------------------------------

    private static readonly SequenceStatistics.PkaScale[] AllPkaScales =
        { SequenceStatistics.PkaScale.Emboss, SequenceStatistics.PkaScale.Bjellqvist };

    private static string DnaComplementOracle(string s) =>
        new(s.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static string RnaComplementOracle(string s) =>
        new(s.Select(c => c switch { 'A' => 'U', 'U' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static string RevCompOracle(string s) =>
        new(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static Arbitrary<string> AlphabetArbitrary(string alphabet, int minLen, int maxLen) =>
        (from n in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements(alphabet.ToCharArray()).ArrayOf(n)
         select new string(chars)).ToArbitrary();

    /// <summary>
    /// F3–F6 (P, MON): the net charge is monotone non-increasing in pH on both pK scales, and bounded by the
    /// number of ionizable groups (termini + D/E/C/Y/H/K/R).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NetCharge_MonotoneNonIncreasingInPh_AndBoundedByGroups()
    {
        return Prop.ForAll(BoundedProteinArbitrary(40), seq =>
        {
            if (seq.Length == 0) return true.ToProperty();
            int acidic = seq.Count(c => c is 'D' or 'E' or 'C' or 'Y') + 1;
            int basic = seq.Count(c => c is 'H' or 'K' or 'R') + 1;
            foreach (var scale in AllPkaScales)
            {
                double previous = double.PositiveInfinity;
                for (int k = 0; k <= 56; k++)
                {
                    double pH = k * 0.25;
                    double q = SequenceStatistics.CalculateNetCharge(seq, pH, scale);
                    if (q > previous + 1e-12)
                        return false.Label($"{scale}: charge rose from {previous} to {q} at pH {pH} for '{seq}'");
                    if (q < -acidic - 1e-12 || q > basic + 1e-12)
                        return false.Label($"{scale}: charge {q} outside [−{acidic}, {basic}] at pH {pH}");
                    previous = q;
                }
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// F5/F6 (P): on both scales pI ∈ [0,14] and the exact net charge changes sign within the rounding
    /// half-width of the returned pI (root bisected to 1e-9, rounded to 0.01) — much tighter than the
    /// oracle-based tolerance above. The 1-argument overload and the default net-charge scale are EMBOSS.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPoint_BothScales_ChargeChangesSignAtPi()
    {
        return Prop.ForAll(BoundedProteinArbitrary(40), seq =>
        {
            if (seq.Length == 0) return true.ToProperty();
            if (SequenceStatistics.CalculateIsoelectricPoint(seq)
                != SequenceStatistics.CalculateIsoelectricPoint(seq, SequenceStatistics.PkaScale.Emboss))
                return false.Label("1-argument pI must equal the EMBOSS scale");
            if (SequenceStatistics.CalculateNetCharge(seq, 7.0)
                != SequenceStatistics.CalculateNetCharge(seq, 7.0, SequenceStatistics.PkaScale.Emboss))
                return false.Label("default net-charge scale must be EMBOSS");

            foreach (var scale in AllPkaScales)
            {
                double pi = SequenceStatistics.CalculateIsoelectricPoint(seq, scale);
                if (pi is < 0.0 or > 14.0)
                    return false.Label($"{scale}: pI {pi} outside [0,14]");
                double below = SequenceStatistics.CalculateNetCharge(seq, pi - 0.0051, scale);
                double above = SequenceStatistics.CalculateNetCharge(seq, pi + 0.0051, scale);
                if (below < 0 || above > 0)
                    return false.Label($"{scale}: charge {below} at pI−0.0051, {above} at pI+0.0051 (pI {pi}, '{seq}')");
            }
            return true.ToProperty();
        });
    }

    /// <summary>F3/F6 (INV): pI and net charge are case-insensitive on both scales.</summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPoint_BothScales_CaseInsensitive()
    {
        return Prop.ForAll(BoundedProteinArbitrary(30), seq =>
            AllPkaScales.All(scale =>
                SequenceStatistics.CalculateIsoelectricPoint(seq.ToLowerInvariant(), scale)
                    == SequenceStatistics.CalculateIsoelectricPoint(seq, scale)
                && SequenceStatistics.CalculateNetCharge(seq.ToLowerInvariant(), 6.3, scale)
                    == SequenceStatistics.CalculateNetCharge(seq, 6.3, scale))
            .Label("lower-case input must give the same pI and charge"));
    }

    /// <summary>
    /// F8 (R + P): for any odd window and edge weight in [0,1] the profile has n−W+1 values in
    /// [−4.5, 4.5] (a normalised weighted mean of KD values, unknown residues contributing 0);
    /// edge = 1 is identical to the unweighted 2-argument call.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HydrophobicityProfile_EdgeWeighted_LengthRangeAndUnitEdgeIdentity()
    {
        var gen = from seq in ProteinArbitrary().Generator
                  from half in Gen.Choose(0, 6)
                  from edgeTenths in Gen.Choose(0, 10)
                  select (seq, w: 2 * half + 1, edge: edgeTenths / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var profile = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge).ToList();
            int expected = Math.Max(0, t.seq.Length - t.w + 1);
            bool lengthOk = profile.Count == expected;
            bool rangeOk = profile.All(v => v is >= -4.5 - 1e-12 and <= 4.5 + 1e-12);
            bool unitOk = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, 1.0)
                .SequenceEqual(SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w));
            return (lengthOk && rangeOk && unitOk)
                .Label($"W={t.w} edge={t.edge}: count {profile.Count}/{expected}, range {rangeOk}, edge-1 identity {unitOk}");
        });
    }

    /// <summary>
    /// F8 (P): the weights are normalised (Σw·kd / Σw), so a homopolymer window returns the residue's
    /// Kyte–Doolittle value for every odd window and every edge weight.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HydrophobicityProfile_EdgeWeighted_HomopolymerIsResidueValue()
    {
        var gen = from aa in Gen.Elements("ARNDCEQGHILKMFPSTWYV".ToCharArray())
                  from half in Gen.Choose(0, 6)
                  from extra in Gen.Choose(0, 10)
                  from edgeTenths in Gen.Choose(0, 10)
                  select (aa, w: 2 * half + 1, extra, edge: edgeTenths / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string seq = new(t.aa, t.w + t.extra);
            double kd = OracleGravy(t.aa.ToString());
            var profile = SequenceStatistics.CalculateHydrophobicityProfile(seq, t.w, t.edge).ToList();
            return (profile.Count == t.extra + 1 && profile.All(v => Math.Abs(v - kd) < 1e-12))
                .Label($"poly-{t.aa} W={t.w} edge={t.edge}: expected {kd}, got [{string.Join(", ", profile)}]");
        });
    }

    /// <summary>
    /// F2 (P): double-stranded MW = single-strand MW of the strand + single-strand MW of its Watson–Crick
    /// complement (DNA and RNA, linear and circular); circular = linear − one water per strand.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideMolecularWeight_DoubleStrandedAndCircular_Decompose()
    {
        var gen = from isDna in Gen.Elements(true, false)
                  from seq in AlphabetArbitrary(isDna ? "ACGT" : "ACGU", 1, 60).Generator
                  select (isDna, seq);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string comp = t.isDna ? DnaComplementOracle(t.seq) : RnaComplementOracle(t.seq);
            bool ok = true;
            var sb = new System.Text.StringBuilder();
            foreach (bool circular in new[] { false, true })
            {
                double ds = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, true, circular);
                double sum = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, false, circular)
                             + SequenceStatistics.CalculateNucleotideMolecularWeight(comp, t.isDna, false, circular);
                if (Math.Abs(ds - sum) > 1e-7) { ok = false; sb.Append($" ds {ds} ≠ ss+comp {sum} (circular {circular});"); }
            }
            foreach (bool ds in new[] { false, true })
            {
                double linear = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, ds, false);
                double circ = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, ds, true);
                double expected = linear - (ds ? 2 : 1) * WaterMassOracle;
                if (Math.Abs(circ - expected) > 1e-7) { ok = false; sb.Append($" circular {circ} ≠ {expected} (ds {ds});"); }
            }
            double legacy = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna);
            if (legacy != SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, false, false))
            { ok = false; sb.Append(" 2-arg overload ≠ ss linear;"); }
            return ok.Label($"{(t.isDna ? "DNA" : "RNA")} '{t.seq}':{sb}");
        });
    }

    /// <summary>
    /// F9/F10 (P + INV): the self-complementary model only adds the symmetry entropy (ΔH unchanged,
    /// ΔS − 1.4 within the 2-dp output rounding); selfComplementary=false equals the 3-argument call;
    /// the RNA/lower-case spelling and inserted N/gaps/whitespace are normalised away (Biopython _check).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Thermodynamics_SymmetryCorrection_AndBiopythonCheckNormalisation()
    {
        return Prop.ForAll(AlphabetArbitrary("ACGT", 2, 40), seq =>
        {
            var plain = SequenceStatistics.CalculateThermodynamics(seq, 0.05, 2.5e-7);
            var nonSelf = SequenceStatistics.CalculateThermodynamics(seq, 0.05, 2.5e-7, selfComplementary: false);
            var self = SequenceStatistics.CalculateThermodynamics(seq, 0.05, 2.5e-7, selfComplementary: true);
            string junked = string.Concat(seq.Select((c, i) => i % 3 == 0 ? "N" + c + " -" : c.ToString()));
            var rna = SequenceStatistics.CalculateThermodynamics(seq.Replace('T', 'U').ToLowerInvariant(), 0.05, 2.5e-7);
            var junk = SequenceStatistics.CalculateThermodynamics(junked, 0.05, 2.5e-7);

            bool overloadOk = nonSelf == plain;
            bool dhOk = self.DeltaH == plain.DeltaH;
            bool dsOk = Math.Abs(self.DeltaS - (plain.DeltaS - 1.4)) <= 0.0100001;
            bool normOk = rna == plain && junk == plain;
            return (overloadOk && dhOk && dsOk && normOk)
                .Label($"'{seq}': overload {overloadOk}, ΔH {self.DeltaH}/{plain.DeltaH}, ΔS {self.DeltaS}/{plain.DeltaS}, normalisation {normOk}");
        });
    }

    /// <summary>
    /// F12/F20 (INV): basic Tm reads U as T and ignores every non-ACGTU symbol, including for the
    /// Wallace/GC switch (it depends on the base count, not the raw length).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_RnaSpellingAndJunk_Invariant()
    {
        return Prop.ForAll(AlphabetArbitrary("ACGT", 1, 30), seq =>
        {
            string junked = string.Concat(seq.Select(c => "N-" + c + " R"));
            bool ok = new[] { true, false }.All(wallace =>
            {
                double tm = SequenceStatistics.CalculateMeltingTemperature(seq, wallace);
                return SequenceStatistics.CalculateMeltingTemperature(seq.Replace('T', 'U').ToLowerInvariant(), wallace) == tm
                       && SequenceStatistics.CalculateMeltingTemperature(junked, wallace) == tm;
            });
            return ok.Label($"'{seq}': RNA spelling or junk insertion changed Tm");
        });
    }

    /// <summary>
    /// F13 (INV + P): Karlin ρ* is computed on the sequence together with its reverse complement, so it is
    /// invariant under reverse complementation and satisfies ρ*_XY = ρ*_revcomp(XY); strandSymmetric=false
    /// equals the 1-argument single-strand ρ.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DinucleotideRatios_StrandSymmetric_RevcompInvariantAndSymmetric()
    {
        return Prop.ForAll(AlphabetArbitrary("ACGTN", 2, 60), seq =>
        {
            var rho = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: true);
            var rhoRc = SequenceStatistics.CalculateDinucleotideRatios(RevCompOracle(seq), strandSymmetric: true);
            bool rcOk = rho.Count == rhoRc.Count
                        && rho.All(kv => rhoRc.TryGetValue(kv.Key, out double v) && Math.Abs(v - kv.Value) <= 1e-12 * Math.Max(1, v));
            bool symOk = rho.All(kv => rho.TryGetValue(RevCompOracle(kv.Key), out double v) && Math.Abs(v - kv.Value) <= 1e-12 * Math.Max(1, v));
            var single = SequenceStatistics.CalculateDinucleotideRatios(seq);
            var single2 = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: false);
            bool overloadOk = single.Count == single2.Count && single.All(kv => single2[kv.Key] == kv.Value);
            return (rcOk && symOk && overloadOk)
                .Label($"'{seq}': revcomp-invariant {rcOk}, ρ*_XY=ρ*_rc(XY) {symOk}, overload {overloadOk}");
        });
    }

    /// <summary>
    /// F19 (P + R): with no S/W in the input, the Biopython Remove mode equals the default profile (both
    /// count A/C/G/T/U case-insensitively and exclude N/gaps); every mode yields values in [0,100]
    /// (percent) and the fraction form is the percent form / 100.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcProfile_AmbiguityModes_RemoveEqualsDefaultAndInRange()
    {
        var gen = from seq in AlphabetArbitrary("ACGTUNacgtun-RY", 1, 60).Generator
                  from w in Gen.Choose(1, 12)
                  from step in Gen.Choose(1, 4)
                  select (seq, w, step);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var def = SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step).ToList();
            var remove = SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, false,
                SequenceExtensions.GcAmbiguityMode.Remove).ToList();
            bool removeOk = def.Count == remove.Count && def.Zip(remove).All(p => Math.Abs(p.First - p.Second) <= 1e-9);
            bool modesOk = Enum.GetValues<SequenceExtensions.GcAmbiguityMode>().All(mode =>
            {
                var pct = SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, false, mode).ToList();
                var frac = SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, true, mode).ToList();
                return pct.Count == def.Count && pct.All(v => v is >= 0.0 and <= 100.0)
                       && pct.Zip(frac).All(p => Math.Abs(p.First / 100.0 - p.Second) <= 1e-12);
            });
            return (removeOk && modesOk).Label($"'{t.seq}' W={t.w} s={t.step}: Remove=default {removeOk}, modes {modesOk}");
        });
    }

    /// <summary>
    /// F21 (P + R): linguistic complexity is the canonical Orlov–Potapov Σ V_k / Σ V_max,k
    /// (SequenceComplexity) for every maxK, and lies in [0,1] over the DNA alphabet.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LinguisticComplexity_EqualsCanonical_InUnitInterval()
    {
        var gen = from seq in AlphabetArbitrary("ACGT", 1, 50).Generator
                  from m in Gen.Choose(1, 10)
                  select (seq, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double lc = SequenceStatistics.CalculateLinguisticComplexity(t.seq, t.m);
            double canonical = SequenceComplexity.CalculateLinguisticComplexity(t.seq, t.m);
            return (lc == canonical && lc is >= 0.0 and <= 1.0)
                .Label($"'{t.seq}' m={t.m}: {lc} vs canonical {canonical}");
        });
    }

    #endregion
    #region FIN-B03 heavy tier — finisher options F27–F51 (MW, pI, hydropathy, NN, profiles, Chou–Fasman, codon, ρ, summary, Tm)

    // -------------------------------------------------------------------------
    // Invariants of the opt-in behaviour added by the B03 finisher (review 2026-09 B03.md F27–F51).
    // Oracles are the definitions quoted in the API docs (Biopython 1.88 molecular_weight / ProtParam /
    // MeltingTemp, EMBOSS 6.6.0 iep / cusp, Karlin ρ, Chou & Fasman 1978), recomputed here.
    // -------------------------------------------------------------------------

    private const string KnownProteinMassAlphabet = "ARNDCEQGHILKMFPSTWYVUO"; // Biopython protein_weights keys
    private const double MonoisotopicWaterOracle = 18.010565;                 // Biopython molecular_weight(monoisotopic=True)

    private static char UpperAscii(char c) => char.IsAsciiLetterLower(c) ? (char)(c - 32) : c;

    private static readonly IReadOnlyDictionary<char, double>[] ShippedScales =
    {
        ProteinScales.KyteDoolittle, ProteinScales.HoppWoods, ProteinScales.Emini, ProteinScales.Janin, ProteinScales.Flexibility,
    };

    /// <summary>
    /// F28 (P + INV): protein MW options. Default / null options equal the 1-argument overload; a circular
    /// (head-to-tail) peptide is the linear one minus exactly one water (average 18.0153 or monoisotopic
    /// 18.010565 Da) whenever a residue has a mass; strict mode throws iff a non-whitespace symbol has no mass
    /// and otherwise equals the lenient value.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeightOptions_Protein_DefaultsCircularStrict()
    {
        return Prop.ForAll(AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVUOacdu BZX*-", 0, 50), seq =>
        {
            double legacy = SequenceStatistics.CalculateMolecularWeight(seq);
            bool defaultsOk = SequenceStatistics.CalculateMolecularWeight(seq, null) == legacy
                              && SequenceStatistics.CalculateMolecularWeight(seq, SequenceStatistics.MolecularWeightOptions.Default) == legacy;

            int residues = seq.Count(c => KnownProteinMassAlphabet.Contains(UpperAscii(c)));
            bool hasUnknown = seq.Any(c => !char.IsWhiteSpace(c) && !KnownProteinMassAlphabet.Contains(UpperAscii(c)));
            var sb = new System.Text.StringBuilder();
            bool ok = defaultsOk;
            foreach (bool mono in new[] { false, true })
            {
                double water = mono ? MonoisotopicWaterOracle : WaterMassOracle;
                double linear = SequenceStatistics.CalculateMolecularWeight(seq, new() { Monoisotopic = mono });
                double circular = SequenceStatistics.CalculateMolecularWeight(seq, new() { Monoisotopic = mono, Circular = true });
                double expected = residues == 0 ? 0 : linear - water;
                if (Math.Abs(circular - expected) > 1e-9) { ok = false; sb.Append($" mono={mono}: circular {circular} ≠ {expected};"); }

                double strict = double.NaN;
                bool threw = false;
                try { strict = SequenceStatistics.CalculateMolecularWeight(seq, new() { Monoisotopic = mono, Strict = true }); }
                catch (ArgumentException) { threw = true; }
                if (threw != hasUnknown) { ok = false; sb.Append($" mono={mono}: strict threw={threw}, unknown={hasUnknown};"); }
                if (!threw && Math.Abs(strict - linear) > 1e-9) { ok = false; sb.Append($" mono={mono}: strict {strict} ≠ lenient {linear};"); }
            }
            return ok.Label($"'{seq}': defaults {defaultsOk};{sb}");
        });
    }

    /// <summary>
    /// F28 (P): for any non-empty peptide over the 20 standard residues and pyrrolysine the monoisotopic MW is
    /// below the average MW (each monoisotopic residue mass is the lightest-isotope mass, below the abundance-
    /// weighted average; the 0.0047 Da per bond water difference is far smaller). Selenocysteine is excluded on
    /// purpose: its monoisotopic mass uses ⁸⁰Se (79.9165), heavier than Se's average 78.971, so Biopython
    /// <c>monoisotopic_protein_weights</c> U = 168.964203 &gt; <c>protein_weights</c> U = 168.0532.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeightOptions_Protein_MonoisotopicBelowAverage()
    {
        return Prop.ForAll(AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVO", 1, 60), seq =>
        {
            double avg = SequenceStatistics.CalculateMolecularWeight(seq, null);
            double mono = SequenceStatistics.CalculateMolecularWeight(seq, new() { Monoisotopic = true });
            double avgCirc = SequenceStatistics.CalculateMolecularWeight(seq, new() { Circular = true });
            double monoCirc = SequenceStatistics.CalculateMolecularWeight(seq, new() { Monoisotopic = true, Circular = true });
            return (mono < avg && monoCirc < avgCirc && mono > 0).Label($"'{seq}': mono {mono} vs avg {avg}; circular {monoCirc} vs {avgCirc}");
        });
    }

    /// <summary>
    /// F28 (P): nucleotide MW options. Circular = linear − one water per strand (2 for double-stranded) with the
    /// mass set's water; monoisotopic &lt; average; default options = the 3-argument overload; strict throws iff a
    /// non-whitespace symbol is outside the DNA (ACGT) / RNA (ACGU) table.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MolecularWeightOptions_Nucleotide_CircularMonoisotopicStrict()
    {
        var gen = from isDna in Gen.Elements(true, false)
                  from seq in AlphabetArbitrary("ACGTUacgtuN -", 0, 50).Generator
                  from ds in Gen.Elements(true, false)
                  select (isDna, seq, ds);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string table = t.isDna ? "ACGT" : "ACGU";
            int monomers = t.seq.Count(c => table.Contains(UpperAscii(c)));
            bool hasUnknown = t.seq.Any(c => !char.IsWhiteSpace(c) && !table.Contains(UpperAscii(c)));
            bool ok = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, (SequenceStatistics.MolecularWeightOptions?)null)
                      == SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds);
            var sb = new System.Text.StringBuilder(ok ? "" : " default ≠ legacy;");
            foreach (bool mono in new[] { false, true })
            {
                double water = mono ? MonoisotopicWaterOracle : WaterMassOracle;
                double linear = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, new SequenceStatistics.MolecularWeightOptions { Monoisotopic = mono });
                double circ = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, new SequenceStatistics.MolecularWeightOptions { Monoisotopic = mono, Circular = true });
                double expected = monomers == 0 ? 0 : linear - (t.ds ? 2 : 1) * water;
                if (Math.Abs(circ - expected) > 1e-8) { ok = false; sb.Append($" mono={mono}: circular {circ} ≠ {expected};"); }
                bool threw = false;
                try { SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, new SequenceStatistics.MolecularWeightOptions { Monoisotopic = mono, Strict = true }); }
                catch (ArgumentException) { threw = true; }
                if (threw != hasUnknown) { ok = false; sb.Append($" mono={mono}: strict threw={threw}, unknown={hasUnknown};"); }
            }
            if (monomers > 0)
            {
                double avg = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, SequenceStatistics.MolecularWeightOptions.Default);
                double mono = SequenceStatistics.CalculateNucleotideMolecularWeight(t.seq, t.isDna, t.ds, new SequenceStatistics.MolecularWeightOptions { Monoisotopic = true });
                if (!(mono < avg)) { ok = false; sb.Append($" mono {mono} ≥ avg {avg};"); }
            }
            return ok.Label($"{(t.isDna ? "DNA" : "RNA")} ds={t.ds} '{t.seq}':{sb}");
        });
    }

    private static Gen<SequenceStatistics.IsoelectricPointOptions> PiOptionsGen(SequenceStatistics.IsoelectricPointWindow window) =>
        from scale in Gen.Elements(AllPkaScales)
        from termini in Gen.Elements(true, false)
        from lys in Gen.Choose(0, 3)
        from dis in Gen.Choose(0, 2)
        select new SequenceStatistics.IsoelectricPointOptions
        {
            Scale = scale, Window = window, Termini = termini, LysinesModified = lys, Disulphides = dis,
        };

    /// <summary>
    /// F43/F44 (P): EMBOSS <c>iep</c> window (<c>embIepIepC</c> / <c>embIepPhConverge</c>). For any options the Try
    /// method returns false iff the net charge has the same strict sign at pH 1 and pH 14 ("Isoelectric Point =
    /// None"; also for empty input); otherwise the pI lies in [1, 14], the charge there is ≥ 0 (lower bracket
    /// bound or exact zero) and ≤ 0 one bracket width (1e-4) above it.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPointOptions_EmbossWindow_NoneIffSameSignAtWindowEnds()
    {
        var gen = from seq in BoundedProteinArbitrary(30).Generator
                  from opts in PiOptionsGen(SequenceStatistics.IsoelectricPointWindow.Emboss)
                  select (seq, opts);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            bool found = SequenceStatistics.TryCalculateIsoelectricPoint(t.seq, t.opts, out double pi);
            if (t.seq.Length == 0)
                return (!found && double.IsNaN(pi)).Label("empty input must report no pI");
            double q1 = SequenceStatistics.CalculateNetCharge(t.seq, 1.0, t.opts);
            double q14 = SequenceStatistics.CalculateNetCharge(t.seq, 14.0, t.opts);
            bool sameSign = (q1 > 0 && q14 > 0) || (q1 < 0 && q14 < 0);
            if (found == sameSign)
                return false.Label($"'{t.seq}' {t.opts}: found={found} but q(1)={q1}, q(14)={q14}");
            if (!found)
                return double.IsNaN(pi).Label("no pI ⇒ NaN");
            double atPi = SequenceStatistics.CalculateNetCharge(t.seq, pi, t.opts);
            double above = SequenceStatistics.CalculateNetCharge(t.seq, pi + 1.0001e-4, t.opts);
            return (pi is >= 1.0 and <= 14.0 && atPi >= -1e-12 && above <= 1e-12)
                .Label($"'{t.seq}' {t.opts}: pI {pi}, q(pI)={atPi}, q(pI+1e-4)={above}");
        });
    }

    /// <summary>
    /// F44 (MON): EMBOSS <c>-lysinemodified N</c> removes positive groups, so charge and pI are non-increasing in
    /// N and saturate at N = #K; <c>-disulphides N</c> removes acidic cysteines, so they are non-decreasing in N and
    /// saturate at 2N ≥ #C (<c>embIepCompC</c>). Library window, both scales, termini on or off.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPointOptions_Modifications_MonotoneAndSaturating()
    {
        var gen = from seq in AlphabetArbitrary("KKCCDEHRYAGLS", 1, 30).Generator
                  from opts in PiOptionsGen(SequenceStatistics.IsoelectricPointWindow.Default)
                  select (seq, baseOpts: opts with { LysinesModified = 0, Disulphides = 0 });
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double Pi(SequenceStatistics.IsoelectricPointOptions o)
            {
                SequenceStatistics.TryCalculateIsoelectricPoint(t.seq, o, out double v);
                return v;
            }
            int nK = t.seq.Count(c => c == 'K'), nC = t.seq.Count(c => c == 'C');
            var sb = new System.Text.StringBuilder();
            for (int n = 0; n <= nK + 1; n++)
            {
                var a = t.baseOpts with { LysinesModified = n };
                var b = t.baseOpts with { LysinesModified = n + 1 };
                if (Pi(b) > Pi(a)) sb.Append($" pI rose with lysinesModified {n}→{n + 1};");
                if (SequenceStatistics.CalculateNetCharge(t.seq, 7.0, b) > SequenceStatistics.CalculateNetCharge(t.seq, 7.0, a) + 1e-12)
                    sb.Append($" charge rose with lysinesModified {n}→{n + 1};");
                if (n >= nK && Pi(b) != Pi(a)) sb.Append($" lysinesModified {n + 1} ≥ #K={nK} must saturate;");
            }
            for (int n = 0; 2 * n <= nC + 2; n++)
            {
                var a = t.baseOpts with { Disulphides = n };
                var b = t.baseOpts with { Disulphides = n + 1 };
                if (Pi(b) < Pi(a)) sb.Append($" pI fell with disulphides {n}→{n + 1};");
                if (SequenceStatistics.CalculateNetCharge(t.seq, 7.0, b) < SequenceStatistics.CalculateNetCharge(t.seq, 7.0, a) - 1e-12)
                    sb.Append($" charge fell with disulphides {n}→{n + 1};");
                if (2 * n >= nC && Pi(b) != Pi(a)) sb.Append($" disulphides {n + 1} (2N ≥ #C={nC}) must saturate;");
            }
            return (sb.Length == 0).Label($"'{t.seq}' {t.baseOpts}:{sb}");
        });
    }

    /// <summary>
    /// F44 (P): EMBOSS <c>-notermini</c> on a sequence without ionizable side chains (no D/E/C/Y/H/K/R/B/Z) has net
    /// charge 0 at every pH, so <c>embIepPhConverge</c>'s first probe ((14 − 1)/2 + 1 = 7.5) is an exact zero and the
    /// EMBOSS-window pI is 7.5.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IsoelectricPointOptions_NoTermini_ChargeFreeSequence_EmbossPiIsSevenPointFive()
    {
        var gen = from seq in AlphabetArbitrary("AGILMFPSTVWNQXagil", 1, 30).Generator
                  from scale in Gen.Elements(AllPkaScales)
                  select (seq, scale);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var opts = new SequenceStatistics.IsoelectricPointOptions
            {
                Scale = t.scale, Termini = false, Window = SequenceStatistics.IsoelectricPointWindow.Emboss,
            };
            bool zero = Enumerable.Range(0, 29).All(k => SequenceStatistics.CalculateNetCharge(t.seq, k * 0.5, opts) == 0.0);
            bool found = SequenceStatistics.TryCalculateIsoelectricPoint(t.seq, opts, out double pi);
            return (zero && found && pi == 7.5).Label($"'{t.seq}' {t.scale}: zero charge {zero}, found {found}, pI {pi}");
        });
    }

    /// <summary>
    /// F29 (P): for an odd window W ≥ 3, any edge weight and residues that all have a scale value, the Biopython
    /// <c>protein_scale</c> convention (paired positions j and W−1−j with the half-list weight, middle weight 1,
    /// divisor 2·Σ half-weights + 1) equals the standard weighted mean (same weights, Σ w = the same divisor).
    /// Checked for every shipped <c>ProtParamData</c> scale.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HydrophobicityOptions_BiopythonConvention_EqualsStandard_OddWindowKnownResidues()
    {
        var gen = from seq in AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVarndc", 0, 40).Generator
                  from half in Gen.Choose(1, 6)
                  from edgeTenths in Gen.Choose(0, 10)
                  from scale in Gen.Elements(ShippedScales)
                  select (seq, w: 2 * half + 1, edge: edgeTenths / 10.0, scale);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var std = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge,
                new SequenceStatistics.HydrophobicityProfileOptions { Scale = t.scale }).ToList();
            var bio = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge,
                new SequenceStatistics.HydrophobicityProfileOptions { Scale = t.scale, Convention = SequenceStatistics.HydropathyWindowConvention.Biopython }).ToList();
            bool ok = std.Count == bio.Count && std.Count == Math.Max(0, t.seq.Length - t.w + 1)
                      && std.Zip(bio).All(p => Math.Abs(p.First - p.Second) <= 1e-12);
            return ok.Label($"'{t.seq}' W={t.w} edge={t.edge}: std [{string.Join(", ", std)}] vs biopython [{string.Join(", ", bio)}]");
        });
    }

    /// <summary>
    /// F30 (R + P): ExPASy ProtScale 0..1 normalisation v' = (v − min)/(max − min). Every profile value lies in
    /// [0, 1] in both conventions (unknown residues add 0); on known-residue input the standard profile is the
    /// affine image (raw − min)/(max − min) of the raw profile (the weighted mean is affine-equivariant).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HydrophobicityOptions_NormalizedScale_InUnitInterval_AndAffine()
    {
        var gen = from seq in AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVXBZ*", 0, 40).Generator
                  from w in Gen.Choose(2, 11)
                  from edgeTenths in Gen.Choose(0, 10)
                  from scale in Gen.Elements(ShippedScales)
                  select (seq, w, edge: edgeTenths / 10.0, scale);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double edge = t.w % 2 == 0 ? 1.0 : t.edge; // the standard convention needs an odd W for edge < 1
            var sb = new System.Text.StringBuilder();
            foreach (var conv in new[] { SequenceStatistics.HydropathyWindowConvention.Standard, SequenceStatistics.HydropathyWindowConvention.Biopython })
            {
                var norm = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, edge,
                    new SequenceStatistics.HydrophobicityProfileOptions { Scale = t.scale, Convention = conv, NormalizeScale = true });
                if (!norm.All(v => v is >= -1e-12 and <= 1.0 + 1e-12)) sb.Append($" {conv}: value outside [0,1];");
            }
            string known = new(t.seq.Where(c => t.scale.ContainsKey(c)).ToArray());
            double min = t.scale.Values.Min(), max = t.scale.Values.Max();
            var raw = SequenceStatistics.CalculateHydrophobicityProfile(known, t.w, edge, new SequenceStatistics.HydrophobicityProfileOptions { Scale = t.scale }).ToList();
            var nrm = SequenceStatistics.CalculateHydrophobicityProfile(known, t.w, edge, new SequenceStatistics.HydrophobicityProfileOptions { Scale = t.scale, NormalizeScale = true }).ToList();
            if (raw.Count != nrm.Count || raw.Zip(nrm).Any(p => Math.Abs((p.First - min) / (max - min) - p.Second) > 1e-12))
                sb.Append(" normalised ≠ affine image of the raw profile;");
            return (sb.Length == 0).Label($"'{t.seq}' W={t.w} edge={edge}:{sb}");
        });
    }

    /// <summary>
    /// F29/F31/F52 (P): the caller-scale overload with the shipped Kyte–Doolittle table, and the options overload
    /// with <c>Scale = KD</c> or default options, reproduce the default profile exactly; strict mode throws iff a
    /// symbol has no value in the scale and otherwise equals the lenient profile.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HydrophobicityOptions_KdCallerScale_ReproducesDefault_StrictIffUnknown()
    {
        var gen = from seq in AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVarnXB-", 0, 40).Generator
                  from half in Gen.Choose(0, 6)
                  from edgeTenths in Gen.Choose(0, 10)
                  select (seq, w: 2 * half + 1, edge: edgeTenths / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var def = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge).ToList();
            bool sameOk = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, ProteinScales.KyteDoolittle, t.w, t.edge).SequenceEqual(def)
                          && SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge, new SequenceStatistics.HydrophobicityProfileOptions { Scale = ProteinScales.KyteDoolittle }).SequenceEqual(def)
                          && SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge, null).SequenceEqual(def);
            bool hasUnknown = t.seq.Any(c => !ProteinScales.KyteDoolittle.ContainsKey(UpperAscii(c)));
            bool threw = false;
            List<double>? strict = null;
            try { strict = SequenceStatistics.CalculateHydrophobicityProfile(t.seq, t.w, t.edge, new SequenceStatistics.HydrophobicityProfileOptions { Strict = true }).ToList(); }
            catch (ArgumentException) { threw = true; }
            bool strictOk = threw == hasUnknown && (threw || strict!.SequenceEqual(def));
            return (sameOk && strictOk).Label($"'{t.seq}' W={t.w} edge={t.edge}: KD scale identity {sameOk}, strict threw {threw} (unknown {hasUnknown})");
        });
    }

    /// <summary>
    /// F27 (P): default / null <c>ThermodynamicsOptions</c> equal the 4-argument overload; for any salt method 0–7 and
    /// ion mix, the unrounded result rounded (ΔH/ΔS/ΔG to 2 dp, Tm to 1 dp) is exactly the rounded result, and both
    /// calls agree on rejecting a salt input Biopython rejects.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ThermodynamicsOptions_DefaultsEqualOverload_UnroundedRoundsToRounded()
    {
        var gen = from seq in AlphabetArbitrary("ACGT", 2, 40).Generator
                  from self in Gen.Elements(true, false)
                  from naMm in Gen.Choose(0, 1000)
                  from method in Gen.Choose(0, 7)
                  from kMm in Gen.Choose(0, 100)
                  from trisMm in Gen.Choose(0, 50)
                  from mgTenthMm in Gen.Choose(0, 100)
                  from dntpTenthMm in Gen.Choose(0, 20)
                  select (seq, self, na: naMm / 1000.0,
                          opts: new SequenceStatistics.ThermodynamicsOptions
                          {
                              SaltCorrection = (NnSaltCorrection)method, PotassiumConcentration = kMm / 1000.0,
                              TrisConcentration = trisMm / 1000.0, MagnesiumConcentration = mgTenthMm / 10000.0,
                              DntpConcentration = dntpTenthMm / 10000.0,
                          });
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            const double ct = 2.5e-7;
            bool defaultsOk = t.na <= 0 || (
                SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self, null) == SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self)
                && SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self, SequenceStatistics.ThermodynamicsOptions.Default) == SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self));

            SequenceStatistics.ThermodynamicProperties r = default, u = default;
            Type? rErr = null, uErr = null;
            try { r = SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self, t.opts); } catch (ArgumentException e) { rErr = e.GetType(); }
            try { u = SequenceStatistics.CalculateThermodynamics(t.seq, t.na, ct, t.self, t.opts with { RoundResults = false }); } catch (ArgumentException e) { uErr = e.GetType(); }
            bool roundOk = rErr == uErr && (rErr != null || (
                Math.Round(u.DeltaH, 2) == r.DeltaH && Math.Round(u.DeltaS, 2) == r.DeltaS
                && Math.Round(u.DeltaG, 2) == r.DeltaG && Math.Round(u.MeltingTemperature, 1) == r.MeltingTemperature
                && double.IsFinite(u.MeltingTemperature)));
            return (defaultsOk && roundOk).Label($"'{t.seq}' self={t.self} Na={t.na} {t.opts}: defaults {defaultsOk}, rounding {roundOk} (errors {rErr?.Name}/{uErr?.Name}; {u} vs {r})");
        });
    }

    private static Gen<(string seq, int w, int step)> ProfileGen(string alphabet) =>
        from seq in AlphabetArbitrary(alphabet, 0, 60).Generator
        from w in Gen.Choose(1, 15)
        from step in Gen.Choose(1, 6)
        select (seq, w, step);

    private static int PartialStartCount(int len, int step) => len == 0 ? 0 : (len + step - 1) / step;

    /// <summary>
    /// F24/F38 (P): Biopython <c>GC_skew</c>-style windowing (<c>range(0, len, step)</c>). With trailing partial windows
    /// the GC and entropy profiles start with exactly the complete-window profile, have ⌈len/step⌉ values (one per
    /// start &lt; len), and every partial value is the statistic of the suffix <c>seq[start:]</c>. Checked for every GC
    /// ambiguity mode and the letter / k = 1 entropy profiles.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcAndEntropyProfiles_PartialWindows_PrefixCountAndSuffixValues()
    {
        return Prop.ForAll(ProfileGen("ACGTNacgt-RS").ToArbitrary(), t =>
        {
            int expectedCount = PartialStartCount(t.seq.Length, t.step);
            var sb = new System.Text.StringBuilder();
            void Check(string name, List<double> complete, List<double> partial, Func<string, double> stat)
            {
                if (partial.Count != expectedCount) sb.Append($" {name}: count {partial.Count} ≠ ⌈{t.seq.Length}/{t.step}⌉;");
                if (!partial.Take(complete.Count).SequenceEqual(complete)) sb.Append($" {name}: prefix ≠ complete profile;");
                for (int k = complete.Count; k < partial.Count; k++)
                {
                    int start = k * t.step;
                    if (start >= t.seq.Length || start + t.w <= t.seq.Length) { sb.Append($" {name}: start {start} is not a partial window;"); break; }
                    if (Math.Abs(partial[k] - stat(t.seq[start..])) > 1e-12) { sb.Append($" {name}: partial {k} ≠ statistic of suffix;"); break; }
                }
            }
            Check("gc", SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, false).ToList(),
                SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, false, includePartialWindow: true).ToList(),
                s => s.AsSpan().CalculateGcFraction() * 100.0);
            foreach (var mode in Enum.GetValues<SequenceExtensions.GcAmbiguityMode>())
                Check($"gc-{mode}", SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, true, mode).ToList(),
                    SequenceStatistics.CalculateGcContentProfile(t.seq, t.w, t.step, true, mode, includePartialWindow: true).ToList(),
                    s => s.AsSpan().CalculateGcFraction(mode));
            Check("entropy", SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step).ToList(),
                SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step, includePartialWindow: true).ToList(),
                SequenceStatistics.CalculateShannonEntropy);
            Check("kmer-entropy", SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step, 2, false).ToList(),
                SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step, 2, true).ToList(),
                s => SequenceComplexity.CalculateKmerEntropy(s, 2));
            return (sb.Length == 0).Label($"'{t.seq}' W={t.w} step={t.step}:{sb}");
        });
    }

    /// <summary>
    /// F38 (R + P): the order-k block-entropy profile of a DNA window of length L holds L − k + 1 overlapping k-mers
    /// over at most 4^k words, so each value lies in [0, log₂ min(4^k, L − k + 1)] (0 when L &lt; k); for k = 1 it
    /// equals the letter entropy profile on A/C/G/T input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerEntropyProfile_BoundedByLog2DistinctWords_K1EqualsLetterProfile()
    {
        var gen = from t in ProfileGen("ACGT")
                  from k in Gen.Choose(1, 5)
                  from partial in Gen.Elements(true, false)
                  select (t.seq, t.w, t.step, k, partial);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var profile = SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step, t.k, t.partial).ToList();
            var starts = new List<int>();
            for (int i = 0; i < t.seq.Length; i += t.step)
                if (i + t.w <= t.seq.Length || t.partial) starts.Add(i);
            if (profile.Count != starts.Count)
                return false.Label($"count {profile.Count} ≠ {starts.Count}");
            for (int j = 0; j < starts.Count; j++)
            {
                int len = Math.Min(t.w, t.seq.Length - starts[j]);
                int words = len - t.k + 1;
                double bound = words <= 1 ? 0.0 : Math.Log2(Math.Min(Math.Pow(4, t.k), words));
                if (profile[j] < -1e-12 || profile[j] > bound + 1e-9)
                    return false.Label($"'{t.seq}' W={t.w} k={t.k}: window {j} (L={len}) entropy {profile[j]} ∉ [0, {bound}]");
            }
            if (t.k == 1)
            {
                var letters = SequenceStatistics.CalculateEntropyProfile(t.seq, t.w, t.step, t.partial).ToList();
                if (!letters.Zip(profile).All(p => Math.Abs(p.First - p.Second) <= 1e-12) || letters.Count != profile.Count)
                    return false.Label($"'{t.seq}' W={t.w}: k=1 block entropy ≠ letter entropy profile");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// F34/F45 (R + INV): the Chou–Fasman 1978 graded assignment returns one state per input character from
    /// {H, E, T, C}; characters without parameters (non-standard, gaps, non-ASCII letters) stay coil; the ASCII
    /// case fold makes lower-case input equivalent.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ChouFasmanGraded1978_LengthAlphabetCoilForNonResidues_CaseFold()
    {
        return Prop.ForAll(AlphabetArbitrary("ARNDCEQGHILKMFPSTWYVAEMLKVIYXBZ-*ſı", 0, 80), seq =>
        {
            string states = SequenceStatistics.PredictSecondaryStructureChouFasman(seq, SequenceStatistics.ChouFasmanRuleSet.Graded1978);
            bool lengthOk = states.Length == seq.Length;
            bool alphabetOk = states.All(c => c is 'H' or 'E' or 'T' or 'C');
            bool coilOk = lengthOk && seq.Select((c, i) => (c, i)).All(p => StandardResidueSet.Contains(p.c) || states[p.i] == 'C');
            string lower = new(seq.Select(c => char.IsAsciiLetterUpper(c) ? (char)(c + 32) : c).ToArray());
            bool caseOk = SequenceStatistics.PredictSecondaryStructureChouFasman(lower, SequenceStatistics.ChouFasmanRuleSet.Graded1978) == states;
            return (lengthOk && alphabetOk && coilOk && caseOk)
                .Label($"'{seq}' → '{states}': length {lengthOk}, alphabet {alphabetOk}, non-residues coil {coilOk}, case {caseOk}");
        });
    }

    private const string StandardResidueSet = "ARNDCEQGHILKMFPSTWYV";

    /// <summary>
    /// F36/F49 (P): EMBOSS <c>cusp</c> normalisations. PerThousand = 1000 × Total on the same keys; PerAminoAcid
    /// ("Fraction") sums to 1 over the observed codons of every synonymous family (amino acid, or '*' for stops)
    /// under the chosen genetic code, each value in (0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonFrequencyNormalizations_PerThousandAndPerAminoAcidFamilies()
    {
        var codes = new[] { GeneticCode.Standard, GeneticCode.VertebrateMitochondrial, GeneticCode.YeastMitochondrial, GeneticCode.BacterialPlastid };
        var gen = from seq in AlphabetArbitrary("ACGTUacgtN", 0, 90).Generator
                  from frame in Gen.Choose(0, 3)
                  from code in Gen.Elements(codes)
                  select (seq, frame, code);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var total = SequenceStatistics.CalculateCodonFrequencies(t.seq, t.frame);
            var totalExplicit = SequenceStatistics.CalculateCodonFrequencies(t.seq, t.frame, SequenceStatistics.CodonFrequencyNormalization.Total, t.code);
            var perK = SequenceStatistics.CalculateCodonFrequencies(t.seq, t.frame, SequenceStatistics.CodonFrequencyNormalization.PerThousand, t.code);
            var perAa = SequenceStatistics.CalculateCodonFrequencies(t.seq, t.frame, SequenceStatistics.CodonFrequencyNormalization.PerAminoAcid, t.code);
            bool totalOk = totalExplicit.Count == total.Count && total.All(kv => totalExplicit[kv.Key] == kv.Value);
            bool perKOk = perK.Count == total.Count && total.All(kv => perK.TryGetValue(kv.Key, out double v) && Math.Abs(v - 1000.0 * kv.Value) <= 1e-9);
            bool keysOk = perAa.Count == total.Count && total.Keys.All(perAa.ContainsKey);
            bool familiesOk = perAa.All(kv => kv.Value is > 0.0 and <= 1.0)
                              && perAa.GroupBy(kv => t.code.Translate(kv.Key)).All(g => Math.Abs(g.Sum(kv => kv.Value) - 1.0) <= 1e-12);
            return (totalOk && perKOk && keysOk && familiesOk)
                .Label($"'{t.seq}' frame {t.frame}: Total {totalOk}, PerThousand {perKOk}, keys {keysOk}, families sum to 1 {familiesOk}");
        });
    }

    /// <summary>
    /// F51 (P): <c>includeAbsent</c> is a superset of the default ratio map with identical values for observed
    /// dinucleotides, every added key has ρ = 0, and the key set is exactly B × B for the set B of counted bases (the
    /// sequence's A/C/G/T, plus their complements for Karlin's strand-symmetric ρ*).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DinucleotideRatios_IncludeAbsent_SupersetWithZerosOverBasePairs()
    {
        var gen = from seq in AlphabetArbitrary("ACGTNacgt-", 0, 60).Generator
                  from sym in Gen.Elements(true, false)
                  select (seq, sym);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var observed = SequenceStatistics.CalculateDinucleotideRatios(t.seq, t.sym);
            var all = SequenceStatistics.CalculateDinucleotideRatios(t.seq, t.sym, includeAbsent: true);
            if (observed.Count == 0)
                return (all.Count == 0).Label($"'{t.seq}': no observed dinucleotide ⇒ empty, got {all.Count}");
            var bases = t.seq.Select(UpperAscii).Where(c => "ACGT".Contains(c)).ToHashSet();
            if (t.sym)
                bases.UnionWith(bases.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToList());
            var expectedKeys = bases.SelectMany(x => bases.Select(y => $"{x}{y}")).ToHashSet();
            bool superset = observed.All(kv => all.TryGetValue(kv.Key, out double v) && v == kv.Value);
            bool zeros = all.Where(kv => !observed.ContainsKey(kv.Key)).All(kv => kv.Value == 0.0);
            bool keys = all.Keys.ToHashSet().SetEquals(expectedKeys);
            return (superset && zeros && keys).Label($"'{t.seq}' ρ*={t.sym}: superset {superset}, zeros {zeros}, keys B×B {keys} ({string.Join(",", all.Keys)})");
        });
    }

    /// <summary>
    /// F33 (P): <c>SequenceSummary.MolecularWeight</c> is the canonical single-stranded linear average MW — RNA table
    /// when the sequence has U and no T, DNA table otherwise — for both Tm weighting flags; and on input without IUPAC
    /// ambiguity codes (only A/C/G/T/U, any case) the IUPAC-weighted Tm (Biopython non-strict) equals the default Tm.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Summary_MolecularWeightCanonical_WeightedTmEqualsDefaultOnPureBases()
    {
        var gen = from seq in AlphabetArbitrary("ACGTUacgtu", 0, 40).Generator
                  from junk in Gen.Elements("", "N", "-", "NN-")
                  select (seq, junk);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string withJunk = t.seq + t.junk;
            bool hasU = withJunk.Any(c => c is 'U' or 'u'), hasT = withJunk.Any(c => c is 'T' or 't');
            double expectedMw = SequenceStatistics.CalculateNucleotideMolecularWeight(withJunk, isDna: !(hasU && !hasT));
            bool mwOk = SequenceStatistics.SummarizeNucleotideSequence(withJunk).MolecularWeight == expectedMw
                        && SequenceStatistics.SummarizeNucleotideSequence(withJunk, weightAmbiguousBasesInTm: true).MolecularWeight == expectedMw;
            bool tmOk = new[] { true, false }.All(wallace =>
                Math.Abs(SequenceStatistics.CalculateMeltingTemperature(t.seq, wallace, weightAmbiguousBases: true)
                         - SequenceStatistics.CalculateMeltingTemperature(t.seq, wallace)) <= 1e-9);
            bool summaryTmOk = Math.Abs(SequenceStatistics.SummarizeNucleotideSequence(t.seq, true).MeltingTemperature
                                        - SequenceStatistics.SummarizeNucleotideSequence(t.seq).MeltingTemperature) <= 1e-9;
            return (mwOk && tmOk && summaryTmOk).Label($"'{withJunk}': MW {mwOk}, weighted Tm = default {tmOk}, summary Tm {summaryTmOk}");
        });
    }

    #endregion
}
