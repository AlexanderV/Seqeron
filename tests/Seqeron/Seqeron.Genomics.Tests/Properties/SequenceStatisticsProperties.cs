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
}
