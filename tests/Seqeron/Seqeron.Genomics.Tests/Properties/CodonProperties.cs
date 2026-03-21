using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for codon tables, translation, and codon usage.
/// Verifies invariants of the genetic code and translation process.
///
/// Test Units: TRANS-CODON-001, TRANS-PROT-001, CODON-USAGE-001, CODON-CAI-001, CODON-OPT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Translation")]
public class CodonProperties
{
    private static Arbitrary<string> CodingDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 3)
            .Select(a => new string(a, 0, a.Length - a.Length % 3)) // trim to codon boundary
            .Where(s => s.Length >= 3)
            .ToArbitrary();

    /// <summary>
    /// The standard genetic code must have exactly 64 codons (4³).
    /// </summary>
    [Test]
    [Category("Property")]
    public void StandardCode_Has64Codons()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.CodonTable.Count, Is.EqualTo(64),
            "Standard genetic code must map all 64 codons");
    }

    /// <summary>
    /// Every codon in the table must be exactly 3 characters.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AllCodons_HaveLength3()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.CodonTable.Keys.All(c => c.Length == 3), Is.True,
            "All codons must be exactly 3 nucleotides");
    }

    /// <summary>
    /// Start and stop codon sets must be non-empty and disjoint.
    /// </summary>
    [Test]
    [Category("Property")]
    public void StartAndStopCodons_AreDisjoint()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.StartCodons.Count, Is.GreaterThan(0));
        Assert.That(code.StopCodons.Count, Is.GreaterThan(0));
        Assert.That(code.StartCodons.Intersect(code.StopCodons).Count(), Is.EqualTo(0),
            "Start and stop codon sets must be disjoint");
    }

    /// <summary>
    /// ATG must be a start codon in the standard code.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ATG_IsStartCodon()
    {
        Assert.That(GeneticCode.Standard.IsStartCodon("ATG"), Is.True);
    }

    /// <summary>
    /// TAA, TAG, TGA must be stop codons in the standard code.
    /// </summary>
    [TestCase("TAA")]
    [TestCase("TAG")]
    [TestCase("TGA")]
    [Category("Property")]
    public void KnownStopCodons_AreStopCodons(string codon)
    {
        Assert.That(GeneticCode.Standard.IsStopCodon(codon), Is.True);
    }

    /// <summary>
    /// Translation preserves length: protein length == DNA length / 3.
    /// (Excluding stop codons which produce '*')
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_OutputLength_MatchesInput()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            int expectedLen = seq.Length / 3;
            return (protein.Sequence.Length == expectedLen)
                .Label($"Protein len={protein.Sequence.Length}, expected={expectedLen}");
        });
    }

    /// <summary>
    /// Translation output only contains valid amino acid characters or '*' (stop).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_ProducesValidAminoAcids()
    {
        const string validAa = "ACDEFGHIKLMNPQRSTVWY*";
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            return protein.Sequence.All(c => validAa.Contains(c))
                .Label($"Invalid amino acid in: {protein.Sequence[..Math.Min(20, protein.Sequence.Length)]}");
        });
    }

    /// <summary>
    /// Codon usage counts sum to total number of codons in the sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonUsage_CountsSumToTotal()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var usage = CodonOptimizer.CalculateCodonUsage(seq);
            int total = usage.Values.Sum();
            int expected = seq.Length / 3;
            return (total == expected)
                .Label($"Usage sum={total}, expected={expected}");
        });
    }

    /// <summary>
    /// CAI is always in (0, 1] for valid coding sequences.
    /// Evidence: CAI = geometric mean of relative adaptiveness values, each in (0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CAI_InRange()
    {
        // Generate coding sequences that are multiples of 3, contain only ACGT, 
        // and have at least one non-stop codon
        var codingGen = Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 3)
            .Select(a => new string(a, 0, a.Length - a.Length % 3))
            .Where(s => s.Length >= 3)
            // Filter out sequences that are all stop codons
            .Where(s =>
            {
                var code = GeneticCode.Standard;
                for (int i = 0; i < s.Length - 2; i += 3)
                {
                    string codon = s.Substring(i, 3);
                    if (!code.IsStopCodon(codon)) return true;
                }
                return false;
            })
            .ToArbitrary();

        return Prop.ForAll(codingGen, seq =>
        {
            double cai = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);
            return (cai >= 0.0 && cai <= 1.0 + 0.0001)
                .Label($"CAI={cai:F4}, must be in [0, 1]");
        });
    }

    /// <summary>
    /// GetCodonsForAminoAcid returns codons that all translate back to the same amino acid.
    /// </summary>
    [TestCase('M')]
    [TestCase('W')]
    [TestCase('L')]
    [TestCase('S')]
    [Category("Property")]
    public void Codons_TranslateBackToSameAminoAcid(char aminoAcid)
    {
        var code = GeneticCode.Standard;
        var codons = code.GetCodonsForAminoAcid(aminoAcid).ToList();
        Assert.That(codons, Is.Not.Empty, $"Must have codons for '{aminoAcid}'");
        foreach (var codon in codons)
            Assert.That(code.Translate(codon), Is.EqualTo(aminoAcid),
                $"Codon {codon} should translate to {aminoAcid}");
    }

    #region CODON-OPT-001: P: optimized translates to same protein; R: only valid codons; D: deterministic

    /// <summary>
    /// INV-1: Optimized sequence translates to the same protein as the original.
    /// Evidence: Codon optimization replaces synonymous codons only — the amino acid
    /// sequence encoded by the DNA must be preserved exactly.
    /// This is the fundamental invariant of codon optimization.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CodonOptimization_PreservesProteinSequence()
    {
        // A known coding sequence (GFP fragment)
        string codingSeq = "AUGGCUAGCAAAGGA";
        var result = CodonOptimizer.OptimizeSequence(codingSeq, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Translate both original and optimized
        string originalProtein = TranslateRna(result.OriginalSequence);
        string optimizedProtein = TranslateRna(result.OptimizedSequence);

        Assert.That(optimizedProtein, Is.EqualTo(originalProtein),
            $"Optimized protein '{optimizedProtein}' ≠ original '{originalProtein}'");
    }

    /// <summary>
    /// INV-2: Optimized sequence translates to same protein for arbitrary coding DNA.
    /// Evidence: Synonymous codon substitution preserves translation (by definition of the genetic code).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_PreservesProtein_Property()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            if (result.OptimizedSequence.Length == 0) return true.ToProperty();

            string origProtein = TranslateRna(result.OriginalSequence);
            string optProtein = TranslateRna(result.OptimizedSequence);
            return (origProtein == optProtein)
                .Label($"Protein mismatch: orig='{origProtein[..Math.Min(10, origProtein.Length)]}' ≠ opt='{optProtein[..Math.Min(10, optProtein.Length)]}'");
        });
    }

    /// <summary>
    /// INV-3: Optimized sequence contains only valid RNA codons (A, C, G, U).
    /// Evidence: CodonOptimizer works with RNA codons internally.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_OnlyValidCodons()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            if (result.OptimizedSequence.Length == 0) return true.ToProperty();

            bool valid = result.OptimizedSequence.All(c => "ACGU".Contains(c));
            return valid.Label($"Invalid chars in optimized: {result.OptimizedSequence[..Math.Min(20, result.OptimizedSequence.Length)]}");
        });
    }

    /// <summary>
    /// INV-4: Optimized sequence length equals original sequence length.
    /// Evidence: Codon optimization replaces codons 1:1, preserving total length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_PreservesLength()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            return (result.OptimizedSequence.Length == result.OriginalSequence.Length)
                .Label($"Length changed: orig={result.OriginalSequence.Length}, opt={result.OptimizedSequence.Length}");
        });
    }

    /// <summary>
    /// INV-5: Codon optimization is deterministic.
    /// Evidence: OptimizeSequence is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonOptimization_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var r1 = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            var r2 = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
            return (r1.OptimizedSequence == r2.OptimizedSequence &&
                    r1.OptimizedCAI == r2.OptimizedCAI)
                .Label("OptimizeSequence must be deterministic");
        });
    }

    #endregion

    #region CODON-CAI-001: R: CAI ∈ [0,1]; M: all optimal codons → CAI close to 1.0; D: deterministic

    /// <summary>
    /// INV-CAI1: A sequence composed entirely of optimal codons produces CAI close to 1.0.
    /// Evidence: CAI = geometric mean of w_i where w_i = freq(codon) / max_freq(AA).
    /// When every codon is the most frequent for its AA, each w_i = 1, so CAI = 1.
    /// Source: Sharp &amp; Li (1987) "The codon adaptation index".
    /// </summary>
    [Test]
    [Category("Property")]
    public void CAI_AllOptimalCodons_HighCAI()
    {
        // Build a sequence using only the most frequent codon for each amino acid in E. coli
        var table = CodonOptimizer.EColiK12;
        var optimalCodons = new List<string>();

        // Select some amino acids and their optimal codons
        string[] aminoAcids = { "L", "S", "A", "G", "V", "T", "P" };
        foreach (var aa in aminoAcids)
        {
            // Find the codon with highest frequency for this AA from the CodonFrequencies
            string best = table.CodonFrequencies
                .Where(kv => table.CodonToAminoAcid.GetValueOrDefault(kv.Key) == aa)
                .OrderByDescending(kv => kv.Value)
                .First().Key;
            optimalCodons.Add(best);
        }

        string optimalSeq = string.Join("", optimalCodons);
        double cai = CodonOptimizer.CalculateCAI(optimalSeq, table);

        Assert.That(cai, Is.GreaterThanOrEqualTo(0.9),
            $"Sequence of optimal codons should have CAI ≈ 1.0, got {cai:F4}");
    }

    /// <summary>
    /// INV-CAI2: CAI calculation is deterministic for any coding sequence.
    /// Evidence: CalculateCAI is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CAI_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            double cai1 = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);
            double cai2 = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);
            return (cai1 == cai2)
                .Label($"CAI must be deterministic: {cai1:F6} vs {cai2:F6}");
        });
    }

    /// <summary>
    /// INV-CAI3: CodonUsageAnalyzer.CalculateCai also returns values in [0, 1].
    /// Evidence: CAI is a geometric mean of ratios in (0, 1], bounded [0, 1].
    /// Source: Sharp &amp; Li (1987).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CAI_CodonUsageAnalyzer_InRange()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            double cai = CodonUsageAnalyzer.CalculateCai(seq, CodonUsageAnalyzer.EColiOptimalCodons);
            return (cai >= 0.0 && cai <= 1.0 + 0.0001)
                .Label($"CAI={cai:F4} must be in [0, 1]");
        });
    }

    #endregion

    #region CODON-USAGE-001: R: usage freqs ≥ 0; P: sum per amino acid = k; D: deterministic

    /// <summary>
    /// INV-U1: All RSCU values are non-negative.
    /// Evidence: RSCU = observed / expected where both are ≥ 0.
    /// Source: Sharp &amp; Li (1986) "An evolutionary perspective on synonymous codon usage".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RSCU_AllValues_NonNegative()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var rscu = CodonUsageAnalyzer.CalculateRscu(seq);
            bool allNonNeg = rscu.Values.All(v => v >= 0.0);
            return allNonNeg.Label($"All RSCU values must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-U2: For each amino acid with observed usage, RSCU values sum to the number
    /// of synonymous codons (degeneracy k).
    /// Evidence: RSCU_i = (observed_i × k) / total_AA. Sum = k × total_AA / total_AA = k.
    /// Source: Sharp &amp; Li (1986).
    /// </summary>
    [Test]
    [Category("Property")]
    public void RSCU_SumPerAminoAcid_EqualsExpected()
    {
        string seq = "ATGATGAAATTCTTACTGCCCCCCACCACCGCTGCTGCAGCAGGCAGCG";
        var rscu = CodonUsageAnalyzer.CalculateRscu(seq);
        var code = GeneticCode.Standard;

        // GeneticCode.Standard uses RNA codons (GCU); RSCU uses DNA codons (GCT)
        foreach (var aaGroup in code.CodonTable.GroupBy(kv => kv.Value))
        {
            var synonymousCodons = aaGroup.Select(kv => kv.Key.Replace('U', 'T')).ToList();
            int k = synonymousCodons.Count;
            double sum = synonymousCodons.Sum(c => rscu.GetValueOrDefault(c, 0.0));

            bool anyUsed = synonymousCodons.Any(c => rscu.GetValueOrDefault(c, 0.0) > 0);
            if (anyUsed)
            {
                Assert.That(sum, Is.EqualTo(k).Within(0.01),
                    $"AA '{aaGroup.Key}': RSCU sum={sum:F3}, expected={k}");
            }
            else
            {
                Assert.That(sum, Is.EqualTo(0).Within(0.01),
                    $"AA '{aaGroup.Key}': unused amino acid should have RSCU sum=0");
            }
        }
    }

    /// <summary>
    /// INV-U3: RSCU calculation is deterministic.
    /// Evidence: CalculateRscu is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RSCU_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var rscu1 = CodonUsageAnalyzer.CalculateRscu(seq);
            var rscu2 = CodonUsageAnalyzer.CalculateRscu(seq);
            bool same = rscu1.Count == rscu2.Count &&
                        rscu1.All(kv => rscu2.ContainsKey(kv.Key) &&
                                        Math.Abs(kv.Value - rscu2[kv.Key]) < 1e-10);
            return same.Label("CalculateRscu must be deterministic");
        });
    }

    /// <summary>
    /// INV-U4: Codon counts are all non-negative and sum to total number of codons.
    /// Evidence: CountCodons tallies occurrences, each ≥ 0; total = seqLen / 3.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CodonCounts_AllNonNegative_SumToTotal()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var counts = CodonUsageAnalyzer.CountCodons(seq);
            bool allNonNeg = counts.Values.All(c => c >= 0);
            int total = counts.Values.Sum();
            int expected = seq.Length / 3;
            return (allNonNeg && total == expected)
                .Label($"counts sum={total}, expected={expected}, allNonNeg={allNonNeg}");
        });
    }

    #endregion

    #region CODON-RARE-001: R: rare codon positions valid; M: lower threshold → more rare codons; D: deterministic

    /// <summary>
    /// INV-RARE1: Rare codon positions are valid nucleotide positions (multiples of 3).
    /// Evidence: FindRareCodons returns position = codon_index × 3.
    /// Source: Sharp et al. (1987).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindRareCodons_Positions_AreMultiplesOf3()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var rare = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, 0.15).ToList();
            return rare.All(r => r.Position % 3 == 0 && r.Position >= 0 && r.Position < seq.Length)
                .Label("All rare codon positions must be multiples of 3 within sequence bounds");
        });
    }

    /// <summary>
    /// INV-RARE2: Rare codon positions are within sequence bounds [0, seqLen - 3].
    /// Evidence: Position refers to the start of a codon triplet.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindRareCodons_Positions_WithinBounds()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var rare = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, 0.15).ToList();
            return rare.All(r => r.Position >= 0 && r.Position + 3 <= seq.Length)
                .Label("Rare codon position + 3 must not exceed sequence length");
        });
    }

    /// <summary>
    /// INV-RARE3: Lower threshold yields more or equal rare codons (monotonicity).
    /// Evidence: threshold₂ &lt; threshold₁ ⇒ all codons below threshold₁ are still below threshold₂,
    /// plus some additional codons with frequency in [threshold₂, threshold₁).
    /// Source: Plotkin &amp; Kudla (2011).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindRareCodons_LowerThreshold_MoreResults()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 6) return true.ToProperty();
            int highCount = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, 0.10).Count();
            int lowCount = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, 0.20).Count();
            return (lowCount >= highCount)
                .Label($"threshold 0.20 → {lowCount} rare codons must be ≥ threshold 0.10 → {highCount}");
        });
    }

    /// <summary>
    /// INV-RARE4: Each reported rare codon has frequency below the threshold.
    /// Evidence: FindRareCodons filters codons where frequency &lt; threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindRareCodons_ReportedFrequency_BelowThreshold()
    {
        const double threshold = 0.15;
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var rare = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, threshold).ToList();
            return rare.All(r => r.Frequency < threshold)
                .Label($"All reported frequencies must be < {threshold}");
        });
    }

    /// <summary>
    /// INV-RARE5: FindRareCodons is deterministic.
    /// Evidence: Pure function — same input always yields same output.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindRareCodons_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var r1 = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12).ToList();
            var r2 = CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12).ToList();
            bool same = r1.Count == r2.Count &&
                        r1.Zip(r2).All(p => p.First.Position == p.Second.Position &&
                                             p.First.Codon == p.Second.Codon);
            return same.Label("FindRareCodons must be deterministic");
        });
    }

    #endregion

    #region TRANS-CODON-001: R: 64 codons mapped; P: start codons → M; P: stop codons → *; D: deterministic

    /// <summary>
    /// INV-TC1: Standard genetic code maps all 64 possible codons.
    /// Evidence: 4³ = 64 possible RNA triplets from {A, U, G, C}.
    /// Source: Wikipedia "Genetic code"; NCBI Translation Tables.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneticCode_Has64Entries()
    {
        var code = GeneticCode.Standard;
        Assert.That(code.CodonTable.Count, Is.EqualTo(64),
            "Standard genetic code must have exactly 64 codon entries");
    }

    /// <summary>
    /// INV-TC2: All start codons (ATG in standard code) translate to Methionine ('M').
    /// Evidence: In the standard code, ATG is both start and Met codon.
    /// Source: Wikipedia "Start codon"; NCBI Translation Tables.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneticCode_StartCodons_TranslateToMet()
    {
        var code = GeneticCode.Standard;
        // ATG is the universal start codon in the standard code
        Assert.That(code.Translate("ATG"), Is.EqualTo('M'),
            "ATG must translate to M (Methionine) in standard code");
    }

    /// <summary>
    /// INV-TC3: All stop codons translate to '*' (termination).
    /// Evidence: Stop codons (UAA, UAG, UGA) signal translation termination.
    /// Source: Wikipedia "Stop codon"; NCBI Translation Tables.
    /// </summary>
    [TestCase("TAA")]
    [TestCase("TAG")]
    [TestCase("TGA")]
    [Category("Property")]
    public void GeneticCode_StopCodons_TranslateToStar(string stopCodon)
    {
        Assert.That(GeneticCode.Standard.Translate(stopCodon), Is.EqualTo('*'),
            $"Stop codon {stopCodon} must translate to '*'");
    }

    /// <summary>
    /// INV-TC4: DNA/RNA equivalence — ATG and AUG produce identical translation.
    /// Evidence: T→U conversion is internal; both notations represent the same codon.
    /// Source: Wikipedia "Genetic code" — DNA/RNA equivalence.
    /// </summary>
    [TestCase("ATG", "AUG")]
    [TestCase("TAA", "UAA")]
    [TestCase("GCT", "GCU")]
    [Category("Property")]
    public void GeneticCode_DnaRna_Equivalence(string dnaCodon, string rnaCodon)
    {
        char dnaResult = GeneticCode.Standard.Translate(dnaCodon);
        char rnaResult = GeneticCode.Standard.Translate(rnaCodon);
        Assert.That(dnaResult, Is.EqualTo(rnaResult),
            $"{dnaCodon}→'{dnaResult}' must equal {rnaCodon}→'{rnaResult}'");
    }

    /// <summary>
    /// INV-TC5: Every codon maps to exactly one deterministic amino acid.
    /// Evidence: The genetic code is a function — each triplet produces one result.
    /// Source: Crick (1968).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneticCode_Translation_IsDeterministic()
    {
        var code = GeneticCode.Standard;
        const string validAa = "ACDEFGHIKLMNPQRSTVWY*";

        foreach (var kv in code.CodonTable)
        {
            char aa1 = code.Translate(kv.Key);
            char aa2 = code.Translate(kv.Key);
            Assert.That(aa1, Is.EqualTo(aa2),
                $"Codon {kv.Key} must translate deterministically");
            Assert.That(validAa.Contains(aa1), Is.True,
                $"Codon {kv.Key} translated to invalid character '{aa1}'");
        }
    }

    /// <summary>
    /// INV-TC6: Degeneracy — each amino acid has at least one codon that translates back to it.
    /// Evidence: The genetic code is degenerate: most amino acids have 2-6 synonymous codons.
    /// Source: Wikipedia "Genetic code" — degeneracy table.
    /// </summary>
    [TestCase('A', 4)]
    [TestCase('L', 6)]
    [TestCase('M', 1)]
    [TestCase('W', 1)]
    [TestCase('S', 6)]
    [Category("Property")]
    public void GeneticCode_Degeneracy_CodonsPerAminoAcid(char aminoAcid, int expectedCount)
    {
        var codons = GeneticCode.Standard.GetCodonsForAminoAcid(aminoAcid).ToList();
        Assert.That(codons.Count, Is.EqualTo(expectedCount),
            $"Amino acid '{aminoAcid}' should have {expectedCount} codons, got {codons.Count}");
    }

    #endregion

    #region TRANS-PROT-001: R: protein len ≤ seqLen/3; P: starts with M if starts with ATG; D: deterministic

    /// <summary>
    /// INV-TP1: Protein length equals seqLen / 3 (exact division, no stop truncation).
    /// Evidence: Translation reads non-overlapping triplets from the entire sequence.
    /// Source: Wikipedia "Translation (biology)".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_ProteinLength_EqualsSeqDivBy3()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            int expected = seq.Length / 3;
            return (protein.Sequence.Length == expected)
                .Label($"Protein len={protein.Sequence.Length}, expected=seqLen/3={expected}");
        });
    }

    /// <summary>
    /// INV-TP2: If sequence starts with ATG, protein starts with 'M' (Methionine).
    /// Evidence: ATG is the universal start codon encoding Met.
    /// Source: Wikipedia "Start codon".
    /// </summary>
    [Test]
    [Category("Property")]
    public void Translation_AtgStart_ProducesMetStart()
    {
        string seq = "ATGGCTGCTGCTGCTGCT"; // ATG + 5 GCT (=Ala) codons
        var protein = Translator.Translate(seq);

        Assert.That(protein.Sequence[0], Is.EqualTo('M'),
            $"Sequence starting with ATG must produce protein starting with M, got '{protein.Sequence[0]}'");
    }

    /// <summary>
    /// INV-TP3: ATG-starting sequences always produce M-starting proteins (FsCheck).
    /// Evidence: Same invariant as INV-TP2, verified with random inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_AtgStart_AlwaysProducesMet()
    {
        // Generate sequences that start with ATG
        var atgArb = Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= 3)
            .Select(a => "ATG" + new string(a, 0, a.Length - a.Length % 3))
            .Where(s => s.Length >= 6)
            .ToArbitrary();

        return Prop.ForAll(atgArb, seq =>
        {
            var protein = Translator.Translate(seq);
            return (protein.Sequence.Length > 0 && protein.Sequence[0] == 'M')
                .Label($"ATG-starting sequence must produce M-starting protein, got '{(protein.Sequence.Length > 0 ? protein.Sequence[0] : '?')}'");
        });
    }

    /// <summary>
    /// INV-TP4: Translation is deterministic — same input always yields same protein.
    /// Evidence: Translator.Translate is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_IsDeterministic()
    {
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var p1 = Translator.Translate(seq);
            var p2 = Translator.Translate(seq);
            return (p1.Sequence == p2.Sequence)
                .Label("Translation must produce identical protein for identical input");
        });
    }

    /// <summary>
    /// INV-TP5: Translation output contains only valid amino acid characters or '*'.
    /// Evidence: The genetic code maps each codon to one of 20 amino acids or a stop signal.
    /// Source: Wikipedia "Genetic code".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Translation_OutputCharacters_AreValidAminoAcids()
    {
        const string validAa = "ACDEFGHIKLMNPQRSTVWY*";
        return Prop.ForAll(CodingDnaArbitrary(), seq =>
        {
            if (seq.Length < 3) return true.ToProperty();
            var protein = Translator.Translate(seq);
            return protein.Sequence.All(c => validAa.Contains(c))
                .Label($"Invalid character in protein: {protein.Sequence[..Math.Min(20, protein.Sequence.Length)]}");
        });
    }

    #endregion

    private static string TranslateRna(string rna)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i + 2 < rna.Length; i += 3)
        {
            string codon = rna.Substring(i, 3);
            // Convert to DNA for GeneticCode
            string dnaCodon = codon.Replace('U', 'T');
            sb.Append(GeneticCode.Standard.Translate(dnaCodon));
        }
        return sb.ToString();
    }
}
