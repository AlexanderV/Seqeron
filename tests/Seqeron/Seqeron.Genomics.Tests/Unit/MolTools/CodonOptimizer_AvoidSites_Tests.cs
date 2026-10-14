namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for CODON-OPT-001: the additive <c>avoidSites</c> pass of
/// <c>CodonOptimizer.OptimizeSequence(...)</c> and the correctly spelled
/// <c>OptimizationStrategy.AvoidRareCodons</c> alias (review 2026-09, B02 F54/F55,
/// audit items A2-13 / A2-14).
///
/// Reference: DNA Chisel 3.2.16 composes codon optimisation with site avoidance in ONE problem —
/// <c>DnaOptimizationProblem(sequence, constraints=[AvoidPattern("CTGCAG"), AvoidPattern("CATATG"),
/// AvoidPattern("CCATGG"), EnforceTranslation()], objectives=[CodonOptimize(species='e_coli',
/// method='use_best_codon')])</c> — and resolves the hard constraints last. Every locked sequence
/// below was fed back into that library: all three <c>AvoidPattern(...).evaluate(problem)</c> score
/// 0 (pass) and <c>EnforceTranslation</c> passes, and for these three genes dnachisel's own solver
/// returns the byte-identical sequence with the same <c>CodonOptimize</c> objective score −1.3892
/// (docs/Evidence/CODON-OPT-001-Evidence.md, review 2026-09 F54).
/// </summary>
[TestFixture]
[Category("Codon")]
[Category("CODON-OPT-001")]
public class CodonOptimizer_AvoidSites_Tests
{
    // PstI, NdeI, NcoI recognition sequences. All three are *created* by E. coli K12 codon
    // optimisation (Leu CUG + Gln CAG = CTGCAG, His CAU + Met AUG = CATATG, Thr ACC + Met AUG =
    // CCATGG) and are absent from the three inputs below, which is exactly the chaining gap:
    // RemoveRestrictionSites had to be run by hand after OptimizeSequence.
    private static readonly string[] Sites = { "CTGCAG", "CATATG", "CCATGG" };

    // MKLQALHMTMAEF* with deliberately sub-optimal codons (42 nt).
    private const string GeneA = "ATGAAGCTACAAGCATTACACATGACAATGGCTGAGTTTTAA";

    // GFP N-terminus (wild-type Aequorea victoria codons) plus a LQHMTMA tail (150 nt).
    private const string GeneB =
        "ATGAGTAAAGGAGAAGAATTTTTATTCACTGGAGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCT" +
        "GTCTCTGGAGAAGGTGAAGGTGATGCAACATACGGAAAACTACAACACATGACAATGGCATAA";

    // MLQHMTMAK* — every site in a 30-nt window (30 nt).
    private const string GeneC = "ATGCTACAACACATGACAATGGCAAAGTAA";

    #region avoidSites — chained site removal

    [Test]
    [Description("A2-13: without avoidSites, MaximizeCAI CREATES the PstI/NdeI/NcoI sites (the defect)")]
    public void OptimizeSequence_WithoutAvoidSites_CanCreateRestrictionSites()
    {
        foreach (string gene in new[] { GeneA, GeneB, GeneC })
        {
            var plain = CodonOptimizer.OptimizeSequence(
                gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
            string dna = plain.OptimizedSequence.Replace('U', 'T');
            Assert.Multiple(() =>
            {
                Assert.That(dna, Does.Contain("CTGCAG"), "PstI site created by codon selection");
                Assert.That(dna, Does.Contain("CATATG"), "NdeI site created by codon selection");
                Assert.That(dna, Does.Contain("CCATGG"), "NcoI site created by codon selection");
            });
        }
    }

    [Test]
    [Description("A2-13: avoidSites clears every site on both strands and preserves the protein")]
    public void OptimizeSequence_AvoidSites_RemovesAllSites_AndPreservesProtein()
    {
        var cases = new[]
        {
            (Gene: GeneA, Protein: "MKLQALHMTMAEF*"),
            (Gene: GeneB, Protein: "MSKGEEFLFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLQHMTMA*"),
            (Gene: GeneC, Protein: "MLQHMTMAK*"),
        };

        foreach (var (gene, protein) in cases)
        {
            foreach (var strategy in new[]
                     {
                         CodonOptimizer.OptimizationStrategy.MaximizeCAI,
                         CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                         CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
                         CodonOptimizer.OptimizationStrategy.HarmonizeExpression,
                     })
            {
                var result = CodonOptimizer.OptimizeSequence(
                    gene, CodonOptimizer.EColiK12, strategy, avoidSites: Sites);
                string dna = result.OptimizedSequence.Replace('U', 'T');
                string reverseComplement = DnaSequence.GetReverseComplementString(dna);

                Assert.Multiple(() =>
                {
                    foreach (string site in Sites)
                    {
                        Assert.That(dna, Does.Not.Contain(site), $"[{strategy}] forward strand keeps {site}");
                        Assert.That(reverseComplement, Does.Not.Contain(site), $"[{strategy}] reverse strand keeps {site}");
                    }

                    Assert.That(result.ProteinSequence, Is.EqualTo(protein), $"[{strategy}] protein");
                    Assert.That(
                        Translator.Translate(result.OptimizedSequence).ToString(), Is.EqualTo(protein),
                        $"[{strategy}] re-translated output");
                    Assert.That(result.OptimizedSequence, Has.Length.EqualTo(result.OriginalSequence.Length));
                });
            }
        }
    }

    [Test]
    [Description("A2-13: locked outputs — identical to dnachisel's own AvoidPattern + CodonOptimize solution")]
    public void OptimizeSequence_AvoidSites_MatchesDnaChiselSolution()
    {
        // dnachisel 3.2.16 DnaOptimizationProblem(..., resolve_constraints() + optimize()) returns
        // exactly these sequences for the same three AvoidPattern constraints, EnforceTranslation
        // and CodonOptimize(species='e_coli', method='use_best_codon'); objective score −1.3892
        // for all three (review 2026-09 F54).
        var expected = new[]
        {
            (Gene: GeneA, Optimized: "AUGAAACUGCAAGCGCUGCACAUGACGAUGGCGGAAUUUUAA", Cai: 0.870294, Changed: 7),
            (Gene: GeneB, Optimized:
                "AUGAGCAAAGGCGAAGAAUUUCUGUUUACCGGCGUGGUGCCGAUUCUGGUGGAACUGGAUGGCGAUGUGAACGGCCAUAAAUUUAGC" +
                "GUGAGCGGCGAAGGCGAAGGCGAUGCGACCUAUGGCAAACUGCAACACAUGACGAUGGCGUAA", Cai: 0.970251, Changed: 30),
            (Gene: GeneC, Optimized: "AUGCUGCAACACAUGACGAUGGCGAAAUAA", Cai: 0.793311, Changed: 4),
        };

        foreach (var (gene, optimized, cai, changed) in expected)
        {
            var result = CodonOptimizer.OptimizeSequence(
                gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
                avoidSites: Sites);

            Assert.Multiple(() =>
            {
                Assert.That(result.OptimizedSequence, Is.EqualTo(optimized));
                Assert.That(result.OptimizedCAI, Is.EqualTo(cai).Within(1e-6), "CAI of the FINAL sequence");
                Assert.That(result.ChangedCodons, Is.EqualTo(changed), "changes describe the final sequence");
                Assert.That(result.Changes, Has.Count.EqualTo(changed));
                Assert.That(
                    result.GcContentOptimized,
                    Is.EqualTo(GcFraction(result.OptimizedSequence)).Within(1e-12),
                    "GC fraction describes the final sequence");
            });
        }

        static double GcFraction(string rna)
        {
            int gc = rna.Count(c => c is 'G' or 'C');
            return rna.Length == 0 ? 0 : (double)gc / rna.Length;
        }
    }

    [Test]
    [Description("A2-13: the site pass costs CAI relative to the unconstrained optimum (hard constraint wins)")]
    public void OptimizeSequence_AvoidSites_ReportsMetricsOfFinalSequence()
    {
        var plain = CodonOptimizer.OptimizeSequence(
            GeneC, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var avoided = CodonOptimizer.OptimizeSequence(
            GeneC, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            avoidSites: Sites);

        Assert.Multiple(() =>
        {
            Assert.That(plain.OptimizedCAI, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(avoided.OptimizedCAI, Is.LessThan(plain.OptimizedCAI));
            Assert.That(avoided.OptimizedCAI, Is.EqualTo(0.793311).Within(1e-6));
            Assert.That(
                avoided.OptimizedCAI,
                Is.EqualTo(CodonOptimizer.CalculateCAI(avoided.OptimizedSequence, CodonOptimizer.EColiK12)).Within(1e-12),
                "the reported CAI is the CAI of the returned sequence");
        });
    }

    [Test]
    [Description("A2-13: null / empty avoidSites is the historical behaviour; the pass is idempotent")]
    public void OptimizeSequence_AvoidSites_DefaultsAreUnchanged()
    {
        var dflt = CodonOptimizer.OptimizeSequence(GeneA, CodonOptimizer.EColiK12);
        var explicitNull = CodonOptimizer.OptimizeSequence(GeneA, CodonOptimizer.EColiK12, avoidSites: null);
        var empty = CodonOptimizer.OptimizeSequence(GeneA, CodonOptimizer.EColiK12, avoidSites: Array.Empty<string>());
        var avoided = CodonOptimizer.OptimizeSequence(GeneA, CodonOptimizer.EColiK12, avoidSites: Sites);
        var twice = CodonOptimizer.OptimizeSequence(
            avoided.OptimizedSequence, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization, avoidSites: Sites);

        Assert.Multiple(() =>
        {
            Assert.That(explicitNull.OptimizedSequence, Is.EqualTo(dflt.OptimizedSequence));
            Assert.That(empty.OptimizedSequence, Is.EqualTo(dflt.OptimizedSequence));
            Assert.That(twice.OptimizedSequence, Is.EqualTo(avoided.OptimizedSequence), "idempotent");
        });
    }

    [Test]
    [Description("A2-13: avoidSites runs with the caller's genetic code and rejects non-IUPAC sites")]
    public void OptimizeSequence_AvoidSites_UsesGeneticCode_AndValidatesSites()
    {
        // Table 2: AGA/AGG are stops, so the Arg family is CGU/CGC/CGA/CGG only — the site pass
        // must pick replacements from that family (protein preserved under table 2).
        const string argGene = "AUGCGUCGUCGUCGUUAA";
        var table2 = CodonOptimizer.OptimizeSequence(
            argGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            geneticCode: GeneticCode.GetByTableNumber(2), avoidSites: new[] { "CGCCGC" });

        Assert.Multiple(() =>
        {
            Assert.That(table2.OptimizedSequence.Replace('U', 'T'), Does.Not.Contain("CGCCGC"));
            Assert.That(
                Translator.Translate(table2.OptimizedSequence, GeneticCode.GetByTableNumber(2)),
                Is.EqualTo(Translator.Translate(argGene, GeneticCode.GetByTableNumber(2))));
            Assert.Throws<ArgumentException>(() => CodonOptimizer.OptimizeSequence(
                GeneA, CodonOptimizer.EColiK12, avoidSites: new[] { "GAATTZ" }));
        });
    }

    #endregion

    #region AvoidRareCodons spelling (A2-14)

    [Test]
    [Description("A2-14: AvoidRareCodons and the historical AvoidRareCodeons are the same member")]
    public void OptimizationStrategy_AvoidRareCodons_IsAliasOfAvoidRareCodeons()
    {
        Assert.Multiple(() =>
        {
            Assert.That(
                CodonOptimizer.OptimizationStrategy.AvoidRareCodons,
                Is.EqualTo(CodonOptimizer.OptimizationStrategy.AvoidRareCodeons));
            Assert.That(
                (int)CodonOptimizer.OptimizationStrategy.AvoidRareCodons,
                Is.EqualTo((int)CodonOptimizer.OptimizationStrategy.AvoidRareCodeons));

            // Both spellings parse, case-insensitively, to that one value (the MCP layer takes the
            // strategy as a JSON string: System.Text.Json's enum converter matches either name).
            foreach (string name in new[]
                     {
                         "AvoidRareCodons", "avoidrarecodons", "AVOIDRARECODONS",
                         "AvoidRareCodeons", "avoidrarecodeons",
                     })
            {
                Assert.That(
                    Enum.Parse<CodonOptimizer.OptimizationStrategy>(name, ignoreCase: true),
                    Is.EqualTo(CodonOptimizer.OptimizationStrategy.AvoidRareCodeons), name);
                Assert.That(
                    System.Text.Json.JsonSerializer.Deserialize<CodonOptimizer.OptimizationStrategy>(
                        $"\"{name}\"",
                        new System.Text.Json.JsonSerializerOptions
                        {
                            Converters = { new System.Text.Json.Serialization.JsonStringEnumConverter() },
                        }),
                    Is.EqualTo(CodonOptimizer.OptimizationStrategy.AvoidRareCodeons), $"JSON {name}");
            }
        });
    }

    [Test]
    [Description("A2-14: both spellings produce identical optimisation output")]
    public void OptimizationStrategy_BothSpellings_BehaveIdentically()
    {
        const string mixed = "AUGCUGCUAAGAGAGUUUUAA"; // CUG common, CUA/AGA/GAG rarer in E. coli

        var legacy = CodonOptimizer.OptimizeSequence(
            mixed, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons);
        var fixedSpelling = CodonOptimizer.OptimizeSequence(
            mixed, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.AvoidRareCodons);

        Assert.Multiple(() =>
        {
            Assert.That(fixedSpelling.OptimizedSequence, Is.EqualTo(legacy.OptimizedSequence));
            Assert.That(fixedSpelling.OptimizedCAI, Is.EqualTo(legacy.OptimizedCAI));
            Assert.That(fixedSpelling.ChangedCodons, Is.EqualTo(legacy.ChangedCodons));
        });
    }

    [Test]
    [Description("A2-14: Enum.GetValues consumers still work — the alias shares a value, it adds no strategy")]
    public void OptimizationStrategy_GetValues_ExposesTheAliasWithoutANewStrategy()
    {
        var values = Enum.GetValues<CodonOptimizer.OptimizationStrategy>();
        var distinct = values.Distinct().ToArray();

        Assert.Multiple(() =>
        {
            // Both fields are enumerated (duplicate values), but there are still five strategies,
            // so a consumer that iterates Enum.GetValues runs the rare-codon strategy twice and
            // never sees an unknown value.
            Assert.That(values, Has.Length.EqualTo(6));
            Assert.That(distinct, Has.Length.EqualTo(5));
            Assert.That(
                distinct,
                Is.EquivalentTo(new[]
                {
                    CodonOptimizer.OptimizationStrategy.MaximizeCAI,
                    CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                    CodonOptimizer.OptimizationStrategy.HarmonizeExpression,
                    CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
                    CodonOptimizer.OptimizationStrategy.AvoidRareCodeons,
                }));
            // ToString keeps the historical name, so no serialized value changes meaning.
            Assert.That(CodonOptimizer.OptimizationStrategy.AvoidRareCodons.ToString(), Is.EqualTo("AvoidRareCodeons"));
        });
    }

    #endregion
}
