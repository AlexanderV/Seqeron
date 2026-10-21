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

    [Test]
    [Description("B08 F71: the pass runs on RestrictionAnalyzer.FindPatternOccurrences but an ambiguous triplet still matches no pattern code (an N in the site is not satisfied by a sequence N)")]
    public void RemoveRestrictionSites_AmbiguousTriplet_NeverCompletesAPatternNWildcard()
    {
        // Site "TCCNGG": the unambiguous twin ATG TCC AGG TTT is an occurrence (TCCAGG) and is removed
        // by a synonymous Ser/Arg swap; with the sequence N in the wildcard slot (ATG TCC NGG TTT)
        // the optimizer, which counts ambiguous triplets as contributing nothing, leaves it alone.
        string clean = CodonOptimizer.RemoveRestrictionSites(
            "AUGUCCAGGUUU", new[] { "TCCNGG" }, CodonOptimizer.EColiK12);
        string ambiguous = CodonOptimizer.RemoveRestrictionSites(
            "AUGUCCNGGUUU", new[] { "TCCNGG" }, CodonOptimizer.EColiK12);

        Assert.Multiple(() =>
        {
            Assert.That(clean.Replace('U', 'T'), Does.Not.Contain("TCCAGG"), "positive control");
            Assert.That(ambiguous, Is.EqualTo("AUGUCCNGGUUU"));
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

    #region GC-constraint-aware site removal (A3-6 / F60)

    // Three constructed genes whose BalancedOptimization output carries exactly one of the three
    // sites, and whose *best-frequency* site-removing synonym flips the GC fraction of a 24-nt
    // window across the 0.60 bound — the chaining defect: the GC pass satisfied
    // EnforceGCContent(0.40, 0.60, window=24), the old site pass broke it again.
    // dnachisel 3.2.16 numbers (docs/Evidence/CODON-OPT-001-Evidence.md, review 2026-09 F60):
    // for each gene the GC-aware output is byte-identical to the sequence dnachisel's own
    // resolve_constraints() + optimize() returns from the same GC-balanced start, and all of
    // AvoidPattern x3 / EnforceGCContent(0.40, 0.60, window=24) / EnforceTranslation pass.
    private const double GcMin = 0.40;
    private const double GcMax = 0.60;
    private const int GcWindow = 24;

    private static readonly (string Gene, string Balanced, string OldChain, double OldBreach, string Expected, string Protein)[] GcCases =
    {
        ("ATGATTAGACTTTGCACGCCAGCCTATGGACGAAGACGCTAA",
         "AUGAUUCGCCUGUGCACUCCAGCAUAUGGUCGUCGUCGCUAA",
         "AUGAUUCGCCUGUGCACUCCAGCAUACGGUCGUCGUCGCUAA", 0.225,
         "AUGAUUCGCCUGUGCACUCCAGCUUAUGGUCGUCGUCGCUAA", "MIRLCTPAYGRRR*"),
        ("ATGCGTATGGCCCTCCGTTCGGCTGTGGGCTGGCTAGTGTAA",
         "AUGCGUAUGGCACUGCGUUCUGCAGUUGGUUGGCUGGUGUAA",
         "AUGCGUAUGGCACUGCGUUCUGCGGUUGGUUGGCUGGUGUAA", 0.225,
         "AUGCGUAUGGCACUGCGUUCUGCUGUUGGUUGGCUGGUGUAA", "MRMALRSAVGWLV*"),
        ("ATGGGTTCCGAAACGTCTGCTTACGGTGCAGCGTTGTGGTAA",
         "AUGGGCAGCGAAACCUCUGCAUAUGGUGCAGCGCUGUGGUAA",
         "AUGGGCAGCGAAACCUCUGCAUACGGUGCAGCGCUGUGGUAA", 0.150,
         "AUGGGCAGCGAAACCUCUGCUUAUGGUGCAGCGCUGUGGUAA", "MGSETSAYGAALW*"),
    };

    [Test]
    [Description("A3-6: the old chain (GC pass, then GC-blind site removal) breaks EnforceGCContent — the defect")]
    public void RemoveRestrictionSites_WithoutGcConstraint_CanBreakTheGcWindows()
    {
        foreach (var (gene, balanced, oldChain, oldBreach, _, _) in GcCases)
        {
            var gcOnly = CodonOptimizer.OptimizeSequence(
                gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                GcMin, GcMax, gcWindow: GcWindow);

            // The public 4-argument RemoveRestrictionSites overload is unchanged: it stays GC-blind.
            string chained = CodonOptimizer.RemoveRestrictionSites(
                gcOnly.OptimizedSequence, Sites, CodonOptimizer.EColiK12);

            Assert.Multiple(() =>
            {
                Assert.That(gcOnly.OptimizedSequence, Is.EqualTo(balanced), "GC-balanced start");
                Assert.That(
                    CodonOptimizer.EvaluateGcContent(balanced, GcMin, GcMax, GcWindow).Score,
                    Is.EqualTo(0).Within(1e-12), "the GC pass satisfied the constraint");
                Assert.That(chained, Is.EqualTo(oldChain), "GC-blind site removal (historical behaviour)");
                Assert.That(
                    -CodonOptimizer.EvaluateGcContent(chained, GcMin, GcMax, GcWindow).Score,
                    Is.EqualTo(oldBreach).Within(1e-6),
                    "dnachisel EnforceGCContent(0.40, 0.60, window=24) breach of the old chain");
            });
        }
    }

    [Test]
    [Description("A3-6: the avoidSites pass keeps the GC windows satisfied — dnachisel's own solution")]
    public void OptimizeSequence_AvoidSites_KeepsTheWindowedGcConstraint()
    {
        foreach (var (gene, _, _, _, expected, protein) in GcCases)
        {
            var result = CodonOptimizer.OptimizeSequence(
                gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                GcMin, GcMax, gcWindow: GcWindow, avoidSites: Sites);

            string dna = result.OptimizedSequence.Replace('U', 'T');
            string reverseComplement = DnaSequence.GetReverseComplementString(dna);

            Assert.Multiple(() =>
            {
                Assert.That(result.OptimizedSequence, Is.EqualTo(expected), "dnachisel-identical output");
                Assert.That(
                    CodonOptimizer.EvaluateGcContent(result.OptimizedSequence, GcMin, GcMax, GcWindow).Score,
                    Is.EqualTo(0).Within(1e-12), "every 24-nt window inside [0.40, 0.60]");
                foreach (string site in Sites)
                {
                    Assert.That(dna, Does.Not.Contain(site), $"forward strand keeps {site}");
                    Assert.That(reverseComplement, Does.Not.Contain(site), $"reverse strand keeps {site}");
                }

                Assert.That(result.ProteinSequence, Is.EqualTo(protein), "protein");
                Assert.That(
                    Translator.Translate(result.OptimizedSequence).ToString(), Is.EqualTo(protein),
                    "re-translated output");
            });
        }
    }

    [Test]
    [Description("A3-6: the global gcTargetMin/Max bounds are kept too (window = null)")]
    public void OptimizeSequence_AvoidSites_KeepsTheGlobalGcBounds()
    {
        // dnachisel: EnforceGCContent(mini=0.40, maxi=0.60) scores −0.019048 on the old chain
        // (GC 0.619048) and 0 on the GC-aware output, which is byte-identical to dnachisel's own
        // resolve_constraints() + optimize() result.
        const string gene = "ATGGTAACTGCGGTACCCATGAGGGGGATGGGCAGACGCTAA";
        const string balanced = "AUGGUUACUGCAGUUCCAAUGCGCGGCAUGGGCCGCCGCUAA";
        const string oldChain = "AUGGUUACCGCAGUUCCAAUGCGCGGCAUGGGCCGCCGCUAA";
        const string expected = "AUGGUUACUGCUGUUCCAAUGCGCGGCAUGGGCCGCCGCUAA";

        var gcOnly = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            GcMin, GcMax);
        string chained = CodonOptimizer.RemoveRestrictionSites(
            gcOnly.OptimizedSequence, Sites, CodonOptimizer.EColiK12);
        var result = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            GcMin, GcMax, avoidSites: Sites);

        Assert.Multiple(() =>
        {
            Assert.That(gcOnly.OptimizedSequence, Is.EqualTo(balanced));
            Assert.That(chained, Is.EqualTo(oldChain), "GC-blind chain");
            Assert.That(
                -CodonOptimizer.EvaluateGcContent(chained, GcMin, GcMax).Score,
                Is.EqualTo(0.019048).Within(1e-6), "dnachisel breach of the old chain");
            Assert.That(result.OptimizedSequence, Is.EqualTo(expected));
            Assert.That(result.GcContentOptimized, Is.InRange(GcMin, GcMax));
            Assert.That(
                CodonOptimizer.EvaluateGcContent(result.OptimizedSequence, GcMin, GcMax).Score,
                Is.EqualTo(0).Within(1e-12));
            Assert.That(result.OptimizedSequence.Replace('U', 'T'), Does.Not.Contain("CTGCAG"));
            Assert.That(Translator.Translate(result.OptimizedSequence).ToString(), Is.EqualTo("MVTAVPMRGMGRR*"));
        });
    }

    [Test]
    [Description("A3-6: unreachable GC bounds — the site is still removed with the least-breach candidate")]
    public void OptimizeSequence_AvoidSites_UnreachableGcBounds_RemovesTheSiteWithLeastBreach()
    {
        // Bounds [0.30, 0.45] over 24-nt windows are unreachable for this protein: dnachisel's
        // resolve_constraints() raises NoSolutionError for the same problem. Our documented
        // deterministic fallback removes the site anyway (AvoidPattern stays the hard constraint)
        // with the candidate that does not increase the pre-existing breach: dnachisel scores the
        // GC-balanced start −0.483333, the GC-blind chain −1.066667 and our output −0.483333.
        const string gene = "ATGGGTTCCGAAACGTCTGCTTACGGTGCAGCGTTGTGGTAA";
        const double tightMin = 0.30;
        const double tightMax = 0.45;

        var gcOnly = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            tightMin, tightMax, gcWindow: GcWindow);
        string chained = CodonOptimizer.RemoveRestrictionSites(
            gcOnly.OptimizedSequence, Sites, CodonOptimizer.EColiK12);
        var result = CodonOptimizer.OptimizeSequence(
            gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            tightMin, tightMax, gcWindow: GcWindow, avoidSites: Sites);

        double Breach(string sequence) =>
            -CodonOptimizer.EvaluateGcContent(sequence, tightMin, tightMax, GcWindow).Score;

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGGGUUCUGAAACUUCUGCUUAUGGUGCAGCACUGUGGUAA"));
            Assert.That(Breach(gcOnly.OptimizedSequence), Is.EqualTo(0.483333).Within(1e-6), "start already breaches");
            Assert.That(Breach(chained), Is.EqualTo(1.066667).Within(1e-6), "GC-blind chain makes it worse");
            Assert.That(Breach(result.OptimizedSequence), Is.EqualTo(0.483333).Within(1e-6), "GC-aware pass does not");
            Assert.That(Breach(result.OptimizedSequence), Is.LessThan(Breach(chained)));
            // The site is gone either way: AvoidPattern is the hard constraint, nothing throws.
            Assert.That(result.OptimizedSequence.Replace('U', 'T'), Does.Not.Contain("CATATG"));
            Assert.That(Translator.Translate(result.OptimizedSequence).ToString(), Is.EqualTo("MGSETSAYGAALW*"));
        });
    }

    [Test]
    [Description("A3-6: strategies without a GC target keep the historical GC-blind site pass")]
    public void OptimizeSequence_AvoidSites_StrategiesWithoutGcTarget_AreUnchanged()
    {
        // MaximizeCAI / HarmonizeExpression declare no GC target, so no GC constraint is threaded
        // and the F54 locked outputs stay byte-identical.
        foreach (var (gene, expected) in new[]
                 {
                     (GeneA, "AUGAAACUGCAAGCGCUGCACAUGACGAUGGCGGAAUUUUAA"),
                     (GeneC, "AUGCUGCAACACAUGACGAUGGCGAAAUAA"),
                 })
        {
            var result = CodonOptimizer.OptimizeSequence(
                gene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
                gcTargetMin: 0.49, gcTargetMax: 0.51, gcWindow: 24, avoidSites: Sites);
            Assert.That(result.OptimizedSequence, Is.EqualTo(expected), gene);
        }
    }

    #endregion
}
