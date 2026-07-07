// ONCO-SIG-004 — Mutational Process Classification (SBS exposure → active processes)
// Evidence: docs/Evidence/ONCO-SIG-004-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-004.md
// Source: Rosenthal R. et al. (2016). deconstructSigs. Genome Biology 17:31. https://doi.org/10.1186/s13059-016-0893-4
//         deconstructSigs whichSignatures.R (signature.cutoff = 0.06). https://github.com/raerose01/deconstructSigs/blob/master/R/whichSignatures.R
//         COSMIC SBS proposed aetiologies. https://cancer.sanger.ac.uk/signatures/sbs/
//         Alexandrov L.B. et al. (2020). Nature 578:94-101. https://doi.org/10.1038/s41586-020-1943-3

using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ClassifyMutationalProcess_Tests
{
    // Hand-derived dataset (Evidence §Test Datasets): raw exposures sum to 100, so normalized = raw/100.
    // SBS2 0.50, SBS13 0.30 -> APOBEC 0.80; SBS1 0.15 -> Aging; SBS4 0.05 < 0.06 -> dropped.
    private static IReadOnlyList<(string, double)> CanonicalDataset() => new (string, double)[]
    {
        ("SBS2", 50), ("SBS13", 30), ("SBS1", 15), ("SBS4", 5),
    };

    private static double ContributionOf(MutationalProcessClassification c, MutationalProcess p) =>
        c.ActiveProcesses.Single(a => a.Process == p).Contribution;

    #region ClassifyMutationalProcess — contributions, cutoff, aggregation, dominance

    // M1 — normalized contribution = exposure / Σ exposure (deconstructSigs: weights normalized 0..1).
    [Test]
    public void ClassifyMutationalProcess_NormalizesContributionsByTotal()
    {
        var result = ClassifyMutationalProcess(CanonicalDataset());

        // APOBEC = SBS2(0.50)+SBS13(0.30) = 0.80; Aging = SBS1(0.15). These are exact raw/100 sums.
        Assert.Multiple(() =>
        {
            Assert.That(ContributionOf(result, MutationalProcess.Apobec), Is.EqualTo(0.80).Within(1e-10),
                "APOBEC = (50+30)/100 = 0.80 per normalized additive contributions (Rosenthal 2016).");
            Assert.That(ContributionOf(result, MutationalProcess.Aging), Is.EqualTo(0.15).Within(1e-10),
                "Aging = SBS1 15/100 = 0.15.");
        });
    }

    // M2 — a signature whose normalized contribution < 0.06 is excluded (deconstructSigs 6% cutoff).
    [Test]
    public void ClassifyMutationalProcess_SubCutoffSignature_IsExcluded()
    {
        var result = ClassifyMutationalProcess(CanonicalDataset());

        // SBS4 = 0.05 < 0.06 -> Tobacco smoking must NOT appear as an active process.
        Assert.That(result.ActiveProcesses.Any(a => a.Process == MutationalProcess.TobaccoSmoking), Is.False,
            "SBS4 at 0.05 is below the 0.06 cutoff and must be dropped (weights[weights<0.06]<-0).");
    }

    // M3 — boundary: a contribution of exactly 0.06 is RETAINED (strict less-than cutoff).
    [Test]
    public void ClassifyMutationalProcess_ContributionExactlyAtCutoff_IsRetained()
    {
        // SBS4 = 6, others = 94: SBS4 normalized = 6/100 = 0.06 exactly; SBS1 = 0.94.
        var exposures = new (string, double)[] { ("SBS1", 94), ("SBS4", 6) };

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses.Any(a => a.Process == MutationalProcess.TobaccoSmoking), Is.True,
                "0.06 is NOT < 0.06, so SBS4 is retained (strict-less-than cutoff, deconstructSigs).");
            Assert.That(ContributionOf(result, MutationalProcess.TobaccoSmoking), Is.EqualTo(0.06).Within(1e-10),
                "Tobacco contribution = 6/100 = 0.06.");
        });
    }

    // M4 — just below the boundary (< 0.06) is excluded.
    [Test]
    public void ClassifyMutationalProcess_ContributionJustBelowCutoff_IsExcluded()
    {
        // SBS4 = 5.999, total = 100 -> 0.05999 < 0.06 -> excluded.
        var exposures = new (string, double)[] { ("SBS1", 94.001), ("SBS4", 5.999) };

        var result = ClassifyMutationalProcess(exposures);

        Assert.That(result.ActiveProcesses.Any(a => a.Process == MutationalProcess.TobaccoSmoking), Is.False,
            "0.05999 < 0.06 so SBS4 is excluded.");
    }

    // M5 — APOBEC aggregation: SBS2 + SBS13 sum into one APOBEC contribution.
    [Test]
    public void ClassifyMutationalProcess_AggregatesApobecSignatures()
    {
        var exposures = new (string, double)[] { ("SBS2", 50), ("SBS13", 30) }; // total 80

        var result = ClassifyMutationalProcess(exposures);

        Assert.That(ContributionOf(result, MutationalProcess.Apobec), Is.EqualTo(1.0).Within(1e-10),
            "Only APOBEC signatures present -> APOBEC = (50+30)/80 = 1.0 (additive weights).");
    }

    // M6 — active set is exactly {APOBEC, Aging} for the canonical dataset.
    [Test]
    public void ClassifyMutationalProcess_ActiveSet_ExcludesSubCutoffProcess()
    {
        var result = ClassifyMutationalProcess(CanonicalDataset());

        var active = result.ActiveProcesses.Select(a => a.Process).ToArray();
        Assert.That(active, Is.EquivalentTo(new[] { MutationalProcess.Apobec, MutationalProcess.Aging }),
            "Active processes are APOBEC and Aging; Tobacco (SBS4 sub-cutoff) is excluded.");
    }

    // M7 — dominant process is the one with the largest aggregated contribution (APOBEC 0.80 > Aging 0.15).
    [Test]
    public void ClassifyMutationalProcess_DominantProcess_IsLargestContribution()
    {
        var result = ClassifyMutationalProcess(CanonicalDataset());

        Assert.Multiple(() =>
        {
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Apobec),
                "APOBEC (0.80) > Aging (0.15) -> dominant = APOBEC.");
            Assert.That(result.ActiveProcesses[0].Process, Is.EqualTo(MutationalProcess.Apobec),
                "Active processes are ordered by descending contribution; APOBEC first.");
        });
    }

    // M8 — MMR aggregation: SBS6/SBS15/SBS20/SBS26 all collapse into one MMR-deficiency process.
    [Test]
    public void ClassifyMutationalProcess_AggregatesMismatchRepairSignatures()
    {
        var exposures = new (string, double)[]
        {
            ("SBS6", 25), ("SBS15", 25), ("SBS20", 25), ("SBS26", 25), // total 100
        };

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses, Has.Count.EqualTo(1),
                "All four signatures map to MMR deficiency -> a single active process.");
            Assert.That(ContributionOf(result, MutationalProcess.MismatchRepairDeficiency),
                Is.EqualTo(1.0).Within(1e-10),
                "MMR = (25+25+25+25)/100 = 1.0 (four MMR signatures summed).");
        });
    }

    // C1 — UV multi-subtype: SBS7a + SBS7b both map to UV and sum.
    [Test]
    public void ClassifyMutationalProcess_AggregatesUvSubtypes()
    {
        var exposures = new (string, double)[] { ("SBS7a", 40), ("SBS7b", 60) }; // total 100

        var result = ClassifyMutationalProcess(exposures);

        Assert.That(ContributionOf(result, MutationalProcess.UltravioletLight), Is.EqualTo(1.0).Within(1e-10),
            "SBS7a (0.40) + SBS7b (0.60) -> UV = 1.0 (all SBS7 subtypes are UV).");
    }

    // S3 — single dominant signature: one signature at 100% -> one active process, it is dominant.
    [Test]
    public void ClassifyMutationalProcess_SingleSignature_IsSoleDominant()
    {
        var exposures = new (string, double)[] { ("SBS4", 200) };

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.TobaccoSmoking),
                "Only Tobacco present -> dominant = Tobacco smoking.");
            Assert.That(ContributionOf(result, MutationalProcess.TobaccoSmoking), Is.EqualTo(1.0).Within(1e-10),
                "Single signature normalized to 1.0.");
        });
    }

    // S4 — custom cutoff override: raising the cutoff to 0.20 drops the 0.15 Aging signature.
    [Test]
    public void ClassifyMutationalProcess_CustomCutoff_DropsSignatureBelowOverride()
    {
        var result = ClassifyMutationalProcess(CanonicalDataset(), contributionCutoff: 0.20);

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses.Any(a => a.Process == MutationalProcess.Aging), Is.False,
                "Aging (SBS1 = 0.15) is below the 0.20 override cutoff and is dropped.");
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Apobec),
                "APOBEC (0.80 >= 0.20) survives and remains dominant.");
        });
    }

    #endregion

    #region Edge cases and invalid input

    // S1 — all-zero exposures: Σ = 0 -> no active processes, dominant = Unknown (INV-05).
    [Test]
    public void ClassifyMutationalProcess_AllZeroExposures_ReturnsEmpty()
    {
        var exposures = new (string, double)[] { ("SBS1", 0), ("SBS2", 0) };

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses, Is.Empty, "Σ exposure = 0 -> no active processes.");
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Unknown),
                "No process active -> dominant = Unknown.");
        });
    }

    // Empty list -> empty result.
    [Test]
    public void ClassifyMutationalProcess_EmptyList_ReturnsEmpty()
    {
        var result = ClassifyMutationalProcess(Array.Empty<(string, double)>());

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses, Is.Empty, "No exposures -> no active processes.");
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Unknown), "No dominant process.");
        });
    }

    // S2 — unmapped signature label contributes to no recognized process.
    [Test]
    public void ClassifyMutationalProcess_UnmappedSignature_ContributesToNoProcess()
    {
        // SBS99 is not in the COSMIC aetiology map; SBS2 is APOBEC.
        var exposures = new (string, double)[] { ("SBS2", 50), ("SBS99", 50) };

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(result.ActiveProcesses, Has.Count.EqualTo(1),
                "Only the mapped SBS2 forms a process; SBS99 (unmapped) contributes to none.");
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Apobec),
                "APOBEC is the only recognized active process.");
        });
    }

    [Test]
    public void ClassifyMutationalProcess_NullExposures_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ClassifyMutationalProcess(null!),
            "Null exposures list is invalid.");
    }

    [Test]
    public void ClassifyMutationalProcess_NegativeExposure_Throws()
    {
        var exposures = new (string, double)[] { ("SBS1", -1.0) };

        Assert.Throws<ArgumentException>(
            () => ClassifyMutationalProcess(exposures),
            "Negative exposures are invalid.");
    }

    // NaN exposure is invalid input (contract §3.3: "a negative or NaN exposure → ArgumentException").
    [Test]
    public void ClassifyMutationalProcess_NaNExposure_Throws()
    {
        var exposures = new (string, double)[] { ("SBS1", double.NaN) };

        Assert.Throws<ArgumentException>(
            () => ClassifyMutationalProcess(exposures),
            "A NaN exposure is invalid.");
    }

    // NaN cutoff is rejected before the [0,1) range check (contract §3.3).
    [Test]
    public void ClassifyMutationalProcess_NaNCutoff_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ClassifyMutationalProcess(CanonicalDataset(), contributionCutoff: double.NaN),
            "A NaN cutoff is outside the valid [0,1) interval.");
    }

    // Tie-break ordering: two processes with equal aggregated contribution are ordered by the
    // MutationalProcess enum value (Aging=1 before Apobec=2), so dominant is the lower enum.
    // Aging (SBS1=0.50) and APOBEC (SBS2=0.50) tie at 0.50 each.
    [Test]
    public void ClassifyMutationalProcess_EqualContributions_OrderedByProcessEnum()
    {
        var exposures = new (string, double)[] { ("SBS1", 50), ("SBS2", 50) }; // total 100

        var result = ClassifyMutationalProcess(exposures);

        Assert.Multiple(() =>
        {
            Assert.That(ContributionOf(result, MutationalProcess.Aging), Is.EqualTo(0.50).Within(1e-10),
                "Aging = SBS1 50/100 = 0.50.");
            Assert.That(ContributionOf(result, MutationalProcess.Apobec), Is.EqualTo(0.50).Within(1e-10),
                "APOBEC = SBS2 50/100 = 0.50.");
            Assert.That(result.ActiveProcesses[0].Process, Is.EqualTo(MutationalProcess.Aging),
                "On equal contribution the lower MutationalProcess enum (Aging) sorts first.");
            Assert.That(result.DominantProcess, Is.EqualTo(MutationalProcess.Aging),
                "Tie is broken deterministically by process enum: Aging precedes APOBEC.");
        });
    }

    [Test]
    public void ClassifyMutationalProcess_CutoffOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ClassifyMutationalProcess(CanonicalDataset(), contributionCutoff: 1.0),
                "Cutoff must be < 1.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ClassifyMutationalProcess(CanonicalDataset(), contributionCutoff: -0.1),
                "Cutoff must be >= 0.");
        });
    }

    #endregion

    #region GetMutationalProcess — COSMIC SBS aetiology lookup

    // M9 — each COSMIC SBS label maps to its proposed-aetiology process.
    [TestCase("SBS1", MutationalProcess.Aging)]
    [TestCase("SBS5", MutationalProcess.Aging)]
    [TestCase("SBS2", MutationalProcess.Apobec)]
    [TestCase("SBS13", MutationalProcess.Apobec)]
    [TestCase("SBS4", MutationalProcess.TobaccoSmoking)]
    [TestCase("SBS7a", MutationalProcess.UltravioletLight)]
    [TestCase("SBS6", MutationalProcess.MismatchRepairDeficiency)]
    [TestCase("SBS15", MutationalProcess.MismatchRepairDeficiency)]
    [TestCase("SBS20", MutationalProcess.MismatchRepairDeficiency)]
    [TestCase("SBS26", MutationalProcess.MismatchRepairDeficiency)]
    public void GetMutationalProcess_MapsLabelToCosmicAetiology(string label, MutationalProcess expected)
    {
        Assert.That(GetMutationalProcess(label), Is.EqualTo(expected),
            $"{label} must map to {expected} per COSMIC proposed aetiology.");
    }

    [Test]
    public void GetMutationalProcess_IsCaseInsensitive()
    {
        Assert.That(GetMutationalProcess("sbs2"), Is.EqualTo(MutationalProcess.Apobec),
            "Label matching is case-insensitive.");
    }

    [Test]
    public void GetMutationalProcess_UnmappedLabel_ReturnsUnknown()
    {
        Assert.That(GetMutationalProcess("SBS99"), Is.EqualTo(MutationalProcess.Unknown),
            "A label outside the COSMIC map resolves to Unknown.");
    }

    [Test]
    public void GetMutationalProcess_NullLabel_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => GetMutationalProcess(null!),
            "Null label is invalid.");
    }

    #endregion
}
