using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Chromosome;
using Seqeron.Genomics.MolTools;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Strict-mode limitation guards (default policy = throw). Each of the 5 runtime-guarded branches must
/// throw <see cref="SeqeronLimitationException"/> in <see cref="LimitationMode.Strict"/> and return the
/// historical best-effort value in <see cref="LimitationMode.Permissive"/>. The assembly default is
/// Permissive (see <c>_LimitationPolicyTestBootstrap</c>); these tests opt back into Strict with a
/// scoped <see cref="LimitationPolicy.UseStrict"/> so they exercise the production default.
/// </summary>
[TestFixture]
public class LimitationPolicy_Strict_Tests
{
    // ── policy plumbing ────────────────────────────────────────────────────────────────────

    [Test]
    public void Enforce_Throws_UnderStrict()
    {
        using (LimitationPolicy.UseStrict())
        {
            Assert.That(LimitationPolicy.IsStrict, Is.True);
            Assert.Throws<SeqeronLimitationException>(() => LimitationPolicy.Enforce("PARSE-FASTQ-001"));
        }
    }

    [Test]
    public void ScopedOverride_NestsAndRestores()
    {
        Assume.That(LimitationPolicy.CurrentMode, Is.EqualTo(LimitationMode.Permissive),
            "assembly default is Permissive");

        using (LimitationPolicy.UseStrict())
        {
            Assert.That(LimitationPolicy.IsStrict, Is.True);
            using (LimitationPolicy.UsePermissive())
                Assert.That(LimitationPolicy.IsStrict, Is.False);
            Assert.That(LimitationPolicy.IsStrict, Is.True, "inner scope restored on dispose");
        }

        Assert.That(LimitationPolicy.CurrentMode, Is.EqualTo(LimitationMode.Permissive),
            "outer scope restored to assembly default");
    }

    [Test]
    public void Enforce_NoOp_UnderPermissive()
        => Assert.DoesNotThrow(() =>
        {
            using (LimitationPolicy.UsePermissive())
                LimitationPolicy.Enforce("PARSE-FASTQ-001");
        });

    // ── catalog integrity ──────────────────────────────────────────────────────────────────

    // Non-ideal-output branches: allowed ONLY under Permissive (throw in Strict AND Moderate).
    private static readonly string[] NonIdealIds =
        { "PARSE-FASTQ-001", "CHROM-CENT-001", "DISORDER-REGION-001", "MIRNA-TARGET-001", "MIRNA-CLEAVAGE-001" };

    // Correct-but-incomplete / narrower-contract branches: allowed under Moderate+Permissive (throw in Strict).
    private static readonly string[] IncompleteIds =
        { "ONCO-MHC-001", "ONCO-IMMUNE-001", "META-BIN-001", "PROBE-DESIGN-001" };

    private static readonly string[] GuardedIds = NonIdealIds.Concat(IncompleteIds).ToArray();

    [Test]
    public void Catalog_HasExactlyTheNineGuardedUnits()
        => Assert.That(LimitationCatalog.Entries.Keys.OrderBy(x => x),
            Is.EqualTo(GuardedIds.OrderBy(x => x)));

    [TestCaseSource(nameof(GuardedIds))]
    public void Catalog_Entry_IsFullyPopulated(string id)
    {
        var info = LimitationCatalog.Get(id);
        Assert.Multiple(() =>
        {
            Assert.That(info.Id, Is.EqualTo(id));
            Assert.That(info.Branch, Is.Not.Empty, "Branch");
            Assert.That(info.Summary, Is.Not.Empty, "Summary");
            Assert.That(info.RelatedTo, Is.Not.Empty, "RelatedTo");
            Assert.That(info.Workaround, Is.Not.Empty, "Workaround");
            Assert.That(info.ReportPath, Is.EqualTo($"docs/Validation/reports/{id}.md"));
            Assert.That(info.MinimumMode, Is.AnyOf(LimitationMode.Moderate, LimitationMode.Permissive));
        });
    }

    // ── three-tier policy semantics ──────────────────────────────────────────────────────────

    [Test]
    public void Modes_AreOrderedStrictModeratePermissive()
        => Assert.That((int)LimitationMode.Strict < (int)LimitationMode.Moderate
                    && (int)LimitationMode.Moderate < (int)LimitationMode.Permissive, Is.True);

    [TestCaseSource(nameof(NonIdealIds))]
    public void NonIdealBranch_MinPermissive_ThrowsUnderStrictAndModerate_AllowedPermissive(string id)
    {
        Assert.That(LimitationCatalog.Get(id).MinimumMode, Is.EqualTo(LimitationMode.Permissive));
        using (LimitationPolicy.UseStrict())
            Assert.Throws<SeqeronLimitationException>(() => LimitationPolicy.Enforce(id));
        using (LimitationPolicy.Use(LimitationMode.Moderate))
            Assert.Throws<SeqeronLimitationException>(() => LimitationPolicy.Enforce(id));
        using (LimitationPolicy.UsePermissive())
            Assert.DoesNotThrow(() => LimitationPolicy.Enforce(id));
    }

    [TestCaseSource(nameof(IncompleteIds))]
    public void IncompleteBranch_MinModerate_ThrowsUnderStrictOnly(string id)
    {
        Assert.That(LimitationCatalog.Get(id).MinimumMode, Is.EqualTo(LimitationMode.Moderate));
        using (LimitationPolicy.UseStrict())
            Assert.Throws<SeqeronLimitationException>(() => LimitationPolicy.Enforce(id));
        using (LimitationPolicy.Use(LimitationMode.Moderate))
            Assert.DoesNotThrow(() => LimitationPolicy.Enforce(id));
        using (LimitationPolicy.UsePermissive())
            Assert.DoesNotThrow(() => LimitationPolicy.Enforce(id));
    }

    [TestCaseSource(nameof(GuardedIds))]
    public void Catalog_ReportFile_Exists(string id)
    {
        var repoRoot = FindRepoRoot();
        if (repoRoot is null)
            Assert.Ignore("repo root not resolvable from test output dir");
        var path = Path.Combine(repoRoot!, LimitationCatalog.Get(id).ReportPath.Replace('/', Path.DirectorySeparatorChar));
        Assert.That(File.Exists(path), Is.True, $"report should exist: {path}");
    }

    [Test]
    public void Catalog_Get_UnknownId_Throws()
        => Assert.Throws<KeyNotFoundException>(() => LimitationCatalog.Get("NOT-A-UNIT-999"));

    [Test]
    public void Exception_Message_CarriesIdRelatedWorkaroundAndPermissiveHint()
    {
        var ex = Assert.Throws<SeqeronLimitationException>(() =>
        {
            using (LimitationPolicy.UseStrict())
                LimitationPolicy.Enforce("DISORDER-REGION-001");
        })!;

        Assert.Multiple(() =>
        {
            Assert.That(ex.LimitationId, Is.EqualTo("DISORDER-REGION-001"));
            Assert.That(ex.Category, Is.EqualTo(LimitationCategory.DataBlocked));
            Assert.That(ex.Message, Does.Contain("DISORDER-REGION-001"));
            Assert.That(ex.Message, Does.Contain("Related to:"));
            Assert.That(ex.Message, Does.Contain("How to obtain the result:"));
            Assert.That(ex.Message, Does.Contain("PredictDisorderRegions"));
            Assert.That(ex.Message, Does.Contain("docs/Validation/reports/DISORDER-REGION-001.md"));
            Assert.That(ex.Message, Does.Contain("Permissive"));
        });
    }

    // ── 1) PARSE-FASTQ-001: ambiguous (overlap-confined) encoding ────────────────────────────

    private const string OverlapOnly = "@ABI"; // ASCII 64,65,66,73 — all in the Phred+33/Phred+64 overlap

    [Test]
    public void Fastq_FileLevelDetect_Ambiguous_ThrowsStrict_DefaultsPhred33Permissive()
    {
        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => QualityScoreAnalyzer.DetectEncoding(new[] { OverlapOnly }))!;
            Assert.That(ex.LimitationId, Is.EqualTo("PARSE-FASTQ-001"));
        }
        using (LimitationPolicy.UsePermissive())
        {
            var r = QualityScoreAnalyzer.DetectEncoding(new[] { OverlapOnly });
            Assert.That(r.Confidence, Is.EqualTo(QualityScoreAnalyzer.EncodingConfidence.Ambiguous));
            Assert.That(r.Encoding, Is.EqualTo(QualityScoreAnalyzer.QualityEncoding.Phred33));
        }
    }

    [Test]
    public void Fastq_PerStringDetectors_Ambiguous_ThrowStrict()
    {
        using (LimitationPolicy.UseStrict())
        {
            Assert.Throws<SeqeronLimitationException>(() => QualityScoreAnalyzer.DetectEncoding(OverlapOnly));
            Assert.Throws<SeqeronLimitationException>(() => FastqParser.DetectEncoding(OverlapOnly));
        }
        using (LimitationPolicy.UsePermissive())
        {
            Assert.That(QualityScoreAnalyzer.DetectEncoding(OverlapOnly),
                Is.EqualTo(QualityScoreAnalyzer.QualityEncoding.Phred33));
            Assert.That(FastqParser.DetectEncoding(OverlapOnly),
                Is.EqualTo(FastqParser.QualityEncoding.Phred33));
        }
    }

    [Test]
    public void Fastq_DefinitivelyResolvable_NeverThrows_EvenStrict()
    {
        // A character below ASCII 64 proves Phred+33; this is the ideal result, not a guess.
        using (LimitationPolicy.UseStrict())
        {
            Assert.DoesNotThrow(() => QualityScoreAnalyzer.DetectEncoding(new[] { "#@AB" }));
            Assert.DoesNotThrow(() => QualityScoreAnalyzer.DetectEncoding("#@AB"));
        }
    }

    // ── 2) CHROM-CENT-001: SF1-vs-SF2 (dimeric) ──────────────────────────────────────────────

    private static string DimericArray()
    {
        var reference = ChromosomeAnalyzer.LoadBundledAlphaSatelliteReference();
        string a = reference.First(r => r.Name == "ALRa").Sequence; // A-type monomer
        string b = reference.First(r => r.Name == "ALRb").Sequence; // B-type monomer
        return string.Concat(Enumerable.Repeat(a + b, 6));          // period-2 A->B dimer => SF1/SF2
    }

    [Test]
    public void Chromosome_Sf1OrSf2Dimeric_ThrowsStrict_ReturnsDimericPermissive()
    {
        string array = DimericArray();

        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => ChromosomeAnalyzer.AssignSuprachromosomalFamily(array))!;
            Assert.That(ex.LimitationId, Is.EqualTo("CHROM-CENT-001"));
        }
        using (LimitationPolicy.UsePermissive())
        {
            var r = ChromosomeAnalyzer.AssignSuprachromosomalFamily(array);
            Assert.That(r.Family, Is.EqualTo(ChromosomeAnalyzer.SuprachromosomalFamily.Sf1OrSf2Dimeric));
        }
    }

    // ── 3) DISORDER-REGION-001: uncalibrated confidence ──────────────────────────────────────

    private const string DisorderSeq =
        "MKSPSEEKQDSPSEEKQRSPSEEKQDSPSEEKQGSPSEEKQDSPSEEKQRSPSEEKQDSPSEEKQ";

    [Test]
    public void Disorder_PredictDisorder_ThrowsStrict()
    {
        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => DisorderPredictor.PredictDisorder(DisorderSeq))!;
            Assert.That(ex.LimitationId, Is.EqualTo("DISORDER-REGION-001"));
        }
    }

    [Test]
    public void Disorder_PredictDisorderRegions_NeverThrows_StripsConfidence_KeepsBoundaries()
    {
        DisorderPredictor.DisorderPredictionResult regionsOnly;
        DisorderPredictor.DisorderPredictionResult full;

        using (LimitationPolicy.UseStrict())
            regionsOnly = DisorderPredictor.PredictDisorderRegions(DisorderSeq); // must not throw

        using (LimitationPolicy.UsePermissive())
            full = DisorderPredictor.PredictDisorder(DisorderSeq);

        Assert.Multiple(() =>
        {
            // Confidence withheld (NaN) in the boundaries-only result.
            Assert.That(regionsOnly.DisorderedRegions.All(r => double.IsNaN(r.Confidence)), Is.True,
                "confidence is withheld (NaN) in PredictDisorderRegions");
            // Boundaries (the validated TOP-IDP result) are identical.
            Assert.That(regionsOnly.DisorderedRegions.Select(r => (r.Start, r.End, r.RegionType)),
                Is.EqualTo(full.DisorderedRegions.Select(r => (r.Start, r.End, r.RegionType))));
            Assert.That(regionsOnly.OverallDisorderContent, Is.EqualTo(full.OverallDisorderContent));
        });
    }

    // ── 4) MIRNA-TARGET-001: partial context++ (omitted features) ────────────────────────────

    private static (string mrna, MiRnaAnalyzer.MiRna miRna, MiRnaAnalyzer.TargetSite site) ContextFixture()
    {
        var let7a = MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
        string mrna = "GGGGG" + "CUACCUC" + "A" + "GGGGG"; // canonical let-7a 8mer site (RC of seed + A)
        var site = MiRnaAnalyzer.FindTargetSites(mrna, let7a, minScore: 0.0).First();
        return (mrna, let7a, site);
    }

    [Test]
    public void ContextPlusPlus_PartialScore_ThrowsStrict_ReturnsPartialPermissive()
    {
        var (mrna, miRna, site) = ContextFixture();

        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, miRna, site))!;
            Assert.That(ex.LimitationId, Is.EqualTo("MIRNA-TARGET-001"));
        }
        using (LimitationPolicy.UsePermissive())
        {
            var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, miRna, site);
            Assert.That(ctx.OmittedFeatures, Is.Not.Empty, "no optional inputs supplied => partial");
        }
    }

    // ── 5) MIRNA-CLEAVAGE-001: approximate 3p (star) span ────────────────────────────────────

    private const string PriMiRna =
        "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

    [Test]
    public void DroshaDicer_ThrowsStrict_ReturnsCleavagePermissive()
    {
        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => MiRnaAnalyzer.PredictDroshaDicerCleavage(PriMiRna, 0))!;
            Assert.That(ex.LimitationId, Is.EqualTo("MIRNA-CLEAVAGE-001"));
        }
        using (LimitationPolicy.UsePermissive())
        {
            var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(PriMiRna, 0);
            Assert.That(cut, Is.Not.Null);
        }
    }

    [Test]
    public void DroshaDicer_InvalidInput_ReturnsNull_DoesNotThrow_EvenStrict()
    {
        using (LimitationPolicy.UseStrict())
        {
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage("", 0), Is.Null);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(PriMiRna, -1), Is.Null);
        }
    }

    // ── incomplete-group behavioural checks (throw only under Strict; allowed under Moderate) ──

    [Test]
    public void Mgb_QualitativeRules_ThrowsStrict_AllowedModerate()
    {
        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACGT"))!;
            Assert.That(ex.LimitationId, Is.EqualTo("PROBE-DESIGN-001"));
        }
        using (LimitationPolicy.Use(LimitationMode.Moderate))
            Assert.DoesNotThrow(() => ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACGT"));
    }

    [Test]
    public void CheckM_DomainLevel_ThrowsStrict_AllowedModerate()
    {
        var markerSets = new[] { new MetagenomicsAnalyzer.MarkerSet(new[] { "m1", "m2" }) };
        var counts = new Dictionary<string, int> { ["m1"] = 1, ["m2"] = 1 };

        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(markerSets, counts))!;
            Assert.That(ex.LimitationId, Is.EqualTo("META-BIN-001"));
        }
        using (LimitationPolicy.Use(LimitationMode.Moderate))
            Assert.DoesNotThrow(() => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(markerSets, counts));
    }

    [Test]
    public void EstimateTumorPurity_ThrowsStrict_AllowedModerate()
    {
        using (LimitationPolicy.UseStrict())
        {
            var ex = Assert.Throws<SeqeronLimitationException>(
                () => ImmuneAnalyzer.EstimateTumorPurity(0.5))!;
            Assert.That(ex.LimitationId, Is.EqualTo("ONCO-IMMUNE-001"));
        }
        using (LimitationPolicy.Use(LimitationMode.Moderate))
            Assert.DoesNotThrow(() => ImmuneAnalyzer.EstimateTumorPurity(0.5));
    }

    // ── helper ───────────────────────────────────────────────────────────────────────────────

    private static string? FindRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir is not null)
        {
            if (Directory.Exists(Path.Combine(dir.FullName, "docs", "Validation", "reports")))
                return dir.FullName;
            dir = dir.Parent;
        }
        return null;
    }
}
