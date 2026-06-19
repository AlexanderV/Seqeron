using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the MiRNA area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-SEED-001 — miRNA seed region / target seed matching (MiRNA).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 74.
///
/// API under test (MiRnaAnalyzer.GetSeedSequence / CreateMiRna / FindTargetSites):
///   The seed is miRNA nucleotides 2–8 (0-based substring [1,8)). Target prediction scans an
///   mRNA for the reverse complement of the seed (the 6-mer core, optionally extended to
///   7mer/8mer); a seed match is therefore necessary for a predicted site.
///
/// Relations (derived from the seed definition, NOT from output):
///   • INV  (3' end changes don't affect seed): the seed is fixed to positions 2–8, so mutating
///          any nucleotide from position 9 onward leaves the seed unchanged.
///   • INV  (seed extraction deterministic): extraction is a pure, case-normalising function,
///          and CreateMiRna additionally maps T→U, so DNA and RNA spellings give the same seed.
///   • COMP (seed match ⊂ target prediction): a planted seed reverse-complement yields a target
///          site whose span contains the 6-mer seed-core match, and an mRNA with no seed match
///          yields no sites — target prediction is gated on the seed.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MiRnaMetamorphicTests
{
    // let-7a-5p (22 nt). Seed (positions 2–8) = "GAGGUAG".
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    #region MIRNA-SEED-001 INV — 3' end changes leave the seed unchanged

    [Test]
    [Description("INV: the seed is positions 2–8, so mutating any nucleotide from position 9 onward (the 3' end) does not change the extracted seed.")]
    public void Seed_3PrimeEndChanges_DoNotAffectSeed()
    {
        string baseSeed = MiRnaAnalyzer.GetSeedSequence(Let7a);
        baseSeed.Should().Be("GAGGUAG", because: "let-7a nucleotides 2–8 spell GAGGUAG");

        // Keep the first 8 nucleotides (which contain the seed) and rewrite the 3' tail.
        foreach (var tail in new[] { "CCCCCCCCCCCCCC", "AAAAAAAAAAAAAA", "GUGUGUGUGUGUGU" })
        {
            string mutated = Let7a.Substring(0, 8) + tail;
            MiRnaAnalyzer.GetSeedSequence(mutated).Should().Be(baseSeed,
                because: "the seed window (positions 2–8) is untouched by 3'-end edits");
            MiRnaAnalyzer.CreateMiRna("m", mutated).SeedSequence.Should().Be(baseSeed,
                because: "CreateMiRna derives the seed from the same fixed window");
        }
    }

    #endregion

    #region MIRNA-SEED-001 INV — seed extraction is deterministic and normalising

    [Test]
    [Description("INV: seed extraction is a pure, case-normalising function; CreateMiRna also maps T→U so DNA and RNA spellings of the same miRNA yield the same seed.")]
    public void Seed_Extraction_IsDeterministicAndNormalising()
    {
        MiRnaAnalyzer.GetSeedSequence(Let7a)
            .Should().Be(MiRnaAnalyzer.GetSeedSequence(Let7a), because: "extraction is deterministic");

        MiRnaAnalyzer.GetSeedSequence(Let7a.ToLowerInvariant())
            .Should().Be(MiRnaAnalyzer.GetSeedSequence(Let7a), because: "extraction upper-cases its input");

        // DNA spelling (T for U) of let-7a must yield the same RNA seed via CreateMiRna.
        string dna = Let7a.Replace('U', 'T');
        MiRnaAnalyzer.CreateMiRna("dna", dna).SeedSequence
            .Should().Be(MiRnaAnalyzer.CreateMiRna("rna", Let7a).SeedSequence,
                because: "CreateMiRna transcribes T→U before extracting the seed");
    }

    #endregion

    #region MIRNA-SEED-001 COMP — a seed match is necessary and contained in a target site

    [Test]
    [Description("COMP: a planted seed reverse-complement produces a target site whose span contains the 6-mer seed core, and an mRNA with no seed match produces no sites.")]
    public void Seed_TargetPrediction_RequiresAndContainsSeedMatch()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence); // RC of positions 2–8 (7 nt)
        string sixmerCore = seedRc.Substring(1, 6);                              // RC of positions 2–7

        // Plant the full seed-RC with an upstream position-8 match and a downstream 'A'
        // (opposite miRNA position 1) so an 8mer site forms and clears the default score gate.
        string mrna = "GGGGG" + seedRc + "A" + "GGGGG";

        var sites = MiRnaAnalyzer.FindTargetSites(mrna, mirna).ToList();

        sites.Should().NotBeEmpty(because: "a seed reverse-complement match yields a predicted target site");
        sites.Should().OnlyContain(s => mrna.Substring(s.Start, s.End - s.Start + 1).Contains(sixmerCore),
            because: "every predicted site's span contains the 6-mer seed-core match it was built around");

        // No seed match anywhere ⇒ no prediction: target prediction is gated on the seed.
        MiRnaAnalyzer.FindTargetSites(new string('A', 40), mirna, minScore: 0.0)
            .Should().BeEmpty(because: "without a seed reverse-complement match there is nothing to predict");
    }

    #endregion
}
