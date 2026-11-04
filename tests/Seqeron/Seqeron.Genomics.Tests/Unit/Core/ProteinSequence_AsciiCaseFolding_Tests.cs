using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Review 2026-09, B02 F31 (TRANS-PROT-001; B01 request from SEQ-VALID-001 F3/F18/F19): B02-owned entry
/// points upper-case ASCII letters only. <c>ToUpperInvariant</c> folds U+017F 'ſ' (LATIN SMALL LETTER
/// LONG S) to 'S' — a valid amino acid and the IUPAC nucleotide code S (G/C) — so "Mſ" was accepted as
/// "MS" and "GGſ" translated as GGS → Gly. U+0131 'ı' (dotless i) is locked as well (Python
/// <c>str.upper</c> folds it to 'I').
/// References: scikit-bio 0.7.4 <c>Protein("Mſ", lowercase=True)</c>, <c>Protein("Mı", lowercase=True)</c>
/// and <c>DNA("GſA", lowercase=True).translate()</c> raise (ASCII-only); Biopython 1.88
/// <c>Seq("GſA").translate()</c> / <c>Seq("AıG").translate()</c> raise; EMBOSS pepstats skips 'ſ'
/// (MSKſ: Residues = 3, MW 364.46, pI 9.15 — identical to MSK).
/// </summary>
[TestFixture]
public class ProteinSequence_AsciiCaseFolding_Tests
{
    private const char LongS = 'ſ';     // ſ → 'S' under ToUpperInvariant
    private const char DotlessI = 'ı';  // ı → 'I' under Python str.upper

    [TestCase("Mſ")]
    [TestCase("Mı")]
    [TestCase("ſ")]
    [TestCase("mkſk")]
    [Description("F31: ProteinSequence rejects non-ASCII letters (scikit-bio Protein raises)")]
    public void Constructor_NonAsciiLetter_Throws(string input)
    {
        var ex = Assert.Throws<ArgumentException>(() => _ = new ProteinSequence(input));
        Assert.Multiple(() =>
        {
            // The offending character itself is reported (not its folded 'S' / 'I').
            Assert.That(ex!.Message, Does.Contain($"'{input.First(c => c > 127)}'"));
            Assert.That(ProteinSequence.TryCreate(input, out var result), Is.False);
            Assert.That(result, Is.Null);
        });
    }

    [Test]
    [Description("F31: ASCII lower case still folds (scikit-bio Protein('ms', lowercase=True) == 'MS')")]
    public void Constructor_AsciiLowerCase_StillUpperCased()
    {
        Assert.That(new ProteinSequence("msik").Sequence, Is.EqualTo("MSIK"));
        Assert.That(ProteinSequence.TryCreate("ms", out var p), Is.True);
        Assert.That(p!.Sequence, Is.EqualTo("MS"));
    }

    [Test]
    [Description("F31: FindMotif folds the pattern ASCII-only — 'ſ' does not match S, 'ı' does not match I")]
    public void FindMotif_NonAsciiPattern_NotFolded()
    {
        var protein = new ProteinSequence("MSKSIL");
        Assert.Multiple(() =>
        {
            Assert.That(protein.FindMotif(LongS.ToString()), Is.Empty);
            Assert.That(protein.FindMotif(DotlessI.ToString()), Is.Empty);
            Assert.That(protein.FindMotif("k" + LongS), Is.Empty);
            Assert.That(protein.FindMotif("s"), Is.EqualTo(new[] { 1, 3 }));
            Assert.That(protein.FindMotif("ks"), Is.EqualTo(new[] { 2 }));
        });
    }

    [TestCase("ATGGGſ")]
    [TestCase("ATGGGı")]
    [Description("F31: Translator.Translate rejects a non-ASCII codon symbol (Biopython/scikit-bio translate raise); " +
                 "before the fix ATGGGſ → \"MG\" (GGS read as Gly)")]
    public void Translator_NonAsciiCodon_Throws(string dna)
    {
        Assert.Throws<ArgumentException>(() => Translator.Translate(dna));
        Assert.That(Translator.Translate("atgggc").Sequence, Is.EqualTo("MG"));
    }

    [Test]
    [Description("F31: GeneticCode.Translate / GetCodonsForAminoAcid fold ASCII only")]
    public void GeneticCode_NonAsciiLetters_NotFolded()
    {
        var code = GeneticCode.Standard;
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => code.Translate("GG" + LongS));
            Assert.Throws<ArgumentException>(() => code.Translate("GG" + DotlessI));
            Assert.That(code.Translate("ggc"), Is.EqualTo('G'));
            Assert.That(code.Translate("ggs"), Is.EqualTo('G'));   // ASCII IUPAC S = G/C: GGG, GGC both Gly
            Assert.That(code.GetCodonsForAminoAcid(LongS), Is.Empty);
            Assert.That(code.GetCodonsForAminoAcid(DotlessI), Is.Empty);
            Assert.That(code.GetCodonsForAminoAcid('s').Count(), Is.EqualTo(6));  // Ser: UCN + AGY
        });
    }

    [Test]
    [Description("F31: CodonOptimizer does not read 'ſ' as IUPAC S — the triplet stays untranslatable and untouched")]
    public void CodonOptimizer_NonAsciiCodon_NotFolded()
    {
        var result = CodonOptimizer.OptimizeSequence("ATGGG" + LongS, CodonOptimizer.EColiK12);
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalSequence, Is.EqualTo("AUGGG" + LongS));
            Assert.That(result.ProteinSequence, Is.EqualTo("MX"));
            Assert.That(result.OptimizedSequence, Does.EndWith("GG" + LongS));
            Assert.Throws<ArgumentException>(() =>
                CodonOptimizer.RemoveRestrictionSites("ATGGAATTC", new[] { "GAA" + LongS + "TC" }, CodonOptimizer.EColiK12));
            Assert.That(CodonOptimizer.CreateCodonUsageTable("t", new Dictionary<string, double> { ["gg" + LongS] = 1.0 })
                .CodonFrequencies.ContainsKey("GGS"), Is.False);
        });
    }

    [Test]
    [Description("F31: ProteinPhysicochemistry skips 'ſ' like EMBOSS pepstats (MSKſ ≡ MSK: Residues 3, MW 364.46, pI 9.15)")]
    public void ProteinPhysicochemistry_NonAsciiLetter_Skipped()
    {
        string withLongS = "MSK" + LongS;
        Assert.Multiple(() =>
        {
            Assert.That(ProteinPhysicochemistry.MolecularWeight(withLongS),
                Is.EqualTo(ProteinPhysicochemistry.MolecularWeight("MSK")));
            Assert.That(ProteinPhysicochemistry.MolecularWeight(withLongS), Is.EqualTo(364.46).Within(0.01));
            Assert.That(ProteinPhysicochemistry.Gravy(withLongS), Is.EqualTo(ProteinPhysicochemistry.Gravy("MSK")));
            Assert.That(ProteinPhysicochemistry.IsoelectricPoint(withLongS),
                Is.EqualTo(ProteinPhysicochemistry.IsoelectricPoint("MSK")));
            Assert.That(ProteinPhysicochemistry.MolecularWeight("msk"),
                Is.EqualTo(ProteinPhysicochemistry.MolecularWeight("MSK")));
        });
    }
}
