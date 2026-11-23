using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-7 (B09 F172): the EMBOSS sequence-input qualifiers <c>-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein</c> of the ported programs, against the installed EMBOSS 6.6.0
/// (TestData/Alignment/emboss_seqopts_oracle.tsv, scripts/alignment/generate_emboss_seqopts_oracle.py: 180 rows per program with random begin / end (inside, at and beyond the ends, negative, zero, reversed
/// ranges), reverse, case, circular and type overrides - also the combinations the program refuses). The rows are evaluated by <see cref="EmbossSeqOptsProducer"/> through the shared reader
/// (<see cref="EmbossInputSequence"/>): what each program REPORTS for ranged and reverse-complemented input (positions, strands, numbering direction, alignment offsets) is compared, not our reading of the C code.
/// </summary>
[TestFixture]
public class EmbossSequenceQualifiers_Oracle_Tests
{
    private static readonly Lazy<List<EmbossSeqOptsProducer.Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_seqopts_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        return EmbossSeqOptsProducer.ReadRows(reader);
    });

    private static void Check(string program, int minimum)
    {
        var rows = Oracle.Value.Where(r => r.Program == program).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum), $"oracle rows of {program}");
        var failures = new List<string>();
        int errors = 0, ranged = 0, reversed = 0;
        foreach (var row in rows)
        {
            string expected = EmbossSeqOptsProducer.Expected(row);
            string got;
            try
            {
                got = EmbossSeqOptsProducer.Produce(row);
            }
            catch (Exception e) when (e is not AssertionException)
            {
                got = $"CRASH {e.GetType().Name}: {e.Message}";
            }

            if (got != expected)
                failures.Add($"{row}\n   expected {expected}\n   actual   {got}");
            if (expected == "ERROR")
                errors++;
            var q = EmbossSeqOptsProducer.Quals(row.Quals1);
            if (q.Begin != 0 || q.End != 0)
                ranged++;
            if (q.Reverse)
                reversed++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(4)));
        Assert.That(ranged, Is.GreaterThanOrEqualTo(minimum / 2), "rows with a begin or an end");
        Assert.That(reversed, Is.GreaterThanOrEqualTo(minimum / 5), "rows with -sreverse");
        Assert.That(errors, Is.LessThan(rows.Count / 2), "the refused combinations stay a minority");
    }

    [Test] public void Needle_QualifiersLikeEmboss_OnAllRows() => Check("needle", 150);
    [Test] public void Water_QualifiersLikeEmboss_OnAllRows() => Check("water", 150);
    [Test] public void Stretcher_QualifiersLikeEmboss_OnAllRows() => Check("stretcher", 150);
    [Test] public void Matcher_QualifiersLikeEmboss_OnAllRows() => Check("matcher", 150);
    [Test] public void Supermatcher_QualifiersLikeEmboss_OnAllRows() => Check("supermatcher", 150);
    [Test] public void Needleall_QualifiersLikeEmboss_OnAllRows() => Check("needleall", 150);
    [Test] public void Merger_QualifiersLikeEmboss_OnAllRows() => Check("merger", 150);
    [Test] public void Megamerger_QualifiersLikeEmboss_OnAllRows() => Check("megamerger", 150);
    [Test] public void Seqmatchall_QualifiersLikeEmboss_OnAllRows() => Check("seqmatchall", 150);
    [Test] public void Fuzznuc_QualifiersLikeEmboss_OnAllRows() => Check("fuzznuc", 150);
    [Test] public void Fuzzpro_QualifiersLikeEmboss_OnAllRows() => Check("fuzzpro", 150);
    [Test] public void Fuzztran_QualifiersLikeEmboss_OnAllRows() => Check("fuzztran", 150);
    [Test] public void Etandem_QualifiersLikeEmboss_OnAllRows() => Check("etandem", 150);
    [Test] public void Equicktandem_QualifiersLikeEmboss_OnAllRows() => Check("equicktandem", 150);
    [Test] public void Dotmatcher_QualifiersLikeEmboss_OnAllRows() => Check("dotmatcher", 150);
    [Test] public void Distmat_QualifiersLikeEmboss_OnAllRows() => Check("distmat", 150);

    [Test]
    public void Read_Range_FollowsAjSeqTrim()
    {
        // ajSeqTrim + ajCvtSposToPosStart measured on needle 6.6.0: negative positions count from the end, both ends are clipped into the sequence, the end never before the start
        var q = new EmbossSequenceQualifiers { Begin = 3, End = 8 };
        var v = EmbossInputSequence.Read("ACGTACGTAC", q);
        Assert.That((v.Residues, v.Offset, v.Offend, v.Reversed), Is.EqualTo(("GTACGT", 2, 2, false)));
        v = EmbossInputSequence.Read("ACGTACGTAC", new EmbossSequenceQualifiers { Begin = -3 });
        Assert.That((v.Residues, v.Offset, v.Offend), Is.EqualTo(("TAC", 7, 0)));
        v = EmbossInputSequence.Read("ACGTACGTAC", new EmbossSequenceQualifiers { Begin = 100 });
        Assert.That((v.Residues, v.Offset), Is.EqualTo(("C", 9)));          // begin beyond the end: the last residue
        v = EmbossInputSequence.Read("ACGTACGTAC", new EmbossSequenceQualifiers { Begin = 8, End = 3 });
        Assert.That((v.Residues, v.Offset), Is.EqualTo(("T", 7)));          // end before begin: one residue at the begin
        v = EmbossInputSequence.Read("ACGTACGTAC", new EmbossSequenceQualifiers { End = -3 });
        Assert.That((v.Residues, v.Offend), Is.EqualTo(("ACGTACGT", 2)));
    }

    [Test]
    public void Read_Reverse_MirrorsTheRangeAndReverseComplementsOnlyNucleotides()
    {
        // needle -sbegin1 3 -send1 8 -sreverse1: the reverse complement of original positions 3..8, reported in original numbering counting down
        var v = EmbossInputSequence.Read("ACGTACGTAC", new EmbossSequenceQualifiers { Begin = 3, End = 8, Reverse = true });
        Assert.That(v.Residues, Is.EqualTo("ACGTAC"));                       // revcomp(GTACGT)
        Assert.That(v.Reversed, Is.True);
        Assert.That(v.ReportedRange(0, 6), Is.EqualTo((8, 3)).Or.EqualTo((3, 8)));
        // a protein is never reversed; -sprotein makes a nucleotide-looking sequence one
        Assert.That(EmbossInputSequence.Read("ACDEFGHIK", new EmbossSequenceQualifiers { Reverse = true }).Reversed, Is.False);
        Assert.That(EmbossInputSequence.Read("ACGTACGT", new EmbossSequenceQualifiers { Reverse = true, Protein = true }).Reversed, Is.False);
        Assert.That(EmbossInputSequence.Read("ACDEFGHIK", new EmbossSequenceQualifiers { Reverse = true, Nucleotide = true }).Reversed, Is.True);
    }

    [Test]
    public void ReverseComplement_IsTheAjBaseCompTable_CasePreserved()
    {
        // seqret -sreverse -snucleotide on A-Z a-z (6.6.0): every IUPAC code complements, U and T give A, X gives N, J I L F E O P Q Z stay
        const string input = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        Assert.That(EmbossInputSequence.ReverseComplement(input), Is.EqualTo("zrnwbaasyqponklmjidcfehgvtZRNWBAASYQPONKLMJIDCFEHGVT"));
    }

    [Test]
    public void Read_ForcedType_RepairsOrRefusesLikeAjSeqTypeCheckIn()
    {
        // fuzznuc (type dna): -snucleotide on protein letters turns every non-nucleotide letter into N (U -> T, X -> N); -sprotein is refused; a protein input refuses -snucleotide
        var q = new EmbossSequenceQualifiers { Nucleotide = true };
        Assert.That(EmbossInputSequence.Read("ACDEFGHIKLMNPQRSTUVWXYZ", q, EmbossSequenceType.Dna).Residues, Is.EqualTo("ACDNNGHNKNMNNNRSTTVWNYN"));
        Assert.Throws<ArgumentException>(() => EmbossInputSequence.Read("ACGT", new EmbossSequenceQualifiers { Protein = true }, EmbossSequenceType.Dna));
        Assert.Throws<ArgumentException>(() => EmbossInputSequence.Read("ACGT", q, EmbossSequenceType.Protein));
        Assert.Throws<ArgumentException>(() => EmbossInputSequence.Read("ACDEFG", null, EmbossSequenceType.Nucleotide));         // letters that are not nucleotides
        Assert.Throws<ArgumentException>(() => EmbossInputSequence.Read("--..", null));                                       // nothing left after the reader
    }

    [Test]
    public void Parse_CommandLineFragment_PerInput()
    {
        var q1 = EmbossSequenceQualifiers.Parse("-sbegin1 3 -send1 -2 -sreverse1 -slower1 -sbegin2 7 -nosreverse2", "1");
        Assert.That((q1.Begin, q1.End, q1.Reverse, q1.Lower), Is.EqualTo((3, -2, true, true)));
        var q2 = EmbossSequenceQualifiers.Parse("-sbegin1 3 -send1 -2 -sreverse1 -slower1 -sbegin2 7 -nosreverse2 -sprotein2 -sformat2 fasta", "2");
        Assert.That((q2.Begin, q2.End, q2.Reverse, q2.Protein), Is.EqualTo((7, 0, false, true)));
        Assert.That(EmbossSequenceQualifiers.Parse("-sbegin 4 -scircular -snucleotide").IsDefault, Is.False);
        Assert.That(EmbossSequenceQualifiers.Parse(null).IsDefault, Is.True);
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.Parse("-sbegin1"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.Parse("-gapopen 10"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.Parse("sbegin1 3"));
    }

    [Test]
    public void ReverseQualifier_NumbersTheOriginalStrandDownwards_ExampleOfTheAudit()
    {
        // audit A5-7 (needle 6.6.0): `needle a.fa b.fa -sbegin1 3 -send1 20 -sreverse2 -slower1` aligns a[3..20] (lower case) with the reverse complement of b, a numbered 3..18, b counting down
        var a = new EmbossSequence("a", "ACGTACGTTAGCCGATTAGGCATCGATCGGATCG", "", new EmbossSequenceQualifiers { Begin = 3, End = 20, Lower = true });
        var b = new EmbossSequence("b", "CGGATCGATCAGGACTATCGGCTAACGTACGTTAGC", "", new EmbossSequenceQualifiers { Reverse = true });
        var result = SequenceAligner.NeedleProgram(a, b);
        var pair = result.Alignments.Single();
        Assert.That(pair.First.Row.Replace("-", ""), Is.EqualTo("gtacgttagccgattagg"));
        Assert.That((pair.First.Offset, pair.First.Reversed), Is.EqualTo((2, false)));
        Assert.That(pair.Second.Reversed, Is.True);
        string text = EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, result.Alignments, new EmbossAlignOptions { Program = "needle", MatrixName = result.MatrixName, Matrix = SubstitutionMatrix.Get(result.MatrixName), Global = true, GapPenalty = 10, ExtendPenalty = 0.5 });
        Assert.That(text, Does.Contain("a                  3 ------gtacgttagccgattagg-------------     20"));
        Assert.That(text, Does.Contain("b                 36 GCTAACGTACGTTAGCCGA-TAGTCCTGATCGATCCG      1"));
    }
}
