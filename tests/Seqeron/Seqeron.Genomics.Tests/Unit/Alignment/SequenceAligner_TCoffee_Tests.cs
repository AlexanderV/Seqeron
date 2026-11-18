using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-13, B09 F58): the progressive stage of the T-Coffee 13.41 program
/// (<see cref="SequenceAligner.MultipleAlignTCoffee"/>) for an explicit primary library and an explicit guide tree.
///
/// Oracle (TestData/Alignment/tcoffee_oracle.tsv, scripts/alignment/generate_tcoffee_oracle.py): the rows written by the installed
/// <c>t_coffee</c> (Debian package t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1) for
/// <c>t_coffee -seq s.fa -lib l.tc_lib -usetree t.dnd -relax_lib {1|0} -output fasta_aln -n_core 1</c> on random DNA sets (3-6 sequences),
/// with libraries built by this repository's global + local alignments ("rand") and by the program's own
/// <c>-method slow_pair,lalign_id_pair</c> ("binary"). Every row of the program must be reproduced exactly.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffee_Tests
{
    private sealed record Row(string Kind, bool Relax, string[] Seqs, string Tree, string Library, string[] Rows);

    private static List<Row> ReadOracle(string file = "tcoffee_oracle.tsv")
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1] == "1", f[2].Split('|'), f[3], f[4], f[5].Split('|')));
        }

        return rows;
    }

    private static TCoffeeLibrary ParseLibrary(string text, int k)
    {
        var lib = new TCoffeeLibrary(k);
        foreach (var e in text.Split(';', StringSplitOptions.RemoveEmptyEntries))
        {
            var t = e.Split(':').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
            lib.Add(t[0], t[1] - 1, t[2], t[3] - 1, t[4]); // the file uses 1-based residues
        }

        return lib;
    }

    private static TCoffeeTreeNode ParseTree(string text)
    {
        int pos = 0;
        TCoffeeTreeNode Node()
        {
            if (text[pos] == '(')
            {
                pos++;
                var l = Node();
                pos++; // ','
                var r = Node();
                pos++; // ')'
                return TCoffeeTreeNode.Join(l, r);
            }

            int start = pos;
            while (pos < text.Length && char.IsDigit(text[pos])) pos++;
            return TCoffeeTreeNode.Leaf(int.Parse(text.AsSpan(start, pos - start), CultureInfo.InvariantCulture));
        }

        return Node();
    }

    [Test]
    public void OracleFile_IsComplete()
    {
        var rows = ReadOracle();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(780));
        Assert.That(rows.Count(r => r.Relax), Is.GreaterThanOrEqualTo(390));
        Assert.That(rows.Count(r => !r.Relax), Is.GreaterThanOrEqualTo(390));
        Assert.That(rows.Select(r => r.Kind).Distinct().Order(), Is.EqualTo(new[] { "binary", "rand" }));
        Assert.That(rows.Select(r => r.Seqs.Length).Distinct().Order(), Is.EqualTo(new[] { 3, 4, 5, 6 }));
    }

    [Test]
    public void MultipleAlignTCoffee_EveryRowEqualsTheProgram()
    {
        var failures = new List<string>();
        var rows = ReadOracle();
        foreach (var o in rows)
        {
            var got = SequenceAligner.MultipleAlignTCoffee(o.Seqs, ParseLibrary(o.Library, o.Seqs.Length), ParseTree(o.Tree), o.Relax);
            if (!got.SequenceEqual(o.Rows))
                failures.Add($"{o.Kind} relax={o.Relax} {string.Join('|', o.Seqs)} {o.Tree}: got {string.Join('|', got)} expected {string.Join('|', o.Rows)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count} sets differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>
    /// Audit A1-14: the progressive stage is alphabet-agnostic (the library carries the weights) - 160 random PROTEIN sets (3-6 sequences of
    /// 3-11 residues), explicit libraries ("rand": identity alignments of this repository; "binary": the program's own BLOSUM-based protein
    /// library) and explicit random trees: every row of <c>t_coffee</c> is reproduced.
    /// </summary>
    [Test]
    public void MultipleAlignTCoffee_Protein_EveryRowEqualsTheProgram()
    {
        var rows = ReadOracle("tcoffee_protein_oracle.tsv");
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(320));
        Assert.That(rows.SelectMany(r => r.Seqs).SelectMany(s => s).Distinct().Count(), Is.GreaterThan(15), "a real protein alphabet");
        var failures = new List<string>();
        foreach (var o in rows)
        {
            var got = SequenceAligner.MultipleAlignTCoffee(o.Seqs, ParseLibrary(o.Library, o.Seqs.Length), ParseTree(o.Tree), o.Relax);
            if (!got.SequenceEqual(o.Rows))
                failures.Add($"{o.Kind} relax={o.Relax} {string.Join('|', o.Seqs)} {o.Tree}: got {string.Join('|', got)} expected {string.Join('|', o.Rows)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count} protein sets differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    [Test]
    public void MultipleAlignTCoffee_RowsAreTheInputsWithGaps()
    {
        foreach (var o in ReadOracle().Take(100))
        {
            var got = SequenceAligner.MultipleAlignTCoffee(o.Seqs, ParseLibrary(o.Library, o.Seqs.Length), ParseTree(o.Tree), o.Relax);
            Assert.That(got.Select(r => r.Replace("-", "")), Is.EqualTo(o.Seqs));
            Assert.That(got.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
        }
    }

    [Test]
    public void MultipleAlignTCoffee_RejectsInvalidInput()
    {
        var lib = new TCoffeeLibrary(2);
        var tree = TCoffeeTreeNode.Join(TCoffeeTreeNode.Leaf(0), TCoffeeTreeNode.Leaf(1));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignTCoffee(null!, lib, tree));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffee(new[] { "ACGT", "ACG", "AC" }, lib, tree));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffee(new[] { "ACGT", "" }, lib, tree));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffee(
            new[] { "ACGT", "ACG" }, lib, TCoffeeTreeNode.Join(TCoffeeTreeNode.Leaf(0), TCoffeeTreeNode.Leaf(0))));
        Assert.Throws<ArgumentException>(() => lib.Add(0, 0, 0, 1, 10));
        lib.Add(0, 9, 1, 0, 10);
        Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffee(new[] { "ACGT", "ACG" }, lib, tree));
    }

    [Test]
    public void ConsistencyEngine_DefaultIsThePaperFormula_AndUnchanged()
    {
        var seqs = new[] { "ACGTTACGT", "ACGACGT", "ACGTACGA" }.Select(s => new DnaSequence(s)).ToList();
        var plain = SequenceAligner.MultipleAlignConsistency(seqs);
        var paper = SequenceAligner.MultipleAlignConsistency(seqs, null, new MsaOptions { ConsistencyEngine = MsaConsistencyEngine.PaperExtension });
        Assert.That(new MsaOptions().ConsistencyEngine, Is.EqualTo(MsaConsistencyEngine.PaperExtension));
        Assert.That(paper.AlignedSequences, Is.EqualTo(plain.AlignedSequences));
    }

    [Test]
    public void ConsistencyEngine_TCoffeeLinked_ProducesAValidAlignmentThroughTheProgramStage()
    {
        var raw = new[] { "ACGTTACGT", "ACGACGT", "ACGTACGA", "TTACGTAG" };
        var seqs = raw.Select(s => new DnaSequence(s)).ToList();
        var result = SequenceAligner.MultipleAlignConsistency(seqs, null, new MsaOptions { ConsistencyEngine = MsaConsistencyEngine.TCoffeeLinked });
        Assert.That(result.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(raw));
        Assert.That(result.AlignedSequences.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
        var again = SequenceAligner.MultipleAlignConsistency(seqs, null, new MsaOptions { ConsistencyEngine = MsaConsistencyEngine.TCoffeeLinked });
        Assert.That(again.AlignedSequences, Is.EqualTo(result.AlignedSequences));
    }
}
