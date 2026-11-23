using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-9 / A2-10, B09): the T-Coffee 13.41 methods that run an EXTERNAL program, reproduced from T-Coffee's side of the interface.
/// Oracle: TestData/Alignment/tcoffee_external_oracle.tsv (scripts/alignment/generate_tcoffee_external_oracle.py): the installed <c>t_coffee</c> on random sets
/// together with what the installed <c>TMalign</c> / <c>mustang</c> (on synthetic CA-only PDB files) or <c>blastall</c> (a local BLAST database of the sequences) wrote for the same inputs.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffeeExternal_Tests
{
    internal sealed record Case(string Feature, string Flags, string[] Seqs, string Rows, JsonElement Extra);

    /// <summary>The parsed oracle, ONE copy per fixture (audit A3-15): the method used to re-read and re-parse the whole embedded TSV (6.0 MB / 3.2 MB as UTF-16 strings, 20 321 / several thousand rows) on every call, i.e. in every test.</summary>
    internal static List<Case> Oracle() => OracleRows.Value;

    private static readonly Lazy<List<Case>> OracleRows = new(ReadOracleRows);

    private static List<Case> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.tcoffee_external_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            list.Add(new Case(f[0], f[1], f[2].Split('|'), f[3], JsonDocument.Parse(f[4]).RootElement.Clone()));
        }

        return list;
    }

    internal static Dictionary<TCoffeePairMethod, IReadOnlyDictionary<(int, int), string>> Outputs(Case c)
    {
        var map = new Dictionary<TCoffeePairMethod, IReadOnlyDictionary<(int, int), string>>();
        foreach (var tool in c.Extra.GetProperty("out").EnumerateObject())
        {
            var byPair = new Dictionary<(int, int), string>();
            foreach (var p in tool.Value.EnumerateObject())
            {
                var ij = p.Name.Split(',');
                byPair[(int.Parse(ij[0]), int.Parse(ij[1]))] = p.Value.GetString()!;
            }

            map[tool.Name == "TMalign_pair" ? TCoffeePairMethod.TMalignPair : TCoffeePairMethod.MustangPair] = byPair;
        }

        return map;
    }

    private static void AssertStructureFeature(string feature, int minSets)
    {
        var cases = Oracle().Where(c => c.Feature == feature).ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(minSets), $"oracle sets of '{feature}'");
        var failures = new List<string>();
        int withStructure = 0, hijacked = 0;
        foreach (var c in cases)
        {
            var pdb = c.Extra.GetProperty("pdb").EnumerateArray().Select(e => e.GetString() is { Length: > 0 } s ? s : null).ToArray();
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with
            {
                StructureTemplates = pdb,
                StructureAlignerOutputs = Outputs(c),
            };
            if (pdb.Count(p => p != null) >= 2) withStructure++; else hijacked++;
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
            if (got != c.Rows) failures.Add($"'{c.Flags}' {string.Join('|', c.Seqs)} pdb {string.Join(',', pdb.Select(p => p ?? "-"))}: got {got} expected {c.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} '{feature}' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
        Assert.That(withStructure, Is.GreaterThanOrEqualTo(minSets / 2), "sets with at least two structures");
    }

    /// <summary>A2-10: <c>-method TMalign_pair</c> (+ slow_pair / proba_pair, relaxation off, UPGMA): pairs of two structures take the alignment TM-align wrote, the others fall back to proba_pair.</summary>
    [Test]
    public void TMalignPair_EveryRowEqualsTheProgram() => AssertStructureFeature("TMalign_pair", 150);

    /// <summary>A2-10: <c>-method mustang_pair</c> (+ slow_pair / proba_pair / TMalign_pair, relaxation off): the alignment MUSTANG wrote (results.afasta) is the library of a pair of two structures.</summary>
    [Test]
    public void MustangPair_EveryRowEqualsTheProgram() => AssertStructureFeature("mustang_pair", 150);

    /// <summary>
    /// A2-12: <c>-method sap_pair</c>: SAP is not installed, so the program runs TM-align instead ("SAP is not installed. TMalign will be used instead"); the rows equal <c>TMalign_pair</c> with the same structures (150 sets).
    /// <c>-mode 3dcoffee</c> (= <c>-in Msap_pair -template_file SELF_P_</c>) is the method list sap_pair ALONE: the 8 committed sets (the sequence names q000, ... are no PDB entries, so no structure is found and every pair is aligned by proba_pair)
    /// equal it; they differ from <c>-method proba_pair,sap_pair</c> on one set because that list has the proba_pair library twice (B09 F311: the "not understood" eighth set of the first run). The runs of the real SAP
    /// are in <see cref="SequenceAligner_TCoffeeSapDaliRuns_Tests"/>.
    /// </summary>
    [Test]
    public void SapPair_IsTMalignWhenSapIsNotInstalled()
    {
        AssertStructureFeature("sap_pair", 140);
        var threeD = Oracle().Where(c => c.Feature == "sap_3dcoffee").ToList();
        Assert.That(threeD.Count, Is.EqualTo(8));
        int equalToSapAlone = 0, equalToProbaPlusSap = 0;
        foreach (var c in threeD)
        {
            var pdb = c.Extra.GetProperty("pdb").EnumerateArray().Select(e => e.GetString() is { Length: > 0 } s ? s : null).ToArray();
            string alone = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, TCoffeeProgramOptions.FromCommandLine(c.Flags)).AlignedSequences);
            if (alone == c.Rows) equalToSapAlone++;
            var options = TCoffeeProgramOptions.FromCommandLine("-method proba_pair,sap_pair") with { StructureTemplates = pdb, StructureAlignerOutputs = Outputs(c) };
            if (string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences) == c.Rows) equalToProbaPlusSap++;
        }

        Assert.That(equalToSapAlone, Is.EqualTo(8), "-mode 3dcoffee sets that equal -method sap_pair alone");
        Assert.That(equalToProbaPlusSap, Is.EqualTo(7), "-mode 3dcoffee sets that equal proba_pair,sap_pair (7 of 8)");
    }

    /// <summary>
    /// A2-12: <c>dali_pair</c> needs DaliLite and <c>fugue_pair</c> needs fugueali; neither is installed. The program stops with exit status 1 at the first pair of two usable structures (dali_pair; pairs without two structures are
    /// aligned by proba_pair, so 4 of the 60 committed sets run through) or when it loads the method (fugue_pair, 60 of 60 sets): both raise <see cref="TCoffeeAbortedException"/> exactly where the program stops.
    /// </summary>
    [Test]
    public void DaliAndFuguePair_TheProgramStopsBecauseTheAlignersAreNotInstalled()
    {
        AssertStructureFeature("dali_pair", 60);
        AssertStructureFeature("fugue_pair", 60);
        Assert.That(Oracle().Where(c => c.Feature == "fugue_pair").All(c => c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.True);
        Assert.That(Oracle().Count(c => c.Feature == "dali_pair" && c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(50));
    }

    private static string BlastXml(Case c)
    {
        byte[] z = Convert.FromBase64String(c.Extra.GetProperty("xmlz").GetString()!);
        using var ms = new MemoryStream(z);
        using var zs = new System.IO.Compression.ZLibStream(ms, System.IO.Compression.CompressionMode.Decompress);
        using var reader = new StreamReader(zs, System.Text.Encoding.Latin1);
        return reader.ReadToEnd();
    }

    /// <summary>
    /// A2-9: <c>-method blastp_msa</c> (alone, with slow_pair, with proba_pair first, without relaxation, UPGMA): the sequences are formatted as a local BLAST database (makeblastdb) and searched against themselves
    /// with <c>blastall -p blastp -m7</c> (ncbi-blast+-legacy); the library is built from the high scoring pairs of that XML. 170 sets; the BLAST search is the external program, the XML it wrote is the input.
    /// </summary>
    [Test]
    public void BlastpMsa_EveryRowEqualsTheProgram()
    {
        var cases = Oracle().Where(c => c.Feature == "blastp_msa").ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(160), "oracle sets of 'blastp_msa'");
        var failures = new List<string>();
        foreach (var c in cases)
        {
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with { BlastXml = BlastXml(c) };
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
            if (got != c.Rows) failures.Add($"'{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} 'blastp_msa' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
    }

    /// <summary>A2-9: the library <c>t_coffee -method blastp_msa -lib_only -out_lib</c> writes (HSP residue pairs, weight = truncated percent identity of the HSP, an entry and its mirror) equals the one built from the BLAST XML.</summary>
    [Test]
    public void BlastpMsa_LibraryEqualsTheProgram()
    {
        var cases = Oracle().Where(c => c.Feature == "blastp_msa" && c.Extra.TryGetProperty("lib", out _)).ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(70), "oracle sets with a library");
        int entries = 0;
        foreach (var c in cases)
        {
            var names = Enumerable.Range(0, c.Seqs.Length).Select(i => $"q{i:D3}").ToArray();
            var got = SequenceAligner.TCoffeeBlastMsaLibraryEntries(names, c.Seqs, BlastXml(c)).Select(e => $"{e.SeqA}:{e.PosA + 1}:{e.SeqB}:{e.PosB + 1}:{e.Weight}:{e.Count}").OrderBy(x => x, StringComparer.Ordinal).ToArray();
            string lib = c.Extra.GetProperty("lib").GetString()!;
            var expected = lib.Length == 0 ? Array.Empty<string>() : lib.Split(';').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            Assert.That(got, Is.EqualTo(expected), $"library of {string.Join('|', c.Seqs)}: only ours {string.Join(',', got.Except(expected).Take(8))}; only the program's {string.Join(',', expected.Except(got).Take(8))}");
            entries += got.Length;
        }

        Assert.That(entries, Is.GreaterThan(1000), "library entries compared");
    }

    /// <summary>
    /// A2-9: <c>-method tblastx_msa</c> (nucleotide sequences, a local database of themselves, <c>blastall -p tblastx -m 7 -S1</c>; every hit of a later sequence gives an entry for each nucleotide of the query range of every HSP, reverse frames
    /// mirrored with the sequence lengths, weight = truncated HSP percent identity): rows of 101 sets (5 sets without any hit stop the program: empty library) and the libraries of the plain-method sets.
    /// </summary>
    [Test]
    public void TblastxMsa_EveryRowAndLibraryEqualsTheProgram()
    {
        var cases = Oracle().Where(c => c.Feature == "tblastx_msa").ToList();
        Assert.That(cases.Count, Is.GreaterThanOrEqualTo(100), "oracle sets of 'tblastx_msa'");
        var failures = new List<string>();
        int libraries = 0, entries = 0;
        foreach (var c in cases)
        {
            var options = TCoffeeProgramOptions.FromCommandLine(c.Flags) with { BlastXml = BlastXml(c) };
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                Assert.Throws<TCoffeeAbortedException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(c.Seqs, options).AlignedSequences);
            if (got != c.Rows) failures.Add($"'{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
            if (!c.Extra.TryGetProperty("lib", out var libProperty)) continue;
            var lib = SequenceAligner.TCoffeeTblastxLibraryEntries(c.Seqs, BlastXml(c)).Select(e => $"{e.SeqA}:{e.PosA + 1}:{e.SeqB}:{e.PosB + 1}:{e.Weight}:{e.Count}").OrderBy(x => x, StringComparer.Ordinal).ToArray();
            string text = libProperty.GetString()!;
            var expected = text.Length == 0 ? Array.Empty<string>() : text.Split(';').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            if (!lib.SequenceEqual(expected)) failures.Add($"library of {string.Join('|', c.Seqs)}: only ours {string.Join(',', lib.Except(expected).Take(5))}; only the program's {string.Join(',', expected.Except(lib).Take(5))}");
            libraries++;
            entries += lib.Length;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {cases.Count} 'tblastx_msa' sets differ from t_coffee; first: {failures.FirstOrDefault()}");
        Assert.That(libraries, Is.GreaterThanOrEqualTo(40), "libraries compared");
        Assert.That(entries, Is.GreaterThan(1000), "library entries compared");
    }

    [Test]
    public void TMalignAlignment_IsTheFirstAndThirdOfTheLastFourLines_EachWithAFinalGap()
    {
        const string output = "TM-score= 0.5\n(You should use TM-score normalized by length of the reference protein)\n\n(\":\" denotes aligned residue pairs of d < 5.0 A)\nAC-DE\n:: ::\nACGDE\n\n";
        var rows = SequenceAligner.TCoffeeTMalignAlignment(output);
        Assert.That(rows, Is.EqualTo(("AC-DE-", "ACGDE-")));
        Assert.That(SequenceAligner.TCoffeeTMalignAlignment("a\nb\n"), Is.Null);
    }

    [Test]
    public void MustangAlignment_IsTheFirstTwoFastaRecords()
    {
        var rows = SequenceAligner.TCoffeeMustangAlignment(">p0.pdb\nAC-DE\nFG\n\n>p1.pdb\nACGDE\n-G\n");
        Assert.That(rows, Is.EqualTo(("AC-DEFG", "ACGDE-G")));
        Assert.That(SequenceAligner.TCoffeeMustangAlignment(">p0.pdb\nAC\n"), Is.Null);
    }
}
