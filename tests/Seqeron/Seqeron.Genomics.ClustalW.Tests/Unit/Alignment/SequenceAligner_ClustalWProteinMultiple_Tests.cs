// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10c part b, B09 F63): <see cref="ClustalWAligner.MultipleAlignClustalW(IEnumerable{string}, ClustalWProteinOptions, MatrixScoring?)"/>,
/// the step-by-step port of ClustalW 2.1's multiple-alignment stage for PROTEIN: the matrix series chosen by the percent identity of the two profiles
/// (GONNET / BLOSUM / PAM / ID), the matrix- and length-scaled gap opening, the residue-specific (Pascarella &amp; Argos), hydrophilic-run and variable-region gap
/// penalties and the gap-separation distance, on top of the guide tree / weights / delayed sequences / Myers–Miller profile alignment shared with the DNA stage.
///
/// Oracle (TestData/Alignment/clustalw_protein_msa_oracle.tsv, scripts/alignment/generate_clustalw_protein_msa_oracle.py): the installed <c>clustalw</c> 2.1
/// <c>-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT</c> with the option variants -MATRIX, -PWMATRIX, -GAPOPEN, -GAPEXT, -PWGAPOPEN, -PWGAPEXT, -MAXDIV, -NOWEIGHTS,
/// -NEGATIVE, -QUICKTREE (-KTUPLE, -WINDOW, -TOPDIAGS, -PAIRGAP), -NOPGAP, -NOHGAP, -NOVGAP, -GAPDIST, -HGAPRESIDUES, -ENDGAPS on 1 000 protein sets (2-25 sequences up to
/// ~320 residues; mutated families of 2-85 % divergence, hydrophilic-run material, low-complexity, repeats, duplicates, tiny sets, B/Z/X/O/U residues). The binary
/// computes its own pairwise distances, guide tree, weights and profile alignments: every output row must be identical, except for the two sets whose only difference
/// is the uninitialised <c>gdist</c> member of ClustalW (see <see cref="KnownUninitialisedGapDistanceSets"/>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWProteinMultiple_Tests
{
    private sealed record Row(string Id, string[] Seqs, string[] Rows, string Options);

    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
        ReadTsv("clustalw_protein_msa_oracle.tsv").Select(f => new Row(f[0], f[1].Split(';'), f[2].Split(';'), f.Length > 3 ? f[3] : "")).ToList());

    /// <summary>
    /// ProfileBase::calcGapCoeff reads the member <c>gdist</c> without ever assigning it when -NOVGAP is set and a profile holds exactly two sequences of more than 60 %
    /// identity (undefined behaviour of the C++: the recycled heap chunk of the profile object still holds whatever was in it - found to be 0, the user's -GAPDIST, or a
    /// stale residue code such as 5 or 20). <see cref="ClustalWProteinOptions.UninitialisedGapDistances"/> emulates the common outcome; these are the sets of the oracle
    /// (id -> gdist of the 1st, 2nd, ... such profile) where the binary's heap held something else. The values were found by scanning gdist per profile against the
    /// binary's <c>-DEBUG=5</c> profile dumps (opening / extension penalties and every profile column); with them the rows equal the binary on every row. Proof of the undefined behaviour (valgrind 3.22.0 on the
    /// installed clustalw 2.1+lgpl-7, Evidence F143): 87 of the 189 <c>-NOVGAP</c> sets of the oracle report "Conditional jump or move depends on uninitialised value(s)", created by a heap allocation (operator new);
    /// 40 sets without <c>-NOVGAP</c> report none. The common outcome (0) is emulated by default; these two sets are the exact list of those where the recycled chunk held something else
    /// (<see cref="MultipleAlignClustalW_Protein_UninitialisedGdist_OnlyTheDocumentedSetsNeedAnOverride"/> asserts that exactly these two differ without the override).
    /// </summary>
    private static readonly Dictionary<string, int[]> KnownUninitialisedGapDistanceSets = new()
    {
        ["368"] = [0, 4, 0, 0],
        ["915"] = [5, 0, 20],
    };

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static ClustalWProfileMatrix ProfileMatrix(string name) => name switch
    {
        "GONNET" => ClustalWProfileMatrix.Gonnet,
        "BLOSUM" => ClustalWProfileMatrix.Blosum,
        "PAM" => ClustalWProfileMatrix.Pam,
        "ID" => ClustalWProfileMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    private static ClustalWPairwiseMatrix PairwiseMatrix(string name) => name switch
    {
        "GONNET" => ClustalWPairwiseMatrix.Gonnet,
        "BLOSUM" => ClustalWPairwiseMatrix.Blosum,
        "PAM" => ClustalWPairwiseMatrix.Pam,
        "ID" => ClustalWPairwiseMatrix.Identity,
        _ => throw new ArgumentException(name),
    };

    /// <summary>The clustalw command line semantics of the oracle's option string (PROTEIN: -KTUPLE sets the pair gap to ktuple + 3, -PAIRGAP overrides).</summary>
    internal static ClustalWProteinOptions Parse(string s)
    {
        var o = new ClustalWProteinOptions();
        if (s.Length == 0)
            return o;
        foreach (var kv in s.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "matrix": o = o with { Matrix = ProfileMatrix(p[1]) }; break;
                case "pwmatrix": o = o with { Pairwise = o.Pairwise with { Matrix = PairwiseMatrix(p[1]) } }; break;
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "pwgapopen": o = o with { Pairwise = o.Pairwise with { GapOpen = D(p[1]) } }; break;
                case "pwgapext": o = o with { Pairwise = o.Pairwise with { GapExtend = D(p[1]) } }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "noweights": o = o with { NoWeights = true }; break;
                case "negative": o = o with { NegativeMatrix = true }; break;
                case "nopgap": o = o with { NoPreferredGapPenalties = true }; break;
                case "nohgap": o = o with { NoHydrophilicGapPenalties = true }; break;
                case "novgap": o = o with { VariableRegionGapPenalties = true }; break;
                case "gapdist": o = o with { GapSeparationDistance = int.Parse(p[1]) }; break;
                case "hgapresidues": o = o with { HydrophilicResidues = p[1] }; break;
                case "endgaps": break;   // the command line flag only writes useEndGaps = false: no effect
                case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                case "ktuple": o = o with { Pairwise = o.Pairwise with { Quick = o.Pairwise.Quick with { KTuple = int.Parse(p[1]), PairGap = int.Parse(p[1]) + 3 } } }; break;
                case "window": o = o with { Pairwise = o.Pairwise with { Quick = o.Pairwise.Quick with { Window = int.Parse(p[1]) } } }; break;
                case "topdiags": o = o with { Pairwise = o.Pairwise with { Quick = o.Pairwise.Quick with { TopDiagonals = int.Parse(p[1]) } } }; break;
                case "pairgap": o = o with { Pairwise = o.Pairwise with { Quick = o.Pairwise.Quick with { PairGap = int.Parse(p[1]) } } }; break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static string[] Align(IEnumerable<string> seqs, ClustalWProteinOptions? o = null) =>
        ClustalWAligner.MultipleAlignClustalW(seqs, o ?? new ClustalWProteinOptions()).AlignedSequences.ToArray();

    private static string[] AlignOracleRow(Row r) =>
        Align(r.Seqs, KnownUninitialisedGapDistanceSets.TryGetValue(r.Id, out var forced) ? Parse(r.Options) with { UninitialisedGapDistances = forced } : Parse(r.Options));

    [Test]
    public void MultipleAlignClustalW_Protein_Rows_EqualClustalW21_OnAllOracleSets()
    {
        var rows = Oracle.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(1000));
        int bad = 0;
        string? first = null;
        foreach (var r in rows)
        {
            string[] got = AlignOracleRow(r);
            if (!got.SequenceEqual(r.Rows) && ++bad == 1)
                first = $"set {r.Id} ({r.Options}): expected {string.Join('/', r.Rows)}, got {string.Join('/', got)}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>
    /// The default emulation of the unassigned <c>gdist</c> reproduces the binary on every set except the documented ones, and on those it really differs (so the list
    /// cannot go stale): 998 of the 1 000 oracle sets - 2 of them with -NOVGAP - are reproduced without any override.
    /// </summary>
    [Test]
    public void MultipleAlignClustalW_Protein_UninitialisedGdist_OnlyTheDocumentedSetsNeedAnOverride()
    {
        var differing = Oracle.Value.Where(r => !Align(r.Seqs, Parse(r.Options)).SequenceEqual(r.Rows)).Select(r => r.Id).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(differing, Is.EquivalentTo(KnownUninitialisedGapDistanceSets.Keys));
            foreach (string id in differing)
                Assert.That(Oracle.Value.First(r => r.Id == id).Options.Split(',').Contains("novgap"), Is.True, $"set {id} uses -NOVGAP");
        });
    }

    [Test]
    public void MultipleAlignClustalW_Protein_OracleCoversEveryOptionAndTheSpecialPaths()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            foreach (string key in new[] { "matrix", "pwmatrix", "gapopen", "gapext", "maxdiv", "noweights", "negative", "quick", "ktuple", "window", "topdiags", "pairgap", "pwgapopen", "pwgapext", "nopgap", "nohgap", "novgap", "gapdist", "hgapresidues", "endgaps" })
                Assert.That(rows.Count(r => r.Options.Split(',').Any(o => o.Split('=')[0] == key)), Is.GreaterThan(20), key);
            foreach (string m in new[] { "GONNET", "BLOSUM", "PAM", "ID" })
                Assert.That(rows.Count(r => r.Options.Contains("matrix=" + m)), Is.GreaterThan(20), "matrix " + m);
            Assert.That(rows.Count(r => r.Options.Length == 0), Is.GreaterThan(200), "default runs");
            Assert.That(rows.Count(r => r.Seqs.Length == 2), Is.GreaterThan(30), "two sequences (rooted two-leaf tree)");
            Assert.That(rows.Count(r => r.Seqs.Length >= 14), Is.GreaterThan(30), "large sets");
            Assert.That(rows.Count(r => r.Seqs.Any(s => s.Length >= 100)), Is.GreaterThan(150), "profiles of >= 100 residues (GONNET250 / 120, 350 / 160 split)");
            Assert.That(rows.Count(r => r.Seqs.Any(s => s.Any(c => "BZXOU".Contains(c)))), Is.GreaterThan(30), "ambiguity codes B Z X and O U");
            Assert.That(rows.Count(r => r.Seqs.Distinct().Count() < r.Seqs.Length), Is.GreaterThan(20), "duplicate sequences (tied distances)");
        });
    }

    /// <summary>
    /// The oracle discriminates every protein mechanism: for each switch there are many sets whose rows change when the mechanism is turned off or varied, so equality
    /// with clustalw is not an accident of easy inputs.
    /// </summary>
    [Test]
    public void MultipleAlignClustalW_Protein_EachMechanism_ChangesRowsOnManyOracleSets()
    {
        var sets = Oracle.Value.Where(r => r.Options.Length == 0 && r.Seqs.Length >= 3).Take(120).ToList();
        int noPref = 0, noHyd = 0, gapDist = 0, vgap = 0, blosum = 0, pam = 0, id = 0, negative = 0, noWeights = 0, delayed = 0, quick = 0, gap = 0, pwMatrix = 0;
        foreach (var r in sets)
        {
            var baseRows = Align(r.Seqs);
            var d = new ClustalWProteinOptions();
            if (!Align(r.Seqs, d with { NoPreferredGapPenalties = true }).SequenceEqual(baseRows)) noPref++;
            if (!Align(r.Seqs, d with { NoHydrophilicGapPenalties = true }).SequenceEqual(baseRows)) noHyd++;
            if (!Align(r.Seqs, d with { GapSeparationDistance = 0 }).SequenceEqual(baseRows)) gapDist++;
            if (!Align(r.Seqs, d with { VariableRegionGapPenalties = true }).SequenceEqual(baseRows)) vgap++;
            if (!Align(r.Seqs, d with { Matrix = ClustalWProfileMatrix.Blosum }).SequenceEqual(baseRows)) blosum++;
            if (!Align(r.Seqs, d with { Matrix = ClustalWProfileMatrix.Pam }).SequenceEqual(baseRows)) pam++;
            if (!Align(r.Seqs, d with { Matrix = ClustalWProfileMatrix.Identity }).SequenceEqual(baseRows)) id++;
            if (!Align(r.Seqs, d with { NegativeMatrix = true }).SequenceEqual(baseRows)) negative++;
            if (!Align(r.Seqs, d with { NoWeights = true }).SequenceEqual(baseRows)) noWeights++;
            if (!Align(r.Seqs, d with { DivergenceCutoff = 0 }).SequenceEqual(baseRows)) delayed++;
            if (!Align(r.Seqs, d with { PairwiseMode = ClustalWPairwiseMode.Quick }).SequenceEqual(baseRows)) quick++;
            if (!Align(r.Seqs, d with { GapOpen = 3, GapExtend = 1 }).SequenceEqual(baseRows)) gap++;
            if (!Align(r.Seqs, d with { Pairwise = d.Pairwise with { Matrix = ClustalWPairwiseMatrix.Blosum } }).SequenceEqual(baseRows)) pwMatrix++;
        }

        Assert.Multiple(() =>
        {
            Assert.That(noPref, Is.GreaterThan(10), "residue-specific gap penalties");
            Assert.That(noHyd, Is.GreaterThan(5), "hydrophilic gap penalties");
            Assert.That(gapDist, Is.GreaterThan(5), "gap separation distance");
            Assert.That(vgap, Is.GreaterThan(3), "variable-region gap penalties");
            Assert.That(blosum, Is.GreaterThan(10), "BLOSUM series");
            Assert.That(pam, Is.GreaterThan(10), "PAM series");
            Assert.That(id, Is.GreaterThan(10), "identity matrix");
            Assert.That(negative, Is.GreaterThan(10), "negative matrix");
            Assert.That(noWeights, Is.GreaterThan(5), "sequence weights");
            Assert.That(delayed, Is.GreaterThan(5), "delayed divergent sequences");
            Assert.That(quick, Is.GreaterThan(3), "fast pairwise stage");
            Assert.That(gap, Is.GreaterThan(20), "gap penalties");
            Assert.That(pwMatrix, Is.GreaterThan(3), "pairwise matrix");
        });
    }

    /// <summary>The Debian clustalw keeps at most MAXHYDRESIDUES = 9 letters of -HGAPRESIDUES.</summary>
    [Test]
    public void MultipleAlignClustalW_Protein_HydrophilicResidues_AreTruncatedToNineLetters()
    {
        var sets = Oracle.Value.Where(r => r.Options.Length == 0 && r.Seqs.Length >= 3).Take(80).ToList();
        int differsFromDefault = 0;
        foreach (var r in sets)
        {
            var twelve = new ClustalWProteinOptions { HydrophilicResidues = "ACDEFGHIKLMN" };
            var nine = new ClustalWProteinOptions { HydrophilicResidues = "ACDEFGHIK" };
            Assert.That(Align(r.Seqs, twelve), Is.EqualTo(Align(r.Seqs, nine)), $"set {r.Id}");
            if (!Align(r.Seqs, nine).SequenceEqual(Align(r.Seqs)))
                differsFromDefault++;
        }

        Assert.That(differsFromDefault, Is.GreaterThan(5));
        // an empty or non-letter list keeps the default GPSNDQEKR
        var s = sets[0].Seqs;
        Assert.That(Align(s, new ClustalWProteinOptions { HydrophilicResidues = "" }), Is.EqualTo(Align(s)));
        Assert.That(Align(s, new ClustalWProteinOptions { HydrophilicResidues = "gpsndqekr" }), Is.EqualTo(Align(s)));
    }

    [Test]
    public void MultipleAlignClustalW_Protein_Invariants_RowsSpellInputs_EqualLength_InputOrder()
    {
        foreach (var r in Oracle.Value.Take(300))
        {
            var rows = AlignOracleRow(r);
            Assert.That(rows.Length, Is.EqualTo(r.Seqs.Length));
            Assert.That(rows.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
            for (int i = 0; i < rows.Length; i++)
                Assert.That(rows[i].Replace("-", ""), Is.EqualTo(r.Seqs[i]), $"set {r.Id} row {i}");
        }
    }

    [Test]
    public void MultipleAlignClustalW_Protein_ProteinSequenceOverload_EqualsStringOverload_AndLowerCase()
    {
        var set = Oracle.Value.First(r => r.Options.Length == 0 && r.Seqs.Length >= 5 && r.Seqs.All(s => s.All(c => "ACDEFGHIKLMNPQRSTVWY".Contains(c))));
        var viaString = ClustalWAligner.MultipleAlignClustalW(set.Seqs, new ClustalWProteinOptions());
        var viaProtein = ClustalWAligner.MultipleAlignClustalW(set.Seqs.Select(s => new ProteinSequence(s)));
        var lower = ClustalWAligner.MultipleAlignClustalW(set.Seqs.Select(s => s.ToLowerInvariant()), new ClustalWProteinOptions());
        var scored = ClustalWAligner.MultipleAlignClustalW(set.Seqs, new ClustalWProteinOptions(), MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5));
        Assert.Multiple(() =>
        {
            Assert.That(viaString.AlignedSequences, Is.EqualTo(set.Rows));
            Assert.That(viaProtein.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(lower.AlignedSequences, Is.EqualTo(viaString.AlignedSequences));
            Assert.That(viaString.ExactTotalScore, Is.EqualTo(SequenceAligner.MatrixSumOfPairsScore(viaString.AlignedSequences, MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2), MsaGapModel.Affine)).Within(1e-9));
            Assert.That(scored.ExactTotalScore, Is.EqualTo(SequenceAligner.MatrixSumOfPairsScore(scored.AlignedSequences, MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5), MsaGapModel.Affine)).Within(1e-9));
            Assert.That(viaString.Consensus.Length, Is.EqualTo(viaString.AlignedSequences[0].Length));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Protein_Contracts()
    {
        string[] two = ["MKVLAAGIVGLLLA", "MKVLSAGIVALLA"];
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(Array.Empty<string>(), new ClustalWProteinOptions()), "A3-4: the program refuses an empty input");
            var one = ClustalWAligner.MultipleAlignClustalW(new[] { "mkvl" }, new ClustalWProteinOptions());
            Assert.That(one.AlignedSequences, Is.EqualTo(new[] { "MKVL" }));
            Assert.That(one.TotalScore, Is.Zero);
            Assert.Throws<ArgumentNullException>(() => ClustalWAligner.MultipleAlignClustalW((IEnumerable<string>)null!, new ClustalWProteinOptions()));
            Assert.Throws<ArgumentNullException>(() => ClustalWAligner.MultipleAlignClustalW(two, (ClustalWProteinOptions)null!));
            Assert.Throws<ArgumentNullException>(() => ClustalWAligner.MultipleAlignClustalW((IEnumerable<ProteinSequence>)null!));
            Assert.Throws<ArgumentNullException>(() => ClustalWAligner.MultipleAlignClustalW(new string[] { "MKVL", null! }, new ClustalWProteinOptions()));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(new[] { "MKVL", "MKJL" }, new ClustalWProteinOptions()), "J is dropped by ClustalW");
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(new[] { "MKVL", "MK*L" }, new ClustalWProteinOptions()), "'*' is dropped by ClustalW");
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(new[] { "MKVL", "" }, new ClustalWProteinOptions()));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { GapOpen = -1 }));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { GapExtend = double.NaN }));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { Pairwise = new ClustalWProteinPairwiseOptions { GapOpen = -2 } }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { Matrix = (ClustalWProfileMatrix)9 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { PairwiseMode = (ClustalWPairwiseMode)9 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { GapSeparationDistance = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ClustalWAligner.MultipleAlignClustalW(two, new ClustalWProteinOptions { Pairwise = new ClustalWProteinPairwiseOptions { Quick = new ClustalWQuickOptions { KTuple = 3 } } }));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Protein_DefaultOptions_AreClustalWProteinDefaults()
    {
        var o = new ClustalWProteinOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.Matrix, Is.EqualTo(ClustalWProfileMatrix.Gonnet));
            Assert.That(o.GapOpen, Is.EqualTo(10.0));
            Assert.That(o.GapExtend, Is.EqualTo(0.2));
            Assert.That(o.Pairwise.Matrix, Is.EqualTo(ClustalWPairwiseMatrix.Gonnet));
            Assert.That(o.Pairwise.GapOpen, Is.EqualTo(10.0));
            Assert.That(o.Pairwise.GapExtend, Is.EqualTo(0.1));
            Assert.That(o.DivergenceCutoff, Is.EqualTo(30));
            Assert.That(o.GapSeparationDistance, Is.EqualTo(4));
            Assert.That(o.HydrophilicResidues, Is.EqualTo("GPSNDQEKR"));
            Assert.That(o.NoWeights || o.NegativeMatrix || o.NoPreferredGapPenalties || o.NoHydrophilicGapPenalties || o.VariableRegionGapPenalties, Is.False);
            Assert.That(o.PairwiseMode, Is.EqualTo(ClustalWPairwiseMode.Full));
            Assert.That(o.UninitialisedGapDistances, Is.Null);
        });
    }

    /// <summary>The DNA overload and the DNA oracle are untouched: a DNA set aligned as before is still the binary's DNA result.</summary>
    [Test]
    public void MultipleAlignClustalW_DnaStringOverload_StillSelectsTheDnaStage()
    {
        string[] seqs = ["ACGTACGTTTGACCAGTAGGCATCGAT", "ACGTACGATTGACCAGTAGCCATCGAT", "GGGGGGAAAATTTTCCCC"];
        var dna = ClustalWAligner.MultipleAlignClustalW(seqs).AlignedSequences;
        var asProtein = ClustalWAligner.MultipleAlignClustalW(seqs, new ClustalWProteinOptions()).AlignedSequences;
        Assert.That(dna.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
        Assert.That(asProtein.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
    }

    /// <summary>
    /// The forced gdist of <see cref="ClustalWProteinOptions.UninitialisedGapDistances"/> only matters on -NOVGAP runs: without that option the list is ignored.
    /// </summary>
    [Test]
    public void MultipleAlignClustalW_Protein_UninitialisedGapDistances_AreIgnoredWithoutVariableRegionPenalties()
    {
        foreach (var r in Oracle.Value.Where(r => !r.Options.Contains("novgap")).Take(60))
        {
            var o = Parse(r.Options);
            Assert.That(Align(r.Seqs, o with { UninitialisedGapDistances = [7, 7, 7, 7, 7] }), Is.EqualTo(Align(r.Seqs, o)), $"set {r.Id}");
        }
    }

    #region matrix series tables

    private static readonly Lazy<List<string[]>> SeriesReference = new(() => ReadTsv("clustalw_protein_series_reference.tsv").ToList());

    private static SubstitutionMatrix Table(string name) => name switch
    {
        "BLOSUM40" => ClustalWMatrices.ClustalWBlosum40,
        "BLOSUM45" => ClustalWMatrices.ClustalWBlosum45,
        "BLOSUM80" => ClustalWMatrices.ClustalWBlosum80,
        "BLOSUM62" => SubstitutionMatrix.Blosum62,
        "PAM20" => ClustalWMatrices.ClustalWPam20,
        "PAM60" => ClustalWMatrices.ClustalWPam60,
        "PAM120" => ClustalWMatrices.ClustalWPam120,
        "PAM350" => ClustalWMatrices.ClustalWPam350,
        _ => throw new ArgumentException(name),
    };

    /// <summary>Every cell of the series tables equals the EMBOSS data file (EBLOSUM40/45/62/80, EPAM20/60/120/350) on all 23 x 23 residue pairs.</summary>
    [Test]
    public void SeriesMatrices_EqualEmbossDataFiles_OnEveryCell()
    {
        var rows = SeriesReference.Value;
        Assert.That(rows.Count, Is.EqualTo(8 * 23 * 23));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            int expected = int.Parse(f[3], CultureInfo.InvariantCulture);
            int got = Table(f[0]).GetScore(f[1][0], f[2][0]);
            if (got != expected && ++bad == 1)
                first = $"{f[0]} {f[1]}/{f[2]}: expected {expected}, got {got}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    [Test]
    public void SeriesMatrices_GonnetTables_AreSymmetricAndOnTheClustalWAlphabet()
    {
        foreach (var m in new[]
        {
            ClustalWMatrices.ClustalWGonnet40, ClustalWMatrices.ClustalWGonnet80, ClustalWMatrices.ClustalWGonnet120, ClustalWMatrices.ClustalWGonnet160,
            ClustalWMatrices.ClustalWGonnet250, ClustalWMatrices.ClustalWGonnet350,
        })
        {
            Assert.That(m.Alphabet, Is.EqualTo(ClustalWMatrices.ClustalWAminoAcidOrder));
            foreach (char a in m.Alphabet)
                foreach (char b in m.Alphabet)
                    Assert.That(m.GetScore(a, b), Is.EqualTo(m.GetScore(b, a)), m.Name);
        }

        // spot values read from matrices.h (tenth-bit tables): gon40mt A-A = 92 and C-C = 163; gon250mt A-A = 24 and W-W = 142
        Assert.Multiple(() =>
        {
            Assert.That(ClustalWMatrices.ClustalWGonnet40.GetScore('A', 'A'), Is.EqualTo(92));
            Assert.That(ClustalWMatrices.ClustalWGonnet40.GetScore('C', 'C'), Is.EqualTo(163));
            Assert.That(ClustalWMatrices.ClustalWGonnet250.GetScore('A', 'A'), Is.EqualTo(24));
            Assert.That(ClustalWMatrices.ClustalWGonnet250.GetScore('W', 'W'), Is.EqualTo(142));
        });
    }

    #endregion
}
