using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>B09 A5-10: the Mash text tools (<c>mash_bounds</c>, <c>mash_dist</c>, <c>mash_triangle</c>); the text is that of the installed mash 2.3 (locked byte for byte by the 660-row oracle of the library).</summary>
[TestFixture]
public class MashToolsTests
{
    private const string G1 = ">g1 first\nACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCGACCCTTTAGGGATCGATTTGACCAGGACTTAGCAATGCCGATTAGGCATCG\n";
    private const string G2 = ">g2 second\nACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCGACCCTTAAGGGATCGATTTGACCAGGACTTAGCAATGCCGATTAGGCATCG\n";

    [Test]
    public void MashBounds_FirstLines()
    {
        var r = AnalysisTools.MashBounds();
        Assert.That(r.Text, Does.Contain("\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0270708\t0.0867606\tinf"));
        // the whole text of `mash bounds -k 9 -p 0.9` (installed mash 2.3)
        Assert.That(AnalysisTools.MashBounds(9, 0.9).Text, Is.EqualTo(
            "\n" +
            "Parameters (run with -h for details):\n" +
            "   k:   9\n" +
            "   p:   0.9\n" +
            "\n" +
            "\tMash distance\n" +
            "Sketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n" +
            "100\t0.0141961\t0.0268375\t0.0501085\t0.0612639\t0.109853\t0.135775\tinf\tinf\n" +
            "500\t0.0061394\t0.0116762\t0.0179134\t0.0259746\t0.0349944\t0.0494808\t0.0657367\t0.0920923\n" +
            "1000\t0.00436992\t0.00811332\t0.0125014\t0.0188267\t0.0248309\t0.0313533\t0.0412718\t0.0603483\n" +
            "5000\t0.00191179\t0.00360949\t0.00558262\t0.007896\t0.0104209\t0.0137499\t0.017764\t0.0233147\n" +
            "10000\t0.00134621\t0.00252462\t0.00386759\t0.00546784\t0.00732082\t0.00960518\t0.0126118\t0.0157367\n" +
            "50000\t0.000596252\t0.0011235\t0.00172091\t0.00240138\t0.00321329\t0.00425859\t0.00547781\t0.00695277\n" +
            "100000\t0.000421678\t0.000789597\t0.00121226\t0.00169611\t0.00227986\t0.0029796\t0.00381786\t0.00486373\n" +
            "500000\t0.000187895\t0.000352542\t0.000539145\t0.000756281\t0.00101426\t0.00132132\t0.00169787\t0.00216651\n" +
            "1000000\t0.00013292\t0.000249078\t0.000381147\t0.000534156\t0.000714634\t0.000934215\t0.00119951\t0.00153067\n" +
            "\n" +
            "\tScreen distance\n" +
            "Sketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n" +
            "100\t0.0142682\t0.02202\t0.0342278\t0.0446972\t0.0726841\t0.100516\t0.65\t0.6\n" +
            "500\t0.00623219\t0.00953706\t0.013749\t0.0174933\t0.0240449\t0.0333684\t0.0505157\t0.0585452\n" +
            "1000\t0.00430608\t0.00676112\t0.00920323\t0.0121109\t0.0157883\t0.0226841\t0.0276793\t0.0449527\n" +
            "5000\t0.00189919\t0.00293418\t0.00408045\t0.00531664\t0.00700955\t0.00925089\t0.0124629\t0.0168392\n" +
            "10000\t0.00133891\t0.00208953\t0.00285862\t0.00379156\t0.00488052\t0.00650786\t0.00849214\t0.0112439\n" +
            "50000\t0.000595956\t0.000928005\t0.00127038\t0.00166802\t0.00216086\t0.00280777\t0.00375355\t0.00504655\n" +
            "100000\t0.000421027\t0.000655067\t0.000899629\t0.00117763\t0.00152698\t0.00199551\t0.00262918\t0.00350552\n" +
            "500000\t0.0001883\t0.000292696\t0.000400747\t0.000525721\t0.000682474\t0.000886472\t0.00115863\t0.00154862\n" +
            "1000000\t0.000133124\t0.000207069\t0.000283342\t0.000372047\t0.00048257\t0.000625956\t0.000819775\t0.00109079\n" +
            "\n"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashBounds(33));
    }

    [Test]
    public void MashDist_And_Triangle_AgreeOnTheDistance()
    {
        var d = AnalysisTools.MashDist([G1], [G2], ["g1.fa"], ["g2.fa"], k: 9, sketchSize: 20);
        var t = AnalysisTools.MashTriangle([G1, G2], ["g1.fa", "g2.fa"], k: 9, sketchSize: 20);
        string distance = d.Text.Split('\t')[2];
        Assert.That(d.Text, Does.StartWith("g1.fa\tg2.fa\t"));
        Assert.That(t.Text, Is.EqualTo($"\t2\ng1.fa\ng2.fa\t{distance}\n"));
        Assert.That(AnalysisTools.MashDist([G1], [G2], k: 9, sketchSize: 20, maxDistance: 0.0).Text, Is.Empty);        // -d 0: only identical pairs pass
        Assert.That(AnalysisTools.MashDist([G1], [G2], k: 9, sketchSize: 20, table: true).Text, Does.StartWith("#query\treference1\nquery1\t"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashDist([], [G2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashTriangle([G1, G2], seed: -1));
    }

    [Test]
    public void MashDist_And_Triangle_RefuseAFileWhoseRecordsAreAllShorterThanK_LikeTheProgram()
    {
        // A6-6: `mash dist -k 11 -s 100 g.fa short.fa` -> "WARNING: All fasta records in input files were shorter than the k-mer size (11).", exit status 1; -i skips the record instead
        const string Short = ">s\nACG\n";
        var ex = Assert.Throws<ArgumentException>(() => AnalysisTools.MashDist([G1], [Short], k: 11, sketchSize: 100))!;
        Assert.That(ex.Message, Does.StartWith("WARNING: All fasta records in input files were shorter than the k-mer size (11)."));
        Assert.That(ex.ParamName, Is.EqualTo("queries"));
        Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.MashDist([Short], [G1], k: 11, sketchSize: 100))!.ParamName, Is.EqualTo("references"));
        Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.MashTriangle([G1, Short], k: 11, sketchSize: 100))!.Message, Does.StartWith("WARNING: All fasta records in input files were shorter than the k-mer size (11)."));
        Assert.That(AnalysisTools.MashDist([G1], [Short], k: 11, sketchSize: 100, individual: true).Text, Is.Empty);
        Assert.That(AnalysisTools.MashTriangle([G1, Short], k: 11, sketchSize: 100, individual: true).Text, Is.EqualTo("\t1\ng1\n"));
    }

    [Test]
    public void MashDist_And_Triangle_TakeACancellationToken()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => AnalysisTools.MashDist([G1], [G2], k: 9, sketchSize: 20, cancellationToken: cts.Token));
        Assert.Throws<OperationCanceledException>(() => AnalysisTools.MashTriangle([G1, G2], k: 9, sketchSize: 20, cancellationToken: cts.Token));
    }

    private const string SketchInput = ">s1 first seq\nACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGTCGATCGATGCTAGCTAGCATCGATCG\n>s2\nTTGACCAGTAGCTAGCTAGGGATCGATTTAGCGCGATATATCGCGCTAGCTAGTTTGGGCCCAAATTTAGC\n";

    [Test]
    public void MashSketch_BloomFilter_AsTheProgram()
    {
        // mash 2.3: mash sketch -k 8 -s 6 -b 100 -o tb t.fa (A5-10b2, F197): the Bloom filter of 100 bytes (one hash function over 800 bits), the counts of the read-set sketch in robin_hood order
        var sketch = AnalysisTools.MashSketch([SketchInput], ["t.fa"], k: 8, sketchSize: 6, bloomFilter: "100");
        Assert.That(sketch.Msh, Is.EqualTo("msh:" + "AAAAACEAAAAAAAAAAwAEAAgAAAAAAAAABgAAAAEAAAAAAAAAAAAAAAwAAAAAAAEAYAAAAAAAAQBlAAAAKgAAAAAAAAAAAAAAAQAAAE8AAAAEAAAAAgAHAAAAAAAAAAAAUgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEQAAACoAAAARAAAA4gAAAB0AAAA0AAAAAAAAAAAAAAAhAAAANAAAAHQuZmEAAAAAWzIgc2Vxc10gczEgZmlyc3Qgc2VxIFsuLi5dAAAAAAAfQEMBkxLfA1WzHAUqyT4Q2GKoEYaxkBICAAAAAgAAAAIAAAACAAAAAwAAAAIAAAABAAAABwAAAAAAAAADAAAAQUNHVAAAAAA="));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashSketch([SketchInput], ["t.fa"], k: 8, sketchSize: 6, bloomFilter: "5X"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashSketch([SketchInput], ["t.fa"], k: 8, sketchSize: 6, minCopies: 2, bloomFilter: "100"));      // -m cannot be used with -b
    }

    [Test]
    public void MashSketch_Info_Paste_AsTheProgram()
    {
        // mash 2.3: mash sketch -k 8 -s 6 -o t t.fa (A5-10b, F193): the file is byte for byte the one the program wrote; the msh: text is accepted by the other tools
        var sketch = AnalysisTools.MashSketch([SketchInput], ["t.fa"], k: 8, sketchSize: 6);
        Assert.That(sketch.Msh, Is.EqualTo("msh:" + "AAAAAB4AAAAAAAAAAwAEAAgAAAAAAAAABgAAAAEAAAAAAAAAAAAAAAwAAAAAAAEAVAAAAAAAAQBZAAAAKgAAAAAAAAAAAAAAAQAAAE8AAAAEAAAAAgAHAAAAAAAAAAAAjAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEQAAACoAAAARAAAA4gAAAB0AAAA0AAAAAAAAAAAAAAAAAAAAAAAAAHQuZmEAAAAAWzIgc2Vxc10gczEgZmlyc3Qgc2VxIFsuLi5dAAAAAACYOngAH0BDAZMS3wNVsxwF9IB0CboKGAsBAAAABwAAAAAAAAADAAAAQUNHVAAAAAA="));
        Assert.That(sketch.Sketches, Is.EqualTo(1));
        Assert.That(AnalysisTools.MashInfo(sketch.Msh).Text, Is.EqualTo("Header:\n  Hash function (seed):          MurmurHash3_x64_128 (42)\n  K-mer size:                    8 (32-bit hashes)\n  Alphabet:                      ACGT (canonical)\n  Target min-hashes per sketch:  6\n  Sketches:                      1\n\nSketches:\n  [Hashes]  [Length]  [ID]  [Comment]\n\n  6         140       t.fa  [2 seqs] s1 first seq [...]\n\n"));
        Assert.That(AnalysisTools.MashInfo(sketch.Msh, tabular: true).Text, Is.EqualTo("#Hashes\tLength\tID\tComment\n6\t140\tt.fa\t[2 seqs] s1 first seq [...]\n"));
        Assert.That(AnalysisTools.MashPaste([sketch.Msh, sketch.Msh]).Sketches, Is.EqualTo(2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashInfo(sketch.Msh, headerOnly: true, tabular: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashInfo("msh:AAAA"));
        // a sketch file is an input of mash_dist (its parameters are inherited) and of mash_triangle
        string dist = AnalysisTools.MashDist([sketch.Msh], [SketchInput], queryNames: ["q"]).Text;
        Assert.That(dist, Is.EqualTo("t.fa\tq\t0\t1.47546e-18\t6/6\n"));
    }
}
