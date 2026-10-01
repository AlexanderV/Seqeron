// MOTIF-CONS-001 — EMBOSS cons residue-type Auto and ragged-row padding (B05 audit group C)
// Evidence: docs/Evidence/MOTIF-CONS-001-Evidence.md
// Source (opened): EMBOSS 6.6.0 (Ubuntu emboss_6.6.0+dfsg.orig.tar.xz) emboss/cons.c,
// nucleus/embcons.c, ajax/core/ajseq.c (ajSeqsetFill, ajSeqsetIsNuc, ajSeqsetIsProt,
// ajSeqIsNuc), ajax/core/ajseqtype.c (ajSeqType, ajSeqSetNuc x/X → n/N, ajSeqTypeGapnucS,
// seqCharNucPure/NucAmbig/Gap, "gapany" ? → X), ajax/core/ajseqread.c (ajSeqsetFromList: set
// type = first sequence's type, Len = longest), ajax/acd/ajacd.c (aligned seqsets → ajSeqsetFill).
// Every expected string is the verbatim output of the EMBOSS 6.6.0 `cons` binary (Ubuntu noble
// package emboss, /usr/lib/emboss/cons) run WITHOUT -snucleotide/-sprotein on the rows written
// to a.fa, for the options in the description (cross-check: 3,000 cases, 0 mismatches).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_EmbossConsensusAutoPad_Tests
{
    private static readonly TestCaseData[] AutoCases =
    {
        new TestCaseData(new[] { "H-UVTGMSHKKY?R", "KNGTG?G~", "?WBDDWAYYHHU", "sryg", "RBNW", ".wta.stasysyvsgr", "C", "DH~B", ".NVA" }, 1.5f, 8, null, "NNNNNNNNNNNNNNNN").SetName("EmbossAuto_000").SetDescription("cons -sequence a.fa -plurality 1.5 -identity 8"),
        new TestCaseData(new[] { "~AAC.", "-", "YFQ", "G", "-CA.AG", "~T", "T" }, 1.0f, 0, null, "nnacag").SetName("EmbossAuto_001").SetDescription("cons -sequence a.fa -plurality 1"),
        new TestCaseData(new[] { ".~SHMBNR", ".YMHZLD-", "KFHDYH*K", "DGLVEAA~", "iwazbnze", "ngrye.rg", "*-RC*PDM", ".C~HKIZW" }, 2.0f, 0, null, "nnnnnnnn").SetName("EmbossAuto_002").SetDescription("cons -sequence a.fa -plurality 2"),
        new TestCaseData(new[] { "MAPNASQPPDEWQE", "KVHRPWNLCAMMRLVQ", "DVMM" }, null, 0, 0.0f, "XVXXxxXxxxxxQxXX").SetName("EmbossAuto_003").SetDescription("cons -sequence a.fa -setcase 0"),
        new TestCaseData(new[] { "tdarddnsedyedyskqhfceasriv", "MYRWINATQSEKTCDQGVIYVGQDLM", "CIVSHLNNGETCDSMCIKAHTLIVAQ" }, null, 0, null, "xxxxxNNSEDxKDxxKxxxYxxxxIM").SetName("EmbossAuto_004").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "E-NF", "CATT", "GML-", "FG.V", "FELP", "T~CG", "ERQH", ".A~~", "t~m." }, 3.0f, 9, null, "XXXX").SetName("EmbossAuto_005").SetDescription("cons -sequence a.fa -plurality 3 -identity 9"),
        new TestCaseData(new[] { "TXMVRUXHDYN", "WNIKYTAWWAF" }, 0.0f, 0, 2.5f, "tnmvrunhdyn").SetName("EmbossAuto_006").SetDescription("cons -sequence a.fa -plurality 0 -setcase 2.5"),
        new TestCaseData(new[] { "WBI", "B*ZISMGNNDQB", "CABSPPPSBXDZFLBQN", "kpwrnrsvqadv", "IMDHSGSIDKMSKEX", "RSNBTQVRQAWFGZBP", "MCBRNEF", "DCDKFBWPNMTP" }, 2.0f, 6, null, "XXXXXXXXXXXXXXXXX").SetName("EmbossAuto_007").SetDescription("cons -sequence a.fa -plurality 2 -identity 6"),
        new TestCaseData(new[] { "~.TTG~", "~tta~g", ".T-AG-" }, null, 3, null, "nNNNNN").SetName("EmbossAuto_008").SetDescription("cons -sequence a.fa -identity 3"),
        new TestCaseData(new[] { "TCTG-AATCA-GTAG", "ATTA-GATTT.AACG", "GGTCTATCTTAATGC", "gatcaggacattccg" }, null, 0, null, "gnTcnaatcanatcG").SetName("EmbossAuto_009").SetDescription("cons -sequence a.fa -identity 0"),
        new TestCaseData(new[] { "CHHMM?AC.~TRBNKRV-KNTV", "X-X-KN.TSAR-~C?C-.SG-T", "XUT.A-BRTWS.Y~HCTVTYS." }, 1.5f, 0, null, "nnnnnnnnnAnnnnnCnnnnnn").SetName("EmbossAuto_010").SetDescription("cons -sequence a.fa -plurality 1.5"),
        new TestCaseData(new[] { "AT~CT", "TCGCG", "TATGA", "aatcg", "~TAAC", "~CGCA", "AATGC", "TTCTA", "C-GTC" }, 5.65f, 8, null, "NNNNN").SetName("EmbossAuto_011").SetDescription("cons -sequence a.fa -plurality 5.65 -identity 8"),
        new TestCaseData(new[] { "-GA..", "CGCAGG" }, null, 0, null, "cGaagg").SetName("EmbossAuto_012").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { ".-X-S", ".Y" }, 1.5f, 0, null, "nnnnn").SetName("EmbossAuto_013").SetDescription("cons -sequence a.fa -plurality 1.5"),
        new TestCaseData(new[] { ".SRA", "--VXH", "CT-~PEK~-", "ni.-..i~" }, -1.0f, 0, null, "-i--hek--").SetName("EmbossAuto_014").SetDescription("cons -sequence a.fa -plurality -1"),
        new TestCaseData(new[] { "BWRYUT?BBKN", "KKHYANMYREQ", "DCSRIMVDDCD", "NGFMSCFIYHC", "ILTFHVMHRGP", "TVFYWMLYHSM", "mtmfdkksleh", "SDFGALMCIIA" }, 2.0f, 7, 0.0f, "NNNNNNNNNNN").SetName("EmbossAuto_015").SetDescription("cons -sequence a.fa -plurality 2 -identity 7 -setcase 0"),
        new TestCaseData(new[] { "TTGA", "cttg", "TCAG", "GCTG", "tagc", "GCGC", "GAAG" }, 0.74f, 6, null, "NNNN").SetName("EmbossAuto_016").SetDescription("cons -sequence a.fa -plurality 0.74 -identity 6"),
        new TestCaseData(new[] { "ACGGGG", "T", "H", "GATT", "A.ACA", "T", "P~RS" }, null, 0, null, "nnnnnn").SetName("EmbossAuto_017").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "~--.GT-TC~ATC", "--G~TAGGGG.T-C", "AGGGAAGCTTTTAT", "..C.C.TA-TTCG.TC", "a~-c.~-gt-gagcaa~~~", "C.AAGAA.-CATA", "~ag-c-c.-a-ac-t", "G.~CA-" }, 1.5f, 0, null, "nngnnannnnntnctnnnn").SetName("EmbossAuto_018").SetDescription("cons -sequence a.fa -plurality 1.5 -identity 0"),
        new TestCaseData(new[] { "YNWKK.LVEN", "TTCCGAAGCTGA" }, null, 0, null, "ynwkkalvenga").SetName("EmbossAuto_019").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "ARI", "rtneyl-tph~gddkaf-qmsmdmt", "~CT.ATCGGTC-" }, null, 0, null, "xxxxxxxxxxxxxxxxxxxxxxxxx").SetName("EmbossAuto_020").SetDescription("cons -sequence a.fa -identity 0"),
        new TestCaseData(new[] { "XTWB~N", "KL.EAELIISHHKA", "R-C.QY*M", "..nsea", "KIG", "-", "-", "A.B..~.N~FX", "lswcfvfzmp.cachh" }, 2.0f, 0, null, "nnnnnnnnnnnnnnnn").SetName("EmbossAuto_021").SetDescription("cons -sequence a.fa -plurality 2"),
        new TestCaseData(new[] { "GRKNNSESVEITY", "RYIDHWEKHPFHR", "hwysqmvwmmqps", "kmcqnmthpdnfm", "AATGGGAGGGCCC", "ATACGCGCCTTTG" }, 3.0f, 4, null, "XXXXXXXXXXXXX").SetName("EmbossAuto_022").SetDescription("cons -sequence a.fa -plurality 3 -identity 4"),
        new TestCaseData(new[] { "-ATTCCATCA", "TAATCTAG.ATGCG", "t-cgtg-tt--tctgta", "G", "TTAGTGTTGCG", "CTTCGGAAATCTGG" }, null, 0, null, "tnnnngatnnnnnnnnn").SetName("EmbossAuto_023").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "CGCCG", "VGD.", "L.KI~QHYT", "ndkn", "LCSSM", "WSC~~GQ", "-Q~", "A.HYHFT" }, null, 0, null, "nnnnnnnnn").SetName("EmbossAuto_024").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "ny-a", "napfh.~sprcfknc", "AQRHD", "LG-VETEFSC", "-", "PMCLIR" }, null, 0, null, "nnnnnnnnnnnnnnn").SetName("EmbossAuto_025").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "-M~BVCG", "bast~-g", "-AUS~KU" }, null, 0, null, "nAnnnnG").SetName("EmbossAuto_026").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "YUAKBRMSWWKWSVRBSXWHKN?TURCBUCMKKBGUNCCDURSKY", "GSXUGYNRUVYMNBRNDWBHRDVC?XKUNATHYKWTDY?K??Y??" }, null, 0, null, "nuaknnnnwnnnnnnnnnnnnnntuncnucnnnngUnCcnunnnn").SetName("EmbossAuto_027").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "M?CMQFTGSXXKYVNIY.", "cs~s-.h-*-.nbeppvs", "bepndz*bzbfqpz-z~r", "QQMEXFBEMIGCCMHZPE", "WGKCF*TSQLVT*XKLCR", "WNTPBIQHV?IEYARKMT" }, null, 0, null, "xxxxxxxxxxxqxxxzxx").SetName("EmbossAuto_028").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "TTGCCGCAGCTGGAG.TCGATG~T", "mpienvig.madte.cq.vkelka", "YDQFTDPYKQ-NFMQTLNCYSGRC", "-RYCHCCLWGQILAMFYVNIGATV" }, null, 0, null, "nnncnncnnnnnnAnntnnnngnn").SetName("EmbossAuto_029").SetDescription("cons -sequence a.fa -identity 0"),
        new TestCaseData(new[] { "T~T", "y", "E" }, null, 0, null, "Tnn").SetName("EmbossAuto_030").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "csiyfsavyfplciacw", "ADRGWHDWWVKPICLCS", "CLCYDDAGVLTDQPNDK", "PFDKVLRIDEEDLQDTV", "idcshpimnimmerkkq", "VMFSPKAQSMYSNKKEE" }, 2.0f, 0, null, "cdcsfxamyIkdqkkcq").SetName("EmbossAuto_031").SetDescription("cons -sequence a.fa -plurality 2"),
        new TestCaseData(new[] { "?~?WBNCAK", "?SUBHG?X~" }, null, 0, null, "nnnnnncan").SetName("EmbossAuto_032").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "DX~VS-A.", "BNB" }, 1.5f, 0, null, "nnnnnnnn").SetName("EmbossAuto_033").SetDescription("cons -sequence a.fa -plurality 1.5"),
        new TestCaseData(new[] { "IPYKWHLTKDLQIYT", "LTDSEPVYDNKMFAN", "VKYGRIYKWINDQNE", "HIGVKQTCHRRDTMT" }, 1.68f, 3, null, "XXXXXXXXXXXXXXX").SetName("EmbossAuto_034").SetDescription("cons -sequence a.fa -plurality 1.68 -identity 3"),
        new TestCaseData(new[] { "wbxrhwktmdyvymyhy", "X?V-NGMB-NH-NYTWK", "WW-D?AVXKAAWT.URB" }, null, 1, null, "nnnnnnnnnnnnnnTnn").SetName("EmbossAuto_035").SetDescription("cons -sequence a.fa -identity 1"),
        new TestCaseData(new[] { "C.GCTTC~GCT", "L~GRPEDISWE", "gg-ctgac~ac", "TCCGAT.TCC~" }, 1.5f, 0, null, "nngcttnnscn").SetName("EmbossAuto_036").SetDescription("cons -sequence a.fa -plurality 1.5"),
        new TestCaseData(new[] { "~~~C.-R", ".*-MFAPCK-K", "k-qecf.", "ICEXPAKP.", "-.g.-kfida", "ZTDZ.-*", "ZZR*KWV~E", "xs~l-l.-ya", "-.ppv" }, null, 0, null, "nnnnnnnnnnn").SetName("EmbossAuto_037").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "?KDHBKT", "EW" }, 1.5f, 0, 0.0f, "nnnnnnN").SetName("EmbossAuto_038").SetDescription("cons -sequence a.fa -plurality 1.5 -identity 0 -setcase 0"),
        new TestCaseData(new[] { "AEA?FBLCVSCMMVCPCYZE*SAWMB", "wfvqg*zdpmhrhz*izg?zkqbmbe" }, null, 0, null, "aeaxfblcvscmmvcpcyzE*sawmB").SetName("EmbossAuto_039").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "~~A~GG.G.GAA-~.-AG.A~~CTA", "AC.C-T-.-T-~AAG~G-G.GGTA~", ".-~~GATCG.GA.-GG~~TAGG~G.", "ATT-A~GCGAA.ACTCGTGCT~T.T" }, 2.46f, 1, null, "nnnnnnnnnnnnnnnnnnnnnnnnn").SetName("EmbossAuto_040").SetDescription("cons -sequence a.fa -plurality 2.46 -identity 1"),
        new TestCaseData(new[] { "YZFPAWAZEPLWTNVNE*", "**KWPQ*DDGIDLPGMLFFVGWFSI*", "vsrbqewwvwfkzwpagerehqw*z" }, 3.0f, 0, null, "xxXxxXxXXxXxxxxxxxxxxxXxxx").SetName("EmbossAuto_041").SetDescription("cons -sequence a.fa -plurality 3 -identity 0"),
        new TestCaseData(new[] { "IBRKAS", "ELKSRTBCRNQCN.SI", "-.rtqtpwcndp", "EX-WWCBSH", "LYSQFLSXBMYK.L*KSIS", "*E*TYZBF", "BVZYTVE*MZNT*BQIQ", "Z.W.R*", "D~PPZN" }, null, 0, null, "Exxxxxxxxxxxxxxxxxx").SetName("EmbossAuto_042").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "BH-~.~~G", "b.g~xy?-mww.", "G--VGA", "HDGUUT~HTXRX", "ahha", "..MA-K.YRY", "BDM.K-C~.", ".-AMRNA" }, null, 0, null, "nnnnnnnnnnnn").SetName("EmbossAuto_043").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "STZKPNV", "qgpbbm*dzxzbr", "BVHK", "ESLIBIIQVNBY", "NIAGXYAX", "DWXFFKH", "HEYMFDRVZH", "LS" }, null, 6, null, "XXXXXXXXXXXXX").SetName("EmbossAuto_044").SetDescription("cons -sequence a.fa -identity 6"),
        new TestCaseData(new[] { "MSBNG*VTBZK", "GW*FDD*ZPRB", "?SFKFTDLED~", "NTR~SYH*KEA" }, null, 0, null, "xSxxxxxxEEx").SetName("EmbossAuto_045").SetDescription("cons -sequence a.fa -identity 0"),
        new TestCaseData(new[] { "xxknh", "VCYHC", "YNEPM", "YGBVZ", "*YIGC", "*X-HV" }, 2.0f, 0, 1.0f, "nnnnC").SetName("EmbossAuto_046").SetDescription("cons -sequence a.fa -plurality 2 -setcase 1"),
        new TestCaseData(new[] { ".G-GAT", "C", "-", "ACCAT", "GC~CC", "~TGATGC", "A~TC", "TT-A", ".GTA" }, null, 4, null, "NNNnNNN").SetName("EmbossAuto_047").SetDescription("cons -sequence a.fa -identity 4"),
        new TestCaseData(new[] { "ysrenq", "HEEVHQWLQA", "SQHIPVRPPLHLPRN", "ARRQ", "FDYGQEERGNSI", "ILMMLPQEMHARTRNY", "WDAFRTK" }, null, 0, null, "YExxxxQxxxxxxxxx").SetName("EmbossAuto_048").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { ".~VTWI", "IC-EKFDE", "L", "TFYMEQG", "dynswlvifkhc" }, null, 0, null, "xxxxxxxxxxxx").SetName("EmbossAuto_049").SetDescription("cons -sequence a.fa -identity 0"),
        new TestCaseData(new[] { "M", "ABAUM", "U?ACS", "ARU", "x?wg?xcauydcwmnu", "byx" }, 0.0f, 0, null, "a-w--ncauydcwmnu").SetName("EmbossAuto_050").SetDescription("cons -sequence a.fa -plurality 0"),
        new TestCaseData(new[] { "HPGAME", "X.THEI", "ZZ~NQT", "HLMABR", "RHZHAP", "*HEXHQ", "MC-EAX", "VKLSPY", "LRS*NX", "CH*XHE", "pbkhmz" }, 1.48f, 0, null, "hhxhqe").SetName("EmbossAuto_051").SetDescription("cons -sequence a.fa -plurality 1.48"),
        new TestCaseData(new[] { "V", "N", ".", "W", "?", "T" }, 3.0f, 5, null, "N").SetName("EmbossAuto_052").SetDescription("cons -sequence a.fa -plurality 3 -identity 5"),
        new TestCaseData(new[] { "ABYYK-", "LRNYB.", "D-.-BB", "MYIIA-", "HV-~YG", "kmfqep", "i~.isz" }, 3.0f, 0, 1.0f, "nnnnnn").SetName("EmbossAuto_053").SetDescription("cons -sequence a.fa -plurality 3 -identity 0 -setcase 1"),
        new TestCaseData(new[] { "GTTCC", "RNC", "Y", "YNGEP", "FMEMQ" }, null, 0, null, "nnnnn").SetName("EmbossAuto_054").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "sk", "IKE", "ra~", "FQ", "M-D", "hymg", "C-S" }, 3.0f, 0, null, "nnnn").SetName("EmbossAuto_055").SetDescription("cons -sequence a.fa -plurality 3"),
        new TestCaseData(new[] { "~~c", "RC~", "DS-", ".WG", "fss" }, 1.3f, 0, null, "nns").SetName("EmbossAuto_056").SetDescription("cons -sequence a.fa -plurality 1.3"),
        new TestCaseData(new[] { "?aryss", "-a.grf", "VTYC~B" }, 0.0f, 0, 1.0f, "-A-y-f").SetName("EmbossAuto_057").SetDescription("cons -sequence a.fa -plurality 0 -identity 0 -setcase 1"),
        new TestCaseData(new[] { "C-~~H~UDVW~?~.YUMB.G.CR-CWC", "RR~SIHI-STZWGCHWSHWPND*~QWQ" }, 1.0f, 2, 2.5f, "NNnNNNNNNNnNNNNNNNNNNNNnNnN").SetName("EmbossAuto_058").SetDescription("cons -sequence a.fa -plurality 1 -identity 2 -setcase 2.5"),
        new TestCaseData(new[] { "D-N~I~MQS~PEP-N-M.CLS", "..CG~AGGTGACGATTAGTTA", "ATATGC~AGC~~AAAGCACT~", "ITM-.ECYW-IMTND.~TTSH" }, 0.69f, 0, 0.0f, "xTAxxxxASxxxAANxxACTS").SetName("EmbossAuto_059").SetDescription("cons -sequence a.fa -plurality 0.69 -setcase 0"),
        new TestCaseData(new[] { "KICLLWW", "GQTIRYQ", "QYLPY", "DDQDEM", "earw", "INICE", "flk", "R" }, 4.14f, 0, null, "xxxxxxx").SetName("EmbossAuto_060").SetDescription("cons -sequence a.fa -plurality 4.14 -identity 0"),
        new TestCaseData(new[] { "qsvymlrvinrr", "AVCGIHIIAWAHH", "WCPDS", "ptpwmc", "ESDSTEH", "QDWVWSA", "REHPYSME", "KVCAYA", "SAVFKMPEWKHWS" }, null, 0, null, "Qxxxxxxxxxxxx").SetName("EmbossAuto_061").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "SPZCSTIH~", "nkrhhbltz", "-HXGV?BDB", "*CCHW~B~P", "ZTHXW?LMF", "MNXP*KSIL" }, null, 1, null, "xxxxxxixx").SetName("EmbossAuto_062").SetDescription("cons -sequence a.fa -identity 1"),
        new TestCaseData(new[] { "ACCCCCAAGCAATACGTTATAA", "HILYIWLSVQFCYNHKNYAASH", "AGGCGGAGATATCTAACGGACG", "ctgccatcctgttgcagatcta" }, 3.0f, 0, null, "nnnCnnnnnnnnYnnnnnnnnn").SetName("EmbossAuto_063").SetDescription("cons -sequence a.fa -plurality 3"),
        new TestCaseData(new[] { "-A", "KP", "F", "MG", "M", "-", "E.", "L" }, null, 0, null, "nn").SetName("EmbossAuto_064").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "DNN~.VWSUR~-~MWAVMWT?M.V", "-C~S~Y.FT-VCHQGQ~.YVIEFS", "cyh~lrvyfwiwleynets.ldq-", "KS~W.-G-~LR~RVH-PFK~VREP" }, null, 0, null, "nCnnnnnnunnnnnnnnnnnnnnn").SetName("EmbossAuto_065").SetDescription("cons -sequence a.fa"),
        new TestCaseData(new[] { "g", "ACCACTTT", "GACCGC", "CCGTCA" }, 0.27f, 0, null, "gccncntt").SetName("EmbossAuto_066").SetDescription("cons -sequence a.fa -plurality 0.27"),
        new TestCaseData(new[] { "-asdbdbxdhyus?rhm", "gk", "XK*DSRB~X*WWWKZ", "MQSKHSYTD~V*ER", "*CP-NR~WH~SBLNWXNRK", "YXKPYCQNP*PQGC" }, -1.0f, 5, null, "NNNNNNNNNNNNNNNNNNN").SetName("EmbossAuto_067").SetDescription("cons -sequence a.fa -plurality -1 -identity 5"),
        new TestCaseData(new[] { "TGCT", "AGPD", "GMWD", "GMTH", "mvwi", "VIWM", "CAVW", "GETQ" }, 1.5f, 0, null, "gnnn").SetName("EmbossAuto_068").SetDescription("cons -sequence a.fa -plurality 1.5"),
        new TestCaseData(new[] { "YCHR~AGAEWFN~FQWGMMKPQ", "ptwcmv.ivsgpqld.serkqd", "NC.PEFIHEL~TM.QYLQRMLC" }, 0.54f, 0, null, "nCxrxvxaElxtqfQWsQRKqq").SetName("EmbossAuto_069").SetDescription("cons -sequence a.fa -plurality 0.54 -identity 0"),
    };

    [TestCaseSource(nameof(AutoCases))]
    public void GenerateEmbossConsensus_AutoPadded_MatchesConsBinary(
        string[] rows, float? plurality, int identity, float? setcase, string expected)
    {
        Assert.That(
            MotifFinder.GenerateEmbossConsensus(rows, true, ConsensusResidueType.Auto, plurality, identity, setcase),
            Is.EqualTo(expected));
    }

    /// <summary>
    /// Sequence typing details, each confirmed with the cons binary (-plurality 0): a
    /// nucleotide-looking row reads X and ? as N (gapany ? → X, then ajSeqSetNuc X → N), a
    /// protein-looking row keeps an unscored X, and an X emitted into a nucleotide set is written
    /// as N (cons.c ajSeqSetNuc on the output sequence).
    /// </summary>
    [Test]
    public void GenerateEmbossConsensus_Auto_PerSequenceTyping_MatchesConsBinary()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "A-", "EX", "EX" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("an"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "A-", "E?", "E?" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("an"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "A-", "E*", "E*" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("a*"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "E-", "AX", "EX" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("E-"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "EXGT", "AXGT", "A?GT" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("ANGT"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(new[] { "A?GT", "E?GT", "E?GT" }, ConsensusResidueType.Auto, 0f), Is.EqualTo("anGT"));
        });
    }

    /// <summary>
    /// Auto reproduces cons where an explicit type cannot: cons takes N/X from the first sequence
    /// (cons -sequence a.fa -identity 0 -setcase 2.15 on V, x, ~ → "n"; the explicit Protein
    /// type gives "x").
    /// </summary>
    [Test]
    public void GenerateEmbossConsensus_Auto_FirstSequenceDecidesType()
    {
        string[] rows = { "V", "x", "~" };
        Assert.That(MotifFinder.GenerateEmbossConsensus(rows, ConsensusResidueType.Auto, setcase: 2.15f), Is.EqualTo("n"));
        Assert.That(MotifFinder.GenerateEmbossConsensus(rows, ConsensusResidueType.Protein, setcase: 2.15f), Is.EqualTo("x"));
        // Protein alignment with a protein-looking first row: Auto = explicit Protein.
        string[] protein = { "MKVLAAGIVG", "MKVLSAGIVA", "MRVLAAG-VG", "MKILTAGLVG" };
        Assert.That(MotifFinder.GenerateEmbossConsensus(protein, ConsensusResidueType.Auto),
            Is.EqualTo(MotifFinder.GenerateEmbossConsensus(protein, ConsensusResidueType.Protein)));
    }

    /// <summary>
    /// cons -snucleotide sets each sequence nucleotide before the gapany ? → X conversion
    /// (ajSeqTypeCheckIn), so '?' stays an unscored X (EDNAFULL code 0) and an emitted X is
    /// written as N; without the flag (Auto) '?' → X → N is scored as N. Binary outputs:
    /// one-column alignments, one row per character.
    /// </summary>
    [Test]
    public void GenerateEmbossConsensus_Nucleotide_QuestionMarkIsUnscored_AsConsSnucleotide()
    {
        static string[] Col(string c) => c.Select(ch => ch.ToString()).ToArray();
        Assert.Multiple(() =>
        {
            // cons -plurality 4.51 -setcase 1.1 -snucleotide → n ; without -snucleotide → N
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("USC?TK"), ConsensusResidueType.Nucleotide, 4.51f, 0, 1.1f), Is.EqualTo("n"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("USC?TK"), ConsensusResidueType.Auto, 4.51f, 0, 1.1f), Is.EqualTo("N"));
            // cons -plurality 1.04 -setcase 3.17 -snucleotide → n ; with N instead of ? → w
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("AWASGd?U"), ConsensusResidueType.Nucleotide, 1.04f, 0, 3.17f), Is.EqualTo("n"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("AWASGdNU"), ConsensusResidueType.Nucleotide, 1.04f, 0, 3.17f), Is.EqualTo("w"));
            // cons -plurality 0 -snucleotide on ???: the emitted X is written as n
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("???"), ConsensusResidueType.Nucleotide, 0f), Is.EqualTo("n"));
            // X itself is read as N under -snucleotide (USCXTK → N)
            Assert.That(MotifFinder.GenerateEmbossConsensus(Col("USCXTK"), ConsensusResidueType.Nucleotide, 4.51f, 0, 1.1f), Is.EqualTo("N"));
        });
    }

    [Test]
    public void GenerateEmbossConsensus_PadRaggedRows_EqualsExplicitGapPadding_AndOriginalOverloadUnchanged()
    {
        string[] ragged = { "ACGTAC", "ACG", "AC" };
        string[] padded = { "ACGTAC", "ACG---", "AC----" };
        Assert.Multiple(() =>
        {
            // cons -sequence a.fa (rows ACGTAC, ACG, AC) → ACGnnn
            Assert.That(MotifFinder.GenerateEmbossConsensus(ragged, true), Is.EqualTo("ACGnnn"));
            Assert.That(MotifFinder.GenerateEmbossConsensus(ragged, true),
                Is.EqualTo(MotifFinder.GenerateEmbossConsensus(padded)));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossConsensus(ragged));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossConsensus(ragged, false));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossConsensus(new[] { "AC", null! }, true));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossConsensus(new[] { "AC" }, true));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.GenerateEmbossConsensus(null!, true));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateEmbossConsensus(padded, true, (ConsensusResidueType)9));
        });
    }
}
