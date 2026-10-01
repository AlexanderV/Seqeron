// MOTIF-GENERATE-001 — Bioconductor DECIPHER ConsensusSequence (B05 audit group C)
// Evidence: docs/Evidence/MOTIF-GENERATE-001-Evidence.md
// Source (opened): DECIPHER 3.9.4 R/ConsensusSequence.R, src/ConsensusSequence.c,
// man/ConsensusSequence.Rd (raw.githubusercontent.com/bioc/DECIPHER/devel).
// Reference: every expected string below is the output of DECIPHER's own R + C source
// (ConsensusSequence.R / ConsensusSequence.c verbatim, built as an R package against
// R 4.3.3 + Biostrings 2.70.2, Ubuntu noble) — 10,000 random cases, 0 mismatches; the
// Decipher_* cases are a stratified sample of that run (seed 20261001).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_DecipherConsensus_Tests
{
    private static readonly TestCaseData[] DecipherCases =
    {
        new TestCaseData(new[] { "naa.", "." }, DecipherSequenceType.Dna, 0.99, false, '+', 0.5, false, false, "-AA-").SetName("Decipher_000"),
        new TestCaseData(new[] { "--NNCCUGUGUAAN", "NGUNACUNGNCGAG", "CAAGNCGCCCCAGU", "NAUCGAAACAUCNA" }, DecipherSequenceType.Rna, 0.99, false, 'A', null, false, true, "CAUSVCUVCVYAAD").SetName("Decipher_001"),
        new TestCaseData(new[] { "---H.B-A-DAA-N.B", "YN.U--WU-YM-.U", "B..", "DCKGGR.A..", "HYY-." }, DecipherSequenceType.Rna, 0.05, false, 'N', null, false, false, "-C-KN--W--AN-U--").SetName("Decipher_002"),
        new TestCaseData(new[] { "----++UUUAAAAUAUUA+U+AAUAUUU", "UAUUUAA+AUAUAAAAUUUUA", "U+AAAUUUAUUAUUUUAUUAUUUAU", "AUA+AUAU+UUUU+UUU+U+UAUAU", "UUAUAUAUUUUUA+AAUUUAUUAUU" }, DecipherSequenceType.Rna, 0.5, true, '-', null, true, true, "UUAUAUAUWUUUAUAUUUUWUWWWUUUU").SetName("Decipher_003"),
        new TestCaseData(new[] { "-acaaa---cccac-aa-aa", "A-AAACAA--A-AAAA-A-C" }, DecipherSequenceType.AminoAcid, 0.05, true, 'N', null, true, false, "AAXAAXAA-CXCAXAAAAAX").SetName("Decipher_004"),
        new TestCaseData(new[] { "CACC" }, DecipherSequenceType.Dna, 0.05, true, '+', 0.8, false, true, "CACC").SetName("Decipher_005"),
        new TestCaseData(new[] { "CACCAC", "----CCCCAAAAAAACC", "caccacccaaccc", "CCACACACACAAC", "CCCACAAACCCCC", "CCCCCACAACACA", "AACCCCACCCACC" }, DecipherSequenceType.Dna, 0.0, false, '-', 1.0, false, false, "MMMMMMMMMMMMMAACC").SetName("Decipher_006"),
        new TestCaseData(new[] { "-AAACAAAAAAACAA", "AAAAAAACAACAAAAAAA" }, DecipherSequenceType.AminoAcid, 0.05, false, '+', 0.5, false, true, "+AAAXAAXAAXAXAAAAA").SetName("Decipher_007"),
        new TestCaseData(new[] { "zzdxbnzzdndzbbjn", "DJBNNNXZJNNXJBZD", "DDDBXDNXNXNDJBB+", "DJNZBXXXZBBBDJXJ", "JXDZBZXN+JX+DBBN" }, DecipherSequenceType.AminoAcid, 0.0, true, 'P', 0.5, true, false, "PXXXXXXXXXXPXXXX").SetName("Decipher_008"),
        new TestCaseData(new[] { "uuaaauaaaauuauaaaa" }, DecipherSequenceType.Rna, 0.25, true, '+', null, false, false, "UUAAAUAAAAUUAUAAAA").SetName("Decipher_009"),
        new TestCaseData(new[] { "aaaaccac", "ACACCAAAACAACAAAAA" }, DecipherSequenceType.Dna, 0.1, true, '-', 0.8, false, false, "AMAMCMAMACAACAAAAA").SetName("Decipher_010"),
        new TestCaseData(new[] { "VWT", "WANGT", "NASWYWRC", "", "--VVDGG" }, DecipherSequenceType.Dna, 0.1006, true, '-', null, false, false, "NWNNNDRC").SetName("Decipher_011"),
        new TestCaseData(new[] { "V", "VHGGA-", "KV.N", "-", "R", "V", "----" }, DecipherSequenceType.Rna, 0.3, true, '+', null, true, false, "RNGNA-").SetName("Decipher_012"),
        new TestCaseData(new[] { "AAAAAAA", "-CAAAA", "---CAAA", "AA-AAAA" }, DecipherSequenceType.AminoAcid, 0.1, true, '+', null, true, false, "AXAXAAA").SetName("Decipher_013"),
        new TestCaseData(new[] { "AAAUUUAAAUAUAAAU", "AUAUUAAAUUAUUAAU" }, DecipherSequenceType.Rna, 0.05, false, '+', null, false, false, "AWAUUWAAWUAUWAAU").SetName("Decipher_014"),
        new TestCaseData(new[] { "GAGCGN", "---AUUAGA" }, DecipherSequenceType.Rna, 0.9, true, '-', null, false, true, "GAGMKUAGA").SetName("Decipher_015"),
        new TestCaseData(new[] { "--ZJJNBJZZ", "JXBZXZXN" }, DecipherSequenceType.AminoAcid, 0.99, true, '+', null, false, false, "JXXXJNBNZZ").SetName("Decipher_016"),
        new TestCaseData(new[] { "HT" }, DecipherSequenceType.Dna, 0.1, false, '+', null, false, false, "-T").SetName("Decipher_017"),
        new TestCaseData(new[] { "++++B++++YV+G++++HA++K+", "--VK+V++++C++++M++++B+N", "++B++R++SMA+++Y+Y+S+++G++D", "GM++", "++++N+R+++++++++++R+", "---+G+++++++++R+S++GM+++W", "DC+RY++++H++" }, DecipherSequenceType.Dna, 0.99, true, '+', null, true, true, "GCSGGRR+SCM+G+NMCHAGCKG+WD").SetName("Decipher_018"),
        new TestCaseData(new[] { "YCBVVNYWBGMAGABTV-W-", "tnk-d-knkcymc-bamhrs" }, DecipherSequenceType.Dna, 0.25, true, 'N', 0.5, true, false, "TNKVNNBWKSHASABWMHDS").SetName("Decipher_019"),
        new TestCaseData(new[] { "O*I*UUIU*LII", "---IUJUOUI*JJ", "---ILI*J*L*OO", "----*", "ojuoloji*u", "JUUJU**L*L", "*IJL*IUIOU" }, DecipherSequenceType.AminoAcid, 0.5, false, '-', null, false, false, "OXUJXXXJ*L*XO").SetName("Decipher_020"),
        new TestCaseData(new[] { "ACCAACCAAA", "aaaacca", "aacaa+a", "aacaaca", "ACA+AAA", "AACAACA", "+CAAAAA" }, DecipherSequenceType.Dna, 0.3, false, 'N', 1.0, false, false, "NMMNNNNAAA").SetName("Decipher_021"),
        new TestCaseData(new[] { "U", "AUAUAU", "uauaau" }, DecipherSequenceType.Rna, 0.25, true, '-', 0.5614, false, true, "WWWWAU").SetName("Decipher_022"),
        new TestCaseData(new[] { "-AAA", "AAAACA", "AAAACA", "aaaaaa", "AAAAAA", "CAAAAA", "----AAAAA", "AAAAAA", "aaaaaa", "A" }, DecipherSequenceType.AminoAcid, 0.0, false, '+', null, false, true, "++++XAAAA").SetName("Decipher_023"),
        new TestCaseData(new[] { "--C-C-C---AAA", "---", "---A-AA----C-CA", "A-A-A" }, DecipherSequenceType.Dna, 0.1, false, '+', null, false, false, "A-M+++M---+M+CA").SetName("Decipher_024"),
        new TestCaseData(new[] { "----ccaaaacaccaaaaacccaaa", "----CCAAACACCAACACCCACCCAAC" }, DecipherSequenceType.Dna, 0.7744, true, '-', 0.5, false, false, "----CCAAAMMMCMAMAMMCMCMMAAC").SetName("Decipher_025"),
        new TestCaseData(new[] { "c-aaanga", "-A--NG--", "ng-g--at", "CTTCTATN" }, DecipherSequenceType.Dna, 0.8909, true, 'W', null, false, true, "CDWVWRDW").SetName("Decipher_026"),
        new TestCaseData(new[] { "RN--AD" }, DecipherSequenceType.AminoAcid, 0.75, true, '+', 0.01, false, true, "RN--AD").SetName("Decipher_027"),
        new TestCaseData(new[] { "+++++A++++++++CA++A++AAA" }, DecipherSequenceType.AminoAcid, 0.05, true, 'N', null, false, true, "+++++A++++++++CA++A++AAA").SetName("Decipher_028"),
        new TestCaseData(new[] { "G++TGTG+++C+TG+AA+C++++", "++++++++ACC+CC++++TA+G+G", "+++++C+C++CT++CG+++G++G" }, DecipherSequenceType.Dna, 0.5, true, '+', null, false, true, "+++++Y++++C+YS+R++YR+++G").SetName("Decipher_029"),
        new TestCaseData(new[] { "-cc--ua----n--g---gagn-", "----NU--A-UG-G------N--" }, DecipherSequenceType.Rna, 0.05, false, '-', 1.0, true, true, "-----U-----G--------G--").SetName("Decipher_030"),
        new TestCaseData(new[] { "TAATAATTAAAATT", "AAATAATATATTAT", "---AATTATAAATAA", "aaaattaatattaa", "AAAAATTATATAAAA" }, DecipherSequenceType.Dna, 0.25, true, 'N', null, false, false, "AAAWAWTATAWWWWA").SetName("Decipher_031"),
        new TestCaseData(new[] { "ITAIEANI", "YHNLR", "MGA", "ayw", "msr", "QIG", "CGP" }, DecipherSequenceType.AminoAcid, 0.1, true, 'P', 1.0, true, false, "XPPJXANI").SetName("Decipher_032"),
        new TestCaseData(new[] { "CU", "nn", "CA", "UG", "AA", "nc", "UA", "U", "NN", "AGCNANG" }, DecipherSequenceType.Rna, 0.75, true, 'A', null, false, false, "UACNANG").SetName("Decipher_033"),
        new TestCaseData(new[] { "*IIUJ*", "IJUL" }, DecipherSequenceType.AminoAcid, 0.75, false, 'W', null, true, true, "XIXXW*").SetName("Decipher_034"),
        new TestCaseData(new[] { "ATTAT", "--TT", "TTTTT", "AAAAT", "AATAT", "AATATA", "AAATA" }, DecipherSequenceType.Dna, 0.99, true, 'N', 1.0, true, false, "NNNNNA").SetName("Decipher_035"),
        new TestCaseData(new[] { "----A", "---+R" }, DecipherSequenceType.AminoAcid, 0.99, true, '-', 0.5, false, false, "---+X").SetName("Decipher_036"),
        new TestCaseData(new[] { "-.-D.", ".D.N--..", "..d-n.db" }, DecipherSequenceType.AminoAcid, 0.5, false, '+', 0.8, false, false, "-D++N-D-").SetName("Decipher_037"),
        new TestCaseData(new[] { ".-+F+E.S.KR", "VR+P+QD..+-", "+M.+ENE++.+", "-++" }, DecipherSequenceType.AminoAcid, 0.0, false, '+', null, false, false, "+++++X+++++").SetName("Decipher_038"),
        new TestCaseData(new[] { "iiliilliliiliiii", "illllililiiiiiii" }, DecipherSequenceType.AminoAcid, 0.25, true, '+', 0.7601, true, false, "IJLJJJLILIIJIIII").SetName("Decipher_039"),
        new TestCaseData(new[] { "a", "A" }, DecipherSequenceType.Rna, 0.9853, true, 'N', null, true, false, "A").SetName("Decipher_040"),
        new TestCaseData(new[] { "-.N.VV+.+.PH-P.KK-.", "----A.+M.E.++--." }, DecipherSequenceType.AminoAcid, 0.0, true, 'N', 0.5, true, true, "NNNNXVNMNEPHNPNKKNN").SetName("Decipher_041"),
        new TestCaseData(new[] { ".UAAU-U..AU.-A-A-.-U-", ".-.U-UAA-UA-A-A.A-.A...U", "--U.--", "-u-uuu.a-uaua-uaa--u.a..." }, DecipherSequenceType.Rna, 0.05, true, '+', null, true, false, "-UWWUUWA-WWUAAWAA--W-A-U-").SetName("Decipher_042"),
        new TestCaseData(new[] { "AA" }, DecipherSequenceType.Dna, 0.9, true, '+', null, true, false, "AA").SetName("Decipher_043"),
        new TestCaseData(new[] { "---.-..-ua-", "-U-U.-.U.", "--.AAA-A-", "a---au-au.u", "u" }, DecipherSequenceType.Rna, 0.25, true, 'N', 0.8, false, false, "WN-NNN-WUNU").SetName("Decipher_044"),
        new TestCaseData(new[] { "", "----", "", "", "----" }, DecipherSequenceType.Dna, 0.75, true, '+', null, true, true, "++++").SetName("Decipher_045"),
        new TestCaseData(new[] { "----j", "---JI", "----*o" }, DecipherSequenceType.AminoAcid, 0.99, false, '-', null, true, false, "----XO").SetName("Decipher_046"),
        new TestCaseData(new[] { "-RNANNDDAN-RADNAR-R-R" }, DecipherSequenceType.AminoAcid, 0.05, false, '+', null, false, true, "-RNANNDDAN-RADNAR-R-R").SetName("Decipher_047"),
        new TestCaseData(new[] { "G.-M.WWSGM", "----a+vr+v-", "NA.N+S.WMT", "--NTBGH.R", "-R-", "k+", "RGV-VT+-SY", ".HMKCVH+BM", "-kyn..+a", "+--.---BHH" }, DecipherSequenceType.Dna, 0.99, true, '+', null, true, false, "GGCTCGAAGC-").SetName("Decipher_048"),
        new TestCaseData(new[] { "----CVUDMANNSHKBDYH", "bvywrcrsmymddsybws", "kcwbrramvuhmmkumrs" }, DecipherSequenceType.Rna, 0.05, false, 'N', null, true, false, "-C--CCW--W----U----").SetName("Decipher_049"),
        new TestCaseData(new[] { "nykubydvwsrnwbhgbuwwc", "ycs+wvmknwvsdnvdhkysc" }, DecipherSequenceType.Rna, 0.3, true, 'N', null, false, false, "YCBNNNNNWNRSWBNDNUHNC").SetName("Decipher_050"),
        new TestCaseData(new[] { "-C" }, DecipherSequenceType.Dna, 0.3, false, '+', 0.8, true, true, "-C").SetName("Decipher_051"),
        new TestCaseData(new[] { "-AAGCTTTTCCG", "CCCGATGATCTAA" }, DecipherSequenceType.Dna, 0.9, true, '+', 0.1721, true, false, "CMMGMTKWTCYRA").SetName("Decipher_052"),
        new TestCaseData(new[] { "D-...+-B-.", "+.-++..-D.", "----+N..++-.-", "---NND.B" }, DecipherSequenceType.AminoAcid, 0.25, true, 'D', null, false, true, "D--DDD-DDD---").SetName("Decipher_053"),
        new TestCaseData(new[] { "---YY+", "+" }, DecipherSequenceType.Dna, 0.0, false, '-', 0.8, true, false, "+----+").SetName("Decipher_054"),
        new TestCaseData(new[] { "-..D-DBD" }, DecipherSequenceType.AminoAcid, 0.1, true, 'N', null, true, true, "NNNDNDBD").SetName("Decipher_055"),
        new TestCaseData(new[] { "-----AR---D-R-A--", "A---NR-RDNNN---" }, DecipherSequenceType.AminoAcid, 0.0, true, '+', 0.1392, false, true, "+---+X++++B++-+--").SetName("Decipher_056"),
        new TestCaseData(new[] { "UUUAUAUAAUAAUUAUUUAUUAAUA", "UAAUUAUAAUAU", "UUAUAAUAUUUAAUAAUAUUAUAUU", "AAUUUAUAUUUAUAUUAUUAUUUUUUU", "----AUAUUUUAUAAAAAAAAUAUAAAA" }, DecipherSequenceType.Rna, 0.1827, false, '+', null, false, true, "++++WWWWWUWWWWWWWWWWWWWUWWWA").SetName("Decipher_057"),
        new TestCaseData(new[] { "--A-", ".TK", "", "SMK" }, DecipherSequenceType.Dna, 0.05, false, '-', 0.01, true, false, "-TA-").SetName("Decipher_058"),
        new TestCaseData(new[] { "--ADANDN", "+ARNDRNNN", "DDRARNADD", "RNNRND+RN" }, DecipherSequenceType.AminoAcid, 0.0, true, '*', null, true, false, "XXXXXXXXB").SetName("Decipher_059"),
        new TestCaseData(new[] { "DDBB" }, DecipherSequenceType.AminoAcid, 0.9, true, '*', 0.5, false, false, "DDBB").SetName("Decipher_060"),
        new TestCaseData(new[] { "", "AA", "-", "", "A", "-AA", "UA" }, DecipherSequenceType.Rna, 0.05, true, 'A', 1.0, true, false, "WAA").SetName("Decipher_061"),
        new TestCaseData(new[] { "RBD-SN", "C-RNM" }, DecipherSequenceType.Dna, 0.0, true, '+', null, false, true, "V+++VN").SetName("Decipher_062"),
        new TestCaseData(new[] { "++a+ac+", "A++AC+A", "+C+AA+CCA" }, DecipherSequenceType.Dna, 0.5, false, 'V', null, true, false, "+++AA+MCA").SetName("Decipher_063"),
        new TestCaseData(new[] { "GGTCCAG", "GCGGTAG", "----ATCCCCA" }, DecipherSequenceType.Dna, 0.99, false, 'N', null, false, false, "GSKSHAGCCCA").SetName("Decipher_064"),
        new TestCaseData(new[] { "+N+T++++NA+NCT+CAG+N+", "C++G+++C++T+N++C+G++", "+gc+c+a++++c+ctc+c+", "+TCN+AGG+++++++++", "-++++G++A++++CC" }, DecipherSequenceType.Dna, 0.1, false, '+', null, false, true, "+++++++++++++++++S+++").SetName("Decipher_065"),
        new TestCaseData(new[] { "BDNDNDDDDD" }, DecipherSequenceType.AminoAcid, 0.9, true, 'S', 1.0, false, false, "BDNDNDDDDD").SetName("Decipher_066"),
        new TestCaseData(new[] { "DDNDANDARDN+ARNADRA", "NRDDANARANANADD+DNN", "ANRRRDRDNRNADR+NAAN", "NANNND+DNAANRDAADAD", "DNDDDA+R+DNRRNNRRDN" }, DecipherSequenceType.AminoAcid, 0.5, true, 'N', 1.0, false, true, "NNNNXNNXNXNNXXNNNXN").SetName("Decipher_067"),
        new TestCaseData(new[] { "AACACCAAACACAAACCC", "CACAACACAAAAACAACC", "CCAAAAAAAACCCCCCCA", "accacccccaaaacccaa" }, DecipherSequenceType.Dna, 0.5, false, 'A', null, false, true, "MMCAMCAMAAAMACMCCM").SetName("Decipher_068"),
        new TestCaseData(new[] { "GGAU---UGGGAAUGACGA-UA", "-CGUACGG-GC-GGAC--ACAA" }, DecipherSequenceType.Rna, 0.99, true, 'N', 0.8, false, false, "GSRUNNNKNGSNRKRMNNANWA").SetName("Decipher_069"),
        new TestCaseData(new[] { "---UANCAAA+", "---u", "---GUNC" }, DecipherSequenceType.Rna, 0.75, true, 'K', null, false, true, "---UWNCAAA+").SetName("Decipher_070"),
        new TestCaseData(new[] { "-AAA---------A-AAA-", "AA-A-A-C-AA-A--A--A", "AA--A--ACA---AA", "-A---------A--A-C" }, DecipherSequenceType.AminoAcid, 0.05, false, 'J', 1.0, true, true, "AAAAAAJXCAAAAAAAXAA").SetName("Decipher_071"),
        new TestCaseData(new[] { "----C..G-N..--.GG.-NGG-N", ".-ua..g." }, DecipherSequenceType.Rna, 0.05, true, '+', 1.0, false, false, "--UA+-+G-N-----GG--NGG-N").SetName("Decipher_072"),
        new TestCaseData(new[] { "-----AACACA", "C-ACA-CCAC", "AACAAAACAC", "ACCCACCCAC", "-acaa-ccca" }, DecipherSequenceType.Dna, 0.5, false, '+', null, false, true, "MMCMAMCCACA").SetName("Decipher_073"),
        new TestCaseData(new[] { "---AAAU--.......-", "----.---A.-.UUAU", ".UA.U-..--A..UA", "---AAAU--UUA.A---", ".-UUUAA-.U--.-.", "---", "-----UA--.U.-U-U.A.", "U.UUA.A--..-", "-.-----AU--.-U.", "-U.AUA-..AU.U--" }, DecipherSequenceType.Rna, 0.1, true, 'N', 0.01, true, false, "UUWWWWWAWWWAUWAU-A-").SetName("Decipher_074"),
        new TestCaseData(new[] { "----a-t-tta", "--at-ta---a", "----ta-----a-a-", "a----t", "T-T-", "--A-TATAAT-", "---t-----a--a-", "T---A-A--A-", "--A--TAA---", "---A-A-T---T---" }, DecipherSequenceType.Dna, 0.3, false, '+', null, true, false, "W-++++++-++++A-").SetName("Decipher_075"),
        new TestCaseData(new[] { "----B", "ndjb", "-BB", "----B", "XXBJ", "-XDXZ", "JJDZ" }, DecipherSequenceType.AminoAcid, 0.99, false, '+', 0.8, false, false, "NDD--").SetName("Decipher_076"),
        new TestCaseData(new[] { "", "--", "-", "-" }, DecipherSequenceType.AminoAcid, 0.0, true, '+', null, false, false, "--").SetName("Decipher_077"),
        new TestCaseData(new[] { "", "AR-", "", "", "", "", "", "", "", "-" }, DecipherSequenceType.AminoAcid, 0.0, true, '+', null, true, true, "AR+").SetName("Decipher_078"),
        new TestCaseData(new[] { "iiil", "", "---L", "--II", "----illi", "LII", "----LLL", "LLLL", "illi", "IIII" }, DecipherSequenceType.AminoAcid, 0.75, true, '-', null, false, false, "IIIJJLLI").SetName("Decipher_079"),
        new TestCaseData(new[] { "--TATATTTTTATAA", "TTAAAATAAAATAATTAT", "TATTTTAAATA", "AAATAATTTATTT" }, DecipherSequenceType.Dna, 0.3, true, 'N', null, false, false, "WWWWWATWWWWWWAWTAT").SetName("Decipher_080"),
        new TestCaseData(new[] { "D+NBNNBDN", "+DB+BDNDB" }, DecipherSequenceType.AminoAcid, 0.25, false, 'N', null, false, false, "NNN+NBNDN").SetName("Decipher_081"),
        new TestCaseData(new[] { "tanatg" }, DecipherSequenceType.Dna, 0.5, false, '-', null, true, false, "TA-ATG").SetName("Decipher_082"),
        new TestCaseData(new[] { "B", "+" }, DecipherSequenceType.AminoAcid, 0.5, true, 'N', null, false, true, "B").SetName("Decipher_083"),
        new TestCaseData(new[] { "--G.U.-", "---+C+C+.", "UU+-+.-AU+A", "A--+--C", "--GCGG+U", "-+cu...", "+.+C+-+", "C+-UC.+", "CGG+U.-", "AC+A-.-" }, DecipherSequenceType.Rna, 0.0, false, 'C', null, false, false, "CCCCCCCCU+A").SetName("Decipher_084"),
        new TestCaseData(new[] { "JUJJO", "+I-O", "-*J**", "J-..-", ".o*jl", "--L*JJ.", "J*", "*j*u+", "+I*L*", "+jl*-" }, DecipherSequenceType.AminoAcid, 0.9, false, '+', 0.8, false, false, "+X+++--").SetName("Decipher_085"),
        new TestCaseData(new[] { "----UUUU" }, DecipherSequenceType.Rna, 0.25, true, '-', 0.8, false, true, "----UUUU").SetName("Decipher_086"),
        new TestCaseData(new[] { "CMK--...SRGT..K.", "BNDKGGHTCK.YNHD-", "DD---Y..VRC.H-A.", "---MKASY-RAATNNWCNW" }, DecipherSequenceType.Dna, 0.1, true, '+', null, true, false, "NNDNKNNYVDVHHNDWCNW").SetName("Decipher_087"),
        new TestCaseData(new[] { "----NDDNBBBDBNBDNDBDDNDDNB", "-BDBDNNDDBNNBNDDDND", "--NDDDDNNBBDNNNNBNBDDDB" }, DecipherSequenceType.AminoAcid, 0.9, true, 'P', null, true, true, "PBBDDDDNBBNDNNBDBNDDDBDDNB").SetName("Decipher_088"),
        new TestCaseData(new[] { "tngnggngcc", "----gtcgagn" }, DecipherSequenceType.Dna, 0.05, true, '+', 1.0, false, false, "TNGNGKNGMSN").SetName("Decipher_089"),
    };

    [TestCaseSource(nameof(DecipherCases))]
    public void GenerateDecipherConsensus_MatchesDecipherSource(
        string[] rows, DecipherSequenceType type, double threshold, bool ambiguity, char noConsensusChar,
        double? minInformation, bool includeNonLetters, bool includeTerminalGaps, string expected)
    {
        Assert.That(
            MotifFinder.GenerateDecipherConsensus(rows, type, threshold, ambiguity, noConsensusChar,
                minInformation, includeNonLetters, includeTerminalGaps),
            Is.EqualTo(expected));
    }

    /// <summary>The examples of DECIPHER man/ConsensusSequence.Rd (stated outputs; reproduced by the R build).</summary>
    [Test]
    public void GenerateDecipherConsensus_ManualExamples()
    {
        string[] aaat = { "A", "A", "A", "T" };
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.GenerateDecipherConsensus(aaat), Is.EqualTo("W"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(aaat, threshold: 0.3), Is.EqualTo("A"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(aaat, threshold: 0.3, minInformation: 0.8), Is.EqualTo("+"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(aaat, threshold: 0.3, minInformation: 0.8, noConsensusChar: 'N'), Is.EqualTo("N"));

            string[] majority = { "GTT", "GAA", "CTG" };
            Assert.That(MotifFinder.GenerateDecipherConsensus(majority), Is.EqualTo("SWD"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(majority, threshold: 0.5), Is.EqualTo("GTD"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(majority, threshold: 0.5, minInformation: 0.75), Is.EqualTo("++D"));

            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "A", "T" }), Is.EqualTo("W"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "A", "T" }, threshold: 0.5), Is.EqualTo("W"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "A", "T" }, DecipherSequenceType.AminoAcid), Is.EqualTo("X"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "A", "T" }, DecipherSequenceType.AminoAcid, 0.5), Is.EqualTo("X"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "I", "L" }, DecipherSequenceType.AminoAcid), Is.EqualTo("J"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "I", "L" }, DecipherSequenceType.AminoAcid, 0.5), Is.EqualTo("J"));

            string[] dna = { "ANGCT-", "-ACCT-" };
            Assert.That(MotifFinder.GenerateDecipherConsensus(dna), Is.EqualTo("ANSCT-"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(dna, includeTerminalGaps: true), Is.EqualTo("+NSCT-"));

            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "ANQIH-", "ADELW." }, DecipherSequenceType.AminoAcid), Is.EqualTo("ABZJX-"));

            string[] nonLetters = { "A-+.A", "AAAAA" };
            Assert.That(MotifFinder.GenerateDecipherConsensus(nonLetters, noConsensusChar: 'N'), Is.EqualTo("ANNNA"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(nonLetters, includeNonLetters: true), Is.EqualTo("AAAAA"));

            string[] degenerate = { "AWNDA", "AAAAA" };
            Assert.That(MotifFinder.GenerateDecipherConsensus(degenerate), Is.EqualTo("AWNDA"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(degenerate, ambiguity: false), Is.EqualTo("AAAAA"));
        });
    }

    /// <summary>Further outputs of the R build: RNA writes U, ragged rows, lower case.</summary>
    [Test]
    public void GenerateDecipherConsensus_RnaRaggedAndCase_MatchDecipher()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "ACGU", "ACGU", "UCGA" }, DecipherSequenceType.Rna), Is.EqualTo("WCGW"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "ACGTAC", "AC", "ACG" }), Is.EqualTo("ACGTAC"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "ACGTAC", "AC", "ACG" }, includeTerminalGaps: true), Is.EqualTo("ACGTAC"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(new[] { "acd", "ACD" }, DecipherSequenceType.AminoAcid), Is.EqualTo("ACD"));
            Assert.That(MotifFinder.GenerateDecipherConsensus(Array.Empty<string>()), Is.EqualTo(""));
        });
    }

    /// <summary>
    /// Hand-derived from makeConsensus (t = 1 − threshold): A 0.6 / C 0.4 → no single base reaches
    /// 0.95, M (A+C = 1, both above G and T) does; A 0.9 / C 0.06 / G 0.04 → M is tested before the
    /// three-base codes (C 0.06 > G 0.04 = excluded maximum).
    /// </summary>
    [Test]
    public void GenerateDecipherConsensus_HandDerived_SourceOrder()
    {
        string[] sixFour = { "A", "A", "A", "C", "C", "A", "A", "C", "C", "A" };
        Assert.That(MotifFinder.GenerateDecipherConsensus(sixFour), Is.EqualTo("M"));

        var rows = Enumerable.Repeat("A", 45).Concat(Enumerable.Repeat("C", 3)).Concat(Enumerable.Repeat("G", 2)).ToArray();
        Assert.That(MotifFinder.GenerateDecipherConsensus(rows), Is.EqualTo("M"));
        // Threshold 0.1 → t = 0.9: A alone (0.9) passes.
        Assert.That(MotifFinder.GenerateDecipherConsensus(rows, threshold: 0.1), Is.EqualTo("A"));
    }

    [Test]
    public void GenerateDecipherConsensus_Guards_AsDecipherR()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.GenerateDecipherConsensus(null!));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "AC", null! }));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, threshold: 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, threshold: -0.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, threshold: double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, minInformation: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, minInformation: 1.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, (DecipherSequenceType)7));
            // noConsensusChar must be in the alphabet (pmatch against DNA_/RNA_/AA_ALPHABET).
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, noConsensusChar: 'Z'));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, DecipherSequenceType.Rna, noConsensusChar: 'T'));
            Assert.DoesNotThrow(() => MotifFinder.GenerateDecipherConsensus(new[] { "A" }, DecipherSequenceType.AminoAcid, noConsensusChar: '*'));
            // Characters outside the Biostrings alphabet.
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "ACZ" }));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "ACT" }, DecipherSequenceType.Rna));
            Assert.Throws<ArgumentException>(() => MotifFinder.GenerateDecipherConsensus(new[] { "AC~" }));
        });
    }
}
