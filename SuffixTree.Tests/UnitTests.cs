using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests
{
    [TestFixture]
    public class UnitTests
    {
        private const int RANDOM_SEED = 42; // Fixed seed for deterministic tests

        [TestCase("abcdefghijklmnpqrst")]
        [TestCase("aaaaaaaaaaaaaaaaaaaaa")]
        [TestCase("bbefefbbaaaaaaaaaaaabbefefbb")]
        [TestCase("bananasbanananananananananassssss")]
        [TestCase("olbafuynhfcxzqhnebecxjrfwfttwrxvgujqxaxuaukbflddcrptlvyoaxuwzlwmoeljnxgmsleapkyzodhtymxuvlchoomsuodicehnzyebqtgsqeplinthhnalituxrisknsyjszuaatwoulznpjbvjmhytqgaqmctqvwgxailhproehwctldlagpjqaawdbialginqmweqrcopiqfnludmjuxkqlsgrydzyhecoojgmspowoykgghnbudhujnmyhqxbkfggxxprgfhraksfylcveevxvlxpzxkcqtkchasarbusvqzimvvfsvredhjykpqyyysyxbzwsuqahpjcroqvhysaynfheehppinszvwmyqlmymyqngrqzuefojczpoqcgbkvkmfpipdoetqxtdigphjhkxuwzieqirlvapypdysohfydtxzppfuufcreorhpsyydvvvsproofmuucwqqtskzieegstlokqkvjbssfythoenpbhlhnnsgknlapaigdwvrvsnyrhxhuzqkzoakldexmvnuvqscxmrysnuumawqrldjbtbmnhytvmmyykdaxuvqifecczafafzewmuplebvkxseatwsxwatbszboybwzhgfdtsjpxckknalqvgwuwwretocfaphnyuoyvnxbtabosfewkfrlbbeiduuidlogxfdacbplkbkpljvthltjjrlxbtejpdqjddnnsfhsljjfvmsjigyxmhjeqfcmmzzqpsxmnkuwlhhvcrtxskfmyieoctweswpkplcnjiqmtjjdloobapntxqmducnkabjcutinyhhekioybfokektjerdojqfvyalkvpqsznlvqvrswhelvburtkzdcceqehyqndhlcvkbieceazmuanqiauhkyhcbcckeydaevunddkwlntezctepnfrchvquxgtsnupoiwneengszjggwxkmahlbiwzsbyryqasufdsaaigulgwjqepccwesmbcfpoymrsjrbqwzjpjmbexpjloxdtwxqbdmggreurdcohfpgbchhrthdopewrsyfindsvrexpkkooxkmzxklsalyfuxscwthbfdbeghnpowbjxcedzogidsrdnjimcybbxmwpdiwnihhgylpsbukpsjtbkktylouakffurdfmpsnndtjcvjkbviezyqdgvhdcllibfbniafffwebrmyvbryjnomzgiglecxjntcvcrngwrvhefqaswhpynyzqwdpvewmjlpndtihwebjqolymkytrtidajqrdyvqzhcsvlvfvqspskkttqjsotdqkcdwzmdxxuevpvcrsijxskruaajrqaqgcarbxfrwerhddeetidequujlxmyaaoriomkhdmqaitbzbvhmnhmuntueqwueagpomwdhturmpwkyszjiwwlucqbhqbxgibuqmghvlrrbypswfsxkhgwjcndjnqblxargeegkzmhlahbahsfecevnpbxqdbuamjffddctbcedlcptoynjiuypvbgeatatnxztxsxvjrihxmoeeqmghwxxdyzrczljthnteqrfrquhvlssswndmdwxcfzrhcszffqdnjmqyjnywrurbsyavdxcwwtjsttcbsnvrpgiqlswqdcqmxjxwoebxjwlhlxbjuxuacdwktlivrfmncnqosxecfccutmikgwkeprlrkdfcinqgeeeompsmpcvxvnopzmrnuvdljcxjurxmliveisyfqsnpxsokkefgdujosxckvrkgeavugntchvztxkdqeiwyluxxgptyuuligmgfjcwcynffbgysjewlaaglqjuujjxytrphnfwncbkgkwswhcvliseqyifouatvszslptxqnhawzjhgfyorphndgksqdeoqohsqvwctwofrvqqpsnfisbcpluhesurrihkxvpugeitmatignbqqqldkdwqzaggxmitqlzobbuqccoeddmsdtjvywnbiiwkbidkjrofmbxjlnzfryzgxjbwgiaxbahchovroigmraoofyuzqheonmrfpskgciitjtxjzbhlpsohvysrwdwviirlxpvemizykpykhipjwhmqxoiwtevhyddyrigooibzrshqmbypvthubgozvhinzmntadmkfplledvglacrbeghcofvsddhokjhyfcqwwhbwjlkafilmaezpwezzgzgajpxhxcgwmcieilzlfrsxjlagjbjryhbrznmsfushtydgfsizclunncsbzpktmkmhmacicjuqhqaozwtihtcokd")]
        public void Build_ShouldNotThrow(string treeContent)
        {
            Assert.DoesNotThrow(() => SuffixTree.Build(treeContent));
        }

        [Test]
        public void Build_WithNullString_ShouldThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SuffixTree.Build(null));
        }

        [Test]
        public void Build_WithNullCharacter_ShouldThrow()
        {
            Assert.Throws<ArgumentException>(() => SuffixTree.Build("abc\0def"));
        }

        [Test]
        public void Build_WithEmptyString_ShouldNotThrow()
        {
            Assert.DoesNotThrow(() => SuffixTree.Build(""));
        }

        [Test]
        public void Contains_EmptyString_ShouldReturnTrue()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.Contains(""), Is.True);
        }

        [Test]
        public void Contains_WithNullString_ShouldThrow()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.Contains(null));
        }

        [Test]
        public void Contains_SingleCharacter_ShouldWork()
        {
            var st = SuffixTree.Build("a");
            Assert.That(st.Contains("a"), Is.True);
            Assert.That(st.Contains("b"), Is.False);
        }

        // ==================== Contains with ReadOnlySpan Tests ====================

        [Test]
        public void Contains_Span_EmptySpan_ShouldReturnTrue()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.Contains(ReadOnlySpan<char>.Empty), Is.True);
        }

        [Test]
        public void Contains_Span_FromStringSlice_ShouldWork()
        {
            var st = SuffixTree.Build("hello world");
            var source = "say hello please";

            // Using span slice to avoid string allocation
            Assert.That(st.Contains(source.AsSpan(4, 5)), Is.True); // "hello"
            Assert.That(st.Contains(source.AsSpan(0, 3)), Is.False); // "say"
        }

        [Test]
        public void Contains_Span_FromCharArray_ShouldWork()
        {
            var st = SuffixTree.Build("abcdef");
            char[] pattern = { 'c', 'd', 'e' };

            Assert.That(st.Contains(pattern.AsSpan()), Is.True);
            Assert.That(st.Contains(pattern.AsSpan(0, 2)), Is.True); // "cd"
        }

        [Test]
        public void Contains_Span_MatchesStringOverload()
        {
            const int CYCLES = 100;
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, 100);
                var st = SuffixTree.Build(s);

                int pos = r.Next(0, s.Length - 5);
                int len = r.Next(1, Math.Min(10, s.Length - pos));
                var pattern = s.Substring(pos, len);

                var stringResult = st.Contains(pattern);
                var spanResult = st.Contains(pattern.AsSpan());

                Assert.That(spanResult, Is.EqualTo(stringResult),
                    $"Mismatch for pattern '{pattern}' in '{s}'");
            }
        }

        [Test]
        public void Contains_Span_WithTerminator_ThrowsArgumentException()
        {
            var st = SuffixTree.Build("hello");
            char[] pattern = { 'h', 'e', '\0', 'l' };

            Assert.Throws<ArgumentException>(() => st.Contains(pattern.AsSpan()));
        }

        [Test]
        public void Contains_OnNonRepeated_ShouldReturnTrue()
        {
            var s = "abcdefghijklmnpqrst";
            var t = s.Substring(3, 6);
            var st = SuffixTree.Build(s);

            Assert.That(st.Contains(t), Is.True);
        }

        [Test]
        public void Contains_OnNonRepeated_ShouldReturnFalse()
        {
            var s = "abcdefghijklmnpqrst";
            var t = "rstb";
            var st = SuffixTree.Build(s);

            Assert.That(st.Contains(t), Is.False);
        }

        [Test]
        public void Contains_Simple_OnLongString_ShouldReturnTrue()
        {
            var s = "olbafuynhfcxzqhnebecxjrfwfttwrxvgujqxaxuaukbflddcrptlvyoaxuwzlwmoeljnxgmsleapkyzodhtymxuvlchoomsuodicehnzyebqtgsqeplinthhnalituxrisknsyjszuaatwoulznpjbvjmhytqgaqmctqvwgxailhproehwctldlagpjqaawdbialginqmweqrcopiqfnludmjuxkqlsgrydzyhecoojgmspowoykgghnbudhujnmyhqxbkfggxxprgfhraksfylcveevxvlxpzxkcqtkchasarbusvqzimvvfsvredhjykpqyyysyxbzwsuqahpjcroqvhysaynfheehppinszvwmyqlmymyqngrqzuefojczpoqcgbkvkmfpipdoetqxtdigphjhkxuwzieqirlvapypdysohfydtxzppfuufcreorhpsyydvvvsproofmuucwqqtskzieegstlokqkvjbssfythoenpbhlhnnsgknlapaigdwvrvsnyrhxhuzqkzoakldexmvnuvqscxmrysnuumawqrldjbtbmnhytvmmyykdaxuvqifecczafafzewmuplebvkxseatwsxwatbszboybwzhgfdtsjpxckknalqvgwuwwretocfaphnyuoyvnxbtabosfewkfrlbbeiduuidlogxfdacbplkbkpljvthltjjrlxbtejpdqjddnnsfhsljjfvmsjigyxmhjeqfcmmzzqpsxmnkuwlhhvcrtxskfmyieoctweswpkplcnjiqmtjjdloobapntxqmducnkabjcutinyhhekioybfokektjerdojqfvyalkvpqsznlvqvrswhelvburtkzdcceqehyqndhlcvkbieceazmuanqiauhkyhcbcckeydaevunddkwlntezctepnfrchvquxgtsnupoiwneengszjggwxkmahlbiwzsbyryqasufdsaaigulgwjqepccwesmbcfpoymrsjrbqwzjpjmbexpjloxdtwxqbdmggreurdcohfpgbchhrthdopewrsyfindsvrexpkkooxkmzxklsalyfuxscwthbfdbeghnpowbjxcedzogidsrdnjimcybbxmwpdiwnihhgylpsbukpsjtbkktylouakffurdfmpsnndtjcvjkbviezyqdgvhdcllibfbniafffwebrmyvbryjnomzgiglecxjntcvcrngwrvhefqaswhpynyzqwdpvewmjlpndtihwebjqolymkytrtidajqrdyvqzhcsvlvfvqspskkttqjsotdqkcdwzmdxxuevpvcrsijxskruaajrqaqgcarbxfrwerhddeetidequujlxmyaaoriomkhdmqaitbzbvhmnhmuntueqwueagpomwdhturmpwkyszjiwwlucqbhqbxgibuqmghvlrrbypswfsxkhgwjcndjnqblxargeegkzmhlahbahsfecevnpbxqdbuamjffddctbcedlcptoynjiuypvbgeatatnxztxsxvjrihxmoeeqmghwxxdyzrczljthnteqrfrquhvlssswndmdwxcfzrhcszffqdnjmqyjnywrurbsyavdxcwwtjsttcbsnvrpgiqlswqdcqmxjxwoebxjwlhlxbjuxuacdwktlivrfmncnqosxecfccutmikgwkeprlrkdfcinqgeeeompsmpcvxvnopzmrnuvdljcxjurxmliveisyfqsnpxsokkefgdujosxckvrkgeavugntchvztxkdqeiwyluxxgptyuuligmgfjcwcynffbgysjewlaaglqjuujjxytrphnfwncbkgkwswhcvliseqyifouatvszslptxqnhawzjhgfyorphndgksqdeoqohsqvwctwofrvqqpsnfisbcpluhesurrihkxvpugeitmatignbqqqldkdwqzaggxmitqlzobbuqccoeddmsdtjvywnbiiwkbidkjrofmbxjlnzfryzgxjbwgiaxbahchovroigmraoofyuzqheonmrfpskgciitjtxjzbhlpsohvysrwdwviirlxpvemizykpykhipjwhmqxoiwtevhyddyrigooibzrshqmbypvthubgozvhinzmntadmkfplledvglacrbeghcofvsddhokjhyfcqwwhbwjlkafilmaezpwezzgzgajpxhxcgwmcieilzlfrsxjlagjbjryhbrznmsfushtydgfsizclunncsbzpktmkmhmacicjuqhqaozwtihtcokd";
            var t = s.Substring(s.Length / 2, s.Length / 50);
            var st = SuffixTree.Build(s);

            Assert.That(st.Contains(t), Is.True);
        }

        [Test]
        public void Contains_DynamicBarrage_OnLongString_ShouldReturnTrue()
        {
            const int CYCLES = 200;
            const int MAXLEN = 200;

            var s = "olbafuynhfcxzqhnebecxjrfwfttwrxvgujqxaxuaukbflddcrptlvyoaxuwzlwmoeljnxgmsleapkyzodhtymxuvlchoomsuodicehnzyebqtgsqeplinthhnalituxrisknsyjszuaatwoulznpjbvjmhytqgaqmctqvwgxailhproehwctldlagpjqaawdbialginqmweqrcopiqfnludmjuxkqlsgrydzyhecoojgmspowoykgghnbudhujnmyhqxbkfggxxprgfhraksfylcveevxvlxpzxkcqtkchasarbusvqzimvvfsvredhjykpqyyysyxbzwsuqahpjcroqvhysaynfheehppinszvwmyqlmymyqngrqzuefojczpoqcgbkvkmfpipdoetqxtdigphjhkxuwzieqirlvapypdysohfydtxzppfuufcreorhpsyydvvvsproofmuucwqqtskzieegstlokqkvjbssfythoenpbhlhnnsgknlapaigdwvrvsnyrhxhuzqkzoakldexmvnuvqscxmrysnuumawqrldjbtbmnhytvmmyykdaxuvqifecczafafzewmuplebvkxseatwsxwatbszboybwzhgfdtsjpxckknalqvgwuwwretocfaphnyuoyvnxbtabosfewkfrlbbeiduuidlogxfdacbplkbkpljvthltjjrlxbtejpdqjddnnsfhsljjfvmsjigyxmhjeqfcmmzzqpsxmnkuwlhhvcrtxskfmyieoctweswpkplcnjiqmtjjdloobapntxqmducnkabjcutinyhhekioybfokektjerdojqfvyalkvpqsznlvqvrswhelvburtkzdcceqehyqndhlcvkbieceazmuanqiauhkyhcbcckeydaevunddkwlntezctepnfrchvquxgtsnupoiwneengszjggwxkmahlbiwzsbyryqasufdsaaigulgwjqepccwesmbcfpoymrsjrbqwzjpjmbexpjloxdtwxqbdmggreurdcohfpgbchhrthdopewrsyfindsvrexpkkooxkmzxklsalyfuxscwthbfdbeghnpowbjxcedzogidsrdnjimcybbxmwpdiwnihhgylpsbukpsjtbkktylouakffurdfmpsnndtjcvjkbviezyqdgvhdcllibfbniafffwebrmyvbryjnomzgiglecxjntcvcrngwrvhefqaswhpynyzqwdpvewmjlpndtihwebjqolymkytrtidajqrdyvqzhcsvlvfvqspskkttqjsotdqkcdwzmdxxuevpvcrsijxskruaajrqaqgcarbxfrwerhddeetidequujlxmyaaoriomkhdmqaitbzbvhmnhmuntueqwueagpomwdhturmpwkyszjiwwlucqbhqbxgibuqmghvlrrbypswfsxkhgwjcndjnqblxargeegkzmhlahbahsfecevnpbxqdbuamjffddctbcedlcptoynjiuypvbgeatatnxztxsxvjrihxmoeeqmghwxxdyzrczljthnteqrfrquhvlssswndmdwxcfzrhcszffqdnjmqyjnywrurbsyavdxcwwtjsttcbsnvrpgiqlswqdcqmxjxwoebxjwlhlxbjuxuacdwktlivrfmncnqosxecfccutmikgwkeprlrkdfcinqgeeeompsmpcvxvnopzmrnuvdljcxjurxmliveisyfqsnpxsokkefgdujosxckvrkgeavugntchvztxkdqeiwyluxxgptyuuligmgfjcwcynffbgysjewlaaglqjuujjxytrphnfwncbkgkwswhcvliseqyifouatvszslptxqnhawzjhgfyorphndgksqdeoqohsqvwctwofrvqqpsnfisbcpluhesurrihkxvpugeitmatignbqqqldkdwqzaggxmitqlzobbuqccoeddmsdtjvywnbiiwkbidkjrofmbxjlnzfryzgxjbwgiaxbahchovroigmraoofyuzqheonmrfpskgciitjtxjzbhlpsohvysrwdwviirlxpvemizykpykhipjwhmqxoiwtevhyddyrigooibzrshqmbypvthubgozvhinzmntadmkfplledvglacrbeghcofvsddhokjhyfcqwwhbwjlkafilmaezpwezzgzgajpxhxcgwmcieilzlfrsxjlagjbjryhbrznmsfushtydgfsizclunncsbzpktmkmhmacicjuqhqaozwtihtcokd";
            var st = SuffixTree.Build(s);
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var pos = r.Next(0, s.Length - 2);
                var len = r.Next(1, Math.Min(s.Length - pos, MAXLEN));
                var ss = s.Substring(pos, len);

                Assert.That(st.Contains(ss), Is.True, $"Failed for substring at pos {pos}, len {len}: '{ss}'");
            }
        }

        [Test]
        public void Contains_InShortRandomString_ShouldReturnTrue()
        {
            const int CYCLES = 500;
            const int STRLEN = 200;
            const int MATCHLEN = 20;

            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, STRLEN);
                var pos = r.Next(0, s.Length - MATCHLEN);
                var ss = s.Substring(pos, MATCHLEN);

                bool res;
                try
                {
                    res = SuffixTree.Build(s).Contains(ss);
                }
                catch (Exception e)
                {
                    throw new Exception($"Cycle {i}: String: '{s}', Substring: '{ss}'", e);
                }

                Assert.That(res, Is.True, $"Cycle {i}: String: '{s}', Substring: '{ss}'");
            }
        }

        [Test]
        public void Contains_RepeatingPattern_ShouldWork()
        {
            var st = SuffixTree.Build("abababab");
            Assert.That(st.Contains("abab"), Is.True);
            Assert.That(st.Contains("baba"), Is.True);
            Assert.That(st.Contains("ab"), Is.True);
            Assert.That(st.Contains("ababa bab"), Is.False);
        }

        [Test]
        public void Contains_FullString_ShouldReturnTrue()
        {
            var s = "teststring";
            var st = SuffixTree.Build(s);
            Assert.That(st.Contains(s), Is.True);
        }

        [Test]
        public void ToString_ShouldReturnMeaningfulString()
        {
            var st = SuffixTree.Build("hello");
            var result = st.ToString();

            Assert.That(result, Does.Contain("SuffixTree"));
            Assert.That(result, Does.Contain("hello"));
        }

        [Test]
        public void PrintTree_ShouldNotThrow()
        {
            var st = SuffixTree.Build("banana");
            Assert.DoesNotThrow(() => st.PrintTree());
        }

        [Test]
        public void Contains_AllSuffixes_ShouldReturnTrue()
        {
            // Verify that all suffixes of the string are found
            var s = "banana";
            var st = SuffixTree.Build(s);

            for (int i = 0; i < s.Length; i++)
            {
                var suffix = s.Substring(i);
                Assert.That(st.Contains(suffix), Is.True, $"Suffix '{suffix}' not found");
            }
        }

        [Test]
        public void Contains_AllSubstrings_ShouldReturnTrue()
        {
            // Verify that ALL possible substrings are found
            var s = "abracadabra";
            var st = SuffixTree.Build(s);

            for (int i = 0; i < s.Length; i++)
            {
                for (int len = 1; len <= s.Length - i; len++)
                {
                    var substr = s.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True, $"Substring '{substr}' at [{i},{i + len}) not found");
                }
            }
        }

        [Test]
        public void Contains_NonExistentSubstrings_ShouldReturnFalse()
        {
            var st = SuffixTree.Build("abcdef");

            Assert.That(st.Contains("xyz"), Is.False);
            Assert.That(st.Contains("abd"), Is.False);  // 'd' follows 'c', not 'b'
            Assert.That(st.Contains("abcdefg"), Is.False);  // longer than original
            Assert.That(st.Contains("aba"), Is.False);
            Assert.That(st.Contains("fgh"), Is.False);
        }

        [Test]
        public void Build_WithRepeatingCharacter_ShouldWork()
        {
            var st = SuffixTree.Build("aaaaaaaaaa");
            Assert.That(st.Contains("aaa"), Is.True);
            Assert.That(st.Contains("aaaaaaaaaa"), Is.True);
            Assert.That(st.Contains("b"), Is.False);
        }

        [Test]
        public void Build_Mississipi_ClassicTest()
        {
            // Classic test case from Ukkonen's paper
            var st = SuffixTree.Build("mississippi");

            // All suffixes
            Assert.That(st.Contains("mississippi"), Is.True);
            Assert.That(st.Contains("ississippi"), Is.True);
            Assert.That(st.Contains("ssissippi"), Is.True);
            Assert.That(st.Contains("sissippi"), Is.True);
            Assert.That(st.Contains("issippi"), Is.True);
            Assert.That(st.Contains("ssippi"), Is.True);
            Assert.That(st.Contains("sippi"), Is.True);
            Assert.That(st.Contains("ippi"), Is.True);
            Assert.That(st.Contains("ppi"), Is.True);
            Assert.That(st.Contains("pi"), Is.True);
            Assert.That(st.Contains("i"), Is.True);

            // Some substrings
            Assert.That(st.Contains("issi"), Is.True);
            Assert.That(st.Contains("sis"), Is.True);
            Assert.That(st.Contains("pp"), Is.True);

            // Non-existent
            Assert.That(st.Contains("spa"), Is.False);
            Assert.That(st.Contains("mississippii"), Is.False);
        }

        [Test]
        public void Contains_EdgeCaseBoundary_ShouldWork()
        {
            // Test matching at edge boundaries
            var st = SuffixTree.Build("abcabxabcd");

            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("abca"), Is.True);
            Assert.That(st.Contains("abcab"), Is.True);
            Assert.That(st.Contains("abx"), Is.True);
            Assert.That(st.Contains("xabc"), Is.True);
            Assert.That(st.Contains("abcd"), Is.True);
            Assert.That(st.Contains("bcd"), Is.True);
        }

        // ==================== FindAllOccurrences Tests ====================

        [Test]
        public void FindAllOccurrences_WithNullPattern_ShouldThrow()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.FindAllOccurrences(null));
        }

        [Test]
        public void FindAllOccurrences_EmptyPattern_ReturnsAllPositions()
        {
            var st = SuffixTree.Build("abc");
            var result = st.FindAllOccurrences("");
            Assert.That(result.Count, Is.EqualTo(3)); // positions 0,1,2 (text length)
            Assert.That(result, Is.EquivalentTo(new[] { 0, 1, 2 }));
        }

        [Test]
        public void FindAllOccurrences_NotFound_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abcdef");
            var result = st.FindAllOccurrences("xyz");
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void FindAllOccurrences_SingleOccurrence_ReturnsCorrectPosition()
        {
            var st = SuffixTree.Build("abcdef");
            var result = st.FindAllOccurrences("cde");
            Assert.That(result, Is.EquivalentTo(new[] { 2 }));
        }

        [Test]
        public void FindAllOccurrences_MultipleOccurrences_ReturnsAll()
        {
            var st = SuffixTree.Build("abcabc");
            var result = st.FindAllOccurrences("abc");
            Assert.That(result, Is.EquivalentTo(new[] { 0, 3 }));
        }

        [Test]
        public void FindAllOccurrences_OverlappingOccurrences_ReturnsAll()
        {
            var st = SuffixTree.Build("aaa");
            var result = st.FindAllOccurrences("aa");
            Assert.That(result, Is.EquivalentTo(new[] { 0, 1 }));
        }

        [Test]
        public void FindAllOccurrences_FullString_ReturnsSinglePosition()
        {
            var st = SuffixTree.Build("hello");
            var result = st.FindAllOccurrences("hello");
            Assert.That(result, Is.EquivalentTo(new[] { 0 }));
        }

        [Test]
        public void FindAllOccurrences_SingleChar_ReturnsAllPositions()
        {
            var st = SuffixTree.Build("abacada");
            var result = st.FindAllOccurrences("a");
            Assert.That(result, Is.EquivalentTo(new[] { 0, 2, 4, 6 }));
        }

        [Test]
        public void FindAllOccurrences_Banana_FindsAllAna()
        {
            var st = SuffixTree.Build("banana");
            var result = st.FindAllOccurrences("ana");
            Assert.That(result, Is.EquivalentTo(new[] { 1, 3 }));
        }

        [Test]
        public void FindAllOccurrences_Mississippi_FindsAllIssi()
        {
            var st = SuffixTree.Build("mississippi");
            var result = st.FindAllOccurrences("issi");
            Assert.That(result, Is.EquivalentTo(new[] { 1, 4 }));
        }

        [Test]
        public void FindAllOccurrences_Stress_VerifyAgainstNaive()
        {
            const int CYCLES = 100;
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, 100);
                var st = SuffixTree.Build(s);

                // Pick a random substring to search
                int pos = r.Next(0, s.Length - 5);
                int len = r.Next(1, Math.Min(10, s.Length - pos));
                var pattern = s.Substring(pos, len);

                var stResult = st.FindAllOccurrences(pattern);
                var naiveResult = NaiveFindAll(s, pattern);

                Assert.That(stResult.OrderBy(x => x), Is.EqualTo(naiveResult.OrderBy(x => x)),
                    $"Mismatch for string '{s}' pattern '{pattern}'");
            }
        }

        private static List<int> NaiveFindAll(string text, string pattern)
        {
            var result = new List<int>();
            int idx = 0;
            while ((idx = text.IndexOf(pattern, idx)) != -1)
            {
                result.Add(idx);
                idx++;
            }
            return result;
        }

        // ==================== CountOccurrences Tests ====================

        [Test]
        public void CountOccurrences_WithNullPattern_ShouldThrow()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.CountOccurrences(null));
        }

        [Test]
        public void CountOccurrences_EmptyPattern_ReturnsTextLength()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.CountOccurrences(""), Is.EqualTo(3)); // positions 0,1,2 (excluding terminator)
        }

        [Test]
        public void CountOccurrences_NotFound_ReturnsZero()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.CountOccurrences("xyz"), Is.EqualTo(0));
        }

        [Test]
        public void CountOccurrences_SingleOccurrence_ReturnsOne()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.CountOccurrences("cde"), Is.EqualTo(1));
        }

        [Test]
        public void CountOccurrences_MultipleOccurrences_ReturnsCorrectCount()
        {
            var st = SuffixTree.Build("abcabc");
            Assert.That(st.CountOccurrences("abc"), Is.EqualTo(2));
        }

        [Test]
        public void CountOccurrences_OverlappingOccurrences_ReturnsCorrectCount()
        {
            var st = SuffixTree.Build("aaaa");
            Assert.That(st.CountOccurrences("aa"), Is.EqualTo(3)); // positions 0,1,2
        }

        [Test]
        public void CountOccurrences_Banana_Ana()
        {
            var st = SuffixTree.Build("banana");
            Assert.That(st.CountOccurrences("ana"), Is.EqualTo(2));
        }

        [Test]
        public void CountOccurrences_Mississippi_Issi()
        {
            var st = SuffixTree.Build("mississippi");
            Assert.That(st.CountOccurrences("issi"), Is.EqualTo(2));
        }

        [Test]
        public void CountOccurrences_MatchesFindAllOccurrences()
        {
            const int CYCLES = 100;
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, 100);
                var st = SuffixTree.Build(s);

                int pos = r.Next(0, s.Length - 5);
                int len = r.Next(1, Math.Min(10, s.Length - pos));
                var pattern = s.Substring(pos, len);

                var findAllCount = st.FindAllOccurrences(pattern).Count;
                var countResult = st.CountOccurrences(pattern);

                Assert.That(countResult, Is.EqualTo(findAllCount),
                    $"Mismatch for string '{s}' pattern '{pattern}'");
            }
        }

        // ==================== LongestRepeatedSubstring Tests ====================

        [Test]
        public void LongestRepeatedSubstring_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_SingleChar_ReturnsEmpty()
        {
            var st = SuffixTree.Build("a");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_AllUnique_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_SimpleRepeat()
        {
            var st = SuffixTree.Build("abcabc");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("abc"));
        }

        [Test]
        public void LongestRepeatedSubstring_OverlappingRepeat()
        {
            var st = SuffixTree.Build("aaa");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("aa"));
        }

        [Test]
        public void LongestRepeatedSubstring_Banana()
        {
            var st = SuffixTree.Build("banana");
            // "ana" appears twice (positions 1 and 3)
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("ana"));
        }

        [Test]
        public void LongestRepeatedSubstring_Mississippi()
        {
            var st = SuffixTree.Build("mississippi");
            // "issi" appears twice
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("issi"));
        }

        [Test]
        public void LongestRepeatedSubstring_ABAB()
        {
            var st = SuffixTree.Build("ABAB");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("AB"));
        }

        [Test]
        public void LongestRepeatedSubstring_GEEKSFORGEEKS()
        {
            var st = SuffixTree.Build("GEEKSFORGEEKS");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("GEEKS"));
        }

        [Test]
        public void LongestRepeatedSubstring_ResultAppearsAtLeastTwice()
        {
            const int CYCLES = 50;
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, 100);
                var st = SuffixTree.Build(s);
                var lrs = st.LongestRepeatedSubstring();

                if (lrs.Length > 0)
                {
                    var count = st.CountOccurrences(lrs);
                    Assert.That(count, Is.GreaterThanOrEqualTo(2),
                        $"LRS '{lrs}' should appear at least twice in '{s}'");
                }
            }
        }

        [Test]
        public void LongestRepeatedSubstring_IsActuallyLongest()
        {
            const int CYCLES = 30;
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < CYCLES; i++)
            {
                var s = MakeRandomString(r, 50);
                var st = SuffixTree.Build(s);
                var lrs = st.LongestRepeatedSubstring();

                // Verify no longer repeated substring exists
                if (lrs.Length > 0 && lrs.Length < s.Length - 1)
                {
                    // Check that lrs+1 characters don't repeat
                    // This is a probabilistic check - we try a few extensions
                    bool foundLonger = false;
                    for (int j = 0; j < s.Length - lrs.Length && !foundLonger; j++)
                    {
                        var candidate = s.Substring(j, lrs.Length + 1);
                        if (st.CountOccurrences(candidate) >= 2)
                            foundLonger = true;
                    }

                    if (foundLonger)
                    {
                        Assert.Fail($"Found longer repeated substring than LRS '{lrs}' in '{s}'");
                    }
                }
            }
        }

        #region GetAllSuffixes Tests

        [Test]
        public void GetAllSuffixes_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.GetAllSuffixes(), Is.Empty);
        }

        [Test]
        public void GetAllSuffixes_SingleChar_ReturnsSingleSuffix()
        {
            var st = SuffixTree.Build("a");
            var suffixes = st.GetAllSuffixes();
            Assert.That(suffixes, Has.Count.EqualTo(1));
            Assert.That(suffixes[0], Is.EqualTo("a"));
        }

        [Test]
        public void GetAllSuffixes_Banana_ReturnsSortedSuffixes()
        {
            var st = SuffixTree.Build("banana");
            var suffixes = st.GetAllSuffixes();

            // All suffixes of "banana": a, ana, anana, banana, na, nana
            var expected = new[] { "a", "ana", "anana", "banana", "na", "nana" };
            Assert.That(suffixes, Is.EqualTo(expected));
        }

        [Test]
        public void GetAllSuffixes_CountMatchesStringLength()
        {
            const string text = "abracadabra";
            var st = SuffixTree.Build(text);
            var suffixes = st.GetAllSuffixes();

            Assert.That(suffixes.Count, Is.EqualTo(text.Length));
        }

        [Test]
        public void GetAllSuffixes_AreSortedLexicographically()
        {
            var st = SuffixTree.Build("mississippi");
            var suffixes = st.GetAllSuffixes();

            var sorted = suffixes.OrderBy(s => s).ToList();
            Assert.That(suffixes, Is.EqualTo(sorted));
        }

        #endregion

        #region LongestCommonSubstring Tests

        [Test]
        public void LongestCommonSubstring_Null_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.LongestCommonSubstring(null));
        }

        [Test]
        public void LongestCommonSubstring_EmptyOther_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.LongestCommonSubstring(""), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.LongestCommonSubstring("abc"), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_NoCommon_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.LongestCommonSubstring("xyz"), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_FullMatch()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.LongestCommonSubstring("abcdef"), Is.EqualTo("abcdef"));
        }

        [Test]
        public void LongestCommonSubstring_PartialMatch()
        {
            var st = SuffixTree.Build("xyzabcdef");
            var lcs = st.LongestCommonSubstring("123abc456");
            Assert.That(lcs, Is.EqualTo("abc"));
        }

        [Test]
        public void LongestCommonSubstring_MultipleMatches_ReturnsLongest()
        {
            var st = SuffixTree.Build("abcdefghijklmnop");
            var lcs = st.LongestCommonSubstring("xxabcxxdefghixx");
            Assert.That(lcs, Is.EqualTo("defghi"));
        }

        [Test]
        public void LongestCommonSubstring_SingleCharMatch()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.LongestCommonSubstring("xay"), Is.EqualTo("a"));
        }

        #endregion

        #region Stress Tests

        [Test]
        public void Contains_SuffixLinkTraversal_Stress()
        {
            var s = string.Concat(Enumerable.Repeat("abcdefgh", 100));
            var st = SuffixTree.Build(s);

            Assert.That(st.Contains("abcdefghabcdefgh"), Is.True);
            Assert.That(st.Contains("efghabcd"), Is.True);
            Assert.That(st.Contains("ghabcdefghabcdef"), Is.True);
            Assert.That(st.Contains(s), Is.True); // Full string
            Assert.That(st.Contains("xyz"), Is.False);
        }

        [Test]
        public void StressTest_BinaryAlphabet()
        {
            // Binary alphabet stress test
            var r = new Random(777);
            var chars = new char[500];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = r.Next(2) == 0 ? 'a' : 'b';

            var s = new string(chars);
            var st = SuffixTree.Build(s);

            // Verify all substrings of length 1-20
            for (int i = 0; i < s.Length - 20; i += 10)
            {
                for (int len = 1; len <= 20; len++)
                {
                    var substr = s.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True,
                        $"Substring '{substr}' at pos {i} not found");
                }
            }
        }

        [Test]
        public void Build_AbcAbxAbc_TrickyCase()
        {
            // This is a known tricky case for Ukkonen's algorithm
            // The split + suffix link handling must be correct
            var st = SuffixTree.Build("abcabxabcd");

            // All suffixes
            string[] suffixes = {
                "abcabxabcd",
                "bcabxabcd",
                "cabxabcd",
                "abxabcd",
                "bxabcd",
                "xabcd",
                "abcd",
                "bcd",
                "cd",
                "d"
            };

            foreach (var suffix in suffixes)
            {
                Assert.That(st.Contains(suffix), Is.True, $"Suffix '{suffix}' not found");
            }

            // Some substrings
            Assert.That(st.Contains("abcab"), Is.True);
            Assert.That(st.Contains("bcabx"), Is.True);
            Assert.That(st.Contains("xab"), Is.True);
            Assert.That(st.Contains("abc"), Is.True);
        }

        [Test]
        public void Contains_PartialMatchAtEdgeBoundary()
        {
            // Test case where match ends exactly at edge boundary
            var st = SuffixTree.Build("aaaaab");

            Assert.That(st.Contains("aaaa"), Is.True);
            Assert.That(st.Contains("aaaaa"), Is.True);
            Assert.That(st.Contains("aaaaab"), Is.True);
            Assert.That(st.Contains("aaaab"), Is.True);
            Assert.That(st.Contains("aaab"), Is.True);
            Assert.That(st.Contains("aab"), Is.True);
            Assert.That(st.Contains("ab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);

            // Non-existent
            Assert.That(st.Contains("ba"), Is.False);
            Assert.That(st.Contains("aaaaaba"), Is.False);
        }

        [Test]
        public void Comprehensive_AllSubstringsUpToLength50()
        {
            // Ultimate test: verify ALL substrings for a medium string
            var r = new Random(999);
            const string alphabet = "abcd";
            var chars = new char[50];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = alphabet[r.Next(alphabet.Length)];

            var s = new string(chars);
            var st = SuffixTree.Build(s);

            int substringCount = 0;
            for (int i = 0; i < s.Length; i++)
            {
                for (int len = 1; len <= s.Length - i; len++)
                {
                    var substr = s.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True,
                        $"String '{s}': Substring '{substr}' at [{i},{i + len}) not found");
                    substringCount++;
                }
            }

            // Sanity check: we tested n*(n+1)/2 substrings
            Assert.That(substringCount, Is.EqualTo(50 * 51 / 2));
        }

        private static string MakeRandomString(Random r, int len)
        {
            const string SET = "abcdefghijklmnopqrstuvwxyz";
            var res = new char[len];
            for (int i = 0; i < len; i++)
                res[i] = SET[r.Next(SET.Length)];

            return new string(res);
        }

        #endregion

        #region Tests for Previous Fixes

        // ===========================================
        // FIX 1: Terminator character added
        // The tree now uses '\0' as terminator to ensure proper suffix tree
        // ===========================================

        [Test]
        public void Fix_TerminatorCharacter_IsPartOfTreeStructure()
        {
            // Terminator '\0' is part of the tree, making all suffixes explicit
            // This ensures the suffix tree (not suffix trie) is correctly formed
            var st = SuffixTree.Build("abc");

            // The terminator is in the tree structure but users typically
            // don't search for it - they search for actual content
            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("bc"), Is.True);
            Assert.That(st.Contains("c"), Is.True);
        }

        [Test]
        public void Fix_TerminatorCharacter_EnsuresProperSuffixTree()
        {
            // Without terminator, "ab" in "aab" would be implicit suffix
            // With terminator, all suffixes are explicit (end at leaf)
            var st = SuffixTree.Build("aab");

            // All suffixes must be found
            Assert.That(st.Contains("aab"), Is.True);
            Assert.That(st.Contains("ab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);

            // The terminator ensures "ab" is not confused with "aab" prefix
            Assert.That(st.Contains("a"), Is.True);
            Assert.That(st.Contains("aa"), Is.True);
        }

        // ===========================================
        // FIX 2: Input validation added
        // ===========================================

        [Test]
        public void Fix_InputValidation_NullThrowsArgumentNullException()
        {
            var ex = Assert.Throws<ArgumentNullException>(() => SuffixTree.Build(null!));
            Assert.That(ex, Is.Not.Null);
            // Parameter name comes from nameof() - just verify exception is thrown
        }

        [Test]
        public void Fix_InputValidation_EmptyStringReturnsValidTree()
        {
            var st = SuffixTree.Build("");

            Assert.That(st, Is.Not.Null);
            Assert.That(st.Contains(""), Is.True, "Empty string should be found in any tree");
            Assert.That(st.Contains("a"), Is.False);
        }

        [Test]
        public void Fix_InputValidation_ContainsNullThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("test");

            var ex = Assert.Throws<ArgumentNullException>(() => st.Contains(null!));
            Assert.That(ex, Is.Not.Null);
            // Parameter name comes from nameof() - just verify exception is thrown
        }

        // ===========================================
        // FIX 3: State reset - multiple Build calls
        // Ensures each Build creates independent tree
        // ===========================================

        [Test]
        public void Fix_StateReset_MultipleBuildCallsAreIndependent()
        {
            // First tree
            var st1 = SuffixTree.Build("hello");
            Assert.That(st1.Contains("hello"), Is.True);
            Assert.That(st1.Contains("world"), Is.False);

            // Second tree - should be completely independent
            var st2 = SuffixTree.Build("world");
            Assert.That(st2.Contains("world"), Is.True);
            Assert.That(st2.Contains("hello"), Is.False);

            // Original tree unchanged
            Assert.That(st1.Contains("hello"), Is.True);
            Assert.That(st1.Contains("world"), Is.False);
        }

        [Test]
        public void Fix_StateReset_BuildResetsActivePoint()
        {
            // Build multiple trees to ensure active point doesn't leak
            for (int i = 0; i < 10; i++)
            {
                var s = new string((char)('a' + i), 5) + "xyz";
                var st = SuffixTree.Build(s);

                Assert.That(st.Contains(s), Is.True, $"Iteration {i}: full string not found");
                Assert.That(st.Contains("xyz"), Is.True, $"Iteration {i}: suffix 'xyz' not found");
            }
        }

        // ===========================================
        // FIX 4: Suffix links stored in Node, not external dictionary
        // ===========================================

        [Test]
        public void Fix_SuffixLinks_ComplexPatternRequiresSuffixLinks()
        {
            // "abcabxabcd" is a classic case requiring correct suffix link handling
            // When splitting "abc" edge at position 2 (to insert "x"), 
            // the new internal node needs suffix link to node after "bc"
            var st = SuffixTree.Build("abcabxabcd");

            // These substrings require proper suffix link traversal during construction
            Assert.That(st.Contains("abcabx"), Is.True);
            Assert.That(st.Contains("bcabxa"), Is.True);
            Assert.That(st.Contains("cabxab"), Is.True);
            Assert.That(st.Contains("abxabc"), Is.True);
        }

        [Test]
        public void Fix_SuffixLinks_RepeatingPatternsStress()
        {
            // Repeating patterns heavily exercise suffix links
            var st = SuffixTree.Build("abababababab");

            // All these require correct suffix link chain
            Assert.That(st.Contains("ababab"), Is.True);
            Assert.That(st.Contains("bababa"), Is.True);
            Assert.That(st.Contains("ababababab"), Is.True);
            Assert.That(st.Contains("abababababab"), Is.True);
        }

        [Test]
        public void Fix_SuffixLinks_MississippiClassicCase()
        {
            // Mississippi is the classic suffix tree test case
            var st = SuffixTree.Build("mississippi");

            // "issi" appears twice - suffix links must handle this
            Assert.That(st.Contains("issi"), Is.True);
            Assert.That(st.Contains("iss"), Is.True);
            Assert.That(st.Contains("ssi"), Is.True);
            Assert.That(st.Contains("si"), Is.True);

            // These require navigating through suffix links
            Assert.That(st.Contains("ississi"), Is.True);
            Assert.That(st.Contains("sissipp"), Is.True);
        }

        // ===========================================
        // FIX 5: Children moved to Node.Children dictionary
        // Original: external dictionary with (Node, char) tuple keys
        // ===========================================

        [Test]
        public void Fix_NodeChildren_LargeAlphabetEfficiency()
        {
            // With dictionary per node, large alphabet should work efficiently
            // Original tuple-key design was O(total_edges) per lookup
            var alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
            var st = SuffixTree.Build(alphabet);

            // Each character should be findable
            foreach (char c in alphabet)
            {
                Assert.That(st.Contains(c.ToString()), Is.True, $"Character '{c}' not found");
            }

            // Full string
            Assert.That(st.Contains(alphabet), Is.True);
        }

        // ===========================================
        // FIX 6: PrintTree now works with any characters, not just a-z
        // ===========================================

        [Test]
        public void Fix_PrintTree_UnicodeCharacters()
        {
            var st = SuffixTree.Build("привет");

            // Should not throw and should find all suffixes
            Assert.DoesNotThrow(() => st.PrintTree());
            Assert.That(st.Contains("привет"), Is.True);
            Assert.That(st.Contains("ривет"), Is.True);
            Assert.That(st.Contains("ивет"), Is.True);
            Assert.That(st.Contains("вет"), Is.True);
            Assert.That(st.Contains("ет"), Is.True);
            Assert.That(st.Contains("т"), Is.True);
        }

        [Test]
        public void Fix_PrintTree_SpecialCharacters()
        {
            var st = SuffixTree.Build("a!@#$%^&*()b");

            Assert.DoesNotThrow(() => st.PrintTree());
            Assert.That(st.Contains("!@#"), Is.True);
            Assert.That(st.Contains("$%^"), Is.True);
            Assert.That(st.Contains("&*()"), Is.True);
            Assert.That(st.Contains("a!@#$%^&*()b"), Is.True);
        }

        [Test]
        public void Fix_PrintTree_DigitsAndMixedCase()
        {
            var st = SuffixTree.Build("Test123ABC");

            Assert.DoesNotThrow(() => st.PrintTree());
            Assert.That(st.Contains("Test"), Is.True);
            Assert.That(st.Contains("123"), Is.True);
            Assert.That(st.Contains("ABC"), Is.True);
            Assert.That(st.Contains("t123A"), Is.True);
        }

        [Test]
        public void Fix_PrintTree_WhitespaceAndNewlines()
        {
            var st = SuffixTree.Build("hello world\ntest");

            Assert.DoesNotThrow(() => st.PrintTree());
            Assert.That(st.Contains("hello world"), Is.True);
            Assert.That(st.Contains(" world\n"), Is.True);
            Assert.That(st.Contains("\ntest"), Is.True);
        }

        // ===========================================
        // FIX 7: ToString() method added
        // ===========================================

        [Test]
        public void Fix_ToString_ContainsOriginalContent()
        {
            var original = "hello world";
            var st = SuffixTree.Build(original);

            // ToString should contain useful information about the tree
            var result = st.ToString();
            Assert.That(result, Does.Contain("hello world"));
            Assert.That(result, Does.Contain("SuffixTree"));
        }

        [Test]
        public void Fix_ToString_EmptyStringShowsEmpty()
        {
            var st = SuffixTree.Build("");
            var result = st.ToString();
            Assert.That(result, Does.Contain("empty"));
        }

        [Test]
        public void Fix_ToString_UnicodeStringShowsContent()
        {
            var original = "日本語テスト";
            var st = SuffixTree.Build(original);
            var result = st.ToString();
            Assert.That(result, Does.Contain("日本語テスト"));
        }

        // ===========================================
        // FIX 8: Random seed in tests for determinism
        // ===========================================

        [Test]
        public void Fix_RandomSeed_Deterministic()
        {
            // Running same random test twice should produce identical results
            var results1 = RunDeterministicRandomTest(12345);
            var results2 = RunDeterministicRandomTest(12345);

            Assert.That(results1, Is.EqualTo(results2),
                "Same seed should produce identical test strings");
        }

        private static string[] RunDeterministicRandomTest(int seed)
        {
            var r = new Random(seed);
            var results = new string[5];
            for (int i = 0; i < 5; i++)
            {
                results[i] = MakeRandomString(r, 20);
            }
            return results;
        }

        // ===========================================
        // FIX 9: Removed dead commented code
        // These tests verify the API is clean and complete
        // ===========================================

        [Test]
        public void Fix_CleanAPI_BuildAndContainsWork()
        {
            // The only public API: Build() and Contains()
            var st = SuffixTree.Build("test");

            Assert.That(st.Contains("test"), Is.True);
            Assert.That(st.Contains("es"), Is.True);
            Assert.That(st.Contains("xyz"), Is.False);
        }

        [Test]
        public void Fix_CleanAPI_PrintTreeDoesNotThrow()
        {
            var st = SuffixTree.Build("abracadabra");

            // PrintTree should work without throwing
            Assert.DoesNotThrow(() => st.PrintTree());
        }

        // ===========================================
        // FIX 10: Edge boundary handling (BOUNDLESS constant)
        // ===========================================

        [Test]
        public void Fix_EdgeBoundary_VeryLongString()
        {
            // Test that BOUNDLESS (-1) handles long strings correctly
            var longString = new string('a', 10000) + "b";
            var st = SuffixTree.Build(longString);

            Assert.That(st.Contains(longString), Is.True);
            Assert.That(st.Contains("aaaaaaaaab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);
        }

        [Test]
        public void Fix_EdgeBoundary_IncrementalConstruction()
        {
            // Verify edges are extended correctly during construction
            var st = SuffixTree.Build("abcdefghijklmnop");

            // Each suffix must be found
            var s = "abcdefghijklmnop";
            for (int i = 0; i < s.Length; i++)
            {
                Assert.That(st.Contains(s.Substring(i)), Is.True,
                    $"Suffix starting at {i} not found");
            }
        }

        #endregion

        #region Additional Algorithm Verification Tests

        [Test]
        public void Algorithm_WalkDown_MultipleEdgeTraversals()
        {
            // Test case that forces multiple walk-down operations
            // When activeLength > edge length, we must walk down
            var st = SuffixTree.Build("abcabcabcabc");

            // All substrings must be found
            var s = "abcabcabcabc";
            for (int i = 0; i < s.Length; i++)
            {
                for (int len = 1; len <= s.Length - i; len++)
                {
                    var substr = s.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True,
                        $"Substring '{substr}' not found");
                }
            }
        }

        [Test]
        public void Algorithm_SuffixLinkChain_DeepNesting()
        {
            // Pattern that creates deep suffix link chains
            var st = SuffixTree.Build("aaabaaabaaab");

            Assert.That(st.Contains("aaabaaab"), Is.True);
            Assert.That(st.Contains("aabaaaba"), Is.True);
            Assert.That(st.Contains("abaaabaa"), Is.True);
            Assert.That(st.Contains("baaabaaab"), Is.True);
        }

        [Test]
        public void Algorithm_EdgeSplit_AtDifferentPositions()
        {
            // Forces splits at various edge positions
            var st = SuffixTree.Build("abcabdabeabf");

            // All these share "ab" prefix but diverge
            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("abd"), Is.True);
            Assert.That(st.Contains("abe"), Is.True);
            Assert.That(st.Contains("abf"), Is.True);
            Assert.That(st.Contains("abcabd"), Is.True);
            Assert.That(st.Contains("abdabe"), Is.True);
            Assert.That(st.Contains("abeabf"), Is.True);
        }

        [Test]
        public void Algorithm_Rule1_ActivePointAtRoot()
        {
            // Test Rule 1: when activeNode is root, decrement activeLength
            var st = SuffixTree.Build("aabbaabb");

            Assert.That(st.Contains("aabb"), Is.True);
            Assert.That(st.Contains("abba"), Is.True);
            Assert.That(st.Contains("bbaa"), Is.True);
            Assert.That(st.Contains("baab"), Is.True);
        }

        [Test]
        public void Algorithm_Rule3_FollowSuffixLink()
        {
            // Test Rule 3: follow suffix link when not at root
            var st = SuffixTree.Build("xabxac");

            // This pattern exercises suffix link following
            Assert.That(st.Contains("xab"), Is.True);
            Assert.That(st.Contains("xac"), Is.True);
            Assert.That(st.Contains("abxa"), Is.True);
            Assert.That(st.Contains("bxac"), Is.True);
        }

        [Test]
        public void Algorithm_ImplicitToExplicit_Terminator()
        {
            // Terminator converts implicit suffixes to explicit
            // In "aa", suffix "a" would be implicit without terminator
            var st = SuffixTree.Build("aa");

            Assert.That(st.Contains("aa"), Is.True);
            Assert.That(st.Contains("a"), Is.True);
        }

        [Test]
        public void Algorithm_ExhaustiveSmallAlphabet()
        {
            // Exhaustive test with binary alphabet
            var r = new Random(555);
            for (int trial = 0; trial < 20; trial++)
            {
                var chars = new char[30];
                for (int i = 0; i < chars.Length; i++)
                    chars[i] = (char)('a' + r.Next(3));

                var s = new string(chars);
                var st = SuffixTree.Build(s);

                // Verify ALL substrings
                for (int i = 0; i < s.Length; i++)
                {
                    for (int len = 1; len <= s.Length - i; len++)
                    {
                        var substr = s.Substring(i, len);
                        Assert.That(st.Contains(substr), Is.True,
                            $"Trial {trial}, string '{s}': Substring '{substr}' not found");
                    }
                }
            }
        }

        [Test]
        public void Algorithm_NegativeTest_NonExistentSubstrings()
        {
            var st = SuffixTree.Build("abcdefgh");

            // These should NOT be found
            Assert.That(st.Contains("xyz"), Is.False);
            Assert.That(st.Contains("abcdefghi"), Is.False); // Too long
            Assert.That(st.Contains("ba"), Is.False);
            Assert.That(st.Contains("hg"), Is.False);
            Assert.That(st.Contains("aba"), Is.False);
        }

        [Test]
        public void Algorithm_AllSuffixesExplicit_MultiplePatterns()
        {
            string[] patterns = {
                "banana",
                "abracadabra",
                "mississippi",
                "abcabxabcd",
                "aaaabaaaabaab"
            };

            foreach (var pattern in patterns)
            {
                var st = SuffixTree.Build(pattern);

                // Every suffix must be found
                for (int i = 0; i < pattern.Length; i++)
                {
                    var suffix = pattern.Substring(i);
                    Assert.That(st.Contains(suffix), Is.True,
                        $"Pattern '{pattern}': Suffix '{suffix}' not found");
                }
            }
        }

        #endregion

        #region Important Additional Tests

        // ===========================================
        // Overlapping patterns - critical for suffix links
        // ===========================================

        [Test]
        public void Contains_OverlappingPatterns_ABAB()
        {
            var st = SuffixTree.Build("ababababab");

            Assert.That(st.Contains("abab"), Is.True);
            Assert.That(st.Contains("baba"), Is.True);
            Assert.That(st.Contains("ababa"), Is.True);
            Assert.That(st.Contains("babab"), Is.True);
            Assert.That(st.Contains("ababababab"), Is.True);
            Assert.That(st.Contains("abababababa"), Is.False); // Too long
        }

        [Test]
        public void Contains_OverlappingPatterns_AAA()
        {
            var st = SuffixTree.Build("aaaaaaaaaa");

            for (int len = 1; len <= 10; len++)
            {
                Assert.That(st.Contains(new string('a', len)), Is.True,
                    $"String of {len} 'a's not found");
            }
            Assert.That(st.Contains(new string('a', 11)), Is.False);
        }

        // ===========================================
        // Prefix/Suffix edge cases
        // ===========================================

        [Test]
        public void Contains_PrefixAndSuffix_SamePattern()
        {
            // String where prefix equals suffix
            var st = SuffixTree.Build("abcabc");

            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("abcabc"), Is.True);
            Assert.That(st.Contains("bca"), Is.True);
            Assert.That(st.Contains("cab"), Is.True);
            Assert.That(st.Contains("abcab"), Is.True);
            Assert.That(st.Contains("bcabc"), Is.True);
        }

        [Test]
        public void Contains_LongCommonPrefix()
        {
            // Multiple suffixes share long common prefix
            var st = SuffixTree.Build("aaaaab");

            Assert.That(st.Contains("aaaaa"), Is.True);
            Assert.That(st.Contains("aaaab"), Is.True);
            Assert.That(st.Contains("aaab"), Is.True);
            Assert.That(st.Contains("aab"), Is.True);
            Assert.That(st.Contains("ab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);
        }

        // ===========================================
        // Single character variations
        // ===========================================

        [Test]
        public void Contains_TwoCharacterAlphabet_Exhaustive()
        {
            var st = SuffixTree.Build("aabbab");

            // All substrings
            string[] expected = { "a", "b", "aa", "ab", "bb", "ba",
                                  "aab", "abb", "bba", "bab",
                                  "aabb", "abba", "bbab",
                                  "aabba", "abbab", "aabbab" };

            foreach (var s in expected)
            {
                Assert.That(st.Contains(s), Is.True, $"'{s}' not found");
            }
        }

        // ===========================================
        // Boundary conditions
        // ===========================================

        [Test]
        public void Contains_FirstAndLastCharacter()
        {
            var st = SuffixTree.Build("xyzabc");

            Assert.That(st.Contains("x"), Is.True);  // First
            Assert.That(st.Contains("c"), Is.True);  // Last
            Assert.That(st.Contains("xy"), Is.True); // First two
            Assert.That(st.Contains("bc"), Is.True); // Last two
        }

        [Test]
        public void Contains_FullString()
        {
            string[] strings = { "a", "ab", "abc", "abcd", "abcdefghij" };

            foreach (var s in strings)
            {
                var st = SuffixTree.Build(s);
                Assert.That(st.Contains(s), Is.True, $"Full string '{s}' not found");
            }
        }

        // ===========================================
        // Repeated substrings
        // ===========================================

        [Test]
        public void Contains_RepeatedSubstring()
        {
            var st = SuffixTree.Build("abcabcabc");

            // "abc" appears 3 times
            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("abcabc"), Is.True);
            Assert.That(st.Contains("abcabcabc"), Is.True);
            Assert.That(st.Contains("bcabca"), Is.True);
            Assert.That(st.Contains("cabcab"), Is.True);
        }

        [Test]
        public void Contains_NoRepeatedCharacters()
        {
            var st = SuffixTree.Build("abcdefghij");

            // Every position should be unique
            for (int i = 0; i < 10; i++)
            {
                char c = (char)('a' + i);
                Assert.That(st.Contains(c.ToString()), Is.True);
            }

            Assert.That(st.Contains("aa"), Is.False);
            Assert.That(st.Contains("ba"), Is.False);
        }

        // ===========================================
        // Walk-down edge cases
        // ===========================================

        [Test]
        public void Algorithm_WalkDown_ExactEdgeLength()
        {
            // Forces walk-down when activeLength exactly equals edge length
            var st = SuffixTree.Build("abcdefabcdef");

            Assert.That(st.Contains("abcdef"), Is.True);
            Assert.That(st.Contains("bcdef"), Is.True);
            Assert.That(st.Contains("cdefab"), Is.True);
            Assert.That(st.Contains("defabc"), Is.True);
        }

        [Test]
        public void Algorithm_WalkDown_MultipleSteps()
        {
            // Long repeating pattern forces multiple walk-down steps
            var st = SuffixTree.Build("abcabcabcabcabc");

            Assert.That(st.Contains("abcabcabcabcabc"), Is.True);
            Assert.That(st.Contains("bcabcabcabcab"), Is.True);
            Assert.That(st.Contains("cabcabcabcabc"), Is.True);
        }

        // ===========================================
        // Special patterns that stress suffix links
        // ===========================================

        [Test]
        public void Algorithm_SuffixLinks_Fibonacci()
        {
            // Fibonacci-like string: each char depends on previous two
            // These strings are known to stress suffix tree construction
            var st = SuffixTree.Build("abaababaabaab");

            Assert.That(st.Contains("abaab"), Is.True);
            Assert.That(st.Contains("baaba"), Is.True);
            Assert.That(st.Contains("aabab"), Is.True);
            Assert.That(st.Contains("abaababaabaab"), Is.True);
        }

        [Test]
        public void Algorithm_SuffixLinks_DeWaarPatterns()
        {
            // de Waard patterns: (ab)^n a - known edge case
            string[] patterns = { "aba", "ababa", "abababa", "ababababa" };

            foreach (var p in patterns)
            {
                var st = SuffixTree.Build(p);

                // Verify all suffixes
                for (int i = 0; i < p.Length; i++)
                {
                    Assert.That(st.Contains(p.Substring(i)), Is.True,
                        $"Pattern '{p}': suffix at {i} not found");
                }
            }
        }

        // ===========================================
        // Performance regression tests
        // ===========================================

        [Test]
        public void Performance_LongString_Build()
        {
            // 10000 characters should build quickly
            var r = new Random(999);
            var chars = new char[10000];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = (char)('a' + r.Next(26));

            var s = new string(chars);

            Assert.DoesNotThrow(() => SuffixTree.Build(s));
        }

        [Test]
        public void Performance_LongString_Contains()
        {
            var r = new Random(888);
            var chars = new char[10000];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = (char)('a' + r.Next(26));

            var s = new string(chars);
            var st = SuffixTree.Build(s);

            // 1000 random searches
            for (int i = 0; i < 1000; i++)
            {
                int start = r.Next(s.Length - 100);
                var substr = s.Substring(start, 50);
                Assert.That(st.Contains(substr), Is.True);
            }
        }

        // ===========================================
        // Unicode and special characters
        // ===========================================

        [Test]
        public void Contains_CyrillicText()
        {
            var st = SuffixTree.Build("привет мир");

            Assert.That(st.Contains("привет"), Is.True);
            Assert.That(st.Contains("мир"), Is.True);
            Assert.That(st.Contains("вет м"), Is.True);
            Assert.That(st.Contains(" "), Is.True);
        }

        [Test]
        public void Contains_ChineseCharacters()
        {
            var st = SuffixTree.Build("你好世界");

            Assert.That(st.Contains("你"), Is.True);
            Assert.That(st.Contains("好"), Is.True);
            Assert.That(st.Contains("世界"), Is.True);
            Assert.That(st.Contains("你好世界"), Is.True);
        }

        [Test]
        public void Contains_Emoji()
        {
            var st = SuffixTree.Build("hello🌍world🌎test");

            Assert.That(st.Contains("hello"), Is.True);
            Assert.That(st.Contains("world"), Is.True);
            Assert.That(st.Contains("🌍"), Is.True);
            Assert.That(st.Contains("🌎"), Is.True);
            Assert.That(st.Contains("🌍world"), Is.True);
        }

        [Test]
        public void Contains_MixedCase()
        {
            var st = SuffixTree.Build("AbCdEfAbCdEf");

            Assert.That(st.Contains("AbCd"), Is.True);
            Assert.That(st.Contains("abcd"), Is.False); // Case sensitive
            Assert.That(st.Contains("ABCD"), Is.False);
            Assert.That(st.Contains("CdEfAb"), Is.True);
        }

        // ===========================================
        // Whitespace handling
        // ===========================================

        [Test]
        public void Contains_Whitespace()
        {
            var st = SuffixTree.Build("a b c a b c");

            Assert.That(st.Contains(" "), Is.True);
            Assert.That(st.Contains("a "), Is.True);
            Assert.That(st.Contains(" b"), Is.True);
            Assert.That(st.Contains("a b c"), Is.True);
            Assert.That(st.Contains(" a "), Is.True);
        }

        [Test]
        public void Contains_TabsAndNewlines()
        {
            var st = SuffixTree.Build("line1\nline2\tdata");

            Assert.That(st.Contains("\n"), Is.True);
            Assert.That(st.Contains("\t"), Is.True);
            Assert.That(st.Contains("line1\n"), Is.True);
            Assert.That(st.Contains("\nline2"), Is.True);
            Assert.That(st.Contains("line2\t"), Is.True);
        }

        // ===========================================
        // Edge cases with numbers
        // ===========================================

        [Test]
        public void Contains_NumericString()
        {
            var st = SuffixTree.Build("123123123");

            Assert.That(st.Contains("123"), Is.True);
            Assert.That(st.Contains("231"), Is.True);
            Assert.That(st.Contains("312"), Is.True);
            Assert.That(st.Contains("123123123"), Is.True);
        }

        [Test]
        public void Contains_AlphanumericMix()
        {
            var st = SuffixTree.Build("a1b2c3a1b2c3");

            Assert.That(st.Contains("a1b2"), Is.True);
            Assert.That(st.Contains("1b2c"), Is.True);
            Assert.That(st.Contains("b2c3a"), Is.True);
            Assert.That(st.Contains("3a1"), Is.True);
        }

        #endregion
    }
}
