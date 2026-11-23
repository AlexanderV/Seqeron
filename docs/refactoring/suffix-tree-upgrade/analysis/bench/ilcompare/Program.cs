using System.Reflection.Metadata;
using System.Reflection.PortableExecutable;
using System.Security.Cryptography;

// ilcompare <old.dll> <new.dll> [type-substring ...]: every method body of both assemblies (the IL bytes, the local variable signature, the exception regions' count) by the full name of its type and its method;
// prints the methods that exist in only one build and the methods whose IL differs, and the number of the methods that are byte for byte the same. With type substrings only the types whose name contains one are listed
// in detail (the count line is always over the whole assembly). Exit code 0 when every method of the named types is identical, 1 otherwise (2 on a usage error).
static class P
{
    static Dictionary<string, string> Bodies(string path)
    {
        using var stream = File.OpenRead(path);
        using var pe = new PEReader(stream);
        MetadataReader md = pe.GetMetadataReader();
        var result = new Dictionary<string, string>();
        foreach (TypeDefinitionHandle th in md.TypeDefinitions)
        {
            TypeDefinition td = md.GetTypeDefinition(th);
            string ns = md.GetString(td.Namespace);
            string name = md.GetString(td.Name);
            string declaring = td.GetDeclaringType().IsNil ? "" : md.GetString(md.GetTypeDefinition(td.GetDeclaringType()).Name) + "+";
            string typeName = (ns.Length > 0 ? ns + "." : "") + declaring + name;
            foreach (MethodDefinitionHandle mh in td.GetMethods())
            {
                MethodDefinition m = md.GetMethodDefinition(mh);
                string key = typeName + "::" + md.GetString(m.Name) + "#" + md.GetBlobBytes(m.Signature).Length;
                string body = "abstract";
                if (m.RelativeVirtualAddress != 0)
                {
                    MethodBodyBlock b = pe.GetMethodBody(m.RelativeVirtualAddress);
                    byte[] il = b.GetILBytes() ?? [];
                    byte[] locals = b.LocalSignature.IsNil ? [] : md.GetBlobBytes(md.GetStandaloneSignature(b.LocalSignature).Signature);
                    body = Convert.ToHexString(SHA256.HashData([.. il, .. locals])) + " il=" + il.Length + " eh=" + b.ExceptionRegions.Length;
                }
                string unique = key;
                for (int i = 2; result.ContainsKey(unique); i++) unique = key + "~" + i;
                result[unique] = body;
            }
        }
        return result;
    }

    static int Main(string[] a)
    {
        if (a.Length < 2) { Console.Error.WriteLine("usage: ilcompare <old.dll> <new.dll> [type-substring ...]"); return 2; }
        Dictionary<string, string> o = Bodies(a[0]), n = Bodies(a[1]);
        string[] focus = a[2..];
        bool Wanted(string k) => focus.Length == 0 || focus.Any(f => k.Contains(f, StringComparison.Ordinal));
        int same = o.Count(kv => n.TryGetValue(kv.Key, out string? v) && v == kv.Value);
        Console.WriteLine($"methods: old {o.Count}, new {n.Count}, byte for byte the same {same}");
        int bad = 0;
        foreach (string k in o.Keys.Union(n.Keys).Order())
        {
            bool inO = o.TryGetValue(k, out string? bo), inN = n.TryGetValue(k, out string? bn);
            string status = !inO ? "ONLY IN NEW" : !inN ? "ONLY IN OLD" : bo == bn ? "" : "DIFFERENT";
            if (status.Length == 0) continue;
            if (Wanted(k)) { Console.WriteLine($"  {status,-12} {k}   old [{bo}] new [{bn}]"); bad++; }
            else Console.WriteLine($"  {status,-12} {k}");
        }
        if (focus.Length > 0) Console.WriteLine($"in the named types ({string.Join(", ", focus)}): {(bad == 0 ? "every method identical" : bad + " method(s) differ")}");
        return bad == 0 ? 0 : 1;
    }
}
