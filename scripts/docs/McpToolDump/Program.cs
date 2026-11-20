using System.ComponentModel;
using System.Reflection;
using System.Runtime.Loader;
using System.Text.Json;
using System.Text.Json.Nodes;

// usage: McpToolDump <server-bin-dir>[:<server-bin-dir>...] [assembly-name ...]   (default assemblies: every Seqeron.Mcp.*.dll of the first directory)
// writes one JSON array to stdout: [{ "server", "name", "description", "method", "parameters": [{ "name", "clr", "json", "items", "hasDefault", "default", "description", "required" }] }]
var binDirs = args[0].Split(':', StringSplitOptions.RemoveEmptyEntries);
AssemblyLoadContext.Default.Resolving += (ctx, name) =>
{
    foreach (var d in binDirs)
    {
        var p = Path.Combine(d, name.Name + ".dll");
        if (File.Exists(p)) return ctx.LoadFromAssemblyPath(p);
    }

    return null;
};

var assemblyNames = args.Length > 1 ? args.Skip(1).ToList() : binDirs.SelectMany(d => Directory.GetFiles(d, "Seqeron.Mcp.*.dll")).Select(Path.GetFileNameWithoutExtension).Distinct().OfType<string>().ToList();
var tools = new JsonArray();
foreach (var an in assemblyNames.OrderBy(x => x, StringComparer.Ordinal))
{
    Assembly asm;
    try { asm = AssemblyLoadContext.Default.LoadFromAssemblyName(new AssemblyName(an)); }
    catch (Exception e) { Console.Error.WriteLine($"cannot load {an}: {e.Message}"); continue; }
    string server = an.Substring("Seqeron.Mcp.".Length);
    foreach (var type in asm.GetExportedTypes().OrderBy(t => t.FullName, StringComparer.Ordinal))
    {
        foreach (var m in type.GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.Instance | BindingFlags.DeclaredOnly))
        {
            var attr = m.GetCustomAttributes().FirstOrDefault(a => a.GetType().Name == "McpServerToolAttribute");
            if (attr is null) continue;
            string? name = attr.GetType().GetProperty("Name")?.GetValue(attr) as string;
            if (string.IsNullOrEmpty(name)) continue;
            var ps = new JsonArray();
            foreach (var p in m.GetParameters())
            {
                var pt = p.ParameterType;
                if (pt == typeof(CancellationToken) || pt.Name is "McpServer" or "IProgress`1" or "RequestContext`1") continue;
                if (pt.Namespace is not null && pt.Namespace.StartsWith("ModelContextProtocol", StringComparison.Ordinal)) continue;
                var (json, items) = JsonType(pt);
                var o = new JsonObject
                {
                    ["name"] = p.Name,
                    ["clr"] = pt.ToString(),
                    ["json"] = json,
                    ["hasDefault"] = p.HasDefaultValue,
                    ["required"] = !p.HasDefaultValue,
                };
                if (items is not null) o["items"] = items;
                if (p.HasDefaultValue && p.DefaultValue is not null && p.DefaultValue != DBNull.Value)
                    o["default"] = DefaultNode(p.DefaultValue);
                o["description"] = p.GetCustomAttribute<DescriptionAttribute>()?.Description;
                ps.Add(o);
            }

            tools.Add(new JsonObject
            {
                ["server"] = server,
                ["name"] = name,
                ["description"] = m.GetCustomAttribute<DescriptionAttribute>()?.Description,
                ["method"] = type.Name + "." + m.Name,
                ["parameters"] = ps,
            });
        }
    }
}

Console.OutputEncoding = new System.Text.UTF8Encoding(false);
Console.Out.Write(tools.ToJsonString(new JsonSerializerOptions { WriteIndented = true, Encoder = System.Text.Encodings.Web.JavaScriptEncoder.UnsafeRelaxedJsonEscaping }));
return 0;

static (string Json, string? Items) JsonType(Type t)
{
    t = Nullable.GetUnderlyingType(t) ?? t;
    if (t == typeof(string) || t == typeof(char) || t.IsEnum) return ("string", null);
    if (t == typeof(bool)) return ("boolean", null);
    if (t == typeof(int) || t == typeof(long) || t == typeof(short) || t == typeof(byte) || t == typeof(sbyte) || t == typeof(uint) || t == typeof(ulong) || t == typeof(ushort)) return ("integer", null);
    if (t == typeof(double) || t == typeof(float) || t == typeof(decimal)) return ("number", null);
    if (t.IsArray) return ("array", JsonType(t.GetElementType()!).Json);
    if (t.IsGenericType && t.GetGenericTypeDefinition() is var g && (g == typeof(List<>) || g == typeof(IEnumerable<>) || g == typeof(IReadOnlyList<>) || g == typeof(IList<>) || g == typeof(IReadOnlyCollection<>) || g == typeof(ICollection<>)))
        return ("array", JsonType(t.GetGenericArguments()[0]).Json);
    return ("object", null);
}

static JsonNode? DefaultNode(object v) => v switch
{
    bool b => JsonValue.Create(b),
    string s => JsonValue.Create(s),
    char c => JsonValue.Create(c.ToString()),
    int i => JsonValue.Create(i),
    long l => JsonValue.Create(l),
    short sh => JsonValue.Create((int)sh),
    byte by => JsonValue.Create((int)by),
    uint ui => JsonValue.Create(ui),
    double d => JsonValue.Create(d),
    float f => JsonValue.Create((double)f),
    decimal m => JsonValue.Create((double)m),
    Enum e => JsonValue.Create(e.ToString()),
    _ => JsonValue.Create(v.ToString()),
};
