# Mutants of the lists of 1b.7d and 1b.7e whose pattern no longer matches the production code at the head exactly once (1b.7e split SaBuilder into the path mode and the section mode and added a second sink to SaSink.cs,
# 1b.9b moved the BUILD_INFO room check into CheckBuildInfoFits): the same mutation written again against the code as it is now. Key: ('d' = list of 1b.7d | 'e' = list of 1b.7e, id).
FIXUPS = {
    ('d', 'Q54'): {'old': 'MappedOutputFile.Create(saPath!, 4 * (n + 1), options.Files)', 'new': 'MappedOutputFile.Create(saPath!, 4 * n, options.Files)'},
    ('d', 'Q56'): {'old': '            saCreated = owned is not null;\n', 'new': ''},
    ('d', 'Q94'): {'old': '                    Trace = options.Trace,\n                }, section is null ?', 'new': '                }, section is null ?'},
    ('d', 'Q96'): {'old': 'new MappedSaSliceSink(saFile, p.WindowBytes, counters)', 'new': 'new MappedSaSliceSink(saFile, p.WindowBytes, null)'},
    ('d', 'Q100'): {'old': '                CheckLength(saPath!, n);\n', 'new': ''},
    ('d', 'Q114'): {'old': 'if (saCreated) DeleteQuietly(saPath!);', 'new': ''},
    ('d', 'Z04'): {'old': 'MemoryMarshal.AsBytes(data));\n\n        public void Dispose() => inner.Dispose();\n    }\n}\n\n/// <summary>The suffix array as a section of a larger file',
                   'new': 'MemoryMarshal.AsBytes(data)[..(data.Length * 2)]);\n\n        public void Dispose() => inner.Dispose();\n    }\n}\n\n/// <summary>The suffix array as a section of a larger file'},
    ('e', 'B116'): {'old': 'BuildInfoText(info, placed.Kind, placed.OverflowCount, options.IncludePsi, childCodec)', 'new': 'BuildInfoText(info, placed.Kind, 0, options.IncludePsi, childCodec)'},
    ('e', 'B118'): {'old': 'if (length > capacity) throw new InvalidOperationException($"BUILD_INFO has', 'new': 'if (length > capacity + 4096) throw new InvalidOperationException($"BUILD_INFO has'},
    # the first form of the mutant is a syntax error (CS1023: the statement of the using is empty)
    ('e', 'G14'): {'old': 'using (var stream = new FileStream(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) stream.SetLength(newLength);', 'new': 'using (var stream = new FileStream(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) _ = stream.Length;'},
}
