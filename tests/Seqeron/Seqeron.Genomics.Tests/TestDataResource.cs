using System.IO.Compression;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Opens an embedded oracle file by the resource name of its PLAIN form (audit R14-20, B09 F299): <c>Seqeron.Genomics.Tests.TestData.Alignment.x_oracle.tsv</c> opens the resource of that name when it exists and
/// otherwise <c>x_oracle.tsv.gz</c> through a <see cref="GZipStream"/>, so a test names the oracle once and the repository stores it compressed (every oracle above 1 MB is gzip; <c>TestData_Size_Tests</c> keeps it so).
/// </summary>
internal static class TestDataResource
{
    /// <summary>The stream of the resource (decompressed when only the <c>.gz</c> form is embedded), or <c>null</c> when neither exists.</summary>
    public static Stream? Open(string name)
    {
        var assembly = typeof(TestDataResource).Assembly;
        var plain = assembly.GetManifestResourceStream(name);
        if (plain is not null)
            return plain;
        var compressed = assembly.GetManifestResourceStream(name + ".gz");
        return compressed is null ? null : new GZipStream(compressed, CompressionMode.Decompress);
    }
}
