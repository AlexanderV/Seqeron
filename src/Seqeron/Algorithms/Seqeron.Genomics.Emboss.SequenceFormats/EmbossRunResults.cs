// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

/// <summary>The run result of a <c>seqret</c>-style output.</summary>
public static class EmbossRunResults
{
    public static EmbossRunResult From(EmbossSeqretOutput o, string stdoutName = "stdout")
    {
        var files = new Dictionary<string, string>();
        string stdout = "";
        if (o.SingleFiles.Count > 0)
        {
            foreach (var (k, v) in o.SingleFiles)
                files[k] = v;
        }
        else if (o.SequenceFileName == stdoutName || o.SequenceFileName.Length == 0)
        {
            stdout = o.SequenceText;
        }
        else if (o.SequenceText.Length > 0 || !o.Died)
        {
            files[o.SequenceFileName] = o.SequenceText;
        }

        return new EmbossRunResult(stdout, EmbossRunResult.Join(o.Messages), o.Died ? 1 : 0, files);
    }
}
