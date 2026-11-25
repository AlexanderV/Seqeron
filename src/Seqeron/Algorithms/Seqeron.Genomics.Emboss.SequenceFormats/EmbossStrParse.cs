// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// <c>ajStrParseWhite</c> / <c>ajStrParseC</c> of <c>ajstr.c</c>: the next token of the string last given (<c>strtok</c> over the white space " \t\n\r"), with the library's static state: a call with no string
/// after the tokens were used up warns "Error in ajStrParseC: NULL argument and not initialised" (B09 WP62, F281).
/// </summary>
public sealed class EmbossStrParse(EmbossFeatureReadContext context)
{
    private static readonly char[] White = [' ', '\t', '\n', '\r'];
    private string _text = "";
    private int _pos;
    private bool _initialised;                       // strParseStr.Ptr != NULL: the last call returned a token

    /// <summary>The next token; <paramref name="text"/> starts a new string, null continues the last one.</summary>
    public string? Next(string? text)
    {
        if (!_initialised && text is null)
        {
            context.Warn("Error in ajStrParseC: NULL argument and not initialised");
            return null;
        }

        if (text is not null)
        {
            _text = text;
            _pos = 0;
        }

        while (_pos < _text.Length && Array.IndexOf(White, _text[_pos]) >= 0)
            _pos++;
        if (_pos >= _text.Length)
        {
            _initialised = false;
            return null;
        }

        int start = _pos;
        while (_pos < _text.Length && Array.IndexOf(White, _text[_pos]) < 0)
            _pos++;
        _initialised = true;
        return _text[start.._pos];
    }
}
