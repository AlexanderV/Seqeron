// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>An EMBOSS 6.6.0 library function stopped the program with <c>ajFatal</c> (B09 WP67, F306): the message is the text of the library, with the source file and the line of the call.</summary>
public sealed class EmbossProgramFatalException(string message) : Exception(message);
