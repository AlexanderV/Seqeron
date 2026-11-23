namespace Seqeron.Genomics.Analysis;

/// <summary>An EMBOSS 6.6.0 library function stopped the program with <c>ajFatal</c> (B09 WP67, F306): the message is the text of the library, with the source file and the line of the call.</summary>
public sealed class EmbossProgramFatalException(string message) : Exception(message);
