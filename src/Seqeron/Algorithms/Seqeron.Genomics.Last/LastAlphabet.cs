// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: Alphabet.cc, ScoreMatrixRow.hh), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Last;

/// <summary>The sequence alphabet of LAST (<c>Alphabet.cc</c>): the letters of the database alphabet get the codes 0 ..., the delimiter (a blank) the next one, then A - Z, '*', a - z and '.'.</summary>
internal sealed class LastAlphabet
{
    /// <summary><c>ALPHABET_CAPACITY</c> of the makefile.</summary>
    public const int Capacity = 66;

    /// <summary>The code of a symbol that is not in the alphabet.</summary>
    public const int DummyCode = Capacity - 1;

    /// <summary><c>Alphabet::dna</c>.</summary>
    public const string Dna = "ACGT";

    public string Letters { get; }

    /// <summary>The number of letters of the alphabet proper (the delimiter has the code <see cref="Size"/>).</summary>
    public int Size { get; private set; }

    public byte[] Encode { get; } = new byte[256];

    public byte[] Decode { get; } = new byte[Capacity];

    public byte[] NumbersToUppercase { get; } = new byte[Capacity];

    public byte[] NumbersToLowercase { get; } = new byte[Capacity];

    public byte[] LettersToUppercase { get; } = new byte[256];

    public byte[] Complement { get; } = new byte[Capacity];

    public LastAlphabet(string letters)
    {
        Letters = letters.ToUpperInvariant();
        Init();
    }

    public bool IsDna => Letters == Dna;

    public bool IsProtein => Letters.Length == 20 && Letters == "ACDEFGHIKLMNPQRSTVWY";

    private void Init()
    {
        Array.Fill(Encode, (byte)DummyCode);
        int code = 0;
        AddLetters(Letters, ref code);
        Size = code;
        AddLetters(" ", ref code);
        if (code - 1 != Letters.Length)
            throw new InvalidOperationException("bad alphabet: " + Letters);
        AddLetters("ABCDEFGHIJKLMNOPQRSTUVWXYZ", ref code);
        AddLetters("*", ref code);
        AddLetters("abcdefghijklmnopqrstuvwxyz", ref code);
        AddLetters(".", ref code);
        for (int i = 0; i < code; i++)
        {
            byte letter = Decode[i];
            NumbersToUppercase[i] = Encode[char.ToUpperInvariant((char)letter) & 0xFF];
            NumbersToLowercase[i] = Encode[char.ToLowerInvariant((char)letter) & 0xFF];
        }

        for (int i = code; i < Capacity; i++)
        {
            NumbersToUppercase[i] = (byte)i;
            NumbersToLowercase[i] = (byte)i;
        }

        if (Letters == Dna)
        {
            // kludge for RNA
            Encode['U'] = Encode['T'];
            Encode['u'] = Encode['t'];
        }

        for (int i = 0; i < 256; i++)
            LettersToUppercase[i] = NumbersToUppercase[Encode[i]];
        MakeComplement();
    }

    private void AddLetters(string letters, ref int code)
    {
        foreach (char c in letters)
        {
            if (Encode[c] == DummyCode)
            {
                Encode[c] = (byte)code;
                Decode[code] = (byte)c;
                code++;
            }
        }
    }

    private void MakeComplement()
    {
        const string x = "ACGTRYSWKMBDHVN";
        const string y = "TGCAYRSWMKVHDBN";
        for (int i = 0; i < Capacity; i++)
            Complement[i] = (byte)i;
        for (int i = 0; i < x.Length; i++)
        {
            Complement[Encode[x[i]]] = Encode[y[i]];
            Complement[Encode[char.ToLowerInvariant(x[i])]] = Encode[char.ToLowerInvariant(y[i])];
        }
    }

    /// <summary><c>Alphabet::tr</c>: the letters of a sequence to codes (a letter that is not in the alphabet stops the program).</summary>
    public void Translate(byte[] seq, int begin, int end, bool keepLowercase)
    {
        byte[] map = keepLowercase ? Encode : LettersToUppercase;
        for (int i = begin; i < end; i++)
        {
            byte code = map[seq[i]];
            if (code == DummyCode)
                throw new LastProgramException("bad symbol in sequence: " + (char)seq[i]);
            seq[i] = code;
        }
    }
}

/// <summary>An error that LAST reports as <c>program: message</c> on the standard error and answers with exit status 1.</summary>
public sealed class LastProgramException(string message) : Exception(message);
