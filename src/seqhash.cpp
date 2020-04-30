#include "myutils.h"
#include "alpha.h"
#include "seqdb.h"
#include "seqhash.h"

uint32 SeqHash32(const byte *Seq, unsigned L)
	{
	unsigned a = 63689;
	unsigned b = 378551;
	uint32 h = 0;

	for (unsigned i = 0; i < L; ++i)
		{
		h = h*a + toupper(Seq[i]);
		a *= b;
		}
	return h;
	}

uint32 SeqHashRC32(const byte *Seq, unsigned L)
	{
	unsigned a = 63689;
	unsigned b = 378551;
	uint32 h = 0;

	for (unsigned k = 0; k < L; ++k)
		{
		unsigned i = L - k - 1;
		h = h*a + toupper(g_CharToCompChar[Seq[i]]);
		a *= b;
		}
	return h;
	}

uint32 SeqHash32_EitherStrand(const byte *Seq, unsigned L)
	{
	uint32 h1 = SeqHash32(Seq, L);
	uint32 h2 = SeqHashRC32(Seq, L);
	uint32 h = h1 ^ h2;
	return h;
	}

bool SeqEq(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2)
	{
	if (L1 != L2)
		return false;
	for (unsigned i = 0; i < L1; ++i)
		if (toupper(Seq1[i]) != toupper(Seq2[i]))
			return false;
	return true;
	}

bool SeqEqRC(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2)
	{
	if (L1 != L2)
		return false;
	for (unsigned i = 0; i < L1; ++i)
		{
		byte c1 = Seq1[i];
		byte c2 = g_CharToCompChar[Seq2[L2 - i - 1]];
		if (toupper(c1) != toupper(c2))
			return false;
		}
	return true;
	}
