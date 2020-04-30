#include "myutils.h"
#include "tenx.h"
#include "alpha.h"
#include <set>

uint32 Tenx_GetIntegerBarcode(const byte *Seq)
	{
	uint32 BC = 0;
	for (unsigned i = 0; i < 16; ++i)
		{
		byte Letter = g_CharToLetterNucleo[Seq[i]];
		BC = BC*4 + Letter;
		}
	return BC;
	}

void Tenx_ReadWhitelist(const string &FileName, set<uint32> &WhiteSet)
	{
	WhiteSet.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	ProgressFile(f, "Whitelist", FileName);
	unsigned Counter = 0;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressLoopTick(Counter++);
		if (Line.empty())
			continue;
		if (SIZE(Line) != 16)
			Die("Bad barcode '%s'", Line.c_str());
		uint32 BC = Tenx_GetIntegerBarcode((const byte *) Line.c_str());
		WhiteSet.insert(BC);
		}
	ProgressDone();
	ProgressLog("%s whitelisted barcodes\n", IntToStrCommas(SIZE(WhiteSet)));
	}

// 012345678901234
// B.4f7c42de.1256
uint32 Tenx_GetIntegerBarcodeFromLabel(const char *Label)
	{
	unsigned n = ustrlen(Label);
	if (n < 12)
		Die("Tenx_GetIntegerBarcodeFromLabel(%s), <12 chars", Label);
	if (Label[0] != 'B' || Label[1] != '.' || Label[10] != '.')
		Die("Tenx_GetIntegerBarcodeFromLabel(%s), must start with 'B.'", Label);

	char s[9];
	for (unsigned i = 0; i < 8; ++i)
		s[i] = Label[i+2];
	s[8] = 0;

	char *p;
	unsigned long ul = strtoul(s, &p, 16);
	asserta(ul <= UINT32_MAX);
	unsigned BC = uint32(ul);
	return BC;
	}

const char *Tenx_GetBarcodeFromSeq(const byte *Seq, string &BC)
	{
	BC.clear();
	for (unsigned i = 0; i < 16; ++i)
		BC += char(Seq[i]);
	return BC.c_str();
	}

const char *Tenx_IntegerBarcodeToStr(uint32 BC, string &s)
	{
	s.clear();
	const byte *B = (const byte *) &BC;
	for (unsigned i = 0; i < 16; ++i)
		{
		unsigned i1 = i/4;
		unsigned i2 = i%4;
		byte Letter = (B[i1] >> i2*2) & 0x3;
		char c = g_LetterToCharNucleo[Letter];
		s = c + s;
		}
	return s.c_str();
	}
