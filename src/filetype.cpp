#include "myutils.h"
#include "filetype.h"

bool FastaFileIsNucleo(FILE *f);

FILE_TYPE GetFileType(const string &FileName, bool *ptrNucleo)
	{
	FILE_TYPE Type = FT_Unknown;
	*ptrNucleo = false;

	FILE *f = OpenStdioFile(FileName);
	uintB FileSize = GetStdioFileSizeB(f);
	if (FileSize == 0)
		Die("Empty file %s", FileName.c_str());

	byte b;
	ReadStdioFile(f, &b, 1);

	if (b == '>')
		{
		Type = FT_FASTA;
		*ptrNucleo = FastaFileIsNucleo(f);
		}
	else if (b == '@')
		{
		Type = FT_FASTQ;
		*ptrNucleo = true;
		}
	CloseStdioFile(f);

	if (Type == FT_Unknown)
		Die("Unknown file format %s", FileName.c_str());

	return Type;
	}
