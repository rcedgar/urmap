#include "myutils.h"
#if 1 // _MSC_VER
#include "zlib.h"
#else
#include <zlib.h>
#endif

FILE *OpenGzipFile(const string &FileName)
	{
	gzFile f = gzopen(FileName.c_str(), "rb");
	if (f == 0)
		Die("Error opening gzip file %s", FileName.c_str());
	return (FILE *) f;
	}

uint32 ReadGzipFile(FILE *f, void *Buff, uint32 MaxBytes)
	{
	int n = gzread(gzFile(f), Buff, MaxBytes);
	if (n < 0)
		Die("Error reading gzip file");
	return unsigned(n);
	}

void CloseGzipFile(FILE *f)
	{
	if (f == 0)
		return;
	gzclose_r(gzFile(f));
	}
