#include "myutils.h"

FILE *g_fSAM;
FILE *g_fTab;
FILE *g_fOut;

void OpenOutputFiles()
	{
	g_fSAM = CreateStdioFile(opt(samout));
	g_fTab = CreateStdioFile(opt(tabbedout));
	g_fOut = CreateStdioFile(opt(output));
	}

void CloseOutputFiles()
	{
	CloseStdioFile(g_fSAM);
	CloseStdioFile(g_fTab);
	CloseStdioFile(g_fOut);
	}
