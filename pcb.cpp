#include "myutils.h"
#include "cmd.h"
#include "seqdb.h"

static char *g_S;
static string g_strQueryFileName;
static string g_strDBFileName;
static const char *g_QueryFileName;
static const char *g_DBFileName;

void SetPCBQueryFileName(const string &FileName)
	{
	g_strQueryFileName = FileName.c_str();
	for (const char *p = g_strQueryFileName.c_str(); *p != 0; ++p)
		{
		if (*p == '/' || *p == '\\')
			g_QueryFileName = p + 1;
		}
	}

void SetPCBDBFileName(const string &FileName)
	{
	g_strDBFileName = FileName.c_str();
	for (const char *p = g_strDBFileName.c_str(); *p != 0; ++p)
		{
		if (*p == '/' || *p == '\\')
			g_DBFileName = p + 1;
		}
	}

static void AllocS()
	{
	if (g_S != 0)
		return;
	g_S = myalloc(char, 128);
	}
