#include "myutils.h"
#include <time.h>

enum PSTATE
	{
	WAITING = 1,
	PARTLINE = 2,
	BAR = 3,
	BUSY = 4,
	};
static PSTATE g_State = WAITING;

static string g_Msg = "Processing";
static FILE *g_fOut = stderr;
static FILE *g_fMon;
static uint64 g_FileSize = UINT64_MAX;
static FN_PROGRESS_CALLBACK g_CB = 0;
static uint32 g_StepCount32;
static uint64 g_StepCount64;
static const uint32 *g_ptrStepIndex32;
static const uint64 *g_ptrStepIndex64;
static time_t g_TimeLastOutput;
static time_t g_TimeStart;
static time_t g_TimeDone;
static double g_MemStart;
static double g_MemDone;
static bool g_Prefix = true;

static string g_ProgressBar;
static string g_CurrentLine;
static string g_BusyStr;

static void ForceNewLine()
	{
	switch (g_State)
		{
	case WAITING:
		return;

	case PARTLINE:
		fputs(g_CurrentLine.c_str(), g_fOut);
		g_CurrentLine.clear();
		fputc('\n', g_fOut);
		time(&g_TimeLastOutput);
		g_State = WAITING;
		return;

	case BAR:
		fputs(g_CurrentLine.c_str(), g_fOut);
		g_CurrentLine.clear();
		fputc('\n', g_fOut);
		time(&g_TimeLastOutput);
		g_State = WAITING;
		return;

	case BUSY:
		fputc('\n', g_fOut);
		g_State = WAITING;
		return;

	default:
		Warning("ForceNewLine state=%d", g_State);
		}
	}

static void Append(char c)
	{
	if (c == '\n')
		{
		if (g_Prefix)
			{
			fputs(GetResourceStr(), g_fOut);
			fputc(' ', g_fOut);
			}
		fputs(g_CurrentLine.c_str(), g_fOut);
		g_CurrentLine.clear();
		fputc('\n', g_fOut);
		time(&g_TimeLastOutput);
		g_State = WAITING;
		}
	else
		{
		g_CurrentLine += c;
		g_State = PARTLINE;
		}
	}

static void Append(const string &s)
	{
	if (s.empty())
		return;
	for (size_t i = 0; i < s.size(); ++i)
		Append(s[i]);
	}

static void UpdateProgressBar(bool Done)
	{
	if (g_CB == 0)
		return;

	if (g_State != BAR)
		ForceNewLine();

	size_t OldL = g_ProgressBar.size();
	string Msg;
	double PctDone;
	if (g_CB == 0)
		return;
	(*g_CB)(Msg, PctDone);

	g_ProgressBar = string(GetResourceStr());
	if (PctDone >= 0.0)
		{
		if (Done)
			g_ProgressBar += " 100%";
		else
			{
			if (PctDone < 1.0)
				Psa(g_ProgressBar, " %.2f%%", PctDone);
			else if (PctDone > 99.9)
				g_ProgressBar += " 99.9%";
			else
				Psa(g_ProgressBar, " %.1f%%", PctDone);
			}
		}
	g_ProgressBar += " ";
	g_ProgressBar += Msg;

	size_t NewL = g_ProgressBar.size();
	fputc('\r', g_fOut);
	fputs(g_ProgressBar.c_str(), g_fOut);
	for (size_t i = NewL; i < OldL; ++i)
		fputc(' ', g_fOut);
	time(&g_TimeLastOutput);
	if (Done)
		{
		Log("%s\n", g_ProgressBar.c_str());
		g_ProgressBar.clear();
		fputc('\n', g_fOut);
		g_State = WAITING;
		}
	else
		g_State = BAR;
	}

static void CB_Step32(string &s, double &PctDone)
	{
	uint32 Index = *g_ptrStepIndex32;
	if (Index >= g_StepCount32+1)
		Index = g_StepCount32;
	PctDone = GetPct(Index, g_StepCount32+1);
	s = g_Msg;
	}

static void CB_Step64(string &s, double &PctDone)
	{
	uint64 Index = *g_ptrStepIndex64;
	if (Index >= g_StepCount64+1)
		Index = g_StepCount64;
	PctDone = GetPct64(Index, g_StepCount64+1);
	s = g_Msg;
	}

static void CB_File(string &s, double &PctDone)
	{
	if (g_fMon == 0)
		return;
	s = g_Msg;

	uint64 Bytes = GetStdioFileBytesRead(g_fMon);
	if (g_FileSize == UINT64_MAX)
		{
		PctDone = -1.0;
		string t = string(MemBytesToStr(Bytes));
		s += string(" (") + t + string(")");
		}
	else
		PctDone = GetPct((double) Bytes, (double) g_FileSize);
	}

static void CB_One(string &s, double &PctDone)
	{
	s = g_Msg;
	PctDone = -1;
	}

void Progress_OnTick()
	{
	if (g_TimeLastOutput == 0)
		return;

	time_t Now = time(0); 
	time_t Delay = Now - g_TimeLastOutput;
	if (Delay == 0)
		return;

	switch (g_State)
		{
	case WAITING:
	case PARTLINE:
		{
		if (Delay < 3)
			return;

		ForceNewLine();
		g_BusyStr = string(GetResourceStr()) + " Busy";
		fputs(g_BusyStr.c_str(), g_fOut);
		g_State = BUSY;
		break;
		}

	case BAR:
		{
		UpdateProgressBar(false);
		break;
		}

	case BUSY:
		{
		if (Delay < 2)
			return;

		unsigned OldLength = SIZE(g_BusyStr);
		g_BusyStr = string(GetResourceStr()) + " Busy";
		unsigned NewLength = SIZE(g_BusyStr);
		fputc('\r', g_fOut);
		fputs(g_BusyStr.c_str(), g_fOut);
		for (unsigned i = OldLength; i < NewLength; ++i)
			fputc(' ', g_fOut);
		break;
		}

	default:
		Warning("Progress_OnTimer state=%d", g_State);
		break;
		}
	}

void ProgressClose()
	{
	if (g_State != WAITING)
		{
		fputc('\n', g_fOut);
		g_State = WAITING;
		}
	}

void ProgressStart(FN_PROGRESS_CALLBACK PCB)
	{
	PSTATE CurrState = g_State;
	switch (CurrState)
		{
	case PARTLINE:
		Warning("ProgressStart() state=PARTLINE");
		ForceNewLine();
		break;

	case WAITING:
		break;

	case BAR:
		Warning("ProgressStart() state=BAR");
		ForceNewLine();
		break;

	case BUSY:
		ForceNewLine();
		break;

	default:
		Warning("ProgressStart() state=%d", CurrState);
		break;
		}

	g_CB = PCB;
	time(&g_TimeStart);
	g_TimeDone = 0;

	g_MemStart = GetMemUseBytes();
	g_MemDone = 0;

	UpdateProgressBar(false);
	}

void Progress(const char *Format, ...)
	{
	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

	Append(Str);
	}

void ProgressOne(const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(g_Msg, Format, ArgList);
	va_end(ArgList);
	g_CB = CB_One;
	UpdateProgressBar(false);
	}

void ProgressFile(FILE *f, const string &Msg, const string &PathName)
	{
	if (f == 0)
		return;
	if (PathName.empty())
		return;
	g_fMon = f;
	int fd = fileno(f);
	string FileName;
	BaseName(PathName, FileName);
	g_FileSize = GetStdioFileSize64_NoFail(f);
	g_Msg = Msg + string(" ") + FileName;
	ProgressStart(CB_File);
	}

void ProgressFileWrite(FILE *f, const string &Msg, const string &PathName)
	{
	if (f == 0)
		return;
	if (PathName.empty())
		return;
	g_fMon = f;
	int fd = fileno(f);
	string FileName;
	BaseName(PathName, FileName);
	g_FileSize = UINT64_MAX;
	g_Msg = Msg + string(" ") + FileName;
	ProgressStart(CB_File);
	}

void ProgressLoop(uint32 Count, uint32 *ptrIndex, const string &Msg)
	{
	g_StepCount32 = Count;
	g_ptrStepIndex32 = ptrIndex;
	*ptrIndex = 0;
	g_Msg = Msg;
	ProgressStart(CB_Step32);
	}

void ProgressLoop64(uint64 Count, uint64 *ptrIndex, const string &Msg)
	{
	g_StepCount64 = Count;
	g_ptrStepIndex64 = ptrIndex;
	*ptrIndex = 0;
	g_Msg = Msg;
	ProgressStart(CB_Step64);
	}

void ProgressDone()
	{
	time(&g_TimeDone);
	UpdateProgressBar(true);
	g_Msg = "Processing";
	g_CB = 0;
	}

void ProgressPrefix(bool On)
	{
	g_Prefix = On;
	}

void ProgressSetMsg(const string &Msg)
	{
	g_Msg = Msg;
	}

void ProgressLog(const char *Format, ...)
	{
	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

	Log("%s", Str.c_str());
	Progress("%s", Str.c_str());
	}

void Progress_Init()
	{
	}

void Progress_OnExit()
	{
	ProgressPrefix(true);
	Progress("Exiting\n");
	fputc('\n', g_fOut);
	ProgressPrefix(false);
	ProgressClose();
	}
