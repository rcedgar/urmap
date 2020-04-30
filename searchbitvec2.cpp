#include "myutils.h"
#include "bitvec.h"
#include "alpha.h"
#include "seqdb.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "omplock.h"
#include <time.h>
#include "fastqseqsource.h"

static const uint MINCOUNT = 1;
static uint32 g_QueryCount;
static uint32 g_FoundCount;
static FILE *g_fOut1;
static FILE *g_fOut2;

bool SearchBitVec1(SeqInfo *Query, const BitVec &BV,
  uint WordLength, uint64 ShiftMask);

static void SearchThread(FASTQSeqSource &FSS1, FASTQSeqSource &FSS2,
  const BitVec &BV, uint WordLength, uint64 ShiftMask)
	{
	SeqInfo *Query1 = ObjMgr::GetSeqInfo();
	SeqInfo *Query2 = ObjMgr::GetSeqInfo();

	uint QueryCount = 0;
	uint FoundCount = 0;
	for (;;)
		{
		Lock();
		bool Ok1 = FSS1.GetNext(Query1);
		bool Ok2 = FSS2.GetNext(Query2);
		Unlock();

		if (!Ok1)
			{
			asserta(!Ok2);
			break;
			}

		++QueryCount;
		if (QueryCount%1000 == 0)
			Progress_OnTick();

		bool Found1 = SearchBitVec1(Query1, BV, WordLength, ShiftMask);
		bool Found2 = SearchBitVec1(Query2, BV, WordLength, ShiftMask);
		if (Found1 || Found2)
			{
			Lock();
			Query1->ToFastq(g_fOut1);
			Query2->ToFastq(g_fOut2);
			Unlock();
			++FoundCount;
			}
		}
	Lock();
	g_QueryCount += QueryCount;
	g_FoundCount += FoundCount;
	Unlock();
	}

void cmd_search_bitvec2()
	{
	const string &QueryFileName1 = opt(search_bitvec2);
	const string &QueryFileName2 = opt(reverse);
	const string &BVFileName = opt(ref);
	g_fOut1 = CreateStdioFile(opt(output1));
	g_fOut2 = CreateStdioFile(opt(output2));

	FILE *fRef = OpenStdioFile(BVFileName);
	uint32 Magic;
	ReadStdioFile(fRef, &Magic, sizeof(Magic));
	if (Magic != BV_MAGIC)
		Die("Invalid .bv file");

	uint32 WordLength;
	ReadStdioFile(fRef, &WordLength, sizeof(WordLength));
	ProgressLog("Word length %u\n", WordLength);

	uint64 BitCount = myipow64(4, WordLength);
	ProgressLog("%" PRIu64 " bits (%s)\n", BitCount, Int64ToStr(BitCount));

	BitVec BV;
	BV.Alloc(BitCount);

	ProgressOne("Reading bitvec");
	ReadStdioFile64(fRef, BV.m_Vec, BitCount/8);
	CloseStdioFile(fRef);
	ProgressDone();

	FASTQSeqSource FSS1;
	FASTQSeqSource FSS2;
	FSS1.Open(QueryFileName1);
	FSS2.Open(QueryFileName2);

	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	time_t t1 = time(0);
	unsigned ThreadCount = GetRequestedThreadCount();
	ProgressFile(FSS1.m_LR.m_f, "Searching", QueryFileName1);
#pragma omp parallel num_threads(ThreadCount)
	{
	SearchThread(FSS1, FSS2, BV, WordLength, ShiftMask);
	}
	ProgressDone();
	CloseStdioFile(g_fOut1);
	CloseStdioFile(g_fOut2);

	time_t t2 = time(0);
	time_t Secs = t2 - t1;
	if (Secs == 0)
		Secs = 1;
	double ReadsPerSec = double(g_QueryCount)/double(Secs);
	ProgressLog("%u / %u found (%.1f%%)\n", 
	  g_FoundCount, g_QueryCount, GetPct(g_FoundCount, g_QueryCount));
	ProgressLog("Search time %s secs", Int64ToStr((uint64) Secs));
	ProgressLog(", %s reads/sec\n",
	  FloatToStr(ReadsPerSec));
	}
