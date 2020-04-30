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
static FILE *g_fOut;

bool SearchBitVec1(SeqInfo *Query, const BitVec &BV,
  uint WordLength, uint64 ShiftMask)
	{
	const uint64 BitCount = BV.m_BitCount;
	const uint L = Query->m_L;
	const byte *Seq = Query->m_Seq;
	for (int Strand = 0; Strand <= 1; ++Strand)
		{
		if (Strand == 1)
			Query->RevCompInPlace();

		uint Count = 0;
		uint64 Word = 0;
		byte K = 0;
		for (uint SeqPos = 0; SeqPos < L - WordLength + 1; ++SeqPos)
			{
			byte c = Seq[SeqPos];
			byte Letter = g_CharToLetterNucleo[c];
			if (Letter == INVALID_LETTER)
				{
				K = 0;
				Word = 0;
				continue;
				}
			if (K < WordLength)
				++K;
			Word = ((Word << uint64(2)) | Letter) & ShiftMask;
			if (K == WordLength)
				{
				if (BV.GetBit(Word))
					{
					uint StartPos = SeqPos - WordLength + 1;
					++Count;
					if (Count >= MINCOUNT)
						return true;
					K = 0;
					Word = 0;
					}
				}
			}
		}
	return false;
	}

static void SearchThread(FASTQSeqSource &FSS, const BitVec &BV,
  uint WordLength, uint64 ShiftMask)
	{
	SeqInfo *Query = ObjMgr::GetSeqInfo();
	uint QueryCount = 0;
	uint FoundCount = 0;
	while (FSS.GetNext(Query))
		{
		if (QueryCount%1000 == 0)
			Progress_OnTick();
		++QueryCount;
		bool Found = SearchBitVec1(Query, BV, WordLength, ShiftMask);
		if (Found)
			{
			Lock();
			Query->ToFastq(g_fOut);
			Unlock();
			++FoundCount;
			}
		}
	Lock();
	g_QueryCount += QueryCount;
	g_FoundCount += FoundCount;
	Unlock();
	}

void cmd_search_bitvec()
	{
	const string &QueryFileName = opt(search_bitvec);
	const string &BVFileName = opt(ref);
	g_fOut = CreateStdioFile(opt(output));

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

	FASTQSeqSource FSS;
	FSS.Open(QueryFileName);

	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	time_t t1 = time(0);
	unsigned ThreadCount = GetRequestedThreadCount();
	ProgressFile(FSS.m_LR.m_f, "Searching", QueryFileName);
#pragma omp parallel num_threads(ThreadCount)
	{
	SearchThread(FSS, BV, WordLength, ShiftMask);
	}
	ProgressDone();
	CloseStdioFile(g_fOut);

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
