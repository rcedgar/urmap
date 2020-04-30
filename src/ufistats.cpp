#include "myutils.h"
#include "ufindex.h"
#include <inttypes.h>

void UFIndex::LogStats()
	{
	if (m_SlotCounts == 0)
		CountSlots();
	if (m_SlotCounts_Minus == 0)
		CountSlots_Minus();

	unsigned FreeCount = 0;
	unsigned SinglePlusCount = 0;
	unsigned SingleBothCount = 0;
	unsigned EndCount = 0;
	unsigned MineCount = 0;
	unsigned OtherCount = 0;
	unsigned TruncCount = 0;
	unsigned LongMineCount = 0;
	unsigned LongOtherCount = 0;
	unsigned Trunc2Count = 0;
	unsigned TotalCount = 0;
	uint32 *PosVec = myalloc(uint32, m_MaxIx);
	vector<unsigned> CountHist(256);
	vector<unsigned> TruncHist(256);
	for (uint64 Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		byte T = GetTally(Slot);
		if (T == TALLY_FREE)
			++FreeCount;
		else if (T == TALLY_PLUS1)
			++SinglePlusCount;
		else if (T == TALLY_BOTH1)
			++SingleBothCount;
		else if (T == TALLY_END)
			++EndCount;
		else if (T == TALLY_NEXT_LONG_MINE)
			++LongMineCount;
		else if (T == TALLY_NEXT_LONG_OTHER)
			++LongOtherCount;
		if (T != TALLY_FREE)
			{
			if (TallyMine(T))
				++MineCount;
			else
				++OtherCount;
			}

		unsigned K = GetRow(Slot, PosVec);
		unsigned n = m_SlotCounts[Slot];
		unsigned n_minus = m_SlotCounts_Minus[Slot];
		++(CountHist[n]);
		TotalCount += n;
		if (n > 0 && K < n && n <= m_MaxIx && n_minus <= m_MaxIx)
			{
			++Trunc2Count;
			TruncCount += n;
			++(TruncHist[n]);
			}
		assert(K <= n);
		}

	unsigned Maxi = 0;
	for (unsigned i = 0; i < 256; ++i)
		if (CountHist[i] > 0)
			Maxi = i;
	if (Maxi > 4)
		Maxi = 4;
	for (unsigned i = 0; i <= Maxi; ++i)
		{
		unsigned CH = CountHist[i];
		Log("[%3u]  %10u", i, CH);
		if (i == 0)
			Log("  %7.7s ", "");
		else
			Log("  %7.2f%%", GetPct(CH, TotalCount));
		unsigned TH = TruncHist[i];
		if ((i > 1 && i <= m_MaxIx) || TH != 0)
			{
			Log("  %10u  ", TH);
			Log("  %7.2f%%", GetPct(TH, CH));
			if (i == 1 && TH > 0)
				Log(" <<< TRUNCATED SINGLES");
			if (TH > 0 && i > m_MaxIx)
				Log(" <<< GT MaxIx %u", m_MaxIx);
			}
		Log("\n");
		}

	unsigned IndexedCount, NotIndexedCount, WildcardCount;
	CountIndexedWords(IndexedCount, NotIndexedCount, WildcardCount);

	uint64 CollisionCount;
	uint64 Indexed2Count;
	GetCollisionCount(CollisionCount, Indexed2Count);

	Log("\n");
	Log("%10" PRIu32 "  Word length\n", m_WordLength);
	Log("%10" PRIu32 "  MaxIx\n", m_MaxIx);
	Log("%10" PRIu32 "  Sequence data (%s)\n", m_SeqDataSize, IntToStr(m_SeqDataSize));
	Log("%10" PRIu64 "  Slots (%s)\n", m_SlotCount, Int64ToStr(m_SlotCount));

#define X(x)	Log("%10" PRIu64 "  " #x "\n", x##Count);
	X(Indexed)
	X(Indexed2)
	X(NotIndexed)
	X(Wildcard)
	X(Free)
	X(Collision)
	X(SingleBoth)
	X(SinglePlus)
	X(End)
	X(Mine)
	X(Other)
	X(Trunc)
	X(Trunc2)
	X(LongMine)
	X(LongOther)
	X(Total)
#undef X
	Log("\n");
	}

void cmd_ufi_stats()
	{
	UFIndex UFI;
	UFI.FromFile(opt(ufi_stats));
	UFI.LogStats();
	}

void cmd_ufi_counts()
	{
	UFIndex UFI;
	UFI.FromFile(opt(ufi_counts));
	UFI.CountSlots();
	FILE *f = CreateStdioFile(opt(output));
	WriteStdioFile64(f, UFI.m_SlotCounts, UFI.m_SlotCount);
	CloseStdioFile(f);
	}

void cmd_ufi_validate()
	{
	UFIndex UFI;
	UFI.FromFile(opt(ufi_validate));
	UFI.Validate();
	}

void cmd_ufi_info()
	{
	const string &FileName = opt(ufi_info);
	FILE *f = OpenStdioFile(FileName);
	uint32 WordLength, MaxIx, SeqDataSize;
	uint64 SlotCount;

	uint32 u;
	ReadStdioFile(f, &u, 4);
	asserta(u == UFI_MAGIC1);

	ReadStdioFile(f, &WordLength, 4);
	ReadStdioFile(f, &MaxIx, 4);
	ReadStdioFile(f, &SeqDataSize, 4);
	ReadStdioFile(f, &SlotCount, 8);

	ProgressLog(" Word length  %u\n", WordLength);
	ProgressLog("       MaxIx  %u\n", MaxIx);
	ProgressLog("     SeqData  %u (%s)\n",
	  SeqDataSize, MemBytesToStr((double) SeqDataSize));
	ProgressLog("       Slots  %" PRIu64 " (%s)\n",
	  SlotCount, MemBytesToStr((double) SlotCount));

	CloseStdioFile(f);
	}

