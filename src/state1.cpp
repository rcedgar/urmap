#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "ufihit.h"
#include "alignresult.h"
#include "sort.h"
#include "cigar.h"
#include "omplock.h"

//////////////
//#include "readsimbench.h"
//bool State1::IsCorrect() const
//	{
//	if (m_TopHit == 0)
//		return false;
//
//	string TrueLabel;
//	uint32 TruePos;
//	bool Strand;
//	GetCoordsFromLabel(m_Query->m_Label, TrueLabel, TruePos, Strand);
//
//	string PredLabel;
//	unsigned PredPos;
//	m_TopHit->GetStartCoord(PredLabel, PredPos);
//
//	int64 d = int64(TruePos) - int64(PredPos);
//	if (d < 0)
//		d = -d;
//
//	bool Correct = (TrueLabel == PredLabel && d <= MAX_LOCUS_D);
//	return Correct;
//	}

const UFIHit *State1::GetSecondHit() const
	{
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		const UFIHit *Hit = m_Hits.Data[HitIndex];
		if (Hit != m_TopHit && Hit->m_Score == m_SecondBestScore)
			return Hit;
		}
	return 0;
	}
///////////

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);

int State1::MISMATCH_SCORE;
int State1::GAP_OPEN_SCORE;
int State1::GAP_EXT_SCORE;
int State1::MIN_HSP_SCORE_PCT;
int State1::TERM_HSP_SCORE_PCT_PHASE3;
int State1::XDROP;
int State1::MAX_PENALTY;
int State1::XPHASE1;
int State1::XPHASE3;
int State1::XPHASE4;
unsigned State1::GLOBAL_BAND_RADIUS;

unsigned State1::m_Method = UINT_MAX;
unsigned State1::m_StartSecs;
unsigned State1::m_QueryCount;
unsigned State1::m_AcceptCount;
unsigned State1::m_RejectCount;
unsigned State1::m_NoHitCount;
unsigned State1::m_Minq = 10;

uint32 GetTrueDBPos(SeqInfo *Query, const UFIndex &UFI)
	{
	string TrueLabel;
	uint32 TrueCoord;
	bool TrueStrand;
	GetCoordsFromLabel(Query->m_Label, TrueLabel, TrueCoord, TrueStrand);
	asserta(TrueCoord > 0);
	uint32 TrueDBPos = UFI.CoordToPos(TrueLabel, TrueCoord-1);
	return TrueDBPos;
	}

void State1::AllocCachedVectors(unsigned QL)
	{
	m_SlotsVec_Plus = alloc_buff(uint64, QL);
	m_SlotsVec_Minus = alloc_buff(uint64, QL);
	m_PosVec = alloc_buff(uint32, m_UFI->m_MaxIx);
	m_BlobVec_Plus = alloc_buff(byte, 5 * QL);
	m_BlobVec_Minus = alloc_buff(byte, 5 * QL);
	m_QPosPendingVec_Plus = alloc_buff(byte, QL);
	m_QPosPendingVec_Minus = alloc_buff(byte, QL);
	m_QPosPendingCount_Plus = 0;
	m_QPosPendingCount_Minus = 0;
	}

void State1::InitPE(SeqInfo *Query)
	{
	Lock();
	++m_QueryCount;
	Unlock();

	m_Query = Query;

	const byte *Q = Query->m_Seq;
	unsigned QL = Query->m_L;
	m_AllocBuffOffset = 0;
	AllocCachedVectors(Query->m_L);
	SetSlotsVec(Q, QL, m_SlotsVec_Plus);

	m_RevCompQuerySeq.Alloc(QL);

	RevCompSeq(Q, QL, m_RevCompQuerySeq.Data);
	const byte *QRC = m_RevCompQuerySeq.Data;

	SetSlotsVec(QRC, QL, m_SlotsVec_Minus);

	m_HitCount = 0;
	m_HSPCount = 0;
	m_TopHit = 0;
	m_SecondHit = 0;
	m_BestScore = 0;
	m_BestHSPScore = 0;
	m_SecondBestScore = 0;
	m_Mapq = -1;
	m_MaxPenalty = MAX_PENALTY;
	m_MappedTargetLabel.clear();
	m_MappedTargetPos = UINT32_MAX;
	}

void State1::SetMappedPos()
	{
	m_MappedTargetLabel.clear();
	m_MappedTargetPos = UINT32_MAX;
	if (m_TopHit == 0)
		return;
	unsigned TargetL = 0;
	m_MappedTargetPos = m_UFI->PosToCoordL(m_TopHit->m_DBStartPos,
	  m_MappedTargetLabel, TargetL);
	if (m_MappedTargetPos + m_Query->m_L > TargetL)
		{
		m_MappedTargetLabel.clear();
		m_MappedTargetPos = UINT32_MAX;
		m_TopHit = 0;
		}
	return;
	}

void State1::SetMethod(unsigned Method)
	{
	m_Method = Method;
	if (Method == 1)
		return;
	if (Method == 6 || Method == 8)
		{
		State1::MISMATCH_SCORE = -3;
		State1::GAP_OPEN_SCORE = -5;
		State1::GAP_EXT_SCORE = -1;
		State1::MIN_HSP_SCORE_PCT = 20;
		State1::TERM_HSP_SCORE_PCT_PHASE3 = 60;
		State1::XDROP = 9;
		State1::MAX_PENALTY = 100;
		State1::XPHASE1 = 1;
		State1::XPHASE3 = 1;
		State1::XPHASE4 = 1;
		State1::GLOBAL_BAND_RADIUS = 12;
		}
	else if (Method == 7)
		{
		State1::MISMATCH_SCORE = -4;
		State1::GAP_OPEN_SCORE = -6;
		State1::GAP_EXT_SCORE = -2;
		State1::MIN_HSP_SCORE_PCT = 35;
		State1::TERM_HSP_SCORE_PCT_PHASE3 = 35;
		State1::XDROP = 12;
		State1::MAX_PENALTY = 75;
		State1::XPHASE1 = 8;
		State1::XPHASE3 = 6;
		State1::XPHASE4 = 5;
		State1::GLOBAL_BAND_RADIUS = 8;
		}
	else
		asserta(false);
	return;
	}

void State1::SetUFI(const UFIndex &UFI)
	{
	m_UFI = &UFI;
	m_WordLength = m_UFI->m_WordLength;
	m_ShiftMask = m_UFI->m_ShiftMask;
	m_SlotCount = m_UFI->m_SlotCount;
	}

void State1::AllocHits(unsigned n)
	{
	unsigned OldSize = m_Hits.MaxSize;
	if (OldSize >= n)
		return;
	unsigned NewSize = n + 64;
	UFIHit **OldHits = m_Hits.Data;
	m_Hits.Alloc(NewSize);
	unsigned NewMaxSize = m_Hits.MaxSize;
	memcpy(m_Hits.Data, OldHits, OldSize*sizeof(UFIHit *));
	for (unsigned i = OldSize; i < NewMaxSize; ++i)
		{
		UFIHit *NewHit = new UFIHit;
		asserta(NewHit != 0);
		m_Hits.Data[i] = NewHit;
		NewHit->m_UD = this;
		}
	}

void State1::AllocHSPs(unsigned n)
	{
	unsigned OldSize = m_HSPs.MaxSize;
	if (OldSize >= n)
		return;
	unsigned NewSize = n + 64;
	UFIHSP **OldHSPs = m_HSPs.Data;
	m_HSPs.Alloc(NewSize);
	unsigned NewMaxSize = m_HSPs.MaxSize;
	memcpy(m_HSPs.Data, OldHSPs, OldSize*sizeof(UFIHSP *));
	for (unsigned i = OldSize; i < NewMaxSize; ++i)
		{
		UFIHSP *NewHit = new UFIHSP;
		asserta(NewHit != 0);
		m_HSPs.Data[i] = NewHit;
		}
	}

UFIHit *State1::OverlapsHit(uint32 DBStartPos, bool Plus)
	{
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		UFIHit *Hit = m_Hits.Data[HitIndex];
		if (DBStartPos/64 == Hit->m_DBStartPos/64)//TODO:better
			return Hit;
		}
	return 0;
	}

unsigned State1::OverlapsHSP(uint32 StartPosQ, uint32 StartPosDB, bool Plus)
	{
	int64 Diag = int64(StartPosDB) - int64(StartPosQ);
	for (unsigned HSPIndex = 0; HSPIndex < m_HSPCount; ++HSPIndex)
		{
		UFIHSP *HSP = m_HSPs.Data[HSPIndex];
		int64 HSPDiag = int64(HSP->m_StartPosDB) - int64(HSP->m_StartPosQ);
		if (Diag == HSPDiag)
			return HSPIndex;
		}
	return UINT_MAX;
	}

void State1::LogPos(uint32 DBPos) const
	{
	PrPos(g_fLog, DBPos);
	}

void State1::PrPos(FILE *f, uint32 DBPos) const
	{
	if (DBPos == UINT32_MAX)
		Pf(f, "*");
	else
		{
		string TargetLabel;
		uint32 TargetPos = m_UFI->PosToCoord(DBPos, TargetLabel);
		Pf(f, "%u=%s(%u)", DBPos, TargetLabel.c_str(), TargetPos+1);
		}
	}

void State1::LogHit(unsigned HitIndex) const
	{
	if (HitIndex == UINT_MAX)
		return;
	const UFIHit *Hit = m_Hits.Data[HitIndex];
	string TargetLabel;
	uint32 TargetPos = m_UFI->PosToCoord(Hit->m_DBStartPos, TargetLabel);
	Log("hit %u dbpos %u score %d(%c) %s(%u)",
	  HitIndex, Hit->m_DBStartPos, Hit->m_Score, pom(Hit->m_Plus),
	  TargetLabel.c_str(), TargetPos+1);
	const unsigned QL = m_Query->m_L;
	string Path = Hit->m_Path;
	if (Path.empty())
		for (unsigned i = 0; i < m_Query->m_L; ++i)
			Path.push_back('M');

	const byte *DBSeq = m_UFI->m_SeqData + Hit->m_DBStartPos;
	bool Plus = Hit->m_Plus;

	const byte *QSeq = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	LogAlnPretty(QSeq, DBSeq, Path.c_str(), false);
	}

void State1::LogHSPSummary() const
	{
	Log("Q>%s  %u HSPs\n", m_Query->m_Label, m_HSPCount);
	if (m_HSPCount == 0)
		return;
	vector<int> Scores;
	vector<unsigned> Order;
	for (unsigned i = 0; i < m_HSPCount; ++i)
		{
		int Score = GetHSP(i)->m_Score;
		Scores.push_back(Score);
		}
	Order.resize(m_HSPCount);
	QuickSortOrderDesc(Scores.data(), m_HSPCount, Order.data());
	for (unsigned k = 0; k < m_HSPCount; ++k)
		{
		unsigned HSPIndex = Order[k];
		const UFIHSP *HSP = GetHSP(HSPIndex);
		string TargetLabel;
		unsigned TargetLo;
		Log("HSP[%2u] ", HSPIndex+1);
		Log("  %6.1f", double(HSP->m_Score));
		if (HSP->m_Aligned)
			Log("  ALN");
		else
			Log("  una");
		TargetLo = m_UFI->PosToCoord(HSP->m_StartPosDB, TargetLabel);
		Log("  %5.5s(%9u)%c", TargetLabel.c_str(), TargetLo+1, pom(HSP->m_Plus));
		Log("\n");
		}
	}

void State1::LogHitsSummary() const
	{
	Log("Q>%s  %u hits\n", m_Query->m_Label, m_HitCount);
	if (m_HitCount == 0)
		return;
	vector<unsigned> Order(m_HitCount);
	GetHitOrder(Order);
	for (unsigned k = 0; k < m_HitCount; ++k)
		{
		unsigned HitIndex = Order[k];
		const UFIHit *Hit = m_Hits.Data[HitIndex];
//		string TargetLabel;
//		unsigned TargetLo;
//		Hit->GetStartCoord(TargetLabel, TargetLo);
		Log("Hit[%2u] ", HitIndex+1);
		Log("  %6.1f", double(Hit->m_Score));
		Log("  %10u", Hit->m_DBStartPos);
//		Log("  %5.5s(%9u)%c", TargetLabel.c_str(), TargetLo+1, pom(Hit->m_Plus));
		if (k == 0)
			Log("  MAPQ %.1f", double(m_Mapq));
		Log("\n");
		}
	}

void State1::LogHits() const
	{
	Log("Q>%s  %u hits\n", m_Query->m_Label, m_HitCount);
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		LogHit(HitIndex);
		Log("\n");
		}
	}

int State1::GetSecondBestDelta() const
	{
	int BestScore = GetBestScore();
	int BestDelta = 99;
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		int s = m_Hits.Data[i]->m_Score;
		if (s == BestScore)
			continue;
		int Delta = BestScore - s;
		if (Delta < BestDelta)
			BestDelta = Delta;
		}
	return BestDelta;
	}

const UFIHit *State1::GetTopHit() const
	{
	if (m_HitCount == 0)
		return 0;

	int BestScore = GetBestScore();
	const UFIHit *TopHit = 0;
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		const UFIHit *Hit = m_Hits.Data[i];
		if (Hit->m_Score == BestScore)
			{
			if (TopHit != 0)
				return 0;
			TopHit = Hit;
			}
		}
	return TopHit;
	}

void State1::SetSlotsVec(const byte *Seq, unsigned L, uint64 *SlotsVec) const
	{
	uint64 Word = 0;
	byte K = 0;
	unsigned SlotCount = 0;
	for (uint32 SeqPos = 0; SeqPos < m_WordLength-1; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		}

	for (uint32 SeqPos = m_WordLength-1; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			SlotsVec[SeqPos-m_WordLength+1] = UINT64_MAX;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			SlotsVec[SeqPos-m_WordLength+1] = Slot;
			}
		else
			SlotsVec[SeqPos-m_WordLength+1] = UINT64_MAX;
		}
	}


void State1::SetSlotsVecRC(const byte *Seq, unsigned L, uint64 *SlotsVec) const
	{
	uint64 Word = 0;
	byte K = 0;
	unsigned SlotCount = 0;
	for (uint32 SeqPos = 0; SeqPos < m_WordLength-1; ++SeqPos)
		{
		byte c = Seq[SeqPos-L-1];
		byte Letter = g_CharToCompLetter[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		}

	for (uint32 SeqPos = m_WordLength-1; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos-L-1];
		byte Letter = g_CharToCompLetter[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			SlotsVec[SeqPos-m_WordLength+1] = Slot;
			}
		else
			SlotsVec[SeqPos-m_WordLength+1] = UINT64_MAX;
		}
	}

void State1::SetTallyVec(const uint64 *Slots, unsigned N, byte *TallyVec) const
	{
	for (unsigned i = 0; i < N; ++i)
		{
		uint64 Slot = Slots[i];
		if (Slot == UINT64_MAX)
			TallyVec[i] = 0;
		else
			TallyVec[i] = m_UFI->GetTally(Slot);
		}
	}

void State1::SetBlobVec(const uint64 *Slots, unsigned N, byte *BlobVec) const
	{
	for (unsigned i = 0; i < N; ++i)
		{
		uint64 Slot = Slots[i];
		if (Slot == UINT64_MAX)
			BlobVec[5*i] = TALLY_FREE;
		else
			m_UFI->GetBlob(Slot, BlobVec + 5*i);
		}
	}

unsigned State1::AddHitX(uint32 StartPosDB, bool Plus,
  int Score, const string &Path)
	{
	if (Score < 10)//TODO:PARAM
		return UINT_MAX;
	asserta(Score > 0);
	UFIHit *Hit = OverlapsHit(StartPosDB, Plus);
	if (Hit != 0)
		return UINT_MAX;

	unsigned QL = m_Query->m_L;
	int Pen = int(QL) - int(Score);
	int MaxPen = Pen - 2*MISMATCH_SCORE;
	if (MaxPen < m_MaxPenalty)
		m_MaxPenalty = MaxPen;

	unsigned HitIndex = m_HitCount;
	AllocHits(m_HitCount+1);
	Hit = m_Hits.Data[HitIndex];
	Hit->m_Score = Score;
	Hit->m_Plus = Plus;
	Hit->m_DBStartPos = StartPosDB;
	Hit->m_Path = Path;

	if (Score > m_BestScore)
		{
		m_SecondBestScore = m_BestScore;
		m_BestScore = Score;
		m_TopHit = Hit;
		}
	else if (Score == m_BestScore)
		m_SecondBestScore = Score;
	else
		{
		asserta(Score < m_BestScore);
		if (Score < m_BestScore - SECONDARY_HIT_MAX_DELTA)
			return UINT_MAX;
		if (Score > m_SecondBestScore)
			m_SecondBestScore = Score;
		}

	++m_HitCount;
	return HitIndex;
	}

void State1::AddHSPX(unsigned StartPosQ, uint32 StartPosDB, bool Plus,
  unsigned Length, int Score)
	{
	asserta(Score > 0);

	if (Score < m_BestScore - 4) //TODO:PARAM
		return;

	unsigned HSPIndex = OverlapsHSP(StartPosQ, StartPosDB, Plus);
	if (HSPIndex != UINT_MAX)
		{
		UFIHSP *HSP = m_HSPs.Data[HSPIndex];
		if (Score > HSP->m_Score)
			{
			HSP->m_StartPosQ = StartPosQ;
			HSP->m_StartPosDB = StartPosDB;
			HSP->m_Score = Score;
			HSP->m_Length = Length;
			HSP->m_Plus = Plus;
			HSP->m_Aligned = false;
			}
		else
			{
			;
			}
		return;
		}

	AllocHSPs(m_HSPCount+1);
	UFIHSP *HSP = m_HSPs.Data[m_HSPCount++];
	HSP->m_StartPosQ = StartPosQ;
	HSP->m_StartPosDB = StartPosDB;
	HSP->m_Score = Score;
	HSP->m_Length = Length;
	HSP->m_Plus = Plus;
	HSP->m_Aligned = false;
	if (Score > m_BestHSPScore)
		m_BestHSPScore = Score;
	}

void State1::HitStats()
	{
	unsigned Secs = GetElapsedSecs() - m_StartSecs;
	double ReadsPerSec = 0.0;
	if (Secs > 0)
		ReadsPerSec = double(m_QueryCount)/Secs;

	unsigned ThreadCount = GetRequestedThreadCount();
	//Log("%u threads, %u secs\n", ThreadCount, Secs);

	double AcceptPct = GetPct(m_AcceptCount, m_QueryCount);
	double RejectPct = GetPct(m_RejectCount, m_QueryCount);
	double NoHitPct = GetPct(m_NoHitCount, m_QueryCount);

	Log("@rps=%.1f\n", ReadsPerSec);

	//Log("%u secs, %s rps", Secs, FloatToStr(ReadsPerSec));
	//Log(", acc %u (%.1f%%)", m_AcceptCount, AcceptPct);
	//Log(", rej %u (%.1f%%)", m_RejectCount, RejectPct);
	//Log(", nohit %u (%.1f%%)\n", m_NoHitCount, NoHitPct);
	//Log("\n");

	ProgressPrefix(false);
	ProgressLog("\n");
	ProgressLog("%16u  Seconds to load index\n", m_StartSecs);
	if (Secs < 180)
		ProgressLog("%16u  Seconds in mapper\n", Secs);
	else if (Secs < 2*60*60)
		ProgressLog("%16.1f  Minutes in mapper\n", Secs/60.0);
	else
		ProgressLog("%16.1f  Hours in mapper\n", Secs/(60.0*60.0));
	ProgressLog("%16s  Reads", IntToStrCommas(m_QueryCount));
	ProgressLog(" (%s)\n", IntToStr(m_QueryCount));
	ProgressLog("%16.16s  Reads/sec. (%u threads)\n", FloatToStr(ReadsPerSec), ThreadCount);
	ProgressLog("%16s  Mapped Q>=%u (%.1f%%)\n", IntToStrCommas(m_AcceptCount), m_Minq, AcceptPct);
	ProgressLog("%16s  Mapped Q< %u (%.1f%%)\n", IntToStrCommas(m_RejectCount), m_Minq, RejectPct);
	ProgressLog("%16s  Unmapped (%.1f%%)\n", IntToStrCommas(m_NoHitCount), NoHitPct);
	ProgressLog("\n");
	ProgressPrefix(true);
	}

const UFIHSP *State1::GetHSP(unsigned Index) const
	{
	asserta(Index < m_HSPCount);
	const UFIHSP *HSP = m_HSPs.Data[Index];
	return HSP;
	}

UFIHSP *State1::GetHSP(unsigned Index)
	{
	asserta(Index < m_HSPCount);
	UFIHSP *HSP = m_HSPs.Data[Index];
	return HSP;
	}

const UFIHit *State1::GetHit(unsigned Index) const
	{
	asserta(Index < m_HitCount);
	const UFIHit *Hit = m_Hits.Data[Index];
	return Hit;
	}

void State1::GetTopHitIndexes(vector<unsigned> &HitIndexes) const
	{
	HitIndexes.clear();
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		if (m_Hits.Data[HitIndex]->m_Score == m_BestScore)
			HitIndexes.push_back(HitIndex);
		}
	}

UFIHit *State1::GetHit(unsigned Index)
	{
	asserta(Index < m_HitCount);
	UFIHit *Hit = m_Hits.Data[Index];
	return Hit;
	}

void State1::GetHitOrder(vector<unsigned> &Order) const
	{
	vector<int> Scores;
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		int Score = GetHit(i)->m_Score;
		Scores.push_back(Score);
		}
	Order.resize(m_HitCount);
	QuickSortOrderDesc(Scores.data(), m_HitCount, Order.data());
	}


void State1::LogState(const string &s) const
	{
	Log("LogState(%s)\n", s.c_str());
#define X(x)	Log("  %6.1f  %s\n", double(x), #x)
	X(MISMATCH_SCORE);
	X(GAP_OPEN_SCORE);
	X(GAP_EXT_SCORE);
	X(MIN_HSP_SCORE_PCT);
	X(TERM_HSP_SCORE_PCT_PHASE3);
	X(XDROP);
	X(MAX_PENALTY);
	X(XPHASE1);
	X(XPHASE3);
	X(XPHASE4);
	X(GLOBAL_BAND_RADIUS);
	X(m_MaxPenalty);
	X(m_BestScore);
	X(m_SecondBestScore);
	X(m_HSPCount);
	X(m_HitCount);
	}

bool State1::PolishCIGAR(string &CIGAR) const
	{
	vector<char> Ops;
	vector<unsigned> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	CIGAROpsFixDanglingMs(Ops, Lengths);
	OpsToCIGAR(Ops, Lengths, CIGAR);
	return true;
	}

void State1::GetCIGAR(string &CIGAR) const
	{
	if (m_TopHit == 0)
		{
		CIGAR = "*";
		return;
		}

	const char *Path = m_TopHit->m_Path.c_str();
	if (*Path == 0)
		{
		unsigned QL = m_Query->m_L;
		Ps(CIGAR, "%uM", QL);
		}
	else
		PathToCIGAR(Path, CIGAR);
	PolishCIGAR(CIGAR);
	}

void State1::WriteSAMHeader(FILE *f, const UFIndex &UFI)
	{
	if (f == 0)
		return;

	unsigned SeqCount = SIZE(UFI.m_SeqLengths);
	asserta(SIZE(UFI.m_Labels) == SeqCount);
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		const char *Label = UFI.m_Labels[i].c_str();
		unsigned Len = UFI.m_SeqLengths[i];
		fprintf(f, "@SQ	SN:%s	LN:%u\n", Label, Len);
		}
	fprintf(f, "@PG	ID:urmap	PN:urmap	VN:%s.%s	CL:",
	  MY_VERSION, SVN_VERSION);
	PrintCmdLine(f);
	}

uint32 State1::GetTopHitPos() const
	{
	if (m_TopHit == 0)
		return UINT32_MAX;
	uint32 Pos = m_TopHit->m_DBStartPos;
	return Pos;
	}

void State1::GetPairLabel(string &Label) const
	{
	Label.clear();
	const char *s = m_Query->m_Label;;
	unsigned n = ustrlen(s);
	if (n > 2 && s[n-2] == '/' &&
	  (s[n-1] == '1' || s[n-1] == '2'))
		n -= 2;

	for (unsigned i = 0; i < n; ++i)
		{
		char c = s[i];
		if (isspace(c))
			break;
		Label += c;
		}
	}

void State1::GetChrPos(const UFIHit *Hit, string &ChrLabel, uint32 &ChrPos) const
	{
	if (Hit == 0)
		{
		ChrLabel = "*";
		ChrPos = 0;
		return;
		}

	ChrPos = m_UFI->PosToCoord(Hit->m_DBStartPos, ChrLabel);
	}
