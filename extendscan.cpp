#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "ufihit.h"

unsigned State1::AddHSPScan(unsigned StartPosQ, uint32 StartPosDB, bool Plus,
  unsigned Length, int Score)
	{
#if COUNTERS
	bool Correct = DBPosIsCorrect(StartPosDB);
#endif
	asserta(Score > 0);
	unsigned HSPIndex = OverlapsHSP(StartPosQ, StartPosDB, Plus);
	if (HSPIndex != UINT_MAX)
		{
		IncCounter(ScanHSPOverlap);
		UFIHSP *HSP = m_HSPs.Data[HSPIndex];
		if (Score > HSP->m_Score)
			{
			IncCounter(ScanHSPReplace);
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
		return HSPIndex;
		}

	HSPIndex = m_HSPCount++;
	AllocHSPs(m_HSPCount);
	UFIHSP *HSP = m_HSPs.Data[HSPIndex];
	HSP->m_StartPosQ = StartPosQ;
	HSP->m_StartPosDB = StartPosDB;
	HSP->m_Score = Score;
	HSP->m_Length = Length;
	HSP->m_Plus = Plus;
	HSP->m_Aligned = false;
	if (Score > m_BestHSPScore)
		m_BestHSPScore = Score;
	return HSPIndex;
	}

unsigned State1::ExtendScan(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus)
	{
	if (SeedPosDB < SeedPosQ)
		return UINT_MAX;

	uint32 DBLo = SeedPosDB - SeedPosQ;
#if COUNTERS
	bool Correct = DBPosIsCorrect(DBLo);
	if (Correct)
		Inc(ScanCorrectExtend);
#endif
	const byte *QSeq = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	const unsigned QL = m_Query->m_L;
	const byte *DBSeq = m_UFI->m_SeqData + DBLo;
	const int MinHSPScore = int(m_WordLength*2);//TODO:PARAM
	int MaxPenalty = 120;//@@TODO

// Assume seed matches
	int Pen = 0;
	int Score = int(m_WordLength);
	int BestScore = 0;
	int EndPos = int(SeedPosQ+m_WordLength) - 1;
	for (int QPos = EndPos + 1; QPos < int(QL); ++QPos)
		{
		byte q = QSeq[QPos];
		byte t = DBSeq[QPos];
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				EndPos = QPos;
				}
			}
		else
			{
			Pen -= MISMATCH_SCORE;
			if (Pen > m_MaxPenalty)
				{
#if COUNTERS
				if (Correct)
					Inc(ScanCorrectExtendFailedMaxPen);
#endif
				return UINT_MAX;
				}
			Score += MISMATCH_SCORE;
			int Drop = (BestScore - Score);
			if (Drop > XDROP)
				break;
			}
		}

	int StartPos = int(SeedPosQ);
	for (int QPos = StartPos-1; QPos >= 0; --QPos)
		{
		byte q = QSeq[QPos];
		byte t = DBSeq[QPos];
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				StartPos = QPos;
				}
			}
		else
			{
			if (Pen > m_MaxPenalty)
				{
#if COUNTERS
				if (Correct)
					Inc(ScanCorrectExtendFailedMaxPen);
#endif
				return UINT_MAX;
				}
			Score += MISMATCH_SCORE;
			int Drop = (BestScore - Score);
			if (Drop > XDROP)
				break;
			}
		}
	if (StartPos == 0 && EndPos == int(QL) -  1)
		{
#if COUNTERS
		{
		if (Correct)
			Inc(ScanCorrectGlobalHSPs);
		}
#endif
		unsigned HitIndex = AddHitX(DBLo, Plus, BestScore, "");
#if COUNTERS
		{
		if (Correct)
			{
			if (HitIndex == UINT_MAX)
				Inc(ScanCorrectGlobalHSPAddFails);
			else
				{
				Inc(ScanCorrectGlobalHSPSuccess);
				Inc(ScanTotalSuccess);
				}
			}
		}
#endif
		return HitIndex;
		}

	if (BestScore < MinHSPScore)
		return UINT_MAX;

	unsigned HSPIndex = AddHSPScan(unsigned(StartPos), DBLo + unsigned(StartPos),
	  Plus, unsigned(EndPos - StartPos + 1), BestScore);
#if COUNTERS
	if (Correct && HSPIndex == UINT_MAX)
		Inc(ScanCorrectHSPAddHSPFail);
#endif
	unsigned HitIndex = UINT_MAX;
	if (HSPIndex != UINT_MAX)
		{
		HitIndex = AlignHSP(HSPIndex);
#if COUNTERS
		if (Correct)
			{
			if (HitIndex == UINT_MAX)
				Inc(ScanCorrectExtendFailedAlign);
			else
				{
				Inc(ScanTotalSuccess);
				Inc(ScanAlignSuccess);
				}
			}
#endif
		}
	return HitIndex;
	}
