#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "ufihit.h"

// Return +score if new hit, -2 if new hsp, -1 otherwise.
int State1::ExtendPen(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus)
	{
	if (SeedPosDB < SeedPosQ)
		return -1;

	uint32 DBLo = SeedPosDB - SeedPosQ;
	UFIHit *Hit = OverlapsHit(DBLo, Plus);
	if (Hit != 0)
		return -1;

	const byte *QSeq = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	const unsigned QL = m_Query->m_L;
	const byte *DBSeq = m_UFI->m_SeqData + DBLo;
	const int MinHSPScore = int(MIN_HSP_SCORE_PCT*QL/100.0);

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
				return -1;
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
			Pen -= MISMATCH_SCORE;
			if (Pen > m_MaxPenalty)
				return -1;
			Score += MISMATCH_SCORE;
			int Drop = (BestScore - Score);
			if (Drop > XDROP)
				break;
			}
		}
	if (StartPos == 0 && EndPos == int(QL) -  1)
		{
		AddHitX(DBLo, Plus, BestScore, "");
		return BestScore;
		}
	else 
		{
		if (BestScore >= MinHSPScore)
			{
			unsigned HSPLength = unsigned(EndPos - StartPos + 1);
			AddHSPX(unsigned(StartPos), DBLo+unsigned(StartPos), Plus,
			  HSPLength, BestScore);
			return -2;
			}
		}
	return -1;
	}

bool State1::ExtendExact(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus)
	{
	if (SeedPosDB < SeedPosQ)
		return false;
	const byte *QSeq = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	const unsigned QL = m_Query->m_L;
	uint32 DBLo = SeedPosDB - SeedPosQ;
	const byte *DBSeq = m_UFI->m_SeqData + DBLo;

// Assume seed matches
	int EndPos = int(SeedPosQ+m_WordLength) - 1;
	for (int QPos = EndPos + 1; QPos < int(QL); ++QPos)
		{
		byte q = QSeq[QPos];
		byte t = DBSeq[QPos];
		if (q != t)
			return false;
		}

	int StartPos = int(SeedPosQ);
	for (int QPos = StartPos-1; QPos >= 0; --QPos)
		{
		byte q = QSeq[QPos];
		byte t = DBSeq[QPos];
		if (q != t)
			return false;
		}

	return true;
	}
