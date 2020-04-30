#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "alignresult.h"
#include "pathinfo.h"
#include "alnparams.h"
#include <inttypes.h>

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);

void State1::Scan(uint32 DBPos, unsigned DBSegLength, bool Plus, bool DoVit)
	{
	int SavedMaxPenalty = m_MaxPenalty;
	unsigned SavedHitCount = m_HitCount;
	m_MaxPenalty = 130;//@@TODO
	ScanSlots(DBPos, DBSegLength, Plus);
	m_MaxPenalty = SavedMaxPenalty;
	if (m_HitCount > SavedHitCount)
		return;
	if (!DoVit)
		return;

	const byte *Q = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	const unsigned QL = m_Query->m_L;
	const byte *T = m_UFI->m_SeqData + DBPos;
	const unsigned TL = DBSegLength;
	float Score = Viterbi(Q, QL, T, TL, true, true, *m_PI);

	if (Score >= QL/3.0) // TODO:PARAM
		{
		unsigned LeftICount = m_PI->TrimLeftIs();
		m_PI->TrimRightIs();
		unsigned StartPosDB = DBPos + LeftICount;
		AddHitX(StartPosDB, Plus, int(Score), m_PI->GetPath());
		}
	}
