#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "alignresult.h"
#include "pathinfo.h"
#include "alnparams.h"
#include <inttypes.h>

#define	TRACE	0
#define BRN		1//TODO:PARAM

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);

unsigned State1::AlignHSP5(unsigned HSPIndex)
	{
	UFIHSP *HSP = m_HSPs.Data[HSPIndex];
	if (HSP->m_Aligned)
		return UINT_MAX;

	HSP->m_Aligned = true;
	dpscore TotalPen = dpscore(HSP->m_Length) - HSP->m_Score;
	dpscore TotalScore = HSP->m_Score;
	if (TotalPen > m_MaxPenalty)
		return UINT_MAX;

	const unsigned StartPosQ = HSP->m_StartPosQ;
	const unsigned StartPosDB = HSP->m_StartPosDB;
	const unsigned HSPLength = HSP->m_Length;
	const bool Plus = HSP->m_Plus;
	const unsigned QL = m_Query->m_L;
	const unsigned TL = m_UFI->m_SeqDataSize;
	unsigned CombinedTLo = StartPosDB;

	const byte *Q = (HSP->m_Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	const byte *T = m_UFI->m_SeqData;
	m_LeftPI->SetEmpty();
	m_RightPI->SetEmpty();

/***
        0    StartPosQ               QL-1
  LeftQLo     LeftQHi     RightQLo   RightQHi
        |       |            |       |
   Q:   LEFTQ   qqqqq HSP qqqq  RIGHTQ
   T:   LEFTT   ttttt HSP tttt  RIGHTT
                |            |
             LeftTHi      RightTLo
            StartPosDB
***/

/////////////////////////////// Left /////////////////////////
// TODO:optimize by doing right first?
	if (StartPosQ > 0)
		{
		if (StartPosDB < StartPosQ)
			return UINT_MAX;

		unsigned LeftQLo = 0;
		unsigned LeftQHi = StartPosQ - 1;
		unsigned LeftQL = LeftQHi - LeftQLo + 1;

		unsigned LeftTHi = StartPosDB - 1;
		unsigned LeftTL = LeftQL + BRN*GLOBAL_BAND_RADIUS;
		if (LeftTL >= LeftTHi)
			return UINT_MAX;
		unsigned LeftTLo = LeftTHi - LeftTL + 1;

		const byte *LeftQ = Q;
		const byte *LeftT = T + LeftTLo;
		dpscore LeftScore = (dpscore) Viterbi(LeftQ, LeftQL, LeftT, LeftTL, true, false, *m_LeftPI);
		unsigned LeftICount = m_LeftPI->TrimLeftIs();
		CombinedTLo = LeftTLo + LeftICount;
		dpscore AllGapScore = GAP_OPEN_SCORE + (LeftQL-1)*GAP_EXT_SCORE;
		if (AllGapScore > LeftScore)
			LeftScore = AllGapScore;
		dpscore BestPossibleScore = dpscore(LeftQL);
		dpscore LeftPen = BestPossibleScore - LeftScore;
		TotalScore += LeftScore;
		TotalPen += LeftPen;
		if (TotalPen > m_MaxPenalty)
			return UINT_MAX;
		} // end Left

/////////////////////////////// Right /////////////////////////
	const unsigned RightQLo = StartPosQ + HSPLength;
	if (RightQLo < QL)
		{
		unsigned RightQHi = QL - 1;
		asserta(RightQHi >= RightQLo);
		unsigned RightQL = RightQHi - RightQLo + 1;

		unsigned RightTLo = StartPosDB + HSPLength;
		unsigned RightTHi = RightTLo + RightQL + BRN*GLOBAL_BAND_RADIUS;
		if (RightTHi >= TL)
			RightTHi = TL - 1;
		unsigned RightTL = RightTHi - RightTLo + 1;

		const byte *RightQ = Q + RightQLo;
		const byte *RightT = T + RightTLo;
		dpscore RightScore = (dpscore) Viterbi(RightQ, RightQL, RightT, RightTL, false, true, *m_RightPI);
		m_RightPI->TrimRightIs();
		dpscore AllGapScore = GAP_OPEN_SCORE + (RightQL-1)*GAP_EXT_SCORE;
		if (AllGapScore > RightScore)
			RightScore = AllGapScore;
		dpscore BestPossibleScore = dpscore(RightQL);
		dpscore RightPen = BestPossibleScore - RightScore;
		TotalScore += RightScore;
		TotalPen += RightPen;
		if (TotalPen > m_MaxPenalty)
			return UINT_MAX;
		} // end Right

	m_PI->SetEmpty();
	m_PI->AppendPath(*m_LeftPI);
	m_PI->AppendMs(HSPLength);
	m_PI->AppendPath(*m_RightPI);

	string Path = string(m_PI->GetPath());
	unsigned HitIndex = AddHitX(CombinedTLo, Plus, TotalScore, Path);
	return HitIndex;
	}
