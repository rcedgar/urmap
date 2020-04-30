#include "myutils.h"
#include "ufindex.h"
#include "state2.h"
#include "seqinfo.h"
#include "omplock.h"
#include "sort.h"

void State2::AdjustTopHitsAndMapqs()
	{
	const unsigned PairCount= SIZE(m_Paired_TotalScoreVec);
	if (PairCount == 0)
		{
		m_UD_Fwd.m_Mapq /= 2;
		m_UD_Rev.m_Mapq /= 2;
		return;
		}
	// BAD IDEA!
	//if (PairCount == 1)
	//	{
	//	m_UD_Fwd.m_Mapq = 30;
	//	m_UD_Rev.m_Mapq = 30;
	//	return;
	//	}

	unsigned QL = m_UD_Fwd.m_Query->m_L + m_UD_Rev.m_Query->m_L;
	double BestPossibleScore = double(QL);
	double Fract = double(m_BestPairScore)/double(BestPossibleScore);
	double Drop = m_BestPairScore - m_SecondBestPairScore;
	if (Drop > 30)
		Drop = 30;
	unsigned mapq = (unsigned)(Drop*Fract*Fract);
	if (mapq > 40)
		mapq = 40;

	if (mapq > m_UD_Fwd.m_Mapq)
		m_UD_Fwd.m_Mapq = mapq;
	if (mapq > m_UD_Rev.m_Mapq)
		m_UD_Rev.m_Mapq = mapq;

	if (m_BestPairIndex != UINT_MAX)
		{
		unsigned fi = m_Paired_HitIndexVec_Fwd[m_BestPairIndex];
		unsigned ri = m_Paired_HitIndexVec_Rev[m_BestPairIndex];

		m_UD_Fwd.m_TopHit = m_UD_Fwd.GetHit(fi);
		m_UD_Rev.m_TopHit = m_UD_Rev.GetHit(ri);
		}

	if (m_SecondPairIndex != UINT_MAX)
		{
		unsigned fi = m_Paired_HitIndexVec_Fwd[m_SecondPairIndex];
		unsigned ri = m_Paired_HitIndexVec_Rev[m_SecondPairIndex];

		m_UD_Fwd.m_SecondHit = m_UD_Fwd.GetHit(fi);
		m_UD_Rev.m_SecondHit = m_UD_Rev.GetHit(ri);
		}
	}

void State2::Search(SeqInfo *Query_Fwd, SeqInfo *Query_Rev)
	{
	switch (m_Method)
		{
	case 4:
		Search4(Query_Fwd, Query_Rev);
		break;
	case 5:
		Search5(Query_Fwd, Query_Rev);
		break;
	default:
		asserta(false);
		}
	Output2();
	}
