#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state2.h"
#include "ufihit.h"
#include "alignresult.h"
#include "sort.h"
#include "omplock.h"

unsigned State2::m_Method;
unsigned State2::m_Minq;

void State2::SetUFI(const UFIndex &UFI)
	{
	m_UD_Fwd.SetUFI(UFI);
	m_UD_Rev.SetUFI(UFI);
	}

void State2::FindPairs()
	{
	m_Paired_HitIndexVec_Fwd.clear();
	m_Paired_HitIndexVec_Rev.clear();
	m_Paired_TotalScoreVec.clear();
	m_Paired_TLVec.clear();

	const unsigned QL2 = (m_UD_Fwd.m_Query->m_L + m_UD_Rev.m_Query->m_L)/2;

	const unsigned HitCount_Fwd = m_UD_Fwd.m_HitCount;
	const unsigned HitCount_Rev = m_UD_Rev.m_HitCount;

	m_BestPairIndex = UINT_MAX;
	m_SecondPairIndex = UINT_MAX;

	m_BestPairScore = -1;
	m_SecondBestPairScore = -1;

	for (unsigned HitIndex_Fwd = 0; HitIndex_Fwd < HitCount_Fwd; ++HitIndex_Fwd)
		{
		const UFIHit *Hit_Fwd = m_UD_Fwd.GetHit(HitIndex_Fwd);
		int Score_Fwd = Hit_Fwd->m_Score;
		int64 DBPos_Fwd = (int64) Hit_Fwd->m_DBStartPos;
		if (Score_Fwd < m_UD_Fwd.m_SecondBestScore - 12)//@@TODO:param
			continue;
		bool Plus_Fwd = Hit_Fwd->m_Plus;
		for (unsigned HitIndex_Rev = 0; HitIndex_Rev < HitCount_Rev; ++HitIndex_Rev)
			{
			const UFIHit *Hit_Rev = m_UD_Rev.GetHit(HitIndex_Rev);
			int Score_Rev = Hit_Rev->m_Score;
			if (Score_Rev < m_UD_Rev.m_SecondBestScore - 12)//@@TODO:param
				continue;
			int64 DBPos_Rev = (int64) Hit_Rev->m_DBStartPos;
			int64 TL = abs(DBPos_Fwd - DBPos_Rev) + int64(QL2);
			if (TL > 1000) // TODO:param/adaptive
				continue;
			bool Plus_Rev = Hit_Rev->m_Plus;
			if (Plus_Rev == Plus_Fwd)
				continue;

			int TotalScore = Score_Fwd + Score_Rev;
			if (TotalScore > m_BestPairScore)
				{
				m_SecondPairIndex = m_BestPairIndex;
				m_SecondBestPairScore = m_BestPairScore;
				m_BestPairScore = TotalScore;
				m_BestPairIndex = SIZE(m_Paired_TotalScoreVec);
				}
			else if (TotalScore == m_BestPairScore)
				{
				m_SecondPairIndex = SIZE(m_Paired_TotalScoreVec);
				m_SecondBestPairScore = m_BestPairScore;
				}
			else if (TotalScore > m_SecondBestPairScore)
				{
				m_SecondPairIndex = m_BestPairIndex;
				m_SecondBestPairScore = TotalScore;
				}

			m_Paired_TotalScoreVec.push_back(TotalScore);
			m_Paired_HitIndexVec_Fwd.push_back(HitIndex_Fwd);
			m_Paired_HitIndexVec_Rev.push_back(HitIndex_Rev);
			m_Paired_TLVec.push_back(unsigned(TL));
			}
		}
	}

void State2::ScanPair()
	{
	const unsigned HitCount_Fwd = m_UD_Fwd.m_HitCount;
	const unsigned HitCount_Rev = m_UD_Rev.m_HitCount;
	const unsigned SCAN_DB_SEG_LENGTH = 1024;//TODO:PARAM

	int SavedMapq_Fwd = m_UD_Fwd.m_Mapq;
	int SavedMapq_Rev = m_UD_Rev.m_Mapq;
	bool DoVit_Fwd = (SavedMapq_Fwd >= 10);
	bool DoVit_Rev = (SavedMapq_Rev >= 10);

	//@@TODO:check underflow
	for (unsigned HitIndex_Fwd = 0; HitIndex_Fwd < HitCount_Fwd; ++HitIndex_Fwd)
		{
		unsigned QL_Rev = m_UD_Fwd.m_Query->m_L;
		const UFIHit *Hit_Fwd = m_UD_Fwd.GetHit(HitIndex_Fwd);
		if (Hit_Fwd->m_Score < m_UD_Fwd.m_SecondBestScore)//@@TODO:param
			continue;
		uint32 DBPos = Hit_Fwd->m_DBStartPos;
		bool Plus = Hit_Fwd->m_Plus;
		if (Plus)
			m_UD_Rev.Scan(DBPos, SCAN_DB_SEG_LENGTH, false, DoVit_Fwd);
		else
			{
			if (DBPos >= SCAN_DB_SEG_LENGTH)
				m_UD_Rev.Scan(DBPos-SCAN_DB_SEG_LENGTH, SCAN_DB_SEG_LENGTH+2*QL_Rev,
				  true, DoVit_Fwd);
			}
		}

	for (unsigned HitIndex_Rev = 0; HitIndex_Rev < HitCount_Rev; ++HitIndex_Rev)
		{
		unsigned QL_Fwd = m_UD_Fwd.m_Query->m_L;
		const UFIHit *Hit_Rev = m_UD_Rev.GetHit(HitIndex_Rev);
		if (Hit_Rev->m_Score < m_UD_Rev.m_SecondBestScore)//@@TODO:param
			continue;
		uint32 DBPos = Hit_Rev->m_DBStartPos;
		bool Plus = Hit_Rev->m_Plus;
		if (Plus)
			m_UD_Fwd.Scan(DBPos, SCAN_DB_SEG_LENGTH, false, DoVit_Rev);
		else
			{
			if (DBPos >= SCAN_DB_SEG_LENGTH)
				m_UD_Fwd.Scan(DBPos-SCAN_DB_SEG_LENGTH, SCAN_DB_SEG_LENGTH+2*QL_Fwd,
				  true, DoVit_Rev);
			}
		}

	m_UD_Fwd.CalcMAPQ6();
	m_UD_Rev.CalcMAPQ6();
	}

void State2::LogPairs() const
	{
	const unsigned PairCount = SIZE(m_Paired_TotalScoreVec);
	Log("\n");
	Log("LogPairs, %u pairs\n", PairCount);
	if (PairCount == 0)
		return;

	vector<unsigned> Order(PairCount);
	QuickSortOrderDesc(m_Paired_TotalScoreVec.data(), PairCount, Order.data());

	Log("Best score %.1f, second %.1f\n",
	  double(m_BestPairScore), double(m_SecondBestPairScore));
	for (unsigned k = 0; k < PairCount; ++k)
		{
		unsigned i = Order[k];
		unsigned fi = m_Paired_HitIndexVec_Fwd[i];
		unsigned ri = m_Paired_HitIndexVec_Rev[i];
		int TotalScore = m_Paired_TotalScoreVec[i];
		const UFIHit *FwdHit =  m_UD_Fwd.GetHit(fi);
		const UFIHit *RevHit =  m_UD_Rev.GetHit(ri);
		int Score_Fwd = FwdHit->m_Score;
		int Score_Rev = RevHit->m_Score;
		unsigned FwdDBPos = FwdHit->m_DBStartPos;
		unsigned RevDBPos = RevHit->m_DBStartPos;
		unsigned Offset = (FwdDBPos >= RevDBPos ?
		  FwdDBPos - RevDBPos : RevDBPos - FwdDBPos);

		//Log("%10u", m_DBPosVec[i]);
		Log("  %6.2f", double(TotalScore));
		Log("  %6.1f", double(Score_Fwd));
		Log("  %6.2f", double(Score_Rev));
		Log("  %7u", Offset);
		Log("  Fwd=");
		m_UD_Fwd.LogPos(FwdDBPos);
		Log(", Rev=");
		m_UD_Rev.LogPos(RevDBPos);
		Log("\n");
		}
	}

void State2::GetChrPos(const UFIHit *Hit, string &ChrLabel, uint32 &ChrPos) const
	{
	m_UD_Fwd.GetChrPos(Hit, ChrLabel, ChrPos);
	}

bool State2::AcceptTemplateLength(int64 TL) const
	{
	return TL <= MAX_TL;
	}
