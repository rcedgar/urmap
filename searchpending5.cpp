#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "omplock.h"
#include <inttypes.h>

void State1::SearchPE_Pending5(unsigned k)
	{
	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	unsigned QWordCount = QL - (m_WordLength - 1);

	m_MaxPenalty = MAX_PENALTY;
	dpscore MinScorePhase1 = dpscore(QL) + XPHASE1*MISMATCH_SCORE;
	dpscore MinScorePhase3 = dpscore(QL) + XPHASE3*MISMATCH_SCORE;
	dpscore MinScorePhase4 = dpscore(QL) + XPHASE4*MISMATCH_SCORE;
	dpscore TermHSPScorePhase3 = (dpscore(QL)*TERM_HSP_SCORE_PCT_PHASE3)/100;

	uint32 QPos;
	bool Plus;
	uint32 DBPos;
	while (k != UINT_MAX)
		k = GetNextBoth1SeedEx(k, QPos, Plus, DBPos);
	if (m_BestScore >= MinScorePhase1)
		{
		m_Mapq = CalcMAPQ6();
		return;
		}
	if (m_BestHSPScore >= TermHSPScorePhase3)
		{
		for (unsigned i = 0; i < m_HSPCount; ++i)
			AlignHSP5(i);
		if (m_BestScore >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}

	const uint64 *SlotsVec_Plus = m_SlotsVec_Plus;
	const uint64 *SlotsVec_Minus = m_SlotsVec_Minus;
	byte *BlobVec_Plus = m_BlobVec_Plus;
	byte *BlobVec_Minus = m_BlobVec_Minus;
	uint32 *PosVec = m_PosVec;

// Pending Plus round 1
	unsigned QPosPendingCount_Plus2 = 0;
	for (unsigned i = 0; i < m_QPosPendingCount_Plus; ++i)
		{
		unsigned QPos = m_QPosPendingVec_Plus[i];
		uint64 Slot = SlotsVec_Plus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Plus + 5*QPos, PosVec);
		if (RowLength > 2)
			{
			m_QPosPendingVec_Plus[QPosPendingCount_Plus2++] = QPos;
			continue;
			}
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, true);
			}
		}

// Pending Minus round 1
	unsigned QPosPendingCount_Minus2 = 0;
	for (unsigned i = 0; i < m_QPosPendingCount_Minus; ++i)
		{
		unsigned QPos = m_QPosPendingVec_Minus[i];
		uint64 Slot = SlotsVec_Minus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Minus + 5*QPos, PosVec);
		if (RowLength > 2)
			{
			m_QPosPendingVec_Minus[QPosPendingCount_Minus2++] = QPos;
			continue;
			}
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, false);
			}
		}

// Pending Plus round 2
	for (unsigned i = 0; i < QPosPendingCount_Plus2; ++i)
		{
		unsigned QPos = m_QPosPendingVec_Plus[i];
		uint64 Slot = SlotsVec_Plus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Plus + 5*QPos, PosVec);
		asserta(RowLength > 2);
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, true);
			}
		}

// Pending Minus round 2
	for (unsigned i = 0; i < QPosPendingCount_Minus2; ++i)
		{
		unsigned QPos = m_QPosPendingVec_Minus[i];
		uint64 Slot = SlotsVec_Minus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Minus + 5*QPos, PosVec);
		asserta(RowLength > 2);
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, false);
			}
		}

	int B = max(m_BestScore, m_BestHSPScore) - 8;
	for (unsigned i = 0; i < m_HSPCount; ++i)
		{
		const UFIHSP *HSP = GetHSP(i);
		if (HSP->m_Score < B)
			continue;
		AlignHSP5(i);
		}

	m_Mapq = CalcMAPQ6();
	}
