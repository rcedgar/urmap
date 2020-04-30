#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "omplock.h"
#include <inttypes.h>

unsigned State1::CalcMAPQ6() const
	{
	if (m_HitCount == 0)
		return 0;
	if (m_BestScore <= 0)
		return 0;

	unsigned QL = m_Query->m_L;
	double BestPossibleScore = double(QL);
	double Second = double(m_SecondBestScore);
	if (Second < BestPossibleScore/2.0)
		{
		Second = BestPossibleScore/2.0;
		if (m_BestScore <= Second)
			return 0;
		}
	double Fract = double(m_BestScore)/double(BestPossibleScore);
	double Drop = m_BestScore - Second;
	if (Drop > 40)
		Drop = 40;
	unsigned mapq = (unsigned) (Drop*Fract*Fract);
	if (mapq > 40)
		mapq = 40;
	return mapq;
	}

void State1::Search_Lo()
	{
	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	unsigned QWordCount = QL - (m_WordLength - 1);

	m_MaxPenalty = MAX_PENALTY;
	dpscore MinScorePhase1 = dpscore(QL) + XPHASE1*MISMATCH_SCORE;
	dpscore MinScorePhase3 = dpscore(QL) + XPHASE3*MISMATCH_SCORE;
	dpscore MinScorePhase4 = dpscore(QL) + XPHASE4*MISMATCH_SCORE;
	dpscore TermHSPScorePhase3 = (dpscore(QL)*TERM_HSP_SCORE_PCT_PHASE3)/100;

	m_BestHSPScore = 0;
	m_AllocBuffOffset = 0;

	uint64 *SlotsVec_Plus = alloc_buff(uint64, QL);
	uint64 *SlotsVec_Minus = alloc_buff(uint64, QL);
	byte *BlobVec_Plus = alloc_buff(byte, 5*QL);
	byte *BlobVec_Minus = alloc_buff(byte, 5*QL);

	SetSlotsVec(Q, QL, SlotsVec_Plus);

	m_RevCompQuerySeq.Alloc(QL);
	RevCompSeq(Q, QL, m_RevCompQuerySeq.Data);
	const byte *QRC = m_RevCompQuerySeq.Data;

	SetSlotsVec(QRC, QL, SlotsVec_Minus);

// ========================== Start Phase 1 =============================
// Extend BOTH1 seeds with stride = word length
// If gapless extension covering whole read with <= 1 diff, done.
// Save HSPs for later gapped extension.
	for (uint32 QPos = 0; QPos < QWordCount; QPos += m_WordLength)
		{
		{/// Plus
		uint64 Slot = SlotsVec_Plus[QPos];
		if (Slot == UINT64_MAX)
			{
			BlobVec_Plus[5*QPos] = TALLY_FREE;
			goto Minus1;
			}
		m_UFI->GetBlob(Slot, BlobVec_Plus + 5*QPos);
		byte T = BlobVec_Plus[5*QPos];
		if (T != TALLY_BOTH1)
			goto Minus1;
		uint32 SeedPosDB = *(uint32 *) (BlobVec_Plus + 5*QPos + 1);
		dpscore Score = ExtendPen(QPos, SeedPosDB, true);
		if (Score >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}
	Minus1:
		{///Minus
		uint64 Slot = SlotsVec_Minus[QPos];
		if (Slot == UINT64_MAX)
			{
			BlobVec_Minus[5*QPos] = TALLY_FREE;
			continue;
			}
		m_UFI->GetBlob(Slot, BlobVec_Minus + 5*QPos);
		byte T = BlobVec_Minus[5*QPos];
		if (T != TALLY_BOTH1)
			continue;
		uint32 SeedPosDB = *(uint32 *) (BlobVec_Minus + 5*QPos + 1);
		dpscore Score = ExtendPen(QPos, SeedPosDB, false);
		if (Score >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}
		}

// ========================== Start Phase 2 =============================
// Extend BOTH1 seeds not tried in Phase 1
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		if (QPos%m_WordLength == 0)
			continue;

		{// Plus
		uint64 Slot = SlotsVec_Plus[QPos];
		if (Slot == UINT64_MAX)
			{
			BlobVec_Plus[5*QPos] = TALLY_FREE;
			goto Minus2;
			}
		m_UFI->GetBlob(Slot, BlobVec_Plus + 5*QPos);
		byte T = BlobVec_Plus[5*QPos];
		if (T != TALLY_BOTH1)
			goto Minus2;
		uint32 SeedPosDB = *(uint32 *) (BlobVec_Plus + 5*QPos + 1);
		dpscore Score = ExtendPen(QPos, SeedPosDB, true);
		if (Score >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}
	Minus2:
		{// Minus
		uint64 Slot = SlotsVec_Minus[QPos];
		if (Slot == UINT64_MAX)
			{
			BlobVec_Minus[5*QPos] = TALLY_FREE;
			continue;
			}
		m_UFI->GetBlob(Slot, BlobVec_Minus + 5*QPos);
		byte T = BlobVec_Minus[5*QPos];
		if (T != TALLY_BOTH1)
			continue;
		uint32 SeedPosDB = *(uint32 *) (BlobVec_Minus + 5*QPos + 1);
		dpscore Score = ExtendPen(QPos, SeedPosDB, false);
		if (Score >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}
		}
// ========================== End Phase 2 =============================

// ========================== Start Phase 3 =============================
// Extend high-scoring HSPs from Phases 1 & 2
//	if (m_BestHSPScore > dpscore(QL)/2)
	if (m_BestHSPScore > TermHSPScorePhase3)
		{
		for (unsigned i = 0; i < m_HSPCount; ++i)
			AlignHSP(i);
		if (m_BestScore >= MinScorePhase1)
			{
			m_Mapq = CalcMAPQ6();
			return;
			}
		}
// ========================== End Phase 3 =============================

// ========================== Start Phase 4 =============================
// Extend seeds with ab<=2 except BOTH1.
// Save HSPs for later gapped extension.
	uint32 *PosVec = alloc_buff(uint32, m_UFI->m_MaxIx);
	unsigned TodoCount_Plus = 0;
	uint32 *QPosTodoVec_Plus = alloc_buff(uint32, QL);
	{//plus
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		byte T = BlobVec_Plus[5*QPos];
		if (T == TALLY_FREE || T == TALLY_BOTH1 || TallyOther(T))
			continue;
		uint64 Slot = SlotsVec_Plus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Plus + 5*QPos, PosVec);
		if (RowLength > 2)
			{
			QPosTodoVec_Plus[TodoCount_Plus] = QPos;
			++TodoCount_Plus;
			continue;
			}
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, true);
			}
		} // end for QPos
	} // end plus

	unsigned TodoCount_Minus = 0;
	uint32 *QPosTodoVec_Minus = alloc_buff(uint32, QL);
	{//minus
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		byte T = BlobVec_Minus[5*QPos];
		if (T == TALLY_FREE || T == TALLY_BOTH1 || TallyOther(T))
			continue;
		uint64 Slot = SlotsVec_Minus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Minus + 5*QPos, PosVec);
		if (RowLength > 2)
			{
			QPosTodoVec_Minus[TodoCount_Minus] = QPos;
			++TodoCount_Minus;
			continue;
			}
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, false);
			}
		} // end for QPos
	} // end minus

	if (m_BestScore >= MinScorePhase3)
		{
		m_Mapq = CalcMAPQ6();
		return;
		}
// ========================== End Phase 4 =============================

// ========================== Start Phase 5 =============================
// Extend all other seeds, i.e. with ab>=2
// Save HSPs for later gapped extension.
	{ // plus
	for (unsigned TodoIndex = 0; TodoIndex < TodoCount_Plus; ++TodoIndex)
		{
		uint32 QPos = QPosTodoVec_Plus[TodoIndex];
		uint64 Slot = SlotsVec_Plus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Plus + 5*QPos, PosVec);
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, true);
			}
		} // end for TodoIndex
	} // end plus

	{ // minus
	for (unsigned TodoIndex = 0; TodoIndex < TodoCount_Minus; ++TodoIndex)
		{
		uint32 QPos = QPosTodoVec_Minus[TodoIndex];
		uint64 Slot = SlotsVec_Minus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Minus + 5*QPos, PosVec);
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			ExtendPen(QPos, SeedPosDB, false);
			}
		} // end for TodoIndex
	} // end minus

// ========================== End Phase 5 =============================
	if (m_BestScore >= MinScorePhase4)
		{
		m_Mapq = CalcMAPQ6();
		return;
		}

// ========================== Start Phase 6 =============================
// Align HSPs found in previous phases
	for (unsigned i = 0; i < m_HSPCount; ++i)
		AlignHSP(i);
// ========================== End Phase 6 =============================
	m_Mapq = CalcMAPQ6();
	}
