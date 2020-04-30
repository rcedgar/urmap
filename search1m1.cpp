#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "omplock.h"
#include <inttypes.h>

void State1::Search1(SeqInfo *Query)
	{
	m_Query = Query;

	m_HitCount = 0;
	m_HSPCount = 0;
	m_TopHit = 0;
	m_BestScore = 0;
	m_SecondBestScore = 0;
	m_Mapq = -1;

	Search1_Lo(Query);
	Output1();
	Unlock();
	}

bool State1::Search1_Lo(SeqInfo *Query)
	{
	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	unsigned QWordCount = QL - (m_WordLength - 1);

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
		bool Yes = ExtendExact(QPos, SeedPosDB, true);
		if (Yes)
			return true;
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
		bool Yes = ExtendExact(QPos, SeedPosDB, false);
		if (Yes)
			return true;
		}
		}

// ========================== Start Phase 2 =============================
// Extend BOTH1 seeds not tried in Phase 1
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		if (QPos%m_WordLength == 0)
			continue;

		{/// Plus
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
		bool Yes = ExtendExact(QPos, SeedPosDB, true);
		if (Yes)
			return true;
		}
	Minus2:
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
		bool Yes = ExtendExact(QPos, SeedPosDB, false);
		if (Yes)
			return true;
		}
	}
// ========================== End Phase 2 =============================

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
			bool Yes = ExtendExact(QPos, SeedPosDB, true);
			if (Yes)
				return true;
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
	//TODO: cache blob
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
			bool Yes = ExtendExact(QPos, SeedPosDB, false);
			if (Yes)
				return true;
			}
		} // end for QPos
	} // end minus

// ========================== Start Phase 5 =============================
// Extend all other seeds, i.e. with ab>=2
	{ // plus
	for (unsigned TodoIndex = 0; TodoIndex < TodoCount_Plus; ++TodoIndex)
		{
		uint32 QPos = QPosTodoVec_Plus[TodoIndex];
		uint64 Slot = SlotsVec_Plus[QPos];
		unsigned RowLength = m_UFI->GetRow_Blob(Slot, BlobVec_Plus + 5*QPos, PosVec);
		for (unsigned RowIndex = 0; RowIndex < RowLength; ++RowIndex)
			{
			uint32 SeedPosDB = PosVec[RowIndex];
			bool Yes = ExtendExact(QPos, SeedPosDB, true);
			if (Yes)
				return true;
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
			bool Yes = ExtendExact(QPos, SeedPosDB, false);
			if (Yes)
				return true;
			}
		} // end for TodoIndex
	} // end minus

	return false;
	}
