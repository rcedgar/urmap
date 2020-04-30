#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "omplock.h"
#include <inttypes.h>

void State1::GetBoth1s(SeqInfo *Query, vector<uint32> &PosVec)
	{
	PosVec.clear();
	m_Query = Query;

	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	unsigned QWordCount = QL - (m_WordLength - 1);

	m_AllocBuffOffset = 0;

	uint64 *SlotsVec_Plus = alloc_buff(uint64, QL);
	uint64 *SlotsVec_Minus = alloc_buff(uint64, QL);
	byte *BlobVec_Plus = alloc_buff(byte, 5*QL);
	byte *BlobVec_Minus = alloc_buff(byte, 5*QL);

// Plus
	SetSlotsVec(Q, QL, SlotsVec_Plus);
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		uint64 Slot = SlotsVec_Plus[QPos];
		if (Slot == UINT64_MAX)
			continue;

		m_UFI->GetBlob(Slot, BlobVec_Plus + 5*QPos);
		byte T = BlobVec_Plus[5*QPos];
		if (T == TALLY_BOTH1)
			{
			uint32 Pos = *(uint32 *) (BlobVec_Plus + 5*QPos + 1);
			PosVec.push_back(Pos);
			}
		}

// Minus
	m_RevCompQuerySeq.Alloc(QL);
	RevCompSeq(Q, QL, m_RevCompQuerySeq.Data);
	const byte *QRC = m_RevCompQuerySeq.Data;
	SetSlotsVec(QRC, QL, SlotsVec_Minus);
	for (uint32 QPos = 0; QPos < QWordCount; ++QPos)
		{
		uint64 Slot = SlotsVec_Minus[QPos];
		if (Slot == UINT64_MAX)
			continue;

		m_UFI->GetBlob(Slot, BlobVec_Minus + 5*QPos);
		byte T = BlobVec_Minus[5*QPos];
		if (T == TALLY_BOTH1)
			{
			uint32 Pos = *(uint32 *) (BlobVec_Minus + 5*QPos + 1);
			PosVec.push_back(Pos);
			}
		}
	}
