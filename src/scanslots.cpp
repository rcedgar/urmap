#include "myutils.h"
#include "alpha.h"
#include "seqinfo.h"
#include "ufindex.h"
#include "state1.h"

void State1::ScanSlots(uint32 DBLo, unsigned DBSegLength, bool Plus)
	{
	unsigned QL = m_Query->m_L;
	unsigned QueryWordCount = QL - (m_WordLength - 1);
	const uint64 *SlotsVec = (Plus ? m_SlotsVec_Plus : m_SlotsVec_Minus);
	const byte *BlobVec = (Plus ? m_BlobVec_Plus : m_BlobVec_Minus);
	if (QL <= m_WordLength*4)
		return;

	const byte *DBSegSeq = m_UFI->m_SeqData + DBLo;
	const byte *Q = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	uint64 Word = 0;
	byte K = 0;
	for (uint32 DBPos = 0; DBPos < m_WordLength-1; ++DBPos)
		{
		byte c = DBSegSeq[DBPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		}

	for (uint32 DBPos = m_WordLength-1; DBPos < DBSegLength; ++DBPos)
		{
		byte c = DBSegSeq[DBPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			for (unsigned k = 0; k < SCANK; ++k)
				{
				unsigned QPos = (k*PRIME_STRIDE)%QueryWordCount;
				if (Slot == SlotsVec[QPos])
					{
					unsigned DBSeedPos = DBLo + DBPos - m_WordLength + 1;
					ExtendScan(QPos, DBSeedPos, Plus);
					}
				}
			}
		}
	}
