#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "omplock.h"
#include <inttypes.h>

unsigned State1::GetFirstBoth1Seed(uint32 &QPos, bool &Plus, uint32 &DBPos)
	{
	unsigned QL = m_Query->m_L;
	unsigned QueryWordCount = QL - (m_WordLength - 1);
	for (unsigned k = 0; k < QueryWordCount; ++k)
		{
		QPos = (k*PRIME_STRIDE)%QueryWordCount;
		{// Plus
		uint64 Slot = m_SlotsVec_Plus[QPos];
		if (Slot == UINT64_MAX)
			goto Minus1;
		m_UFI->GetBlob(Slot, m_BlobVec_Plus + 5*QPos);
		byte T = m_BlobVec_Plus[5*QPos];
		if (TallyOther(T))
			goto Minus1;
		if (T != TALLY_BOTH1)
			{
			m_QPosPendingVec_Plus[m_QPosPendingCount_Plus++] = QPos;
			goto Minus1;
			}
		DBPos = *(uint32 *) (m_BlobVec_Plus + 5*QPos + 1);
		Plus = true;
		return k;
		}
	Minus1:
		{// Minus
		uint64 Slot = m_SlotsVec_Minus[QPos];
		if (Slot == UINT64_MAX)
			continue;
		m_UFI->GetBlob(Slot, m_BlobVec_Minus + 5*QPos);
		byte T = m_BlobVec_Minus[5*QPos];
		if (TallyOther(T))
			continue;
		if (T != TALLY_BOTH1)
			{
			m_QPosPendingVec_Minus[m_QPosPendingCount_Minus++] = QPos;
			continue;
			}
		DBPos = *(uint32 *) (m_BlobVec_Minus + 5*QPos + 1);
		Plus = false;
		return k;
		}
		}

	return UINT_MAX;
	}

unsigned State1::GetNextBoth1Seed(unsigned ak,
  uint32 &aQPos, bool &Plus, uint32 &DBPos)
	{
	unsigned QL = m_Query->m_L;
	unsigned QueryWordCount = QL - (m_WordLength - 1);
	{
// Special case: if last BOTH1 seed was Plus=true, need to check 
// minus strand at same position before incrementing k.
	if (!Plus)
		goto Skip;
	unsigned QPos = (ak*PRIME_STRIDE)%QueryWordCount;
	uint64 Slot = m_SlotsVec_Minus[QPos];
	if (Slot == UINT64_MAX)
		goto Skip;
	m_UFI->GetBlob(Slot, m_BlobVec_Minus + 5*QPos);
	byte T = m_BlobVec_Minus[5*QPos];
	if (!TallyMine(T))
		goto Skip;
	if (T == TALLY_BOTH1)
		{
		uint32 NewDBPos = *(uint32 *) (m_BlobVec_Minus + 5*QPos + 1);
		if (NewDBPos - QPos != DBPos - aQPos)
			{
			DBPos = NewDBPos;
			aQPos = QPos;
			Plus = false;
			return ak;
			}
		else
			m_QPosPendingVec_Minus[m_QPosPendingCount_Minus++] = QPos;
		}
	}
Skip:
	for (unsigned k = ak + 1; k < QueryWordCount; ++k)
		{
		unsigned QPos = (k*PRIME_STRIDE)%QueryWordCount;
		{// Plus
		uint64 Slot = m_SlotsVec_Plus[QPos];
		if (Slot == UINT64_MAX)
			goto Minus1;
		m_UFI->GetBlob(Slot, m_BlobVec_Plus + 5*QPos);
		byte T = m_BlobVec_Plus[5*QPos];
		if (TallyOther(T))
			goto Minus1;
		if (T != TALLY_BOTH1)
			{
			m_QPosPendingVec_Plus[m_QPosPendingCount_Plus++] = QPos;
			goto Minus1;
			}
		uint32 NewDBPos = *(uint32 *) (m_BlobVec_Plus + 5*QPos + 1);
		if (NewDBPos - QPos == DBPos - aQPos)
			goto Minus1;
		DBPos = NewDBPos;
		aQPos = QPos;
		Plus = true;
		return k;
		}
	Minus1:
		{// Minus
		uint64 Slot = m_SlotsVec_Minus[QPos];
		if (Slot == UINT64_MAX)
			continue;
		m_UFI->GetBlob(Slot, m_BlobVec_Minus + 5*QPos);
		byte T = m_BlobVec_Minus[5*QPos];
		if (TallyOther(T))
			continue;
		if (T != TALLY_BOTH1)
			{
			m_QPosPendingVec_Minus[m_QPosPendingCount_Minus++] = QPos;
			continue;
			}
		uint32 NewDBPos = *(uint32 *) (m_BlobVec_Minus + 5*QPos + 1);
		if (NewDBPos - QPos == DBPos - aQPos)
			continue;
		DBPos = NewDBPos;
		aQPos = QPos;
		Plus = false;
		return k;
		}
		}

	return UINT32_MAX;
	}

void State1::GetBoth1Seeds(vector<unsigned> &QPosVec,
  vector<uint32> &DBPosVec, vector<bool> &PlusVec,  vector<uint64> &DiagVec)
	{
	QPosVec.clear();
	DBPosVec.clear();
	PlusVec.clear();
	DiagVec.clear();

	unsigned QL = m_Query->m_L;
	unsigned QueryWordCount = QL - (m_WordLength - 1);

	for (unsigned QPos = 0; QPos < QueryWordCount; ++QPos)
		{
		for (int iPlus = 0; iPlus < 2; ++iPlus)
			{
			bool Plus = (iPlus == 0);
			const uint64 *SlotsVec = (Plus ? m_SlotsVec_Plus : m_SlotsVec_Minus);
			uint64 Slot = SlotsVec[QPos];
			if (Slot == UINT64_MAX)
				continue;
			byte T = m_UFI->GetTally(Slot);
			if (T != TALLY_BOTH1)
				continue;
			unsigned RowLength = m_UFI->GetRow(Slot, m_PosVec);
			asserta(RowLength == 1);
			uint32 DBPos = m_PosVec[0];
			uint64 Diag = uint64(DBPos) - uint64(QPos);
			bool Found = false;
			for (unsigned i = 0; i < SIZE(DiagVec); ++i)
				if (DiagVec[i] == Diag)
					{
					Found = true;
					break;
					}
			if (!Found)
				DiagVec.push_back(Diag);

			QPosVec.push_back(QPos);
			DBPosVec.push_back(DBPos);
			PlusVec.push_back(Plus);
			}
		}
	}

unsigned State1::GetFirstBoth1SeedEx(uint32 &QPos,
  bool &Plus, uint32 &DBPos)
	{
	unsigned k = GetFirstBoth1Seed(QPos, Plus, DBPos);
	while (k != UINT_MAX)
		{
		int Score = ExtendPen(QPos, DBPos, Plus);
		if (Score != -1)
			return k;
		k = GetNextBoth1Seed(k, QPos, Plus, DBPos);
		}
	return UINT_MAX;
	}

unsigned State1::GetNextBoth1SeedEx(unsigned k, uint32 &QPos,
  bool &Plus, uint32 &DBPos)
	{
	uint64 LastDiag = uint64(DBPos) - uint64(QPos);
	for (;;)
		{
		k = GetNextBoth1Seed(k, QPos, Plus, DBPos);
		if (k == UINT_MAX)
			return UINT_MAX;
		uint64 Diag = uint64(DBPos) - uint64(QPos);
		if (Diag == LastDiag)
			continue;
		int Score = ExtendPen(QPos, DBPos, Plus);
		if (Score != -1)
			return k;
		}
	return UINT_MAX;
	}
