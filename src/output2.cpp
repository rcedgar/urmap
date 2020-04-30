#include "myutils.h"
#include "seqinfo.h"
#include "state1.h"
#include "state2.h"
#include "omplock.h"
#include "outfiles.h"

static int64 g_PrevSeqIndex;

void State2::Output2()
	{
	OutputSAM2();
	OutputTab2();
	m_UD_Fwd.UpdateHitStats();
	m_UD_Rev.UpdateHitStats();
	}

uint32 State2::GetPairedFlags(bool First, bool RevComp, bool MateRevComp, bool MateUnmapped)
	{
	uint32 Flags = 0;

	if (First)
		Flags |= 0x41;
	else
		Flags |= 0x81;

	if (RevComp)
		Flags |= 0x10;

	if (MateUnmapped)
		Flags |= 0x08;
	else if (MateRevComp)
		Flags |= 0x20;

	return Flags;
	}

void State2::OutputSAM2()
	{
	if (g_fSAM == 0)
		return;
	SetSAM2();
	Lock();
	fputs(m_UD_Fwd.m_SAM.Data, g_fSAM);
	fputs(m_UD_Rev.m_SAM.Data, g_fSAM);
	Unlock();
	}

unsigned State2::GetTemplateLength(const UFIHit *Hit1, const UFIHit *Hit2) const
	{
	if (Hit1 == 0 || Hit2 == 0)
		return UINT_MAX;

	unsigned L1 = m_UD_Fwd.m_Query->m_L;
	unsigned L2 = m_UD_Rev.m_Query->m_L;

	unsigned Pos1 = Hit1->m_DBStartPos;
	unsigned Pos2 = Hit2->m_DBStartPos;

	int iTL;
	if (Pos1 <= Pos2)
		iTL = int(Pos2 + m_UD_Rev.m_Query->m_L) - int(Pos1);
	else
		iTL = int(Pos1 + m_UD_Fwd.m_Query->m_L) - int(Pos2);
	if (iTL < 0 || iTL > 1000)
		iTL = 0;

	return unsigned(iTL);
	}

void State2::SetSAM2()
	{
	m_UD_Fwd.SetMappedPos();
	m_UD_Rev.SetMappedPos();

	const string &TargetLabel1 = m_UD_Fwd.m_MappedTargetLabel;
	const string &TargetLabel2 = m_UD_Rev.m_MappedTargetLabel;

	uint32 TargetPos1 = m_UD_Fwd.m_MappedTargetPos;
	uint32 TargetPos2 = m_UD_Rev.m_MappedTargetPos;

	int TLEN1 = 0;
	int TLEN2 = 0;

	bool Mapped1 = (m_UD_Fwd.m_MappedTargetPos != UINT32_MAX);
	bool Mapped2 = (m_UD_Rev.m_MappedTargetPos != UINT32_MAX);

	const UFIHit *Hit1 = m_UD_Fwd.m_TopHit;
	const UFIHit *Hit2 = m_UD_Rev.m_TopHit;

	bool Plus1 = (Hit1 != 0 && Hit1->m_Plus);
	bool Plus2 = (Hit2 != 0 && Hit2->m_Plus);
	bool StrandsConsistent = (Hit1 != 0 && Hit2 != 0 && (Plus1 != Plus2));

	bool CorrectlyPaired = false;
	if (Mapped1 && Mapped2)
		{
		if (TargetPos1 <= TargetPos2)
			{
			TLEN1 = int(TargetPos2 + m_UD_Rev.m_Query->m_L) - int(TargetPos1);
			if (TLEN1 > 0 && TLEN1 < 1000 && StrandsConsistent)
				CorrectlyPaired = true;
			if (TLEN1 > 1000)
				TLEN1 = 0;
			TLEN2 = -TLEN1;
			}
		else
			{
			TLEN2 = int(TargetPos1 + m_UD_Fwd.m_Query->m_L) - int(TargetPos2);
			if (TLEN2 > 0 && TLEN2 < 1000 && StrandsConsistent)
				CorrectlyPaired = true;
			if (TLEN2 > 1000)
				TLEN2 = 0;
			TLEN1 = -TLEN2;
			}
		}

	bool RevComp1 = (Mapped1 && Hit1 != 0 && !Hit1->m_Plus);
	bool RevComp2 = (Mapped2 && Hit2 != 0 && !Hit2->m_Plus);

	uint32 Flags1 = GetPairedFlags(true, RevComp1, RevComp2, !Mapped2);
	uint32 Flags2 = GetPairedFlags(false, RevComp2, RevComp1, !Mapped1);

	if (CorrectlyPaired)
		{
		Flags1 |= 0x02;
		Flags2 |= 0x02;
		}

	m_UD_Fwd.SetSAM(Flags1, TargetLabel2, TargetPos2, TLEN1);
	m_UD_Rev.SetSAM(Flags2, TargetLabel1, TargetPos1, TLEN2);
	}
