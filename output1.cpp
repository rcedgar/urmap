#include "myutils.h"
#include "ufindex.h"
#include "state1.h"
#include "seqinfo.h"
#include "omplock.h"
#include "outfiles.h"

void State1::Output1()
	{
	Lock();
	if (g_fSAM != 0)
		{
		SetSAM(0, "*", UINT32_MAX, 0);
		fputs(m_SAM.Data, g_fSAM);
		}
	Unlock();
	UpdateHitStats();
	}

void State1::UpdateHitStats()
	{
	Lock();
	if (m_TopHit == 0)
		++m_NoHitCount;
	else if (m_Mapq >= int(m_Minq))
		++m_AcceptCount;
	else
		++m_RejectCount;
	Unlock();
	}
