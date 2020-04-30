#include "myutils.h"
#include "ufindex.h"
#include "state1.h"
#include "seqinfo.h"
#include "omplock.h"

void State1::Search(SeqInfo *Query)
	{
	m_Query = Query;
	Lock();
	++m_QueryCount;
	Unlock();

	m_HitCount = 0;
	m_HSPCount = 0;
	m_TopHit = 0;
	m_BestScore = 0;
	m_SecondBestScore = 0;
	m_Mapq = -1;

	Search_Lo();
	SetMappedPos();
//	Output1();
	}

void State1::Search_PE(SeqInfo *Query)
	{
	m_Query = Query;

	m_HitCount = 0;
	m_HSPCount = 0;
	m_TopHit = 0;
	m_BestScore = 0;
	m_SecondBestScore = 0;
	m_Mapq = -1;
	Search_Lo();
	}
