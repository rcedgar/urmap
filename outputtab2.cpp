#include "myutils.h"
#include "state2.h"
#include "outfiles.h"
#include "omplock.h"

void State2::GetInfoStr(string &s) const
	{
	s.clear();

	if (GetTopHit1() == 0 || GetTopHit2() == 0 || GetSecondHit1() == 0 || GetSecondHit2() == 0)
		return;

	unsigned TopTL = GetTemplateLength(GetTopHit1(), GetTopHit2());
	unsigned SecondTL = GetTemplateLength(GetSecondHit1(), GetSecondHit2());

	if (TopTL == SecondTL)
		Psasc(s, "TL=%u", TopTL, SecondTL);
	else
		Psasc(s, "TL/%u,%u", TopTL, SecondTL);

	int TopScore = GetTopHit1()->m_Score + GetTopHit2()->m_Score;
	int SecondScore = GetSecondHit1()->m_Score + GetSecondHit2()->m_Score;
	if (TopScore == SecondScore)
		Psasc(s, "Score=%d", TopScore);
	else
		Psasc(s, "Score/%d,%d", TopScore, SecondScore);
	}

void State2::GetPairPosStr1(const UFIHit *Hit, bool Fwd, string &s) const
	{
	asserta(Hit != 0);

	string ChrLabel;
	uint32 ChrPos;
	GetChrPos(Hit, ChrLabel, ChrPos);

	char OneOrTwo = (Fwd ? '1' : '2');
	char Strand = pom(Hit->m_Plus);

	Ps(s, "%s:%u(%c)/%c", ChrLabel.c_str(), ChrPos+1, Strand, OneOrTwo);
	}

void State2::GetPairPosStr(const UFIHit *Hit1, const UFIHit *Hit2, string &s) const
	{
	if (Hit1 == 0 && Hit2 == 0)
		{
		s = "*";
		return;
		}

	if (Hit1 != 0 && Hit2 == 0)
		{
		GetPairPosStr1(Hit1, true, s);
		return;
		}

	if (Hit1 == 0 && Hit2 != 0)
		{
		GetPairPosStr1(Hit2, false, s);
		return;
		}

	string ChrLabel1;
	string ChrLabel2;

	uint32 ChrPos1;
	uint32 ChrPos2;

	GetChrPos(Hit1, ChrLabel1, ChrPos1);
	GetChrPos(Hit2, ChrLabel2, ChrPos2);
	if (ChrLabel1 == ChrLabel2 && Hit1->m_Plus != Hit2->m_Plus)
		{
		char Strand = pom(Hit1->m_Plus);
		Ps(s, "%s:%u-%u", ChrLabel1.c_str(), ChrPos1+1, ChrPos2+1);
		return;
		}

	string s1;
	string s2;
	GetPairPosStr1(Hit1, true, s1);
	GetPairPosStr1(Hit2, false, s2);
	s = s1 + "," + s2;
	}

void State2::OutputTab2()
	{
	if (g_fTab == 0)
		return;

	string PairLabel;
	m_UD_Fwd.GetPairLabel(PairLabel);

	string sTop;
	GetPairPosStr(GetTopHit1(), GetTopHit2(), sTop);

	string sSecond = "*";
	string InfoStr;
	if (GetSecondHit1() != 0)
		{
		asserta(GetSecondHit2() != 0);
		GetPairPosStr(GetSecondHit1(), GetSecondHit2(), sSecond);
		}
	GetInfoStr(InfoStr);

	Lock();
	fputs(PairLabel.c_str(), g_fTab);
	fputc('\t', g_fTab);
	fputs(sTop.c_str(), g_fTab);
	fputc('\t', g_fTab);
	fprintf(g_fTab, "%u,%u", m_UD_Fwd.m_Mapq, m_UD_Rev.m_Mapq);
	fputc('\t', g_fTab);
	fputs(sSecond.c_str(), g_fTab);
	if (!InfoStr.empty())
		{
		fputc('\t', g_fTab);
		fputs(InfoStr.c_str(), g_fTab);
		}
	fputc('\n', g_fTab);
	Unlock();
	}
