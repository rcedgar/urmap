#include "myutils.h"
#include "ufindex.h"
#include "seqinfo.h"
#include "alpha.h"
#include "state1.h"
#include "ufihit.h"
#include "alignresult.h"
#include "sort.h"
#include "omplock.h"
#include "cigar.h"

void State1::SetSAM_Unmapped(uint32 aFlags)
	{
	uint32 Flags = 0x04; // unmapped

	if (aFlags & 0x01)	// multiseg
		Flags |= 0x01;

	if (aFlags & 0x40)	// first
		Flags |= 0x40;
	else if (aFlags & 0x80)	// second
		Flags |= 0x80;

	if (aFlags & 0x08)	// mate unmapped
		Flags |= 0x08;
	else if (aFlags & 0x20)	// mate revcomp'd
		Flags |= 0x20;

	m_Mapq = 0;

	unsigned QL = m_Query->m_L;
	m_SAM.Alloc(3*max(QL, 200u));
	char *p = m_SAM.Data;
	const char *Label = m_Query->m_Label;
	unsigned n = ustrlen(Label);
	if (n > 2 && Label[n-2] == '/' && (Label[n-1] == '1' || Label[n-1] == '2'))
		n -= 2;
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Label[i];
		if (c == ' ' || c == '\t')
			break;
		*p++ = c;
		}
	*p++ = '\t';

	char Tmp[16];
	sprintf(Tmp, "%u", Flags);
	n = ustrlen(Tmp);
	for (unsigned i = 0; i < n; ++i)
		*p++ = Tmp[i];

	const char s[] = "	*	0	0	*	*	0	0	";
	n = sizeof(s) - 1;
	memcpy(p, s, n);
	p += n;

	memcpy(p, m_Query->m_Seq, QL);
	p += QL;
	*p++ = '\t';

	const char *Qual = m_Query->m_Qual;
	if (Qual == 0)
		*p++ = '*';
	else
		{
		memcpy(p, Qual, QL);
		p += QL;
		}
	*p++ = '\n';
	*p = 0;
	m_SAM.Size = unsigned(p - m_SAM.Data);
	}

void State1::SetSAM(uint32 Flags, const string &MateTargetLabel,
  const uint32 MateTargetPos, int TLEN)
	{
	if (m_TopHit == 0)
		{
		SetSAM_Unmapped(Flags);
		return;
		}

	unsigned QL = m_Query->m_L;
	m_SAM.Alloc(3*max(QL, 200u));
	char *p = m_SAM.Data;
	const char *Label = m_Query->m_Label;
	unsigned n = ustrlen(Label);
	if (n > 2 && Label[n-2] == '/' && (Label[n-1] == '1' || Label[n-1] == '2'))
		n -= 2;
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Label[i];
		if (c == ' ' || c == '\t')
			break;
		*p++ = c;
		}
	*p++ = '\t';

	bool Plus = m_TopHit->m_Plus;
	char Tmp[16];
	sprintf(Tmp, "%u", Flags);
	n = ustrlen(Tmp);
	for (unsigned i = 0; i < n; ++i)
		*p++ = Tmp[i];
	*p++ = '\t';

	//string TargetLabel;
	//unsigned TargetLo;
	//m_TopHit->GetStartCoord(TargetLabel, TargetLo);
	if (m_MappedTargetPos == UINT32_MAX)
		{
		SetSAM_Unmapped(Flags);
		return;
		}

	n = SIZE(m_MappedTargetLabel);
	memcpy(p, m_MappedTargetLabel.c_str(), n);
	p += n;
	*p++ = '\t';

	sprintf(Tmp, "%u", m_MappedTargetPos+1);
	n = ustrlen(Tmp);
	memcpy(p, Tmp, n);
	p += n;
	*p++ = '\t';

	sprintf(Tmp, "%u", m_Mapq);
	n = ustrlen(Tmp);
	memcpy(p, Tmp, n);
	p += n;
	*p++ = '\t';

	if (m_TopHit->m_Path.empty())
		{
		sprintf(Tmp, "%d", QL);
		n = ustrlen(Tmp);
		memcpy(p, Tmp, n);
		p += n;
		*p++ = 'M';
		}
	else
		{
		string CIGAR;
		GetCIGAR(CIGAR);
		n = ustrlen(CIGAR);
		memcpy(p, CIGAR.c_str(), n);
		p += n;
		}
	*p++ = '\t';

	if (MateTargetLabel == "" || MateTargetLabel == "*")
		*p++ = '*';
	else if (MateTargetLabel == m_MappedTargetLabel)
		*p++ = '=';
	else
		{
		n = ustrlen(MateTargetLabel);
		memcpy(p, MateTargetLabel.c_str(), n);
		p += n;
		}
	*p++ = '\t';

	if (MateTargetPos == 0 || MateTargetPos == UINT32_MAX)
		*p++ = '0';
	else
		{
		sprintf(Tmp, "%u", MateTargetPos+1);
		n = ustrlen(Tmp);
		memcpy(p, Tmp, n);
		p += n;
		}
	*p++ = '\t';

	sprintf(Tmp, "%d", TLEN);
	n = ustrlen(Tmp);
	memcpy(p, Tmp, n);
	p += n;
	*p++ = '\t';

	const byte *Seq = (Plus ? m_Query->m_Seq : m_RevCompQuerySeq.Data);
	memcpy(p, Seq, QL);
	p += QL;
	*p++ = '\t';

	const char *Qual = m_Query->m_Qual;
	if (Qual ==0)
		*p++ = '*';
	else
		{
		if (Plus)
			{
			memcpy(p, Qual, QL);
			p += QL;
			}
		else
			{
			for (unsigned i = 1; i <= QL; ++i)
				*p++ = Qual[QL-i];
			}
		}

	*p++ = '\n';
	*p = 0;

	m_SAM.Size = unsigned(p - m_SAM.Data);
	}
