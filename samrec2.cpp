#include "myutils.h"
#include "seqinfo.h"
#include "samrec2.h"

void SAMRec2::FromLine(const string &Line)
	{
	Split(Line, m_Fields, '\t');
	asserta(SIZE(m_Fields) >= 11);
	m_Flags = StrToUint(Fld(FLAG));
	m_Pos = StrToUint(Fld(POS));
	m_Mapq  = StrToUint(Fld(MAPQ));
	}

void SAMRec2::ToSeqInfo(SeqInfo *SI) const
	{
	const string &Label = Fld(QNAME);
	const string &Seq = Fld(SEQ);
	const string &Qual = Fld(QUAL);

	unsigned LabelBytes = SIZE(Label) + 1;
	unsigned L = SIZE(Seq);

	SI->AllocSeq(L);
	SI->AllocQual(L);
	SI->AllocLabel(LabelBytes);

	memcpy(SI->m_SeqBuffer, Seq.c_str(), L);
	memcpy(SI->m_QualBuffer, Qual.c_str(), L);
	memcpy(SI->m_LabelBuffer, Label.c_str(), LabelBytes);

	if (IsRevComp())
		RevCompSeq(SI->m_SeqBuffer, L, SI->m_SeqBuffer);

	SI->m_Seq = SI->m_SeqBuffer;
	SI->m_Label = SI->m_LabelBuffer;
	SI->m_L = L;
	}

void SAMRec2::ToFile(FILE *f) const
	{
	if (f == 0)
		return;
	const unsigned N = SIZE(m_Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		if (i > 0)
			fputc('\t', f);
		fputs(m_Fields[i].c_str(), f);
		}
	fputc('\n', f);
	}
