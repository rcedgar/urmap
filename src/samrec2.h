#ifndef samrec2_h
#define samrec2_h

#include "seqinfo.h"

enum SAM_FIELDS
	{
	F_QNAME = 0,
	F_FLAG = 1,
	F_RNAME = 2,
	F_POS = 3,
	F_MAPQ = 4,
	F_CIGAR = 5,
	F_RNEXT = 6,
	F_PNEXT = 7,
	F_TLEN = 8,
	F_SEQ = 9,
	F_QUAL = 10
	};

enum SAM_BITFLAGS
	{
	B_PAIRED = 0x1,				//	1		template having multiple segments in sequencing
	B_PROPER = 0x2,				//	2		each segment properly aligned according to the aligner
	B_UNMAPPED = 0x4,			//	4		segment unmapped
	B_MATE_UNMAPPED = 0x8,		//	8		next segment in the template unmapped
	B_SEQ_REVCOMP = 0x10,		//	16		SEQ being reverse complemented
	B_SEQ_MATE_REVCOMP = 0x20,	//	32		SEQ of the next segment in the template being reversed
	B_R1 = 0x40,				//	64		the first segment in the template
	B_R2 = 0x80,				//	128		the last segment in the template
	B_SECONDARY = 0x100,		//	256		secondary alignment (BWA uses 0x800 instead)
	B_QCFAIL = 0x200,			//	512		not passing quality control
	B_DUPE = 0x400,				//	1024	PCR or optical duplicate
	B_SECONDARY2 = 0x800		//	2048	secondary alignment (BWA only, not in SAM spec).
	};

#define Fld(x)	m_Fields[F_##x]

class SAMRec2
	{
public:
	vector<string> m_Fields;
	unsigned m_Flags;
	unsigned m_Pos;
	unsigned m_Mapq;

public:
	void FromLine(const string &Line);

	const string &GetReadLabel() const
		{
		return Fld(QNAME);
		}

	const string &GetTargetLabel() const
		{
		return Fld(RNAME);
		}

	bool IsRevComp() const
		{
		return (m_Flags & B_SEQ_REVCOMP) != 0;
		}

	bool IsSec() const
		{
		return (m_Flags & B_SECONDARY) != 0;
		}

	void ToSeqInfo(SeqInfo *SI) const;
	void ToFile(FILE *f) const;
	};

#endif // samrec2_h
