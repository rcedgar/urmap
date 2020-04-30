#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "linereader.h"
#include "samrec.h"
#include "cigar.h"

static unsigned WriteAln(FILE *f, const byte *ReadSeq, unsigned ReadL,
  const byte *RefSeq, unsigned RefSeqL, const string &CIGAR)
	{
	if (f == 0)
		return 0;

	if (CIGAR.empty() || CIGAR == "*")
		Die("CIGAR='%s'", CIGAR.c_str());

	vector<char> Ops;
	vector<unsigned> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);

	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);

	unsigned ReadPos = 0;
	unsigned RefPos = 0;
	string ReadRow;
	string RefRow;
	string AnnotRow;
	unsigned IdCount = 0;
	unsigned SubCount = 0;
	unsigned InsCount = 0;
	unsigned DelCount = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		unsigned n = Lengths[i];
		switch (Op)
			{
		case 'S':
			for (unsigned j = 0; j < n; ++j)
				{
				byte q = toupper(ReadSeq[ReadPos++]);
				ReadRow += q;
				RefRow += ' ';
				AnnotRow += 'S';
				}
			break;

		case 'M':
			for (unsigned j = 0; j < n; ++j)
				{
				byte q = toupper(ReadSeq[ReadPos++]);
				byte t = toupper(RefSeq[RefPos++]);
				ReadRow += q;
				RefRow += t;
				if (q == t)
					{
					++IdCount;
					AnnotRow += '|';
					}
				else
					{
					++SubCount;
					AnnotRow += 'x';
					}
				}
			break;

		case 'D':
			for (unsigned j = 0; j < n; ++j)
				{
				byte t = toupper(RefSeq[RefPos++]);
				++DelCount;
				ReadRow += '-';
				RefRow += t;
				AnnotRow += ' ';
				}
			break;

		case 'I':
			for (unsigned j = 0; j < n; ++j)
				{
				++InsCount;
				byte q = toupper(ReadSeq[ReadPos++]);
				ReadRow += q;
				RefRow +=  '-';
				AnnotRow += ' ';
				}
			break;

		case 'H':
			break;

		default:
			Die("Invalid op '%c' in CIGAR '%s'", Op, CIGAR.c_str());
			}
		}

	fprintf(f, "%s\n", ReadRow.c_str());
	fprintf(f, "%s\n", AnnotRow.c_str());
	fprintf(f, "%s\n", RefRow.c_str());

	return IdCount;
	}

void cmd_sam2aln()
	{
	const string SamFileName = opt(sam2aln);
	const string &DBFileName = opt(ref);
	FILE *fOut = CreateStdioFile(opt(output));

	SeqDB DB;
	DB.FromFasta(DBFileName);

	t_LineBuff LB;
	LineReader LR;
	LR.Open(SamFileName);

	ProgressFile(LR.m_f, "sam2aln", SamFileName);
	for (;;)
		{
		bool Ok = LR.ReadLine(LB);
		if (!Ok)
			break;
		if (LB.Size == 0 || LB.Data[0] == '@')
			continue;

		string Line = string(LB.Data);
		
		SAMRec Rec;
		Rec.FromLine(Line);

		if (Rec.IsUnmapped())
			continue;

		bool Plus = Rec.GetStrand();

		const char *ReadSeq = Rec.m_ReadSeq;
		unsigned ReadL = Rec.m_ReadSeqLength;

		const string &RefLabel = Rec.m_TargetLabel;
		unsigned RefPos = Rec.m_TargetLo;
		asserta(RefPos > 0);
		const string &CIGAR = Rec.m_Cigar;

		unsigned RefSeqIndex = DB.GetSeqIndex(RefLabel);
		const byte *RefSeq = DB.GetSeq(RefSeqIndex) + RefPos - 1;
		unsigned RefSeqL = DB.GetSeqLength(RefSeqIndex);

		int Mapq = Rec.m_Mapq;
		unsigned BitFlags = Rec.m_BitFlags;

		fprintf(fOut, "\n");
		fprintf(fOut, "@%s", Rec.m_ReadLabel);
		if (BitFlags & 0x40)
			fprintf(fOut, " (R1)");
		else if (BitFlags & 0x80)
			fprintf(fOut, " (R2)");
		fprintf(fOut, "\n");
		unsigned IdCount = WriteAln(fOut, (const byte *) ReadSeq, ReadL, RefSeq, RefSeqL, CIGAR);
		double PctId = GetPct(IdCount, ReadL);
		fprintf(fOut, "%s:%u(%c) Q%d, %.1f%% id\n",
		  RefLabel.c_str(), RefPos, pom(Plus), Mapq, PctId);
		}
	ProgressDone();
	LR.Close();
	CloseStdioFile(fOut);
	}
