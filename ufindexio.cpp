#include "myutils.h"
#include "ufindex.h"

void UFIndex::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	ProgressOne("Writing %s", FileName.c_str());
	FILE *f = CreateStdioFile(FileName);
	ToFile(f);
	ProgressDone();
	CloseStdioFile(f);
	}

void UFIndex::ToFile(FILE *f) const
	{
	uint32 u;
	u = UFI_MAGIC1;
	WriteStdioFile(f, &u, 4);

	WriteStdioFile(f, &m_WordLength, 4);
	WriteStdioFile(f, &m_MaxIx, 4);
	WriteStdioFile(f, &m_SeqDataSize, 4);
	WriteStdioFile(f, &m_SlotCount, 8);

	const uint32 SeqCount = SIZE(m_Labels);
	WriteStdioFile(f, &SeqCount, 4);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		WriteStdioFile(f, &m_SeqLengths[SeqIndex], 4);
		WriteStdioFile(f, &m_Offsets[SeqIndex], 4);

		const string &Label = m_Labels[SeqIndex];
		u = SIZE(Label);
		WriteStdioFile(f, &u, 4);
		WriteStdioFile(f, Label.c_str(), u);
		}

	u = UFI_MAGIC2;
	WriteStdioFile(f, &u, 4);
	WriteStdioFile64(f, m_Blob, 5*m_SlotCount);

	u = UFI_MAGIC3;
	WriteStdioFile(f, &u, 4);
	WriteStdioFile(f, m_SeqData, m_SeqDataSize);

	u = UFI_MAGIC5;
	WriteStdioFile(f, &u, 4);
	}

void UFIndex::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	ProgressFile(f, "Reading index", FileName);
	FromFile(f);
	ProgressDone();
	CloseStdioFile(f);
	}

void UFIndex::FromFile(FILE *f)
	{
	uint32 u;
	ReadStdioFile(f, &u, 4);
	asserta(u == UFI_MAGIC1);

	ReadStdioFile(f, &m_WordLength, 4);
	ReadStdioFile(f, &m_MaxIx, 4);
	ReadStdioFile(f, &m_SeqDataSize, 4);
	ReadStdioFile(f, &m_SlotCount, 8);
	m_StartPos = 0;
	m_EndPos = m_SeqDataSize - 1;
	SetShiftMask();

	uint32 SeqCount;
	ReadStdioFile(f, &SeqCount, 4);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint32 SeqLength;
		uint32 Offset;
		ReadStdioFile(f, &SeqLength, 4);
		ReadStdioFile(f, &Offset, 4);

		m_SeqLengths.push_back(SeqLength);
		m_Offsets.push_back(Offset);
		m_GenomeSize += SeqLength;

		ReadStdioFile(f, &u, 4);
		char *s = myalloc(char, u+1);
		ReadStdioFile(f, s, u);
		s[u] = 0;

		string Label;
		Label = string(s);
		m_Labels.push_back(Label);
		}

	ReadStdioFile(f, &u, 4);
	asserta(u == UFI_MAGIC2);

	m_Blob = (byte *) malloc(5*m_SlotCount);
	if (m_Blob == 0)
		Die("Out of memory %s:%u", __FILE__, __LINE__);
	ReadStdioFile64(f, m_Blob, 5*m_SlotCount);

	ReadStdioFile(f, &u, 4);
	asserta(u == UFI_MAGIC3);
//	m_SeqData = myalloc(byte, m_SeqDataSize);
	m_SeqData = (byte *) malloc(m_SeqDataSize);
	if (m_SeqData == 0)
		Die("Out of memory %s:%u", __FILE__, __LINE__);
	ReadStdioFile(f, m_SeqData, m_SeqDataSize);

	ReadStdioFile(f, &u, 4);
	asserta(u == UFI_MAGIC5);
	}

void cmd_make_ufi()
	{
	const string &RefFileName = opt(make_ufi);
	double LoadFactor = opt(load_factor);
	unsigned WordLength = 24;
	unsigned MaxIx = 32;
	if (!optset_notrunclabels)
		{
		opt_trunclabels = true;
		optset_trunclabels = true;
		optused_trunclabels = true;
		}

	if (optset_wordlength)
		WordLength = opt(wordlength);

	if (opt(veryfast))
		MaxIx = 3;
	if (optset_maxix)
		MaxIx = opt(maxix);

	FILE *f = OpenStdioFile(RefFileName);
	int64 GenomeSize = GetStdioFileSize64(f);
	CloseStdioFile(f);

	uint64 SlotCount;
	if (optset_slots)
		SlotCount = StrToUint64(opt(slots));
	else
		{
		int64 iMinSlotCount = int64(GenomeSize/LoadFactor);
		uint64 GetPrime(uint64 n);
		SlotCount = GetPrime(iMinSlotCount);
		}

	if (GenomeSize > UINT32_MAX - 100000)
		Die("Genome too big (%s), email robert@drive5.com for update",
		  Int64ToStr(GenomeSize));

	Progress("\n");
	Progress("  Genome size  %s\n", Int64ToStr(GenomeSize));
	Progress("        Slots  %" PRIu64 " (%s)\n", SlotCount, Int64ToStr(SlotCount));
	Progress("  Load factor  %.2f\n", GenomeSize/double(SlotCount));
	Progress("  Word length  %u\n", WordLength);
	Progress("   Max abund.  %u\n", MaxIx);
	Progress("\n");

	UFIndex UFI;
	UFI.m_WordLength = WordLength;
	UFI.m_MaxIx = MaxIx;
	UFI.m_SlotCount = SlotCount;
	UFI.SetShiftMask();
	UFI.ReadSeqData(RefFileName);
	UFI.MakeIndex();
	if (opt(validate))
		UFI.Validate();
	if (opt(logstats))
		{
		UFI.LogMe();
		UFI.LogStats();
		}
	UFI.ToFile(opt(output));
	}

void cmd_suggest_slots()
	{
	double GenomeSize = StrToFloat(opt(suggest_slots))*1e6;
	double LoadFactor = opt(load_factor);

	Progress("\n");
	Progress("Genome size %s\n", FloatToStr(GenomeSize));
	Progress("Load factor %.2f\n", LoadFactor);

	int64 iMinSlotCount = int64(GenomeSize/LoadFactor);

	uint64 GetPrime(uint64 n);
	uint64 SlotCount = GetPrime(iMinSlotCount);

	Progress("Prime slots %" PRIu64 " (%s)\n", SlotCount, Int64ToStr(SlotCount));
	Progress("Prime load  %.2f\n", GenomeSize/double(SlotCount));
	Progress("\n");
	}
