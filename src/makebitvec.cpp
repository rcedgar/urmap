#include "myutils.h"
#include "bitvec.h"
#include "alpha.h"
#include "seqdb.h"

static void Scan(BitVec &BV, SeqDB &Ref, uint WordLength, bool Exclude)
	{
	uint64 BitCount = BV.m_BitCount;
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 IncludedCount = 0;
	uint64 ExcludedCount = 0;
	const uint SeqCount = Ref.GetSeqCount();
	ProgressOne(Exclude ? "Scanning exclude" : "Scanning ref");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		Progress_OnTick();
		byte *Seq = Ref.GetSeq(SeqIndex);
		const uint L = Ref.GetSeqLength(SeqIndex);
		for (int Strand = 0; Strand <= 1; ++Strand)
			{
			if (Strand == 1)
				RevCompSeq(Seq, L , Seq);
			uint64 Word = 0;
			byte K = 0;
			for (uint SeqPos = 0; SeqPos < L - WordLength + 1; ++SeqPos)
				{
				byte c = Seq[SeqPos];
				byte Letter = g_CharToLetterNucleo[c];
				if (Letter == INVALID_LETTER)
					{
					K = 0;
					Word = 0;
					continue;
					}
				if (K < WordLength)
					++K;
				Word = ((Word << uint64(2)) | Letter) & ShiftMask;
				if (K == WordLength)
					{
					if (Exclude)
						{
						if (BV.GetBit(Word))
							{
							++ExcludedCount;
							BV.ClearBit(Word);
							}
						}
					else
						{
						if (!BV.GetBit(Word))
							{
							++IncludedCount;
							BV.SetBit(Word);
							}
						}
					}
				}
			}
		}
	ProgressDone();
	if (Exclude)
		ProgressLog("%u words excluded (%s)\n",
		  ExcludedCount, IntToStr(ExcludedCount));
	else
		ProgressLog("%u words included (%s)\n",
		  IncludedCount, IntToStr(IncludedCount));
	}

void cmd_make_bitvec()
	{
	const string &RefFileName = opt(make_bitvec);
	const string &OutputFileName = opt(output);
	const string &Input2FileName = opt(input2);

	asserta(optset_wordlength);
	uint WordLength = opt(wordlength);

	uint64 BitCount = myipow64(4, WordLength);
	ProgressLog("%" PRIu64 " bits (%s)\n", BitCount, Int64ToStr(BitCount));

	FILE *fOut = CreateStdioFile(OutputFileName);

	BitVec BV;
	BV.Alloc(BitCount);

	SeqDB Ref;
	Ref.FromFasta(RefFileName);

	SeqDB Input2;
	Input2.FromFasta(Input2FileName);

	Scan(BV, Ref, WordLength, false);
	Scan(BV, Input2, WordLength, true);

	ProgressFile(fOut, "Writing bitvec", OutputFileName);
	WriteStdioFile(fOut, &BV_MAGIC, sizeof(BV_MAGIC));
	WriteStdioFile(fOut, &WordLength, sizeof(WordLength));
	uint64 Bytes = BV.m_BitCount/8;
	WriteStdioFile64(fOut, BV.m_Vec, Bytes);
	ProgressDone();
	CloseStdioFile(fOut);
	}
