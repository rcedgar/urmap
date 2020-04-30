#include "myutils.h"
#include "ufindex.h"
#include "alpha.h"
#include <inttypes.h>

void UFIndex::LogSeqDict() const
	{
	unsigned n = SIZE(m_Labels);
	Log("%u seqs\n", n);
	for (unsigned i = 0; i < n; ++i)
		{
		Log("[%3u]", i);
		Log("  %10u", m_SeqLengths[i]);
		Log("  %10u", m_Offsets[i]);
		Log("  %s", m_Labels[i].c_str());
		Log("\n");
		}
	}

void UFIndex::SetShiftMask()
	{
	m_ShiftMask = 0;
	for (unsigned i = 0; i < 2u*m_WordLength; ++i)
		m_ShiftMask |= (uint64(1) << i);

	for (unsigned i = 0; i < m_WordLength; ++i)
		m_LowerCaseShiftMask |= (uint64(1) << i);
	}

const char *UFIndex::GetStr(uint32 Pos, unsigned n, string &s) const
	{
	s.clear();
	for (unsigned i = 0; i < n; ++i)
		{
		byte c = m_SeqData[Pos+i];
		s += c;
		}
	return s.c_str();
	}

const char *UFIndex::WordToStr(uint64 Word, string &s) const
	{
	s.clear();
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte Letter = Word & 3;
		byte c = g_LetterToCharNucleo[Letter];
		s = char(c) + s;
		Word >>= 2;
		}
	return s.c_str();
	}

uint64 UFIndex::SeqToWord(const byte *Seq) const
	{
	uint64 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte c = Seq[i];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter > 3)
			return UINT64_MAX;
		Word = (Word << uint64(2)) | Letter;
		}
	return Word;
	}

uint64 UFIndex::GetWord(uint32 Pos) const
	{
	asserta(Pos < m_SeqDataSize);
	uint64 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte c = m_SeqData[Pos+i];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter > 3)
			return UINT64_MAX;
		Word = (Word << uint64(2)) | Letter;
		}
	return Word;
	}

void UFIndex::MakeIndex()
	{
	if (m_StartPos == UINT_MAX)
		{
		m_StartPos = 0;
		m_EndPos = m_SeqDataSize;
		}

	m_TruncatedCount = 0;
	ProgressOne("Init vectors");
	uint64 Size = 5*m_SlotCount;
	m_Blob = myalloc64(byte, Size);
	memset_zero(m_Blob, Size);
	for (uint64 Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		SetPos(Slot, UINT32_MAX);
		SetTally(Slot, TALLY_FREE);
		}
	ProgressDone();

	CountSlots();
	CountSlots_Minus();
//	LogCountHist();

	uint64 Word = 0;
	byte K = 0;
	uint32 SeqPos;
	ProgressOne("Index");
	for (SeqPos = m_StartPos; SeqPos < m_EndPos; ++SeqPos)
		{
		ProgressLoopTick(SeqPos);
		byte c = m_SeqData[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint32 StartPos = SeqPos - (m_WordLength-1);
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
#if TSLOT
			if (Slot == TSLOT)
				{
				uint64 Word2 = GetWord(StartPos);
				asserta(Word2 == (Word & m_ShiftMask));
				string s;
				GetStr(StartPos, m_WordLength, s);
				Log("\n");
				Log("============================================\n");
				Log("Update\n");
				Log("Slot   %" PRIx64 "\n", Slot);
				Log("Word   %" PRIx64 "\n", Word & m_ShiftMask);
				Log("Word2  %" PRIx64 "\n", Word2);
				Log("Seq    %s\n", s.c_str());
				LogList(Slot);
				}
#endif
			UpdateSlot(Slot, StartPos);
			}
		}
	ProgressDone();
	Progress("%u slots truncated\n", m_TruncatedCount);
	}

void UFIndex::TruncateSlot(uint64 Slot)
	{
	++m_TruncatedCount;

	byte T  = GetTally(Slot);
	asserta(!TallyOther(T));

	uint32 Pos = GetPos(Slot);
	uint64 Slot2 = Slot;
	unsigned K = 0;
	for (;;)
		{
		byte T  = GetTally(Slot2);
		uint32 Pos = GetPos(Slot2);
		SetTally(Slot2, TALLY_FREE);
		SetPos(Slot2, UINT32_MAX);
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			assert(Slot2 == Slot);
			return;
			}

		if (T == TALLY_END)
			return;
		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	Die("TruncateSlot(%" PRIx64 ")", Slot);
	}

void UFIndex::UpdateSlot(uint64 Slot, uint32 Pos)
	{
	byte n = m_SlotCounts[Slot];
	byte n_minus = m_SlotCounts_Minus[Slot];
#if TSLOT
	if (Slot == TSLOT)
		{
		Log("\n");
		Log("_______________________________________\n");
		Log("\nUpdateSlot(%" PRIx64 ", %u) nplus %u, nminus %u\n", Slot, Pos, n, n_minus);
		Log("Before update... ");
		LogList(Slot);
		}
#endif
	if (n > m_MaxIx || n_minus > m_MaxIx)
		{
#if TSLOT
		if (Slot == TSLOT)
			Log("n=%u, n_minus=%u > m_MaxIx=%u\n", n, n_minus, m_MaxIx);
#endif
		return;
		}
	byte T = GetTally(Slot);
	if (T == TALLY_FREE)
		{
		assert(GetPos(Slot) == UINT32_MAX);
		SetPos(Slot, Pos);
		if (n == 1 && n_minus == 0)
			SetTally(Slot, TALLY_BOTH1);
		else
			SetTally(Slot, TALLY_PLUS1);
#if TSLOT
		if (Slot == TSLOT)
			{
			Log("TALLY_FREE (n=%u): ", n);
			LogRow(Slot);
			Log("\n");
			}
#endif
		return;
		}
	uint64 EOLSlot = FindEndOfList(Slot);
#if TSLOT
	if (Slot == TSLOT)
		Log("EOLSlot=%" PRIx64 "\n", EOLSlot);
#endif
	unsigned Step = FindFreeSlot(EOLSlot);
#if TSLOT
	if (Slot == TSLOT)
		Log("Step=%u\n", Step);
#endif
	if (Step == UINT_MAX)
		{
		TruncateSlot(Slot);
		return;
		}

	uint64 FreeSlot = (EOLSlot + Step)%m_SlotCount;
#if TSLOT
	if (Slot == TSLOT)
		Log("FreeSlot=%" PRIx64 "\n", FreeSlot);
#endif
	if (Step > TALLY_MAX_NEXT)
		{
		unsigned Step2 = FindFreeSlot(FreeSlot);
#if TSLOT
		if (Slot == TSLOT)
			Log("Step2=%u\n", Step2);
#endif
		if (Step2 == UINT_MAX)
			{
			TruncateSlot(Slot);
			return;
			}
		uint32 EOLPos = GetPos(EOLSlot);
#if TSLOT
		if (Slot == TSLOT)
			Log("EOLPos=%" PRIu64 "\n", EOLPos);
#endif

		uint64 FreeSlot2 = (FreeSlot + Step2)%m_SlotCount;
#if TSLOT
		if (Slot == TSLOT)
			Log("FreeSlot2=%" PRIx64 "\n", FreeSlot2);
#endif

		if (EOLSlot == Slot)
			SetTallyNext(EOLSlot, TALLY_NEXT_LONG_MINE);
		else
			SetTallyNext(EOLSlot, TALLY_NEXT_LONG_OTHER);
		SetPos(EOLSlot, Step | (Step2 << 16));

		SetTally(FreeSlot, TALLY_NEXT_LONG_OTHER);
		SetPos(FreeSlot, EOLPos);

		SetTally(FreeSlot2, TALLY_END);
		SetPos(FreeSlot2, Pos);

#if TSLOT
		if (Slot == TSLOT)
			{
			Log("After update long...");
			LogList(Slot);
			}
#endif
		return;
		}

#if TSLOT
	if (Slot == TSLOT)
		Log("FreeSlot=%" PRIx64 "\n");
#endif

// EOL slot
	byte Next = byte(Step);
	SetTallyNext(EOLSlot, Next);

	assert(GetTally(FreeSlot) == TALLY_FREE);
	SetTally(FreeSlot, TALLY_END);
	SetPos(FreeSlot, Pos);
#if TSLOT
	if (Slot == TSLOT)
		{
		Log(" After update, row: ");
		LogRow(Slot);
		Log("\n");
		}
#endif
	}

void UFIndex::FreeTempBuildData()
	{
	if (m_SlotCounts != 0)
		{
		myfree(m_SlotCounts);
		m_SlotCounts = 0;
		}
	if (m_SlotCounts_Minus != 0)
		{
		myfree(m_SlotCounts_Minus);
		m_SlotCounts_Minus = 0;
		}
	}

void UFIndex::CountSlots()
	{
	ProgressOne("Init slot counts");
	m_SlotCounts = myalloc64(byte, m_SlotCount);
	memset_zero(m_SlotCounts, m_SlotCount);
	ProgressDone();

	uint64 Word = 0;
	byte K = 0;
	uint32 SeqPos;
	ProgressOne("Count slots");
	for (SeqPos = m_StartPos; SeqPos < m_EndPos; ++SeqPos)
		{
		ProgressLoopTick(SeqPos);
		byte c = m_SeqData[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			if (m_SlotCounts[Slot] < 255)
				++(m_SlotCounts[Slot]);
			}
		}
	ProgressDone();
	}

void UFIndex::CountSlots_Minus()
	{
	ProgressOne("Init minus slot counts");
	m_SlotCounts_Minus = myalloc64(byte, m_SlotCount);
	memset_zero(m_SlotCounts_Minus, m_SlotCount);
	ProgressDone();

	uint64 Word = 0;
	byte K = 0;
	uint32 nnn;
	uint32 SeqPos = m_EndPos;
	ProgressOne("Count slots minus");
	for (nnn = m_StartPos; nnn < m_EndPos; ++nnn)
		{
		ProgressLoopTick(nnn);
		byte c = m_SeqData[--SeqPos];
		byte Letter = g_CharToCompLetter[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			if (m_SlotCounts_Minus[Slot] < 255)
				++(m_SlotCounts_Minus[Slot]);
			}
		}
	asserta(SeqPos == m_StartPos);
	ProgressDone();
	}

void UFIndex::CountIndexedWords(unsigned &IndexedCount,
  unsigned &NotIndexedCount, unsigned &WildcardCount)
	{
	IndexedCount = 0;
	NotIndexedCount = 0;
	WildcardCount = 0;

	uint64 Word = 0;
	byte K = 0;
	uint32 *PosVec = myalloc(uint32, m_MaxIx);
	ProgressOne("Count indexed words");
	for (uint32 SeqPos = m_StartPos; SeqPos < m_EndPos; ++SeqPos)
		{
		byte c = m_SeqData[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			++WildcardCount;
			continue;
			}
		if (K < m_WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == m_WordLength)
			{
			uint64 Slot = WordToSlot(Word & m_ShiftMask, m_SlotCount);
			unsigned RowLen = GetRow(Slot, PosVec);
			bool Found = false;
			uint32 StartPos = SeqPos - (m_WordLength-1);
			for (unsigned k = 0; k < RowLen; ++k)
				{
				uint32 Pos_k = PosVec[k];
				if (Pos_k == StartPos)
					{
					Found = true;
					break;
					}
				}
			if (Found)
				++IndexedCount;
			else
				++NotIndexedCount;
			}
		else
			++WildcardCount;
		}
	ProgressDone();
	myfree(PosVec);
	}

void UFIndex::ReadSeqData(const string &FastaFileName)
	{
	SeqDB DB;
	DB.FromFasta(FastaFileName);
	DB.ToUpper();
	const unsigned SeqCount = DB.GetSeqCount();
	ProgressOne("Initialize seqs");
	m_SeqDataSize = 0;
	bool AltWarningDone = false;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Label = DB.GetLabel(SeqIndex);
		if (EndsWith(string(Label), "_alt"))
			{
			Warning("\nURMAP is not alt-aware, recommend removing alt chromosomes\n");
			AltWarningDone = true;
			}

		unsigned L = DB.GetSeqLength(SeqIndex);

		m_Labels.push_back(string(Label));
		m_SeqLengths.push_back(L);
		m_Offsets.push_back(m_SeqDataSize);

		m_SeqDataSize += L;
		m_GenomeSize += L;
		if (SeqIndex + 1 != SeqCount)
			m_SeqDataSize += PADGAP;
		}
	m_SeqData = myalloc(byte, m_SeqDataSize);

	unsigned Offset = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Label = DB.GetLabel(SeqIndex);
		const byte *Seq = DB.GetSeq(SeqIndex);
		unsigned L = DB.GetSeqLength(SeqIndex);
		memcpy(m_SeqData + Offset, Seq, L);
		Offset += L;
		if (SeqIndex + 1 != SeqCount)
			{
			for (unsigned i = 0; i < PADGAP; ++i)
				m_SeqData[Offset+i] = '-';
			Offset += PADGAP;
			}
		}
	asserta(Offset == m_SeqDataSize);
	DB.Free();
	ProgressDone();
	}

void UFIndex::LogMe() const
	{
	Log("\n");
	Log("Word %u\n", m_WordLength);
	Log("m_MaxIx %u\n", m_MaxIx);
	LogSlots(0, 100);
	Log("\n\n");
	LogSlots(1000, 2000);
	Log("\n\n");
	LogSlots(m_SlotCount - 101, 100);
	}

void UFIndex::LogSlots(uint64 SlotLo, uint64 N) const
	{
	for (uint64 i = 0; i < N; ++i)
		{
		uint64 Slot = SlotLo + i;
		byte T = GetTally(Slot);
		bool Mine = TallyMine(T);
		uint32 Pos = GetPos(Slot);
		Log("%10" PRIx64, Slot);
		if (m_SlotCounts != 0)
			Log("  %3u", m_SlotCounts[Slot]);
		Log("  %02x", T);
		if (Pos == UINT32_MAX)
			Log("  ..........");
		else if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			Log("  steps %u, %u", (Pos & 0xffff), Pos >> 16);
		else
			Log("  %10u", Pos);
		Log("  ");
		LogTally(T);
		if (TallyMine(T))
			LogRow(Slot);
		Log("\n");
		}
	}

void UFIndex::GetCollisionCount(uint64 &CollisionCount, uint64 &IndexedCount) const
	{
	CollisionCount = 0;
	IndexedCount = 0;
	string s;
	uint32 *PosVec = myalloc(uint32, m_MaxIx);
	uint64 Slot = 0;
	ProgressLoop64(m_SlotCount, &Slot, "Collisions");
	for (Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		unsigned K = GetRow(Slot, PosVec);
		asserta(K < 256);
		IndexedCount += K;
		const byte *Str0 = m_SeqData + PosVec[0];
		for (unsigned k = 1; k < K; ++k)
			{
			uint32 Pos = PosVec[k];
			asserta(Pos < m_SeqDataSize);
			const byte *Str = m_SeqData + PosVec[k];
			if (memcmp(Str, Str0, m_WordLength) != 0)
				++CollisionCount;
			}
		}
	ProgressDone();
	myfree(PosVec);
	}

void UFIndex::LogSlot(uint64 Slot) const
	{
	Log("\n");
	Log("LogSlot(%" PRIx64 ")\n", Slot);
	uint32 *PosVec = myalloc(uint32, m_MaxIx);
	unsigned K = GetRow(Slot, PosVec);
	Log(" K=%u\n", K);
	asserta(K <= m_MaxIx);
	for (unsigned k = 0; k < K; ++k)
		{
		uint32 Pos = PosVec[k];
		asserta(Pos < m_SeqDataSize);

		Log("\n");
		Log(" k          %u\n", k);
		Log(" Pos        %u\n", Pos);

		uint64 Word = GetWord(Pos);
		uint64 Slot2 = WordToSlot(Word, m_SlotCount);

		string s;
		GetStr(Pos, 24, s);
		uint64 Word2 = SeqToWord((const byte *) s.c_str());

		Log(" Str        %s\n", s.c_str());
		Log(" Word       %" PRIx64 "\n", Word);
		Log(" Word2      %" PRIx64 "\n", Word2);
		Log(" Slot       %" PRIx64 "\n", Slot);
		Log(" Slot2      %" PRIx64 "\n", Slot2);
		}
	myfree(PosVec);
	}

void UFIndex::ValidateSlot(uint64 Slot) const
	{
	uint32 *PosVec = myalloc(uint32, m_MaxIx);
	unsigned K = GetRow_Validate(Slot, PosVec);
	asserta(K <= m_MaxIx);
	for (unsigned k = 0; k < K; ++k)
		{
		uint32 Pos = PosVec[k];
		asserta(Pos < m_SeqDataSize);
		uint64 Word = GetWord(Pos);
		uint64 Slot2 = WordToSlot(Word, m_SlotCount);
		if (Slot2 != Slot)
			{
			Log("\n");
			Log("============== validate fail 0x%" PRIx64 "===============\n", Slot);
			if (Slot >= 2)
				LogSlots(Slot-2, 256);
			string s;
			GetStr(Pos, 24, s);
			uint64 Word2 = SeqToWord((const byte *) s.c_str());

			Log("\n");
			Log(" k          %u\n", k);
			Log(" K          %u\n", K);
			Log(" Pos        %u\n", Pos);
			Log(" Str        %s\n", s.c_str());
			Log(" Word       %" PRIx64 "\n", Word);
			Log(" Word2      %" PRIx64 "\n", Word2);
			Log(" Slot       %" PRIx64 "\n", Slot);
			Log(" Slot2      %" PRIx64 "\n", Slot2);
			Log(" SlotCount  %" PRIx64 "\n", m_SlotCount);
			Die("WordToSlot != Slot");
			}
		}
	myfree(PosVec);
	}

void UFIndex::Validate() const
	{
	uint64 Slot;
	ProgressLoop64(m_SlotCount, &Slot, "Validate");
	for (Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		ProgressLoopTick(Slot);
		ValidateSlot(Slot);
		}
	ProgressDone();
	}

const byte *UFIndex::GetSeq(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_Offsets));
	unsigned Offset = m_Offsets[SeqIndex];
	return m_SeqData + Offset;
	}

unsigned UFIndex::GetSeqIndex(const string &Label) const
	{
	for (unsigned i = 0; i < SIZE(m_Labels); ++i)
		if (m_Labels[i] == Label)
			return i;
	Die("Label not found >%s", Label.c_str());
	return UINT_MAX;
	}

uint32 UFIndex::CoordToPos(const string &Label, uint32 Coord) const
	{
	unsigned SeqIndex = GetSeqIndex(Label);
	uint32 Offset = m_Offsets[SeqIndex];
	uint32 Pos = Offset + Coord;
	return Pos;
	}

void UFIndex::GetUserPosStr(uint32 Pos, string &s) const
	{
	string Label;
	unsigned Coord = PosToCoord(Pos, Label);
	if (Coord == UINT32_MAX)
		{
		s = "*";
		return;
		}

	if (Label[0] == 'c' && Label[1] == 'h' && Label[2] == 'r')
		s = Label.substr(3, string::npos);
	else
		s = Label;
	Psa(s, ":%u", Coord + 1);
	}

uint32 UFIndex::PosToCoord(uint32 Pos, string &Label) const
	{
	const unsigned SeqCount = SIZE(m_Labels);

	unsigned Lo = 0;
	unsigned Hi = SeqCount - 1;
	while (Lo <= Hi)
		{
		unsigned SeqIndex = (Lo + Hi)/2;
		uint32 Offset = m_Offsets[SeqIndex];
		uint32 SeqLength = m_SeqLengths[SeqIndex];
		if (Pos >= Offset && Pos < Offset + SeqLength)
			{
			Label = m_Labels[SeqIndex];
			uint32 Coord = Pos - Offset;
			return Coord;
			}
		else if (Pos > Offset)
			Lo = SeqIndex + 1;
		else
			Hi = SeqIndex - 1;
		}

// Anomaly -- this point reached (rarely)
// when alignment extends into inter-chr padding
	return UINT32_MAX;
	}

uint32 UFIndex::PosToCoordL(uint32 Pos, string &Label, unsigned &L) const
	{
	const unsigned SeqCount = SIZE(m_Labels);
	unsigned Lo = 0;
	unsigned Hi = SeqCount - 1;
	while (Lo <= Hi)
		{
		unsigned SeqIndex = (Lo + Hi)/2;
		uint32 Offset = m_Offsets[SeqIndex];
		uint32 SeqLength = m_SeqLengths[SeqIndex];
		if (Pos >= Offset && Pos < Offset + SeqLength)
			{
			Label = m_Labels[SeqIndex];
			uint32 Coord = Pos - Offset;
			L = m_SeqLengths[SeqIndex];
			return Coord;
			}
		else if (Pos > Offset)
			Lo = SeqIndex + 1;
		else
			Hi = SeqIndex - 1;
		}

// Anomaly -- this point is reached (rarely)
// when alignment extends into inter-chr padding
	return UINT32_MAX;
	}

unsigned UFIndex::GetSeqLength(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_SeqLengths));
	unsigned L = m_SeqLengths[SeqIndex];
	return L;
	}

unsigned UFIndex::GetSeqLength(const string &Label) const
	{
	const unsigned SeqCount = SIZE(m_Labels);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (Label == m_Labels[SeqIndex])
			return m_SeqLengths[SeqIndex];
		}
	Die("GetSeqLength(%s)", Label.c_str());
	return 0;
	}

unsigned UFIndex::GetRow(uint64 Slot, uint32 *PosVec) const
	{
	byte T  = GetTally(Slot);
	if (TallyOther(T))
		return 0;

	uint32 Pos = GetPos(Slot);
	uint64 Slot2 = Slot;
	unsigned K = 0;
	for (;;)
		{
		byte T  = GetTally(Slot2);
		uint32 Pos = GetPos(Slot2);
		PosVec[K++] = Pos;
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			assert(Slot2 == Slot);
			return K;
			}
		if (K == m_MaxIx)
			return K;

#if DEBUG
		{
		if (Slot2 == Slot)
			assert(TallyMine(T));
		else
			assert(TallyOther(T));
		}
#endif
		if (T == TALLY_END)
			return K;
		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			uint32 PosA = GetPos(SlotA);
			PosVec[K-1] = PosA;
#if DEBUG
			{
			byte TA = GetTally(SlotA);
			assert(TA == TALLY_NEXT_LONG_OTHER);
			}
#endif
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	Die("GetRow(%" PRIx64 ")", Slot);
	return 0;
	}

unsigned UFIndex::GetRow_Validate(uint64 Slot, uint32 *PosVec) const
	{
	byte T  = GetTally(Slot);
	if (TallyOther(T))
		return 0;

	uint32 Pos = GetPos(Slot);
	uint64 Slot2 = Slot;
	unsigned K = 0;
	for (;;)
		{
		byte T  = GetTally(Slot2);
		uint32 Pos = GetPos(Slot2);
		PosVec[K++] = Pos;
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			assert(Slot2 == Slot);
			return K;
			}
		if (K == m_MaxIx)
			return K;
		if (Slot2 == Slot)
			asserta(TallyMine(T));
		else
			asserta(TallyOther(T));
		if (T == TALLY_END)
			return K;
		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			uint32 PosA = GetPos(SlotA);
			PosVec[K-1] = PosA;
			byte TA = GetTally(SlotA);
			asserta(TA == TALLY_NEXT_LONG_OTHER);
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	Die("GetRow(%" PRIx64 ")", Slot);
	return 0;
	}

unsigned UFIndex::GetRow_Blob(uint64 Slot, const byte *ptrBlob,
  uint32 *PosVec) const
	{
	byte T = *ptrBlob;
	if (TallyOther(T))
		return 0;
	uint64 Slot2 = Slot;
	uint32 Pos = *(uint32 *) (ptrBlob + 1);
	unsigned K = 0;
	for (;;)
		{
		if (K > 0)
			{
			T  = GetTally(Slot2);
			Pos = GetPos(Slot2);
			}
		assert(Pos < m_SeqDataSize);
		PosVec[K++] = Pos;
		if (K == m_MaxIx)
			return K;
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			return 1;
			}

#if DEBUG
		{
		if (K == 1)
			assert(TallyMine(T));
		else
			assert(TallyOther(T));
		}
#endif
		if (T == TALLY_END)
			return K;

		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			uint32 PosA = GetPos(SlotA);
			PosVec[K-1] = PosA;
#if DEBUG
			{
			byte TA = GetTally(SlotA);
			assert(TA == TALLY_NEXT_LONG_OTHER);
			}
#endif
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	Die("GetRow_Blob(%" PRIx64 ")", Slot);
	return 0;
	}

uint64 UFIndex::FindEndOfList(uint64 Slot) const
	{
	byte T  = GetTally(Slot);
	uint32 Pos = GetPos(Slot);
	assert(TallyMine(T));

	uint64 Slot2 = Slot;
	unsigned PosCount = 0;
	for (;;)
		{
		byte T  = GetTally(Slot2);
		uint32 Pos = GetPos(Slot2);
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			assert(Slot2 == Slot);
			return Slot2;
			}

		if (Slot2 == Slot)
			assert(TallyMine(T));
		else
			assert(TallyOther(T));
		if (T == TALLY_END)
			return Slot2;
		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	Die("FindEndOfList failed");
	return UINT64_MAX;
	}

unsigned UFIndex::FindFreeSlot(uint64 Slot)
	{
	for (unsigned i = 1; i < MAX_LINK_STEP; ++i)
		{
		uint64 Slot2 = (Slot + i)%m_SlotCount;
		byte n = m_SlotCounts[Slot2];
		if (n > 0 && n <= m_MaxIx)
			continue;
		byte T = GetTally(Slot2);
		if (T == TALLY_FREE)
			return i;
		}
	return UINT_MAX;
	}

void UFIndex::LogRow(uint64 Slot) const
	{
	uint32 PosVec[128];
	unsigned K = GetRow(Slot, PosVec);
	for (unsigned k = 0; k < K; ++k)
		{
		uint32 Pos = PosVec[k];
		string s;
		GetStr(Pos, m_WordLength, s);
		string Label;
		unsigned Coord = PosToCoord(Pos, Label);
		Log(" %s(%u=%s:%u)", s.c_str(), Pos, Label.c_str(), Coord+1);
		}
	}

void UFIndex::LogCountHist() const
	{
	asserta(m_SlotCounts != 0);
	Log("\n");
	Log("CountHist\n");
	vector<unsigned> CountHist(256);
	unsigned Total = 0;
	for (uint64 Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		unsigned n = m_SlotCounts[Slot];
		Total += n;
		++(CountHist[n]);
		}
	unsigned Maxi = 0;
	for (unsigned i = 0; i < 256; ++i)
		if (CountHist[i] > 0)
			Maxi = i;
	for (unsigned i = 0; i <= Maxi; ++i)
		{
		unsigned CH = CountHist[i];
		Log("[%3u]  %10u", i, CH);
		if (i > 0)
			Log("  %7.2f%%", GetPct(CH, Total));
		Log("\n");
		}
	Log("Total %u\n", Total);
	}

void UFIndex::LogTally(byte T) const
	{
	Log("T=%02x/", T);
	if (T == TALLY_FREE)
		Log("Free");
	else if (T == TALLY_END)
		Log("End");
	else if (T == TALLY_BOTH1)
		Log("Both1");
	else if (T == TALLY_PLUS1)
		Log("Plus1");
	else if (T == TALLY_NEXT_LONG_MINE)
		Log("LongM");
	else if (T == TALLY_NEXT_LONG_OTHER)
		Log("LongO");
	else
		Log("%u", GetTallyNext(T));
	}

void UFIndex::LogList(uint64 Slot) const
	{
	Log("\n");
	Log("LogList(%x)\n", Slot);
	byte T  = GetTally(Slot);
	uint32 Pos = GetPos(Slot);
	if (TallyOther(T))
		{
		Log("%10" PRIx64, Slot);
		Log("  %10u", Pos);
		Log("  %02X", T);
		LogTally(T);
		Log("\n");
		return;
		}

	uint64 Slot2 = Slot;
	unsigned PosCount = 0;
	for (;;)
		{
		byte T  = GetTally(Slot2);
		uint32 Pos = GetPos(Slot2);
		Log("%10" PRIx64, Slot2);
		Log("  %10u", Pos);
		Log("  %02X  ", T);
		LogTally(T);
		if (T != TALLY_NEXT_LONG_MINE && T != TALLY_NEXT_LONG_OTHER)
			{
			assert(Pos < m_SeqDataSize);
			string s;
			GetStr(Pos, m_WordLength, s);
			Log("  [%08" PRIx64 "]%u=%s", Slot, Pos, s.c_str());
			Log("\n");
			}
		if (T == TALLY_PLUS1 || T == TALLY_BOTH1)
			{
			assert(Slot2 == Slot);
			return;
			}

		if (Slot2 == Slot)
			assert(TallyMine(T));
		else
			assert(TallyOther(T));
		if (T == TALLY_END)
			return;
		if (T == TALLY_NEXT_LONG_MINE || T == TALLY_NEXT_LONG_OTHER)
			{
			uint32 StepA = Pos & 0xffff;
			uint32 StepB = (Pos >> 16);
			uint64 SlotA = (Slot2 + StepA)%m_SlotCount;
			Slot2 = (SlotA + StepB)%m_SlotCount;
			Log(" ...StepA=%u, StepB=%u, SlotA %" PRIx64 ", Slot2 -> %" PRIx64 "\n", StepA, StepB, SlotA, Slot2);
			Log("SlotA:\n");
			uint32 PosA = GetPos(SlotA);
			byte TA = GetTally(SlotA);
			Log("%10" PRIx64, SlotA);
			Log("  %10u", PosA);
			Log("  %02X  ", TA);
			LogTally(TA);
			string s;
			GetStr(PosA, m_WordLength, s);
			Log("  [%08" PRIx64 "]%u=%s", SlotA, PosA, s.c_str());
			Log("\n");
			}
		else
			{
			byte Next = GetTallyNext(T);
			assert(Next > 0 && Next <= TALLY_MAX_NEXT);
			Slot2 = (Slot2 + Next)%m_SlotCount;
			}
		}
	}
