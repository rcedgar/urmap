#ifndef ufindex_h
#define ufindex_h

#include "seqdb.h"
#include "readsimbench.h"

#define TSLOT		0

class SeqInfo;

const unsigned UFI_MAGIC1 = MAGIC('U', 'F', 'I', '1');
const unsigned UFI_MAGIC2 = MAGIC('U', 'F', 'I', '2');
const unsigned UFI_MAGIC3 = MAGIC('U', 'F', 'I', '3');
const unsigned UFI_MAGIC4 = MAGIC('U', 'F', 'I', '4');
const unsigned UFI_MAGIC5 = MAGIC('U', 'F', 'I', '5');

/***
Tally = Mnnnnnnn
	M		My bit, 1=mine, 0=other.
	nnnnn	7-bit offset of next in list, 1111111=end.
***/

static const byte TALLY_FREE			= 0;		// 0b00000000;  // 0
static const byte TALLY_END				= 127;		// 0b01111111;  // 127

static const byte TALLY_MY_BIT			= 128;		// 0b10000000;  // 128
static const byte TALLY_PLUS1			= 254;		// 0b11111110;  // 254
static const byte TALLY_BOTH1			= 255;		// 0b11111111;  // 255
static const byte TALLY_MAX_NEXT		= 124;		// 0b01111100;  // 124
static const byte TALLY_NEXT_LONG		= 125;		// 0b01111101;  // 125
static const byte TALLY_NEXT_MASK		= 127;		// 0b01111111;  // 127

static const byte TALLY_NEXT_LONG_MINE	= 253;		// 0b11111101;  // 253
static const byte TALLY_NEXT_LONG_OTHER	= 125;		// 0b01111101;  // 125

static const unsigned MAX_LINK_STEP = 0xffff;

//// http://xoshiro.di.unimi.it/splitmix64.c
//uint64_t
//splittable64(uint64_t x)
//{
//    x ^= x >> 30;
//    x *= UINT64_C(0xbf58476d1ce4e5b9);
//    x ^= x >> 27;
//    x *= UINT64_C(0x94d049bb133111eb);
//    x ^= x >> 31;
//    return x;
//}

static inline uint64 murmur64(uint64 h)
	{
	h ^= (h >> 33);
	h *= 0xff51afd7ed558ccdL;
	h ^= (h >> 33);
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= (h >> 33);
	return h;
	}

static inline uint64 WordToSlot(uint64 Word, uint64 SlotCount)
	{
	uint64 h = murmur64(Word);
	uint64 Slot = h%SlotCount;
	return Slot;
	}

static inline bool TallyEnd(byte T)
	{
	return T == TALLY_END;
	}

static inline bool TallyMine(byte T)
	{
	return (T & TALLY_MY_BIT) != 0;
	}

static inline bool TallyHasPos(byte T)
	{
	return T != TALLY_FREE && T != TALLY_NEXT_LONG_MINE;
	}

static inline byte TallyOther(byte T)
	{
	return (T & TALLY_MY_BIT) == 0;
	}

static inline byte GetTallyNext(byte T)
	{
	assert(T != TALLY_END);
	byte Next = (T & TALLY_NEXT_MASK);
	assert(Next != 0);
	return Next;
	}

const unsigned PADGAP = 32;

class UFIndex
	{
public:
	uint32 m_WordLength;
	uint32 m_MaxIx;
	byte *m_Blob;
	uint64 m_SlotCount;

// Misc helpers
	uint64 m_ShiftMask;
	uint64 m_LowerCaseShiftMask;

// Sequence data
	byte *m_SeqData;
	uint32 m_SeqDataSize;

// Sequence directory
	vector<string> m_Labels;
	vector<uint32> m_SeqLengths;
	vector<uint32> m_Offsets;

// Temp data needed for build
	byte *m_SlotCounts;
	byte *m_SlotCounts_Minus;

// Build stats
	uint64 m_GenomeSize;
	uint64 m_IndexedPosCount;
	unsigned m_TruncatedCount;

	unsigned m_StartPos;
	unsigned m_EndPos;

public:
	UFIndex()
		{
		m_WordLength = 0;
		m_MaxIx = 0;
		m_Blob = 0;
		m_SlotCount = 0;
		m_ShiftMask = 0;
		m_SeqData = 0;
		m_GenomeSize = 0;
		m_SeqDataSize = 0;
		m_IndexedPosCount = 0;
		m_SlotCounts = 0;
		m_SlotCounts_Minus = 0;
		m_TruncatedCount = 0;

		m_StartPos = UINT_MAX;
		m_EndPos = UINT_MAX;
		}

	void SetTally(uint64 Slot, byte T)
		{
#if TSLOT
		if (Slot == TSLOT)
			{
			Log("SetTally(%" PRIx64 ") ", Slot, T);
			LogTally(T);
			Log("\n");
			}
#endif
		m_Blob[5*Slot] = T;
		}

	byte GetTally(uint64 Slot) const
		{
		byte T = m_Blob[5*Slot];
		return T;
		}

	void SetPos(uint64 Slot, uint32 Pos)
		{
#if TSLOT
		if (Slot == TSLOT)
			Log("SetPos(%" PRIx64 ", %x)\n", Slot, Pos);
#endif
		*(uint32 *) (m_Blob + 5*Slot + 1) = Pos;
		}

	uint32 GetPos(uint64 Slot) const
		{
		uint32 Pos = *(uint32 *) (m_Blob + 5*Slot + 1);
		return Pos;
		}

	void GetBlob(uint64 Slot, byte *ptrBlob) const
		{
		memcpy(ptrBlob, m_Blob + 5*Slot, 5);
		}

	void SetTallyNext(uint64 Slot, byte Next)
		{
		byte T = GetTally(Slot);
		T = (T & TALLY_MY_BIT) | Next;
		SetTally(Slot, T);
		}

	void LogMe() const;
	void LogSeqDict() const;
	void LogSlots(uint64 SlotLo, uint64 N) const;
	void ToFile(FILE *f) const;
	void ToFile(const string &FileName) const;
	void FromFile(FILE *f);
	void FromFile(const string &FileName);
	void Validate() const;
	void FreeHashTable() { myfree(m_Blob), m_Blob = 0; }
	void FreeTempBuildData();
	void LogSlot(uint64 Slot) const;
	void ValidateSlot(uint64 Slot) const;
	void GetCollisionCount(uint64 &n, uint64 &N) const;
	const char *GetStr(uint32 Pos, unsigned n, string &s) const;
	const char *WordToStr(uint64 Word, string &s) const;
	uint64 SeqToWord(const byte *Seq) const;
	uint64 GetWord(uint32 Pos) const;
	void CountSlots();
	void CountSlots_Minus();
	//void IncCountSlots_Minus();
	void MakeIndex();
	void LogCountHist() const;
	void CountIndexedWords(unsigned &IndexedCount,
	  unsigned &NotIndexedCount, unsigned &WildcardCount);
	unsigned GetSeqIndex(const string &Label) const;
	const byte *GetSeq(unsigned SeqIndex) const;
	uint32 PosToCoord(uint32 Pos, string &Label) const;
	uint32 PosToCoordL(uint32 Pos, string &Label, unsigned &L) const;
	uint32 CoordToPos(const string &Label, uint32 Coord) const;
	void GetUserPosStr(uint32 Pos, string &s) const;
	unsigned GetSeqLength(const string &Label) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;
	void SetShiftMask();
	void ReadSeqData(const string &FastaFileName);
	void UpdateSlot(uint64 Slot, uint32 Pos);
	unsigned GetRow(uint64 Slot, uint32 *PosVec) const;
	unsigned GetRow_Validate(uint64 Slot, uint32 *PosVec) const;
	unsigned GetRow_Blob(uint64 Slot, const byte *ptrBlob, uint32 *PosVec) const;
	void LogStats();
	uint64 FindEndOfList(uint64 Slot) const;
	unsigned FindFreeSlot(uint64 Slot);
	void LogRow(uint64 Slot) const;
	void TruncateSlot(uint64 Slot);
	void LogTally(byte T) const;
	void LogList(uint64 Slot) const;
	};

#endif // ufindex_h
