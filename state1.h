#ifndef state1_h
#define state1_h

typedef int dpscore;

#include "gobuff.h"
#include "objmgr.h"
#include "ufihit.h"
#include "ufihsp.h"
#include "xdpmem.h"

// const int MATCH_SCORE = 1;
const unsigned MAX_LOCUS_D = 64;
const unsigned LOCUS_DEGLITCH_BIG = 2000;
const int ACCEPT_MIN_MAPQ = 10;
const int SECONDARY_HIT_MAX_DELTA = 12;
const unsigned PRIME_STRIDE = 27;
const unsigned SCANK = 4;
const int MAX_TL = 1000;

class UFIndex;
class SeqInfo;
class AlignResult;

class State1
	{
public:
	static int MISMATCH_SCORE;
	static int GAP_OPEN_SCORE;
	static int GAP_EXT_SCORE;
	static int MIN_HSP_SCORE_PCT;
	static int TERM_HSP_SCORE_PCT_PHASE3;
	static int XDROP;
	static int MAX_PENALTY;
	static int XPHASE1;
	static int XPHASE3;
	static int XPHASE4;
	static unsigned GLOBAL_BAND_RADIUS;
	static unsigned m_Minq;

private:
	static unsigned m_Method;

public:
	UFIHit *m_TopHit;
	UFIHit *m_SecondHit;

	uint32 m_WordLength;
	uint64 m_SlotCount;
	uint64 m_ShiftMask;
	SeqInfo *m_Query;
	GoBuff<UFIHit *> m_Hits;
	GoBuff<UFIHSP *> m_HSPs;
	GoBuff<char> m_SAM;
	GoBuff<char> m_SUM;
	unsigned m_HitCount;
	unsigned m_HSPCount;
	int m_MaxPenalty;
	int m_BestScore;
	int m_SecondBestScore;
	int m_BestHSPScore;
	unsigned m_Mapq;

	string m_MappedTargetLabel;
	uint32 m_MappedTargetPos;

	GoBuff<byte> m_RevCompQuerySeq;

	byte *m_AllocBuff;
	unsigned m_AllocBuffSize;
	unsigned m_AllocBuffOffset;

	bool m_QLabel;
	XDPMem m_XDPMem;
	PathInfo *m_PI;
	PathInfo *m_LeftPI;
	PathInfo *m_RightPI;

/////////////////////////////////////
// Cached vectors
	uint64 *m_SlotsVec_Plus;		// QL
	uint64 *m_SlotsVec_Minus;		// QL
	uint32 *m_PosVec;				// MaxIx
	byte *m_BlobVec_Plus;			// 5*QL
	byte *m_BlobVec_Minus;			// 5*QL
	byte *m_QPosPendingVec_Plus;	// QL
	byte *m_QPosPendingVec_Minus;	// QL
	unsigned m_QPosPendingCount_Plus;
	unsigned m_QPosPendingCount_Minus;
/////////////////////////////////////

public:
	const UFIndex *m_UFI;

	static unsigned m_StartSecs;
	static unsigned m_QueryCount;
	static unsigned m_AcceptCount;
	static unsigned m_RejectCount;
	static unsigned m_NoHitCount;

public:
	State1()
		{
		m_TopHit = 0;
		m_Query = 0;
		m_WordLength = 0;
		m_SlotCount = 0;
		m_ShiftMask = 0;
		m_HitCount = 0;
		m_MaxPenalty = -1;

		m_QLabel = opt(qlabel);
		m_AllocBuffSize = 1024*1024;
		m_AllocBuff = myalloc(byte, m_AllocBuffSize);
		m_AllocBuffOffset = 0;

		m_PI = ObjMgr::GetPathInfo();
		m_LeftPI = ObjMgr::GetPathInfo();
		m_RightPI = ObjMgr::GetPathInfo();

	///////////////////////////////////////////////
		m_SlotsVec_Plus = 0;
		m_SlotsVec_Minus = 0;
		m_PosVec = 0;
		m_BlobVec_Plus = 0;
		m_BlobVec_Minus = 0;
		m_QPosPendingVec_Plus = 0;
		m_QPosPendingVec_Minus = 0;
		m_QPosPendingCount_Plus = 0;
		m_QPosPendingCount_Minus = 0;
	///////////////////////////////////////////////

		}

	bool IsCorrect() const;
	const UFIHit *GetSecondHit() const;

	void AllocHits(unsigned n);
	void AllocHSPs(unsigned n);
	void SetUFI(const UFIndex &UFI);

	void Search(SeqInfo *Query);
	void Search_PE(SeqInfo *Query);
	void Search_Lo();
	void SearchPE_Pending(unsigned k);
	void SearchPE_Pending5(unsigned k);

	void Search1(SeqInfo *Query);
	bool Search1_Lo(SeqInfo *Query);

	void Output1();
	void UpdateHitStats();

	unsigned ExtendScan(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus);
	int ExtendPen(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus);
	bool ExtendExact(uint32 SeedPosQ, uint32 SeedPosDB, bool Plus);

	unsigned CalcMAPQ6() const;

	unsigned OverlapsHSP(uint32 StartPosQ, uint32 StartPosDB, bool Plus);
	UFIHit *OverlapsHit(uint32 DBStartPos, bool Plus);

	int GetBestScore() const { return m_BestScore; }
	int GetSecondBestScore() const { return m_SecondBestScore; }
	int GetSecondBestDelta() const;
	const UFIHit *GetTopHit() const;
	void SetSlotsVec(const byte *Seq, unsigned L, uint64 *SlotsVec) const;
	void SetSlotsVecRC(const byte *Seq, unsigned L, uint64 *SlotsVec) const;
	void SetTallyVec(const uint64 *Slots, unsigned N, byte *TallyVec) const;
	void SetBlobVec(const uint64 *Slots, unsigned N, byte *BlobVec) const;
	unsigned AddHitX(uint32 StartPosDB, bool Plus, int Score, const string &Path);
	void AddHSPX(unsigned SeedPosQ, uint32 SeedPosDB, bool Plus,
	  unsigned Length, int Score);
	unsigned AddHSPScan(unsigned SeedPosQ, uint32 SeedPosDB, bool Plus,
	  unsigned Length, int Score);
	void SetSAM(uint32 Flags, const string &MateTargetLabel,
	  const uint32 MateTargetPos, int TLEN);
	void SetSUM();
	void SetSAM_Unmapped(uint32 Flags);
	unsigned AlignHSP(unsigned HSPIndex);
	unsigned AlignHSP5(unsigned HSPIndex);
	void GetHitOrder(vector<unsigned> &Order) const;
	const UFIHit *GetHit(unsigned Index) const;
	const UFIHSP *GetHSP(unsigned Index) const;
	UFIHit *GetHit(unsigned Index);
	UFIHSP *GetHSP(unsigned Index);
	float Viterbi(const byte *A, unsigned LA, const byte *B, unsigned LB,
	  bool Left, bool Right, PathInfo &PI);
	void GetTopHitIndexes(vector<unsigned> &HitIndexes) const;
	void GetBoth1s(SeqInfo *Query, vector<uint32> &PosVec);
	uint32 GetTopHitPos() const;

	void Scan(uint32 DBLo, unsigned DBSegLength, bool Plus, bool DoVit);
	void ScanSlots(uint32 DBLo, unsigned DBSegLength, bool Plus);
	void GetCIGAR(string &CIGAR) const;
	bool PolishCIGAR(string &CIGAR) const;

	void LogHit(unsigned HitIndex) const;
	void LogHSP(const UFIHSP *HSP) const;
	void LogHits() const;
	void LogHitsSummary() const;
	void LogHSPSummary() const;
	void LogPos(uint32 DBPos) const;
	void PrPos(FILE *f, uint32 DBPos) const;
	void LogState(const string &s = "") const;
	void GetPairLabel(string &Label) const;

	void *AllocBuff(unsigned Bytes)
		{
		void *p = m_AllocBuff + m_AllocBuffOffset;
		m_AllocBuffOffset += Bytes;
		asserta(m_AllocBuffOffset <= m_AllocBuffSize);
		return p;
		}

	void GetChrPos(const UFIHit *Hit, string &ChrLabel,
	  uint32 &ChrPos) const;

/////////////////////////////////////////
	void InitPE(SeqInfo *Query);
	void SetMappedPos();
	//void ResetPE();
	void AllocCachedVectors(unsigned QL);
	unsigned GetFirstBoth1Seed(uint32 &QPos, bool &Plus, uint32 &DBPos);
	unsigned GetNextBoth1Seed(unsigned k, uint32 &QPos, bool &Plus, uint32 &DBPos);
	void GetBoth1Seeds(vector<unsigned> &QPosVec, vector<uint32> &DBPosVec,
	  vector<bool> &PlusVec, vector<uint64> &DiagVec);
	unsigned GetFirstBoth1SeedEx(uint32 &QPos, bool &Plus, uint32 &DBPos);
	unsigned GetNextBoth1SeedEx(unsigned k, uint32 &QPos, bool &Plus, uint32 &DBPos);
/////////////////////////////////////////

public:
	static void HitStats();
	static void SetMethod(unsigned Method);
	static void WriteSAMHeader(FILE *f, const UFIndex &UFI);
	};

#define alloc_buff(t, n) (t *) AllocBuff((n)*sizeof(t))

#endif // state1_h
