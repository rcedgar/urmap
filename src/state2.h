#ifndef state2_h
#define state2_h

#include "state1.h"

class State2
	{
public:
	static unsigned m_Minq;

public:
	State1 m_UD_Fwd;
	State1 m_UD_Rev;
	unsigned m_k_Fwd;
	unsigned m_k_Rev;

	vector<unsigned> m_Paired_HitIndexVec_Fwd;
	vector<unsigned> m_Paired_HitIndexVec_Rev;
	vector<int> m_Paired_TotalScoreVec;
	vector<unsigned> m_Paired_TLVec;

	int m_BestPairScore;
	int m_SecondBestPairScore;
	unsigned m_BestPairIndex;
	unsigned m_SecondPairIndex;

	GoBuff<uint32> m_B1_QPosVecf;
	GoBuff<uint32> m_B1_QPosVecr;
	GoBuff<bool> m_B1_PlusVecf;
	GoBuff<bool> m_B1_PlusVecr;
	GoBuff<uint32> m_B1_DBPosVecf;
	GoBuff<uint32> m_B1_DBPosVecr;
	unsigned m_B1Sizef;
	unsigned m_B1Sizer;
	int m_TermPairScorePhase1;

public:
	static unsigned m_Method;

public:
	State2()
		{
		m_B1Sizef = 0;
		m_B1Sizer = 0;
		}

public:
	void SetUFI(const UFIndex &UFI);
	void Search(SeqInfo *Query_Fwd, SeqInfo *Query_Rev);

	void Search4(SeqInfo *Query_Fwd, SeqInfo *Query_Rev);
	bool ExtendBoth1Pair4(uint32 QPosf, uint32 DBPosf, bool Plusf,
	  uint32 QPosr, uint32 DBPosr);

	void Search5(SeqInfo *Query_Fwd, SeqInfo *Query_Rev);
	bool ExtendBoth1Pair5(uint32 QPosf, uint32 DBPosf, bool Plusf,
	  uint32 QPosr, uint32 DBPosr);

	void AdjustTopHitsAndMapqs();
	void FindPairs();
	void ScanPair();
	void LogPairs() const;
	bool AcceptTemplateLength(int64 TL) const;

	void LogPos(uint32 DBPos) const { m_UD_Fwd.LogPos(DBPos); }

	void GetChrPos(const UFIHit *Hit, string &ChrLabel,
	  uint32 &ChrPos) const;
	void GetPairPosStr(const UFIHit *Hit1, const UFIHit *Hit2,
	  string &s) const;
	void GetPairPosStr1(const UFIHit *Hit, bool Fwd, string &s) const;
	void GetInfoStr(string &s) const;
	unsigned GetTemplateLength(const UFIHit *Hit1, const UFIHit *Hit2) const;

	const UFIHit *GetTopHit1() const { return m_UD_Fwd.m_TopHit; }
	const UFIHit *GetTopHit2() const { return m_UD_Rev.m_TopHit; }
	const UFIHit *GetSecondHit1() const  { return m_UD_Fwd.m_SecondHit; }
	const UFIHit *GetSecondHit2() const  { return m_UD_Rev.m_SecondHit; }

	void SetSAM2();
	void Output2();
	void OutputSAM2();
	void OutputTab2();

public:
	static uint32 GetPairedFlags(bool First, bool RevComp, bool MateRevComp, bool MateUnmapped);
	};

#endif // state2_h
