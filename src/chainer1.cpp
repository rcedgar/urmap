#include "myutils.h"
#include "chainer1.h"
#include "sort.h"

#define TRACE	0
#define TEST	0

#if	BRUTE
#include "combo.h"
#endif

static const float MINUS_INFINITY = -9e9f;

void SortBPVecInPlace(BPData *BPVec, unsigned N);

void Chainer1::Clear()
	{
	m_BPs.Free();
	m_ChainIndexes.Free();
	m_ChainScores.Free();
	m_TB.Free();
	}

const unsigned *Chainer1::Chain(const unsigned *Los, const unsigned *His,
  float *Scores, unsigned N, unsigned &ChainLength)
	{
	if (N == 0)
		{
		ChainLength = 0;
		return 0;
		}

#if	TRACE
	Log("Chainer1::Chain(N=%u)\n", N);
#endif
	m_BPs.Alloc(2*N);

	BPData *BPVec = m_BPs.Data;
	BPData *BP = BPVec;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Lo = Los[i];
		unsigned Hi = His[i];
		asserta(Hi >= Lo);

		BP->Index = i;
		BP->IsLo = true;
		BP->Pos = Lo;
		BP++;

		BP->Index = i;
		BP->IsLo = false;
		BP->Pos = Hi;
		BP++;
		}
#if	0 // TRACE
	{
	Log("BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (unsigned i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	SortBPVecInPlace(BPVec, 2*N);

#if	TRACE
	{
	Log("Sorted BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (unsigned i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	m_TB.Alloc(N);
	m_ChainScores.Alloc(N);
	unsigned *TB = m_TB.Data;
	float *ChainScores = m_ChainScores.Data;

	asserta(BPVec[0].IsLo);
	unsigned Index0 = BPVec[0].Index;
	unsigned BestChainEnd = UINT_MAX;
	TB[0] = UINT_MAX;

	ChainScores[0] = MINUS_INFINITY;

	for (unsigned i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];

		assert(BP.Index < N);
		float Score = Scores[BP.Index];

		if (BP.IsLo)
			{
			TB[BP.Index] = BestChainEnd;
			if (BestChainEnd == UINT_MAX)
				ChainScores[BP.Index] = Score;
			else
				ChainScores[BP.Index] = ChainScores[BestChainEnd] + Score;
			}
		else
			{
			if (BestChainEnd == UINT_MAX || ChainScores[BP.Index] > ChainScores[BestChainEnd])
				BestChainEnd = BP.Index;
			}
		}

	asserta(BestChainEnd < N);

#if	TRACE
	{
	Log("\n");
	Log("BestChainEnd %u, Score %.1f\n", BestChainEnd, ChainScores[BestChainEnd]);
	Log("Index  ChainScore     TB\n");
	Log("-----  ----------  -----\n");
	for (unsigned i = 0; i < N; ++i)
		{
		Log("%5u", i);
		float Score = ChainScores[i];
		if (Score == MINUS_INFINITY)
			Log("  %10.10s", "*");
		else
			Log("  %10.1f", Score);
		unsigned t = TB[i];
		if (t == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", t);
		Log("\n");
		}
	}
#endif

	m_ChainIndexes.Alloc(N);
	unsigned *ChainIndexes = m_ChainIndexes.Data;
	ChainLength = 0;
	unsigned Index = BestChainEnd;
	for (;;)
		{
		asserta(ChainLength < N);
		ChainIndexes[N - ++ChainLength] = Index;
		asserta(Index < N);
		Index = TB[Index];
		if (Index == UINT_MAX)
			break;
		}
	const unsigned *ChainPtr = ChainIndexes + N - ChainLength;

#if	TRACE
	{
	Log("\n");
	Log("Chain:\n");
	Log("Index     Lo     Hi   Score\n");
	Log("-----  -----  -----  ------\n");
	float Sum = 0.0;
	for (unsigned i = 0; i < ChainLength; ++i)
		{
		unsigned Index = ChainPtr[i];
		asserta(Index < N);
		Log("%5u", Index);
		Log("  %5u", Los[Index]);
		Log("  %5u", His[Index]);
		Log("  %6.1f", Scores[Index]);
		Sum += Scores[Index];
		Log("\n");
		}
	Log("Sum %.1f\n", Sum);
	}
#endif

	assert(IsValidChain(Los, His, N, ChainPtr, ChainLength));
	return ChainPtr;
	}

bool Chainer1::IsValidChain(const unsigned *Los, const unsigned *His, unsigned N,
  const unsigned *Chain, unsigned ChainLength)
	{
	asserta(ChainLength > 0);
	for (unsigned i = 0; i < ChainLength; ++i)
		{
		unsigned Index = Chain[i];
		asserta(Index < N);
		asserta(Los[Index] <= His[Index]);
		if (i > 0)
			{
			unsigned PrevIndex = Chain[i-1];
			if (Los[Index] <= His[PrevIndex])
				return false;
			}
		}
	return true;
	}

float Chainer1::GetChainScore(const unsigned *Los, const unsigned *His,
  const float *Scores, unsigned N, const unsigned *Chain, unsigned ChainLength)
	{
	float Sum = 0.0;
	for (unsigned i = 0; i < ChainLength; ++i)
		{
		unsigned Index = Chain[i];
		assert(Index < N);
		Sum += Scores[Index];
		}
	return Sum;
	}

#if	BRUTE

static vector<unsigned> g_BestChain;
static const unsigned *g_Los;
static const unsigned *g_His;
static unsigned g_N;
static const float *g_Scores;
static float g_BestScore;

static void OnPerm(const vector<unsigned> &v)
	{
	unsigned n = SIZE(v);
	if (!Chainer1::IsValidChain(g_Los, g_His, g_N, v.data(), n))
		return;

	float Sum = 0.0;
	for (unsigned i = 0; i < n; ++i)
		{
		unsigned k = v[i];
		Sum += g_Scores[k];
		}

	if (Sum > g_BestScore)
		{
		g_BestScore = Sum;
		g_BestChain = v;
		}
	}

const unsigned *Chainer1::ChainBrute(const unsigned *Los, const unsigned *His,
  float *Scores, unsigned N, unsigned &ChainLength)
	{
	g_BestChain.clear();
	g_Los = Los;
	g_His = His;
	g_Scores = Scores;
	g_N = N;
	g_BestScore = MINUS_INFINITY;
	EnumPowerSetPerms(N, OnPerm);
	ChainLength = SIZE(g_BestChain);
	const unsigned *ChainPtr = g_BestChain.data();
#if TRACE
	{
	Log("\n");
	Log("ChainBrute:\n");
	Log("Index     Lo     Hi   Score\n");
	Log("-----  -----  -----  ------\n");
	float Sum = 0.0;
	for (unsigned i = 0; i < ChainLength; ++i)
		{
		unsigned Index = ChainPtr[i];
		asserta(Index < N);
		Log("%5u", Index);
		Log("  %5u", Los[Index]);
		Log("  %5u", His[Index]);
		Log("  %6.1f", Scores[Index]);
		Sum += Scores[Index];
		Log("\n");
		}
	Log("Sum %.1f\n", Sum);
	}
#endif
	asserta(IsValidChain(Los, His, N, ChainPtr, ChainLength));
	return ChainPtr;
	}
#endif

#if	TEST
const unsigned MaxTries = 100;
const unsigned MinCount = 1;
const unsigned MaxCount = 8;
const unsigned MinLen = 1;
const unsigned MaxLen = 100;
const unsigned MaxPos = 100;
const unsigned MinScore = 1;
const unsigned MaxScore = 100;
const unsigned RandSeed = 0;

static void GetRandomLoHi(unsigned MaxPos, unsigned MinLen, unsigned MaxLen,
  unsigned MinScore, unsigned MaxScore, unsigned &Lo, unsigned &Hi, float &Score)
	{
	asserta(MinLen <= MaxLen);
	asserta(MinScore <= MaxScore);

	Lo = unsigned(rand()%MaxPos);
	unsigned Length = MinLen + unsigned(rand()%(MaxLen - MinLen + 1));
	Hi = Lo + Length - 1;
	Score = float(MinScore + unsigned(rand()%(MaxScore - MinScore + 1)));
	}

static unsigned GetRandomLoHis(
  unsigned MinCount, unsigned MaxCount,
  unsigned MaxPos,
  unsigned MinLen, unsigned MaxLen,
  unsigned MinScore, unsigned MaxScore,
  unsigned *Los, unsigned *His, float *Scores)
	{
	asserta(MinCount <= MaxCount);
	unsigned Count = MinCount + unsigned(rand()%(MaxCount - MinCount + 1));
	for (unsigned i = 0; i < Count; ++i)
		GetRandomLoHi(MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los[i], His[i], Scores[i]);
	return Count;
	}

void cmd_test()
	{
	srand(RandSeed);

	Chainer1 C;

	unsigned *Los = myalloc(unsigned, MaxCount);
	unsigned *His = myalloc(unsigned, MaxCount);
	float *Scores = myalloc(float, MaxCount);

	for (unsigned Try = 0; Try < MaxTries; ++Try)
		{
		ProgressStep(Try, MaxTries, "Testing");
		unsigned N = GetRandomLoHis(MinCount, MaxCount, MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los, His, Scores);

		unsigned ChainLength;
		const unsigned *Chain = C.Chain(Los, His, Scores, N, ChainLength);
		float Score = Chainer1::GetChainScore(Los, His, Scores, N, Chain, ChainLength);

#if	BRUTE
		unsigned ChainLengthBrute;
		const unsigned *ChainBrute = C.ChainBrute(Los, His, Scores, N, ChainLengthBrute);
		float BruteScore = Chainer1::GetChainScore(Los, His, Scores, N, ChainBrute, ChainLengthBrute);
		asserta(feq(Score, BruteScore));

		Log("N %u, chain %u, brute chain %u, Score %.1f, brute %.1f\n",
		  N, ChainLength, ChainLengthBrute, Score, BruteScore);
#else
		Log("N %u, chain %u, Score %.1f\n", N, ChainLength, Score);
#endif
		}
	}
#endif // TEST
