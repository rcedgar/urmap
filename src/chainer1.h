#ifndef chainer1_h
#define chainer1_h

#include "hsp.h" // for BPData
#include "gobuff.h"

#define BRUTE	0

class Chainer1
	{
public:
	GoBuff<BPData> m_BPs;
	GoBuff<unsigned> m_ChainIndexes;
	GoBuff<float> m_ChainScores;
	GoBuff<unsigned> m_TB;

public:
	const unsigned *Chain(const unsigned *Los, const unsigned *His,
	  float *Scores, unsigned N, unsigned &ChainLength);
	void Clear();

	static bool IsValidChain(const unsigned *Los, const unsigned *His, unsigned N,
	  const unsigned *Chain, unsigned ChainLength);
	static float GetChainScore(const unsigned *Los, const unsigned *His, const float *Scores,
	  unsigned N, const unsigned *Chain, unsigned ChainLength);

#if	BRUTE
	const unsigned *ChainBrute(const unsigned *Los, const unsigned *His,
	  float *Scores, unsigned N, unsigned &ChainLength);
#endif
	};

#endif // chainer1_h
