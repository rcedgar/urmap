#include "myutils.h"
#include "diagbox.h"
#include "alnparams.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "state1.h"
#include "xdpmem.h"

void TraceBackBitMem(XDPMem &m_XDPMem, unsigned LA, unsigned LB, char State, PathInfo &PI);

float State1::Viterbi(const byte *A, unsigned LA,
  const byte *B, unsigned LB, bool Left, bool Right, PathInfo &PI)
	{
	if (LA == 0 || LB == 0)
		{
		if (LA == 0 && LB == 0)
			{
			PI.SetEmpty();
			return 0.0f;
			}
		if (LA == 0 && LB > 0)
			{
			PI.Alloc(LB+1);
			PI.AppendIs(LB);
			float Score = float(GAP_OPEN_SCORE + (LB-1)*GAP_EXT_SCORE);
			return Score;
			}
		if (LA > 0 && LB == 0)
			{
			PI.Alloc(LA+1);
			PI.AppendDs(LA);
			float Score = float(GAP_OPEN_SCORE + (LA-1)*GAP_EXT_SCORE);
			return Score;
			}
		Die("Viterbi LA=%u LB=%u", LA, LB);
		}

// Main diagonal
// d = LA - i + j = 1 .. LA+LB-1
// Left term: i=0,j=0 d=LA
// Right term: i=LA-1,j=LB-1 d=LA-(LA-1)+(LB-1) = LB
	unsigned DiagLo = min(LA, LB);
	unsigned DiagHi = max(LA, LB);
	asserta(DiagLo <= DiagHi);
	if (DiagLo > GLOBAL_BAND_RADIUS)
		DiagLo -= GLOBAL_BAND_RADIUS;
	else
		DiagLo = 1;
	DiagHi += GLOBAL_BAND_RADIUS;
	unsigned MaxDiag = LA + LB - 1;
	if (DiagHi > MaxDiag)
		DiagHi = MaxDiag;

// Verify diagonal range includes terminals
	DiagBox Box(LA, LB, DiagLo, DiagHi);
	asserta(Box.InBox(0, 0));
	asserta(Box.InBox(LA-1, LB-1));

	m_XDPMem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);

	float OpenA = float(GAP_OPEN_SCORE);
	float ExtA = float(GAP_EXT_SCORE);
	if (Left)
		{
		OpenA = 0;
		ExtA = 0;
		}

	float *Mrow = m_XDPMem.GetDPRow1();
	float *Drow = m_XDPMem.GetDPRow2();
	byte **TB = m_XDPMem.GetTBBit();

// Use Mrow[j-1] when j=0, so...
	Mrow[-1] = MINUS_INFINITY;

// ? Surely don't need to initialize all entries in vector
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}

// Main loop
	for (unsigned i = 0; i < LA; ++i)
		{
		unsigned Startj, Endj;
		Box.GetRange_j(i, Startj, Endj);
		if (Endj == 0)
			{
			static bool WarningDone = false;
			if (!WarningDone)
				{
				Warning("Endj==0");
				WarningDone = true;
				}
			continue;
			}

		//float OpenB = Startj == 0 ? AP.LOpenB : AP.OpenB;
		//float ExtB = Startj == 0 ? AP.LExtB : AP.ExtB;
		float OpenB = float((Startj == 0 && Left) ? 0 : GAP_OPEN_SCORE);
		float ExtB = float((Startj == 0 && Left) ? 0 : GAP_EXT_SCORE);

		byte a = A[i];
		float I0 = MINUS_INFINITY;
		float M0;
		if (i == 0)
			M0 = 0;
		else
			{
			if (Startj == 0)
				M0 = MINUS_INFINITY;
			else
				M0 = Mrow[int(Startj)-1];
			}

		byte *TBrow = TB[i];
		if (Startj > 0)
			TBrow[int(Startj)-1] = TRACEBITS_IM;

		for (unsigned j = Startj; j < Endj; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
			
			float xM = M0;
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				}
			M0 = Mrow[j];
			Mrow[j] = xM + (a == b ? 1 : MISMATCH_SCORE);
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
		// I0 = DPI[i][j+1]
			}
			
			OpenB = float(GAP_OPEN_SCORE);
			ExtB = float(GAP_EXT_SCORE);
			
			TBrow[j] = TraceBits;
			}

	// Special case for end of Drow[]
		{
	// M0 = DPM[i][LB]
	// Drow[LB] = DPD[i][LB]
		
		TBrow[LB] = 0;
		float md = M0 + GAP_OPEN_SCORE;
		Drow[LB] += GAP_EXT_SCORE;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
	// Drow[LB] = DPD[i+1][LB]
		}
		
		M0 = MINUS_INFINITY;
		OpenA = float(GAP_OPEN_SCORE);
		ExtA = float(GAP_EXT_SCORE);
		}
	
	unsigned Startj, Endj;
	Box.GetRange_j(LA-1, Startj, Endj);
	asserta(Endj == LB);

// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	Mrow[int(Startj)-1] = MINUS_INFINITY;
	float GapOp = float(GAP_OPEN_SCORE);
	float GapEx = float(GAP_EXT_SCORE);
	if (Right)
		{
		GapOp = 0;
		GapEx = 0;
		}
	for (unsigned j = Startj; j < Endj; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		TBrow[j] = 0;
//		float mi = Mrow[int(j)-1] + AP.ROpenA;
//		I1 += AP.RExtA;
		float mi = Mrow[int(j)-1] + GapOp;
		I1 += GapEx;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	TraceBackBitMem(m_XDPMem, LA, LB, State, PI);

	return Score;
	}

#if 0
static void Test(const char *Q, const char *T, bool Left)
	{
	State1 UD;
	PathInfo *PI = ObjMgr::GetPathInfo();

	void LogAlnPretty(const byte *A, const byte *B, const char *Path,
	  bool StripTermGaps);

	float Score = UD.Viterbi((const byte *) Q, ustrlen(Q),
	  (const byte *) T, ustrlen(T), Left, *PI);

	unsigned M, D, I;
	PI->GetCounts(M, D, I);

	Log("============================= %s %.1f =====================\n",
	  Left ? "Left" : "Right", Score);
	Log("\n");
	Log("%s\n", PI->GetPath());
	Log("M %d, D %u, I %u\n", M, D, I);
	LogAlnPretty((const byte *) Q, (const byte *) T, PI->GetPath(), false);
	}

void cmd_test()
	{
	opt(test);
	const char *Q = "GGGGATTAC";
	const char *T = "GGGGATTACA";
	Test(Q, T, true);
	Test(Q, T, false);

	Q = "GGGGATT";
	T = "GGGGATTACA";
	Test(Q, T, true);
	Test(Q, T, false);

	Q =   "GGATTACA";
	T = "GGGGATTACA";
	Test(Q, T, true);
	Test(Q, T, false);
	}
#endif // 0
