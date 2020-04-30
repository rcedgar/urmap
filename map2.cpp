#include "myutils.h"
#include "ufindex.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastqseqsource.h"
#include <inttypes.h>
#include "state2.h"
#include "omplock.h"
#include "outfiles.h"

static void UFIMapThread(const UFIndex &UFI, FASTQSeqSource &SS1, 
  FASTQSeqSource &SS2)
	{
	State2 UD;
	UD.m_UD_Fwd.SetMethod(6);
	UD.m_UD_Rev.SetMethod(6);
	if (UD.m_Method == 5)
		{
		UD.m_UD_Fwd.GLOBAL_BAND_RADIUS = 4;
		UD.m_UD_Rev.GLOBAL_BAND_RADIUS = 4;
		}
	UD.SetUFI(UFI);
	SeqInfo *Query1 = ObjMgr::GetSeqInfo();
	SeqInfo *Query2 = ObjMgr::GetSeqInfo();
	for (;;)
		{
		Lock();
		bool Ok1 = SS1.GetNext(Query1);
		bool Ok2 = SS2.GetNext(Query2);
		Unlock();
		if (Ok1 != Ok2)
			Die("Premature end of file in FASTQ%c", Ok1 ? '2' : '1');
		if (!Ok1)
			break;
		UD.Search(Query1, Query2);
		}
	}

void cmd_map2()
	{
	const string &FwdFileName = opt(map2);
	if (!optset_reverse)
		Die("-reverse required");
	const string &RevFileName = opt(reverse);

	State2::m_Method = 4;
	if (opt(veryfast))
		State2::m_Method = 5;

	void InitGlobals(bool Nucleo);
	InitGlobals(true);

	FASTQSeqSource SS1;
	FASTQSeqSource SS2;
	SS1.Open(FwdFileName);
	SS2.Open(RevFileName);

	UFIndex UFI;
	UFI.FromFile(opt(ufi));

	//{
	//byte Blob[1024];
	//uint32 PosVec[64];
	//UFI.GetRow_Blob(18446744073709551615, Blob, PosVec);
	//return;
	//}

	State1::WriteSAMHeader(g_fSAM, UFI);

	if (opt(veryfast) && UFI.m_MaxIx > 3)
		Warning("Genome index not optimal for -veryfast");

	unsigned ThreadCount = GetRequestedThreadCount();

	State1::m_StartSecs = GetElapsedSecs();
	State1::m_Minq = opt(minq);
	State2::m_Minq = opt(minq);

	if (SS1.m_LR.m_gz)
		Progress("Mapping paired %s\n", FwdFileName.c_str());
	else
		ProgressFile(SS1.m_LR.m_f, "Mapping (paired)", FwdFileName);
#pragma omp parallel num_threads(ThreadCount)
	{
	UFIMapThread(UFI, SS1, SS2);
	}
	ProgressDone();

	State1::HitStats();
	}
