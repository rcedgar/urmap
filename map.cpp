#include "myutils.h"
#include "ufindex.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastqseqsource.h"
#include <inttypes.h>
#include "state1.h"
#include "sort.h"
#include "outfiles.h"

static void MapThread(const UFIndex &UFI, FASTQSeqSource &FSS)
	{
	State1 UD;
	UD.SetUFI(UFI);
	SeqInfo *Query = ObjMgr::GetSeqInfo();

	for (;;)
		{
		bool Ok = FSS.GetNext(Query);
		if (!Ok)
			break;
		UD.Search(Query);
		UD.Output1();
		}
	}

void cmd_map()
	{
	const string &FileName = opt(map);

	void InitGlobals(bool Nucleo);
	InitGlobals(true);

	unsigned Method = 6;
	if (opt(veryfast))
		Method = 7;
	State1::SetMethod(Method);

	FASTQSeqSource FSS;
	FSS.Open(FileName);

	UFIndex UFI;
	UFI.FromFile(opt(ufi));
	State1::WriteSAMHeader(g_fSAM, UFI);

	if (opt(veryfast) && UFI.m_MaxIx > 3)
		Warning("index not optimal for -veryfast");

	unsigned ThreadCount = GetRequestedThreadCount();
	Progress("%u threads\n", ThreadCount);
	Progress("Mapping %s\n", FileName.c_str());
	State1::m_StartSecs = GetElapsedSecs();

	if (FSS.m_LR.m_gz)
		Progress("Mapping unpaired %s\n", FileName.c_str());
	else
		ProgressFile(FSS.m_LR.m_f, "Mapping unpaired ", FileName);
#pragma omp parallel num_threads(ThreadCount)
	{
	MapThread(UFI, FSS);
	}
	if (!FSS.m_LR.m_gz)
		ProgressDone();
	ProgressSetMsg("Closing");

	State1::HitStats();
	}
