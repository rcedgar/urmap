#include "myutils.h"
#include "fastaseqsource.h"
#include "fastqseqsource.h"
#include "seqdbseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "cmd.h"
#include "filetype.h"

#define TIME_LOCKS	0

#if TIME_LOCKS
#include "getticks.h"
static TICKS g_tLocks;
static TICKS g_tUnLocks;
#endif

SeqSource::SeqSource()
	{
	m_SI = ObjMgr::GetSeqInfo();
	m_SeqCount = 0;
	m_DoGetLock = true;
	}

SeqSource::~SeqSource()
	{
	ObjMgr::Down(m_SI);
	}

bool SeqSource::GetNext(SeqInfo *SI)
	{
	if (m_DoGetLock)
		{
#if	TIME_LOCKS
		TICKS t1 = GetClockTicks();
#endif
		LOCK_CLASS();
#if	TIME_LOCKS
		TICKS t2 = GetClockTicks();
		g_tLocks += (t2 - t1);
#endif
		}
	bool Ok = GetNextLo(SI);
	if (m_DoGetLock)
		{
#if	TIME_LOCKS
		TICKS t1 = GetClockTicks();
#endif
		UNLOCK_CLASS();
#if	TIME_LOCKS
		TICKS t2 = GetClockTicks();
		g_tUnLocks += (t2 - t1);
#endif
		}

	if (!Ok)
		{
#if	TIME_LOCKS
		Log("SeqSource locks %.3e, unlocks %.3e\n", double(g_tLocks), double(g_tUnLocks));
#endif
		return false;
		}

	++m_SeqCount;
	return true;
	}

SeqSource *MakeSeqSource(const string &FileName, DB_SORT SortOrder)
	{
	bool Nucleo;
	FILE_TYPE FileType = GetFileType(FileName, &Nucleo);

	SeqSource *SS = 0;
	if (SortOrder != DBS_None)
		{
		SeqDB *DB = new SeqDB;
		SeqSource *SS2 = MakeSeqSource(FileName, DBS_None);
		DB->FromSS(*SS2);
		DB->Sort(SortOrder);

		SeqDBSeqSource *DBSS = new SeqDBSeqSource;
		DBSS->Init(DB);
		SS = DBSS;
		}
	else
		{
		if (FileType == FT_FASTA)
			{
			FASTASeqSource *FSS = new FASTASeqSource;
			FSS->Open(FileName);
			SS = FSS;
			}
		else if (FileType == FT_FASTQ)
			{
			FASTQSeqSource *FSS = new FASTQSeqSource;
			FSS->Open(FileName);
			SS = FSS;
			}
		else
			Die("Input file format not supported by %s: %s", 
			  CmdToStr(GetCmd()), FileName.c_str());
		}

	asserta(SS != 0);
	return SS;
	}
