#ifndef seqsource_h
#define seqsource_h

#include "lockobj.h"
#include "filetype.h"

class SeqInfo;
class ObjMgr;

enum DB_SORT
	{
	DBS_None,
	DBS_Size,
	DBS_Length,
	DBS_Label,
	};

class SeqSource
	{
	LOCKABLE(SeqSource)

public:
	bool m_DoGetLock;

protected:
	SeqInfo *m_SI;
	unsigned m_SeqCount;

public:
	SeqSource();
	virtual ~SeqSource();

public:
	virtual bool GetIsNucleo() = 0;
	virtual const char *GetFileNameC() const = 0;

public:
	virtual bool GetNextLo(SeqInfo *SI) = 0;

public:
	bool GetNext(SeqInfo *SI);
	};

SeqSource *MakeSeqSource(const string &FileName, DB_SORT SortOrder);

#endif // seqsource_h
