#ifndef ufihit_h
#define ufihit_h

class State1;
class AlignResult;

class UFIHit
	{
public:
	State1 *m_UD;
	unsigned m_DBStartPos;
	bool m_Plus;
	int m_Score;
	string m_Path;
	};

#endif // ufihit_h
