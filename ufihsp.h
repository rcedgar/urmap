#ifndef ufihsp_h
#define ufihsp_h

struct UFIHSP
	{
	uint32 m_StartPosQ;
	uint32 m_StartPosDB;
	uint32 m_Length;
	int m_Score;
	bool m_Plus;
	bool m_Aligned;
	};

#endif //  ufihsp_h
