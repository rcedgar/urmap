#ifndef bitvec_h
#define bitvec_h

class BitVec
	{
public:
	byte *m_Vec;
	uint64 m_BitCount;

public:
	BitVec();
	virtual ~BitVec();
	void Alloc(uint64 BitCount);
	void Free();
	bool GetBit(uint64 n) const;
	void SetBit(uint64 n);
	void ClearBit(uint64 n);
	void Zero();
	};

static const uint32 BV_MAGIC = MAGIC('B', 'V', '1', '0');

#endif // bitvec_h
