#include "myutils.h"
#include "bitvec.h"

BitVec::BitVec()
	{
	m_Vec = 0;
	m_BitCount = 0;
	}

BitVec::~BitVec()
	{
	Free();
	}

void BitVec::Alloc(uint64 BitCount)
	{
	asserta(m_Vec == 0);
	uint64 Bytes = BitCount/8 + 1;
	m_Vec = myalloc64(byte, Bytes);
	memset_zero(m_Vec, Bytes);
	m_BitCount = BitCount;
	}

void BitVec::Zero()
	{
	uint64 Bytes = m_BitCount/8 + 1;
	memset_zero(m_Vec, Bytes);
	}

void BitVec::Free()
	{
	if (m_Vec != 0)
		{
		myfree(m_Vec);
		m_Vec = 0;
		m_BitCount = 0;
		}
	}

bool BitVec::GetBit(uint64 n) const
	{
	byte Byte = m_Vec[n/8];
	return (Byte & (1 << n%8)) != 0;
	}

void BitVec::SetBit(uint64 n)
	{
	m_Vec[n/8] |= (1 << (n%8));
	}

void BitVec::ClearBit(uint64 n)
	{
	m_Vec[n/8] &= ~(1 << (n%8));
	}
