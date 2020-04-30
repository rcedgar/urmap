#include "myutils.h"

#define TEST	0

static uint64 g_Primes[] =
	{
#include "primes.h"
	};
static unsigned g_PrimeCount = sizeof(g_Primes)/sizeof(g_Primes[0]);

uint64 GetPrime(uint64 n)
	{
	for (unsigned i = 0; i < g_PrimeCount; ++i)
		{
		uint64 Prime = g_Primes[i];
		if (Prime >= n)
			return Prime;
		}
	Die("GetPrime(%.3g) overflow", double(n));
	return 0;
	}
