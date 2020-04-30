#ifndef lockobj_h
#define lockobj_h

#define L(x)	extern omp_lock_t g_Lock##x;
#include "lockobjs.h"

#define LOCKABLE(x)  \
public:  \
	static void LOCK_CLASS() { omp_set_lock(&g_Lock##x); }  \
	static void UNLOCK_CLASS() { omp_unset_lock(&g_Lock##x); }

#endif // lockobj_h
