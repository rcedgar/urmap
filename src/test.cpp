#include "myutils.h"
#include "tenx.h"

void cmd_test()
	{
	unsigned n = count_t(~1);
	ProgressLog("count_t(~1) = %u\n", n);
	}
