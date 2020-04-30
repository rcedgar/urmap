#include "myutils.h"
#include "cmd.h"

#define IS_TRUE(x)	case CMD_##x: return true;
//#define ____(x)	case CMD_##x: return false;

const char *CmdToStr(CMD Cmd)
	{
	switch (Cmd)
		{
#define A(x)	case CMD_##x: return #x;
#include "cmds.h"
	default:
		asserta(false);
		}
	return 0;
	}

CMD StrToCmd(const char *Str)
	{
#define A(x)	if (!strcmp(Str, #x)) return CMD_##x;
	Die("Invalid cmd '%s'", Str);
	return CMD_none;
	}
