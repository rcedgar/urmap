#include "myutils.h"
#include "cmd.h"
#include "alpha.h"
#include "outfiles.h"

int main(int argc, char **argv)
	{
//	setbuf(stdout, 0);
	setbuf(stderr, 0);

	MyCmdLine(argc, argv);

	if (!opt(quiet))
		{
		PrintProgramInfo(stderr);
		PrintCopyright(stderr);
		}

	SetLogFileName(opt(log));
	LogProgramInfoAndCmdLine();
	Progress_Init();

	InitAlpha();

	OpenOutputFiles();

	CMD Cmd = GetCmd();
	switch (Cmd)
		{
#define A(x)	case CMD_##x: { void cmd_##x(); cmd_##x(); break; }
#include "cmds.h"
	default:
		asserta(false);
		}
	ProgressSetMsg("Closing");
	CloseOutputFiles();
	CheckUsedOpts(opt_log_used_opts);
	LogElapsedTimeAndRAM();
	Progress_OnExit();
	return 0;
	}
