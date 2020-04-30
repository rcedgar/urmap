#ifndef outfiles_h
#define outfiles_h

extern FILE *g_fOut;
extern FILE *g_fSAM;
extern FILE *g_fTab;

void OpenOutputFiles();
void CloseOutputFiles();

#endif // outfiles_h
