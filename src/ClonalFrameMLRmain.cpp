#include "ClonalFrameMLRmain.h"

void ClonalFrameMLRmain(int* argv_nrow, int* argv_ncol, char** argv_in, int* EXIT_STATUS)
{
	const int argc = *argv_nrow;
	const char** argv = (const char**)argv_in;
	*EXIT_STATUS = main(argc,argv);
	return;
}
