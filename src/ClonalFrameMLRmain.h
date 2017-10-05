#ifndef _CLONALFRAMEML_RMAIN_H_
#define _CLONALFRAMEML_RMAIN_H_
#include <R.h>

extern "C" {
	void ClonalFrameMLRmain(int* argv_nrow, int* argv_ncol, char** argv, int* EXIT_STATUS);
}

int main (const int argc, const char* argv[]);

#endif

