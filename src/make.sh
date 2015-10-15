echo "#define ClonalFrameML_GITRevision \"`git describe --tags`\"" > version.h
g++ main.cpp -o ClonalFrameML -I ./ -I ./myutils -I ./coalesce -O3

