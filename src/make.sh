echo -n "#define ClonalFrameML_SVNRevision " > version.h; svn info | grep "Revision: " | cut -f 2 -d' ' >> version.h
g++ main.cpp -o ClonalFrameML -I ./ -I ./myutils -I ./coalesce -O3

