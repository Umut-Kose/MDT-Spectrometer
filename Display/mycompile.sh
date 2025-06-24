rm -rf eve_display.dSYM eve_display MyDisplayDict_rdict.pcm MyDisplayDict.cxx

rootcling -f MyDisplayDict.cxx -c MyDisplay.h LinkDef.h

g++ -std=c++17 -g -O0 -Wall -Wextra \
    MyDisplayDict.cxx MyDisplay.cpp main.cpp -o eve_display \
    `root-config --cflags --glibs` -lEve -lRGL -lGeom -lGed -lGLEW
