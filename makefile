
ROOTFLAGS=`root-config --cflags --libs` 
CXXDEPFLAGS = -MMD -MP

all: main

#main:  analyzer_base.o analyzer_signal.o analyzer_AODjets.o main.C 
main:  QCD.o main.C 
	g++ --std=c++11 $^ -o QCD.exe ${ROOTFLAGS} ${CXXDEPFLAGS}

%.o: %.C  
	g++ --std=c++11 -c $^ -o $@ ${ROOTFLAGS} ${CXXDEPFLAGS} 

#-include *.d

#g++ -Wno-deprecated main.C -o theplot.exe -I$ROOTSYS/include -L$ROOTSYS/lib -I$CMS_PATH/$SCRAM_ARCH/lcg/roofit/5.32.03-cms9/include -L$CMS_PATH/$SCRAM_ARCH/lcg/roofit/5.32.03-cms9/lib -lRooFit -lRooFitCore `root-config --cflags` `root-config --libs`
