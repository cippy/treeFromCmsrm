CC = g++ 
CCFLAGS = -Wno-deprecated -g -O3
OBJ_MAIN = tmp/objects/main.o
OBJ_MY_FUNCTIONS = tmp/objects/functionsForAnalysis.o
OBJECTS =  $(OBJ_MAIN) $(OBJ_MY_FUNCTIONS)
ROOTLIB = $(shell $(ROOTSYS)/bin/root-config --cflags --libs)
EXE = tmp/Analysis

Analisys : Analysis.cc $(OBJECTS)
	$(CC) $(CCFLAGS) -o  $(EXE) Analysis.cc $(OBJECTS) $(ROOTLIB)

Analisys.cc : $(OBJECTS)
	$(CC) $(CCFLAGS) -o $(EXE) Analysis.cc $(OBJECTS) $(ROOTLIB)

functionsForAnalysis.o : functionsForAnalysis.h functionsForAnalysis.cc
	$(CC)  -c functionsForAnalysis.cc $(ROOTLIB)

functionsForAnalysis.cc : functionsForAnalysis.h
	$(CC)  -c functionsForAnalysis.cc -o $(OBJ_MY_FUNCTIONS) $(ROOTLIB)

main.o : main.cc
	$(CC) -c main.cc -o $(OBJ_MAIN) $(ROOTLIB)

.PHONY : clean
clean: 
	-rm $(OBJECTS)


