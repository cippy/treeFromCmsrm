#ifndef _FUNCTIONS_FOR_ANALYSIS
#define _FUNCTIONS_FOR_ANALYSIS

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

void myDashes();

void myDashesInFile(ofstream &file);

Double_t myGetPhi(Double_t, Double_t);

void myFillOverflowBin(TH1D*, const Double_t);

void myFillOverflowBinW(TH1D*, const Double_t, const Double_t);

Double_t myRejectionFactor(Int_t, Int_t);

Double_t myCrystalBall(double* , double* );

Double_t my2sideCrystalBall(double* , double* );

void myAddDefaultPackages(FILE *, const char* );

char myAskToSaveFile(const char*);

Int_t myGetBin(Double_t, Double_t*, Int_t);

#endif
