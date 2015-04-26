#ifndef _FUNCTIONS_FOR_ANALYSIS
#define _FUNCTIONS_FOR_ANALYSIS

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

void myDashes(ostream &);

Double_t myGetPhi(Double_t, Double_t);

void myFillOverflowBin(TH1D*, const Double_t);

void myFillOverflowBinW(TH1D*, const Double_t, const Double_t);

void myDrawOverflow(const TH1D *, const char*); 
void myDrawOverflow(const TH1D *, const char*, const Int_t); 
TH1D * myOverflowInLastBin(const TH1D *); 
void myOverflowInLastBin2(TH1D *, const TH1D *); 
void myAddOverflowInLastBin(TH1D *);

Double_t myRejectionFactor(Int_t, Int_t);

Double_t myCrystalBall(double* , double* );

Double_t my2sideCrystalBall(double* , double* );

Double_t myResolutionFunction(double* , double* );

Double_t myResolutionFunctionNoB(double* , double* );

void myAddDefaultPackages(FILE *, const char* );

char myAskToSaveFile(const char*);

Int_t myGetBin(Double_t, Double_t*, Int_t);

void myPrintEventYields(ostream &, const Double_t, const Int_t, const Double_t *);

#endif
