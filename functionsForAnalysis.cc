#include "TMath.h"
#include <TAxis.h>
#include <TH1D.h>
#include <TF1.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <iostream>
#include <string>
#include <iomanip> //for input/output manipulators
#include <fstream>
 
#include "functionsForAnalysis.h"

using namespace std;

void myDashes() {

    cout<<"-----------------------------------------------------------------------------------"<<endl; 
    //cout<<"adding this line"<<endl;
}

void myDashesInFile(ofstream &file) {

    file<<"-----------------------------------------------------------------------------------"<<endl; 

}

Double_t myGetPhi(const Double_t x_component, const Double_t y_component) {

  Double_t phi = 0;
  //some useful variables
   Double_t myPi = TMath::Pi(); //just not to call TMath::Pi() every time
   Double_t twoPi = 2*TMath::Pi(); //just not to calculate 2*TMath::Pi() every time

     //while evaluating phi =atan(y/x), we require x!=0. If x=0, then 
     //phi= +/- Pi()/2 with the same sign as the y 
     //for x!=0, since atan(y/x) yields a value in [-Pi()/2, +Pi()/2]
     //we must sum or subtract Pi() depending on the dial

     if (x_component!=0) {

       if ( (x_component<0) && (y_component<0) ) {
         phi = atan(y_component/x_component)-myPi;
       } else if ( (x_component<0) && (y_component>0) ) {
         phi = atan(y_component/x_component)+myPi;
       } else {
         phi = atan(y_component/x_component);
       }

     } else {

      if ( (y_component<0) ) {
	 phi = -myPi/2.;
       } else { 
	 phi = myPi/2.;
      }
     }
     return phi;

}

//this function put the overflows in the last bin of an histogram
void myFillOverflowBin(TH1D *histo, const Double_t value) {
  
  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value>(axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()));
  } else {
    histo->Fill(value);
  }

}

void myFillOverflowBinW(TH1D *histo, const Double_t value, const Double_t weight) {
  
  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value>(axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()),weight);
  } else {
    histo->Fill(value, weight);
  }

}

//calculation of rejection factors
Double_t myRejectionFactor(Int_t selected, Int_t total) {
 
  return (Double_t) total/selected;

}

Double_t myCrystalBall(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t >= -absAlpha)  {
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A/TMath::Power(B-t,n);
  }

}


Double_t my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

void myAddDefaultPackages(FILE *fp, const char* filename) {

  if (fp == NULL) {
    cout<<"Error loading default packages in file \""<<filename<<"\".\nThese will not be added to the file."<<endl;
  } else {
    fprintf(fp,"\\documentclass[11pt,a4paper]{article}\n");
    fprintf(fp,"\\usepackage[T1]{fontenc}\n");
    fprintf(fp,"\\usepackage[utf8]{inputenc}\n");
    fprintf(fp,"\\usepackage[english,italian]{babel}\n");
    fprintf(fp,"\\usepackage{graphicx}\n");
    fprintf(fp,"\\usepackage{wrapfig}\n");
    fprintf(fp,"\\usepackage{sidecap}\n");
    fprintf(fp,"\\usepackage{booktabs}\n");
    fprintf(fp,"\\usepackage{amsmath}\n");
    fprintf(fp,"\\usepackage{amsfonts}\n");
    fprintf(fp,"\\usepackage{amssymb}\n");
    fprintf(fp,"\\usepackage{multirow}\n");
    fprintf(fp,"\\usepackage{subfig}\n");
    fprintf(fp,"\\usepackage{siunitx}\n");
    fprintf(fp,"\n");
  }
  
}

char myAskToSaveFile(const char* filename) {

  char answer = '\0';

  do {
     cout<<"Save event yields on file \""<<filename<<"\"?(y/n) : ";
     cin>>answer;
       if (cin.fail()) {
	 cout<<"Error: cin failed. Data will not be written on \""<<filename<<"\""<<endl;
	 answer = 'n';
       }
       cout<<"answer is \""<<answer<<"\""<<endl;
   } while ( (answer != 'y') && (answer != 'n') );

  return answer;

}

Int_t myGetBin(Double_t x, Double_t *inputArray, Int_t nBins) {
  
// given an array containing a set of edges defining some bins, this function returns the bin which the value x belongs to. 
// e.g.:  array = {1.2, 3.4, 4.5} has 3 elements and defines 2 bins, so that nBins = 2. If x = 3.7, the function returns 1

  if ( (x < *(inputArray)) || (x > *(inputArray + nBins)) ) {
    cout<<"Warning in function myGetBin(): bin not found! End of programme!"<<endl;
    exit(EXIT_FAILURE);
  }

  Int_t correctBin = 0;
  for (Int_t bin = 0; bin < nBins; bin++) {
    if (x > *(inputArray+bin) )  correctBin = bin;  
  } 

  return correctBin;
 
}
