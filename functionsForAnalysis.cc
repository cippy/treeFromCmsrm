#include "TMath.h"
#include <TAxis.h>
#include <TH1D.h>
#include <TF1.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip> //for input/output manipulators
#include <fstream>
 
#include "functionsForAnalysis.h"

using namespace std;

void myDashes(ostream &myOutStream = cout) {

    myOutStream<<"-----------------------------------------------------------------------------------"<<endl; 

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

void myFillOverflowBin(TH1D *histo, const Double_t value) {
  
  //this function put the overflows in the last bin of an histogram

  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value>(axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()));
  } else {
    histo->Fill(value);
  }

}

void myFillOverflowBinW(TH1D *histo, const Double_t value, const Double_t weight) {

  // filling last bin with overflow for histograms with weights

  TAxis *axis = 0;
  axis = histo->GetXaxis();
  if(value > (axis->GetBinUpEdge(axis->GetLast()))) {
    histo->Fill(axis->GetBinCenter(axis->GetLast()),weight);
  } else {
    histo->Fill(value, weight);
  }

}

void myDrawOverflow(const TH1D *oldHisto, const char* myOptions) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D *htemp = new TH1D("htemp","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp->SetLineColor(oldHisto->GetLineColor());
  htemp->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  htemp->Draw(myOptions);

}

void myDrawOverflow(const TH1D *oldHisto, const char* myOptions, const Int_t mySetStats = 1) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D *htemp = new TH1D("htemp","",newNBins,newLowerEdge,newUpperEdge);
  // TH1D histo("histo","",newNBins,newLowerEdge,newUpperEdge);
  // TH1D *htemp = &histo;

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp->SetLineColor(oldHisto->GetLineColor());
  htemp->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  if( !mySetStats ) htemp->SetStats(kFALSE);
  htemp->Draw(myOptions);

}

TH1D *myOverflowInLastBin(const TH1D *oldHisto) {

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  static TH1D *htemp1 = new TH1D("htemp1","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp1->SetBinContent(i,oldHisto->GetBinContent(i));
    htemp1->SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp1->SetBinContent(0,oldHisto->GetBinContent(0));
  htemp1->SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp1->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp1->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp1->SetLineColor(oldHisto->GetLineColor());
  htemp1->SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp1->SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  return htemp1;

}

void myOverflowInLastBin2(TH1D *newHisto, const TH1D *oldHisto) {

  // this function, when called with a histo defined as "TH1D *histo = 0" as newHisto, produces segmentation fault!!!

  // this function create a new histogram with one more bin than oldHisto. This additional bin will contain the overflow of oldHisto. 
  // The width of this additional bin is set equal to the width of the last bin of oldHisto. This is not relevant for histograms where the
  // width is always the same for each bin, but care should be taken for histograms whose bin's width depends on the bit
  // By definition, the overflow bin of newHisto is empty
  // The strategy is creating a temporary array htemp which, before being deleted (since it's created with "new" it' better to delete it at the end
  // of the function), will be copied into newHisto

  TH1::SetDefaultSumw2();            // all the following histograms will automatically call TH1::Sumw2() 
  TH1::StatOverflows();                 // enable use of underflows and overflows for statistics computation 

  TAxis *axis = 0;
  axis = oldHisto->GetXaxis();
  Double_t overflowBinWidth = axis->GetBinWidth(axis->GetLast());
  Double_t newUpperEdge = axis->GetBinUpEdge(axis->GetLast()) + overflowBinWidth;
  Double_t newLowerEdge = axis->GetBinLowEdge(axis->GetFirst());
  Int_t newNBins = oldHisto->GetNbinsX() + 1;   //  +1 because I want to add an additional bin for the overflow

  // some checks to keep things under control 
  // cout<<"inside myOverflowInLastBin()"<<endl;
  // cout<<"newNBins = "<<newNBins<<endl;
  // cout<<"newUpperEdge = "<<newUpperEdge<<endl;
  // cout<<"oldHisto->GetBinContent(newNBins) = "<<oldHisto->GetBinContent(newNBins)<<endl;

  TH1D htemp2("htemp2","",newNBins,newLowerEdge,newUpperEdge);

  // setting bin content and error
  for (Int_t i = 1; i <= newNBins; i++) {
    htemp2.SetBinContent(i,oldHisto->GetBinContent(i));
    htemp2.SetBinError(i,oldHisto->GetBinError(i));
  }
  // filling underflow bin 
  htemp2.SetBinContent(0,oldHisto->GetBinContent(0));
  htemp2.SetBinError(0,oldHisto->GetBinError(0));
  // setting line color and axis titles
  htemp2.GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
  htemp2.GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
  htemp2.SetLineColor(oldHisto->GetLineColor());
  htemp2.SetFillColor(oldHisto->GetFillColor());
  // setting number of entries to effective entries (same as entries for unweighted histograms) + content 
  // of overflow bin (which is not included by Get(Effective)Entries())
  htemp2.SetEntries(oldHisto->GetBinContent(newNBins) + oldHisto->GetEffectiveEntries());

  *newHisto = htemp2;

}

void myAddOverflowInLastBin(TH1D *h) {

  // to avoid problems regarding memory leak for not deleting htemp in previous function, I sum directly the content of overflow bin in last bin

  Int_t lastBinNumber = h->GetNbinsX();
  Int_t overflowBinNumber = 1 + lastBinNumber;
  Double_t lastBinContent = h->GetBinContent(lastBinNumber);
  Double_t overflowBinContent = h->GetBinContent(overflowBinNumber);
  Double_t lastBinError = h->GetBinError(lastBinNumber);
  Double_t overflowBinError = h->GetBinError(overflowBinNumber);

  // add content of overflow bin in last bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)
  h->SetBinContent(lastBinNumber, lastBinContent + overflowBinContent);
  h->SetBinError(lastBinNumber, sqrt(lastBinError * lastBinError + overflowBinError * overflowBinError));

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

Double_t myResolutionFunction(double* x, double* par) {

  // resolution function dE/E is given by the sum in quadrature of 3 terms: a / sqrt(E) ; b / E ; c
  // thus, sigma(E) is obtained by multiplying dE / E by E  -->  sigma(E) = a*sqrt(E) + b + c*E where + is the sum in quadrature
  // a, b and c are the stochastic, noise and constant term respectively

  Double_t E = x[0];      // energy, pT or whatever
  Double_t a = par[0];
  Double_t b = par[1];
  Double_t c = par[2];

  return sqrt( a*a * E + b*b + c*c * E*E );

}

Double_t myResolutionFunctionNoB(double* x, double* par) {

  // resolution function dE/E is given by the sum in quadrature of 3 terms: a / sqrt(E) ; b / E ; c
  // thus, sigma(E) is obtained by multiplying dE / E by E  -->  sigma(E) = a*sqrt(E) + b + c*E where + is the sum in quadrature
  // a, b and c are the stochastic, noise and constant term respectively

  // the term with b could be neglected at big values of E

  Double_t E = x[0];      // energy, pT or whatever
  Double_t a = par[0];
  Double_t c = par[1];

  return sqrt( a*a * E + c*c * E*E );

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


void myPrintEventYields(ostream & myOutStream, const Double_t lumi, const Int_t cutSteps, const Double_t * eventsInStep) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        EVENT YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"trigger cuts disabled"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<eventsInStep[0]<<endl;
  myOutStream<<"n:               number of events after i-th cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/n(1)"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<setw(3)<<"cut"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;
  for (Int_t i = 1; i <= cutSteps; i++) {
    if (i == 1) {	  
      myOutStream<<setw(3)<<i<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<1.0<<setw(8)<<1.0<<endl;     
    } else {
      myOutStream<<setw(3)<<i<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<eventsInStep[i]/eventsInStep[1]<<
	setw(8)<<eventsInStep[i]/eventsInStep[i-1]<<endl;  
    }
  }
  myOutStream<<"--------------------------------------------------------"<<endl;
  myOutStream<<endl;

}
