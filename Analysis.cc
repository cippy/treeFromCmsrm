#define Analysis_cxx
#include "Analysis.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TVirtualFitter.h>

//my headers
#include "functionsForAnalysis.h"

using namespace std;

//path of directory where saved canvases are stored
#define MY_PDF_PATH "./pdfsFromAnalysis/"
#define MY_TEX_DIR "texsFromAnalysis/"


#define LUMI 5.                                 //this number refers to the integrated luminosity used in multiple of 1fb^-1
#define K_FACTOR 1.27                    //the weight wgt happens to be 1.27 times smaller than what it should be, it's just a k-factor 
#define NVTXS  23                           // # of points for study of mumetX and mumetY vs # of reconstructed vertices nvtx
#define FIRST_NVTX  6                    // starting number of vertices for met study                  
#define NCUTS 11                             //number of cuts for the analysis (see CUTS FLOW table below)
#define ENABLE_TRIGGER_CUTS false  //when "true", trigger cuts are enabled for analysis: this means that the initial number of entries is the number of entries surviving trigger cuts. These cuts are described in the CUT FLOW table below
//with the following lines, when the ENABLE_TRIGGER_CUTS directive is set to "true", trigger thresholds are defined properly (one can change them to desired values), otherwise, if it's set to "false", thresholds are set to -1, which means no thresholds at all (positive quantities by definition). 
#if ENABLE_TRIGGER_CUTS
#define JET_TRIGGER 80                        //trigger threshold for centraljet
#define PFMET_TRIGGER 130                 //trigger threshold for particle flow met (pfmet)
#define MHT_TRIGGER 130                     //trigger threshold for mht (met computed from sum of jets)
#else
#define JET_TRIGGER (-1)          
#define PFMET_TRIGGER (-1)   
#define MHT_TRIGGER (-1)
#endif
//thesholds for cuts
#define NJETS 2
#define J1PT 110
#define J1ETA 2.4
#define J1NHEF 0.7
#define J1CHEF 0.2
#define J1NEEF 0.7
#define J2PT 30
#define J2ETA 4.5
#define J2NHEF 0.7
#define J2NEEF 0.9
#define J1J2DPHI 2.4
#define NMUONS 0
#define NELECTRS 0
#define NTAUS 0
#define MUMET 200

//----------------------------------------------------------------------------------
//implementation of methods of class cut declared at the end of Analysis.h

//static data members
Int_t cut::nCuts_ = 0;
vector<cut*> cut::listOfCuts; 

void cut::printCutFlow(ostream & myOutStream, const Int_t cutSteps, const UInt_t *singleCutMask) {

  //this function prints the cut flow on the appropriate ofstream (can be a file or cout). Since it's a member function of class cut, it only needs to get the number of cut steps. For now it needs the array of masks which is not a class data member yet

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*          CUTS FLOW          *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<"Printing list of cuts applied at each step"<<endl;
  for (Int_t i = 0; i < cutSteps; i++) {
    myOutStream<<"-----------------------------------"<<endl;
    myOutStream<<setw(2)<<(i+1)<<endl;
    for (Int_t j = 0; j < cut::getNCuts(); j++) {
      if ((singleCutMask[i] >> j) & 1) {
	myOutStream<<"      ";
	cut::listOfCuts[j]->print(myOutStream,1); 
      }
    }
  }
  myOutStream<<"-----------------------------------"<<endl;  

}

void cut::printActiveCuts(ostream & myOutStream) {

  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  myOutStream<<"Printing list of activated cuts"<<endl;
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  for (Int_t i = 0; i < cut::getNCuts(); i++ ) {
    if ( cut::listOfCuts[i]->isActive() ) {
      cut::listOfCuts[i]->printAllInfo(myOutStream);
    }
  }
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;

}

cut::cut(Bool_t flag, const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const string comment) {
  flag_ = flag;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = threshold;
  comment_ = comment;
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(Bool_t flag, const char *cut_name, const char *var_name, const char *condition, const Double_t threshold) {
  flag_ = flag;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = threshold;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(Bool_t flag, const char *cut_name, const char *var_name, const char *condition) {
  flag_ =flag;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = 0.;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(Bool_t flag, const char *cut_name, const char *var_name) {
  flag_ =flag;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = " ";
  thr_ = 0.;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::~cut() {
  cout<<"~cut() called for "<<cut_<<endl;
  nCuts_--;
}

void cut::printAllInfo(ostream & myOutStream) const {

  myOutStream<<setw(18)<<left<<cut_<<": id = "<<setw(2)<<right<<id_<<" | "<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" ";
  myOutStream<<setw(4)<<right<<thr_<<"   "<<left<<comment_<<endl; 

 }

void cut::print(ostream & myOutStream = cout, Bool_t addComment = 0) const {

  //print cut definition on the right ofstream (cout is default), if addComment = 1, cut  comment is printed if any
  if (addComment) {
    myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<thr_<<"   "<<left<<comment_<<endl; 
  } else {
    myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<thr_<<endl; 
  }

 }

Bool_t cut::isPassed(Double_t input) {

  if (cond_ == "<") return (input < thr_) ? true : false;
  else if (cond_== ">") return (input > thr_) ? true : false;
  else if (cond_== "=") return (input == thr_) ? true : false;
  else if (cond_ == "<=") return (input <= thr_) ? true : false;
  else if (cond_ == ">=") return (input >= thr_) ? true : false;
  else {
    cout<<"Error in cut::isPassed()! No condition fulfilled."<<endl;
    cout<<"End of programme"<<endl;
    exit(EXIT_FAILURE);
  }

}


//----------------------------------------------------------------------------------

void Analysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Analysis.C
//      Root > Analysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  

   fChain->SetBranchStatus("nvtx",1);

   fChain->SetBranchStatus("njets",1);         // # of jets with pt>30GeV
   fChain->SetBranchStatus("nmuons",1);    //for muon veto, # of muons passing loose cuts
   fChain->SetBranchStatus("nelectronsnew",1);    //for electron veto, # of electrons passing loose cuts
   fChain->SetBranchStatus("ntaus",1);    //for tau veto, # of taus (I guess they are taus passing veto but Adish didn't specify 

   fChain->SetBranchStatus("wgt",1);                                //weight for the event: it has to be multiplied by a k-factor which is equal to 1.27 (look at the "define" above)
   fChain->SetBranchStatus("signaljetCHfrac",1);              //charged H fraction
   fChain->SetBranchStatus("signaljetNHfrac",1);             //neutral H fraction
   fChain->SetBranchStatus("signaljetEMfrac",1);             //neutral EM fraction
   fChain->SetBranchStatus("secondjetNHfrac",1);           //neutral H fraction
   fChain->SetBranchStatus("secondjetEMfrac",1);            //neutral EM fraction
   
   fChain->SetBranchStatus("signaljetpt",1);
   fChain->SetBranchStatus("signaljetphi",1);
   fChain->SetBranchStatus("signaljeteta",1);
   fChain->SetBranchStatus("secondjetpt",1);
   fChain->SetBranchStatus("secondjetphi",1);
   fChain->SetBranchStatus("secondjeteta",1);
   fChain->SetBranchStatus("thirdjetpt",1);
   fChain->SetBranchStatus("thirdjetphi",1);
   fChain->SetBranchStatus("thirdjeteta",1);

   fChain->SetBranchStatus("jetjetdphi",1); //delta phi in [-PI,PI] between first and second leading jet

   fChain->SetBranchStatus("l1pt",1);
   fChain->SetBranchStatus("l1phi",1);
   fChain->SetBranchStatus("l1eta",1);   
   fChain->SetBranchStatus("l2pt",1);
   fChain->SetBranchStatus("l2phi",1);
   fChain->SetBranchStatus("l2eta",1);   

   fChain->SetBranchStatus("mht",1);  // necessary for trigger cuts
   fChain->SetBranchStatus("pfmet",1);  // necessary for trigger cuts
   fChain->SetBranchStatus("mumet",1); 
   fChain->SetBranchStatus("mumetphi",1); 
   fChain->SetBranchStatus("t1pfmet",1); 
   fChain->SetBranchStatus("t1pfmetphi",1); 

   fChain->SetBranchStatus("wzmass",1);   
   fChain->SetBranchStatus("wzpt",1);  
   fChain->SetBranchStatus("wzeta",1);  
   fChain->SetBranchStatus("wzphi",1);   

 /*

*******************************
*       GENERAL COMMENTS      *
*******************************

In the following, cuts due to trigger are not considered, that is, I assumed that the ntuple was made out of data passing the trigger selection: this assumption is corroborated  by the evidence of a sharp cut at 200 GeV in the mumet variable stored in the tree. However their contribution can reduce the initial number of entries and it would be useful to have the possibility to take them into account. 

---------------- TRIGGER CUTS -------------------  taken from Nadir's slides-> https://indico.cern.ch/event/368971/contribution/6/material/slides/0.pdf

1)signaljetpt>140 GeV & pfmet>100 GeV & mht>140 GeV       (variables are stored in tree.root)
2)signaljetpt>80 GeV & pfmet>130 GeV & mht>130 GeV  
Use one or the other trigger, the second should provide bigger acceptance (mainly more events)

--------------  CUTS FLOW  -----------------------    the number of cuts is equal to NCUTS defined above. 
cuts are cumulative ( e.g.: 3) also includes 1) and 2) ). 
Note that the cut flow can change depending on the analysis, the following is just a guideline: refer to globalCutMask[] 

1) MET > 250 GeV (I use mumet)
2) signaljetpt>110 GeV; |signaljeteta|<2.5; Noise cleaning on first jet ( NHEF<0.7; CHEF>0.2; NEEF<0.7 )
3) secondjetpt>30 GeV; |secondjeteta|<4.5,  NHEF<0.7; NEEF<0.9 (only for events with 2 jets)
4) dphi(jet1,jet2)<2.5;
5) njets<= 2 (njets is the # of jets with pt>30)
6) muon veto: reject if pt>10GeV and |eta|<2.4
7) electron veto: reject if pt>10GeV and |eta|<2.5
8) tau veto: reject if pt>20GeV and |eta|<2.3 (in this analysis I reject events with taus because i don't have any variable for taus but their number)
9) MET > 300
10) MET > 400
11) MET > 500

--------------------------------------------------------------
 
Cuts are applied using a mask. I built a class cut (see bottom of Analysis.h) which will be represented by a specific bit in a mask. Each bit says if the selection based on its related variable was passed or not.
In these way, there will be a mask for each event recording which selections will be passed in that event. Then, I will be free to fill histograms with the desired variable according to a specific step of the cut flow. As an example, if my mask for variables A, B, C, D is M=1011 (in binary representation) for a given event, than all selections except for C (starts from the rightest bit) were passed. If I'm interested in the histogram of variable x after cuts A and D, the corresponding mask is m=1001.
I will fill the histogram if ((M & m) == m). that is to say, only if bit equal to 1 in m are also 1 in M (independently of the fact that M as other bits equal to 1).
The biggest advantage is that I don't need nested "if" to implement the selection: I' m just labeling the event, and the labels merely categorize events depending on the selections that would be passed by those events.

The following is the list of variables, reflecting in the bits numbering (if Mask is UInt_t, I can have only 32 variables). I added variables that seemed reasonable to be used in the future (names taken from Adish's file). Actually the masks will be updated each time the programme run, so the numbering is more a guideline than an actual rule (it will depend on the order that cuts are declared, which is absolutely arbitrary).

1)  njets
2)  signaljetpt
3)  signaljeteta
4)  signaljetphi
5)  signaljetCHfrac
6)  signaljetNHfrac
7)  signaljetEMfrac
8)  signaljetCEMfrac
9)    secondjetpt
10)  secondjeteta
11)  secondjetphi
12)  secondjetCHfrac
13)  secondjetNHfrac
14)  secondjetEMfrac
15)  secondjetCEMfrac
16)  jetjetdphi
17)  mumet
18)  nmuons
19)  nelectronsnew
20)  ntaus
21)  thirdjetpt
22)  thirdjeteta
23)  thirdjetphi
24) mumet300
25) mumet400
26) mumet500

    */

   //I use the following convention: all classes of type cut are named with a final C (C stands for "Cut") 
   cut njetsC(true,"njetsC","njets","<=",NJETS);
   cut jet1ptC(true,"jet1ptC","signaljetpt",">",J1PT);
   cut jet1etaC(true,"jet1etaC","|signaljeteta|","<",J1ETA);
   //cut jet1phiC(false,"jet1phiC","signaljetphi");
   cut jet1CHfracC(true,"jet1CHfracC","signaljetCHfrac",">",J1CHEF);
   cut jet1NHfracC(true,"jet1NHfracC","signaljetNHfrac","<",J1NHEF);
   cut jet1EMfracC(true,"jet1EMfracC","signaljetEMfrac","<",J1NEEF);
   //cut jet1CEMfracC(false,"jet1CEMfracC","signaljetCEMfrac"); 
   //cut jet2ptC(true,"jet2ptC","secondjetpt",">",J2PT,"only if njets = 2");
   cut jet2etaC(true,"jet2etaC","|secondjeteta|","<",J2ETA,"only if njets = 2");
   //cut jet2phiC(false,"jet2phiC","secondjetphi");
   //cut jet2CHfracC(false,"jet2CHfracC","secondjetCHfrac");
   cut jet2NHfracC(true,"jet2NHfracC","secondjetNHfrac","<",J2NHEF,"only if njets = 2");
   cut jet2EMfracC(true,"jet2EMfracC","secondjetEMfrac","<",J2NEEF,"only if njets = 2");
   //cut jet2CEMfracC(false,"jet2CEMfracC","secondjetCEMfrac"); 
   cut jet1jet2dphiC(true,"jet1jet2dphiC","|jetjetdphi|","<",J1J2DPHI,"only if njets = 2");
   //cut mumetC(true,"mumetC","mumet",">",MUMET);
   cut nmuonsC(true,"nmuonsC","nmuons","=",NMUONS);
   cut nelectronsnewC(true,"nelectronsnewC","nelectronsnew","=",NELECTRS);
   cut ntausC(true,"ntausC","ntaus","=",NTAUS);
   cut mumet250C(true,"mumet250C","mumet",">",250);
   cut mumet300C(true,"mumet300C","mumet",">",300);
   //cut mumet350C(true,"mumet350C","mumet",">",350);
   cut mumet400C(true,"mumet400C","mumet",">",400);
   //cut mumet450C(true,"mumet450C","mumet",">",450);
   cut mumet500C(true,"mumet500C","mumet",">",500);
   //cut mumet550C(true,"mumet550C","mumet",">",550);

   if (cut::getNCuts() > 8*sizeof(UInt_t)) {
     cout<<"Warning: not enough bits in the mask to accomodate all "<<cut::getNCuts()<<" cuts (max is "<<8*sizeof(UInt_t)<<").\n End of programme."<<endl;
     exit(EXIT_FAILURE);
   }

   cut::printActiveCuts(cout);

   //building masks for our analysis. If a different cut flow is needed, changes must be done here
   UInt_t singleCutMask[NCUTS] = {};   //mask with cuts in a specific step
   singleCutMask[0] = mumet250C.get2ToId();
   singleCutMask[1] = njetsC.get2ToId();
   singleCutMask[2] = jet1CHfracC.get2ToId() + jet1NHfracC.get2ToId() + jet1EMfracC.get2ToId();
   singleCutMask[3] = jet1ptC.get2ToId() + jet1etaC.get2ToId();
   singleCutMask[4] = jet1jet2dphiC.get2ToId() + /*jet2ptC.get2ToId() +*/ jet2etaC.get2ToId() + jet2NHfracC.get2ToId() + jet2EMfracC.get2ToId();
   singleCutMask[5] = nmuonsC.get2ToId();
   singleCutMask[6] = nelectronsnewC.get2ToId();
   singleCutMask[7] = ntausC.get2ToId();
   singleCutMask[8] = mumet300C.get2ToId();
   singleCutMask[9] = mumet400C.get2ToId();
   singleCutMask[10] = mumet500C.get2ToId();

   // singleCutMask[0] = mumetC.get2ToId();
   // singleCutMask[1] = njetsC.get2ToId();
   // singleCutMask[2] = jet1ptC.get2ToId() + jet1etaC.get2ToId() + jet1CHfracC.get2ToId() + jet1NHfracC.get2ToId() + jet1EMfracC.get2ToId();
   // singleCutMask[3] = jet2ptC.get2ToId() + jet2etaC.get2ToId() + jet2NHfracC.get2ToId() + jet2EMfracC.get2ToId();
   // singleCutMask[4] = jet1jet2dphiC.get2ToId(); 
   // singleCutMask[5] = nmuonsC.get2ToId();
   // singleCutMask[6] = nelectronsnewC.get2ToId();
   // singleCutMask[7] = ntausC.get2ToId();
   // singleCutMask[8] = mumet300C.get2ToId();
   // singleCutMask[9] = mumet400C.get2ToId();
   // singleCutMask[10] = mumet500C.get2ToId();

   cut::printCutFlow(cout,NCUTS,singleCutMask);

   UInt_t globalCutMask[NCUTS] ;   //mask adding to a specific step all previous cuts 
   for (Int_t i = 0; i < NCUTS; i++ ) {
     globalCutMask[i] = 0.0;
   }
   //the following loop' s aim is to have incapsulated conditions (globalMask of step i-th must include all previous steps)  
   for (Int_t i = 1; i < NCUTS; i++) {
     if (i == 0) globalCutMask[0] = singleCutMask[0];
     else globalCutMask[i] = singleCutMask[i] + globalCutMask[i-1];
   }
   cout<<endl;

   Double_t nWeightedEvents[NCUTS+1];  // number of weighted events (it's a Double_t because weights are non-integer numbers)
   for (Int_t i = 0; i <= NCUTS; i++ ) {
     nWeightedEvents[i] = 0.0;
   }

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   TH1D *Hjet1ptTriggerC = new TH1D("Hjet1ptTriggerC","",30,0,1500);  //this is jet1 pt after trigger cuts (alternative to Hjet1pt[0] defined just after it)
   
   TH1D *Hjet1pt[NCUTS+1];     
   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};                     //Hjet1pt[0] is the histogram without any cut
   for (Int_t i=0; i<=NCUTS; i++) {
     Hjet1pt[i] = new TH1D(Form("Hjet1pt[%i]",i),"",30,0,1500);
     Hjet1pt[i]->SetLineColor(Hcolor[i]);
     Hjet1pt[i]->SetFillColor(Hcolor[i]);
   } 
   //Int_t ncut; //index for the cut considered: it goes from 0 (no cuts) to NCUTS (see the CUTS FLOW table above) and changes in the loop

   TH1D *Hjet1etanoC = new TH1D("Hjet1etanoC","",100,-5.,5.);
   TH1D *Hjet1etaallCMet250 = new TH1D("Hjet1etaallCMet250","",50,-2.5,2.5);
   TH1D *Hjet1etaallCMet500 = new TH1D("Hjet1etaallCMet500","",50,-2.5,2.5); //events in this histogram must have passed the cut |eta|<J1ETA
   TH1D *Hjet2ptnoC = new TH1D("Hjet2ptnoC","",30,0,1500);
   TH1D *Hjet2ptallCMet250 = new TH1D("Hjet2ptallCMet250","",30,0,1500);
   TH1D *Hjet2ptallCMet300 = new TH1D("Hjet2ptallCMet300","",30,0,1500);
   TH1D *Hjet2ptallCMet400 = new TH1D("Hjet2ptallCMet400","",30,0,1500);
   TH1D *Hjet2ptallCMet500 = new TH1D("Hjet2ptallCMet500","",30,0,1500);
   TH1D *Hjet2etanoC = new TH1D("Hjet2etanoC","",100,-5.,5.);
   TH1D *Hjet2etaallCMet250 = new TH1D("Hjet2etaallCMet250","",100,-5.,5.);
   TH1D *Hjet2etaallCMet500 = new TH1D("Hjet2etaallCMet500","",100,-5.,5.); //events in this histogram must have passed the cut |eta|<J2ETA
   TH1D *Hjet1jet2dphinoC = new TH1D("Hjet1jet2dphinoC","",80,-4.,4.);
   TH1D *Hjet1jet2dphiallCMet250 = new TH1D("Hjet1jet2dphiallCMet250","",80,-4.,4.);
   TH1D *Hjet1jet2dphiallCMet500 = new TH1D("Hjet1jet2dphiallCMet500","",80,-4.,4.);
   TH1D *HnjetsnoC = new TH1D("HnjetsnoC","",11,-0.5,10.5);
   TH1D *HnjetsallCMet250 = new TH1D("HnjetsallCMet250","",11,-0.5,10.5);
   TH1D *HnjetsallCMet500 = new TH1D("HnjetsallCMet500","",11,-0.5,10.5);
   TH1D *HmumetnoC = new TH1D("HmumetnoC","",30,0,1500);
   TH1D *HmumetallCMet250 = new TH1D("HmumetallCMet250","",30,0,1500);   
   //when dealing with leptons, allC refers to all cuts including those on leptons for background estimation
   TH1D *HlepptnoC = new TH1D("HlepptnoC","",150,0,1500);
   TH1D *HlepptallCnoIso = new TH1D("HlepptallCnoIso","",150,0,1500);      
   TH1D *HlepptallC = new TH1D("HlepptallC","",150,0,1500);
   TH1D *HlepmassnoC = new TH1D("HlepmassnoC","",350,0.,700);
   TH1D *HlepmassallC = new TH1D("HlepmassallC","",30,60,120);
   TH1D *HZptnoC = new TH1D("HZptnoC","",150,0,1500);
   TH1D *HZptallC = new TH1D("HZptallC","",150,0,1500);
   TH1D *HZmassnoC = new TH1D("HZmassnoC","",350,0.,700);
   TH1D *HZmassallC = new TH1D("HZmassallC","",30,60,120);

   // TH1D *HtestNevents = new TH1D("HtestNevents","",1,0,1);
   // TH1D *HtestNeventsBis = new TH1D("HtestNeventsBis","",1,0,6000);

   TH1D *HmumetOrtZvsNvtx[NVTXS];
   TH1D *HmumetParZvsNvtx[NVTXS];
   TH1D *HmumetX[NVTXS];
   TH1D *HmumetY[NVTXS]; 
   for (Int_t i = 0; i < NVTXS; i++) {
     Int_t color = i+1;
     HmumetOrtZvsNvtx[i] = new TH1D(Form("HmumetOrtZvsNvtx[%i]",i),"",80,-200,200);  // 5 GeV bins (was 10 before)
     //HmumetOrtZvsNvtx[i]->SetLineColor(i);
     HmumetParZvsNvtx[i] = new TH1D(Form("HmumetParZvsNvtx[%i]",i),"",80,-200,200);  // 5 GeV bins (was 10 before)
     //HmumetParZvsNvtx[i]->SetLineColor(i);
     HmumetX[i] = new TH1D(Form("HmumetX[%i]",i),"",500,-1250,1250); 
     HmumetY[i] = new TH1D(Form("HmumetY[%i]",i),"",500,-1250,1250); 
   }
   //in a older version I used t1pfmet as well to compare with mumet and (as expected) there was no difference except at first and last bin
   // TH1D *Ht1pfmetOrtZvsNvtx[NVTXS];
   // TH1D *Ht1pfmetParZvsNvtx[NVTXS];
   // for (Int_t i = 0; i < NVTXS; i++) {
   //   Ht1pfmetOrtZvsNvtx[i] = new TH1D(Form("Ht1pfmetOrtZvsNvtx[%i]",i),"",40,-200,200);  // 10 GeV bins
   //   Ht1pfmetOrtZvsNvtx[i]->SetLineColor(i);
   //   Ht1pfmetParZvsNvtx[i] = new TH1D(Form("Ht1pfmetParZvsNvtx[%i]",i),"",40,-200,200);  // 10 GeV bins
   //   Ht1pfmetParZvsNvtx[i]->SetLineColor(i);
   // }
   
   Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 390., 410., 430., 450., 470., 500., 530., 560, 600., 640., 700., 800.};
   Int_t nBinsForResponse = sizeof(ZptBinEdges)/sizeof(Double_t) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements
   TH1D *HZptBinned[nBinsForResponse];
   //the following histograms will give the distribution of met|| / wzpt. The mean value will be used to create the response curve, that is (<met|| / wzpt>) vs wzpt
   // for each point, wzpt will be taken as the average wzpt in the range considered
   TH1D *HmumetParZratio[nBinsForResponse];
   TH1D *HmumetOrtZvsZpt[nBinsForResponse];
   TH1D *HmumetParZvsZpt[nBinsForResponse]; 
   //cout<<"nbins = "<<nBinsForResponse<<endl;
   for (Int_t i = 0; i < nBinsForResponse; i++) {   
     //HZptBinned[i] are histograms with 5 bins in the range given by ZptBinEdges[i] and ZptBinEdges[i+1]
     // the mean wzpt in each bin will be computed as the histogram's mean
     HZptBinned[i] = new TH1D(Form("HZptBinned[%i]",i),"",5,ZptBinEdges[i],ZptBinEdges[i+1]); 
     HmumetParZratio[i] = new TH1D(Form("HmumetParZratio[%i]",i),"",100,0.0,2.0); 
     HmumetOrtZvsZpt[i] = new TH1D(Form("HmumetOrtZvsZpt[%i]",i),"",250,200,700); 
     HmumetParZvsZpt[i] = new TH1D(Form("HmumetParZvsZpt[%i]",i),"",250,200,700); 
   }

   Double_t nev1j=0, nev2j=0, nev3j=0;                                           //# of events with 1,2,3 jets
   Double_t nev1jmore30=0, nev2jmore30=0, nev3jmore30=0;        //# of events with jet 1,2,3 with pt >30 GeV 
   Double_t nevj1pteq0=0, nevj2pteq0=0, nevj3pteq0=0;                 //# of events with jet 1,2,3 with pt = 30 GeV 

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     /**************************************************/
     // variables used in the loop or in the following
     /**************************************************/

     // Double_t wgtTest = wgt*K_FACTOR;
     // HtestNevents->Fill(wgtTest);
     // HtestNeventsBis->Fill(signaljetpt,wgtTest);

     UInt_t eventMask = 0;    //mask that characterize the event. It's initialized to 0 for each event and is built inside the cut loop 

     Double_t newwgt = wgt*K_FACTOR*LUMI;  //the correct weight to use is newwgt, LUMI is useful to compare to Adish's data
     //Double_t newwgt = 1.;
     nWeightedEvents[0] += newwgt;  // this element will represent the total number of weighted events (no selection applied so far)

     TLorentzVector jet1, jet2, jet3;
     TLorentzVector nu1, nu2, nu1plusnu2, Z;
     TLorentzVector LVt1pfmet, LVmumet;
     Double_t lepInvariantMass=0, leppt=0;
     jet1.SetPtEtaPhiM(signaljetpt,signaljeteta,signaljetphi,0);
     jet2.SetPtEtaPhiM(secondjetpt,secondjeteta,secondjetphi,0);
     jet3.SetPtEtaPhiM(thirdjetpt,thirdjeteta,thirdjetphi,0);
     LVt1pfmet.SetPtEtaPhiM(t1pfmet,0,t1pfmetphi,0);
     LVmumet.SetPtEtaPhiM(mumet,0,mumetphi,0);
     Z.SetPtEtaPhiM(wzpt,wzeta,wzphi,wzmass);
     nu1.SetPtEtaPhiM(l1pt,l1eta,l1phi,0);
     nu2.SetPtEtaPhiM(l2pt,l2eta,l2phi,0);

     nu1plusnu2 = nu1+nu2;
     // it's useful to store the following quantities so that there's no need to recompute them every time I need them
     lepInvariantMass = nu1plusnu2.Mag();
     leppt = nu1plusnu2.Pt();

     /**************************************************/
     // computing met resolutions
     /**************************************************/

     Double_t dphiZmumet_pi_pi = LVmumet.DeltaPhi(Z);
     Double_t mumetPar = mumet*TMath::Cos(dphiZmumet_pi_pi);
     //cout<<"mumetPar = "<<mumetPar<<endl;
     Double_t mumetOrt = mumet*TMath::Sin(dphiZmumet_pi_pi);
     // Double_t dphiZt1pfmet_pi_pi = LVt1pfmet.DeltaPhi(Z);
     // Double_t t1pfmetPar = t1pfmet*TMath::Cos(dphiZt1pfmet_pi_pi);
     // Double_t t1pfmetOrt = t1pfmet*TMath::Sin(dphiZt1pfmet_pi_pi);

     if (wzpt > 250) {   //this corresponds to trigger efficiency plateaux (actually it would be mumet, but mumet recoils against wzpt so their pt should be of the same order of magnitude). In any case, note that in our tree is always mumet > 200 

       Int_t nvtxBin = nvtx-FIRST_NVTX;
  
       if ((nvtxBin >= 0) && (nvtx < (NVTXS + FIRST_NVTX))) {

	 if (wzpt < 500) {                        // (met||-wzpt) distribution's width depends on pt, thus I use this range
	   HmumetParZvsNvtx[nvtxBin]->Fill(mumetPar-wzpt,newwgt);
	   //Ht1pfmetParZvsNvtx[nvtxBin]->Fill(t1pfmetPar-wzpt,newwgt);
	 }
     	 HmumetOrtZvsNvtx[nvtxBin]->Fill(mumetOrt,newwgt);       
	 //Ht1pfmetOrtZvsNvtx[nvtxBin]->Fill(t1pfmetOrt,newwgt); 
	 HmumetX[nvtxBin]->Fill(mumet*TMath::Cos(mumetphi),newwgt);  
	 HmumetY[nvtxBin]->Fill(mumet*TMath::Sin(mumetphi),newwgt);  
    
       }

     }
     
     /**************************************************/
     // computing met responses
     /**************************************************/

     // first of all I make sure that wzpt is in the appropriate range
     if ((wzpt > ZptBinEdges[0]) && (wzpt < ZptBinEdges[nBinsForResponse])) {

       Int_t bin = myGetBin(wzpt,ZptBinEdges,nBinsForResponse);
       //cout<<"bin = "<<bin<<endl;
       HZptBinned[bin]->Fill(wzpt,newwgt);        
       HmumetParZratio[bin]->Fill(mumetPar/wzpt,newwgt);     //the mean value of this histogram is the response
       HmumetOrtZvsZpt[bin]->Fill(mumetOrt,newwgt);
       HmumetParZvsZpt[bin]->Fill(mumetPar-wzpt,newwgt);

     }

     /**************************************************/
     // counting number of events with exactly 1, 2 ,3
     /**************************************************/

     if (njets == 1) nev1j += newwgt;
     if (njets == 2) nev2j += newwgt;
     if (njets == 3) nev3j += newwgt;
     if (signaljetpt >30) nev1jmore30 += newwgt;
     if (signaljetpt == 0) nevj1pteq0 += newwgt;
     if (secondjetpt >30) nev2jmore30 += newwgt;
     if (secondjetpt == 0) nevj2pteq0 += newwgt;     
     if (thirdjetpt >30) nev3jmore30 += newwgt;
     if (thirdjetpt == 0) nevj3pteq0 += newwgt;  

     /**************************************************/
     // filling histograms without any cut applied
     /**************************************************/
     /*---------------------------------- jets ---------------------------------------*/
 
     //ncut = 0 means no cuts for filling the histogram
     HnjetsnoC->Fill(njets,newwgt);
     Hjet1pt[0]->Fill(signaljetpt,newwgt);  
     Hjet1etanoC->Fill(signaljeteta,newwgt);
     Hjet2ptnoC->Fill(secondjetpt,newwgt);
     if (secondjetpt != 0) {
       Hjet2etanoC->Fill(secondjeteta,newwgt);
       Hjet1jet2dphinoC->Fill(jetjetdphi,newwgt);
     }

     /*---------------------------------- Z boson ---------------------------------------*/
     HZptnoC->Fill(wzpt,newwgt);
     HZmassnoC->Fill(wzmass,newwgt);

     /*---------------------------------- leptons ---------------------------------------*/
     HlepptnoC->Fill(leppt,newwgt);
     HlepmassnoC->Fill(lepInvariantMass,newwgt);

     /*---------------------------------- mumet ---------------------------------------*/
     HmumetnoC->Fill(mumet,newwgt);

     //***********************************
     //         STARTING WITH CUTS
     //***********************************  
     //cuts on 2nd jet must be applied only for 2 jet events (I could have 1 jet but the current implementation would apply a cut on the 2nd jet nonetheless): thus, if njets == 2 cuts are actually applied and the corresponding bits of eventMask are set to 1 or 0 depending on the selection being passed or not respectively, otherwise those bits are automatically set to 1 (which is like saying that that selection is irrelevant). The problem and the reason for doing so is that the globalCutMask is defined before the cut loop and I cannot choose whether to set its 2nd jet bits to 0 or 1 according to the effective number of jets before looping on events: I must set them to 1 in advance, thus declaring that those cuts will be applied when I have to

     eventMask += njetsC.addToMask(njets);
     eventMask += jet1ptC.addToMask(signaljetpt);
     eventMask += jet1etaC.addToMask(fabs(signaljeteta));
     eventMask += jet1CHfracC.addToMask(signaljetCHfrac);
     eventMask += jet1NHfracC.addToMask(signaljetNHfrac);
     eventMask += jet1EMfracC.addToMask(signaljetEMfrac);
     if (njets == 2) {
       //eventMask += jet2ptC.addToMask(secondjetpt);
       eventMask += jet2etaC.addToMask(fabs(secondjeteta));
       eventMask += jet2NHfracC.addToMask(secondjetNHfrac);
       eventMask += jet2EMfracC.addToMask(secondjetEMfrac);
       eventMask += jet1jet2dphiC.addToMask(fabs(jetjetdphi));
     } else if (njets == 1) {
       //eventMask += jet2ptC.get2ToId();
       eventMask += jet2etaC.get2ToId();
       eventMask += jet2NHfracC.get2ToId();
       eventMask += jet2EMfracC.get2ToId();
       eventMask += jet1jet2dphiC.get2ToId();
     }
     eventMask += nmuonsC.addToMask(nmuons);
     eventMask += nelectronsnewC.addToMask(nelectronsnew);
     eventMask += ntausC.addToMask(ntaus);
     eventMask += mumet250C.addToMask(mumet);
     eventMask += mumet300C.addToMask(mumet);
     eventMask += mumet400C.addToMask(mumet);
     eventMask += mumet500C.addToMask(mumet);
     //cout<<"event = "<<jentry<<"       eventMask = "<<eventMask<<endl;

     //warning: mask's index starts from 0 (first cut) while Hjet1pt[] starts with index 1 because 0 stands for no cuts. Maybe I should change definition of histogram...
     if ( (eventMask & globalCutMask[0]) == globalCutMask[0]) {
       nWeightedEvents[1] += newwgt;
       Hjet1pt[1]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[1]) == globalCutMask[1]) {
       nWeightedEvents[2] += newwgt;
       Hjet1pt[2]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[2]) == globalCutMask[2]) {
       nWeightedEvents[3] += newwgt;
       Hjet1pt[3]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[3]) == globalCutMask[3]) {
       nWeightedEvents[4] += newwgt;
       Hjet1pt[4]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[4]) == globalCutMask[4]) {
       nWeightedEvents[5] += newwgt;
       Hjet1pt[5]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[5]) == globalCutMask[5]) {
       nWeightedEvents[6] += newwgt;
       Hjet1pt[6]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[6]) == globalCutMask[6]) {
       nWeightedEvents[7] += newwgt;
       Hjet1pt[7]->Fill(signaljetpt, newwgt);
     }
     if ( (eventMask & globalCutMask[7]) == globalCutMask[7]) {
       //here tau veto is fulfilled
       nWeightedEvents[8] += newwgt;       
       Hjet1pt[8]->Fill(signaljetpt, newwgt);
       //cut on mumet > 250 is the first cut: "allCMet250" means all cut flow with cut on met >250 (the following cuts are met> 300 and so on)
       HnjetsallCMet250->Fill(njets,newwgt);
       Hjet1etaallCMet250->Fill(signaljeteta,newwgt);
       Hjet2ptallCMet250->Fill(secondjetpt,newwgt);
       if (secondjetpt != 0) {
	 Hjet2etaallCMet250->Fill(secondjeteta,newwgt);
	 Hjet1jet2dphiallCMet250->Fill(jetjetdphi,newwgt);
       }
       HmumetallCMet250->Fill(mumet,newwgt);
     }
     if ( (eventMask & globalCutMask[8]) == globalCutMask[8]) {
       nWeightedEvents[9] += newwgt;      
       Hjet1pt[9]->Fill(signaljetpt, newwgt);
       Hjet2ptallCMet300->Fill(secondjetpt,newwgt);       
     }  
     if ( (eventMask & globalCutMask[9]) == globalCutMask[9]) {
        nWeightedEvents[10] += newwgt;
	Hjet1pt[10]->Fill(signaljetpt, newwgt);
	Hjet2ptallCMet400->Fill(secondjetpt,newwgt);
     }
     if ( (eventMask & globalCutMask[10]) == globalCutMask[10]) {
       nWeightedEvents[11] += newwgt;
       Hjet1pt[11]->Fill(signaljetpt, newwgt);
       HnjetsallCMet500->Fill(njets,newwgt);
       Hjet1etaallCMet500->Fill(signaljeteta,newwgt);
       Hjet2ptallCMet500->Fill(secondjetpt,newwgt);
       if (secondjetpt != 0) {
	 Hjet2etaallCMet500->Fill(secondjeteta,newwgt);
	 Hjet1jet2dphiallCMet500->Fill(jetjetdphi,newwgt);
       }
     }
     

   ///***********************************
     //       END OF CUTS
     //***********************************

     //--------------------------------------------------------------------------------//
   //   END OF  LOOP
   //--------------------------------------------------------------------------------//

   }

   Double_t nwentries = nWeightedEvents[0];                //storage for clearer code in the following 
   Double_t nwentriesCut1 = nWeightedEvents[1];
   Double_t nevmore1j = nwentries-nev1j;
   Double_t nevmore2j = nevmore1j-nev2j;
   Double_t nevmore3j = nevmore2j-nev3j;
   Double_t nevj1ptneq0 = nwentries-nevj1pteq0;
   Double_t nevj2ptneq0 = nwentries-nevj2pteq0;
   Double_t nevj3ptneq0 = nwentries-nevj3pteq0;

   cout<<endl;   
   cout<<"-----------------------------------------------------------------------------------------------"<<endl;  
   cout<<"total # of weighted events: "<<nwentries<<"     Hjet1pt[0]->GetSumOfWeights = "<<Hjet1pt[0]->GetSumOfWeights()<<endl;
   cout<<"-----------------------------------------------------------------------------------------------"<<endl;  
   for (Int_t i = 1; i <= NCUTS; i++) {
     cout<<"cut "<<setw(2)<<i<<" : nWeightedEvents = "<<nWeightedEvents[i]<<"     Hjet1pt->GetSumOfWeights = "<<Hjet1pt[i]->GetSumOfWeights()<<endl;
   }
   cout<<endl;   
   myPrintEventYields(cout,LUMI,NCUTS,nWeightedEvents);
   cout<<"Following numbers concerning jets refers to jets with pt>30 GeV"<<endl;
   cout<<setw(25)<<left<<"event type"<<"fraction"<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"exactly 1 jet "<<nev1j/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"exactly 2 jets "<<nev2j/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"exactly 3 jets "<<nev3j/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"more than 1 jet "<<nevmore1j/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"more than 2 jet "<<nevmore2j/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"more than 3 jet "<<nevmore3j/nwentries<<endl;
   cout<<"--------------------------------------------------------"<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"signaljetpt>30GeV "<<nev1jmore30/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"secondjetpt>30GeV "<<nev2jmore30/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"thirdjetpt>30GeV "<<nev3jmore30/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"signaljetpt = 0 "<<nevj1pteq0/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"secondjetpt = 0 "<<nevj2pteq0/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"thirdjetpt = 0 "<<nevj3pteq0/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"signaljetpt != 0 "<<nevj1ptneq0/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"secondjetpt != 0 "<<nevj2ptneq0/nwentries<<endl;
   cout<<setprecision(4)<<setw(25)<<left<<"thirdjetpt != 0 "<<nevj3ptneq0/nwentries<<endl;
   cout<<"--------------------------------------------------------"<<endl;
   cout<<"-----------------------------------------------------------------------------------"<<endl;  
   cout<<endl;  

  //------------------------------------------------
   // beginning of file
  //------------------------------------------------

   //the following file is to record the last sequence of cuts used 

   ofstream recoFile("lastCuts.txt",ios::out);
     if ( !recoFile.is_open() ) {
       cout<<"Error: unable to open file "<<"lastCuts.txt"<<" !"<<endl;
     } else {
       cout<<"Saving latest results of cut based analysis in file lastCuts.txt ..."<<endl;
       //when writing on file, the date is printed as well unless an error occurs
       if ( system("date>>lastCuts.txt") != 0) cout<<"Error during \"system(\"date>>lastCuts.txt\")\" call"<<endl;   
       cut::printCutFlow(recoFile, NCUTS, singleCutMask);
       myPrintEventYields(recoFile,LUMI,NCUTS,nWeightedEvents);
       recoFile.close();
     } 

   char answer = '\0';
   const char* filename = "cutFlowData.txt";   //this file is like a logbook

   answer = myAskToSaveFile(filename);

   if (answer == 'y') {

     ofstream myfile(filename,ios::app);

     if ( !myfile.is_open() ) {
       cout<<"Error: unable to open file "<<filename<<" !"<<endl;

     } else {
       //when writing on file, the date is printed as well unless an error occurs
       if ( system("date>>cutFlowData.txt") != 0) cout<<"Error during \"system(\"date>>cutFlowData.txt\")\" call"<<endl;  
       cut::printCutFlow(myfile, NCUTS, singleCutMask);
       myPrintEventYields(myfile,LUMI,NCUTS,nWeightedEvents);
       myfile.close();

     }
 
   }
   //------------------------------------------------
   //end of file
   //------------------------------------------------

   const char *texfname = "table.tex";
   
   answer = myAskToSaveFile(texfname);

   if (answer == 'y') {

     //creating a .tex file to build tables with data
     FILE *fp;	 
     //cout<<"name given to file: "<<texfname<<endl;
     if ( (fp=fopen(texfname,"w")) == NULL) {
       cout<<"Error: '"<<texfname<<"' not opened"<<endl;
     } else {
       cout<<"creating file '"<<texfname<<"' ..."<<endl;
       myAddDefaultPackages(fp,texfname);
       fprintf(fp,"\\begin{document}\n");
       fprintf(fp,"\n");
       //building table
       fprintf(fp,"\\begin{table}\n");
       fprintf(fp,"\\caption{\\emph{Entries after each cut, normalised to %4.2lf $fb^{-1}$; n is the number of entries; aA is the ratio between n and entries of first cut;"
	       " aR is the ratio between n and n-1 (for n = 1, it's set to 1.0). Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV. "
               "}}\n",(Double_t)LUMI);
       fprintf(fp,"\\[\n");
       fprintf(fp,"\\begin{array}{rrrrr}\n");
       fprintf(fp,"\\toprule\n");
       fprintf(fp,"\\ \\textbf{cut}& \\textbf{definition}  & n  & aA & aR \\\\ \n");
       fprintf(fp,"\\midrule\n");
       for (Int_t i = 1; i <= NCUTS; i++) {
	 if (i == 1) {
	   fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i,nWeightedEvents[i],1.0,1.0);
	 } else {
	   fprintf(fp," %i  & & %6.0lf & %1.4lf & %1.4lf \\\\\n",i,nWeightedEvents[i],nWeightedEvents[i]/nwentriesCut1,nWeightedEvents[i]/nWeightedEvents[i-1]);
	 }
	 for (Int_t j = 0; j < cut::getNCuts(); j++) {
	   if ((singleCutMask[i-1] >> j) & 1) {
	     string str = cut::listOfCuts[j]->getCutDefinition();
	     fprintf(fp,"& %s & & & \\\\\n",str.c_str());
	   }
	 }
	 fprintf(fp,"\\midrule\n");
       }
       fprintf(fp,"\\bottomrule\n");
       fprintf(fp,"\\end{array}\n");
       fprintf(fp,"\\]\n");
       fprintf(fp,"\\end{table}\n");
       //end of table
       fprintf(fp,"\n");
       fprintf(fp,"\\end{document}\n");      

       fclose(fp);

     }

   }

   //saving histograms on file .root
   TFile *rootFile = TFile::Open("histograms.root","RECREATE");

   for (Int_t i = 0; i <= NCUTS; i++) {
     Hjet1pt[i]->Write();
   }

   Hjet1etanoC->Write();
   Hjet1etaallCMet250->Write();
   Hjet1etaallCMet500->Write();
   Hjet2ptnoC->Write();
   Hjet2ptallCMet250->Write();
   Hjet2ptallCMet300->Write();
   Hjet2ptallCMet400->Write();
   Hjet2ptallCMet500->Write();
   Hjet2etanoC->Write();
   Hjet2etaallCMet250->Write();
   Hjet2etaallCMet500->Write();
   Hjet1jet2dphinoC->Write();
   Hjet1jet2dphiallCMet250->Write();
   Hjet1jet2dphiallCMet500->Write();
   HnjetsnoC->Write();
   HnjetsallCMet250->Write(); 
   HnjetsallCMet500->Write();
   HmumetnoC->Write();
   HmumetallCMet250->Write(); 

   for (Int_t i = 0; i < NVTXS; i++) {
     HmumetOrtZvsNvtx[i]->Write(); 
     HmumetParZvsNvtx[i]->Write(); 
     HmumetX[i]->Write(); 
     HmumetY[i]->Write();  
   }

   for (Int_t i = 0;  i < nBinsForResponse; i++) {
     HZptBinned[i]->Write(); 
     HmumetParZratio[i]->Write(); 
     HmumetOrtZvsZpt[i]->Write(); 
     HmumetParZvsZpt[i]->Write(); 
   }

   rootFile->Close();

   /*
GetEffectiveEntries() returns ((Sum w)^2)/(Sum(w^2)) where w means weight and Sum is the sum over the arguments. This number is the number of events that would be needed by an unweighted histogram to have the same statistical power as the weighted one
As an example, if I have 5 events with weight 1 and other 5 with weight 0.2, I get 6.9 effective events instead of 10. 

To get the number of entries taking weights into account (that is to say, the bin content or integral of the histogram) I must use GetSumOfWeights() 
    */

     /********************************************************/
     //  ||||||||||||||||||||||||||       CANVASES             ||||||||||||||||||||||||
     /********************************************************/

/************************************/
 //                     JET1 PT, TEST FOR MASKS
 /************************************/

   // TCanvas *Cjet1pt = new TCanvas("Cjet1pt","signaljet pt");
   // Cjet1pt->SetLogy();
   // Cjet1pt->SetGridx();
   // Cjet1pt->SetGridy();
   // Hjet1pt[1]->Draw("HE");
   // Hjet1pt[1]->SetStats(kFALSE);
   // Cjet1pt->SaveAs(MY_PDF_PATH"jet1ptMaskTest.pdf");

 /************************************/
 //                     JET1 PT
 /************************************/

   TCanvas *Cjet1pt = new TCanvas("Cjet1pt","signaljet pt");
   TLegend *LegCjet1pt = new TLegend(0.7,0.5,0.9,0.9);
   Cjet1pt->SetLogy();
   Cjet1pt->SetGridy();
   Cjet1pt->SetGridx();
   Hjet1pt[0]->Draw("HE");
   Hjet1pt[0]->SetStats(kFALSE);
   LegCjet1pt->AddEntry(Hjet1pt[0],"no cuts","lf");
   Hjet1pt[0]->GetXaxis()->SetTitle("signaljet pt [GeV]");
   Hjet1pt[0]->GetYaxis()->SetTitle("events / 50 GeV");
   Hjet1pt[0]->GetYaxis()->SetTitleOffset(1.4);
   if (ENABLE_TRIGGER_CUTS) {
     Hjet1ptTriggerC->Draw("SAME HE");
     Hjet1ptTriggerC->SetStats(kFALSE);
     Hjet1ptTriggerC->SetLineColor(14);
     LegCjet1pt->AddEntry(Hjet1ptTriggerC,"trigger cuts","lf");
   }
   for ( Int_t i=1; i<=NCUTS; i++ ) {
     Hjet1pt[i]->Draw("SAME HE");
     Hjet1pt[i]->SetStats(kFALSE);
     LegCjet1pt->AddEntry(Hjet1pt[i],Form("cut %i",i),"lf");
   }
   LegCjet1pt->Draw(); 
   LegCjet1pt->SetMargin(0.5); 
   gPad->RedrawAxis();
   gPad->RedrawAxis("G");    //without this, the filled area covers axis' s ticks, same for grid 
   if (ENABLE_TRIGGER_CUTS) {
     TPaveText *PaveCjet1pt = new TPaveText(0.55,0.7,0.7,0.9,"NB NDC"); //NB abort borders drawing (no shadows), NDC is to use NDC coordinates
     PaveCjet1pt->AddText("Trigger");
     PaveCjet1pt->AddLine(0.0,0.75,1.0,0.75);
     PaveCjet1pt->AddText(Form("jet pt > %i GeV",JET_TRIGGER));
     PaveCjet1pt->AddText(Form("pfmet > %i GeV",PFMET_TRIGGER));
     PaveCjet1pt->AddText(Form("pfmet > %i GeV",MHT_TRIGGER));
     PaveCjet1pt->Draw();
     Cjet1pt->SaveAs(Form(MY_PDF_PATH"jet1ptCutFlowTrig_%i_%i_%i.pdf",JET_TRIGGER,PFMET_TRIGGER,MHT_TRIGGER));  
   } else {
     Cjet1pt->SaveAs(MY_PDF_PATH"jet1ptFill.pdf");
   }

   TCanvas *Cjet1ptBis = new TCanvas("Cjet1ptBis","signaljet pt");
   TLegend *LegCjet1ptBis = new TLegend(0.6,0.70,0.9,0.9);
   Cjet1ptBis->SetLogy();
   Hjet1pt[0]->UseCurrentStyle();
   Hjet1pt[8]->UseCurrentStyle();
   Hjet1pt[11]->UseCurrentStyle();
   Hjet1pt[0]->Draw("HE");
   Hjet1pt[0]->SetStats(kFALSE);
   Hjet1pt[0]->SetLineColor(1);
   LegCjet1ptBis->AddEntry(Hjet1pt[0],"no cuts","l");
   Hjet1pt[0]->GetXaxis()->SetTitle("signaljet pt [GeV]");
   Hjet1pt[0]->GetYaxis()->SetTitle("events / 50 GeV");
   Hjet1pt[0]->GetYaxis()->SetTitleOffset(1.4);
   LegCjet1ptBis->AddEntry((TObject *)0,"all cuts and MET","");
   Hjet1pt[8]->Draw("SAME HE");
   Hjet1pt[8]->SetStats(kFALSE);
   Hjet1pt[8]->SetLineColor(2);
   LegCjet1ptBis->AddEntry(Hjet1pt[8]," > 250 GeV ","l");
   Hjet1pt[11]->Draw("SAME HE");
   Hjet1pt[11]->SetStats(kFALSE);
   Hjet1pt[11]->SetLineColor(4);
   LegCjet1ptBis->AddEntry(Hjet1pt[11]," > 500 GeV ","l");
   LegCjet1ptBis->Draw(); 
   LegCjet1ptBis->SetMargin(0.5); 
   Cjet1ptBis->SaveAs(MY_PDF_PATH"jet1ptBis.pdf");

/************************************/
 //                     JET1 PT (TO TRIG OR NOT TO TRIG?)
 /************************************/

   // if (ENABLE_TRIGGER_CUTS) {
   //   TCanvas *Cjet1ptTr_comp = new TCanvas("Cjet1ptTr_comp","signaljet pt");
   //   TLegend *LegCjet1ptTr_comp = new TLegend(0.7,0.7,0.9,0.9);
   //   Cjet1ptTr_comp->SetLogy();
   //   Cjet1ptTr_comp->SetGridy();
   //   Cjet1ptTr_comp->SetGridx();
   //   Hjet1pt[0]->Draw("HE");
   //   Hjet1pt[0]->SetStats(kFALSE);
   //   LegCjet1ptTr_comp->AddEntry(Hjet1pt[0],"no cuts","l");
   //   Hjet1pt[0]->GetXaxis()->SetTitle("signaljet pt [GeV]");
   //   Hjet1pt[0]->GetXaxis()->SetRangeUser(-25.,1500.);
   //   Hjet1pt[0]->GetYaxis()->SetTitle("events / 50 GeV");
   //   Hjet1pt[0]->GetYaxis()->SetTitleOffset(1.4);
   //   Hjet1ptTriggerC->Draw("SAME HE");
   //   Hjet1ptTriggerC->SetStats(kFALSE);
   //   Hjet1ptTriggerC->SetLineColor(kRed);
   //   LegCjet1ptTr_comp->AddEntry(Hjet1ptTriggerC,"trigger cuts","l");
   //   LegCjet1ptTr_comp->Draw(); 
   //   TPaveText *PaveCjet1ptTr_comp = new TPaveText(0.55,0.7,0.7,0.9,"NB NDC");
   //   PaveCjet1ptTr_comp->AddText("Trigger");
   //   PaveCjet1ptTr_comp->AddLine(0.0,0.75,1.0,0.75);
   //   PaveCjet1ptTr_comp->AddText(Form("jet pt > %i GeV",JET_TRIGGER));
   //   PaveCjet1ptTr_comp->AddText(Form("pfmet > %i GeV",PFMET_TRIGGER));
   //   PaveCjet1ptTr_comp->AddText(Form("pfmet > %i GeV",MHT_TRIGGER));
   //   PaveCjet1ptTr_comp->Draw();
   //   LegCjet1ptTr_comp->SetMargin(0.5); 
   //   Cjet1ptTr_comp->SaveAs(Form(MY_PDF_PATH"jet1ptTrigEffect_%i_%i_%i.pdf",JET_TRIGGER,PFMET_TRIGGER,MHT_TRIGGER));
   // }

// /************************************/
//  //                     JET1 ETA
//  /************************************/

   TCanvas *Cjet1eta = new TCanvas("Cjet1eta","signaljet eta");
   TLegend *LegCjet1eta = new TLegend(0.60,0.70,0.9,0.9);
   Cjet1eta->SetLogy();
   Hjet1etanoC->Draw("HE");
   Hjet1etanoC->SetStats(kFALSE);
   Hjet1etanoC->SetLineColor(1);
   LegCjet1eta->AddEntry(Hjet1etanoC,"no cuts","l");
   Hjet1etanoC->GetXaxis()->SetTitle("signaljet eta");
   Hjet1etanoC->GetYaxis()->SetTitle("events / 0.1");
   Hjet1etanoC->GetYaxis()->SetTitleOffset(1.4);
   LegCjet1eta->AddEntry((TObject *)0,"all cuts and MET","");
   Hjet1etaallCMet250->Draw("SAME HE");
   Hjet1etaallCMet250->SetStats(kFALSE);
   Hjet1etaallCMet250->SetLineColor(2);
   LegCjet1eta->AddEntry(Hjet1etaallCMet250," > 250 GeV ","l");
   Hjet1etaallCMet500->Draw("SAME HE");
   Hjet1etaallCMet500->SetStats(kFALSE);
   Hjet1etaallCMet500->SetLineColor(4);
   LegCjet1eta->AddEntry(Hjet1etaallCMet500," > 500 GeV ","l");
   LegCjet1eta->Draw(); 
   LegCjet1eta->SetMargin(0.5); 
   Cjet1eta->SaveAs(MY_PDF_PATH"jet1eta.pdf");  

// /************************************/
//  //                     JET2 PT
//  /************************************/

   TCanvas *Cjet2pt = new TCanvas("Cjet2pt","secondjet pt");
   TLegend *LegCjet2pt = new TLegend(0.60,0.65,0.9,0.9);
   Cjet2pt->SetLogy();
   Hjet2ptnoC->Draw("HE");
   Hjet2ptnoC->SetStats(kFALSE);
   Hjet2ptnoC->SetLineColor(1);
   LegCjet2pt->AddEntry(Hjet2ptnoC,"no cuts","l");
   Hjet2ptnoC->GetXaxis()->SetTitle("secondjet pt [GeV]");
   Hjet2ptnoC->GetYaxis()->SetTitle("events / 50 GeV");
   Hjet2ptnoC->GetYaxis()->SetTitleOffset(1.4);
   LegCjet2pt->AddEntry((TObject *)0,"all cuts and MET","");
   Hjet2ptallCMet250->Draw("SAME HE");
   Hjet2ptallCMet250->SetStats(kFALSE);
   Hjet2ptallCMet250->SetLineColor(2);
   LegCjet2pt->AddEntry(Hjet2ptallCMet250," > 250 GeV ","l");
   Hjet2ptallCMet300->Draw("SAME HE");
   Hjet2ptallCMet300->SetStats(kFALSE);
   Hjet2ptallCMet300->SetLineColor(3);
   LegCjet2pt->AddEntry(Hjet2ptallCMet300," > 300 GeV ","l");
   Hjet2ptallCMet400->Draw("SAME HE");
   Hjet2ptallCMet400->SetStats(kFALSE);
   Hjet2ptallCMet400->SetLineColor(4);
   LegCjet2pt->AddEntry(Hjet2ptallCMet400," > 400 GeV ","l");
   Hjet2ptallCMet500->Draw("SAME HE");
   Hjet2ptallCMet500->SetStats(kFALSE);
   Hjet2ptallCMet500->SetLineColor(6);
   LegCjet2pt->AddEntry(Hjet2ptallCMet500," > 500 GeV ","l");
   LegCjet2pt->Draw(); 
   LegCjet2pt->SetMargin(0.5); 
   Cjet2pt->SaveAs(MY_PDF_PATH"jet2pt.pdf");

// /************************************/
//  //                     JET2 ETA
//  /************************************/

   TCanvas *Cjet2eta = new TCanvas("Cjet2eta","secondjet eta");
   TLegend *LegCjet2eta = new TLegend(0.60,0.70,0.9,0.9);
   Cjet2eta->SetLogy();
   Hjet2etanoC->Draw("HE");
   Hjet2etanoC->SetStats(kFALSE);
   Hjet2etanoC->SetLineColor(1);
   Hjet2etanoC->GetYaxis()->SetRangeUser(2.0, 6000.0);
   LegCjet2eta->AddEntry(Hjet2etanoC,"no cuts","l");
   Hjet2etanoC->GetXaxis()->SetTitle("secondjet eta");
   Hjet2etanoC->GetYaxis()->SetTitle("events / 0.1");
   Hjet2etanoC->GetYaxis()->SetTitleOffset(1.4);
   LegCjet2eta->AddEntry((TObject *)0,"all cuts and MET","");
   Hjet2etaallCMet250->Draw("SAME HE");
   Hjet2etaallCMet250->SetStats(kFALSE);
   Hjet2etaallCMet250->SetLineColor(2);
   LegCjet2eta->AddEntry(Hjet2etaallCMet250," > 250 GeV ","l");
   Hjet2etaallCMet500->Draw("SAME HE");
   Hjet2etaallCMet500->SetStats(kFALSE);
   Hjet2etaallCMet500->SetLineColor(4);
   LegCjet2eta->AddEntry(Hjet2etaallCMet500," > 500 GeV ","l");
   LegCjet2eta->Draw(); 
   LegCjet2eta->SetMargin(0.5); 
   Cjet2eta->SaveAs(MY_PDF_PATH"jet2eta.pdf"); 

   TCanvas *Cjet2etaallCMet500 = new TCanvas("Cjet2etaallCMet500","secondjet eta");
   Cjet2etaallCMet500->SetLogy();
   Hjet2etaallCMet500->Draw(" HE");
   //Hjet2etaallCMet500->SetStats(kFALSE);
   Hjet2etaallCMet500->SetLineColor(4);
   Hjet2etaallCMet500->GetXaxis()->SetTitle("secondjet eta");
   Hjet2etaallCMet500->GetYaxis()->SetTitle("events / 0.1");
   Hjet2etaallCMet500->GetYaxis()->SetTitleOffset(1.4);
   Cjet2etaallCMet500->SaveAs(MY_PDF_PATH"jet2etaallCMet500.pdf"); 

   // /************************************/
//  //                   NJETS
//  /************************************/

   TCanvas *Cnjets = new TCanvas("Cnjets","njets");
   TLegend *LegCnjets = new TLegend(0.60,0.70,0.9,0.9);
   Cnjets->SetLogy();
   HnjetsnoC->Draw("HE");
   HnjetsnoC->SetStats(kFALSE);
   HnjetsnoC->SetLineColor(1);
   LegCnjets->AddEntry(HnjetsnoC,"no cuts","l");
   HnjetsnoC->GetXaxis()->SetTitle("# of jets");
   HnjetsnoC->GetYaxis()->SetTitle("events");
   HnjetsnoC->GetYaxis()->SetTitleOffset(1.4);
   LegCnjets->AddEntry((TObject *)0,"all cuts and MET","");
   HnjetsallCMet250->Draw("SAME HE");
   HnjetsallCMet250->SetStats(kFALSE);
   HnjetsallCMet250->SetLineColor(2);
   LegCnjets->AddEntry(HnjetsallCMet250," > 250 GeV ","l");
   HnjetsallCMet500->Draw("SAME HE");
   HnjetsallCMet500->SetStats(kFALSE);
   HnjetsallCMet500->SetLineColor(4);
   LegCnjets->AddEntry(HnjetsallCMet500," > 500 GeV ","l");
   LegCnjets->Draw(); 
   LegCnjets->SetMargin(0.5); 
   Cnjets->SaveAs(MY_PDF_PATH"njets.pdf"); 

//   /************************************/
//  //                 JET1JET2DPHI
//  /************************************/

   TCanvas *Cjet1jet2dphi = new TCanvas("Cjet1jet2dphi","jet1jet2dphi after all cuts");
   TLegend *LegCjet1jet2dphi = new TLegend(0.6,0.70,0.9,0.9);
   Cjet1jet2dphi->SetLogy();
   Hjet1jet2dphinoC->Draw("HE");
   Hjet1jet2dphinoC->SetStats(kFALSE);
   Hjet1jet2dphinoC->SetLineColor(1);
   Hjet1jet2dphinoC->GetXaxis()->SetTitle("jet1jet2dphi");
   Hjet1jet2dphinoC->GetYaxis()->SetTitle("events / 0.1");
   Hjet1jet2dphinoC->GetYaxis()->SetTitleOffset(1.4); 
   Hjet1jet2dphinoC->GetYaxis()->SetRangeUser(8.0,5000.0); 
   Hjet1jet2dphiallCMet250->Draw("SAME HE");
   Hjet1jet2dphiallCMet250->SetStats(kFALSE);
   Hjet1jet2dphiallCMet250->SetLineColor(2);
   Hjet1jet2dphiallCMet500->Draw("SAME HE");
   Hjet1jet2dphiallCMet500->SetStats(kFALSE);
   Hjet1jet2dphiallCMet500->SetLineColor(4);
   LegCjet1jet2dphi->AddEntry(Hjet1jet2dphinoC,"no cuts","l");
   LegCjet1jet2dphi->AddEntry((TObject *)0,"all cuts and","");
   LegCjet1jet2dphi->AddEntry(Hjet1jet2dphiallCMet250,"MET > 250 GeV","l");
   LegCjet1jet2dphi->AddEntry(Hjet1jet2dphiallCMet500,"MET > 500 GeV","l");
   LegCjet1jet2dphi->Draw(); 
   LegCjet1jet2dphi->SetMargin(0.5); 
   Cjet1jet2dphi->SaveAs(MY_PDF_PATH"jet1jet2dphi.pdf");

//   /************************************/
//  //                    MUMET
//  /************************************/

   TCanvas *Cmumet = new TCanvas("Cmumet","mumet after all cuts");
   TLegend *LegCmumet = new TLegend(0.55,0.75,0.9,0.9);
   Cmumet->SetLogy();
   HmumetnoC->Draw("HE");
   HmumetnoC->SetStats(kFALSE);
   HmumetnoC->SetLineColor(2);
   HmumetnoC->GetXaxis()->SetTitle("mumet [GeV]");
   HmumetnoC->GetYaxis()->SetTitle("events / 50 GeV");
   HmumetnoC->GetYaxis()->SetTitleOffset(1.4);
   HmumetallCMet250->Draw("SAME HE");
   HmumetallCMet250->SetStats(kFALSE);
   HmumetallCMet250->SetLineColor(4);
   LegCmumet->AddEntry(HmumetnoC,"no cuts","l");
   LegCmumet->AddEntry(HmumetallCMet250,"all cuts, MET > 250 GeV","l");
   LegCmumet->Draw(); 
   LegCmumet->SetMargin(0.5); 
   Cmumet->SaveAs(MY_PDF_PATH"mumet.pdf");

//   /****************************************************/
//  //                    Z AND LEPTONS INVARIANT MASS (COMPARISON)
//    /****************************************************/

//   TCanvas *CZlepMass_comp = new TCanvas("CZlepMass_comp","leptons' and Z's invariant mass");
//    CZlepMass_comp->SetLogy();
//    HlepmassnoC->ResetStats();
//    HlepmassnoC->GetXaxis()->SetRangeUser(0.,600.);
//    HlepmassnoC->Draw("HE");
//    HlepmassnoC->SetFillColor(20);
//    HlepmassnoC->SetLineColor(kBlue);
//    HlepmassnoC->SetTitle("invariant mass");
//    HlepmassnoC->GetXaxis()->SetTitle("m [GeV]");
//    HlepmassnoC->GetYaxis()->SetTitle("events / 2 GeV");
//    HlepmassnoC->GetYaxis()->SetTitleOffset(1.5); 
//    HlepmassnoC->SetStats(kFALSE);
//    HZmassnoC->Draw("SAME HE");
//    HZmassnoC->SetStats(kFALSE);
//    HZmassnoC->SetLineColor(kRed);
//    TLegend *LegCZlepMass_comp = new TLegend(0.65,0.75,0.9,0.90);
//    LegCZlepMass_comp->AddEntry(HlepmassnoC,"leptons","lf");
//    LegCZlepMass_comp->AddEntry(HZmassnoC,"Z","l");
//    LegCZlepMass_comp->Draw(); 
//    LegCZlepMass_comp->SetMargin(0.5); 
//    CZlepMass_comp->SaveAs(MY_PDF_PATH"ZlepMass_comp.pdf");

//    /*---------------- as the previous canvas but logaritmic x axis and axis range upto 700 ---------*/  
//    TCanvas *CZlepMass_comp_2 = new TCanvas("CZlepMass_comp_2","leptons' and Z's invariant mass");
//    CZlepMass_comp_2->SetLogy();
//    CZlepMass_comp_2->SetLogx();
//    HlepmassnoC->ResetStats();
//    HlepmassnoC->Draw("HE");
//    HlepmassnoC->GetXaxis()->SetRangeUser(0.,700.);
//    HlepmassnoC->SetFillColor(20);
//    HlepmassnoC->SetLineColor(kBlue);
//    //without the following line the histogram inherits changes made on it in the previous canvas
//    //HlepmassnoC->UseCurrentStyle(); 
//    HlepmassnoC->SetTitle("invariant mass");
//    HlepmassnoC->GetXaxis()->SetTitle("m [GeV]");
//    HlepmassnoC->GetYaxis()->SetTitle("events / 2 GeV");
//    HlepmassnoC->GetYaxis()->SetTitleOffset(1.5); 
//    HlepmassnoC->SetStats(kFALSE);
//    HZmassnoC->Draw("SAME HE");
//    HZmassnoC->SetStats(kFALSE);
//    HZmassnoC->SetLineColor(kRed);
//    TLegend *LegCZlepMass_comp_2 = new TLegend(0.65,0.75,0.9,0.90);
//    LegCZlepMass_comp_2->AddEntry(HlepmassnoC,"leptons","lf");
//    LegCZlepMass_comp_2->AddEntry(HZmassnoC,"Z","l");
//    LegCZlepMass_comp_2->Draw(); 
//    LegCZlepMass_comp_2->SetMargin(0.5); 
//    CZlepMass_comp_2->SaveAs(MY_PDF_PATH"ZlepMass_comp_2.pdf");

//    TCanvas *CZlepMassallC_comp = new TCanvas("CZlepMassallC_comp","leptons' and Z's invariant mass");
//    CZlepMassallC_comp->SetLogy();
//    HlepmassallC->ResetStats();
//    HlepmassallC->Draw("HE");
//    HlepmassallC->SetLineColor(kBlue);
//    HlepmassallC->SetFillColor(20);
//    HlepmassallC->SetTitle("invariant mass");
//    HlepmassallC->GetXaxis()->SetTitle("m [GeV]");
//    HlepmassallC->GetYaxis()->SetTitle("events / 2 GeV");
//    HlepmassallC->GetYaxis()->SetTitleOffset(1.5); 
//    HlepmassallC->SetStats(kFALSE);
//    HZmassallC->Draw("SAME HE");
//    HZmassallC->SetStats(kFALSE);
//    HZmassallC->SetLineColor(kRed);
//    TLegend *LegCZlepMassallC_comp = new TLegend(0.65,0.75,0.9,0.90);
//    LegCZlepMassallC_comp->AddEntry(HlepmassallC,"leptons","lf");
//    LegCZlepMassallC_comp->AddEntry(HZmassallC,"Z","l");
//    LegCZlepMassallC_comp->Draw(); 
//    LegCZlepMassallC_comp->SetMargin(0.5); 
//    CZlepMassallC_comp->SaveAs(MY_PDF_PATH"ZlepMass_allC_comp.pdf");

//    /******************************************************/
//    //                    LEPTONS PT (NEUTRINOS FROM Z)
//    /******************************************************/

//    TCanvas *Cleppt = new TCanvas("Cleppt","leptons pt");
//    TLegend *LegCleppt = new TLegend(0.6,0.7,0.9,0.9);
//    Cleppt->SetLogy();
//    HlepptnoC->Draw("HE");
//    HlepptnoC->SetStats(kFALSE);
//    HlepptnoC->SetLineColor(1);
//    HlepptnoC->GetXaxis()->SetTitle("leptons pt [GeV]");
//    HlepptnoC->GetYaxis()->SetTitle("events / 10 GeV");
//    HlepptnoC->GetYaxis()->SetTitleOffset(1.4);
//    HlepptallCnoIso->Draw("SAME HE");
//    HlepptallCnoIso->SetStats(kFALSE);
//    HlepptallCnoIso->SetLineColor(2);
//    HlepptallCnoIso->SetFillColor(18);
//    HlepptallC->Draw("SAME HE");
//    HlepptallC->SetStats(kFALSE);
//    HlepptallC->SetLineColor(4);
//    LegCleppt->AddEntry(HlepptnoC,"no cuts","l");
//    LegCleppt->AddEntry(HlepptallCnoIso,"all cuts, no isol","lf");
//    LegCleppt->AddEntry(HlepptallC,"all cuts","l");
//    LegCleppt->Draw(); 
//    LegCleppt->SetMargin(0.5); 
//    Cleppt->SaveAs(MY_PDF_PATH"leppt.pdf");

//    /******************************************************/
//    //                    LEPTONS AND Z PT (COMPARISON)
//    /******************************************************/

//    TCanvas *CZleppt_comp = new TCanvas("CZleppt_comp","leptons and Z pt");
//    TLegend *LegCZleppt_comp = new TLegend(0.6,0.6,0.9,0.9);
//    CZleppt_comp->SetLogy();
//    HlepptnoC->Draw("HE");
//    HlepptnoC->SetStats(kFALSE);
//    HlepptnoC->SetFillColor(19);
//    HlepptnoC->SetLineColor(1);
//    HlepptnoC->GetXaxis()->SetTitle("pt [GeV]");
//    HlepptnoC->GetYaxis()->SetTitle("events / 10 GeV");
//    HlepptnoC->GetYaxis()->SetTitleOffset(1.4);
//    HZptnoC->Draw("SAME HE");
//    HZptnoC->SetStats(kFALSE);
//    HZptnoC->SetLineColor(4);
//    HlepptallC->Draw("SAME HE");
//    HlepptallC->SetStats(kFALSE);
//    HlepptallC->SetFillColor(17);
//    HlepptallC->SetLineColor(8);
//    HZptallC->Draw("SAME HE");
//    HZptallC->SetStats(kFALSE);
//    HZptallC->SetLineColor(2);
//    LegCZleppt_comp->AddEntry((TObject *)0," -- neutrinos --","");
//    LegCZleppt_comp->AddEntry(HlepptnoC,"no cuts","lf");
//    LegCZleppt_comp->AddEntry(HlepptallC,"all cuts","lf");
//    LegCZleppt_comp->AddEntry((TObject *)0," -- Z boson --","");
//    LegCZleppt_comp->AddEntry(HZptnoC,"no cuts","l");
//    LegCZleppt_comp->AddEntry(HZptallC,"all cuts","l");
//    LegCZleppt_comp->Draw(); 
//    LegCZleppt_comp->SetMargin(0.5); 
//    CZleppt_comp->SaveAs(MY_PDF_PATH"Zleppt_comp.pdf");

   /****************************************************/
   // DOUBLE-SIDE CRYSTAL BALL (DSCB)
   /****************************************************/ 

   // TCanvas *dscbtrial = new TCanvas("dscbtrial","trial");
   // the last parameter (5) is the number of parameter to define the function
   // TF1 *dscb = new TF1("dscb",&my2sideCrystalBall,-10,10,7);
   // parameters are alphasx, nsx, mean (called mu), sigma, normalization (called N), alphadx, ndx
   // dscb->SetParameters(1,2,2,1,10,1,2);
   // dscb->Draw();
   // dscbtrial->SaveAs(MY_PDF_PATH"doubleSideCrystalBall.pdf");

   /****************************************************/
   // FIT TO MET|| AND MET_|_ USING MUMET (HmumetParZratio[1])
   /****************************************************/ 

   TCanvas *CmumetParZFitCB = new TCanvas("CmumetParZFitCB","(mumet|| -Zpt)/Zpt");
   //CmumetParZFitCB->SetLogy();
   //HmumetParZFitCB->ResetStats();
   HmumetParZratio[1]->Draw("HE");
   HmumetParZratio[1]->GetXaxis()->SetTitle("mumet|| / Zpt ");
   HmumetParZratio[1]->GetYaxis()->SetTitle("events / 0.02");
   HmumetParZratio[1]->GetYaxis()->SetTitleOffset(1.4); 
   gStyle->SetOptStat(1000001111); 
   //before fitting with a Crystal Ball I use a gaussian to get the main parameters
   Double_t gaussEdgeL = 0.88;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
   Double_t gaussEdgeR = 1.08;  //right side ...
   HmumetParZratio[1]->Fit("gaus","L 0","",gaussEdgeL,gaussEdgeR);
   Double_t gaussNorm = HmumetParZratio[1]->GetFunction("gaus")->GetParameter(0);
   Double_t gaussMean = HmumetParZratio[1]->GetFunction("gaus")->GetParameter(1);
   Double_t gaussSigma = HmumetParZratio[1]->GetFunction("gaus")->GetParameter(2);
   // fit with Double-side Crystal Ball
   Double_t fncRangeL = 0.8;
   Double_t fncRangeR = 1.35;
   TF1 *cb1 = new TF1("cb1",&my2sideCrystalBall,fncRangeL,fncRangeR,7);
   cb1->SetParNames("alphaL","nL","mu","sigma","N","alphaR","nR");  
   cb1->SetParLimits(cb1->GetParNumber("nL"),0.1,15);
   cb1->SetParLimits(cb1->GetParNumber("nR"),0.1,15);
   //cb1->FixParameter(cb1->GetParNumber("nR"),10.);
   cb1->SetParLimits(cb1->GetParNumber("alphaL"),-10.,-0.01); 
   cb1->SetParLimits(cb1->GetParNumber("alphaR"),0.01,10);
   cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm,(gaussEdgeR-gaussMean)/gaussSigma,5);
   //cb1->SetParameters(1,5,HmumetParZratio[1]->GetBinCenter(HmumetParZratio[1]->GetMaximumBin()),HmumetParZratio[1]->GetRMS(),HmumetParZratio[1]->GetBinContent(HmumetParZratio[1]->GetMaximumBin()),1,5);
   TFitResultPtr frp1 = HmumetParZratio[1]->Fit(cb1,"L I S Q B","HE",fncRangeL,fncRangeR);
   frp1->Print("V");
   TMatrixDSym corrMatrix1 = frp1->GetCorrelationMatrix();
   corrMatrix1.Print(); 
   gStyle->SetOptFit(1100);
   CmumetParZFitCB->SaveAs(MY_PDF_PATH"mumetParZratioFitCB.pdf");

   //I use HmumetOrtZvsNvtx[9] (PU = 15)

   TCanvas *CmumetOrtZFitCB = new TCanvas("CmumetOrtZFitCB","mumet_|_");
   //CmumetOrtZFitCB->SetLogy();
   //HmumetOrtZFitCB->ResetStats();
   HmumetOrtZvsNvtx[9]->Draw("HE");
   HmumetOrtZvsNvtx[9]->GetXaxis()->SetTitle("mumet_|_ [GeV]");
   HmumetOrtZvsNvtx[9]->GetYaxis()->SetTitle("events / 10 GeV");
   HmumetOrtZvsNvtx[9]->GetYaxis()->SetTitleOffset(1.4); 
   gStyle->SetOptStat(1000001111); 
   //before fitting with a Crystal Ball I use a gaussian to get the main parameters
   gaussEdgeL = -30.;
   gaussEdgeR = 30.;
   HmumetOrtZvsNvtx[9]->Fit("gaus","L","",gaussEdgeL,gaussEdgeR);
   gaussNorm = HmumetOrtZvsNvtx[9]->GetFunction("gaus")->GetParameter(0);
   gaussMean = HmumetOrtZvsNvtx[9]->GetFunction("gaus")->GetParameter(1);
   gaussSigma = HmumetOrtZvsNvtx[9]->GetFunction("gaus")->GetParameter(2);
  //fit with Double-side Crystal Ball
   // TF1 *cb2 = new TF1("cb2",&my2sideCrystalBall,-50,50,7);
   // cb2->SetParNames("alphaL","nL","mu","sigma","N","alphaR","nR");
   // cb2->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm,(gaussEdgeR-gaussMean)/gaussSigma,5);
   // //cb2->SetParameters(1,5,HmumetOrtZvsNvtx[9]->GetBinCenter(HmumetOrtZvsNvtx[9]->GetMaximumBin()),HmumetOrtZvsNvtx[9]->GetRMS(),HmumetOrtZvsNvtx[9]->GetBinContent(HmumetOrtZvsNvtx[9]->GetMaximumBin()),1,5);
   // TFitResultPtr frp2 = HmumetOrtZvsNvtx[9]->Fit(cb2,"L I Q S R","HE"); 
   // frp2->Print();
   // TMatrixDSym corrMatrix2 = frp2->GetCorrelationMatrix();
   // corrMatrix2.Print(); 
   gStyle->SetOptFit(1100);
   CmumetOrtZFitCB->SaveAs(MY_PDF_PATH"mumetOrtZFitCB.pdf");

   /************************************/
   //                    MET|| & MET_|_ VS NVTX (using mumet or t1pfmet) AND MET_X & MET_Y
   /************************************/

   Double_t xValues[NVTXS];
   Double_t yValues[NVTXS];
   Double_t yValuesErr[NVTXS];
   for (Int_t i = 0; i < NVTXS; i++) {
     xValues[i] = i + FIRST_NVTX;
     yValues[i] = HmumetParZvsNvtx[i]->GetRMS();
     yValuesErr[i] = HmumetParZvsNvtx[i]->GetRMSError();
   }

   TCanvas *CresolutionMumetParZvsNvtx = new TCanvas("CresolutionMumetParZvsNvtx","resolution ||");
   TGraphErrors *GresolutionMumetParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMumetParZvsNvtx->SetTitle("resolution || from histogram's RMS");
   GresolutionMumetParZvsNvtx->Draw("AP");
   GresolutionMumetParZvsNvtx->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMumetParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMumetParZvsNvtx->GetYaxis()->SetTitle("mumet|| resolution [GeV]");
   GresolutionMumetParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   //resolution vs nvtx should be parametrized as sqrt( sigmaC^2 + nvtx * sigmaPU^2 / 0.7)  (see article about met reconstruction performance for reference)
   //fit to || component
   TF1 *fitResoParVsNvtx = new TF1("fitResoParVsNvtx","sqrt([0] + [1] * x / 0.7)",6.,28.);  // [0] is sigmaC^2 and [1] is sigmaPU^2
   fitResoParVsNvtx->SetParameter(0,200.);
   fitResoParVsNvtx->SetParameter(1,15.);
   TFitResultPtr frpPar = GresolutionMumetParZvsNvtx->Fit(fitResoParVsNvtx,"I S Q B R","");
   frpPar->Print("V");
   Double_t sigmaCPar = TMath::Sqrt(GresolutionMumetParZvsNvtx->GetFunction("fitResoParVsNvtx")->GetParameter(0));
   Double_t sigmaCErrPar = 0.5 * (TMath::Sqrt(GresolutionMumetParZvsNvtx->GetFunction("fitResoParVsNvtx")->GetParError(0))) / TMath::Sqrt(sigmaCPar);
   Double_t sigmaPUPar = TMath::Sqrt(GresolutionMumetParZvsNvtx->GetFunction("fitResoParVsNvtx")->GetParameter(1));
   Double_t sigmaPUErrPar = 0.5 * (TMath::Sqrt(GresolutionMumetParZvsNvtx->GetFunction("fitResoParVsNvtx")->GetParError(1))) / TMath::Sqrt(sigmaPUPar);
   gStyle->SetOptFit(10100);
   CresolutionMumetParZvsNvtx->SaveAs(MY_PDF_PATH"resolutionMumetParZvsNvtx.pdf");

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = HmumetOrtZvsNvtx[i]->GetRMS();
     yValuesErr[i] = HmumetOrtZvsNvtx[i]->GetRMSError();
   }

   TCanvas *CresolutionMumetOrtZvsNvtx = new TCanvas("CresolutionMumetOrtZvsNvtx","resolution _|_");
   TGraphErrors *GresolutionMumetOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMumetOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMumetOrtZvsNvtx->Draw("AP");
   GresolutionMumetOrtZvsNvtx->SetMarkerStyle(7);
   GresolutionMumetOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMumetOrtZvsNvtx->GetYaxis()->SetTitle("mumet_|_ resolution [GeV]");
   GresolutionMumetOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   //resolution vs nvtx should be parametrized as sqrt( sigmaC^2 + nvtx * sigmaPU^2 / 0.7)  (see article about met reconstruction performance for reference)
   //fit to _|_ component
   TF1 *fitResoOrtVsNvtx = new TF1("fitResoOrtVsNvtx","sqrt([0] + [1] * x / 0.7)",6.,28.);  // [0] is sigmaC^2 and [1] is sigmaPU^2
   fitResoOrtVsNvtx->SetParameter(0,200.);
   fitResoOrtVsNvtx->SetParameter(1,15.);
   TFitResultPtr frpOrt = GresolutionMumetOrtZvsNvtx->Fit(fitResoOrtVsNvtx,"I S Q B R","");
   frpOrt->Print("V");
   //given the parameter [i], sigma = sqrt( [i] ) and the uncertainty is given by 0.5 * sigma([i]) / sqrt([i]) according to error propagation
   Double_t sigmaCOrt = TMath::Sqrt(GresolutionMumetOrtZvsNvtx->GetFunction("fitResoOrtVsNvtx")->GetParameter(0));
   Double_t sigmaCErrOrt = 0.5 * (TMath::Sqrt(GresolutionMumetOrtZvsNvtx->GetFunction("fitResoOrtVsNvtx")->GetParError(0))) / TMath::Sqrt(sigmaCOrt);
   Double_t sigmaPUOrt = TMath::Sqrt(GresolutionMumetOrtZvsNvtx->GetFunction("fitResoOrtVsNvtx")->GetParameter(1));
   Double_t sigmaPUErrOrt = 0.5 * (TMath::Sqrt(GresolutionMumetOrtZvsNvtx->GetFunction("fitResoOrtVsNvtx")->GetParError(1))) / TMath::Sqrt(sigmaPUOrt);
   gStyle->SetOptFit(10100);
   CresolutionMumetOrtZvsNvtx->SaveAs(MY_PDF_PATH"resolutionMumetOrtZvsNvtx.pdf");   
  
   cout<<"------------------------------------------------------"<<endl;
   cout<<"Fit to met_|_ resolution"<<endl;
   cout<<"sigmaC = "<<sigmaCOrt<<" +/- " <<sigmaCErrOrt<<endl;
   cout<<"sigmaPU = "<<sigmaPUOrt<<" +/- " <<sigmaPUErrOrt<<endl;
   cout<<"------------------------------------------------------"<<endl;
   cout<<"Fit to met|| resolution"<<endl;
   cout<<"sigmaC = "<<sigmaCPar<<" +/- " <<sigmaCErrPar<<endl;
   cout<<"sigmaPU = "<<sigmaPUPar<<" +/- " <<sigmaPUErrPar<<endl;
   cout<<"------------------------------------------------------"<<endl;

   //now I use t1pfmet
   // for (Int_t i = 0; i < NVTXS; i++) {
   //   yValues[i] = Ht1pfmetParZvsNvtx[i]->GetRMS();
   //   yValuesErr[i] = Ht1pfmetParZvsNvtx[i]->GetRMSError();
   // }

   // TCanvas *CresolutionT1pfmetParZvsNvtx = new TCanvas("CresolutionT1pfmetParZvsNvtx","resolution ||");
   // TGraphErrors *GresolutionT1pfmetParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   // GresolutionT1pfmetParZvsNvtx->SetTitle("resolution || from histogram's RMS");
   // GresolutionT1pfmetParZvsNvtx->Draw("AP");
   // GresolutionT1pfmetParZvsNvtx->SetMarkerStyle(21);
   // GresolutionT1pfmetParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   // GresolutionT1pfmetParZvsNvtx->GetYaxis()->SetTitle("t1pfmet|| resolution [GeV]");
   // GresolutionT1pfmetParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   // CresolutionT1pfmetParZvsNvtx->SaveAs(MY_PDF_PATH"resolutionT1pfmetParZvsNvtx.pdf");

   // for (Int_t i = 0; i < NVTXS; i++) {
   //   yValues[i] = Ht1pfmetOrtZvsNvtx[i]->GetRMS();
   //   yValuesErr[i] = Ht1pfmetOrtZvsNvtx[i]->GetRMSError();
   // }

   // TCanvas *CresolutionT1pfmetOrtZvsNvtx = new TCanvas("CresolutionT1pfmetOrtZvsNvtx","resolution _|_");
   // TGraphErrors *GresolutionT1pfmetOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   // GresolutionT1pfmetOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   // GresolutionT1pfmetOrtZvsNvtx->Draw("AP");
   // GresolutionT1pfmetOrtZvsNvtx->SetMarkerStyle(21);
   // GresolutionT1pfmetOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   // GresolutionT1pfmetOrtZvsNvtx->GetYaxis()->SetTitle("t1pfmet_|_ resolution [GeV]");
   // GresolutionT1pfmetOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   // CresolutionT1pfmetOrtZvsNvtx->SaveAs(MY_PDF_PATH"resolutionT1pfmetOrtZvsNvtx.pdf");

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = HmumetX[i]->GetMean();
     yValuesErr[i] = HmumetX[i]->GetMeanError();
   }

   TCanvas *CmumetX = new TCanvas("CmumetX","mumet along X axis");
   TGraphErrors *GmumetX = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GmumetX->SetTitle("<Met_X> vs Nvtx");
   GmumetX->Draw("AP");
   GmumetX->SetMarkerStyle(7);
   GmumetX->GetXaxis()->SetTitle("nvtx");
   GmumetX->GetYaxis()->SetTitle("<MetX> [GeV]");
   GmumetX->GetYaxis()->SetTitleOffset(1.4); 
   TF1 *fitmumetX = new TF1("fitmumetX","[0] * x + [1] ",6.,28.); 
   fitmumetX->SetParameter(0,-1.);
   fitmumetX->SetParameter(1,0.);
   TFitResultPtr frpMumetX = GmumetX->Fit(fitmumetX,"I S Q B R","");
   frpMumetX->Print("V");
   gStyle->SetOptFit(1111);
   CmumetX->SaveAs(MY_PDF_PATH"mumetX.pdf");   

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = HmumetY[i]->GetMean();
     yValuesErr[i] = HmumetY[i]->GetMeanError();
   }

   TCanvas *CmumetY = new TCanvas("CmumetY","mumet along Y axis");
   TGraphErrors *GmumetY = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GmumetY->SetTitle("<Met_Y> vs Nvtx");
   GmumetY->Draw("AP");
   GmumetY->SetMarkerStyle(7);
   GmumetY->GetXaxis()->SetTitle("nvtx");
   GmumetY->GetYaxis()->SetTitle("<MetY> [GeV]");
   GmumetY->GetYaxis()->SetTitleOffset(1.4); 
   TF1 *fitmumetY = new TF1("fitmumetY","[0] * x + [1] ",6.,28.); 
   fitmumetY->SetParameter(0,-1.);
   fitmumetY->SetParameter(1,0.);
   TFitResultPtr frpMumetY = GmumetY->Fit(fitmumetY,"I S Q B R","");
   frpMumetY->Print("V");
   gStyle->SetOptFit(1111);
   CmumetY->SaveAs(MY_PDF_PATH"mumetY.pdf");    

   /************************************/
   //                    MET RESPONSE (using mumet)
   /************************************/

   Double_t response[nBinsForResponse];
   Double_t responseErr[nBinsForResponse];
   Double_t meanZpt[nBinsForResponse];
   Double_t meanZptErr[nBinsForResponse];
   for (Int_t i = 0; i < nBinsForResponse; i++) {
     meanZpt[i] = HZptBinned[i]->GetMean();
     meanZptErr[i] = HZptBinned[i]->GetMeanError();
     response[i] = HmumetParZratio[i]->GetMean();
     responseErr[i] = HmumetParZratio[i]->GetMeanError();
     //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;
   }

   TCanvas *CresponseCurve = new TCanvas("CresponseCurve","response curve using < met|| / wzpt >");
   TGraphErrors *GresponseCurve = new TGraphErrors(nBinsForResponse,meanZpt,response,0,responseErr);
   GresponseCurve->SetTitle("response curve");
   GresponseCurve->Draw("AP");
   GresponseCurve->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve->GetXaxis()->SetTitle("Z_pt [GeV]");
   GresponseCurve->GetYaxis()->SetTitle(" < Met|| / wzpt >");
   GresponseCurve->GetYaxis()->SetRangeUser(0.5, 1.1);
   GresponseCurve->GetYaxis()->SetTitleOffset(1.4); 
   CresponseCurve->SaveAs(MY_PDF_PATH"responseCurve.pdf");

   /************************************/
   //                    MET|| & MET_|_ VS Z_PT
   /************************************/

   Double_t resoMumetParZvsZpt[nBinsForResponse];
   Double_t resoMumetParZvsZptErr[nBinsForResponse];
   Double_t resoMumetOrtZvsZpt[nBinsForResponse];
   Double_t resoMumetOrtZvsZptErr[nBinsForResponse];
   for (Int_t i = 0; i < nBinsForResponse; i++) {
     resoMumetParZvsZpt[i] = HmumetParZvsZpt[i]->GetRMS();
     resoMumetParZvsZptErr[i] = HmumetParZvsZpt[i]->GetRMSError();
     resoMumetOrtZvsZpt[i] = HmumetOrtZvsZpt[i]->GetRMS();
     resoMumetOrtZvsZptErr[i] = HmumetOrtZvsZpt[i]->GetRMSError();
   }

   TCanvas *CresolutionMumetParZvsZpt = new TCanvas("CresolutionMumetParZvsZpt","resolution ||");
   TGraphErrors *GresolutionMumetParZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMumetParZvsZpt,0,resoMumetParZvsZptErr);
   GresolutionMumetParZvsZpt->SetTitle("resolution || from histogram's RMS");
   GresolutionMumetParZvsZpt->Draw("AP");
   GresolutionMumetParZvsZpt->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMumetParZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMumetParZvsZpt->GetYaxis()->SetTitle("mumet|| resolution [GeV]");
   GresolutionMumetParZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   CresolutionMumetParZvsZpt->SaveAs(MY_PDF_PATH"resolutionMumetParZvsZpt.pdf");

   TCanvas *CresolutionMumetOrtZvsZpt = new TCanvas("CresolutionMumetOrtZvsZpt","resolution _|_");
   TGraphErrors *GresolutionMumetOrtZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMumetOrtZvsZpt,0,resoMumetOrtZvsZptErr);
   GresolutionMumetOrtZvsZpt->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMumetOrtZvsZpt->Draw("AP");
   GresolutionMumetOrtZvsZpt->SetMarkerStyle(7);
   GresolutionMumetOrtZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMumetOrtZvsZpt->GetYaxis()->SetTitle("mumet_|_ resolution [GeV]");
   GresolutionMumetOrtZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   CresolutionMumetOrtZvsZpt->SaveAs(MY_PDF_PATH"resolutionMumetOrtZvsZpt.pdf");  

   /****************************************************/
   // END OF CANVASES
   /****************************************************/ 

}      

   /****************************************************/
   // END OF ANALYSIS.LOOP()
   /****************************************************/ 




// END OF PROGRAMME.

