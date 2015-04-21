
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>

#include "Analysis.h"

using namespace std;



int main(int argc, char* argv[]) {


      //================ Creating chain 
 

  TFile* fin= new TFile("~/Documents/myRootWorkDir/provevarie/treeFromCmsrm/tree.root");
  
  TChain* chain =new TChain("tree/tree");

  chain->Add("~/Documents/myRootWorkDir/provevarie/treeFromCmsrm/tree.root/tree/tree");  


  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis
  Analysis tree( chain );
  tree.Loop();
  
  delete fin;
  delete chain;
  return 0;
}
