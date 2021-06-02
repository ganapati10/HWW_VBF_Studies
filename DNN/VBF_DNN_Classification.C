//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  root -l VBF_TMVAClassification.C\(\"DNN\"\)
//  root -l VBF_TMVAClassification.C\(\"DNN\"\)
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjString.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Tools.h"

// User defined function
//#include "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/mlljj.C" 
//#include "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/mlj.C"
//#include "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/RecoMELA_VBF.C"
//#include "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/RecoMELA_Phi.C"


void VBF_DNN_Classification(TString myMethodList = "") 
{
  // Load the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;  

  gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so");
  gSystem->Load("libJHUGenMELAMELA.so");

  gROOT->ProcessLineSync(".L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/RecoMELA_VBF.C+");
  gROOT->ProcessLineSync(".L mlljj.C+");
  //gROOT->ProcessLineSync(".L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v7/TMVA/RecoMELA_Phi.C+");
    

  Use["DNN"]   = 0;  // Uses Adaptive Boost

  Use["BDT"] = 0;

  if (myMethodList != "") {

    for (std::map<std::string,int>::iterator it=Use.begin(); it!=Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');

    for (UInt_t i=0; i<mlist.size(); i++) {

      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {

	std::cout << " Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;

	for (std::map<std::string,int>::iterator it=Use.begin(); it!=Use.end(); it++) std::cout << it->first << " ";

	std::cout << std::endl;

	return;
      }

      Use[regMethod] = 1;
    }
  }


  // Input and output files
  //----------------------------------------------------------------------------
  TString workdir16 = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer16_102X_nAODv7_Full2016v7/MCl1loose2016v7__MCCorr2016v7__l2loose__l2tightOR2016v7/";
  TString workdir17 = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Fall2017_102X_nAODv7_Full2017v7/MCl1loose2017v7__MCCorr2017v7__l2loose__l2tightOR2017v7/";
  TString workdir18 = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Autumn18_102X_nAODv7_Full2018v7/MCl1loose2018v7__MCCorr2018v7__l2loose__l2tightOR2018v7/";

  TString outfileName("VBF_DNN_TMVAClassification.root");

  TFile* outputFile = TFile::Open(outfileName, "recreate");


  // Create the factory object. The first argument is the base of the name of all the weight files
  //----------------------------------------------------------------------------
  TString factoryName("VBF_DNN_TMVAClassification");

  //TMVA::Factory* factory = new TMVA::Factory(factoryName, outputFile,
  //					     "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification",
					      outputFile,
					      "AnalysisType=Classification" );

  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

  dataloader->AddVariable("mjj", 'F');
  //dataloader->AddVariable("log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)", 'F');
  dataloader->AddVariable("Jet_qgl[0]", 'F');
  dataloader->AddVariable("Jet_qgl[1]", 'F');
  dataloader->AddVariable("detajj", 'F');
  dataloader->AddVariable("Lepton_eta[0]-Lepton_eta[1]", 'F');
  dataloader->AddVariable("sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])", 'F');
  dataloader->AddVariable("dphill", 'F');
  dataloader->AddVariable("dphijjmet", 'F');
  dataloader->AddVariable("dphilljetjet", 'F');
  dataloader->AddVariable("drll", 'F');
  dataloader->AddVariable("Lepton_eta[0]", 'F');
  dataloader->AddVariable("Lepton_eta[1]", 'F');
  dataloader->AddVariable("Lepton_pt[0]", 'F');
  dataloader->AddVariable("Lepton_pt[1]", 'F');
  dataloader->AddVariable("Lepton_phi[0]", 'F');
  dataloader->AddVariable("Lepton_phi[1]", 'F');
  dataloader->AddVariable("CleanJet_eta[0]", 'F');
  dataloader->AddVariable("CleanJet_eta[1]", 'F');
  dataloader->AddVariable("CleanJet_phi[0]", 'F');
  dataloader->AddVariable("CleanJet_phi[1]", 'F');
  dataloader->AddVariable("CleanJet_pt[0]", 'F');
  dataloader->AddVariable("CleanJet_pt[1]", 'F');
  dataloader->AddVariable("PuppiMET_pt", 'F');
  dataloader->AddVariable("PuppiMET_phi", 'F');
  dataloader->AddVariable("mth", 'F');
  dataloader->AddVariable("mTi", 'F');
  dataloader->AddVariable("mtw2", 'F');
  dataloader->AddVariable("Lepton_eta[0]-CleanJet_eta[0]", 'F');
  dataloader->AddVariable("Lepton_eta[0]-CleanJet_eta[1]", 'F');
  dataloader->AddVariable("Lepton_eta[1]-CleanJet_eta[0]", 'F');
  dataloader->AddVariable("Lepton_eta[1]-CleanJet_eta[1]", 'F');
  dataloader->AddVariable("ptll", 'F');
  //dataloader->AddVariable("mlljj:=mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])", 'F');
  //dataloader->AddVariable("ml1j1:=mlj(CleanJet_pt[0], CleanJet_phi[0], CleanJet_eta[0], Lepton_pt[0],  Lepton_phi[0], Lepton_eta[0])", 'F');
  //dataloader->AddVariable("ml2j1:=mlj(CleanJet_pt[0], CleanJet_phi[0], CleanJet_eta[0], Lepton_pt[1],  Lepton_phi[1], Lepton_eta[1])", 'F');
  //dataloader->AddVariable("ml1j2:=mlj(CleanJet_pt[1], CleanJet_phi[1], CleanJet_eta[1], Lepton_pt[0],  Lepton_phi[0], Lepton_eta[0])", 'F');
  //dataloader->AddVariable("ml2j2:=mlj(CleanJet_pt[1], CleanJet_phi[1], CleanJet_eta[1], Lepton_pt[1],  Lepton_phi[1], Lepton_eta[1])", 'F');
  //dataloader->AddVariable("Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 'F');
  //dataloader->AddVariable("mll", 'F');
  //dataloader->AddVariable("RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])", 'F');
  

  // Input files
  //----------------------------------------------------------------------------
  std::vector<TFile*> InputFiles_signal;
  std::vector<TFile*> InputFiles_background;

  InputFiles_signal.clear();
  InputFiles_background.clear();

  InputFiles_signal.push_back(TFile::Open(Form("%s/nanoLatino_VBFHToWWTo2L2Nu_M125__part0.root", workdir16.Data())));
  InputFiles_signal.push_back(TFile::Open(Form("%s/nanoLatino_VBFHToWWTo2L2Nu_M125__part1.root", workdir16.Data())));

  for (UInt_t k=0; k<22; k++){
    InputFiles_signal.push_back(TFile::Open(Form("%s/nanoLatino_VBFHToWWTo2L2Nu_M125__part%d.root", workdir17.Data(), k)));
  }

  for (UInt_t k=0; k<12; k++){
    InputFiles_signal.push_back(TFile::Open(Form("%s/nanoLatino_VBFHToWWTo2L2Nu_M125__part%d.root", workdir18.Data(), k)));
  }

  // GGH

  for (UInt_t k=0; k<21; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part%d.root", workdir16.Data(), k)));
  }

  for (UInt_t k=0; k<37; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GluGluHToWWTo2L2Nu_M125__part%d.root", workdir17.Data(), k)));
  }

  for (UInt_t k=0; k<44; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part%d.root", workdir17.Data(), k)));
  }

  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GluGluHToWWTo2L2Nu_M125__part0.root", workdir18.Data())));
  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GluGluHToWWTo2L2Nu_M125__part1.root", workdir18.Data())));

  for (UInt_t k=0; k<57; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part%d.root", workdir18.Data(), k)));
  }


  //TOP

  for (UInt_t k=0; k<72; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_TTTo2L2Nu__part%d.root", workdir16.Data(), k)));
  }

  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_s-channel__part0.root", workdir16.Data())));
  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_s-channel__part1.root", workdir16.Data())));
  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_s-channel__part2.root", workdir16.Data())));

  //for (UInt_t k=0; k<35; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_antitop__part%d.root", workdir16.Data(), k)));
    //}

  //for (UInt_t k=0; k<143; k++){
    // InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_top__part%d.root", workdir16.Data(), k)));
    //}

  //for (UInt_t k=0; k<9; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_antitop__part%d.root", workdir16.Data(), k)));
    //}

  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part0.root", workdir16.Data())));
  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part1.root", workdir16.Data())));
  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part2.root", workdir16.Data())));
  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part3.root", workdir16.Data())));
  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part4.root", workdir16.Data())));
  //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part5.root", workdir16.Data())));
  ///InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part6.root", workdir16.Data())));


  for (UInt_t k=0; k<75; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_TTTo2L2Nu_PSWeights__part%d.root", workdir17.Data(), k)));
  }

  //for (UInt_t k=0; k<10; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_s-channel__part%d.root", workdir17.Data(), k)));
    //}

  //for (UInt_t k=0; k<17; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_antitop__part%d.root", workdir17.Data(), k)));
    //}

  //for (UInt_t k=0; k<7; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_top__part%d.root", workdir17.Data(), k)));
    //}

  //for (UInt_t k=0; k<11; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_antitop__part%d.root", workdir17.Data(), k)));
    //}

  //for (UInt_t k=0; k<7; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top__part%d.root", workdir17.Data(), k)));
    //}

  for (UInt_t k=0; k<59; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_TTTo2L2Nu__part%d.root", workdir18.Data(), k)));
  }

  //for (UInt_t k=0; k<19; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_s-channel_ext1__part%d.root", workdir18.Data(), k)));
    //}

  //for (UInt_t k=0; k<60; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_antitop__part%d.root", workdir18.Data(), k)));
    //}

  //for (UInt_t k=0; k<102; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_t-channel_top__part%d.root", workdir18.Data(), k)));
    //}

  //for (UInt_t k=0; k<8; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_antitop_ext1__part%d.root", workdir18.Data(), k)));
    //}

  //for (UInt_t k=0; k<9; k++){
    //InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_ST_tW_top_ext1__part%d.root", workdir18.Data(), k)));
    //}


  // WW

  
  for (UInt_t k=0; k<4; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WWTo2L2Nu__part%d.root", workdir16.Data(), k)));
  }

  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WpWmJJ_EWK__part0.root", workdir16.Data())));
  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WpWmJJ_EWK__part1.root", workdir16.Data())));

  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WWTo2L2Nu__part0.root", workdir17.Data())));
  InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WWTo2L2Nu__part1.root", workdir17.Data())));

  for (UInt_t k=0; k<16; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WpWmJJ_EWK__part%d.root", workdir17.Data(), k)));
  }

  for (UInt_t k=0; k<10; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WWTo2L2Nu__part%d.root", workdir18.Data(), k)));
  }

  for (UInt_t k=0; k<30; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_WpWmJJ_EWK__part%d.root", workdir18.Data(), k)));
  }


  // DY


  for (UInt_t k=0; k<21; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part%d.root", workdir16.Data(), k)));
  }

  for (UInt_t k=0; k<26; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToLL_M-10to50-LO__part%d.root", workdir16.Data(), k)));
  }

  for (UInt_t k=0; k<73; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToTT_MuEle_M-50__part%d.root", workdir17.Data(), k)));
  }

  for (UInt_t k=0; k<76; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part%d.root", workdir17.Data(), k)));
  }

  for (UInt_t k=0; k<101; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToTT_MuEle_M-50__part%d.root", workdir18.Data(), k)));
  }

  for (UInt_t k=0; k<78; k++){
    InputFiles_background.push_back(TFile::Open(Form("%s/nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part%d.root", workdir18.Data(), k)));
  }





  // Apply cuts on the signal and background samples (can be different)
  //----------------------------------------------------------------------------
  TCut mycut = "Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) &&  Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && (mth > 40 && mth < 125) && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3 && RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1]) > 0.8";

  TCut mycuts = "";

  TCut test = "mth > 40"; 
  
  // Create factory for signal and background samples
  double tmpWeight = 1.;

  for (UInt_t i=0; i<InputFiles_signal.size(); ++i) {
    TTree* tmpsTree = (TTree*)InputFiles_signal.at(i)->Get("Events");
    dataloader->AddSignalTree(tmpsTree, tmpWeight);
  }

  for (UInt_t k=0; k<InputFiles_background.size(); ++k) {
    TTree* tmpbTree = (TTree*)InputFiles_background.at(k)->Get("Events");
    dataloader->AddBackgroundTree(tmpbTree, tmpWeight);
  }


  //TString dataString = "nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=1000:nTest_Background=1000:SplitMode=Random:NormMode=NumEvents:!V";

  dataloader->PrepareTrainingAndTestTree("", "", "nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=1000:nTest_Background=1000:SplitMode=Random::SplitSeed=10:NormMode=EqualNumEvents");

//dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random::SplitSeed=10:NormMode=None:!V");
//dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=2404:nTrain_Background=11671:SplitMode=Block::SplitSeed=10:NormMode=EqualNumEvents");
//dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=Alternate::SplitSeed=10:NormMode=EqualNumEvents");
//dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=Random:NormMode=NumEvents:!V");
//dataloader->PrepareTrainingAndTestTree(mycut, mycut, "SplitMode=Random::SplitSeed=10:NormMode=EqualNumEvents");

  
  TString configString = "!H:V";
  configString += ":VarTransform=N";

  configString += ":ErrorStrategy=CROSSENTROPY";

  configString += ":WeightInitialization=XAVIERUNIFORM";

  TString inputLayoutString = "InputLayout=1|1|41";

  TString layoutString = "Layout=RELU|70,RELU|50,RELU|30,RELU|30,RELU|20,RELU|10,SIGMOID";

  TString trainingString1 = "TrainingStrategy="
    "LearningRate=0.001,"
    "Repetitions=1,"
    "ConvergenceSteps=10,"
    "MaxEpoch=200,"
    "Optimizer=RMSPROP,"
    "BatchSize=512,"
    "WeightDecay=0.001,"
    "Multithreading=True";

  TString trainingString = "TrainingStrategy=LearningRate=1e-1,MaxEpoch=500";

  configString += ":" + inputLayoutString + ":" + layoutString + ":" + trainingString1;

  // Book MVA methods
  //----------------------------------------------------------------------------

  printf("Start configuration of the DNN\n");
  if (Use["DNN"]) factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_GPU", configString);
  printf("Start training \n");
  

  if (Use["BDT"]) factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
				      "!H:!V:NTrees=250:MinNodeSize=0.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=20");
  
  // Now you can tell the factory to train, test, and evaluate the MVAs
  //----------------------------------------------------------------------------
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();    


  // Save the output
  //----------------------------------------------------------------------------
  outputFile->Close();

  delete factory;
  delete dataloader;

  
  // Launch the GUI for the root macros
  //----------------------------------------------------------------------------
  if (!gROOT->IsBatch()) TMVA::TMVAGui(outfileName);
}
