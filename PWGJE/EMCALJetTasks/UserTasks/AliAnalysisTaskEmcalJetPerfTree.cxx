/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskEmcalJetPerfTree.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetPerfTree);
/// \endcond

namespace PWGJE {
  namespace EMCALJetTasks {

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetPerfTree::AliAnalysisTaskEmcalJetPerfTree() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
  fPartLevelResp(true),
  fMinFractionShared(0.),
  fCreateTree(false),
  fMultiplicity(0.)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetPerfTree::AliAnalysisTaskEmcalJetPerfTree(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
  fPartLevelResp(true),
  fMinFractionShared(0.),
  fCreateTree(false),
  fMultiplicity(0.)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetPerfTree::~AliAnalysisTaskEmcalJetPerfTree()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetPerfTree::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateJetHistograms();
  if (fCreateTree) AllocateJetTree();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  PostData(2, fTree); // Post data for ALL output slots > 0 here.
}


/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetPerfTree::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
      }

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      }
    }
  }
  //now initialize the performance histograms for the matched jets 
  groupname = "JetPerformance";
  fHistManager.CreateHistoGroup(groupname);
  for (Int_t cent = 0; cent < fNcentBins; cent++) {
    histname = TString::Format("%s/histJetPtresp_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;#it{p}_{T,jet1} (GeV/#it{c});#it{p}_{T,jet2} (GeV/#it{c}); JES; counts", histname.Data());
    fHistManager.CreateTH3(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt, 200, 0, 5);
  }
}

  void AliAnalysisTaskEmcalJetPerfTree::AllocateJetTree()
  {
    //fill this in once I figure out how to do trees 
    //make sure that the info for looking at the jet overlap is in here
    //also make sure the info about the pThard bins is in here

    fTree = new TTree("fTree", "A tree of variables for matched jets that fill a response");
    fTree->Branch("Jet1_Pt", &fJet1_Pt);
    fTree->Branch("Jet1_Phi", &fJet1_Phi);
    fTree->Branch("Jet1_Eta", &fJet1_Eta);
    fTree->Branch("Jet1_Area", &fJet1_Area);
    fTree->Branch("Jet2_Pt", &fJet2_Pt);
    fTree->Branch("Jet2_Phi", &fJet2_Phi);
    fTree->Branch("Jet2_Eta", &fJet2_Eta);
    fTree->Branch("Jet2_Area", &fJet2_Area);

    fTree->Branch("Event_BackgroundDensity",&fEvent_BackgroundDensity);
    fTree->Branch("Event_Vertex_X",&fEvent_Vertex_X);
    fTree->Branch("Event_Vertex_Y",&fEvent_Vertex_Y);
    fTree->Branch("Event_Vertex_Z",&fEvent_Vertex_Z);
    fTree->Branch("Event_Centrality",&fEvent_Centrality);
    fTree->Branch("Event_Multiplicity",&fEvent_Multiplicity);
    fTree->Branch("Event_ID",&fEvent_ID);

    fTree->Branch("Event_PtHard",&fEvent_PtHard);
    fTree->Branch("Event_Weight",&fEvent_Weight);
    fTree->Branch("Event_ImpactParameter",&fEvent_ImpactParameter);
  }

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetPerfTree::FillHistograms()
{
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetPerfTree::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }


}


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetPerfTree::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetPerfTree::Run()
{
  SetGlobalVariables();
  DoResponse();
  return kTRUE;
}

/**
 *make sure all the global variables are set to fill the tree later
 */
  
void AliAnalysisTaskEmcalJetPerfTree::SetGlobalVariables()
{
  feventID = 0;
  const AliVVertex* myVertex = InputEvent()->GetPrimaryVertex();
  if(!myVertex && MCEvent())
    myVertex = MCEvent()->GetPrimaryVertex();
  Double_t vtx[3] = {0, 0, 0};
  if(myVertex)
    {
      vtx[0] = myVertex->GetX(); vtx[1] = myVertex->GetY(); vtx[2] = myVertex->GetZ();
    }
  fvtx_X                        = vtx ? vtx[0] : 0;
  fvtx_Y                        = vtx ? vtx[1] : 0;
  fvtx_Z                        = vtx ? vtx[2] : 0;


  AliVHeader* eventIDHeader = InputEvent()->GetHeader();
  if(eventIDHeader)
    feventID = eventIDHeader->GetEventIdAsLong();


  fMultiplicity = 0;
  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++)
    fMultiplicity += GetParticleContainer(iCont)->GetNAcceptedParticles();
  
  frho = GetJetContainer(0)->GetRhoVal();

}

/**
 *separate function to find all the matching for the jets and fill the response histogram and tree 
 */

 void AliAnalysisTaskEmcalJetPerfTree::DoResponse()
 {
   TString histname;
   TString groupname;

   AliJetContainer * jetsHybrid = GetJetContainer("hybridLevelJets");
   AliJetContainer * jetsDetLevel = GetJetContainer("detLevelJets");
   AliJetContainer * jetsPartLevel = GetJetContainer("partLevelJets");
   if (!jetsHybrid) {
     AliErrorStream() << "Could not retrieve hybrid jet collection.\n";
     return;
   }
   if (!jetsDetLevel) {
     AliErrorStream() << "Could not retrieve det level jet collection.\n";
     return;
   }
   if (fPartLevelResp && !jetsPartLevel) {
     AliErrorStream() << "Could not retrieve part level jet collection.\n";
     return;
   }
   
   // Handle matching of jets.                                                                                             
   for (auto jet1 : jetsHybrid->accepted())
     {
       AliDebugStream(4) << "jet1: " << jet1->toString() << "\n";
       AliDebugStream(4) << "jet1 address: " << jet1 << "\n";

       // Get jet the det level jet from the hybrid jet                                                                  
       AliEmcalJet * jet2 = jet1->ClosestJet();
       if(!jet2) continue;

       AliDebugStream(4) << "jet2: " << jet2->toString() << "\n";
       AliDebugStream(4) << "jet2 address: " << jet2 << "\n";

       // Check shared fraction                                                                                           
       double sharedFraction = jetsHybrid->GetFractionSharedPt(jet1);
       //       fHistManager.FillTH1("fHistFractionSharedPt", sharedFraction);
       if (sharedFraction < fMinFractionShared) {
	 AliDebugStream(4) << "Rejecting jet due to momentum fraction of " << sharedFraction << ", smaller than the minimum.\n";
	 continue;
       }
       else {
	 AliDebugStream(4) << "Jet passed momentum fraction cut with value of " << sharedFraction << "\n";
       }
     
       // Apply additional selection to jet 2                                                                                     

       // Get MC level jet                                                                                                      
       AliEmcalJet * jetToPass = 0;
       if (fPartLevelResp) {
	 AliEmcalJet * jet3 = jet2->ClosestJet();
     
	 // Accept jet 3                                                                                                              
	 UInt_t rejectionReason = 0;
	 if (!jetsPartLevel->AcceptJet(jet3, rejectionReason)) {
	   // TODO: Store rejection reasons                                                                                     
	   //fHistRejectionReason2->Fill(jets2->GetRejectionReasonBitPosition(rejectionReason), jet2->Pt());                    
	   continue;
	 }
     
	 AliDebugStream(4) << "jet3: " << jet3->toString() << "\n";
	 AliDebugStream(4) << "jet3 address: " << jet3 << "\n";
     
	 // Use for the response                                                                                               
	 AliDebugStream(4) << "Using part level jet for response\n";
	 jetToPass = jet3;
       }
       else {
	 // Use for the response                                                                                              
	 AliDebugStream(4) << "Using det level jet for response\n";
	 jetToPass = jet2;
       }
       
       // Fill response                                                                                                       
       FillResponseMatrix(jet1, jetToPass);
       //Fill tree
       FillResponseTree(jet1, jetToPass);
     }
   //what about the not matched truth jets? can we do something with those?
 }

 void AliAnalysisTaskEmcalJetPerfTree::FillResponseMatrix(AliEmcalJet * jet1, AliEmcalJet * jet2)
 {

   if (!jet1 || !jet2) {
     AliErrorStream() << "Null jet passed to fill response matrix";
   }
   
   //fill response histogram
   TString histname;
   TString groupname;

   groupname = "JetPerformance";
   histname = TString::Format("%s/histJetPtresp_%d", groupname.Data(), fCentBin);
   fHistManager.FillTH3(histname, jet2->Pt(), jet1->Pt(), jet1->Pt()/jet2->Pt());
 }

 void AliAnalysisTaskEmcalJetPerfTree::FillResponseTree(AliEmcalJet * jet1, AliEmcalJet * jet2)
 {
   //fill the tree once I figure out how to do this

   Jet1_Pt                                   = jet1->Pt() - rho*jet1->Area();
   Jet1_Phi                                   = jet1->Phi();
   Jet1_Eta                                   = jet1->Eta();
   Jet1_Area                                  = jet1->Area();
   Jet2_Pt                                   = jet2->Pt() - rho*jet2->Area();
   Jet2_Phi                                   = jet2->Phi();
   Jet2_Eta                                   = jet2->Eta();
   Jet2_Area                                  = jet2->Area();

   
   // Set event properties                                                                                                                
   fBuffer_Event_BackgroundDensity               = rho;
   fBuffer_Event_Vertex_X                        = vertex ? vertex[0] : 0;
   fBuffer_Event_Vertex_Y                        = vertex ? vertex[1] : 0;
   fBuffer_Event_Vertex_Z                        = vertex ? vertex[2] : 0;
   fBuffer_Event_Centrality                      = fCent;
   fBuffer_Event_Multiplicity                    = fMultiplicity;
   fBuffer_Event_ID                              = eventID;
 }


/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetPerfTree::Terminate(Option_t *) 
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
 AliAnalysisTaskEmcalJetPerfTree * AliAnalysisTaskEmcalJetPerfTree::AddTaskEmcalJetPerfTree(const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetPerfTree", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetPerfTree", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  /*  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
    }*/

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  // Setup task name                                                                                                               
  std::string taskName = "AliAnalysisTaskEmcalJetPerfTree";
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  AliAnalysisTaskEmcalJetPerfTree* task = new AliAnalysisTaskEmcalJetPerfTree(taskName.c_str());
  //  sampleTask->SetCaloCellsName(cellName);
  task->SetVzRange(-10,10);
  // Set a few general default.                                                                                                    
  task->SetNCentBins(5);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(taskName);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );

  return task;
}
  }
}
