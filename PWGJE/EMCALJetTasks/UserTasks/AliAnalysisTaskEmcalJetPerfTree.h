#ifndef ALIANALYSISTASKEMCALJETPERFTREE_H
#define ALIANALYSISTASKEMCALJETPERFTREE_H
/**
 * \file AliAnalysisTaskEmcalJetPerfTree.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetPerfTree
 *
 * In this header file the class AliAnalysisTaskEmcalJetPerfTree is declared.
 * This is  task to extract info about an MC jet to fill a response (also fills response and JES)
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

class AliJetContainer;
class AliEmcalJet;
class AliVCaloCells;




/**
 * \class AliAnalysisTaskEmcalJetPerfTree
 * \brief Implementation of a sample jet analysis task.
 * Classs to evaluate the performance jets
 * It derives from AliAnalysisTaskEmcalJet.
 */

//namespace PWGJE { namespace EMCALJetTasks { class AliAnalysisTaskEmcalJetHPerformance; } }

namespace PWGJE {
  namespace EMCALJetTasks {



class AliAnalysisTaskEmcalJetPerfTree : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetPerfTree()                                               ;
  AliAnalysisTaskEmcalJetPerfTree(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalJetPerfTree()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskEmcalJetPerfTree* AddTaskEmcalJetPerfTree( const char *suffix             = "");

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateJetTree()                                 ;

  void                        DoJetLoop()                                       ;
  void                        DoResponse()                                      ;

  void                        FillResponseMatrix(AliEmcalJet * jet1, AliEmcalJet * jet2)                              ;
  void                        FillResponseTree(AliEmcalJet * jet1, AliEmcalJet * jet2)                                ;

  THistManager                fHistManager                                      ;///< Histogram manager
  TTree   *fTree;    //!<! Output tree                                                                                                           

 private:
  AliAnalysisTaskEmcalJetPerfTree(const AliAnalysisTaskEmcalJetPerfTree&)           ; // not implemented
  AliAnalysisTaskEmcalJetPerfTree &operator=(const AliAnalysisTaskEmcalJetPerfTree&); // not implemented

  bool fPartLevelResp; ///<  If true, create the response between hybrid and particle level and if false create response between hybrid and detector level
  double fMinFractionShared;             ///<  Minimum fraction of shared jet pt required for matching a hybrid jet to detector level     
  bool fCreateTree; ///< If true, also fill a tree with MC variables in addition to the response and JES

  Float_t         fJet1_Pt;                        //!<! array buffer                                                                      
  Float_t         fJet1_Eta;                       //!<! array buffer                                                                      
  Float_t         fJet1_Phi;                       //!<! array buffer                                                                      
  Float_t         fJet1_Area;                      //!<! array buffer  
  Float_t         fJet2_Pt;                        //!<! array buffer                                                                      
  Float_t         fJet2_Eta;                       //!<! array buffer                                                                      
  Float_t         fJet2_Phi;                       //!<! array buffer                                                                      
  Float_t         fJet2_Area;                      //!<! array buffer  
  Float_t         fEvent_BackgroundDensity;      //!<! array buffer                                                                      
  Float_t         fEvent_Vertex_X;               //!<! array buffer                                                                      
  Float_t         fEvent_Vertex_Y;               //!<! array buffer                                                                      
  Float_t         fEvent_Vertex_Z;               //!<! array buffer                                                                      
  Float_t         fEvent_Centrality;             //!<! array buffer                                                                      
  Int_t           fEvent_Multiplicity;           //!<! array buffer                                                                      
  Long64_t        fEvent_ID;                     //!<! array buffer                                                                      
  Float_t         fEvent_PtHard;                 //!<! array buffer                                                                      
  Float_t         fEvent_Weight;                 //!<! array buffer                                                                      
  Float_t         fEvent_ImpactParameter;        //!<! array buffer         

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetPerfTree, 1);
  /// \endcond
};

  }/* namespace EMCALJetTasks */
} /* namespace PWGJE */
#endif
