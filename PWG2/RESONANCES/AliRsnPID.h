/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id: AliRsnPID.h,v 1.5 2007/02/21 14:33:25 pulvir Exp $ */

//-------------------------------------------------------------------------
//                      Class AliRsnPID
//  Simple collection of reconstructed tracks, selected from an ESD event
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPID_H
#define ALIRSNPID_H

#include <TNamed.h>

class AliRsnDaughter;
class AliRsnEvent;

class AliRsnPID : public TNamed
{
public:

    // types enum
    enum EType {
        kElectron = 0,
        kMuon,
        kPion,
        kKaon,
        kProton,
        kUnknown,
        kSpecies = 5
    };

    AliRsnPID();
    virtual ~AliRsnPID() {}

    // conversions from PDG code to local type
    static EType        InternalType(Int_t pdgCode);
    
    // retrieve particle informations from internal type
    static Int_t        PDGCode(EType pid);
    static const char*  ParticleName(EType pid, Bool_t shortName = kTRUE);
    static const char*  ParticleNameLatex(EType pid);
    static Double_t     ParticleMass(EType pid);

    // identification routines
    Bool_t Process(AliRsnEvent *e);
    Bool_t ComputeProbs(AliRsnDaughter *d);
    Bool_t IdentifiedAs(AliRsnDaughter *d, EType type, Short_t charge = 0);
    EType  TrackType(AliRsnDaughter *d);

    // data members
    void     SetPriorProbability(EType type, Double_t p);
    void     SetMinProb(Double_t p) {fMinProb = p;}
    void     SetMaxPt(Double_t p) {fMaxPt = p;}
    Double_t GetPriorProbability(EType type) {return fPrior[(Int_t)type];}
    Double_t GetMinProb() {return fMinProb;}
    Double_t GetMaxPt() {return fMaxPt;}

    // other
    void DumpPriors();

private:

    Double_t  fPrior[kSpecies]; // prior probabilities
    Double_t  fMaxPt;           // pt threshold for realistic PID
    Double_t  fMinProb;         // threshold on acceptable largest probability

    static const Double_t  fgkParticleMass[kSpecies + 1];      // PDG particle mass
    static const char*     fgkParticleNameShort[kSpecies + 1]; // short particle name
    static const char*     fgkParticleNameLong[kSpecies + 1];  // long particle name
    static const char*     fgkParticleNameLatex[kSpecies + 1]; // latex particle name
    static const Int_t     fgkParticlePDG[kSpecies + 1];       // PDG code of particle

    ClassDef(AliRsnPID,1);
};

#endif
