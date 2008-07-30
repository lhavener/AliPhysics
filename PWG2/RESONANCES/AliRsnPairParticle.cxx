/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//
// Class AliRsnPairParticle
//
// Implementation of a pair of tracks, for several purposes
// - computing the total 4-momentum & inv. mass for output histos filling
// - evaluating cut checks on the pair of particles
//
// author: Martin Vala (martin.vala@cern.ch)
// revised by: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliRsnPairParticle.h"

ClassImp (AliRsnPairParticle)

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle()
{
//
// Constructor.
// Initializes all variables to meaningless values.
//

    Int_t i, j;

    for (i = 0; i < 3; i++) {
        fPTot[i] = 0.0;
        fPTotMC[i] = 0.0;
        if (i < 2) {
            fMotherLabel[i] = -1;
            fMotherPDG[i] = 0;
            fDaughter[i] = 0x0;
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = 0.0;
            fPTrackMC[j][i] = 0.0;
        }
    }
}

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle(const AliRsnPairParticle &obj) :
  TObject(obj)
{
//
// Copy constructor.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

    Int_t i, j;
    for (i = 0; i < 3; i++) {
        fPTot[i] = obj.fPTot[i];
        fPTotMC[i] = obj.fPTotMC[i];
        if (i < 2) {
            fMotherLabel[i] = obj.fMotherLabel[i];
            fMotherPDG[i] = obj.fMotherPDG[i];
            fDaughter[i] = obj.fDaughter[i];
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = obj.fPTrack[j][i];
            fPTrackMC[j][i] = obj.fPTrackMC[j][i];
        }
    }
}

//_____________________________________________________________________________
AliRsnPairParticle& AliRsnPairParticle::operator=(const AliRsnPairParticle &obj)
{
//
// Assignment operator.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

    Int_t i, j;
    for (i = 0; i < 3; i++) {
        fPTot[i] = obj.fPTot[i];
        fPTotMC[i] = obj.fPTotMC[i];
        if (i < 2) {
            fMotherLabel[i] = obj.fMotherLabel[i];
            fMotherPDG[i] = obj.fMotherPDG[i];
            fDaughter[i] = obj.fDaughter[i];
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = obj.fPTrack[j][i];
            fPTrackMC[j][i] = obj.fPTrackMC[j][i];
        }
    }

    return (*this);
}

//_____________________________________________________________________________
AliRsnPairParticle::~AliRsnPairParticle()
{
//
// Desctructor.
// Does nothing.
//
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMass(Double_t mass0, Double_t mass1)
{
//
// Compute invariant mass using reconstructed values.
// Mass in argument #1 is assigned to first track in the pair (fDaughter[0]),
// mass in argument #2 is assigned to second track in the pair (fDaughter[1])
// Then, the invariant mass of the pair is computed by using their total momentum
// and the sum of their energies as they result from assigned masses.
//

    if (!fDaughter[0] || !fDaughter[1]) {
        AliError("One of the two tracks is NULL in this pair!");
        return -1000.0;
    }

    // compute track energies using the shortcut method defined in AliRsnDaughter
    Double_t etot = 0.0;
    etot += fDaughter[0]->E(mass0);
    etot += fDaughter[1]->E(mass1);

    // compute & return invariant mass
    return  TMath::Sqrt (etot * etot - GetP2());
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMassMC(Double_t mass0, Double_t mass1)
{
//
// Compute invariant mass using MC values.
// Mass in argument #1 is assigned to first track in the pair (fDaughter[0]),
// mass in argument #2 is assigned to second track in the pair (fDaughter[1])
// Then, the invariant mass of the pair is computed by using their total momentum
// and the sum of their energies as they result from assigned masses.
//

    if (!fDaughter[0] || !fDaughter[1]) {
        AliError("One of the two tracks is NULL in this pair!");
        return -1000.0;
    }
    if (!fDaughter[0]->GetMCInfo() || !fDaughter[1]->GetMCInfo()) {
        AliError("One of the two tracks has a NULL MCInfo in this pair!");
        return -1000.0;
    }

    // compute track energies using the shortcut method defined in AliRsnDaughter
    Double_t etot = 0.0;
    etot += fDaughter[0]->GetMCInfo()->E(mass0);
    etot += fDaughter[1]->GetMCInfo()->E(mass1);

    // compute & return invariant mass
    return  TMath::Sqrt (etot * etot - GetP2());
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetAngle() const
{
//
// Returns the relative angle between the vector momenta of the tracks
// Return value is in DEGREES.
//

    Double_t dotProd = 0.0;
    dotProd += fDaughter[0]->Px() * fDaughter[1]->Px();
    dotProd += fDaughter[0]->Py() * fDaughter[1]->Pz();
    dotProd += fDaughter[0]->Pz() * fDaughter[1]->Pz();
    
    Double_t cosAngle = dotProd / (fDaughter[0]->P() * fDaughter[1]->P());
    
    return TMath::ACos(cosAngle) * TMath::RadToDeg();
}

//_____________________________________________________________________________
Bool_t AliRsnPairParticle::IsTruePair(Int_t refPDG)
{
//
// Checks if the two tracks in the pair come from the same resonance.
// This can be known if MC info is present, looking at the GEANT label of mother
// (which should be the same).
// If the argument is 0, the answer is kTRUE whenever the labels of mothers of the
// two tracks is the same. When the argument is not zero, the answer is kTRUE only
// if the mother is the same and its PDG code is equal to the argument.
//

    // if MC info is not available, the pairs is not true by default
    if (!fDaughter[0]->GetMCInfo() || !fDaughter[1]->GetMCInfo()) {
        return kFALSE;
    }

    // check that labels are the same
    if (fDaughter[0]->GetMCInfo()->Mother() != fDaughter[1]->GetMCInfo()->Mother()) {
        return kFALSE;
    }

    // if we reach this point, the two tracks have the same mother
    // let's check now the PDG code of this common mother
    Int_t motherPDG = TMath::Abs(fDaughter[0]->GetMCInfo()->MotherPDG());
    if (refPDG == 0) return kTRUE;
    else return (motherPDG == refPDG);
}

//_____________________________________________________________________________
void AliRsnPairParticle::SetPair(AliRsnDaughter *daughter1, AliRsnDaughter *daughter2)
{
//
// Accepts two AliRsnDaughter's which are the two tracks in the pair,
// fills all data-members which contain their momenta & info,
// and computes the total momentum for REC data and MC if available
//

    Int_t i;

    fDaughter[0] = daughter1;
    fDaughter[1] = daughter2;

    // copy MC info (if available)
    if (fDaughter[0]->GetMCInfo() && fDaughter[1]->GetMCInfo()) {
        for (i = 0; i < 2; i++) {
            fPTrackMC[i][0] = fDaughter[i]->GetMCInfo()->Px();
            fPTrackMC[i][1] = fDaughter[i]->GetMCInfo()->Py();
            fPTrackMC[i][2] = fDaughter[i]->GetMCInfo()->Pz();
            fMotherPDG[i] = fDaughter[i]->GetMCInfo()->MotherPDG();
        }
        for (i = 0; i < 3; i++) fPTotMC[i] = fPTrackMC[0][i] + fPTrackMC[1][i];
    }

    // copy reconstructed info (always available)
    for (i = 0; i < 2; i++) {
        fPTrack[i][0] = fDaughter[i]->Px();
        fPTrack[i][1] = fDaughter[i]->Py();
        fPTrack[i][2] = fDaughter[i]->Pz();
    }
    for (i = 0; i < 3; i++) fPTot[i] = fPTrack[0][i] + fPTrack[1][i];
}

//_____________________________________________________________________________
void AliRsnPairParticle::PrintInfo (const Option_t *option)
{
//
// Print some info of the pair. 
// The options are passed to the AliRsnDaughter::Print() method
//

    AliInfo("======== BEGIN PAIR INFO ===========");
    AliInfo("Track #1");
    fDaughter[0]->Print(option);
    AliInfo("Track #2");
    fDaughter[1]->Print(option);
    AliInfo ("========= END PAIR INFO ===========");
}
