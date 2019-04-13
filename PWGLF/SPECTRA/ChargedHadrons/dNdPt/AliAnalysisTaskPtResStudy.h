/// \class AliAnalysisTaskPtResStudy
/// \brief Analsysis task to study pT resolution correction
///
/// to study pt Resolution in MC and Data
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Apr 8, 2019

#ifndef AliAnalysisTaskPtResStudy_H
#define AliAnalysisTaskPtResStudy_H

#include "AliAnalysisTaskMKBase.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskPtResStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskPtResStudy();
                                AliAnalysisTaskPtResStudy(const char *name);
        virtual                 ~AliAnalysisTaskPtResStudy();

        virtual void            AddOutput(); //called at the beginning
        virtual void            AnaTrack();  //called once for every track
        virtual void            AnaEvent();  //called once for every event        
        
        static AliAnalysisTaskPtResStudy* AddTaskPtResStudy(const char* name = "TaskPtResStudy", const char* outfile = 0);

    protected:    
        THnSparseD*             fHistPtResCov;     //-> pt resolution from covariance matrix
        THnSparseD*             fHistPtResMC;      //-> pt resolution from mc
        THnSparseD*             fHistPtRes;        //-> pt resolution, combination of covariance and track fit
        
    private:
        AliAnalysisTaskPtResStudy(const AliAnalysisTaskPtResStudy&); // not implemented
        AliAnalysisTaskPtResStudy& operator=(const AliAnalysisTaskPtResStudy&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskPtResStudy, 1);
    /// \endcond        
};

#endif