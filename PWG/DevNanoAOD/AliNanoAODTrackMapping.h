/// \class AliNanoAODTrackMapping

#ifndef _ALINANOAODTRACKMAPPING_H_
#define _ALINANOAODTRACKMAPPING_H_

#include "TObject.h"
#include "TFile.h"
#include "AliLog.h"
#include "TSystem.h"
#include "TTree.h"
#include "TDirectory.h"

class AliNanoAODTrackMapping : public TObject
{
public:
  AliNanoAODTrackMapping();
  AliNanoAODTrackMapping(const char * mappingString);
  virtual ~AliNanoAODTrackMapping(){;}

  void Print(const Option_t * opt = "") const;

  static AliNanoAODTrackMapping * GetInstance(const char * vars = 0) {
    if(!fInstance) {
      if (vars) fInstance = new AliNanoAODTrackMapping(vars);
      else LoadInstance() ;
    }
    
    return fInstance;
  }  
  

  const char * GetVarName(Int_t index) const;
  Int_t GetVarIndex(TString varName); // cannot be const (uses stl map)

  //TODO: implement custom variables

  // Getters
  //  Internal
  Int_t GetSize()             const { return fSize;             }  
  //  Kin vars
  Int_t GetPt()               const { return fPt;               }
  Int_t GetPhi()              const { return fPhi;              }
  Int_t GetTheta()            const { return fTheta;            }
  Int_t GetChi2PerNDF()       const { return fChi2PerNDF;       }
  Int_t GetPosX()             const { return fPosX;             }
  Int_t GetPosY()             const { return fPosY;             }
  Int_t GetPosZ()             const { return fPosZ;             }
  Int_t GetPDCAX()            const { return fPDCAX;            }
  Int_t GetPDCAY()            const { return fPDCAY;            }
  Int_t GetPDCAZ()            const { return fPDCAZ;            }
  Int_t GetPosDCAx()          const { return fPosDCAx;          }
  Int_t GetPosDCAy()          const { return fPosDCAy;          }
  Int_t GetRAtAbsorberEnd()   const { return fRAtAbsorberEnd;   }
  Int_t GetTPCncls()          const { return fTPCncls;          }
  Int_t GetID()               const { return fID;               }
  Int_t GetTPCnclsF()         const { return fTPCnclsF;         }
  Int_t GetTPCNCrossedRows()  const { return fTPCNCrossedRows;  }
  Int_t GetTrackPhiOnEMCal()  const { return fTrackPhiOnEMCal;  }
  Int_t GetTrackEtaOnEMCal()  const { return fTrackEtaOnEMCal;  }
  Int_t GetTrackPtOnEMCal()   const { return fTrackPtOnEMCal;   }
  Int_t GetITSsignal()        const { return fITSsignal;        }
  Int_t GetTPCsignal()        const { return fTPCsignal;        }
  Int_t GetTPCsignalTuned()   const { return fTPCsignalTuned;   }
  Int_t GetTPCsignalN()       const { return fTPCsignalN;       }
  Int_t GetTPCmomentum()      const { return fTPCmomentum;      }
  Int_t GetTPCTgl()           const { return fTPCTgl;           }
  Int_t GetTOFsignal()        const { return fTOFsignal;        }
  Int_t GetintegratedLenght() const { return fintegratedLenght; }
  Int_t GetTOFsignalTuned()   const { return fTOFsignalTuned;   }
  Int_t GetHMPIDsignal()      const { return fHMPIDsignal;      }
  Int_t GetHMPIDoccupancy()   const { return fHMPIDoccupancy;   }
  Int_t GetTRDsignal()        const { return fTRDsignal;        }
  Int_t GetTRDChi2()          const { return fTRDChi2;          }
  Int_t GetTRDnSlices()       const { return fTRDnSlices;       }
  Int_t GetIsMuonTrack()      const { return fIsMuonTrack;      }
  Int_t GetTPCnclsS()         const { return fTPCnclsS;         }
  Int_t GetFilterMap()        const { return fFilterMap;        }
  Int_t GetTOFBunchCrossing() const { return fTOFBunchCrossing; }
  Int_t GetCovMat(Int_t index)  const { return fcovmat[index];  }

  


private:

  static void  LoadInstance() ;
  
  Int_t fSize; ///< Number of variables actually allocated
  void  SetSize (Int_t var) { fSize = var;}
  // FIXME: should this be static?
  Int_t fPt;      	  ///< Mapping variable
  Int_t fPhi;		  ///< Mapping variable
  Int_t fTheta;		  ///< Mapping variable
  Int_t fChi2PerNDF;	  ///< Mapping variable
  Int_t fPosX;		  ///< Mapping variable
  Int_t fPosY;		  ///< Mapping variable
  Int_t fPosZ;		  ///< Mapping variable
  Int_t fPDCAX;		  ///< Mapping variable
  Int_t fPDCAY;		  ///< Mapping variable
  Int_t fPDCAZ;		  ///< Mapping variable
  Int_t fPosDCAx;	  ///< Mapping variable
  Int_t fPosDCAy;	  ///< Mapping variable
  Int_t fRAtAbsorberEnd;  ///< Mapping variable
  Int_t fTPCncls;	  ///< Mapping variable
  Int_t fID;		  ///< Mapping variable
  Int_t fTPCnclsF;	  ///< Mapping variable
  Int_t fTPCNCrossedRows; ///< Mapping variable
  Int_t fTrackPhiOnEMCal; ///< Mapping variable
  Int_t fTrackEtaOnEMCal; ///< Mapping variable
  Int_t fTrackPtOnEMCal; ///< Mapping variable
  Int_t fITSsignal;	  ///< Mapping variable
  Int_t fTPCsignal;	  ///< Mapping variable
  Int_t fTPCsignalTuned;  ///< Mapping variable
  Int_t fTPCsignalN;	  ///< Mapping variable
  Int_t fTPCmomentum;	  ///< Mapping variable
  Int_t fTPCTgl;	  ///< Mapping variable
  Int_t fTOFsignal;	  ///< Mapping variable
  Int_t fintegratedLenght;///< Mapping variable
  Int_t fTOFsignalTuned;  ///< Mapping variable
  Int_t fHMPIDsignal;	  ///< Mapping variable
  Int_t fHMPIDoccupancy;  ///< Mapping variable
  Int_t fTRDsignal;	  ///< Mapping variable
  Int_t fTRDChi2;	  ///< Mapping variable
  Int_t fTRDnSlices;	  ///< Mapping variable
  Int_t fIsMuonTrack;     ///< Mapping variable
  Int_t fTPCnclsS;    ///< Mapping variable
  Int_t fFilterMap;      ///< Mapping variable
  Int_t fTOFBunchCrossing;      ///< Mapping variable
  Int_t fcovmat[21];          ///< Mapping variables

  static AliNanoAODTrackMapping * fInstance; ///< instance, needed for the singleton implementation
  static TString fMappingString; ///< the string which this class was initialized with
  std::map<TString,int> fMapCstVar;// Map of indexes of custom variables: CACHE THIS TO CONST INTs IN YOUR TASK TO AVOID CONTINUOUS STRING COMPARISONS
  ClassDef(AliNanoAODTrackMapping, 2)
  
};

#endif /* _ALINANOAODTRACKMAPPING_H_ */
