/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BRPCHITPAR_H
#define R3BRPCHITPAR_H

#include <TObjString.h>
#include <TVector3.h>

#include "FairParGenericSet.h"
#include "FairParamList.h"

class R3BRpcHitPar : public FairParGenericSet
{
  public:
    R3BRpcHitPar(const char* name = "R3BRpcHitPar",
                 const char* title = "RPC Hit Finder Parameters",
                 const char* context = "TestDefaultContext");
    ~R3BRpcHitPar(void){};
    void clear(void){};
    void putParams(FairParamList* list);
    Bool_t getParams(FairParamList* list);

    void Print(Option_t* option = "") const;
    /** Accessor functions **/
    const Double_t GetExample() { return fExample; }

    void SetExample(Double_t value) { fExample = value; }

  private:
    // Whatever
    Double_t fExample;

    ClassDef(R3BRpcHitPar, 1); //
};

#endif /* !R3BRPCHITPAR_H*/