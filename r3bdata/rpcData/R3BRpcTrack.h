/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// -----------------------------------------------------------------------------
// -----                              R3BRpcTrack                             -----
// -----                 Created on 09.03.2020 by M.Heil                   -----
// -----------------------------------------------------------------------------

#ifndef R3BRPCTRACK_H
#define R3BRPCTRACK_H

#include "TObject.h"
#include "TVector3.h"

class R3BRpcTrack : public TObject
{
  public:
    R3BRpcTrack();
    R3BRpcTrack(Double_t X,
             Double_t Y,
             Double_t Z,
             Double_t dx,
             Double_t dy,
             Double_t TX,
             Double_t TY,
             TVector3 PoQ,
             Double_t Beta,
             Double_t Gamma,
             Double_t FlighPath,
             Double_t AoZ,
             Double_t mul
             );
    virtual ~R3BRpcTrack();

    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Double_t& GetZ() const { return fZ; }
    inline const Double_t& Getdx() const { return fdx; }
    inline const Double_t& Getdy() const { return fdy; }
    inline const Double_t& GetTX() const { return fTX; }
    inline const Double_t& GetTY() const { return fTY; }
    inline const TVector3& GetPoQ() const { return fPoQ; }
    inline const Double_t& GetBeta() const { return fBeta; }
    inline const Double_t& GetGamma() const { return fGamma; }
    inline const Double_t& GetFlightPath() const { return fFlightPath; }
    inline const Double_t& GetAoZ() const { return fAoZ; }
    inline const Double_t& GetMul() const { return fMul; }

  protected:
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fdx;
    Double_t fdy;
    Double_t fTX;
    Double_t fTY;
    TVector3 fPoQ;
    Double_t fGamma;
    Double_t fBeta;
    Double_t fFlightPath;
    Double_t fAoZ;
    Double_t fMul;

  public:
    ClassDef(R3BRpcTrack, 1)
};

#endif
