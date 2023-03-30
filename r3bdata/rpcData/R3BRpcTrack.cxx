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

#include "R3BRpcTrack.h"

R3BRpcTrack::R3BRpcTrack()
    : fX(0.)
    , fY(0.)
    , fZ(0.)
    , fdx(0.)
    , fdy(0.)
    , fTX(0.)
    , fTY(0.)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fBeta(0)
    , fGamma(0.)
    , fFlightPath(0.)
{
}

R3BRpcTrack::R3BRpcTrack(Double_t x,
                   Double_t y,
                   Double_t z,
                   Double_t dx,
                   Double_t dy,
                   Double_t TX,
                   Double_t TY,
                   Double_t px,
                   Double_t py,
                   Double_t pz,
                   Double_t Beta,
                   Double_t Gamma,
                   Double_t FlightPath)
    : fX(x)
    , fY(y)
    , fZ(z)
    , fdx(dx)
    , fdy(dy)
    , fTX(TX)
    , fTY(TY)
    , fPx(px)
    , fPy(py)
    , fPz(pz)
    , fBeta(Beta)
    , fGamma(Gamma)
    , fFlightPath(FlightPath)
{
}

R3BRpcTrack::~R3BRpcTrack() {}

ClassImp(R3BRpcTrack)
