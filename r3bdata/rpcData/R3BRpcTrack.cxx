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
// -----                 Created on 27.01.2022 by M.Xarepe                   -----
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
	, fPoQ(0.,0.,0.)
	, fBeta(0)
	, fGamma(0.)
	, fFlightPath(0.)
	, fAoZ(0.)
	  ,fMul(0.)
{
}

R3BRpcTrack::R3BRpcTrack(Double_t x,
		Double_t y,
		Double_t z,
		Double_t dx,
		Double_t dy,
		Double_t TX,
		Double_t TY,
		TVector3 PoQ,
		Double_t Beta,
		Double_t Gamma,
		Double_t FlightPath,
		Double_t AoZ,
		Double_t Mul)
	: fX(x)
	, fY(y)
	, fZ(z)
	, fdx(dx)
	, fdy(dy)
	, fTX(TX)
	, fTY(TY)
	, fPoQ(PoQ)
	, fBeta(Beta)
	, fGamma(Gamma)
	, fFlightPath(FlightPath)
	, fAoZ(AoZ)
	  , fMul(Mul)
{
}

R3BRpcTrack::~R3BRpcTrack() {}

ClassImp(R3BRpcTrack)