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
// -----                              R3BTrack                             -----
// -----                 Created on 12.03.2023 by A.Lagni for S522  -----
// -----------------------------------------------------------------------------

#include "R3BTrackS522.h"

R3BTrackS522::R3BTrackS522()
	: fF1(0.,0.,0.)
	, fF2(0.,0.,0.)
	, fF15(0.,0.,0.)
	, fF16(0.,0.,0.)
	, fTR(0.)
	, fTL(0.)
	, fPR(0.)
	, fPL(0.)
	, fER(0.)
	, fEL(0.)
	, fX(0.)
	, fY(0.)
	, fZ(0.)
	, fPx(0.)
	, fPy(0.)
	, fPz(0.)
	, fVx(0.)
	, fVy(0.)
	, fVz(0.)
	, fOpa(0)
	, fQ(0)
	, fAoZ(0.)
	, fMul(0.)
	, fChix(0.)
	, fToFm(0.)
	, fFlightPathm(0.)
	, fChiy(0.)
	  , fQuality(0.)
{
}

R3BTrackS522::R3BTrackS522(TVector3 F1,
		TVector3 F2,
		TVector3 F15,
		TVector3 F16,
		Double_t tr,
		Double_t tl,
		Double_t pr,
		Double_t pl,
		Double_t er,
		Double_t el,
		Double_t x,
		Double_t y,
		Double_t z,
		Double_t px,
		Double_t py,
		Double_t pz,
		Double_t vx,
		Double_t vy,
		Double_t vz,
		Double_t opa,
		Double_t q,
		Double_t AoZ,
		Int_t Mul,
		Double_t ToFm,
		Double_t FlightPathm,
		Double_t chix,
		Double_t chiy,
		Int_t quality)
	: fF1(F1)
	, fF2(F2)
	, fF15(F15)
	, fF16(F16)
	, fTR(tr)
	, fTL(tl)
	, fPR(pr)
	, fPL(pl)
	, fER(er)
	, fEL(el)
	, fX(x)
	, fY(y)
	, fZ(z)
	, fPx(px)
	, fPy(py)
	, fPz(pz)
	, fVx(vx)
	, fVy(vy)
	, fVz(vz)
	, fOpa(opa)
	, fQ(q)
	, fAoZ(AoZ)
	, fMul(Mul)
	, fToFm(ToFm)
	, fFlightPathm(FlightPathm)
	, fChix(chix)
	, fChiy(chiy)
	  , fQuality(quality)
{
}

R3BTrackS522::~R3BTrackS522() {}

ClassImp(R3BTrackS522)
