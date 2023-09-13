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
// -----                              R3BTrackS522                             -----
// -----                 Created on 22.03.2023 by A.Lagni for S522     -----
// -----------------------------------------------------------------------------

#ifndef R3BTrackS522_H
#define R3BTrackS522_H

#include "TObject.h"

class R3BTrackS522 : public TObject
{
	public:
		R3BTrackS522();
		R3BTrackS522(Double_t tr,
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
				Double_t Flight_Pathm,
				Double_t chix,
				Double_t chiy,
				Int_t quality);
		virtual ~R3BTrackS522();

		inline const Double_t& GetTr() const { return fTR; }
		inline const Double_t& GetTl() const { return fTL; }
		inline const Double_t& GetPr() const { return fPR; }
		inline const Double_t& GetPl() const { return fPL; }
		inline const Double_t& GetEr() const { return fER; }
		inline const Double_t& GetEl() const { return fEL; }
		inline const Double_t& GetX() const { return fX; }
		inline const Double_t& GetY() const { return fY; }
		inline const Double_t& GetZ() const { return fZ; }
		inline const Double_t& GetPx() const { return fPx; }
		inline const Double_t& GetPy() const { return fPy; }
		inline const Double_t& GetPz() const { return fPz; }
		inline const Double_t& GetVx() const { return fVx; }
		inline const Double_t& GetVy() const { return fVy; }
		inline const Double_t& GetVz() const { return fVz; }
		inline const Double_t& GetOpa() const { return fOpa; }
		inline const Double_t& GetToF() const { return fToFm; }
		inline const Double_t& GetFlightPath() const { return fFlightPathm; }
		inline const Double_t& GetQ() const { return fQ; }
		inline const Double_t& GetAoZ() const { return fAoZ; }
		inline const Int_t& GetMul() const { return fMul; }
		inline const Double_t& GetChix() const { return fChix; }
		inline const Double_t& GetChiy() const { return fChiy; }
		inline const Int_t& GetQuality() const { return fQuality; }

	protected:
		Double_t fTR;
		Double_t fTL;
		Double_t fPR;
		Double_t fPL;
		Double_t fER;
		Double_t fEL;
		Double_t fX;
		Double_t fY;
		Double_t fZ;
		Double_t fPx;
		Double_t fPy;
		Double_t fPz;
		Double_t fVx;
		Double_t fVy;
		Double_t fVz;
		Double_t fOpa;
		Double_t fQ;
		Double_t fAoZ;
		Int_t fMul;
		Double_t fToFm;
		Double_t fFlightPathm;
		Double_t fChix;
		Double_t fChiy;
		Int_t fQuality;

	public:
		ClassDef(R3BTrackS522, 1)
};

#endif
