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

// Created on 23/22/2022 by V.Panin

#ifndef R3BTRACKINGSPROTONS
#define R3BTRACKINGSPROTONS

#include "FairTask.h"
#include "R3BEventHeader.h"
#include "R3BFootHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMDFWrapper.h"
#include "R3BRpcTrack.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TRotation.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TTree.h"
class R3BTrackingProtons : public FairTask
{
	public:
		R3BTrackingProtons();
		R3BTrackingProtons(const char* name, Int_t iVerbose = 1);
		virtual ~R3BTrackingProtons();
		virtual InitStatus Init();
		virtual void Exec(Option_t* option);
		virtual void FinishEvent();
		virtual void FinishTask();

		// Set MDF functions from the steering macro
		void Set_MDF_PoQ(TString name)             { MDF_PoQ_filename = name; }
		void Set_MDF_FlightPath(TString name)      { MDF_FlightPath_filename = name; }
		void Set_MDF_TX0(TString name)             { MDF_TX0_filename = name; }
		void Set_MDF_TY0(TString name)             { MDF_TY0_filename = name; }
		void SetTrigger(Int_t trigger)             { fTrigger = trigger; }
		void SetTpat(Int_t tpat)                   { fTpat = tpat; }
		void SetMaxEvent(Int_t nev)                { maxevent = nev; }

		// Set lab positions and angles of the detectors from the steering macro
		// MWPCs
		void SetPositionMwpc0(double x, double y, double z)    { m0_position.SetXYZ(x, y, z);  }// cm
		void SetPositionMwpc1(double x, double y, double z)    { m1_position.SetXYZ(x, y, z);  }// cm
		void SetEulerAnglesMwpc0(double x, double y, double z) { m0_angles.SetXYZ(x, y, z);    }// rad
		void SetEulerAnglesMwpc1(double x, double y, double z) { m1_angles.SetXYZ(x, y, z);    }// rad

		// FOOTs
		void SetPositionFoot1(double x, double y, double z)    { f1_position.SetXYZ(x, y, z);  }// cm
		void SetPositionFoot2(double x, double y, double z)    { f2_position.SetXYZ(x, y, z);  }// cm
		void SetEulerAnglesFoot1(double x, double y, double z) { f1_angles.SetXYZ(x, y, z);    }// rad
		void SetEulerAnglesFoot2(double x, double y, double z) { f2_angles.SetXYZ(x, y, z);    }// rad

		//RPC
		void SetPositionRpc(double x, double y, double z)     { rpc_position.SetXYZ(x, y, z);}// cm
		void SetEulerAnglesRpc(double x, double y, double z)  { rpc_angles.SetXYZ(x, y, z);  }// rad

		// Get lab positions and angles (needed by alignment function)
		inline TVector3 GetPositionMwpc0()     { return m0_position;   } // cm
		inline TVector3 GetPositionMwpc1()     { return m1_position;   } // cm
		inline TVector3 GetPositionFoot1()     { return f1_position;   } // cm
		inline TVector3 GetPositionFoot2()     { return f2_position;   } // cm
		inline TVector3 GetPositionRpc()       { return rpc_position;  } // cm

		inline TVector3 GetEulerAnglesMwpc0()  { return m0_angles;     } // cm
		inline TVector3 GetEulerAnglesMwpc1()  { return m1_angles;     } // cm
		inline TVector3 GetEulerAnglesFoot1()  { return f1_angles;     } // cm
		inline TVector3 GetEulerAnglesFoot2()  { return f2_angles;     } // cm
		inline TVector3 GetEulerAnglesRpc()    { return rpc_angles;    } // rad

		R3BMDFWrapper* Get_MDF_PoQ() { return MDF_PoQ; }

		void SetTofOffset(double offset) { tof_offset = offset; } // ns

		// Transofrming input detector hit (point) into laboratory system
		void TransformPoint  (TVector3& point, TVector3* rotation, TVector3* translation);
		// Setup energy cuts in foot and fibers 
		void SetFootEnergyMinMax(double min, double max);

		//Storing indices of hits in TCA for potential track candidates
		struct Track{
			double f1_x={};
			double f1_z={};
			double f2_y={};
			double f2_z={};
			double f1_Q={};
			double f2_Q={};
			double vertex_mwpc_X={};
			double vertex_mwpc_Y={};
			double rpc_x={};
			double rpc_y={};
			double rpc_z={};
			double rpc_tof={};
		};
		std::vector<Track> tracks_in; //track candidates in FOOT
		std::vector<Track> tracks_out;//track candidates in FOOT

		//Storage of hit indices
		vector<UInt_t> f1_hits={};
		vector<UInt_t> f2_hits={};
		TVector3 m0_point, m1_point, f1_point, f2_point, rpc_point;
		TVector3 m0_point_i, m1_point_i, f1_point_i, f2_point_i, rpc_point_i;
		double rpc_tof={0},rpc_tof_i={0};
		bool IsGoodFootHit(R3BFootHitData* fhit);
		bool SortFootData();
		bool MakeIncomingTracks();
		bool MakeOutgoingTracks();

	private:
		// Input hit data from the TClonesArray
		// do not change the order, add new det in the end
		enum DetectorInstances
		{
			FOOT_HITDATA,
			MWPC0_HITDATA,
			MWPC1_HITDATA,
			FRS_DATA,
			DET_RPC,
			DET_LOS,
			MUSLI,
			DET_MAX
		};

		// Names of essential branches in the input tree
		// do not change the order! add new data in the end
		const char* fDetectorNames[DET_MAX + 1] = {"FootHitData", "Mwpc0HitData", "Mwpc1HitData", "FrsData", "R3BRpcHitData", "LosHit", "MusliHitData", NULL };

		R3BEventHeader* fHeader;
		std::vector<TClonesArray*> fDataItems; // input data
		TClonesArray* fTrackItems;             // output data

		bool is_good_event;

		TRotation r;
		TVector3 v3_localX;
		TVector3 v3_localZ;

		TVector3 m0_position;
		TVector3 m1_position;
		TVector3 f1_position;
		TVector3 f2_position;

		TVector3 m0_angles;
		TVector3 m1_angles;
		TVector3 f1_angles;
		TVector3 f2_angles;

		TVector3 rpc_position;
		TVector3 rpc_angles;

		R3BMDFWrapper* MDF_FlightPath;
		R3BMDFWrapper* MDF_PoQ;
		R3BMDFWrapper* MDF_TX0;
		R3BMDFWrapper* MDF_TY0;

		TString MDF_FlightPath_filename;
		TString MDF_PoQ_filename;
		TString MDF_TX0_filename;
		TString MDF_TY0_filename;

		Double_t mdf_data[8];   // data container for the MDF function
		unsigned long fNEvents; // Event counter
		Int_t fTrigger;
		Int_t fTpat;
		Int_t maxevent;
		Double_t tof_offset; // ns

		// Energy range in FOOT 
		Double_t FootEnergyMin;
		Double_t FootEnergyMax;


		TTree tree_out;
		UInt_t N_glob_tracks={};
		UInt_t N_in_tracks={};
		UInt_t N_out_tracks={};
		static constexpr UInt_t N_glob_tracks_max = 10000000;
		UInt_t good_ev={};
		UInt_t mul_los={};
		UInt_t mul_m0={};
		UInt_t mul_m1={};
		UInt_t mul_f1={};
		UInt_t mul_f2={};
		UInt_t mul_rpc={};
		UInt_t mul_foot={};
		vector <double> vtx;
                vector <double> vty;
                vector <double> vtz;
                vector <double> distm;
                vector <double> opan;

		Int_t Tpat={};
		Double_t PoQ[N_glob_tracks_max]={};
		Double_t FlightPath[N_glob_tracks_max]={};
		Double_t TX0[N_glob_tracks_max]={};
		Double_t TY0[N_glob_tracks_max]={};
		Double_t Beta[N_glob_tracks_max]={};
		Double_t Gamma[N_glob_tracks_max]={};
		Double_t mdf_AoZ[N_glob_tracks_max]={};

		Float_t  m0_X={};
		Float_t  m0_Y={};
		Float_t  m0_Z={};

		Float_t  m1_X={};
		Float_t  m1_Y={};
		Float_t  m1_Z={};

		Float_t vertex_foot_X[N_glob_tracks_max]={};
		Float_t vertex_foot_Y[N_glob_tracks_max]={};

		Float_t vertex_mwpc_X[N_glob_tracks_max]={};
		Float_t vertex_mwpc_Y[N_glob_tracks_max]={};

		Float_t dx_vertex[N_glob_tracks_max]={};
		Float_t dy_vertex[N_glob_tracks_max]={};

		Float_t f1_X[N_glob_tracks_max]={};
		Float_t f1_Y[N_glob_tracks_max]={};
		Float_t f1_Z[N_glob_tracks_max]={};
		Float_t f1_Q[N_glob_tracks_max]={};
		Float_t f1_T[N_glob_tracks_max]={};

		Float_t f2_X[N_glob_tracks_max]={};
		Float_t f2_Y[N_glob_tracks_max]={};
		Float_t f2_Z[N_glob_tracks_max]={};
		Float_t f2_Q[N_glob_tracks_max]={};
		Float_t f2_T[N_glob_tracks_max]={};

		Float_t rpc_X[N_glob_tracks_max]={};
		Float_t rpc_Y[N_glob_tracks_max]={};
		Float_t rpc_Z[N_glob_tracks_max]={};
		Float_t rpc_T[N_glob_tracks_max]={};

		// Essential constants
		const Double_t SPEED_OF_LIGHT = 29.9792458; // cm/ns
		const Double_t AMU = 0.9314940038;          // GeV/c2

		// Private method to fill output track data
		R3BRpcTrack* AddTrackData(TVector3 mw, Double_t dx, Double_t dy, Double_t TX, Double_t TY, TVector3 poq, Double_t beta, Double_t gamma, Double_t flightPath, Double_t AoZ, Double_t mul);

	public:
		ClassDef(R3BTrackingProtons, 1)
};

#endif
