/******************************************************************************
 *   Copyright (C) 2024 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2024 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// Created on 09/02/2024 by M.Xarepe

#ifndef R3BTRACKINGG249
#define R3BTRACKINGG249

#include "FairTask.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BNeulandHit.h"
#include "R3BTofdHitData.h"
#include "R3BFootHitData.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"
#include "TCanvas.h"
#include "TLatex.h"
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

class TH2F;
class R3BTrackingG249 : public FairTask
{
	public:
		R3BTrackingG249();
		R3BTrackingG249(const char* name, int iVerbose = 1);
		virtual ~R3BTrackingG249();
		virtual InitStatus Init();
		virtual void Exec(Option_t* option);
		virtual void FinishEvent();
		virtual void FinishTask();

		void Reset_Tracker_Histo();
		// GLAD current which was used for training MDF function
		// Should  be set form the steering macro
		void SetGladReferenceCurrent(double cur) { GladReferenceCurrent = cur; }

		// GLAD current in the run being anlaysed, set from steering macro
		void SetGladCurrent(double cur) { GladCurrent = cur; }
		bool function_type;
		void SetVertexfunction(bool tmp) { function_type = true; }

		// Set MDF functions from the steering macro
		void Set_MDF_PoQ(TString name) { MDF_PoQ_filename = name; }
		void Set_MDF_FlightPath(TString name) { MDF_FlightPath_filename = name; }
		void Set_MDF_TX0(TString name) { MDF_TX0_filename = name; }
		void Set_MDF_TY0(TString name) { MDF_TY0_filename = name; }
		void Set_MDF_TX1(TString name) { MDF_TX1_filename = name; }
		void Set_MDF_TY1(TString name) { MDF_TY1_filename = name; }

		void SetTrigger(int trigger) { fTrigger = trigger; }
		void SetTpat(int tpat) { fTpat = tpat; }
		void SetMaxEvent(int nev) { maxevent = nev; }

		// Set lab positions and angles of the detectors from the steering macro
		// MWPCs
		void SetPositionMwpc0(double x, double y, double z) {  m0_position.SetXYZ(x, y, z); }// cm
		void SetPositionMwpc1(double x, double y, double z) {  m1_position.SetXYZ(x, y, z); }// cm
		void SetEulerAnglesMwpc0(double x, double y, double z) {  m0_angles.SetXYZ(x, y, z); }// rad
		void SetEulerAnglesMwpc1(double x, double y, double z) {  m1_angles.SetXYZ(x, y, z); }// rad

		//Fibers
		void SetPositionFiber30(double x, double y, double z) { f30_position.SetXYZ(x, y, z); }// cm
		void SetPositionFiber31(double x, double y, double z) { f31_position.SetXYZ(x, y, z); }// cm
		void SetPositionFiber32(double x, double y, double z) { f32_position.SetXYZ(x, y, z); }// cm
		void SetPositionFiber33(double x, double y, double z) { f33_position.SetXYZ(x, y, z); }// cm
		void SetEulerAnglesFiber30(double x, double y, double z) { f30_angles.SetXYZ(x, y, z); }// rad
		void SetEulerAnglesFiber31(double x, double y, double z) { f31_angles.SetXYZ(x, y, z); }// rad
		void SetEulerAnglesFiber32(double x, double y, double z) { f32_angles.SetXYZ(x, y, z); }// rad
		void SetEulerAnglesFiber33(double x, double y, double z) { f33_angles.SetXYZ(x, y, z); }// rad

		// Get lab positions and angles (needed by alignment function)
		inline TVector3 GetPositionMwpc0() { return m0_position; }   // cm
		inline TVector3 GetPositionMwpc1() { return m1_position; }   // cm
		inline TVector3 GetPositionFiber30() { return f30_position; }   // cm
		inline TVector3 GetPositionFiber31() { return f31_position; }   // cm
		inline TVector3 GetPositionFiber32() { return f32_position; }   // cm
		inline TVector3 GetPositionFiber33() { return f33_position; }   // cm

		inline TVector3 GetEulerAnglesMwpc0() { return m0_angles; }   // cm
		inline TVector3 GetEulerAnglesMwpc1() { return m1_angles; }   // cm
		inline TVector3 GetEulerAnglesFiber30() { return f30_angles; }   // cm
		inline TVector3 GetEulerAnglesFiber31() { return f31_angles; }   // cm
		inline TVector3 GetEulerAnglesFiber32() { return f32_angles; }   // cm
		inline TVector3 GetEulerAnglesFiber33() { return f33_angles; }   // cm

		R3BMDFWrapper* Get_MDF_PoQ() { return MDF_PoQ; }
		R3BMDFWrapper* Get_MDF_TX0() { return MDF_TX0; }
		R3BMDFWrapper* Get_MDF_TY0() { return MDF_TY0; }

		void SetTofOffset(double offset) { tof_offset = offset; } // ns

		// Transofrming input detector hit (point) into laboratory system
		void TransformPoint(TVector3& point, TVector3* rotation, TVector3* translation);
		void TransformPoint1(TVector3& point1,  TVector3 rotation1, TVector3 translation1);
		// Setup energy cuts in foot and fibers 
		void SetFiberEnergyMinMax(double min, double max){FiberEnergyMin = min; FiberEnergyMax=max;};
		// Setup energy cuts in foot and fibers 
		void SetFiberTimeMinMax(double min, double max){FiberTimeMin = min; FiberTimeMax=max;};
		void SetFootEnergyMinMax(double min, double max){FootEnergyMin = min; FootEnergyMax=max;};

		// Setters for the alignment procedure
		void SetReferencePoQ(double val) { reference_PoQ = val; }
		double GetReferencePoQ() { return reference_PoQ; }

		//Setters in case no FRS
		void SetReference_frs_beta(double val) { FRS_BETA = val; }
		double GetReference_frs_beta() { return FRS_BETA; }

		void Alignment();
		static double AlignmentErrorG249(const double* par);

		//Storing indices of hits in TCA for potential track candidates

		struct Track
		{
			double f1_x;
			double f1_z;
			double f1_Q;
			double f2_y;
			double f2_z;
			double f2_Q;
			double f15_x;
			double f15_z;
			double f15_Q;
			double f16_y;
			double f16_z;
			double f16_Q;
			double tx0;
			double ty0;
			double f32_x;
			double f32_z;
			double f30_x;
			double f30_y;
			double f30_z;
			double last_x;
			double last_z;
			double mw_ax;
			double mw_ay;
			bool fiber;
		};

		std::vector<Track> tracks_in; //track candidates in FOOT
		std::vector<Track> tracks_out;//track candidates in FOOT

		TTree *vt;
		bool IsGoodFiberHit(R3BFiberMAPMTHitData* fhit);
		bool IsGoodFootHit(R3BFootHitData* fhit);
		bool SortFootData();
		bool MakeIncomingTracks();
		bool MakeOutgoingTracks();
		// Data containesr needed only for the Alignment() function
		struct det_points
		{
			TVector3 mw1;
			TVector3 mw0;
			TVector3 f30;
			TVector3 f32;
			TVector3 flast;
		};
		vector<int> f1_hits;
		vector<int> f2_hits;
		vector<int> f15_hits;
		vector<int> f16_hits;

		std::vector<det_points> det_points_vec;
		TVector3 m0_point, m1_point, f1_point, f2_point, f15_point,f16_point, f30_point, f32_point, flast_point;
		TVector3 m0_point_i, m1_point_i, f1_point_i, f2_point_i, f15_point_i,f16_point_i, f30_point_i, f32_point_i, flast_point_i;
		TVector3 m0_point_i1, m1_point_i1, f1_point_i1, f2_point_i1, f15_point_i1,f16_point_i1, f30_point_i1, f32_point_i1, flast_point_i1;

	private:
		// Input hit data from the TClonesArray
		// do not change the order, add new det in the end
		R3BCoarseTimeStitch* fTimeStitch;
		enum DetectorInstances
		{
			DET_FI_FIRST,
			DET_FI30 = DET_FI_FIRST,
			DET_FI31,
			DET_FI32,
			DET_FI33,
			DET_FI_LAST = DET_FI33,
			DET_TOFD,
			FOOT_HITDATA,
			MWPC0_HITDATA,
			MWPC1_HITDATA,
			FRS_DATA,
			LOS_DATA,
			/*NEULAND_DATA,*/
			DET_MAX
		};

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

		// Names of essential branches in the input tree
		// do not change the order! add new data in the end
		const char* fDetectorNames[DET_MAX + 1] = { "Fi30Hit", "Fi31Hit", "Fi32Hit", "Fi33Hit",
			"TofdHit", "FootHitData", "Mwpc0HitData", "Mwpc1HitData", "FrsSciTcalData", "LosCal",/* "NeulandHits",*/ NULL };

		R3BEventHeader* fHeader;
		std::vector<TClonesArray*> fDataItems; // input data
		TClonesArray* fTrackItems;             // output data

		bool is_good_event;

		TRotation r;
		TVector3 v3_localX;
		TVector3 v3_localZ;

		TVector3 m0_position;
		TVector3 m1_position;
		TVector3 f30_position;
		TVector3 f31_position;
		TVector3 f32_position;
		TVector3 f33_position;
		TVector3 f1_position;
		TVector3 f2_position;
		TVector3 f15_position;
		TVector3 f16_position;

		TVector3 m0_angles;
		TVector3 m1_angles;
		TVector3 f30_angles;
		TVector3 f31_angles;
		TVector3 f32_angles;
		TVector3 f33_angles;
		TVector3 f1_angles;
		TVector3 f2_angles;
		TVector3 f15_angles;
		TVector3 f16_angles;

		R3BMDFWrapper* MDF_FlightPath;
		R3BMDFWrapper* MDF_PoQ;
		R3BMDFWrapper* MDF_TX0;
		R3BMDFWrapper* MDF_TX1;
		R3BMDFWrapper* MDF_TY0;
		R3BMDFWrapper* MDF_TY1;

		TString MDF_FlightPath_filename;
		TString MDF_PoQ_filename;
		TString MDF_TX0_filename;
		TString MDF_TY0_filename;
		TString MDF_TX1_filename;
		TString MDF_TY1_filename;

		double mdf_data[8];   // data container for the MDF function
		unsigned long fNEvents=0; // Event counter
		int fTrigger;
		int fTpat;
		int maxevent;
		double GladCurrent;
		double GladReferenceCurrent;
		double reference_PoQ ;
		double FRS_BETA = 0.714549;// 12C
		//double FRS_BETA = 0.722; // 10C
		//double FRS_BETA = 0.699804; // 10C
		//double FRS_BETA = 0.714; // 16C
		//double FRS_BETA = 0.714034; // 12C 2nd
		//double FRS_BETA = 0.706721; // 14C
		Bool_t DoAlignment;
		//double tof_offset = 27.431851; // ns 12C
		//double tof_offset = 26.997; // ns 16C
		double tof_offset = 27.285005; // ns 14C
		//double tof_offset = 27.010147; // ns 12C 2nd

		// Cut on fiber hit energy set by SetFiberEnergyMinMax():
		double FiberEnergyMin = 0;
		double FiberEnergyMax = 0;
		double FiberTimeMin = 0;
		double FiberTimeMax = 0;
		double FootEnergyMin = 0;
		double FootEnergyMax = 0;

		static constexpr int N_glob_tracks_max = 10000000;
		int mul_los=-999;
		int mul_m0=-999;
		int mul_m1=-999;
		int mul_f1=-999;
		int mul_f2=-999;
		int mul_f15=-999;
		int mul_f16=-999;
		int mul_f30=-999;
		int mul_f31=-999;
		int mul_f32=-999;
		int mul_f33=-999;
		int mul_tofd=-999;
		int mul_foot=-999;
		int mul_frsi=-999;
		int mul_neuland=-999;
		int Tpat = -999;
		bool cond=false;
		int add_track_counter = 0;
		double a_counter=0;
		double b_counter=0;
		double c_counter=0;
		double d_counter=0;
		double e_counter=0;
		double f_counter=0;
		double g_counter=0;
		double h_counter=0;
		double hh_counter=0;
		double hhh_counter=0;

		double tx0 = -999;
		double ty0 = -999;
		double tx1 = -999;
		double ty1 = -999;
		double beta = -999;
		double gamma = -999;
		double beta_tof = -999;
		double gamma_tof = -999;
		double poq = -999;
		double tofdq = -999;
		double flight_p = -999;
		double tof = -999;
		double maoz = -999;
		double maoz_tof = -999;

		// Essential constants
		const double SPEED_OF_LIGHT = 29.9792458; // cm/ns
		const double AMU = 0.9314940038;          // GeV/c2

		/* ----- tracker Canvases ----- */
		TCanvas* trackerCanvas;
		TCanvas* frs_betaCanvas;
		TCanvas* AoQ_vs_pos_both_fibCanvas;
		TCanvas* Q_vs_pos_Canvas;
		TCanvas* trackerCanvas_tpat;
		TCanvas* trackerAnglesCanvas;
		TCanvas* neutronAnglesCanvas;
		TCanvas* betaCanvas;
		TCanvas* e_relCanvas;
		TCanvas* trackerCanvas_Vs_tpat;

		/* ----- tracker Histograms ----- */
		TH2F* AoQ_vs_fib_pos_q_6_aoz_2;

		TH2F* AoQ_vs_TOFD_pos_q_1;
		TH2F* AoQ_vs_TOFD_pos_q_2;
		TH2F* AoQ_vs_TOFD_pos_q_3;
		TH2F* AoQ_vs_TOFD_pos_q_4;
		TH2F* AoQ_vs_TOFD_pos_q_5;
		TH2F* AoQ_vs_TOFD_pos_q_6;

		TH2F* Q_TOFD_vs_TOFD_pos;
		TH2F* AoQ_Vs_Q_TOFD;
		TH2F* AoQ_vs_TOFD_pos;
		TH2F* AoQ_Vs_Q_TOFD_tpat;
		TH2F* AoQ_tof_Vs_Q_TOFD_tpat;
		TH2F* AoQ_tof_Vs_Q_TOFD;
		TH2F* AoQ_Vs_tpat;
		TH2F* TX0_Vs_fib_pos_Q_5_AoZ_2;
		TH2F* TY0_Vs_fib_pos_Q_5_AoZ_2;

		TH1F* TX_pos_Q_5_AoZ_2;
		TH1F* TY_pos_Q_5_AoZ_2;
		TH1F* neutron_beta;
		TH1F* frs_beta;
		TH1F* neutron_e_rel;
		// Private method to fill output track data
		R3BTrack* AddTrackData(double x, double y, double z, TVector3 poq_vec, Double_t charge, Double_t aoz);

	public:
		ClassDef(R3BTrackingG249, 1)
};

#endif
