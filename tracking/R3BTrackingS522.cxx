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

#include "R3BTrackingS522.h"
#include "R3BFootHitData.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BRpcHitData.h"
#include "R3BFrsData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BRpcTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <TRandomGen.h>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
R3BTrackingS522* gMDFTrackerS522;

double counter_align = 0;

R3BTrackingS522::R3BTrackingS522()
	: R3BTrackingS522("TrackingS522", 1)
{
}

R3BTrackingS522::R3BTrackingS522(const char* name, Int_t iVerbose)
	: FairTask(name, iVerbose)
	, fTrigger(-1)
	, fTpat(-1)
	, fNEvents(0)
	, maxevent(0)
	  , DoAlignment(false)
	, fTrackItems(new TClonesArray("R3BRpcTrack"))
	, reference_PoQ(0.)
	, GladCurrent(-1)
	, GladReferenceCurrent(-1)
	, FootEnergyMin(-1)
	, FootEnergyMax(-1)
	  , fHeader(nullptr)
{
}

R3BTrackingS522::~R3BTrackingS522()
{
	if (fTrackItems)
		delete fTrackItems;
}

InitStatus R3BTrackingS522::Init()
{
	LOG(info) << "R3BTrackingS522::Init()";
	FairRootManager* mgr = FairRootManager::Instance();
	if (NULL == mgr)
	{
		LOG(fatal) << "FairRootManager not found";
	}
	fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
	if (!fHeader)
	{
		LOG(warn) << "R3BTrackingS522::Init() EventHeader. not found";
	}
	// Reading all detector branches
	cout << "\nDET_MAX = " << DET_MAX << endl;
	assert(DET_MAX + 1 >= sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
	for (int det = 0; det < DET_MAX; det++)
	{
		fDataItems.push_back((TClonesArray*)mgr->GetObject(Form("%s", fDetectorNames[det])));
		if (NULL == fDataItems.at(det))
		{
			R3BLOG(fatal, Form("\n\n Cannot find tree branch %s \n\n", fDetectorNames[det]));
		}
		//mgr->Register(Form("%s", fDetectorNames[det]), std::to_string(det).c_str(), fDataItems[det], true);
	}

	// check if all cuts are properly set
	if (GladCurrent < 0 || GladReferenceCurrent < 0 || 
			FootEnergyMin < 0 || FootEnergyMax < 0) 
	{
		R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
	}
	// Initializing all MDF functions
	LOG(info) << "Reading MDF function for FlightPath";
	MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

	LOG(info) << "Reading MDF function for PoQ";
	MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());

	LOG(info) << "Reading MDF function for TX0";
	MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

	LOG(info) << "Reading MDF function for TY0";
	MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());
	// linking to global pointer (needed by alignment)
	gMDFTrackerS522 = this;

	tree_out.SetName("tree_out");
	tree_out.Branch("N_glob_tracks", &N_glob_tracks , "N_glob_tracks/i");
	tree_out.Branch("N_in_tracks"  , &N_in_tracks   , "N_in_tracks/i");
	tree_out.Branch("N_out_tracks" , &N_out_tracks  , "N_out_tracks/i");

	tree_out.Branch("good_ev"      , &good_ev       , "good_ev/i");

	tree_out.Branch("mul_m0"       , &mul_m0        , "mul_m0/i");
	tree_out.Branch("mul_m1"       , &mul_m1        ,  "mul_m1/i");
	tree_out.Branch("mul_foot"     , &mul_foot      , "mul_foot/i");
	tree_out.Branch("mul_f1"       , &mul_f1        , "mul_f1/i");
	tree_out.Branch("mul_f2"       , &mul_f2        , "mul_f2/i");
	tree_out.Branch("mul_f15"       , &mul_f15        , "mul_f15/i");
	tree_out.Branch("mul_f16"       , &mul_f16        , "mul_f16/i");
	tree_out.Branch("mul_tofd"     , &mul_tofd      , "mul_tofd/i");
	tree_out.Branch("mul_los"      , &mul_los       , "mul_los/i");

	tree_out.Branch("Tpat"         , &Tpat          , "Tpat/I");
	tree_out.Branch("PoQ"          , PoQ            , "PoQ[N_glob_tracks]/D");
	tree_out.Branch("FlightPath"   , FlightPath     , "FlightPath[N_glob_tracks]/D");
	tree_out.Branch("TX0"          , TX0            , "TX0[N_glob_tracks]/D");
	tree_out.Branch("TY0"          , TY0            , "TY0[N_glob_tracks]/D");
	tree_out.Branch("Beta"         , Beta           , "Beta[N_glob_tracks]/D");
	tree_out.Branch("Gamma"        , Gamma          , "Gamma[N_glob_tracks]/D");
	tree_out.Branch("mdf_AoZ"      , mdf_AoZ        , "mdf_AoZ[N_glob_tracks]/D");

	tree_out.Branch("m0_X"         , &m0_X          , "m0_X/F");
	tree_out.Branch("m0_Y"         , &m0_Y          , "m0_Y/F");
	tree_out.Branch("m0_Z"         , &m0_Z          , "m0_Z/F");

	tree_out.Branch("m1_X"         , &m1_X          , "m1_X/F");
	tree_out.Branch("m1_Y"         , &m1_Y          , "m1_Y/F");
	tree_out.Branch("m1_Z"         , &m1_Z          , "m1_Z/F");

	tree_out.Branch("frs_Brho"         , &frs_Brho          , "frs_Brho/F");
	tree_out.Branch("frs_A"         , &frs_A          , "frs_A/F");

	tree_out.Branch("vertex_mwpc_X", vertex_mwpc_X  , "vertex_mwpc_X[N_in_tracks]/F");
	tree_out.Branch("vertex_mwpc_Y", vertex_mwpc_Y  , "vertex_mwpc_Y[N_in_tracks]/F");

	tree_out.Branch("vertex_foot_X", vertex_foot_X  , "vertex_foot_X[N_in_tracks]/F");
	tree_out.Branch("vertex_foot_Y", vertex_foot_Y  , "vertex_foot_Y[N_in_tracks]/F");
	tree_out.Branch("dx_vertex"    , dx_vertex      , "dx_vertex[N_in_tracks]/F");
	tree_out.Branch("dy_vertex"    , dy_vertex      , "dy_vertex[N_in_tracks]/F");

	tree_out.Branch("f1_Q"          , f1_Q            , "f1_Q[N_glob_tracks]/F");
	tree_out.Branch("f2_Q"          , f2_Q            , "f2_Q[N_glob_tracks]/F");
	tree_out.Branch("f15_Q"         , f15_Q           , "f15_Q[N_glob_tracks]/F");
	tree_out.Branch("f16_Q"         , f16_Q           , "f16_Q[N_glob_tracks]/F");

	tree_out.Branch("f1_X_i"         , f1_X_i            , "f1_X_i[N_glob_tracks]/F");
	tree_out.Branch("f1_Y_i"         , f1_Y_i            , "f1_Y_i[N_glob_tracks]/F");
	tree_out.Branch("f1_Z_i"         , f1_Z_i            , "f1_Z_i[N_glob_tracks]/F");

	tree_out.Branch("f2_X_i"         , f2_X_i            , "f2_X_i[N_glob_tracks]/F");
	tree_out.Branch("f2_Y_i"         , f2_Y_i            , "f2_Y_i[N_glob_tracks]/F");
	tree_out.Branch("f2_Z_i"         , f2_Z_i            , "f2_Z_i[N_glob_tracks]/F");

	tree_out.Branch("f15_X_i"        , f15_X_i           , "f15_X_i[N_glob_tracks]/F");
	tree_out.Branch("f15_Y_i"        , f15_Y_i           , "f15_Y_i[N_glob_tracks]/F");
	tree_out.Branch("f15_Z_i"        , f15_Z_i           , "f15_Z_i[N_glob_tracks]/F");

	tree_out.Branch("f16_X_i"        , f16_X_i           , "f16_X_i[N_glob_tracks]/F");
	tree_out.Branch("f16_Y_i"        , f16_Y_i           , "f16_Y_i[N_glob_tracks]/F");
	tree_out.Branch("f16_Z_i"        , f16_Z_i           , "f16_Z_i[N_glob_tracks]/F");

	tree_out.Branch("rpc_X_i"        , rpc_X_i          , "rpc_X_i[N_glob_tracks]/F");
	tree_out.Branch("rpc_Y_i"        , rpc_Y_i          , "rpc_Y_i[N_glob_tracks]/F");
	tree_out.Branch("rpc_Z_i"        , rpc_Z_i          , "rpc_Z_i[N_glob_tracks]/F");
	tree_out.Branch("rpc_T_i"        , rpc_T_i          , "rpc_T_i[N_glob_tracks]/F");

	tree_out.Branch("f1_X"           , f1_X           , "f1_X[N_glob_tracks]/F");
	tree_out.Branch("f1_Y"           , f1_Y           , "f1_Y[N_glob_tracks]/F");
	tree_out.Branch("f1_Z"           , f1_Z           , "f1_Z[N_glob_tracks]/F");

	tree_out.Branch("f2_X"           , f2_X           , "f2_X[N_glob_tracks]/F");
	tree_out.Branch("f2_Y"           , f2_Y           , "f2_Y[N_glob_tracks]/F");
	tree_out.Branch("f2_Z"           , f2_Z           , "f2_Z[N_glob_tracks]/F");

	tree_out.Branch("f15_X"         , f15_X           , "f15_X[N_glob_tracks]/F");
	tree_out.Branch("f15_Y"         , f15_Y           , "f15_Y[N_glob_tracks]/F");
	tree_out.Branch("f15_Z"         , f15_Z           , "f15_Z[N_glob_tracks]/F");

	tree_out.Branch("f16_X"         , f16_X           , "f16_X[N_glob_tracks]/F");
	tree_out.Branch("f16_Y"         , f16_Y           , "f16_Y[N_glob_tracks]/F");
	tree_out.Branch("f16_Z"         , f16_Z           , "f16_Z[N_glob_tracks]/F");

	tree_out.Branch("rpc_X"        , rpc_X          , "rpc_X[N_glob_tracks]/F");
	tree_out.Branch("rpc_Y"        , rpc_Y          , "rpc_Y[N_glob_tracks]/F");
	tree_out.Branch("rpc_Z"        , rpc_Z          , "rpc_Z[N_glob_tracks]/F");
	tree_out.Branch("rpc_T"        , rpc_T          , "rpc_T[N_glob_tracks]/F");

	tree_out.Branch("tofd_Q"       , tofd_Q         , "tofd_Q[N_glob_tracks]/F");

	// Request storage of R3BTrack data in the output tree
	mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);
	return kSUCCESS;
}

void R3BTrackingS522::Exec(Option_t* option){

	if (fNEvents / 1000. == (int)fNEvents / 1000)
		std::cout << "\rEvents: " << fNEvents \
			<< " / " << maxevent << " ("  \
			<< (int)(fNEvents * 100. / maxevent)
			<< " %) " << std::flush;
	fNEvents += 1;
	N_glob_tracks=0;
	N_in_tracks=0;
	N_out_tracks=0;
	good_ev=0;

	Tpat = fHeader->GetTpat();//vairable in the output tree
	if(Tpat==0 || Tpat>64) return;//if Tpat is not set
	mul_m0   = fDataItems[MWPC0_HITDATA]->GetEntries();
	mul_m1   = fDataItems[MWPC1_HITDATA]->GetEntries();
	mul_los  = fDataItems[DET_LOS]->GetEntries();
	mul_foot = fDataItems[FOOT_HITDATA]->GetEntries();
	mul_rpc  = fDataItems[DET_RPC]->GetEntries();
	mul_tofd = fDataItems[DET_TOFD]->GetEntries();

	// cout << "\n mul_m0 : " << mul_m0; 
	// cout << "\n mul_m1 : " << mul_m1; 
	// cout << "\n mul_los : " << mul_los; 
	// cout << "\n mul_foot : " << mul_foot; 
	// cout << "\n mul_rpc : " << mul_rpc; 
	// cout << "\n mul_tofd : " << mul_tofd; 


	if(mul_m0!=1 || mul_m1!=1 || mul_foot<1 || mul_los>1 || mul_rpc==0) return;//for now take only mul=1 in mwpcs

	// cout << "\n mul_m0 after : " << mul_m0 << endl; 
	// cout << "\n mul_m1 after : " << mul_m1 << endl; 
	// cout << "\n mul_los after : " << mul_los << endl; 
	// cout << "\n mul_foot after : " << mul_foot << endl; 
	// cout << "\n mul_rpc after : " << mul_rpc << endl; 
	// cout << "\n mul_tofd after : " << mul_tofd << endl;    

	//FRS data
	auto frs_DataItems = fDataItems.at(FRS_DATA);
	if(frs_DataItems->GetEntries() < 1) return; 
	auto frs_data = static_cast<R3BFrsData*>(frs_DataItems->At(0));
	if(frs_data->GetBrho()<6.224 || frs_data->GetBrho()>6.25) return;
	if(frs_data->GetZ()<8.7 || frs_data->GetZ()>9) return;
	frs_Brho = frs_data->GetBrho();
	frs_A = frs_data->GetZ();
	////------ Get TOFD data 
	//R3BTofdHitData* tofd_hit{};
	//bool is_good_tofd = false;
	//for (auto i = 0; i < fDataItems[DET_TOFD]->GetEntriesFast(); ++i)
	//{
	// tofd_hit = static_cast<R3BTofdHitData*>(fDataItems[DET_TOFD]->At(i));
	// if (tofd_hit->GetDetId() == 1) // only hits from first plane
	// {
	//  is_good_tofd = true;
	//  break;
	// }
	//}
	//if(!is_good_tofd){ 
	// //cout<<"Bad_tofd"<<endl; 
	// return;
	//}
	if(!MakeIncomingTracks()){ 
		//cout<<"Bad_incoming"<<endl; 
		return;//at least one good track candidate in FOOT
	}
	if(!MakeOutgoingTracks()){ 
		//cout<<"Bad_outgoing"<<endl; 
		return;//at least one good track candidate in Rpc
	}
	if(mul_f1>1 || mul_f2>1 || mul_f15>1 || mul_f16>1) return;
	good_ev=1;
	double a=0,b=0,c=0;
	double delta_TY0, delta_TX0;
	//cout << "new event" << endl;
	for (auto & tin : tracks_in){
		a++;
		b=0;
		for (auto & tout : tracks_out){
			b++;
			//preserve the order, it is expected by the MDF function!
			mdf_data[0] = tin.f16_y;
			mdf_data[1] = tin.f16_z;
			mdf_data[2] = tin.f15_x;
			mdf_data[3] = tin.f15_z;
			mdf_data[4] = tout.rpc_x;
			mdf_data[5] = tout.rpc_y;
			mdf_data[6] = tout.rpc_z;
			mdf_data[7] = tout.rpc_tof;

			// Calculate all required MDF values

			TVector3 vertex_foot;//projection to the center of the target (0,0,0)

			vertex_foot.SetX(tin.f15_x - MDF_TX0->MDF(mdf_data) * tin.f15_z);
			vertex_foot.SetY(tin.f16_y - MDF_TY0->MDF(mdf_data) * tin.f16_z);
			vertex_foot.SetZ(0);

			double dx = vertex_foot.X() - tin.vertex_mwpc_X;
			double dy = vertex_foot.Y() - tin.vertex_mwpc_Y;
			//   cout << " delta : " << dx << " " << dy << endl;
			//   cout << " counters : " << a << " " << b << " " << c << endl;
			//   cout << " Foot : "  << tin.f16_y << " " << tin.f16_z << " " << tin.f15_x << " " << tin.f15_z << endl;
			//   cout << " RPC : " << tout.rpc_x << " " << tout.rpc_y << " " << tout.rpc_z << " " << tout.rpc_tof << endl;
			//   cout << " vertex_foot : "  <<  vertex_foot.X() << " " << vertex_foot.Y() << endl;
			//   cout << " vertex_mwpc : "  <<  tin.vertex_mwpc_X << " " << tin.vertex_mwpc_Y << endl;
			//if(dy <-.2 || dy>.4){continue;}
			//if(dx <-0.8 || dx>1.3){continue;}
			vertex_mwpc_X[N_glob_tracks] = tin.vertex_mwpc_X;
			vertex_mwpc_Y[N_glob_tracks] = tin.vertex_mwpc_Y;
			vertex_foot_X[N_glob_tracks] = vertex_foot.X();
			vertex_foot_Y[N_glob_tracks] = vertex_foot.Y();
			dx_vertex[N_glob_tracks] = dx;
			dy_vertex[N_glob_tracks] = dy;

			PoQ[N_glob_tracks] = MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent;
			FlightPath[N_glob_tracks] = MDF_FlightPath->MDF(mdf_data);
			TX0[N_glob_tracks] = MDF_TX0->MDF(mdf_data);
			TY0[N_glob_tracks] = MDF_TY0->MDF(mdf_data);
			Beta[N_glob_tracks] = FlightPath[N_glob_tracks] / rpc_T[N_glob_tracks] / SPEED_OF_LIGHT;
			Gamma[N_glob_tracks] = 1. / sqrt(1 - pow(Beta[N_glob_tracks], 2));
			mdf_AoZ[N_glob_tracks] = PoQ[N_glob_tracks] / Beta[N_glob_tracks] / Gamma[N_glob_tracks] / AMU;

			f15_X_i[N_glob_tracks]   = f15_point_i.X();
			f15_Z_i[N_glob_tracks]   = f15_point_i.Z();

			f16_Y_i[N_glob_tracks]   = f16_point_i.Y();
			f16_Z_i[N_glob_tracks]   = f16_point_i.Z();

			rpc_X_i[N_glob_tracks]   = rpc_point_i.X();
			rpc_Y_i[N_glob_tracks]   = rpc_point_i.Y();
			rpc_Z_i[N_glob_tracks]   = rpc_point_i.Z();
			rpc_T_i[N_glob_tracks]   = rpc_tof_i;
	
			f1_X[N_glob_tracks]   = tin.f1_x;
			f1_Z[N_glob_tracks]   = tin.f1_z;

			f2_Y[N_glob_tracks]   = tin.f2_y;
			f2_Z[N_glob_tracks]   = tin.f2_z;

			f15_X[N_glob_tracks]   = tin.f15_x;
			f15_Z[N_glob_tracks]   = tin.f15_z;

			f16_Y[N_glob_tracks]   = tin.f16_y;
			f16_Z[N_glob_tracks]   = tin.f16_z;

			f15_Q[N_glob_tracks]   = tin.f15_Q;
			f16_Q[N_glob_tracks]   = tin.f16_Q;

			rpc_X[N_glob_tracks]   = tout.rpc_x;
			rpc_Y[N_glob_tracks]   = tout.rpc_y;
			rpc_Z[N_glob_tracks]   = tout.rpc_z;
			rpc_T[N_glob_tracks]   = tout.rpc_tof;
			//tofd_Q[N_glob_tracks] = tofd_hit->GetEloss();

			TVector3 rpc(rpc_X[N_glob_tracks], rpc_Y[N_glob_tracks], rpc_Z[N_glob_tracks]);
			TVector3 vec_PoQ(TX0[N_glob_tracks], TY0[N_glob_tracks], 1);
			vec_PoQ.SetMag(PoQ[N_glob_tracks]);

			if (DoAlignment && mul_f15==1 && mul_f16==1 && mul_rpc==1 && PoQ[N_glob_tracks]>1.744 && PoQ[N_glob_tracks]<1.750 && counter_align<100)
			{
				counter_align++;
				det_points align_data;
				align_data.f15.SetXYZ(f15_point_i.X(),0, f15_point_i.Z());
				align_data.f16.SetXYZ(0, f16_point_i.Y(), f16_point_i.Z());
				align_data.frpc.SetXYZ(rpc_point_i.X(), rpc_point_i.Y(), rpc_point_i.Z());
				align_data.rpc_tof = rpc_tof_i;
				det_points_vec.push_back(align_data);
			}
			N_glob_tracks++;
			if(N_glob_tracks == N_glob_tracks_max) return;
		}
	}
	return;
}

void R3BTrackingS522::FinishEvent()
{
	for (auto& DataItem : fDataItems)
	{
		DataItem->Clear();
	}
	fTrackItems->Clear();
	f15_hits.clear();
	f16_hits.clear();
	tracks_in.clear();
	tracks_out.clear();
	if(good_ev==1){tree_out.Fill();}
}

void R3BTrackingS522::FinishTask()
{
	LOG(info) << "Processed " << fNEvents << " events\n\n";
	if(DoAlignment){
		Alignment();
	}
	tree_out.Write();
	//cout << "\n\n------- Statisitcs summary --------- ";
}

bool R3BTrackingS522::IsGoodFootHit(R3BFootHitData* fhit)
{
	if(fhit->GetEnergy()>FootEnergyMin && fhit->GetEnergy()<FootEnergyMax)
		return true;
	else 
		return false;
}

bool R3BTrackingS522::SortFootData()
{
	f15_hits.clear(); f16_hits.clear(); f2_hits.clear();f1_hits.clear();
	mul_foot = fDataItems[FOOT_HITDATA]->GetEntriesFast();
	if(mul_foot==0) return false;
	for (auto f=0; f<mul_foot; ++f)
	{
		auto foot_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f));
		if(!IsGoodFootHit(foot_hit)) continue;
		switch(foot_hit->GetDetId())
		{
			case 1 ://Second FOOT1 for (Y)
				f1_hits.push_back(f);
				break;
			case 2 ://First FOOT2 (X)
				f2_hits.push_back(f);
				break;
			case 15 ://Second FOOT15 for (Y)
				f15_hits.push_back(f);
				break;
			case 16 ://First FOOT16 (X)
				f16_hits.push_back(f);
				break;
		}
	}
	if(f1_hits.empty() || f2_hits.empty() || f15_hits.empty() || f16_hits.empty())     
		return false; 
	mul_f1 = f1_hits.size(); mul_f2 = f2_hits.size(); mul_f15 = f15_hits.size(); mul_f16  = f16_hits.size();
	return true;
}

bool R3BTrackingS522::MakeOutgoingTracks(){
	if(fDataItems[DET_RPC]->GetEntriesFast()==0 || fDataItems[DET_RPC]->GetEntriesFast()==0){
		return false;
	}
	tracks_out.clear();
	Track tr;
	N_out_tracks=0;
	for (auto i=0; i<fDataItems[DET_RPC]->GetEntriesFast(); ++i){
		auto rpc_hit = static_cast<R3BRpcHitData*>(fDataItems[DET_RPC]->At(i));
		if(rpc_hit->GetDetId()!=0 || rpc_hit->GetChannelId()!=21) continue;
		rpc_point.SetXYZ(-1 * (rpc_hit->GetPos()*0.1 -75) ,rpc_hit->GetChannelId()*3. - 21.*3., 0); //cm
		rpc_point_i = rpc_point;
		TransformPoint(rpc_point, &rpc_angles, &rpc_position);
		tr.rpc_x = rpc_point.X();
		tr.rpc_y = rpc_point.Y();
		tr.rpc_z = rpc_point.Z();
		tr.rpc_tof = (rpc_hit->GetTof() + 28.1768 + 29.1868 );
		rpc_tof_i = tr.rpc_tof;
		//tr.rpc_tof = (rpc_hit->GetTof() + 42.2 + 22.8272);
		N_out_tracks++;
		tracks_out.push_back(tr);
	}
	if(tracks_out.empty()) return false;
	return true;
}

bool R3BTrackingS522::MakeIncomingTracks()
{
	if(fDataItems[MWPC0_HITDATA]->GetEntriesFast()==0 || fDataItems[MWPC1_HITDATA]->GetEntriesFast()==0){ 
		return false;
	}
	tracks_in.clear();
	N_in_tracks=0;
	if(!SortFootData()) return false;//at least 1 hit in every FOOT
	TVector3 vertex_mwpc;//projection to the center of the target (0,0,0)
	double tx_in={}, ty_in={};
	Track tr;
	//Get MWPC hits, for now only first hit
	auto m0_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC0_HITDATA]->At(0));
	auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
	m0_point.SetXYZ(m0_hit->GetX()*0.1, m0_hit->GetY()*0.1, 0.);//cm
	m1_point.SetXYZ(m1_hit->GetX()*0.1, m1_hit->GetY()*0.1, 0.);//cm
	TransformPoint(m0_point, &m0_angles, &m0_position);//lab
	TransformPoint(m1_point, &m1_angles, &m1_position);//lab
	//Fill output tree variables
	m0_X = m0_point.X();  m0_Y = m0_point.Y();  m0_Z = m0_point.Z();
	m1_X = m1_point.X();  m1_Y = m1_point.Y();  m1_Z = m1_point.Z();
	//cout<<"m0_x--> "<<m0_X<<endl;
	//------- Project mwpc track to the center of the target
	tx_in = (m0_point.X() - m1_point.X())/(m0_point.Z() - m1_point.Z());
	ty_in = (m0_point.Y() - m1_point.Y())/(m0_point.Z() - m1_point.Z()); 
	vertex_mwpc.SetX(m1_point.X() - tx_in * m1_point.Z());
	vertex_mwpc.SetY(m1_point.Y() - ty_in * m1_point.Z());
	vertex_mwpc.SetZ(0);
	if(abs(vertex_mwpc.X())>50 || abs(vertex_mwpc.Y())>50){return false;}
	//----- Make track candidates in FOOT and project them to the center of the target
	for (auto & f2 : f2_hits){
		auto f2_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f2));
		f2_point.SetXYZ(0, f2_hit->GetPosLab().Y() * 0.1, 0); //cm
		TransformPoint(f2_point,  &f2_angles,  &f2_position);
		for (auto & f1 : f1_hits){
			auto f1_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f1));
			f1_point.SetXYZ((-1.) * f1_hit->GetPosLab().X() * 0.1, 0, 0); //cm
			TransformPoint(f1_point,  &f1_angles,  &f1_position);
			for (auto & f16 : f16_hits){
				auto f16_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f16));
				f16_point.SetXYZ(0, f16_hit->GetPosLab().Y() * 0.1, 0); //cm
				f16_point_i=f16_point;
				TransformPoint(f16_point,  &f16_angles,  &f16_position);
				for (auto & f15 : f15_hits){
					auto f15_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f15));
					f15_point.SetXYZ((-1.) * f15_hit->GetPosLab().X() * 0.1, 0, 0); //cm
					f15_point_i=f15_point;
					TransformPoint(f15_point,  &f15_angles,  &f15_position);       
					//project foot track to the center of the target
					tr.f1_x    = f1_point.X();
					tr.f1_z    = f1_point.Z();
					tr.f2_y    = f2_point.Y();
					tr.f2_z    = f2_point.Z();
					tr.f15_x   = f15_point.X();
					tr.f15_z   = f15_point.Z();
					tr.f16_y   = f16_point.Y();
					tr.f16_z   = f16_point.Z();
					tr.vertex_mwpc_X = vertex_mwpc.X();
					tr.vertex_mwpc_Y = vertex_mwpc.Y();
					tr.f15_Q  = f15_hit->GetEnergy();
					tr.f16_Q  = f16_hit->GetEnergy();
					tracks_in.push_back(tr);
					N_in_tracks++;
					if(N_in_tracks==N_glob_tracks_max/2) return true;
				}
			}
		}
	}
	if(tracks_in.empty()) return false;
	else return true;
}


void R3BTrackingS522::TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
{
	r.SetToIdentity();
	// First Euler rotation around Y axis
	r.RotateY(rot->Y());
	// get local X axis after first rotation
	v3_localX.SetMagThetaPhi(1, r.ThetaX(), r.PhiX());
	// Second Euler rotation around local X axis
	r.Rotate(rot->X(), v3_localX);
	// get local Z axis after second rotation
	v3_localZ.SetMagThetaPhi(1, r.ThetaZ(), r.PhiZ());
	// final rotation around local Z axis
	r.Rotate(rot->Z(), v3_localZ);
	point.Transform(r);
	point += (*trans);
	return;
}

void R3BTrackingS522::TransformPoint1(TVector3& point1, TVector3 rot1, TVector3 trans1)
{
	r1.SetToIdentity();
	// First Euler rotation around Y axis
	r1.RotateY(rot1.Y());
	// get local X axis after first rotation
	v3_localX1.SetMagThetaPhi(1, r1.ThetaX(), r1.PhiX());
	// Second Euler rotation around local X axis
	r1.Rotate(rot1.X(), v3_localX1);
	// get local Z axis after second rotation	
	v3_localZ1.SetMagThetaPhi(1, r1.ThetaZ(), r1.PhiZ());
	// final rotation around local Z axis
	r1.Rotate(rot1.Z(), v3_localZ1);
	point1.Transform(r1);
	point1 += (trans1);
	return;
}

// Setup energy cuts in foot and fibers 
void R3BTrackingS522::SetFootEnergyMinMax(double min, double max)
{
	FootEnergyMin = min;
	FootEnergyMax = max;
	return;
}

void R3BTrackingS522::Alignment()
{
	int run_flag = 0;
	std::cout << "\n\n----------------------------------------------------";
	std::cout << "\n Ready for detector alignment using the following data set:";
	std::cout << "\n\tNumber of reference tracks: " << det_points_vec.size();
	std::cout << "\n\tP/Z reference: " << reference_PoQ << "GeV/c";
	std::cout << "\n\n Do you want to continue? (0 = no, 1 = yes):  ";
	std::cin >> run_flag;
	if (run_flag == 0){
		return;
	}
	// Now define minimizer
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	const UInt_t num_det = 3;//number of detector
	const UInt_t off_p_det = 3;//offsets
	const UInt_t angles_p_det = 1;//angles
	const UInt_t NVarsFunctor = 10; //(angles_p_det + off_p_det)*num_det; // number of the alignment offsets: (number of  angles + number of offsets) x number of detectors
	const double* xs={0};
	double step[NVarsFunctor]={0};
	double offset[NVarsFunctor]={0};
	Double_t mean_offset[NVarsFunctor]={0};
	Double_t max_offset[NVarsFunctor]={0};
	TH1F* hist[NVarsFunctor]={0};
	char hname[100];

	// For every detector: 3 Euler angles and 3 shifts
	for (auto d = 0; d < num_det; ++d)
	{
		if (d == 0) {
			for (auto a = 0; a < angles_p_det; ++a) // angle shifts in rad
			{ 
				mean_offset[a] = 0;
				max_offset[a] = 0.1;
				step[a] = 0.001;
				cout << "here 0 : " << a << endl;
			}

			for (auto o = 0; o < 2/*off_p_det*/; ++o) // position shifts in cm
			{
				mean_offset[angles_p_det + o] = 0;
				max_offset[angles_p_det + o] = 20;
				//if(o==1){mean_offset[angles_p_det + o] = -5;max_offset[angles_p_det + o]= 10;}
				step[angles_p_det + o] = 0.01; //Valerii
				cout << "here 1 : " << angles_p_det + o << endl;
			}
			mean_offset[angles_p_det + off_p_det -1] = 0;
			max_offset[angles_p_det + off_p_det -1] = 10;
			step[angles_p_det + off_p_det -1] = 0.001; //Valerii
			cout << "here 2 : " << angles_p_det + off_p_det -1 << endl;
		}
		else {
			for (auto o = 0; o < off_p_det; ++o) // angle shifts in rad
			{ 
				mean_offset[d * num_det + angles_p_det + o] = 0;
				max_offset[d * num_det + angles_p_det + o] = 20;
				step[d * num_det + angles_p_det + o] = 0.01;
				cout << "here 3 : " << d * num_det + angles_p_det + o << endl;
			}
		}
	}
	for (UInt_t i = 0; i < NVarsFunctor; i++)
	{
		sprintf(hname, "par%d", i);
		hist[i] = new TH1F(hname, hname, 100, mean_offset[i] - max_offset[i], mean_offset[i] + max_offset[i]);
	}
	std::cout << "\n\n-- Perfroming minimization for the detector alignment. Please wait... \n";

	// Setting up Minimizer function parameters
	Double_t precision = 1e-10; // 0 - default precision will be automaticalle determined
	Double_t tolerance = 0.01;
	minimizer->SetMaxFunctionCalls(1000000000); // for Minuit/Minuit2
	minimizer->SetMaxIterations(2000);           // for GSL
	minimizer->SetTolerance(tolerance);
	minimizer->SetPrecision(precision);
	minimizer->SetPrintLevel(2);
	minimizer->SetStrategy(0); // to run faster
	ROOT::Math::Functor f(&R3BTrackingS522::AlignmentErrorS522, NVarsFunctor);
	minimizer->SetFunction(f);
	minimizer->SetPrintLevel(0);
	// Getting experimental data and doing minimization
	UInt_t i;
	Int_t bs = 0;
	bool good_minimum = false;
	while (bs < 100) // running several minimizations
	{	
		minimizer->Clear();
		for (i = 0; i < NVarsFunctor; i++) // sampling +=50% from limits
		{
			offset[i] = gRandom->Uniform(mean_offset[i] - max_offset[i] , mean_offset[i] + max_offset[i]);
			minimizer->SetVariable(i, Form("par%d", i), offset[i], step[i]);
			minimizer->SetVariableLimits(i, mean_offset[i] - max_offset[i], mean_offset[i] + max_offset[i]);	
		}
		minimizer->Minimize();
		xs = minimizer->X();
		//cout<<"Min status--> "<<minimizer->Status()<<endl;
		//if(minimizer->Status() !=0) continue; //valid minimum
		// Check if all paramters are "far" from limits
		for (i = 0; i < NVarsFunctor; i++)
		{
			if (fabs((xs[i] - (mean_offset[i] - max_offset[i])) / (mean_offset[i] - max_offset[i])) < 0.3 ||
					fabs(((mean_offset[i] + max_offset[i]) - xs[i]) / (mean_offset[i] + max_offset[i])) < 0.3)
			{
				break;
			}
		}
		if (i != NVarsFunctor)
		{
			//cout << "parameter is to close to the limit" << endl;
			continue;
		}
		else
		{
			std::cout << "\n\n-- Good parameters are found!!";
			good_minimum = true;
		}
		std::cout << "\n-- Minimization No. " << bs << "\n-- Minimum: f(";
		for (i = 0; i < NVarsFunctor; i++)
		{
			std::cout << xs[i];
			hist[i]->Fill(xs[i]);
			if (i != NVarsFunctor - 1)
				std::cout << ",";
		}
		std::cout << ") = " << minimizer->MinValue();
		std::cout << "\n-- Minimizer status: " << minimizer->Status() << std::endl;
		cout<<"BS---> "<<bs<<endl;	
		bs++;
	}
	// Outputing final info and histograms
	if (minimizer->MinValue() < tolerance && f(xs) < tolerance)
		std::cout << "-- Minimizer "
			<< "   converged to the right minimum" << std::endl;
	else
	{
		std::cout << "-- Minimizer "
			<< "   failed to converge !!!" << std::endl;
	}
	TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
	c1->Divide(3, 4);
	for (UInt_t v = 0; v < NVarsFunctor; v++)
	{
		c1->cd(v + 1);
		hist[v]->Draw();
	}
	return;
}


double R3BTrackingS522::AlignmentErrorS522(const double* par)
{
	gMDFTrackerS522->f15_pos_offset.SetXYZ(par[4], par[5], par[6]);
	//	gMDFTrackerS522->f15_ang_offset.SetXYZ(par[3], 0 , 0);
	//
	gMDFTrackerS522->f16_pos_offset.SetXYZ(par[7], par[8], par[9]);
	//	gMDFTrackerS522->f16_ang_offset.SetXYZ(0, par[7], 0);

	gMDFTrackerS522->rpc_ang_offset.SetXYZ(0, par[0],0);
	//RPC frame X and Y
	gMDFTrackerS522->rpc_pos_offset.SetXYZ(par[1], 0, par[2]);
	gMDFTrackerS522->rpc_tof_offset = par[3];

	Double_t mdf_input[8]; // data container for the MDF function

	double v1=0,v2=0,v3=0;
	double v = 0;
	int counter = 0;
	for (auto& d : (gMDFTrackerS522->det_points_vec))
	{
		gMDFTrackerS522->f15_point_i = d.f15;
		gMDFTrackerS522->f16_point_i = d.f16;
		gMDFTrackerS522->rpc_point_i = d.frpc;
		gMDFTrackerS522->rpc_tof_i = d.rpc_tof;
		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f15_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot15() + gMDFTrackerS522->f15_ang_offset,
				gMDFTrackerS522->GetPositionFoot15() + gMDFTrackerS522->f15_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f16_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot16() + gMDFTrackerS522->f16_ang_offset,
				gMDFTrackerS522->GetPositionFoot16() + gMDFTrackerS522->f16_pos_offset);

		//	 	should apply before rotating
		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->rpc_point_i,
				gMDFTrackerS522->GetEulerAnglesRpc() + gMDFTrackerS522->rpc_ang_offset,
				gMDFTrackerS522->GetPositionRpc() + gMDFTrackerS522->rpc_pos_offset);

		mdf_input[0] = gMDFTrackerS522->f16_point_i.Y();
		mdf_input[1] = gMDFTrackerS522->f16_point_i.Z();
		mdf_input[2] = gMDFTrackerS522->f15_point_i.X();
		mdf_input[3] = gMDFTrackerS522->f15_point_i.Z();
		mdf_input[4] = gMDFTrackerS522->rpc_point_i.X();
		mdf_input[5] = gMDFTrackerS522->rpc_point_i.Y();
		mdf_input[6] = gMDFTrackerS522->rpc_point_i.Z();
		mdf_input[7] = gMDFTrackerS522->rpc_tof_i + gMDFTrackerS522->rpc_tof_offset;

		TVector3 vec_PoQ_a(gMDFTrackerS522->Get_MDF_TX0()->MDF(mdf_input),gMDFTrackerS522->Get_MDF_TY0()->MDF(mdf_input), 1);
		vec_PoQ_a.SetMag(gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input));

		//Difference between TXO and TYO and MWPC angles		

		v1 += pow((gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input) - gMDFTrackerS522->GetReferencePoQ()), 2);
		v2 += pow((vec_PoQ_a.X() - 0), 2);
		v3 += pow((vec_PoQ_a.Y() - 0), 2);

		counter++;
	}
	v1 /= counter;
	v2 /= counter;
	v3 /= counter;
	//v = sqrt(v1);
	v = sqrt(v1 + v2 + v3);
	return v;
}

ClassImp(R3BTrackingS522);
