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

#include "R3BTrackingProtons.h"
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
R3BTrackingProtons* gMDFTrackerProtons;

R3BTrackingProtons::R3BTrackingProtons()
	: R3BTrackingProtons("TrackingProtons", 1)
{
}

R3BTrackingProtons::R3BTrackingProtons(const char* name, Int_t iVerbose)
	: FairTask(name, iVerbose)
	, fTrigger(-1)
	, fTpat(-1)
	, fNEvents(0)
	, maxevent(0)
	, fTrackItems(new TClonesArray("R3BRpcTrack"))
	, FootEnergyMin(-1)
	, FootEnergyMax(-1)
	, fHeader(nullptr)
{
}

R3BTrackingProtons::~R3BTrackingProtons()
{
	if (fTrackItems)
		delete fTrackItems;
}

InitStatus R3BTrackingProtons::Init()
{
	LOG(info) << "R3BTrackingProtons::Init()";
	FairRootManager* mgr = FairRootManager::Instance();
	if (NULL == mgr)
	{
		LOG(fatal) << "FairRootManager not found";
	}
	fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
	if (!fHeader)
	{
		LOG(warn) << "R3BTrackingProtons::Init() EventHeader. not found";
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
	}

	// check if all cuts are properly set
	if (FootEnergyMin < 0 || FootEnergyMax < 0) 
	{
		R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
	}
	// Initializing all MDF functions
	LOG(info) << "Reading MDF function for FlightPath";
	MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

	LOG(info) << "Reading MDF function for TX0";
	MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

	LOG(info) << "Reading MDF function for TY0";
	MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());

	LOG(info) << "Reading MDF function for PoQ";
	MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());

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

	tree_out.Branch("vertex_mwpc_X", vertex_mwpc_X  , "vertex_mwpc_X[N_in_tracks]/F");
	tree_out.Branch("vertex_mwpc_Y", vertex_mwpc_Y  , "vertex_mwpc_Y[N_in_tracks]/F");

	tree_out.Branch("vertex_foot_X", vertex_foot_X  , "vertex_foot_X[N_in_tracks]/F");
	tree_out.Branch("vertex_foot_Y", vertex_foot_Y  , "vertex_foot_Y[N_in_tracks]/F");
	tree_out.Branch("dx_vertex"    , dx_vertex      , "dx_vertex[N_in_tracks]/F");
	tree_out.Branch("dy_vertex"    , dy_vertex      , "dy_vertex[N_in_tracks]/F");

	tree_out.Branch("f1_Q"         , f1_Q           , "f1_Q[N_glob_tracks]/F");
	tree_out.Branch("f2_Q"         , f2_Q           , "f2_Q[N_glob_tracks]/F");

	tree_out.Branch("f1_X"         , f1_X           , "f1_X[N_glob_tracks]/F");
	tree_out.Branch("f1_Y"         , f1_Y           , "f1_Y[N_glob_tracks]/F");
	tree_out.Branch("f1_Z"         , f1_Z           , "f1_Z[N_glob_tracks]/F");

	tree_out.Branch("f2_X"         , f2_X           , "f2_X[N_glob_tracks]/F");
	tree_out.Branch("f2_Y"         , f2_Y           , "f2_Y[N_glob_tracks]/F");
	tree_out.Branch("f2_Z"         , f2_Z           , "f2_Z[N_glob_tracks]/F");

	tree_out.Branch("rpc_X"        , rpc_X          , "rpc_X[N_glob_tracks]/F");
	tree_out.Branch("rpc_Y"        , rpc_Y          , "rpc_Y[N_glob_tracks]/F");
	tree_out.Branch("rpc_Z"        , rpc_Z          , "rpc_Z[N_glob_tracks]/F");
	tree_out.Branch("rpc_T"        , rpc_T          , "rpc_T[N_glob_tracks]/F");

	mgr->Register("MDFRpcTracks", "MDFRpcTracks data", fTrackItems, kTRUE);
	return kSUCCESS;
}

void R3BTrackingProtons::Exec(Option_t* option){

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

	if(mul_m0!=1 || mul_m1!=1 || mul_foot<1 || mul_los != 1 || mul_rpc==0) return;//for now take only mul=1 in mwpcs

	//FRS data
	//auto frs_DataItems = fDataItems.at(FRS_DATA);
	//if(frs_DataItems->GetEntries() < 1) return; 
	//auto frs_data = static_cast<R3BFrsData*>(frs_DataItems->At(0));
	//cout << frs_data->GetBrho() << " " << frs_data->GetZ() << endl;
	//if(frs_data->GetBrho()<17.5 || frs_data->GetBrho()>17.7) return;
	//if(frs_data->GetZ()<5.2 || frs_data->GetZ()>6.7) return;
	if(!MakeOutgoingTracks()){ 
	//	cout<<"Bad_outgoing"<<endl; 
		return;//at least one good track candidate in Rpc
	}
	if(!MakeIncomingTracks()){ 
	//	cout<<"Bad_incoming"<<endl; 
		return;//at least one good track candidate in FOOT
	}
	//if(mul_f1>5 || mul_f2>5) return;
	good_ev=1;
	for (auto & tin : tracks_in){
		for (auto & tout : tracks_out){
			mdf_data[0] = tin.f2_y;
			mdf_data[1] = tin.f2_z;
			mdf_data[2] = tin.f1_x;
			mdf_data[3] = tin.f1_z;
			mdf_data[4] = tout.rpc_x;
			mdf_data[5] = tout.rpc_y;
			mdf_data[6] = tout.rpc_z;
			mdf_data[7] = tout.rpc_tof;
			// Calculate all required MDF values

			TVector3 vertex_foot;//projection to the center of the target (0,0,0)

			vertex_foot.SetX(tin.f1_x - MDF_TX0->MDF(mdf_data) * tin.f1_z);
			vertex_foot.SetY(tin.f2_y - MDF_TY0->MDF(mdf_data) * tin.f2_z);
			vertex_foot.SetZ(0);

			double dx = vertex_foot.X() - tin.vertex_mwpc_X;
			double dy = vertex_foot.Y() - tin.vertex_mwpc_Y;
			vertex_mwpc_X[N_glob_tracks] = tin.vertex_mwpc_X;
			vertex_mwpc_Y[N_glob_tracks] = tin.vertex_mwpc_Y;
			vertex_foot_X[N_glob_tracks] = vertex_foot.X();
			vertex_foot_Y[N_glob_tracks] = vertex_foot.Y();
			dx_vertex[N_glob_tracks] = dx;
			dy_vertex[N_glob_tracks] = dy;

			PoQ[N_glob_tracks] = MDF_PoQ->MDF(mdf_data);
			FlightPath[N_glob_tracks] = MDF_FlightPath->MDF(mdf_data);
			TX0[N_glob_tracks] = MDF_TX0->MDF(mdf_data);
			TY0[N_glob_tracks] = MDF_TY0->MDF(mdf_data);
			Beta[N_glob_tracks] = FlightPath[N_glob_tracks] / rpc_T[N_glob_tracks] / SPEED_OF_LIGHT;
			Gamma[N_glob_tracks] = 1. / sqrt(1 - pow(Beta[N_glob_tracks], 2));
			mdf_AoZ[N_glob_tracks] = PoQ[N_glob_tracks] / Beta[N_glob_tracks] / Gamma[N_glob_tracks] / AMU;

			f1_X[N_glob_tracks]   = tin.f1_x;
			f1_Z[N_glob_tracks]   = tin.f1_z;

			f2_Y[N_glob_tracks]   = tin.f2_y;
			f2_Z[N_glob_tracks]   = tin.f2_z;

			f1_Q[N_glob_tracks]   = tin.f1_Q;
			f2_Q[N_glob_tracks]   = tin.f2_Q;

			rpc_X[N_glob_tracks]   = tout.rpc_x;
			rpc_Y[N_glob_tracks]   = tout.rpc_y;
			rpc_Z[N_glob_tracks]   = tout.rpc_z;
			rpc_T[N_glob_tracks]   = tout.rpc_tof;

			TVector3 rpc(rpc_X[N_glob_tracks], rpc_Y[N_glob_tracks], rpc_Z[N_glob_tracks]);
			TVector3 vec_PoQ(TX0[N_glob_tracks], TY0[N_glob_tracks], 1);
			vec_PoQ.SetMag(PoQ[N_glob_tracks]);

			AddTrackData(rpc, dx, dy, TX0[N_glob_tracks], TY0[N_glob_tracks], vec_PoQ, Beta[N_glob_tracks],  Gamma[N_glob_tracks],  FlightPath[N_glob_tracks], mdf_AoZ[N_glob_tracks], N_glob_tracks);

			N_glob_tracks++;
			if(N_glob_tracks == N_glob_tracks_max) return;
		}
	}
	return;
}

void R3BTrackingProtons::FinishEvent()
{
	for (auto& DataItem : fDataItems)
	{
		DataItem->Clear();
	}
	fTrackItems->Clear();
	f1_hits.clear();
	f2_hits.clear();
	tracks_in.clear();
	tracks_out.clear();
	if(good_ev==1){tree_out.Fill();}
}

void R3BTrackingProtons::FinishTask()
{
	LOG(info) << "Processed " << fNEvents << " events\n\n";
	tree_out.Write();
	//cout << "\n\n------- Statisitcs summary --------- ";
}

bool R3BTrackingProtons::IsGoodFootHit(R3BFootHitData* fhit)
{
	if(fhit->GetEnergy()>FootEnergyMin && fhit->GetEnergy()<FootEnergyMax)
		return true;
	else 
		return false;
}

bool R3BTrackingProtons::SortFootData()
{
	f1_hits.clear(); f2_hits.clear(); 
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
		}
	}
	if(f1_hits.empty() || f2_hits.empty())     
		return false; 
	mul_f1  = f1_hits.size(); mul_f2  = f2_hits.size();
	return true;
}

bool R3BTrackingProtons::MakeOutgoingTracks(){
	tracks_out.clear();
	Track tr;
	N_out_tracks=0;
	for (auto i=0; i<fDataItems[DET_RPC]->GetEntriesFast(); ++i){
		auto rpc_hit = static_cast<R3BRpcHitData*>(fDataItems[DET_RPC]->At(i));
		if(rpc_hit->GetDetId()!=0) continue;
		rpc_point.SetXYZ(-1 * (rpc_hit->GetPos()*0.1 -75) ,rpc_hit->GetChannelId()*3. - 21.*3., 0); //cm
		TransformPoint(rpc_point, &rpc_angles, &rpc_position);
		tr.rpc_x = rpc_point.X();
		tr.rpc_y = rpc_point.Y();
		tr.rpc_z = rpc_point.Z();

		tr.rpc_tof = (rpc_hit->GetTof() + 42.2 + 23.1578) + tof_offset;
		N_out_tracks++;
		tracks_out.push_back(tr);
	}
	if(tracks_out.empty()) return false;
	return true;
}

bool R3BTrackingProtons::MakeIncomingTracks()
{
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
	//------- Project mwpc track to the center of the target
	tx_in = (m0_point.X() - m1_point.X())/(m0_point.Z() - m1_point.Z());
	ty_in = (m0_point.Y() - m1_point.Y())/(m0_point.Z() - m1_point.Z()); 
	vertex_mwpc.SetX(m1_point.X() - tx_in * m1_point.Z());
	vertex_mwpc.SetY(m1_point.Y() - ty_in * m1_point.Z());
	vertex_mwpc.SetZ(0);
	//----- Make track candidates in FOOT and project them to the center of the target
	for (auto & f2 : f2_hits){
		auto f2_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f2));
		f2_point.SetXYZ(0, f2_hit->GetPosLab().Y() * 0.1, 0); //cm
		f2_point_i=f2_point;
		TransformPoint(f2_point,  &f2_angles,  &f2_position);
		for (auto & f1 : f1_hits){
			auto f1_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f1));
			f1_point.SetXYZ(f1_hit->GetPosLab().X() * 0.1, 0, 0); //cm
			f1_point_i=f1_point;
			TransformPoint(f1_point,  &f1_angles,  &f1_position);       
			//project foot track to the center of the target
			tr.f1_x   = f1_point.X();
			tr.f1_z   = f1_point.Z();
			tr.f2_y   = f2_point.Y();
			tr.f2_z   = f2_point.Z();
			tr.vertex_mwpc_X = vertex_mwpc.X();
			tr.vertex_mwpc_Y = vertex_mwpc.Y();
			tr.f1_Q  = f1_hit->GetEnergy();
			tr.f2_Q  = f2_hit->GetEnergy();
			tracks_in.push_back(tr);
			N_in_tracks++;
			if(N_in_tracks==N_glob_tracks_max/2) return true;
		}
	}
	if(tracks_in.empty()) return false;
	//if(abs(vertex_mwpc.X()>50) || abs(vertex_mwpc.Y()>50)) return false;
	else return true;
}


void R3BTrackingProtons::TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
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

R3BRpcTrack* R3BTrackingProtons::AddTrackData(TVector3 mw, Double_t dx, Double_t dy, Double_t TX, Double_t TY, TVector3 poq, Double_t beta, Double_t gamma, Double_t flightPath, Double_t AoZ, Double_t mul)
{
	// Filling output track info
	TClonesArray& clref = *fTrackItems;
	Int_t size = clref.GetEntriesFast();
	return new (clref[size]) R3BRpcTrack(mw.X(), mw.Y(), mw.Z(), dx, dy, TX, TY, poq, beta, gamma, flightPath, AoZ, mul);
}

// Setup energy cuts in foot and fibers 
void R3BTrackingProtons::SetFootEnergyMinMax(double min, double max)
{
	FootEnergyMin = min;
	FootEnergyMax = max;
	return;
}

ClassImp(R3BTrackingProtons);
