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
#include "R3BFiberMAPMTHitData.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BFrsData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BTrackS522.h"

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
#include <R3BCoarseTimeStitch.h>

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
	, fTrackItems(new TClonesArray("R3BTrackS522"))
	, reference_PoQ(0.)
	, GladCurrent(-1)
	, GladReferenceCurrent(-1)
	, FootEnergyMin(-1)
	, FootEnergyMax(-1)
	, FiberEnergyMin(-1)
	, FiberEnergyMax(-1)
	  , fHeader(nullptr)
{
}

R3BTrackingS522::~R3BTrackingS522()
{
	if (fTrackItems){
		delete fTrackItems;}
}

InitStatus R3BTrackingS522::Init()
{
	LOG(info) << "R3BTrackingS522::Init()";
	FairRootManager* mgr = FairRootManager::Instance();
	if (NULL == mgr)
	{
		LOG(fatal) << "FairRootManager not found";
	}
	fTimeStitch = new R3BCoarseTimeStitch();
	fTimeStitch->SetClockTDC150();
	fHeader = (R3BEventHeader*)mgr->GetObject("EventHeader.");
	//fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
	if (!fHeader)
	{
		LOG(warn) << "R3BTrackingS522::Init() EventHeader. not found";
	}
	// Reading all detector branches
	cout << "\nDET_MAX = " << DET_MAX << endl;
	assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
	LOG(info) << "Reading " << NOF_FIB_DET << " fiber detectors";
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
			FootEnergyMin < 0 || FootEnergyMax < 0 || 
			FiberEnergyMin < 0 || FiberEnergyMax < 0)
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

	LOG(info) << "Reading MDF function for TX1";
	MDF_TX1 = new R3BMDFWrapper(MDF_TX1_filename.Data());

	LOG(info) << "Reading MDF function for TY1";
	MDF_TY1 = new R3BMDFWrapper(MDF_TY1_filename.Data());

	//Read output from the vertex macro
	// linking to global pointer (needed by alignment)
	gMDFTrackerS522 = this;
	// Request storage of R3BTrack data in the output tree
	mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);
	return kSUCCESS; 
}

void R3BTrackingS522::Exec(Option_t* option)
{
	if (fNEvents / 1000. == (int)fNEvents / 1000)
		std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
			<< " %) " << std::flush;
	//	FairRootManager* mgr = FairRootManager::Instance();
	//	R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");
	N_glob_tracks=0;
	fNEvents += 1;
	N_in_tracks=0;
	N_out_tracks=0;
	is_good_event = false;

	Tpat = fHeader->GetTpat();//vairable in the output tree
	if(Tpat>=64 || Tpat==0) return;//if Tpat is not set	

	mul_los=-999;
	mul_m0=-999;
	mul_m1=-999;
	mul_foot=-999;
	mul_f32=-999;
	mul_f30=-999;
	mul_f31=-999;
	mul_f33=-999;
	mul_tofd=-999;
	cond=false;

	mul_los   = fDataItems[LOS_DATA]->GetEntriesFast();
	mul_m0   = fDataItems[MWPC0_HITDATA]->GetEntriesFast();
	mul_m1   = fDataItems[MWPC1_HITDATA]->GetEntriesFast();
	mul_foot = fDataItems[FOOT_HITDATA]->GetEntriesFast();
	mul_f32  = fDataItems[DET_FI32]->GetEntriesFast();
	mul_f30  = fDataItems[DET_FI30]->GetEntriesFast();
	mul_f31  = fDataItems[DET_FI31]->GetEntriesFast();
	mul_f33  = fDataItems[DET_FI33]->GetEntriesFast();
	mul_tofd = fDataItems[DET_TOFD]->GetEntriesFast();

	if(mul_m0!=1 || mul_m1!=1 || mul_foot<1) return;//for now take only mul=1 in mwpcs
	if(mul_f32<1 || mul_f30<1 || (mul_f31==0 && mul_f33==0) ||  mul_tofd<1) return;
	if(mul_los!=1) return;
	//FRS data
	//auto frs_DataItems = fDataItems.at(FRS_DATA);
	//if(frs_DataItems->GetEntriesFast() < 1) return; 
	//auto frs_data = (R3BFrsData*)frs_DataItems->At(0);
	//if(frs_data->GetBrho()<17 || frs_data->GetBrho()>18) return;
	//if(frs_data->GetAq()<2.675 || frs_data->GetAq()>2.694) return;
	//if(frs_data->GetZ()<5.2 || frs_data->GetZ()>6.7) return;
	//cout << "\nGood event!\n";
	//------ Get TOFD data 
	dd++;
	R3BTofdHitData* tofd_hit{};
	int mul_tofd1=0;
	double tofdq_temp =0;
	bool is_good_tofd = false;
	for (auto i = 0; i < fDataItems[DET_TOFD]->GetEntriesFast(); ++i)
	{
		tofd_hit = static_cast<R3BTofdHitData*>(fDataItems[DET_TOFD]->At(i));
		if (tofd_hit->GetDetId() == 1 && tofd_hit->GetTof() > 22 && tofd_hit->GetTof() < 32) // only hits from first plane, add Z later
		{
			is_good_tofd = true;
			mul_tofd1++;
			tofdq_temp=tofd_hit->GetEloss();
			//break;
		}
	}
	if(!is_good_tofd){ 
		return;
	}
	if(!MakeIncomingTracks()){ 
		return;//at least one good track candidate in FOOT
	}
	ee++;
	if(!MakeOutgoingTracks()){ 
		return;//at least one good track candidate in Fibers
	}
	//cout << "\nGood event\n";
	is_good_event = true;
	ff++;
	cond=true;
	double delta_TX1, delta_TX0, delta_TY0;
	for (auto & tin : tracks_in){
		for (auto & tout : tracks_out){
			if(vertex_mwpc.X()>-15 && vertex_mwpc.X()<15 && vertex_mwpc.Y()<10 && vertex_mwpc.Y()>-10 && mul_tofd1==1){
				//preserve the order, it is expected by the MDF function!
				mdf_data[0] = tin.f2_y;
				mdf_data[1] = tin.f2_z;
				mdf_data[2] = tin.f1_x;
				mdf_data[3] = tin.f1_z;
				mdf_data[4] = tout.f32_x;
				mdf_data[5] = tout.f32_z;
				mdf_data[6] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
				mdf_data[7] = (tout.f30_y  - tin.f2_y)/(tout.f30_z - tin.f2_z);
				// Calculate all required MDF values
				PoQ[N_glob_tracks] = MDF_PoQ->MDF(mdf_data) * (GladCurrent / GladReferenceCurrent);	
				poqm.push_back(MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent);
				FlightPath[N_glob_tracks] = MDF_FlightPath->MDF(mdf_data);
				flight_p.push_back(MDF_FlightPath->MDF(mdf_data));
				TX0[N_glob_tracks] = MDF_TX0->MDF(mdf_data);
				TX1[N_glob_tracks] = MDF_TX1->MDF(mdf_data);
				TY0[N_glob_tracks] = MDF_TY0->MDF(mdf_data);
				TY1[N_glob_tracks] = MDF_TY1->MDF(mdf_data);
				ToF[N_glob_tracks] = FlightPath[N_glob_tracks] / 0.904252 / SPEED_OF_LIGHT;
				Beta[N_glob_tracks] = FlightPath[N_glob_tracks] / ToF[N_glob_tracks] / SPEED_OF_LIGHT;
				Gamma[N_glob_tracks] = 1. / sqrt(1 - pow(Beta[N_glob_tracks], 2));
				mdf_AoZ[N_glob_tracks] = PoQ[N_glob_tracks] / Beta[N_glob_tracks] / Gamma[N_glob_tracks] / AMU;
				tof.push_back(FlightPath[N_glob_tracks] / 0.904252 / SPEED_OF_LIGHT);

				delta_TX1 = TX1[N_glob_tracks] - (tout.last_x-tout.f32_x)/(tout.last_z - tout.f32_z);//fib31
				//delta_TX1 = TX1[N_glob_tracks] - (tout.f33_x-tout.f32_x)/(tout.f33_z +tout.f32_z);//fib33
				//		if(delta_TX1>(-0.00014) && delta_TX1<(-0.00013)){
				//	if(delta_TY0<(-0.1) || delta_TY0>(0.1))
				//if(delta_TY0<(0.014) || delta_TY0>(0.02))//f31 run 159
				//		continue;

				//delta_TX0 = TX0[N_glob_tracks] - (vtx.at(fNEvents-1)*0.1 - m0_point.X())/((vtz.at(fNEvents-1)-441)*0.1 - m0_point.Z());
				//cout<<delta_TX0<<endl;
				//if(delta_TX0<(-0.029) || delta_TX0>(-0.023))//f31
				//	if(delta_TX0<(-0.16) || delta_TX0>(-0.03))//f33 vertex
				//	if(delta_TX0<(-0.02) || delta_TX0>(-0.012))//f33
				//if(delta_TX0<(-0.022) || delta_TX0>(-0.012))//f31 run 159
				//		continue;

				//		delta_TY0 = TY0[N_glob_tracks] - (tin.f16_y-tin.f2_y)/(tin.f16_z - tin.f2_z);
				//if(delta_TY0>(-0.02) && delta_TY0<(0.02)){
				//if(delta_TY0<(0.014) || delta_TY0>(0.02))//f31 run 159
				//continue;

				//		delta_TX0 = TX0[N_glob_tracks] - (tin.f15_x - tin.f1_x)/(tin.f15_z-tin.f1_z);
				//if(delta_TX0>(-0.045) && delta_TX0<(0.)){ //f31
				//if(delta_TX0<(-0.02) || delta_TX0>(-0.012))//f33
				//if(delta_TX0<(-0.022) || delta_TX0>(-0.012))//f31 run 159
				//continue;
				/*
				   f1_X[N_glob_tracks]   = tin.f1_x;
				   f1_Z[N_glob_tracks]   = tin.f1_z;

				   f2_Y[N_glob_tracks]   = tin.f2_y;
				   f2_Z[N_glob_tracks]   = tin.f2_z;

				   f15_X[N_glob_tracks]   = tin.f15_x;
				   f15_Z[N_glob_tracks]   = tin.f15_z;

				   f16_Y[N_glob_tracks]   = tin.f16_y;
				   f16_Z[N_glob_tracks]   = tin.f16_z;
				   */	   
				f30_Y[N_glob_tracks]  = tout.f30_y;
				f30_Z[N_glob_tracks]  = tout.f30_z;

				f32_X[N_glob_tracks]  = tout.f32_x;
				f32_Z[N_glob_tracks]  = tout.f32_z;

				last_X[N_glob_tracks] = tout.last_x;
				last_Z[N_glob_tracks] = tout.last_z;

				tofd_Q[N_glob_tracks] = tofdq_temp;

				TVector3 foot1(tin.f1_x, tin.f1_z,tin.f1_Q);
				TVector3 foot2(tin.f2_y, tin.f2_z,tin.f2_Q);
				TVector3 foot15(tin.f15_x, tin.f15_z,tin.f15_Q);
				TVector3 foot16(tin.f16_y, tin.f16_z,tin.f16_Q);

				TVector3 a1_foot(tin.f15_x, tin.f15_z,0.);
				TVector3 a2_foot(tin.f16_y, tin.f16_z,0.);
				TVector3 vec_foot(tin.tx0, tin.ty0,1);
				TVector3 vec_PoQ(TX0[N_glob_tracks], TY0[N_glob_tracks], 1);
				TVector3 mul_vec(N_in_tracks,N_out_tracks,N_glob_tracks);
				vec_PoQ.SetMag(PoQ[N_glob_tracks]);
				TVector3 vert_foot(vertex_mwpc.X(), vertex_mwpc.Y(), 1.);
				double opan_ang=f16_Z[N_glob_tracks];
				AddTrackData(foot1, foot2, foot15, foot16, mul_vec, a2_foot, vec_foot, vec_PoQ, vert_foot, opan_ang, tofd_Q[N_glob_tracks] , mdf_AoZ[N_glob_tracks], N_glob_tracks, ToF[N_glob_tracks], FlightPath[N_glob_tracks]);
				if(DoAlignment && PoQ[N_glob_tracks]>0 && mul_f32==1 &&  mul_f30==1 && mul_f33==1 && mul_tofd1==1 && mul_foot<250 && tofd_Q[N_glob_tracks]>5.6 && tofd_Q[N_glob_tracks]<6.4 && vec_PoQ.Z()>4.22 && vec_PoQ.Z()<4.44 && N_glob_tracks==0) //f33
				//if (DoAlignment && PoQ[N_glob_tracks]>0 && mul_f32==1 &&  mul_f30==1 && mul_f33==1 && mul_tofd==4 && mul_foot<200 && tofd_Q[N_glob_tracks]>5.5 && tofd_Q[N_glob_tracks]<6.5 && vec_PoQ.Z()>3.79 && vec_PoQ.Z()<3.85 && N_glob_tracks>0 && N_glob_tracks==1) //f33
				{
					cout<<"Ref"<<endl;
					det_points align_data;
					//align_data.mw.SetXYZ(vertex_mwpc.X(),vertex_mwpc.Y(), 0.);
					align_data.f1.SetXYZ(f1_point_i.X() ,0., f1_point_i.Z());
					align_data.f2.SetXYZ(0. ,f2_point_i.Y(), f2_point_i.Z());
					align_data.f15.SetXYZ(f15_point_i.X(), 0., f15_point_i.Z());
					align_data.f16.SetXYZ(0., f16_point_i.Y(), f16_point_i.Z());
					align_data.f30.SetXYZ(f30_point_i.X(), f30_point_i.Y(), f30_point_i.Z());
					align_data.f32.SetXYZ(f32_point_i.X(),f32_point_i.Y(), f32_point_i.Z());
					align_data.flast.SetXYZ(flast_point_i.X(),flast_point_i.Y(), flast_point_i.Z());
					det_points_vec.push_back(align_data);
				}
				N_glob_tracks++;
			}
		}
	}
	nglt.push_back(N_glob_tracks);	
	return;
}

void R3BTrackingS522::FinishEvent()
{
	for (auto& DataItem : fDataItems)
	{
		DataItem->Clear();
	}
	ver_Z.clear();
	ver_Y.clear();
	ver_X.clear();
	fTrackItems->Clear();
	f1_hits.clear();
	f2_hits.clear();
	f15_hits.clear();
	f16_hits.clear();
	tracks_in.clear();
	tracks_out.clear();
	tof.clear();
	poqm.clear();
	flight_p.clear();
	tofdq.clear();	
	nint.clear();
	nogt.clear();
	nglt.clear();	
}

void R3BTrackingS522::FinishTask()
{
	//LOG(info) << "Processed " << fNEvents << " events\n\n";
	if (DoAlignment) Alignment();
	//tree_out.Write();
	vtx.clear();
	vty.clear();
	vtz.clear();
	opan.clear();
	distm.clear();
	tpr.clear();
	tpl.clear();
	ppr.clear();
	ppl.clear();
	ecr.clear();
	ecl.clear();
	//cout<<"WRITE"<<endl;
	//cout << "\n\n------- Statisitcs summary --------- ";
}

bool R3BTrackingS522::IsGoodFootHit(R3BFootHitData* fhit)
{
	if(fhit->GetEnergy()>FootEnergyMin && fhit->GetEnergy()<FootEnergyMax)
		return true;
	else 
		return false;
}

bool R3BTrackingS522::IsGoodFiberHit(R3BFiberMAPMTHitData* fhit)
{
	if((fhit->GetEloss() > FiberEnergyMin) && (fhit->GetEloss() < FiberEnergyMax) && 
			(fhit->GetTime() < 20000 && fhit->GetTime()>(-20000) )
	  )
		return true;
	else 
		return false;
}

bool R3BTrackingS522::SortFootData()
{
	f1_hits.clear(); f2_hits.clear(); f15_hits.clear(); f16_hits.clear();
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
			case 15 ://Last FOOT15 (X) 
				f15_hits.push_back(f);
				break;
			case 16 ://3rd FOOT16 (Y)
				f16_hits.push_back(f);
				break;
		}
	}
	if(f1_hits.empty() || f2_hits.empty() || f15_hits.empty() || f16_hits.empty())     
		return false; 
	mul_f1  = f1_hits.size(); mul_f2  = f2_hits.size();
	mul_f15 = f15_hits.size(); mul_f16 = f16_hits.size();
	return true;
}

bool R3BTrackingS522::MakeIncomingTracks()
{
	if(fDataItems[MWPC0_HITDATA]->GetEntriesFast()==0 || fDataItems[MWPC1_HITDATA]->GetEntriesFast()==0) 
		return false;
	tracks_in.clear();
	N_in_tracks=0;
	if(!SortFootData()) return false;//at least 1 hit in every FOOT
	//TVector3 vertex_mwpc, vertex_foot;//projection to the center of the target (0,0,0)
	double tx_in = -999;
	double ty_in = -999; 
	double dx_vertex = -999;
	double dy_vertex = -999;
	double dx_angle = -999;
	double dy_angle = -999;	
	Track tr;
	//Get MWPC hits, for now only first hit
	auto m0_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC0_HITDATA]->At(0));
	auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
	m0_point.SetXYZ(m0_hit->GetX()*0.1, m0_hit->GetY()*0.1, 0.);//cm
	m1_point.SetXYZ(m1_hit->GetX()*0.1, m1_hit->GetY()*0.1, 0.);//cm
	m0_point_i1=m0_point;
	m1_point_i1=m1_point;	

	TransformPoint(m0_point, &m0_angles, &m0_position);//lab
	TransformPoint(m1_point, &m1_angles, &m1_position);//lab

	//Fill output tree variables
	m0_X = m0_point.X();  m0_Y = m0_point.Y();  m0_Z = m0_point.Z();
	m1_X = m1_point.X();  m1_Y = m1_point.Y();  m1_Z = m1_point.Z();
	//cout<<"m0_x--> "<<m0_X<<endl;
	//------- Project mwpc track to the center of the target
	tx_in = (m0_point.X() - m1_point.X())/(m0_point.Z() - m1_point.Z());
	ty_in = (m0_point.Y() - m1_point.Y())/(m0_point.Z() - m1_point.Z()); 
	vertex_mwpc.SetX((m1_point.X() - tx_in * m1_point.Z())/*-1.111*/);
	vertex_mwpc.SetY((m1_point.Y() - ty_in * m1_point.Z())/*-0.91*/);
	vertex_mwpc.SetZ(0);

	//----- Make track candidates in FOOT and project them to the center of the target
	for (auto & f2 : f2_hits){
		auto f2_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f2));
		f2_point.SetXYZ(0, f2_hit->GetPosLab().Y() * 0.1, 0.); //cm
		f2_point_i=f2_point;
		TransformPoint(f2_point,  &f2_angles,  &f2_position);

		for (auto & f1 : f1_hits){
			auto f1_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f1));
			f1_point.SetXYZ(f1_hit->GetPosLab().X() * 0.1 , 0, 0.); //cm
			f1_point_i=f1_point;
			TransformPoint(f1_point,  &f1_angles,  &f1_position);       
			for (auto & f16 : f16_hits){
				auto f16_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f16));
				f16_point.SetXYZ(0, f16_hit->GetPosLab().Y() * 0.1, 0); //cm
				f16_point_i=f16_point;
				TransformPoint(f16_point, &f16_angles, &f16_position);

				for (auto & f15 : f15_hits){
					auto f15_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f15));
					f15_point.SetXYZ(f15_hit->GetPosLab().X() * 0.1, 0, 0); //cm
					f15_point_i=f15_point;
					TransformPoint(f15_point, &f15_angles, &f15_position);

					//project foot track to the center of the target
					tr.tx0 = (f15_point.X() - f1_point.X())/(f15_point.Z() - f1_point.Z());
					tr.ty0 = (f16_point.Y() - f2_point.Y())/(f16_point.Z() - f2_point.Z());
					vertex_foot.SetX(f1_point.X() - tr.tx0 * f1_point.Z());
					vertex_foot.SetY(f2_point.Y() - tr.ty0 * f2_point.Z());
					vertex_foot.SetZ(0);
					//Condition on the vertex matching
					dx_vertex = vertex_foot.X() - vertex_mwpc.X();
					dy_vertex = vertex_foot.Y() - vertex_mwpc.Y();
					dx_angle = tr.tx0 - tx_in;
					dy_angle = tr.ty0 - ty_in;
					//s522
					if((dx_vertex>-0.6 && dx_vertex<0.6) && (dy_vertex>-0.7 && dy_vertex<0.7)){
						if(abs(dx_angle)<0.02 && abs(dy_angle)<0.04){
							//s509
							tr.f1_x   = f1_point.X();
							tr.f1_z   = f1_point.Z();
							tr.f2_y   = f2_point.Y();
							tr.f2_z   = f2_point.Z();
							tr.f15_x  = f15_point.X();
							tr.f15_z  = f15_point.Z();
							tr.f16_y  = f16_point.Y();
							tr.f16_z  = f16_point.Z();

							tr.f1_Q   = f1_hit->GetEnergy();
							tr.f2_Q   = f2_hit->GetEnergy();
							tr.f15_Q   = f15_hit->GetEnergy();
							tr.f16_Q   = f16_hit->GetEnergy();

							tr.mw_ax  = tx_in;
							tr.mw_ay  = ty_in;
							tracks_in.push_back(tr);
							N_in_tracks++;
						}
					}
				}
			}
		}
	}
	if((vertex_mwpc.X()>50 || vertex_mwpc.X()<-50) || (vertex_mwpc.Y()>50 || vertex_mwpc.Y()<-50)) return false;
	else return true;
}

bool R3BTrackingS522::MakeOutgoingTracks()
{
	if(fDataItems[DET_FI32]->GetEntriesFast() == 0 || fDataItems[DET_FI30]->GetEntriesFast() == 0 || 
			(fDataItems[DET_FI33]->GetEntriesFast() == 0 && fDataItems[DET_FI31]->GetEntriesFast()==0) ) 
		return false;
	tracks_out.clear();
	//tracks_out2.clear();
	Track tr;
	N_out_tracks=0;
	double angle_out, f30_slope, f30_offset, track_slope, track_offset;
	std::vector<double> f32x;
	std::vector<double> f30y;
	std::vector<double> flastx;
	std::vector<double> flast2x;
	std::vector<double> f32e;
	std::vector<double> f30e;
	std::vector<double> flaste;
	std::vector<double> flast2e;
	TVector3 f30_edge[2];//to extract z and x in f30
	for (auto i=0; i<fDataItems[DET_FI32]->GetEntriesFast(); ++i)
	{
		auto f32 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI32]->At(i));
		if(!IsGoodFiberHit(f32)) continue;
		double fitime = fTimeStitch->GetTime(f32->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -1200 && fitime < 15000))
			//if((fitime > 13460 && fitime < 13490))
		{	f32x.push_back(f32->GetX());
			f32e.push_back(f32->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI30]->GetEntriesFast(); ++i)
	{
		auto f30 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI30]->At(i));
		if(!IsGoodFiberHit(f30)) continue;
		double fitime = fTimeStitch->GetTime(f30->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -1200 && fitime < 15000))
			//if((fitime > 13460 && fitime < 13490))
		{
			f30y.push_back(f30->GetY());
			f30e.push_back(f30->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI33]->GetEntriesFast(); ++i)
	{
		auto f33 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI33]->At(i));
		if(!IsGoodFiberHit(f33)) continue;

		double fitime = fTimeStitch->GetTime(f33->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -1200 && fitime < 15000))
			//if((fitime > 13460 && fitime < 13490))
		{
			flastx.push_back(f33->GetX());
			flaste.push_back(f33->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI31]->GetEntriesFast(); ++i)
	{
		auto f31 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI31]->At(i));
		if(!IsGoodFiberHit(f31)) continue;
		double fitime = fTimeStitch->GetTime(f31->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -1200 && fitime < 15000))
			//if((fitime > 13460 && fitime < 13490))
		{
			flast2x.push_back(f31->GetX());
			flast2e.push_back(f31->GetEloss());
		}
	}
	if(f32x.size() < 1 || f30y.size() < 1 || (flastx.size() < 1 && flast2x.size() < 1))
		return false;
	Double_t cluster2[f32x.size()][f32x.size()][2];
	Bool_t set2[f32x.size()];
	Int_t num_hit2[f32x.size()];
	Int_t hit_clust_id2[f32x.size()];
	Int_t num_clust2 = 0;
	for(Int_t i = 0; i < f32x.size(); i++)
	{
		set2[i] = false;
		num_hit2[i] = 0;
		hit_clust_id2[i] = 0;
		for(Int_t j = 0; j < f32x.size(); j++)
		{
			cluster2[i][j][0] = 0./0.;
			cluster2[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < f32x.size(); i++)
	{
		if(!set2[i])
		{
			cluster2[num_clust2][num_hit2[num_clust2]][0] = f32x[i];
			cluster2[num_clust2][num_hit2[num_clust2]][1] = f32e[i];
			num_hit2[num_clust2]++;
			set2[i] = true;
			hit_clust_id2[i] = num_clust2;
			num_clust2++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit2[hit_clust_id2[i]]; k++)
		{
			for(Int_t j = 0; j < f32x.size(); j++)
			{
				if(set2[j])
					continue;

				if(fabs(f32x[j]-cluster2[hit_clust_id2[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id2[i];
					cluster2[id][num_hit2[id]][0] = f32x[j];
					cluster2[id][num_hit2[id]][1] = f32e[j];
					num_hit2[id]++;
					set2[j] = true;
				}
			}
		}
	}
	Double_t cluster0[f30y.size()][f30y.size()][2];
	Bool_t set0[f30y.size()];
	Int_t num_hit0[f30y.size()];
	Int_t hit_clust_id0[f30y.size()];
	Int_t num_clust0 = 0;
	for(Int_t i = 0; i < f30y.size(); i++)
	{
		set0[i] = false;
		num_hit0[i] = 0;
		hit_clust_id0[i] = 0;
		for(Int_t j = 0; j < f30y.size(); j++)
		{
			cluster0[i][j][0] = 0./0.;
			cluster0[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < f30y.size(); i++)
	{
		if(!set0[i])
		{
			cluster0[num_clust0][num_hit0[num_clust0]][0] = f30y[i];
			cluster0[num_clust0][num_hit0[num_clust0]][1] = f30e[i];
			num_hit0[num_clust0]++;
			set0[i] = true;
			hit_clust_id0[i] = num_clust0;
			num_clust0++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit0[hit_clust_id0[i]]; k++)
		{
			for(Int_t j = 0; j < f30y.size(); j++)
			{
				if(set0[j])
					continue;

				if(fabs(f30y[j]-cluster0[hit_clust_id0[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id0[i];
					cluster0[id][num_hit0[id]][0] = f30y[j];
					cluster0[id][num_hit0[id]][1] = f30e[j];
					num_hit0[id]++;
					set0[j] = true;
				}
			}
		}
	}
	Double_t cluster3[flastx.size()][flastx.size()][2];
	Bool_t set3[flastx.size()];
	Int_t num_hit3[flastx.size()];
	Int_t hit_clust_id3[flastx.size()];
	Int_t num_clust3 = 0;
	for(Int_t i = 0; i < flastx.size(); i++)
	{
		set3[i] = false;
		num_hit3[i] = 0;
		hit_clust_id3[i] = 0;
		for(Int_t j = 0; j < flastx.size(); j++)
		{
			cluster3[i][j][0] = 0./0.;
			cluster3[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < flastx.size(); i++)
	{
		if(!set3[i])
		{
			cluster3[num_clust3][num_hit3[num_clust3]][0] = flastx[i];
			cluster3[num_clust3][num_hit3[num_clust3]][1] = flaste[i];
			num_hit3[num_clust3]++;
			set3[i] = true;
			hit_clust_id3[i] = num_clust3;
			num_clust3++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit3[hit_clust_id3[i]]; k++)
		{
			for(Int_t j = 0; j < flastx.size(); j++)
			{
				if(set3[j])
					continue;

				if(fabs(flastx[j]-cluster3[hit_clust_id3[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id3[i];
					cluster3[id][num_hit3[id]][0] = flastx[j];
					cluster3[id][num_hit3[id]][1] = flaste[j];
					num_hit3[id]++;
					set3[j] = true;
				}
			}
		}
	}
	Double_t cluster1[flast2x.size()][flast2x.size()][2];
	Bool_t set1[flast2x.size()];
	Int_t num_hit1[flast2x.size()];
	Int_t hit_clust_id1[flast2x.size()];
	Int_t num_clust1 = 0;
	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		set1[i] = false;
		num_hit1[i] = 0;
		hit_clust_id1[i] = 0;
		for(Int_t j = 0; j < flast2x.size(); j++)
		{
			cluster1[i][j][0] = 0./0.;
			cluster1[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		if(!set1[i])
		{
			cluster1[num_clust1][num_hit1[num_clust1]][0] = flast2x[i];
			cluster1[num_clust1][num_hit1[num_clust1]][1] = flast2e[i];
			num_hit1[num_clust1]++;
			set1[i] = true;
			hit_clust_id1[i] = num_clust1;
			num_clust1++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit1[hit_clust_id1[i]]; k++)
		{
			for(Int_t j = 0; j < flast2x.size(); j++)
			{
				if(set1[j])
					continue;

				if(fabs(flast2x[j]-cluster1[hit_clust_id1[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id1[i];
					cluster1[id][num_hit1[id]][0] = flast2x[j];
					cluster1[id][num_hit1[id]][1] = flast2e[j];
					num_hit1[id]++;
					set1[j] = true;
				}
			}
		}
	}
	Double_t fi32xcl[num_clust2];
	Double_t fi30ycl[num_clust0];
	Double_t filastxcl[num_clust3];
	Double_t filast2xcl[num_clust1];

	for(Int_t i = 0; i < num_clust2; i++)
	{
		fi32xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit2[i]; j++)
		{
			fi32xcl[i] += cluster2[i][j][0]*cluster2[i][j][1];
			sum_ener +=cluster2[i][j][1];
		}
		if(sum_ener > 0.)
			fi32xcl[i] = (double)fi32xcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust0; i++)
	{
		fi30ycl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit0[i]; j++)
		{
			fi30ycl[i] += cluster0[i][j][0]*cluster0[i][j][1];
			sum_ener +=cluster0[i][j][1];
		}
		if(sum_ener > 0.)
			fi30ycl[i] = (double)fi30ycl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust3; i++)
	{
		filastxcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit3[i]; j++)
		{
			filastxcl[i] += cluster3[i][j][0]*cluster3[i][j][1];
			sum_ener +=cluster3[i][j][1];
		}
		if(sum_ener > 0.)
			filastxcl[i] = (double)filastxcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust1; i++)
	{
		filast2xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit1[i]; j++)
		{
			filast2xcl[i] += cluster1[i][j][0]*cluster1[i][j][1];
			sum_ener +=cluster1[i][j][1];
		}
		if(sum_ener > 0.)
			filast2xcl[i] = (double)filast2xcl[i]/sum_ener;
	}
	for(auto i = 0; i < num_clust2; i ++)
	{
		f32_point.SetXYZ(fi32xcl[i], 0, 0); //cm
		f32_point_i=f32_point;
		TransformPoint(f32_point, &f32_angles, &f32_position);
		tr.f32_x = f32_point.X();
		tr.f32_z = f32_point.Z();

		for (auto j=0; j<num_clust0; ++j)
		{
			f30_point.SetXYZ(0,fi30ycl[j],0); //cm
			f30_point_i=f30_point;
			TransformPoint(f30_point, &f30_angles, &f30_position);
			tr.f30_y = f30_point.Y();
			//make combination with every hit in fibers 33 and 31

			for (auto k = 0; k<num_clust3; ++k)
			{
				flast_point.SetXYZ(filastxcl[k], 0, 0); //cm
				flast_point_i=flast_point;
				TransformPoint(flast_point, &f33_angles, &f33_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//cout<<"a_33-> "<<angle_out<<endl;
				if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &f30_angles, &f30_position);
				TransformPoint(f30_edge[1], &f30_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				N_out_tracks++;
				tracks_out.push_back(tr);
				//	if(N_out_tracks==N_glob_tracks_max/2) return true;
			}

			for (auto k = 0; k<num_clust1; ++k)
			{
				flast_point.SetXYZ(filast2xcl[k], 0, 0); //cm
				flast_point_i=flast_point;
				TransformPoint(flast_point, &f31_angles, &f31_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//cout<<"a_31-> "<<angle_out<<endl;
				if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &f30_angles, &f30_position);
				TransformPoint(f30_edge[1], &f30_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				N_out_tracks++;
				tracks_out.push_back(tr);
				//	if(N_out_tracks==N_glob_tracks_max/2) return true;
			}
		}
	}
	if(tracks_out.empty()) return false;
	return true;
	nogt.push_back(N_out_tracks);	
}





/*

   bool R3BTrackingS522::MakeOutgoingTracks()
   {
   if(fDataItems[DET_FI32]->GetEntriesFast() == 0 || fDataItems[DET_FI30]->GetEntriesFast() == 0 || 
   (fDataItems[DET_FI33]->GetEntriesFast() == 0 && fDataItems[DET_FI31]->GetEntriesFast()==0) ) 
   return false;
   tracks_out.clear();
//id.clear();
Track tr;
N_out_tracks=0;
double angle_out = -999;
double f30_slope = -999; 
double f30_offset = -999;
double track_slope = -999;
double track_offset = -999;
TVector3 f30_edge[2];//to extract z and x in f30

//if (tofd_hit->GetToF < 20 && tofd_hit->GetToF>35 ) continue;
// only hits from first plane

for (auto i=0; i<fDataItems[DET_FI32]->GetEntriesFast(); ++i)
{
auto f32 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI32]->At(i));
//tofd_hit->GetEloss();
if(!IsGoodFiberHit(f32)) continue;
// TOFF if(f32->GetTime()<7990 || f32->GetTime()>8040) continue;  
f32_point.SetXYZ(f32->GetX(), 0, 0); //cm
f32_point_i=f32_point;
TransformPoint(f32_point, &f32_angles, &f32_position);
tr.f32_x = f32_point.X();
tr.f32_z = f32_point.Z();

for (auto j=0; j<fDataItems[DET_FI30]->GetEntriesFast(); ++j)
{
auto f30 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI30]->At(j));
bb++;
if(!IsGoodFiberHit(f30)) continue;
b1++;
f30_point.SetXYZ(0,f30->GetY(), 0); //cm
f30_point_i=f30_point;
TransformPoint(f30_point, &f30_angles, &f30_position);
tr.f30_y = f30_point.Y();
tr.f30_z = f30_point.Z();

//if(fabs(f32->GetTime_ns() - f30->GetTime_ns())>30) continue;
// TOFF if(f30->GetTime()<7980 || f30->GetTime()>8040) continue;
//make combination with every hit in fibers 33 and 31
for (auto k = 0; k<fDataItems[DET_FI33]->GetEntriesFast(); ++k)
{
auto f33 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI33]->At(k));
if(!IsGoodFiberHit(f33)) continue; //Messel side
//  TOFF if(f33->GetTime()<7980 || f33->GetTime()>8040) continue;
//if((fabs(f32->GetTime_ns() - f33->GetTime_ns())>30) && (fabs(f32->GetTime_ns() - f33->GetTime_ns())<-5)) continue;
flast_point.SetXYZ(f33->GetX(), 0, 0); //cm
flast_point_i1=flast_point;
TransformPoint(flast_point, &f33_angles, &f33_position);
tr.last_x = flast_point.X();
tr.last_z = flast_point.Z();
angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
if(angle_out>(-10.) || angle_out<(-18.)) continue;
// We need to extrapolate Z position in f30 because it was used for Y measurement
// Define two (X,Z) points on the f30 plane:
//Now track every combination of upstream and downstream tracks 
f30_edge[0].SetXYZ(-1, 0, 0);
f30_edge[1].SetXYZ(1, 0, 0);
TransformPoint(f30_edge[0], &f30_angles, &f30_position);
TransformPoint(f30_edge[1], &f30_angles, &f30_position);
// Parameterize f30 plane
f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
track_offset = (tr.last_x - track_slope * tr.last_z);
// Extrapolate final X and Z position in f30
tr.f30_z = (track_offset - f30_offset) / (f30_slope - track_slope);// extrapolated
tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
N_out_tracks++;
tracks_out.push_back(tr);
id.push_back(33);
//if(N_out_tracks==N_glob_tracks/2) return true;
}

//make combination with every hit in fibers 33 and 31
for (auto k = 0; k<fDataItems[DET_FI31]->GetEntriesFast(); ++k)
{
	auto f31 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI31]->At(k));
	if(!IsGoodFiberHit(f31)) continue; //Messel side
	as1++;
	// TOFFif(f31->GetTime()<7990 && f31->GetTime()>8040) continue;
	//if((f32->GetTime_ns() - f31->GetTime_ns())>20 ||(f32->GetTime_ns() - f31->GetTime_ns())<(-10) ) continue;
	as++;
	//cout<<"Ciao"<<endl;
	flast_point.SetXYZ(f31->GetX(), 0, 0); //cm
	flast_point_i=flast_point;
	TransformPoint(flast_point, &f31_angles, &f31_position);
	tr.last_x = flast_point.X();
	tr.last_z = flast_point.Z();
	angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();

	//cout<<"Ciao1"<<endl;
	if(angle_out>(-10.) || angle_out<(-18.)) continue;
	// We need to extrapolate Z position in f30 because it was used for Y measurement
	// Define two (X,Z) points on the f30 plane:
	//Now track every combination of upstream and downstream tracks 
	f30_edge[0].SetXYZ(-1, 0, 0);
	f30_edge[1].SetXYZ(1, 0, 0);
	TransformPoint(f30_edge[0], &f30_angles, &f30_position);
	TransformPoint(f30_edge[1], &f30_angles, &f30_position);
	// Parameterize f30 plane
	f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
	f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
	track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
	track_offset = (tr.last_x - track_slope * tr.last_z);
	// Extrapolate final X and Z position in f30
	tr.f30_z = (track_offset - f30_offset) / (f30_slope - track_slope);// extrapolated
	tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated

	N_out_tracks++;
	tracks_out.push_back(tr);
	//cout<<"3"<<endl;
	//id.push_back(31);
	//cout<<"4"<<endl;
	//if(N_out_tracks==N_glob_tracks/2) return true;
}
}
}
if(tracks_out.empty()) return false;
//if(id.empty()) return false;
return true;
nogt.push_back(N_out_tracks);	
}
	*/	
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

	/*	cout<<"pointx "<<point1.X()<<endl;
		cout<<"pointy "<<point1.Y()<<endl;
		cout<<"pointz "<<point1.Z()<<endl;

		cout<<"rotx "<<rot1.X()<<endl;
		cout<<"roty "<<rot1.Y()<<endl;
		cout<<"rotz "<<rot1.Z()<<endl;

		cout<<"transx "<<trans1.X()<<endl;
		cout<<"transy "<<trans1.Y()<<endl;
		cout<<"transz "<<trans1.Z()<<endl;
		*/
	r1.SetToIdentity();
	// First Euler rotation around Y axis
	r1.RotateY(rot1.Y());
	// get local X axis after first rotation
	v3_localX1.SetMagThetaPhi(1, r1.ThetaX(), r1.PhiX());
	//cout<<"R1tx "<<r1.ThetaX()<<endl;
	//cout<<"R1px "<<r1.PhiX()<<endl;
	//cout<<"rotz "<<rot.Z()<<endl;


	// Second Euler rotation around local X axis
	r1.Rotate(rot1.X(), v3_localX1);
	// get local Z axis after second rotation	
	v3_localZ1.SetMagThetaPhi(1, r1.ThetaZ(), r1.PhiZ());
	//cout<<"R1tZ "<<r1.ThetaZ()<<endl;
	//cout<<"R1pZ "<<r1.PhiZ()<<endl;

	// final rotation around local Z axis
	r1.Rotate(rot1.Z(), v3_localZ1);
	point1.Transform(r1);
	point1 += (trans1);
	//cout<<"pointx_a  "<<point1.X()<<endl;
	//cout<<"pointy_a  "<<point1.Y()<<endl;
	//cout<<"pointz_a  "<<point1.Z()<<endl;


	return;
}

R3BTrackS522* R3BTrackingS522::AddTrackData(TVector3 F1, TVector3 F2, TVector3 F15, TVector3 F16,TVector3 a1, TVector3 a2, TVector3 mw, TVector3 poq, TVector3 vrt, double opag,  double charge, double aoz, int Mult, double ToFm, double FlightPathm)
{
	// Filling output track info
	TClonesArray& clref = *fTrackItems;
	Int_t size = clref.GetEntriesFast();
	return new (clref[size]) R3BTrackS522(F1, F2, F15, F16, a1.X(), a1.Y(), a1.Z(), a2.X(), a2.Y(), a2.Z(), mw.X(), mw.Y(), mw.Z(), poq.X(), poq.Y(), poq.Z(),vrt.X(), vrt.Y(), vrt.Z(), opag,  charge, aoz, Mult, ToFm, FlightPathm, 0., 0., 0.);
}

// Setup energy cuts in foot and fibers 
void R3BTrackingS522::SetFootEnergyMinMax(double min, double max)
{
	FootEnergyMin = min;
	FootEnergyMax = max;
	return;
}
void R3BTrackingS522::SetFiberEnergyMinMax(double min, double max)
{
	FiberEnergyMin = min;
	FiberEnergyMax = max;
	return;
}
void R3BTrackingS522::Alignment()
{
	int run_flag = 0;
	std::cout << "\n\n----------------------------------------------------";
	std::cout << "\n Ready for detector alignment using the following data set:";
	std::cout << "\n\tNumber of reference tracks: " << det_points_vec.size();
	//std::cout << "\n\tFRS Brho: min = " << FrsBrhoMin << ", max = " << FrsBrhoMax;
	//std::cout << "\n\tFRS PoQ: min = " << FrsBrhoMin / 3.3356 << ", max = " << FrsBrhoMax / 3.3356;
	std::cout << "\n\tP/Z reference: " << reference_PoQ << "GeV/c";
	std::cout << "\n\n Do you want to continue? (0 = no, 1 = yes):  ";
	std::cin >> run_flag;
	if (run_flag == 0)
		return;

	// Now define minimizer
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	const Int_t NVarsFunctor = 3; // number of the alignment offsets: (2 angles + 3 offsets) x 4 detectors
	const double* xs={0};
	double step[NVarsFunctor]={0};
	double offset[NVarsFunctor]={0};
	double min_offset[NVarsFunctor]={0};
	double max_offset[NVarsFunctor]={0};
	TH1F* hist[NVarsFunctor]={0};
	TH2D* hist1[NVarsFunctor]={0};
	char hname[100];

	// For every detector: 3 Euler angles and 3 shifts
	for (auto d = 0; d < 3; ++d)
	{
		/*	for (auto a = 0; a < 0; ++a) // angle shifts in rad
			{
			min_offset[d * 4 + a] = -0.1;
			max_offset[d * 4 + a] = 0.1;
			step[d * 4 + a] = 0.001;
			}*/
		for (auto o = 0; o < 2; ++o) // position shifts in cm
		{

			if(d==0){
				if(o==0){
					min_offset[d * 1 + 0 + o] = -30;
					max_offset[d * 1 + 0 + o] = 30;
					step[d * 1 + 0 + o] = 0.01; //Valerii
					//step[d * 5 + 2 + o] = 0.01;
				}

				/*if(o==1){
				  min_offset[d * 2 + 0 + o] = -60;
				  max_offset[d * 2 + 0 + o] = 60;
				  step[d * 2 + 0 + o] = 0.01; //Valerii
				//step[d * 5 + 2 + o] = 0.01;
				}

				if(o==2){
				min_offset[d * 3 + 0 + o] = -10;
				max_offset[d * 3 + 0 + o] = 10;
				step[d * 3 + 0 + o] = 0.01; //Valerii
				//step[d * 5 + 2 + o] = 0.01;
				}*/
			}

			/*	if(d==0){

				min_offset[d * 3 + 0 + o] = -2;
				max_offset[d * 3 + 0 + o] = 2;
				step[d * 3 + 0 + o] = 0.001; //Valerii
			//step[d * 5 + 2 + o] = 0.01;
			}
			*/
			if(d==1){
				if(o==0){
					min_offset[d * 1 + 0 + o] = -30;
					max_offset[d * 1 + 0 + o] =30;
					step[d * 1 + 0 + o] = 0.01; //Valerii
					//step[d * 5 + 2 + o] = 0.01;
				}
				/*if(o==1){
				  min_offset[d * 2 + 0 + o] = -60;
				  max_offset[d * 2 + 0 + o] = 60;
				  step[d * 2 + 0 + o] = 0.01; //Valerii

				  }

				  if(o==2){
				  min_offset[d * 3 + 0 + o] = -10;
				  max_offset[d * 3 + 0 + o] = 10;
				  step[d * 3 + 0 + o] = 0.01; //Valerii

				  }*/
			}

			if(d==2){

				if(o==0){
					min_offset[d * 1 + 0 + o] = -30;
					max_offset[d * 1 + 0 + o] = 30;
					step[d * 1 + 0 + o] = 0.01; //Valerii
					//step[d * 5 + 2 + o] = 0.01;
				}

				/*if(o==1){
				  min_offset[d * 2 + 0 + o] = -60;
				  max_offset[d * 2 + 0 + o] = 60;
				  step[d * 2 + 0 + o] = 0.01; //Valerii

				  }

				  if(o==2){
				  min_offset[d * 3 + 0 + o] = -10;
				  max_offset[d * 3 + 0 + o] = 10;
				  step[d * 3 + 0 + o] = 0.01; //Valerii

				  }*/
			}

			if(d==3){

				if(o==0){
					min_offset[d * 1 + 0 + o] = -12;
					max_offset[d * 1 + 0 + o] = 12;
					step[d * 1 + 0 + o] = 0.01; //Valerii
					//step[d * 5 + 2 + o] = 0.01;
				}
				/*if(o==1){
				  min_offset[d * 2 + 0 + o] = -15;
				  max_offset[d * 2 + 0 + o] = 15;
				  step[d * 2 + 0 + o] = 0.01; //Valerii

				  }

				  if(o==2){
				  min_offset[d * 3 + 0 + o] = -10;
				  max_offset[d * 3 + 0 + o] = 10;
				  step[d * 3 + 0 + o] = 0.01; //Valerii

				  }*/
			}

			if(d==4){
				if(o==0){
					min_offset[d * 1 + 0 + o] = -10;
					max_offset[d * 1 + 0 + o] = 10;
					step[d * 1 + 0 + o] = 0.01; //Valerii
					//step[d * 5 + 2 + o] = 0.1;
				}
				/*if(o==1){
				  min_offset[d * 2 + 0 + o] = -15;
				  max_offset[d * 2 + 0 + o] = 15;
				  step[d * 2 + 0 + o] = 0.01; //Valerii

				  }*/
			}
			/*	if(d==5){
				if(o==0){
				min_offset[d * 1 + 0 + o] = -8;
				max_offset[d * 1 + 0 + o] = 8;
				step[d * 1 + 0 + o] = 0.01; //Valerii
			//step[d * 5 + 2 + o] = 0.1;
			}
			if(o==1){
			min_offset[d * 2 + 0 + o] = -8;
			max_offset[d * 2 + 0 + o] = 8;
			step[d * 2 + 0 + o] = 0.01; //Valerii

			}
			}
			if(d==6){
			if(o==0){
			min_offset[d * 1 + 0 + o] = -8;
			max_offset[d * 1 + 0 + o] = 8;
			step[d * 1 + 0 + o] = 0.01; //Valerii
			//step[d * 5 + 2 + o] = 0.1;
			}
			if(o==1){
			min_offset[d * 2 + 0 + o] = -8;
			max_offset[d * 2 + 0 + o] = 8;
			step[d * 2 + 0 + o] = 0.01; //Valerii

			}
			}*/
		}
	}
	//}
for (Int_t i = 0; i < NVarsFunctor; i++)
{
	sprintf(hname, "par%d", i);
	hist[i] = new TH1F(hname, hname, 6000, min_offset[i], max_offset[i]);
	hist1[i] = new TH2D(hname, hname, 1000, -30,30, 1000, 0.001, 0.08);
}
std::cout << "\n\n-- Perfroming minimization for the detector alignment. Please wait... \n";

// Setting up Minimizer function parameters
double precision = 1e-10; // 0 - default precision will be automaticalle determined
double tolerance = 0.001;
minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
minimizer->SetMaxIterations(1000);           // for GSL
minimizer->SetTolerance(tolerance);
minimizer->SetPrecision(precision);
minimizer->SetPrintLevel(2);
minimizer->SetStrategy(0); // to run faster
ROOT::Math::Functor f(&R3BTrackingS522::AlignmentErrorS522, NVarsFunctor);
minimizer->SetFunction(f);
minimizer->SetPrintLevel(0);
// Getting experimental data and doing minimization
Int_t i = 0;
Int_t bs = 0;
bool good_minimum = false;
while (bs <1000) // running several minimizations
{	
	minimizer->Clear();
	for (i = 0; i < NVarsFunctor; i++) // sampling +=50% from limits
	{


		if(i==0){offset[i]=3.924;}
		if(i==1){offset[i]=-0.4542;}
		//		   if(i==2){offset[i]=-1.669;}
		//		   if(i==3){offset[i]=2.542;}
		//	   if(i==4){offset[i]=-0.165392;}
		//	   if(i==5){offset[i]=1.47196;}
		//	   if(i==6){offset[i] = 1.88823;}
		//	   if(i==7){offset[i]=1.05481;}
		//if(i==8){offset[i] = 0.913645;}
		//	   if(i==9){offset[i]=-1.70050;}
		//	   if(i==10){offset[i]=0.66066;}
		//	   if(i==11){offset[i]=0.308029;}
		//	   if(i==12){offset[i]=0.739507;}

		if(i==0 || i==1){
			step[i]=0;
		}
		else {

			offset[i] = gRandom->Uniform(min_offset[i]*0.5 , max_offset[i]*0.5 );
			//cout<<"Off--> "<<offset[i]<<endl;
		}
		minimizer->SetVariable(i, Form("par%d", i), offset[i], step[i]);
		minimizer->SetVariableLimits(i, min_offset[i], max_offset[i]);	

	}
	minimizer->Minimize();
	xs = minimizer->X();
	//cout<<"Min status--> "<<minimizer->Status()<<endl;
	// if(minimizer->Status() !=0) continue; //valid minimum
	// Check if all paramters are "far" from limits
	for (i = 0; i < NVarsFunctor; i++)
	{
		//cout<<"/n cond1--> "<<fabs((xs[i] - min_offset[i]) / min_offset[i])<<" less then 0.1 "<<min_offset[i]<<" "<<xs[i]<<endl;

		cout<<"cond2--> "<</*fabs((max_offset[i] -xs[i]) / max_offset[i])<<" less then 0.1 "<<max_offset[i]<<" "<<*/xs[i]<<"--> "<<i<<endl;

		if (fabs((xs[i] - min_offset[i]) / min_offset[i]) < 0.1 ||
				fabs((max_offset[i] - xs[i]) / max_offset[i]) < 0.1)
		{
			break;
		}
	}
	if (i != NVarsFunctor)
	{
		std::cout << "\n\n -- Parameter is too close to the limit!! ";
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
		hist1[i]->Fill(xs[i],minimizer->MinValue());
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
c1->Divide(3, 1);
for (Int_t v = 0; v < NVarsFunctor; v++)
{
	c1->cd(v + 1);
	hist[v]->Draw();
}

TCanvas* c2 = new TCanvas("c2", "c2", 1000, 1000);
c2->Divide(3, 1);
for (Int_t v = 0; v < NVarsFunctor; v++)
{
	c2->cd(v + 1);
	hist1[v]->Draw();
}
return;
}


double R3BTrackingS522::AlignmentErrorS522(const double* par)
{




	//gMDFTrackerS522->f1_ang_offset.SetXYZ(par[0], par[1], par[2]);
	//gMDFTrackerS522->f1_pos_offset.SetXYZ(par[0], 0  , 0);
	//gMDFTrackerS522->f2_pos_offset.SetXYZ(0, par[1]  , 0);
	//gMDFTrackerS522->f15_pos_offset.SetXYZ(par[2], 0  , 0);
	//gMDFTrackerS522->f16_pos_offset.SetXYZ(0, par[3]  , 0);

	//gMDFTrackerS522->f1_ang_offset.SetXYZ(par[0], par[1], par[2]);
	//gMDFTrackerS522->m0_pos_offset.SetXYZ(par[2], par[3]  , 0);	

	//gMDFTrackerS522->f1_ang_offset.SetXYZ(par[0], par[1], par[2]);
	//gMDFTrackerS522->mw_pos_offset.SetXYZ(par[0], par[1]  , 0);

	//gMDFTrackerS522->f2_ang_offset.SetXYZ(par[4], par[5], par[6]);
	//gMDFTrackerS522->f2_pos_offset.SetXYZ(0, par[7], 0);

	//gMDFTrackerS522->f32_ang_offset.SetXYZ(par[8], par[9], par[10]);
	gMDFTrackerS522->f32_pos_offset.SetXYZ(par[0], 0 , 0);

	//gMDFTrackerS522->f30_ang_offset.SetXYZ(par[12], par[13], par[14]);
	gMDFTrackerS522->f30_pos_offset.SetXYZ( 0, par[1], 0);

	//gMDFTrackerS522->flast_ang_offset.SetXYZ(par[16], par[17], par[18]);
	gMDFTrackerS522->flast_pos_offset.SetXYZ(par[2], 0, 0);


	double mdf_input[8]; // data container for the MDF function

	double tx_in_i = -999;
	double ty_in_i = -999; 
	double v1 = 0;
	double v1_1 = 0;
	double v2 = 0;
	double v3 = 0;
	double v4 = 0;
	double v3_1 = 0;
	double v4_1 = 0;
	double v = 0;
	double v_tot = 0;
	int counter = 0;
	for (auto& d : (gMDFTrackerS522->det_points_vec))
	{
		//gMDFTrackerS522->vertex_mwpc = d.mw;
		gMDFTrackerS522->f1_point_i = d.f1;
		gMDFTrackerS522->f2_point_i = d.f2;
		gMDFTrackerS522->f15_point_i = d.f15;
		gMDFTrackerS522->f16_point_i = d.f16;
		gMDFTrackerS522->f32_point_i = d.f32;
		gMDFTrackerS522->f30_point_i = d.f30;
		gMDFTrackerS522->flast_point_i = d.flast;	

		TVector3 rotv(0,0,0);	
		//cout<<"MDF_MWx-> "<<gMDFTrackerS522->vertex_mwpc.X()<<endl;
		//cout<<"BT-> "<<gMDFTrackerS522->f1_point_i.X()<<endl;
		// This will transform "det_point" vectors into lab frame
		/*
		   cout<<"f1x_posoff_x-> "<<gMDFTrackerS522->f1_point_i.X()<<endl;
		   cout<<"f1y_posoff_z-> "<<gMDFTrackerS522->f1_point_i.Y()<<endl;
		   cout<<"f1z_posoff_y-> "<<gMDFTrackerS522->f1_point_i.Z()<<endl;
		   cout<<"f2x_posoff_z-> "<<gMDFTrackerS522->f2_point_i.X()<<endl;
		   cout<<"f2y_ang_x-> "<<gMDFTrackerS522->f2_point_i.Y()<<endl;
		   cout<<"f2z_ang_y-> "<<gMDFTrackerS522->f2_point_i.Z()<<endl;

*/
		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f1_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot1()/* + gMDFTrackerS522->mw_ang_offset*/,
				gMDFTrackerS522->GetPositionFoot1()/* +  gMDFTrackerS522->f1_pos_offset*/);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f2_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot2()/* + gMDFTrackerS522->mw_ang_offset*/,
				gMDFTrackerS522->GetPositionFoot2()/* +  gMDFTrackerS522->f2_pos_offset*/);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f15_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot15()/* + gMDFTrackerS522->mw_ang_offset*/,
				gMDFTrackerS522->GetPositionFoot15()/* +  gMDFTrackerS522->f15_pos_offset*/);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f16_point_i,
				gMDFTrackerS522->GetEulerAnglesFoot16()/* + gMDFTrackerS522->mw_ang_offset*/,
				gMDFTrackerS522->GetPositionFoot16()/* +  gMDFTrackerS522->f16_pos_offset*/);



		//gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->m0_point_i,
		//		gMDFTrackerS522->GetEulerAnglesMwpc0(),/*, + gMDFTrackerS522->mw_ang_offset,
		//		gMDFTrackerS522->GetPositionMwpc0() + gMDFTrackerS522->m0_pos_offset);


		/*gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f2_point_i,
		  gMDFTrackerS522->GetEulerAnglesFoot2() + gMDFTrackerS522->f2_ang_offset,
		  gMDFTrackerS522->GetPositionFoot2() + gMDFTrackerS522->f2_pos_offset);
		  */
		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f32_point_i,
				gMDFTrackerS522->GetEulerAnglesFiber32(), /*+ gMDFTrackerS522->f32_ang_offset,*/
				gMDFTrackerS522->GetPositionFiber32() + gMDFTrackerS522->f32_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f30_point_i,
				gMDFTrackerS522->GetEulerAnglesFiber30(), /*+ gMDFTrackerS522->f30_ang_offset,*/
				gMDFTrackerS522->GetPositionFiber30() + gMDFTrackerS522->f30_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->flast_point_i,
				gMDFTrackerS522->GetEulerAnglesFiber33(), /*+ gMDFTrackerS522->flast_ang_offset,*/
				gMDFTrackerS522->GetPositionFiber33() + gMDFTrackerS522->flast_pos_offset);

		/*cout<<"f1_posoff_x-> "<<gMDFTrackerS522->f1_pos_offset.X()<<endl;
		  cout<<"f1_posoff_z-> "<<gMDFTrackerS522->f1_pos_offset.Z()<<endl;
		  cout<<"f2_posoff_y-> "<<gMDFTrackerS522->f2_pos_offset.Y()<<endl;
		  cout<<"f2_posoff_z-> "<<gMDFTrackerS522->f2_pos_offset.Z()<<endl;
		  cout<<"f1_ang_x-> "<<gMDFTrackerS522->f1_ang_offset.X()<<endl;
		  cout<<"f1_ang_x-> "<<gMDFTrackerS522->f1_ang_offset.Y()<<endl;
		  */




		/*	mdf_input[0] = 0.;
			mdf_input[1] = 0.;
			mdf_input[2] = 0.;
			mdf_input[3] = gMDFTrackerS522->f32_point_i.X();
			mdf_input[4] = gMDFTrackerS522->f32_point_i.Z();
			mdf_input[5] = (gMDFTrackerS522->flast_point_i.X() - mdf_input[3]) / (gMDFTrackerS522->flast_point_i.Z() - mdf_input[4]);
			mdf_input[6] = (gMDFTrackerS522->f30_point_i.Y() - mdf_input[1]) / (gMDFTrackerS522->f30_point_i.Z() - mdf_input[2]);
			*/







		//IN beam alig	
		mdf_input[0] = gMDFTrackerS522->f2_point_i.Y();
		mdf_input[1] = gMDFTrackerS522->f2_point_i.Z();
		mdf_input[2] = gMDFTrackerS522->f1_point_i.X();
		mdf_input[3] = gMDFTrackerS522->f1_point_i.Z();
		mdf_input[4] = gMDFTrackerS522->f32_point_i.X();
		mdf_input[5] = gMDFTrackerS522->f32_point_i.Z();
		mdf_input[6] = (gMDFTrackerS522->flast_point_i.X() - mdf_input[4]) / (gMDFTrackerS522->flast_point_i.Z() - mdf_input[5]);
		mdf_input[7] = (gMDFTrackerS522->f30_point_i.Y() - mdf_input[0]) / (gMDFTrackerS522->f30_point_i.Z() - mdf_input[1]);

		TVector3 foot_angles_i;	
		foot_angles_i[0]=(gMDFTrackerS522->f15_point_i.X()-gMDFTrackerS522->f1_point_i.X())/(gMDFTrackerS522->f15_point_i.Z()-gMDFTrackerS522->f1_point_i.Z());
		foot_angles_i[1]=(gMDFTrackerS522->f16_point_i.Y()-gMDFTrackerS522->f2_point_i.Y())/(gMDFTrackerS522->f16_point_i.Z()-gMDFTrackerS522->f2_point_i.Z());
		foot_angles_i[2]=1.;
		double GladCurrent =2545.;
		double GladReferenceCurrent =2400.;
		foot_angles_i.SetMag(gMDFTrackerS522->GetReferencePoQ());
		//TVector3 vec_PoQ_a(foot_angles_i[0],foot_angles_i[1], foot_angles_i[2]);
		TVector3 vec_PoQ_a(gMDFTrackerS522->Get_MDF_TX0()->MDF(mdf_input),gMDFTrackerS522->Get_MDF_TY0()->MDF(mdf_input), 1);
		vec_PoQ_a.SetMag(gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input)* (GladCurrent / GladReferenceCurrent));

		//Difference between TXO and TYO and MWPC angles		

		v2 += pow((gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input) * (GladCurrent / GladReferenceCurrent) - gMDFTrackerS522->GetReferencePoQ()), 2);
		//v1 += pow(vec_PoQ_a.Z() - gMDFTrackerS522->GetReferencePoQ(), 2);
		//v2 += pow(vec_PoQ_a.Z() - foot_angles_i.Z(), 2);
		//v3 += pow(vec_PoQ_a.X() - foot_angles_i.X(), 2);
		//v4 += pow(vec_PoQ_a.Y() - foot_angles_i.Y(), 2);

		v3 += pow((vec_PoQ_a.X()/vec_PoQ_a.Z()) - (0.003445), 2);
		v4 += pow((vec_PoQ_a.Y()/vec_PoQ_a.Z()) - (-0.004355), 2);

		//cout<<"MDF_f1x-> "<<gMDFTrackerS522->f1_point_i.X()<<endl;
		//cout<<"MDF_f1z-> "<<gMDFTrackerS522->f1_point_i.Z()<<endl;
		cout<<"ref_PoQ-> "<<gMDFTrackerS522->GetReferencePoQ()<<endl;
		//cout<<"V2-> "<<v2<<endl;
		//cout<<"MDF_PoQ "<<gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input)<<endl;
		counter++;
	}
	v1 /=counter;
	v2 /= counter;
	v3 /= counter;
	v4 /= counter;
	v = sqrt(v2);
	v1_1 = sqrt(v1);
	v3_1 = sqrt(v3);
	v4_1 = sqrt(v4);
	v_tot = v+v3_1+v4_1;
	cout<<"v1-> "<<v1_1<<endl;	
	cout<<"v3-> "<<v3_1<<endl;	
	cout<<"v4-> "<<v4_1<<endl;	
	cout<<"v-> "<<v<<endl;	
	cout<<"vtot-> "<<v_tot<<endl;	

	std::cout << "\nReturning error: " << v_tot;
	//	for(int i=0; i<25; i++){cout<<"Par-> "<<par[i]<<endl;}
	return v_tot;

}

ClassImp(R3BTrackingS522);
