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

#include "R3BTrackingS091.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BFrsData.h"
#include "R3BTttxHitData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "R3BShared.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TFile.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <R3BCoarseTimeStitch.h>
#include "THttpServer.h"

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
R3BTrackingS091* gMDFTrackerS522;

R3BTrackingS091::R3BTrackingS091()
	: R3BTrackingS091("TrackingS091", 1)
{
}

R3BTrackingS091::R3BTrackingS091(const char* name, Int_t iVerbose)
	: FairTask(name, iVerbose)
	, fTrigger(-1)
	, fTpat(-1)
	, fNEvents(0)
	, maxevent(0)
	  , DoAlignment(false)
	, fTrackItems(new TClonesArray("R3BTrack"))
	, reference_PoQ(0.)
	, GladCurrent(-1)
	, GladReferenceCurrent(-1)
	, FiberTimeMin(-1)
	, FiberTimeMax(-1)
	, FiberEnergyMin(-1)
	, FiberEnergyMax(-1)
	  , fHeader(nullptr)
{
}

R3BTrackingS091::~R3BTrackingS091()
{
	if (fTrackItems){
		delete fTrackItems;}
}

InitStatus R3BTrackingS091::Init()
{
	LOG(info) << "R3BTrackingS091::Init()";
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
		LOG(warn) << "R3BTrackingS091::Init() EventHeader. not found";
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
	}
	// check if all cuts are properly set
	if (GladCurrent < 0 || GladReferenceCurrent < 0 || 
			FiberEnergyMin < 0 || FiberEnergyMax < 0)
	{
		R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
	}
	// Initializing all MDF functions
	LOG(info) << "Reading MDF function for TX0";
	MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

	LOG(info) << "Reading MDF function for FlightPath";
	MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

	LOG(info) << "Reading MDF function for TY0";
	MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());

	LOG(info) << "Reading MDF function for TX1";
	MDF_TX1 = new R3BMDFWrapper(MDF_TX1_filename.Data());

	LOG(info) << "Reading MDF function for TY1";
	MDF_TY1 = new R3BMDFWrapper(MDF_TY1_filename.Data());

	LOG(info) << "Reading MDF function for PoQ";
	MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());


	//Read output from the vertex macro
	// linking to global pointer (needed by alignment)
	gMDFTrackerS522 = this;
	// Request storage of R3BTrack data in the output tree
	mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);

	//online

	TFolder* mainfol = new TFolder("tracker", "tracker_info");
	FairRunOnline* run = FairRunOnline::Instance();
	run->GetHttpServer()->Register("", this);
	run->AddObject(mainfol);

	run->GetHttpServer()->RegisterCommand("Reset_tracker", Form("/Objects/%s/->Reset_Tracker_Histo()", GetName()));

	// inititallize canvas and histos
	trackerCanvas = new TCanvas("tracker_Canvas", "trackerCanvas");
	trackerCanvas->Divide(2,1);

	AoQ_Vs_Q_TOFD = R3B::root_owned<TH2F>("AoQ_Vs_Q_TOFD", "Using constant Beta", 500, 0., 5, 400, 0, 8);

	AoQ_Vs_Q_TOFD->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD->GetYaxis()->SetTitle("TOFD Q");

	AoQ_tof_Vs_Q_TOFD = R3B::root_owned<TH2F>("AoQ_tof_Vs_Q_TOFD", " Calculating beta with TOFD", 500, 0., 5, 400, 0, 8);

	AoQ_tof_Vs_Q_TOFD->GetXaxis()->SetTitle("AoQ");
	AoQ_tof_Vs_Q_TOFD->GetYaxis()->SetTitle("TOFD Q");

	trackerCanvas->cd(1);
	AoQ_Vs_Q_TOFD->Draw("COLZ");
	trackerCanvas->cd(2);
	AoQ_tof_Vs_Q_TOFD->Draw("COLZ");

	mainfol->Add(trackerCanvas);

	//////////// with tpat

	trackerCanvas_tpat = new TCanvas("tracker_Canvas_tpat", "trackerCanvas_tpat");
	trackerCanvas_tpat->Divide(2,1);

	AoQ_Vs_Q_TOFD_tpat = R3B::root_owned<TH2F>("AoQ_Vs_Q_TOFD_tpat", " Using constant Beta + reaction tpat", 500, 0., 5., 400, 0, 8);

	AoQ_Vs_Q_TOFD_tpat->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD_tpat->GetYaxis()->SetTitle("TOFD Q");

	AoQ_tof_Vs_Q_TOFD_tpat = R3B::root_owned<TH2F>("AoQ_tof_Vs_Q_TOFD_tpat", "Calculating beta with TOFD + reaction tpat", 500, 0., 5., 400, 0, 8);

	AoQ_tof_Vs_Q_TOFD_tpat->GetXaxis()->SetTitle("AoQ");
	AoQ_tof_Vs_Q_TOFD_tpat->GetYaxis()->SetTitle("TOFD Q");

	trackerCanvas_tpat->cd(1);
	AoQ_Vs_Q_TOFD_tpat->Draw("COLZ");
	trackerCanvas_tpat->cd(2);
	AoQ_tof_Vs_Q_TOFD_tpat->Draw("COLZ");

	mainfol->Add(trackerCanvas_tpat);
	/////////////////////////////////

	trackerCanvas_Vs_tpat = new TCanvas("tracker_Canvas_Vs_tpat", "trackerCanvas_Vs_tpat");

	AoQ_Vs_tpat = R3B::root_owned<TH2F>("AoQ_Vs_tpat", " Using constant Beta Vs reaction tpat", 16, 1., 16., 500, 1, 5);

	AoQ_Vs_tpat->GetXaxis()->SetTitle("Tpat");
	AoQ_Vs_tpat->GetYaxis()->SetTitle("AoQ");

	trackerCanvas_Vs_tpat->cd();
	AoQ_Vs_tpat->Draw("COLZ");

	mainfol->Add(trackerCanvas_Vs_tpat);
	//////////////////////
	// Folder for mapped data
	AoQ_vs_pos_both_fibCanvas = new TCanvas("AoQ_vs_pos_both_fib", "AoQ_vs_pos_both_fib");
	AoQ_vs_pos_both_fibCanvas->Divide(3,2);
	// Q = 6
	AoQ_vs_TOFD_pos_q_6 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_6", "AoQ_vs_TOFD_pos_q_6", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_6->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_6->GetYaxis()->SetTitle("AoQ");

	// Q = 5
	AoQ_vs_TOFD_pos_q_5 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_5", "AoQ_vs_TOFD_pos_q_5", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_5->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_5->GetYaxis()->SetTitle("AoQ");

	// Q = 4
	AoQ_vs_TOFD_pos_q_4 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_4", "AoQ_vs_TOFD_pos_q_4", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_4->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_4->GetYaxis()->SetTitle("AoQ");

	// Q = 3
	AoQ_vs_TOFD_pos_q_3 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_3", "AoQ_vs_TOFD_pos_q_3", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_3->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_3->GetYaxis()->SetTitle("AoQ");

	// Q = 2
	AoQ_vs_TOFD_pos_q_2 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_2", "AoQ_vs_TOFD_pos_q_2", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_2->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_2->GetYaxis()->SetTitle("AoQ");

	// Q = 1
	AoQ_vs_TOFD_pos_q_1 = R3B::root_owned<TH2F>("AoQ_vs_TOFD_pos_q_1", "AoQ_vs_TOFD_pos_q_1", 500, -100, 100, 500, 0, 5);

	AoQ_vs_TOFD_pos_q_1->GetXaxis()->SetTitle("TOFD_X");
	AoQ_vs_TOFD_pos_q_1->GetYaxis()->SetTitle("AoQ");


	AoQ_vs_pos_both_fibCanvas->cd(1);
	AoQ_vs_TOFD_pos_q_1->Draw("COLZ");
	AoQ_vs_pos_both_fibCanvas->cd(2);
	AoQ_vs_TOFD_pos_q_2->Draw("COLZ");
	AoQ_vs_pos_both_fibCanvas->cd(3);
	AoQ_vs_TOFD_pos_q_3->Draw("COLZ");
	AoQ_vs_pos_both_fibCanvas->cd(4);
	AoQ_vs_TOFD_pos_q_4->Draw("COLZ");
	AoQ_vs_pos_both_fibCanvas->cd(5);
	AoQ_vs_TOFD_pos_q_5->Draw("COLZ");
	AoQ_vs_pos_both_fibCanvas->cd(6);
	AoQ_vs_TOFD_pos_q_6->Draw("COLZ");

	mainfol->Add(AoQ_vs_pos_both_fibCanvas);
	///////////////////////
	
	trackerAnglesCanvas = new TCanvas("tracker_AnglesCanvas", "trackerAnglesCanvas");
	trackerAnglesCanvas->Divide(2,1);

	TX0_Vs_fib_pos_Q_5_AoZ_2 = R3B::root_owned<TH2F>("TX0_Vs_fib_pos_Q_5_AoQ_2", " TX0_Vs_fib_pos_Q_5_AoQ_2", 500, -500, -50, 500, -0.1, 0.1);
	TY0_Vs_fib_pos_Q_5_AoZ_2 = R3B::root_owned<TH2F>("TY0_Vs_fib_pos_Q_5_AoQ_2", " TY0_Vs_fib_pos_Q_5_AoQ_2", 500, -500, -50, 500, -0.1, 0.1);

	TX0_Vs_fib_pos_Q_5_AoZ_2->GetXaxis()->SetTitle("fib_pos");
	TX0_Vs_fib_pos_Q_5_AoZ_2->GetYaxis()->SetTitle("TX0");
	TY0_Vs_fib_pos_Q_5_AoZ_2->GetXaxis()->SetTitle("fib_pos");
	TY0_Vs_fib_pos_Q_5_AoZ_2->GetYaxis()->SetTitle("TY0");

	trackerAnglesCanvas->cd(1);
	TX0_Vs_fib_pos_Q_5_AoZ_2->Draw("colz");
	trackerAnglesCanvas->cd(2);
	TY0_Vs_fib_pos_Q_5_AoZ_2->Draw("colz");

	mainfol->Add(trackerAnglesCanvas);

	//// Neutron
	TFolder* neutronfol = new TFolder("neutron", "neutron");
	mainfol->Add(neutronfol);

	betaCanvas = new TCanvas("BetasCanvas", "BetaCanvas");

	neutron_beta = R3B::root_owned<TH1F>("neutron_beta", "neutron_beta", 500, 0, 1);

	neutron_beta->GetXaxis()->SetTitle("beta");

	betaCanvas->cd();
	neutron_beta->Draw("hist");

	neutronfol->Add(betaCanvas);

	neutronAnglesCanvas = new TCanvas("neutron_AnglesCanvas", "neutronAnglesCanvas");
	neutronAnglesCanvas->Divide(2,1);

	TX_pos_Q_5_AoZ_2 = R3B::root_owned<TH1F>("TX_pos_Q_5_AoQ_2", " TX_pos_Q_5_AoQ_2", 500, -0.1, 0.1);
	TY_pos_Q_5_AoZ_2 = R3B::root_owned<TH1F>("TY_pos_Q_5_AoQ_2", " TY_pos_Q_5_AoQ_2", 500, -0.1, 0.1);

	TX_pos_Q_5_AoZ_2->GetXaxis()->SetTitle("TX0");
	TY_pos_Q_5_AoZ_2->GetXaxis()->SetTitle("TY0");

	neutronAnglesCanvas->cd(1);
	TX_pos_Q_5_AoZ_2->Draw("hist");
	neutronAnglesCanvas->cd(2);
	TY_pos_Q_5_AoZ_2->Draw("hist");

	neutronfol->Add(neutronAnglesCanvas);

	e_relCanvas = new TCanvas("e_relsCanvas", "e_relCanvas");

	neutron_e_rel = R3B::root_owned<TH1F>("neutron_e_rel", "neutron_e_rel", 500, 0, 15);

	neutron_e_rel->GetXaxis()->SetTitle("e_rel");

	e_relCanvas->cd();
	neutron_e_rel->Draw("hist");

	neutronfol->Add(e_relCanvas);

	return kSUCCESS; 
}

void R3BTrackingS091::Reset_Tracker_Histo(){

	AoQ_Vs_Q_TOFD->Reset();
	AoQ_Vs_Q_TOFD->Reset();
	AoQ_tof_Vs_Q_TOFD->Reset();
	AoQ_Vs_Q_TOFD_tpat->Reset();
	AoQ_tof_Vs_Q_TOFD_tpat->Reset();
	AoQ_Vs_tpat->Reset();
	AoQ_vs_TOFD_pos_q_1->Reset();
	AoQ_vs_TOFD_pos_q_2->Reset();
	AoQ_vs_TOFD_pos_q_3->Reset();
	AoQ_vs_TOFD_pos_q_4->Reset();
	AoQ_vs_TOFD_pos_q_5->Reset();
	AoQ_vs_TOFD_pos_q_6->Reset();
	TX_pos_Q_5_AoZ_2->Reset();
	TY_pos_Q_5_AoZ_2->Reset();
	neutron_beta->Reset();
	neutron_e_rel->Reset();
	return;
}

void R3BTrackingS091::Exec(Option_t* option)
{
	if (fNEvents / 1000. == (int)fNEvents / 1000)
		std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
			<< " %) " << std::flush;
	//	FairRootManager* mgr = FairRootManager::Instance();
	//	R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");
	fNEvents += 1;
	is_good_event = false;

	h++;
	Tpat = fHeader->GetTpat();//vairiable in the output tree

	if(Tpat & 0xf000){
		return;
	}

	hh++;

	if(Tpat == 0){
		return;
	}
	Int_t tpatbin = 0;
	std::vector<int> tpatindex;
	for (int i = 0; i < 16; i++)
	{
		tpatbin = (Tpat & (1 << i));
		if (tpatbin != 0)
			tpatindex.push_back(i + 1);
	}
	mul_los=-999;
	mul_m0=-999;
	mul_m1=-999;
	mul_f32=-999;
	mul_f30=-999;
	mul_f31=-999;
	mul_f33=-999;
	mul_tofd=-999;
	mul_frsi=-999;
	cond=false;

	mul_los  = fDataItems[LOS_DATA]->GetEntriesFast();
	mul_m0   = fDataItems[MWPC0_HITDATA]->GetEntriesFast();
	mul_m1   = fDataItems[MWPC1_HITDATA]->GetEntriesFast();
	mul_tttx = fDataItems[TTTX_HITDATA]->GetEntriesFast();
	mul_f32  = fDataItems[DET_FI32]->GetEntriesFast();
	mul_f30  = fDataItems[DET_FI30]->GetEntriesFast();
	mul_f31  = fDataItems[DET_FI31]->GetEntriesFast();
	mul_f33  = fDataItems[DET_FI33]->GetEntriesFast();
	mul_tofd = fDataItems[DET_TOFD]->GetEntriesFast();
	mul_frsi = fDataItems[FRS_DATA]->GetEntriesFast();
	mul_neuland = fDataItems[NEULAND_DATA]->GetEntriesFast();
	hhh++;
	if(mul_los!=1){
		return;
	}
	a++;
	if(mul_tofd<1){
		return;
	}
	b++;
	//if(mul_m0!=1)return;
	if(mul_m1!=1){
		return;
	}
	c++;
	if(mul_f32 < 1 || mul_f30 < 1 || (mul_f31 < 1 && mul_f33 < 1)){
		return;
	}
	d++;
	//	auto frs_DataItems = fDataItems.at(FRS_DATA);
	//	if(frs_DataItems->GetEntriesFast() < 1) return; 
	//	hhh++;
	//	auto frs_data = (R3BFrsData*)frs_DataItems->At(0);
	//	frs_beta->Fill(frs_data->GetBeta());
	//	//if(frs_data->GetBeta()<0.715 || frs_data->GetBeta()>0.705) return;
	//	//if(frs_data->GetZ()<5.2 || frs_data->GetZ()>6.7) return;

	double mintime = 2000000; 
	double minE = -1;
	TVector3 minvec;
	bool n_true = false;
	R3BNeulandHit* hit{};
	for (auto i = 0; i < fDataItems[NEULAND_DATA]->GetEntriesFast(); ++i)
	{
		hit = static_cast<R3BNeulandHit*>(fDataItems[NEULAND_DATA]->At(i));
		const Double_t tcorr = hit->GetT() - (hit->GetPosition().Mag() - 1564.) / 29.9792458;
		if(hit->GetE() > 0. && tcorr > 59.5 && tcorr < 96.)
		{
			if(hit->GetT() < mintime)
			{
				mintime = hit->GetT();
				minvec = hit->GetPosition();
				minE = hit->GetE();
				n_true = true;
			}
		}
	}
	//------ Get TOFD data 
	R3BTofdHitData* tofd_hit{};
	int mul_tofd1=0;
	double tofdq_temp =0;
	double tttxq_temp =0;
	double tofd_tof_temp =0;
	double tofd_pos =0;
	bool is_good_tofd = false;
	bool is_good_tttx = false;

	for (auto i = 0; i < fDataItems[DET_TOFD]->GetEntriesFast(); ++i)
	{
		tofd_hit = static_cast<R3BTofdHitData*>(fDataItems[DET_TOFD]->At(i));
		if (tofd_hit->GetDetId() == 1 && tofd_hit->GetTof() < 120 && tofd_hit->GetTof() > 100) // only hits from first plane, add Z later
		{
			is_good_tofd = true;
			mul_tofd1++;
			tofdq_temp=tofd_hit->GetEloss() - 0.2;
			tofd_pos=tofd_hit->GetX();
			tofd_tof_temp=tofd_hit->GetTof() - 104.2 + tof_offset;
			//break;
		}
	}
	if(!is_good_tofd || mul_tofd1!=1){ 
		return;
	}
	e++;
	if(!MakeIncomingTracks()){ 
		return;//at least one good track candidate in FOOT
	}
	f++;
	if(!MakeOutgoingTracks()){ 
		return;//at least one good track candidate in Fibers
	}
	g++;
	//cout << "\nGood event\n";
	is_good_event = true;
	cond=true;
	int counter = 0;
	double aoq_offset =0;
	double TX0_offset =0;
	double TY0_offset =0;

	if(tracks_out.size()==2){
		if(tracks_out[0].fiber == tracks_out[1].fiber){
			return;
		}
	}
	if(tracks_in.size()>1 || tracks_out.size()>2){
		return;
	}
	double delta_TX1, delta_TX0, delta_TY0;
	for (auto & tin : tracks_in){
		for (auto & tout : tracks_out){
			counter++;
			if(tracks_in.size() == 2 && !tout.fiber){
				continue;
			}
			//preserve the order, it is expected by the MDF function!
			mdf_data[0] = tin.mw1_x;
			mdf_data[1] = tin.mw1_y;
			mdf_data[2] = tin.mw1_z;
			mdf_data[3] = tout.f32_x;
			mdf_data[4] = tout.f32_z;
			mdf_data[5] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
			mdf_data[6] = (tout.f30_y  - tin.mw1_y)/(tout.f30_z - tin.mw1_z);
			// Calculate all required MDF values

			flight_p = MDF_FlightPath->MDF(mdf_data);
			poq = MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent;
			tx0 = MDF_TX0->MDF(mdf_data);
			ty0 = MDF_TY0->MDF(mdf_data);
			//		tx1 = MDF_TX1->MDF(mdf_data);
			//		ty1 = MDF_TY1->MDF(mdf_data);
			tof = flight_p / FRS_BETA / SPEED_OF_LIGHT;

			beta = FRS_BETA;
			gamma = 1. / sqrt(1 - pow(beta, 2));
			maoz = poq / beta / gamma / AMU;

			beta_tof = flight_p / tofd_tof_temp / SPEED_OF_LIGHT;
			gamma_tof = 1. / sqrt(1 - pow(beta_tof, 2));
			maoz_tof = poq / beta_tof / gamma_tof / AMU;

			if(counter == 1){
				TY0_offset = .00432;//10C
				if(tout.fiber){
					aoq_offset = -.1;//10C
					TX0_offset = .008 + 0.0055;//10C
					//aoq_offset = -.13;//16C
				}
				else{
					//aoq_offset = -.29;//12C
					TX0_offset = .02 + 0.008;//10C
					aoq_offset = -.25;//10C
					//aoq_offset = -.39;//16C
				}
				maoz = maoz + aoq_offset;
				maoz_tof = maoz_tof + aoq_offset;
				AoQ_Vs_Q_TOFD->Fill(maoz,tofdq_temp);
				AoQ_tof_Vs_Q_TOFD->Fill(maoz_tof,tofdq_temp);
				if((Tpat & 16) == 16 || (Tpat & 32) == 32){
					AoQ_Vs_Q_TOFD_tpat->Fill(maoz,tofdq_temp);
					AoQ_tof_Vs_Q_TOFD_tpat->Fill(maoz_tof,tofdq_temp);
				}
				for (const auto& itpat : tpatindex){
					AoQ_Vs_tpat->Fill(itpat,maoz);  
				}
				if(fabs(tofdq_temp - 1) < 0.5){
					AoQ_vs_TOFD_pos_q_1->Fill(tofd_pos,maoz);
				}
				if(fabs(tofdq_temp - 2) < 0.5){
					AoQ_vs_TOFD_pos_q_2->Fill(tofd_pos,maoz);
				}
				if(fabs(tofdq_temp - 3) < 0.5){
					AoQ_vs_TOFD_pos_q_3->Fill(tofd_pos,maoz);
				}
				if(fabs(tofdq_temp - 4) < 0.5){
					AoQ_vs_TOFD_pos_q_4->Fill(tofd_pos,maoz);
				}
				if(fabs(tofdq_temp - 5) < 0.5){
					AoQ_vs_TOFD_pos_q_5->Fill(tofd_pos,maoz);
					if(fabs(maoz-2.8)<0.1){
						TX0_Vs_fib_pos_Q_5_AoZ_2->Fill(tout.f32_x,tx0 + TX0_offset);
						TY0_Vs_fib_pos_Q_5_AoZ_2->Fill(tout.f32_x,ty0 + TY0_offset);
						if(n_true)
						{
							Double_t m_neut = 939.565;
							Double_t dalt = 931.494;
							Double_t m_f = dalt*(14.025404 - 0.0027429) ; //AME20 //14B
							Double_t nx_corr = 0., ny_corr = 0., nz_corr = 0.;
							Double_t dx_neu = 0., dy_neu = 0., dz_neu = 0.;
							Double_t fx_corr = 0., fy_corr = 0., fz_corr = 0.;
							Double_t cos_ang = 0.;
							Double_t e_rel =0;
							nz_corr = minvec.Z()*10.;
							dx_neu = (minvec.X()*10.)/nz_corr;
							dx_neu = dx_neu + 0.003365;
							dy_neu = (minvec.Y()*10.)/nz_corr;
							dy_neu = dy_neu - 0.003572;
							fx_corr = tx0 + TX0_offset;
							fy_corr = ty0 + TY0_offset;
							fz_corr = 1.;
				
							cos_ang = ((dx_neu*fx_corr) + (dy_neu*fy_corr) + 1.)/(sqrt(dx_neu*dx_neu + dy_neu*dy_neu + 1.)*sqrt(fx_corr*fx_corr + fy_corr*fy_corr + 1.));
							Double_t r_neu = nz_corr*sqrt(pow(dx_neu,2) + pow(dy_neu,2) + 1.);
							Double_t beta_neu = (r_neu/1000.)/(mintime);
							Double_t beta_frag = 0.7127;
							beta_neu = beta_neu/(TMath::C() / pow(10,9));
							Double_t gamma_neu = 1./sqrt(1. - pow(beta_neu,2));
							Double_t gamma_frag = 1./sqrt(1. - pow(beta_frag,2));
							e_rel = sqrt(m_f*m_f + m_neut*m_neut + 2*gamma_neu*gamma_frag*m_f*m_neut*(1-beta_neu*beta_frag*cos_ang)) - m_f - m_neut;

							neutron_beta->Fill(beta_neu);
							neutron_e_rel->Fill(e_rel);
							TX_pos_Q_5_AoZ_2->Fill(dx_neu);
							TY_pos_Q_5_AoZ_2->Fill(dy_neu);

						}
					}
				}
				if(fabs(tofdq_temp - 6) < 0.5){
						AoQ_vs_TOFD_pos_q_6->Fill(tofd_pos,maoz);
				}
			}
			//TVector3 vec_PoQ(0, 0, 1);
			TVector3 vec_PoQ(tx0 + TX0_offset, ty0, 1);
			vec_PoQ.SetMag(poq);
			AddTrackData(tin.mw1_x, tin.mw1_y, tin.mw1_z, vec_PoQ, tofdq_temp, maoz); // chix, chiy, quality
		}
	}
	return;
}

void R3BTrackingS091::FinishEvent()
{
	for (auto& DataItem : fDataItems)
	{
		DataItem->Clear();
	}
    	fTrackItems->Clear();
	if (fNEvents / 10000. == (int)fNEvents / 10000)
		cout << " \n finish event " 
			<< " \n before any cuts " << h 
			<< " \n survives offspil tpat => " << hh 
			<< " \n survives tpat 0 => " << hhh 
			//		<< " \n survives frs != 0 => " << hhh 
			<< " \n survives los==1 => " << a 
			<< " \n survives mul!=0 in tofd => " << b 
			<< " \n survives mul!=0 in mwpc1 => " << c 
			<< " \n survives fiber mul in atleast 3 fib => " << d 
			<< " \n survives there is at least one hit in tofd plane 1 => " << e 
			<< " \n survives mwpc1 is not nan => " << f 
			<< " \n survives survived the clustering in fibs => " << g 
			<< endl;
}

void R3BTrackingS091::FinishTask()
{
	LOG(info) << "Processed " << fNEvents << " events\n\n";
	AoQ_Vs_Q_TOFD->Write();
	AoQ_vs_TOFD_pos_q_1->Write();
	AoQ_vs_TOFD_pos_q_2->Write();
	AoQ_vs_TOFD_pos_q_3->Write();
	AoQ_vs_TOFD_pos_q_4->Write();
	AoQ_vs_TOFD_pos_q_5->Write();
	AoQ_vs_TOFD_pos_q_6->Write();
	AoQ_Vs_tpat->Write();
	AoQ_Vs_Q_TOFD_tpat->Write();
	AoQ_tof_Vs_Q_TOFD_tpat->Write();
	TX0_Vs_fib_pos_Q_5_AoZ_2->Write();
	TY0_Vs_fib_pos_Q_5_AoZ_2->Write();
	TX_pos_Q_5_AoZ_2->Write();
	TY_pos_Q_5_AoZ_2->Write();
	neutron_beta->Write();
	neutron_e_rel->Write();
	//cout<<"WRITE"<<endl;
}


bool R3BTrackingS091::MakeIncomingTracks()
{
	tracks_in.clear();
	TVector3 vertex_mwpc;
	Track tr;
	//Get MWPC hits, for now only first hit
	auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
	m1_point.SetXYZ(m1_hit->GetX()*0.1, m1_hit->GetY()*0.1, 0.);//cm

	TransformPoint(m1_point, &m1_angles, &m1_position);//lab
	tr.mw1_x   = m1_point.X();
	tr.mw1_y   = m1_point.Y();
	tr.mw1_z   = m1_point.Z();
	if(isnan(m1_point.X()) || isnan(m1_point.Y()) || isnan(m1_point.Z())){
		return false;
	}
	tracks_in.push_back(tr);
	return true;
}

void R3BTrackingS091::TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
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

R3BTrack* R3BTrackingS091::AddTrackData(double x, double y, double z, TVector3 poq_vec, Double_t charge, Double_t aoz)
{
	// Filling output track info
	add_track_counter++;
	TClonesArray& clref = *fTrackItems;
	Int_t size = clref.GetEntriesFast();
	return new (clref[size]) R3BTrack(x, y, z, poq_vec.X(),  poq_vec.Y(), poq_vec.Z(), charge, aoz, 0., h, add_track_counter);
}

bool R3BTrackingS091::IsGoodFiberHit(R3BFiberMAPMTHitData* fhit)
{
	if((fhit->GetEloss() > FiberEnergyMin) && (fhit->GetEloss() < FiberEnergyMax) && 
			(fhit->GetTime() < 20000 && fhit->GetTime()>(-20000) )){
		return true;
	}
	else{
		return false;
	}
}
bool R3BTrackingS091::MakeOutgoingTracks()
{
	if(fDataItems[DET_FI32]->GetEntriesFast() > 3 || fDataItems[DET_FI30]->GetEntriesFast() > 3 || (fDataItems[DET_FI33]->GetEntriesFast() > 3 || fDataItems[DET_FI31]->GetEntriesFast() > 3)){
		return false;
	}

	tracks_out.clear();
	Track tr;
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
		if(!IsGoodFiberHit(f32)){
		       	continue;
		}
		double fitime = f32->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{	f32x.push_back(f32->GetX());
			f32e.push_back(f32->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI30]->GetEntriesFast(); ++i)
	{
		auto f30 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI30]->At(i));
		if(!IsGoodFiberHit(f30)){
		       	continue;
		}
		double fitime = f30->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			f30y.push_back(f30->GetY());
			f30e.push_back(f30->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI33]->GetEntriesFast(); ++i)
	{
		auto f33 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI33]->At(i));
		if(!IsGoodFiberHit(f33)){
		       	continue;
		}

		double fitime = f33->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			flastx.push_back(f33->GetX());
			flaste.push_back(f33->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI31]->GetEntriesFast(); ++i)
	{
		auto f31 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI31]->At(i));
		if(!IsGoodFiberHit(f31)){
		       	continue;
		}
		double fitime = f31->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			flast2x.push_back(f31->GetX());
			flast2e.push_back(f31->GetEloss());
		}
	}
	if(f32x.size() < 1 || f30y.size() < 1 || (flastx.size() < 1 && flast2x.size() < 1)){
		return false;
	}
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
		else{
			continue;
		}
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
		else{
			continue;
		}
		for(Int_t k = 0; k < num_hit0[hit_clust_id0[i]]; k++)
		{
			for(Int_t j = 0; j < f30y.size(); j++)
			{
				if(set0[j]){
					continue;
				}
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
		else{
			continue;
		}
		for(Int_t k = 0; k < num_hit3[hit_clust_id3[i]]; k++)
		{
			for(Int_t j = 0; j < flastx.size(); j++)
			{
				if(set3[j]){
					continue;
				}

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
		else{
			continue;
		}
		for(Int_t k = 0; k < num_hit1[hit_clust_id1[i]]; k++)
		{
			for(Int_t j = 0; j < flast2x.size(); j++)
			{
				if(set1[j]){
					continue;
				}

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
		if(sum_ener > 0.){
			fi32xcl[i] = (double)fi32xcl[i]/sum_ener;
		}
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
		if(sum_ener > 0.){
			fi30ycl[i] = (double)fi30ycl[i]/sum_ener;
		}
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
		if(sum_ener > 0.){
			filastxcl[i] = (double)filastxcl[i]/sum_ener;
		}
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
		if(sum_ener > 0.){
			filast2xcl[i] = (double)filast2xcl[i]/sum_ener;
		}
	}
	for(auto i = 0; i < num_clust2; i ++)
	{
		f32_point.SetXYZ(fi32xcl[i], 0, 0); //cm
		TransformPoint(f32_point, &f32_angles, &f32_position);
		tr.f32_x = f32_point.X();
		tr.f32_z = f32_point.Z();

		for (auto j=0; j<num_clust0; ++j)
		{
			f30_point.SetXYZ(0,fi30ycl[j],0); //cm
			TransformPoint(f30_point, &f30_angles, &f30_position);
			tr.f30_y = f30_point.Y();
			//make combination with every hit in fibers 33 and 31

			for (auto k = 0; k<num_clust3; ++k)
			{
				flast_point.SetXYZ(filastxcl[k], 0, 0); //cm
				TransformPoint(flast_point, &f33_angles, &f33_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//if(angle_out>(-10.) || angle_out<(-18.)) continue;
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
				tr.fiber = true;
				tracks_out.push_back(tr);
			}

			for (auto k = 0; k<num_clust1; ++k)
			{
				flast_point.SetXYZ(filast2xcl[k], 0, 0); //cm
				TransformPoint(flast_point, &f31_angles, &f31_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//if(angle_out>(-10.) || angle_out<(-18.)) continue;
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
				tr.fiber = false;
				tracks_out.push_back(tr);
			}
		}
	}
	if(tracks_out.empty()){
	       	return false;
	}
	return true;
}

ClassImp(R3BTrackingS091);
