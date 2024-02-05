/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "TClonesArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"

#include "R3BRpcCal2HitPar.h"
#include "R3BRpcCal2Hit.h"
#include "R3BRpcCalData.h"
#include "R3BEventHeader.h"

#include "R3BRpcHitPar.h"

#include <iostream>
#include <stdlib.h>

R3BRpcCal2HitPar::R3BRpcCal2HitPar()
    : R3BRpcCal2HitPar("R3B Rpc Tot Calibration Parameters Finder ", 1)
{
}

R3BRpcCal2HitPar::R3BRpcCal2HitPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitPar(NULL)
    , fCalDataCA(NULL)
    , fNumChannels(64)
    , fR3BEventHeader(NULL)
    , fDebugMode(false)
{
 for (Int_t i = 0; i < N_STRIP_NB; i++)
 {
  fhPos[i] = NULL;
  fhPar[i] = NULL;
  fhTime[i] = NULL;
 }
}

R3BRpcCal2HitPar::~R3BRpcCal2HitPar()
{
    LOG(info) << "R3BRpcCal2HitPar: Delete instance";
    if (fCalDataCA)
        delete fCalDataCA;
}

void R3BRpcCal2HitPar::SetParContainers()
{
 // Parameter Container
 // Reading RPCHitPar from FairRuntimeDb
 FairRuntimeDb* rtdb = FairRuntimeDb::instance();
 if (!rtdb)
 {
     LOG(error) << "FairRuntimeDb not opened!";
 }
}

void R3BRpcCal2HitPar::SetParameter() {}

InitStatus R3BRpcCal2HitPar::Init()
{
 LOG(info) << "R3BRpcCal2HitPar::Init()";

 FairRootManager* rootManager = FairRootManager::Instance();
 if (!rootManager)
 {
     LOG(error) << "R3BRpcCal2HitPar::Init() FairRootManager not found";
     return kFATAL;
 }

 fCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcCalData"));
 if (!fCalDataCA)
 {
     LOG(error) << "R3BRpcCal2HitPar::Init() R3BRpcCalData not found";
     return kFATAL;
 }

 fR3BEventHeader = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
 if (!fR3BEventHeader)
 {
     LOG(error) << "R3BRpcCal2HitPar::Init() EventHeader. not found";
     return kFATAL;
 }

 FairRuntimeDb* rtdb = FairRuntimeDb::instance();
 if (!rtdb)
 {
     LOG(error) << "R3BRpcCal2HitPar::Init() FairRuntimeDb not found";
     return kFATAL;
 }

 fHitPar = dynamic_cast<R3BRpcHitPar*>(rtdb->getContainer("RpcHitPar"));
 if (!fHitPar)
 {
     LOG(error) << "R3BRpcCal2HitPar::Init() Couldn't get handle on RPCHitPar container";
     return kFATAL;
 }

 // Set container with mapping parameters
 SetParameter();

 return kSUCCESS;
}

InitStatus R3BRpcCal2HitPar::ReInit()
{
 SetParContainers();
 SetParameter();
 return kSUCCESS;
}

void R3BRpcCal2HitPar::Exec(Option_t* opt)
{
 //loop over strip data
 Int_t nHits = fCalDataCA->GetEntries();
 UInt_t iDetector = 0;
 double time_rpc=0;
 int chn_rpc= 0;
 double pos_rpc= 0;
 int chn_scint= 0;
 double pos_scint= 0;
 double time_scint = 0;
 bool rpc_hit=false;
 bool scint_hit=false;
 bool tpat1 = fR3BEventHeader->GetTpat() & 0xf000;
 bool tpat2 = fR3BEventHeader->GetTpat() & 0x0fff;
 UInt_t inum=0 ;
  for (Int_t i = 0; i < nHits; i++){
   auto map1 = dynamic_cast<R3BRpcCalData*>(fCalDataCA->At(i));
   iDetector=map1->GetDetId();
   inum = iDetector * 41 + map1->GetChannelId() -1;
   if(iDetector==0){
    time_rpc = (map1->GetTimeR_B()+map1->GetTimeL_T())/2.;
    chn_rpc = map1->GetChannelId();
    pos_rpc = (map1->GetTimeL_T() - map1->GetTimeR_B());
    rpc_hit = true;
    if (tpat2>0){ 
     if (NULL == fhPos[inum]){
      char strName[255];
      char strName_par[255];
      sprintf(strName, "%s_poscaldata_%d", fHitPar->GetName(),inum);
      sprintf(strName_par, "%s_parcaldata_%d", fHitPar->GetName(),inum);
      fhPos[inum] = new TH1F(strName, "", 4000,-2000.,2000.);
      fhPar[inum] = new TH1F(strName_par, "", 4000,0.,4000.);
     }
     fhPos[inum]->Fill(pos_rpc*CSTRIP/2.);
    }
    if (tpat2>0){ 
     if (NULL == fhTime[inum]){
      char strName[255];
      sprintf(strName, "%s_timecaldata_%d", fHitPar->GetName(),inum);
      fhTime[inum] = new TH1F(strName, "", 500,-40,-15);
     }
    } 
   }
   if(iDetector==1 && tpat2>0){
    time_scint = (map1->GetTimeR_B()+map1->GetTimeL_T())/2.;
    chn_scint = map1->GetChannelId();
    pos_scint = (map1->GetTimeL_T() - map1->GetTimeR_B());
    if(chn_scint == 1){
     scint_hit = true;
    }
    if (NULL == fhPos[inum]){
     char strName[255];
     sprintf(strName, "%s_poscaldata_%d", fHitPar->GetName(),inum);
     fhPos[inum] = new TH1F(strName, "", 5000,-2500.,2500.);
    }
    fhPos[inum]->Fill(pos_scint*CSCINT/2.);
   }
  }
  if(rpc_hit && scint_hit){
   fhTime[chn_rpc -1]->Fill(time_rpc-time_scint);
  }
 return;
}

void R3BRpcCal2HitPar::Reset() {}

void R3BRpcCal2HitPar::FinishEvent() {}

void R3BRpcCal2HitPar::FinishTask() {

 CalculateParsPosStrip();
 CalculateParsTimeStrip();
 fHitPar->setChanged();
 fHitPar->printParams();
 fHitPar->Write();
}

void R3BRpcCal2HitPar::CalculateParsPosStrip(){
 for (int t = 0; t < N_STRIP_NB; t++){
  if (NULL == fhPos[t] || t>=41 || NULL == fhPar[t]){continue;}
  double sum = 0; 
  for(int i = 1; i <=4000 ; i++){
   sum += fhPos[t]->GetBinContent(i);
   fhPar[t]->Fill(i,sum);
  }
  fhPar[t]->Scale(fhPar[t]->GetXaxis()->GetBinWidth(1)/(sum));
  fHitPar->SetCalParams1(fhPar[t]->FindFirstBinAbove(0.01,1),t);
 }
}

void R3BRpcCal2HitPar::CalculateParsTimeStrip(){
 for (int t = 0; t < 41; t++){
  if (NULL == fhTime[t]){continue;}
  if(fhTime[t]->GetEntries()==0){continue;}
  fhTime[t]->Fit("gaus");
  TF1 *g6 = (TF1*)fhTime[t]->GetListOfFunctions()->FindObject("gaus");
  fHitPar->SetCalParams2(g6->GetParameter(1),t);
 }
}
ClassImp(R3BRpcCal2HitPar)
