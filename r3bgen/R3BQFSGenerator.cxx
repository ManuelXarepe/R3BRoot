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

/* A + i = (B+a) + i = B + (a->i) = B + a + i 
 *
 *
 * 
 */

#include "R3BQFSGenerator.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <FairLogger.h>
#include <FairPrimaryGenerator.h>

#include <TMath.h>

/////////////////Kai's Modification//////////////////
#include "FairIon.h"
#include "FairRunSim.h"




using namespace std;

R3BQFSGenerator::R3BQFSGenerator()
{
	this->SetValues(0, 0, 0, 0, 0, false, false, false);
//	T_LIMIT = 10000000.;
	fRandom.SetSeed(0);
	//    GammaEnergy = 0;
	cout<<"\n ##Starting QFS Gen##"<<endl;
}

Bool_t R3BQFSGenerator::Init()
{
	fRandom.SetSeed(0);
	//foo_Cugnon_DXS = new TF1("foo_Cugnon_DXS","");
	return kTRUE;
}

///////////////////////////////////////Kai's Modification/////////////////////////////////////


/////THE OLD CODE/////

// void R3BQFSGenerator::SetHeavyNucleus(int A, double M_def, double M_res)
// {
//     A_BEAM = A;
//     double UNIT = 931.494061;
//     MA = A * Ungle between two vectorNIT + M_def;
//     MB = (A - 1) * UNIT + M_res;
//     return;
// }

/////THE OLD CODE/////

void R3BQFSGenerator::SetAddGamma(double energy)
{
	AddGamma = true;
	GammaEnergy = energy;
	SetExcitation(GammaEnergy);
	return;
}

void R3BQFSGenerator::SetProjectile(int A, int Z, double E)
{
	fIonBeam = new FairIon(TString::Format("FairIon_%d_%d_%d", Z, A, Z), Z, A, Z);
	MA = fIonBeam->GetMass() * 1000.;//GeV to MeV
	A_BEAM = A;
	ENERGY = E;
	//cout << "\nProjectile mass = " << MA << "\n\n";
	return;
}

void R3BQFSGenerator::SetHeavyFragment(int A, int Z)
{
	fIonFragment = new FairIon(TString::Format("FairIon_%d_%d_%d", Z, A, Z), Z, A, Z);
	MB = fIonFragment->GetMass() * 1000.; //Gev to MeV
	//cout << "\nFragment mass = " << MB << "\n\n";
	//auto run = FairRunSim::Instance();
	//run->AddNewIon(fIonFragment);
	return;
}

///////////////////////////////////////Kai's Modification/////////////////////////////////////

void R3BQFSGenerator::SetLightNucleus(double ma, double mi, bool is_pn)
{
	IS_PN = is_pn;
	Ma = ma;
	Mi = mi;
	return;
}

void R3BQFSGenerator::SetExcitation(double exe)
{
	MB += exe;
	return;
}

void R3BQFSGenerator::SetInverse(bool inverse)
{
	INVERSE = inverse;
	return;
}
void R3BQFSGenerator::SetIsotropic(bool isotropic)
{
	ISOTROPIC = isotropic;
	return;
}

void R3BQFSGenerator::SetMomDistrib(double mom)
{
	MOM_SIGMA = mom;
	return;
}

void R3BQFSGenerator::SetVertexDistrib(TVector3 mean, TVector3 sigma)
{
	mean_vertex = mean;
	sigma_vertex = sigma;
	return;
}

/////////////////////////////////////////////08.08.2024////////////////////////////////////////////////
//Create Gamma Emission in the cm frame
std::vector<Double_t> R3BQFSGenerator::GetGamma()
{
//	cout << "energyofGamma  " << GammaEnergy << endl; 
	Double_t phi = fRandom.Uniform(0.,2.*TMath::Pi());
	Double_t costheta = fRandom.Uniform(-1,1);
	Double_t theta = TMath::ACos(costheta);

	Double_t energy = GammaEnergy ;
	Double_t px =energy * TMath::Sin(theta)*TMath::Cos(phi);
	Double_t py =energy * TMath::Sin(theta)*TMath::Sin(phi);
	Double_t pz =energy * TMath::Cos(theta);
	//SetExcitation(energy);
	return {px, py, pz, energy};
}

void R3BQFSGenerator::SetLorentzBoost(TLorentzVector& gammaMomentum, const TLorentzVector& fragmentMomentum)
{
	TVector3 boostVector = fragmentMomentum.BoostVector();
	gammaMomentum.Boost(boostVector);
	return;

}


///////////////////////////////////////////////08.08.2024///////////////////////////////////////////////


////////////////////////////////////////////////Kai's Modification///////////////////////////////////////////

void R3BQFSGenerator::Print()
{
	cout << "***** Print generator values *****" << endl;
	cout << "Energy: \t" << ENERGY << endl;
	cout << "A Beam: \t" << A_BEAM << endl;
	cout << "Mass A: \t" << MA << endl;
	cout << "Mass B: \t" << MB << endl;
	cout << "Mass a: \t" << Ma << endl;
	cout << "Mass i: \t" << Mi << endl;
//	cout << "T_Limit: \t" << T_LIMIT << endl;
	cout << "Sigma Fermi: \t" << MOM_SIGMA << endl;
	cout << "WQ Isotropic: \t" << ISOTROPIC << endl;
	cout << "Inverse kinematics: \t" << INVERSE << endl;
	cout << "Enable Gamma: \t" << AddGamma << endl; //<<<<<<<<<<<<<<<<<<<<<<<<<<
	cout << "Gamma Energy: \t" << GammaEnergy << endl;
}

//void R3BQFSGenerator::get_ascii_input(std::vector <R3BTrackingParticle> *for_ascii)
std::vector <R3BTrackingParticle> R3BQFSGenerator::get_ascii_input()
{
	std::vector <R3BTrackingParticle> for_ascii;
	double E = ENERGY;
	// GenerateEvent(ENERGY,primGen);

	if (E <= 0)
	{
		LOG(error) << "R3BQFSGenerator: E < 0!";
		return for_ascii;
	}
	fRandom.SetSeed(0);
	//  	gRandom = new TRandom3(); //default
	//  	gRandom->SetSeed(0);//using computer time
	//  	TRandom2 r1;
	//  	r1.SetSeed(0);
	//R3BQFSGenerator::Print();
	bool evt = false;


	//-------------- The reaction 11C(p,pn)10C --------------------
	double Tkin = E * A_BEAM; // Total kinetic energy (MeV) of the projectile nucleus for Inverse Kinematics

	if(!INVERSE){Tkin = E;} //For Direct Kinematics

	const double sigma = MOM_SIGMA;
	TLorentzVector LVa; // Target
	TLorentzVector LVi; // Proton Beam 
	TLorentzVector* LVg = nullptr; // For Gamma

	double* gA = nullptr;
	double gi;
	double* bA = nullptr;
	double bi;
	double PA, Pi;
	double EA, Ei;
	double S_first;

	//Note: MA for inverse-kinematics beam, MB = fragment, Mi and Ma are proton or neutron
	if(INVERSE) //Beam Parameters for Inverse Kinematics
	{
		//--------------------- Beam parameters -----------------------------
		PA = sqrt(Tkin * (Tkin + 2 * MA));         // Total 3-momentum of the beam
		EA = sqrt(MA * MA + PA * PA);              // Total energy of the beam
		bA = new double(-PA / EA);                             // Beta of the beam
		gA = new double(1/TMath::Sqrt(1-*bA* *bA));                // Gamma of the beam
		S_first = (EA + Mi) * (EA + Mi) - PA * PA; // Invariant mass (Mandelstam S-variable)
		//	TLorentzVector LVa; // Proton Target
		//	TLorentzVector LVi; // Beam
		LVi.SetPxPyPzE(0., 0.,0.,Mi); //Lorentz Vector of the Beam
		//	TLorentzVector* LVg = nullptr; // For Gamma

		LOG(debug) << "\n****** Beam Parameters for Inverse Kinematics ********";
		LOG(debug) << "\nMA:\t" << MA << " MeV";
		LOG(debug) << "\nTotal momentum:\t" << PA << " MeV";
		LOG(debug) << "\nTotal energy:\t" << EA << " MeV";
		LOG(debug) << "\nBeta (beam):\t" << (-*bA) << "\nGamma (beam):\t" << *gA << "\n\n";
	}

	else //For Direct Kinematics
	{
		//Using the same variables' names here, so the following code is shorter
		Pi = sqrt(Tkin*(Tkin+2*Mi)); // Total 3 momentum of the proton beam
		Ei = sqrt(Mi*Mi + Pi*Pi);    // Total energy of the proton beam
		bi = -Pi/Ei;		  // Beta of the proton beam
		gi = 1/sqrt(1-bi*bi);	  // Gamma of the proton beam
		S_first = (Ei+MA)*(Ei+MA)-Pi*Pi; //Madelstam S-invariant of the proton beam

		//	TLorentzVector LVa; // Target
		//	TLorentzVector LVi; // Proton Beam 
		LVi.SetPxPyPzE(0.,0.,Pi, Ei); //of the Bram

		//	TLorentzVector* LVg = nullptr; //For Gamma Radition from the fragment

		LOG(debug) << "\n******Beam Parameters for Direct Kinematics******"; 
		LOG(debug) << "\nMi: \t" << Mi << "MeV"; //Mass of the proton beam
		LOG(debug) << "\nMA: \t" << MA << "MeV"; //Mass of A target
		LOG(debug) << "\nTotal Momentum:\t" << Pi << "MeV"; //Total Momentum of the proton beam
		LOG(debug) << "\nBeta (beam):\t" << (-bi) << "MeV"; //Beta of the proton beami
	}

	while (!evt) // eventloop
	{

		//------------ Internal momentum of a cluster -------------------
		TVector3 Pa;
		/*Pa.SetX(input_xy->GetRandom()); ### it is for momentum distribution from carlos calculations
		 *input_filename Pa.SetY(input_xy->GetRandom()); Pa.SetZ(input_z->GetRandom()); ### */
		Pa.SetX(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		Pa.SetY(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		Pa.SetZ(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		//------------ Internal momentum of the residual-----------------
		TVector3 PB;
		PB.SetX(-Pa.X());
		PB.SetY(-Pa.Y());
		PB.SetZ(-Pa.Z());


		// From the energy conservation in the virtual dissociation A->B+a
		double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());
		if (rrtt <= 0)
		{
			// cout<<"\nerror off-shell mass!!";//non-zero and real off-shell mass
			// cout<<"\nP:"<< PBx << "\t" << PBy << "\t" << PBz_rf << "\n";//non-zero and real off-shell mass
			continue;
		}

		// Off-shell mass of the bound cluster
		double Ma_off = sqrt(rrtt);

		// Total energies of "a" and "B" in the restframe of "A"
		double EaL = sqrt(Ma_off * Ma_off + Pa.Mag2());
		double EBL = sqrt(MB * MB + PB.Mag2());



		if(INVERSE)
		{
			//------- Lorentz transformations into laboratory system ---------
			std::pair<double, double> lora = R3BQFSGenerator::Lorentz(*gA,*bA, EaL, Pa.Z());
			EaL = lora.first; // cluster energy in lab
			Pa.SetZ(lora.second);    // cluster Pz in lab
		//	LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaLL); //<<<<<<<<<<<<<<<<<<<<<<<<<<		

			std::pair<double, double> lorB = R3BQFSGenerator::Lorentz(*gA, *bA, EBL, PB.Z());
			EBL = lorB.first; // energy of the residual B in lab
			PB.SetZ(lorB.second);    // Pz of the residual B in lab
			
		}

		/////////////////////////////////////////08.08.2024////////////////////////////////////////
		// Gamma radition from fragment
		if(AddGamma)
		{		
//			cout << "Gamma Energy Check:  " << GammaEnergy << endl;
			std::vector<Double_t> Gamma_vec = GetGamma();
			LVg = new TLorentzVector(Gamma_vec[0], Gamma_vec[1], Gamma_vec[2], Gamma_vec[3]);
		}
		/////////////////////////////////////////08.08.2024/////////////////////////////////////////

		//		std::pair<double, double> lorB = R3BQFSGenerator::Lorentz(gA, bA, EB, PB.Z());
		//		double EBL = lorB.first; // energy of the residual B in lab
		//		PB.SetZ(lorB.second);    // Pz of the residual B in lab
		
	//	Ea = EaL - Ma_off;
		LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaL); 
		//---------- Generating CM scattering process ----------
		//            double S = Ma_off * Ma_off + Mi * Mi + 2 * Mi * EaL; // Mandelstam invariant  //<<<<<<<<<<
		TLorentzVector LVStart = LVa + LVi;
		double S = LVStart.Mag2(); // Manndelstam invariant

		// Now generate CM scattering kinematics
		CM_values CM = R3BQFSGenerator::CENMASS(S, Ma_off, Mi, Ma, ISOTROPIC);
		if (!CM.good)
			continue; // non-physical output

		TVector3 P1cm(0., 0., 1.), P2cm(0., 0., 1.);
		double phi_rand = fRandom.Uniform(-PI, PI);

		P2cm.SetMag(CM.p_clust);
		P2cm.SetTheta(CM.theta_clust);
		P2cm.SetPhi(phi_rand);

		P1cm.SetX(-P2cm.X());
		P1cm.SetY(-P2cm.Y());
		P1cm.SetZ(-P2cm.Z());

		//------- Calculate realtive to the direction of the quasi-particle (cluster) --------
		//            double beta_cm = -Pa.Mag() / (EaL + Mi); // <<<<<<<<<<<<<
		double beta_cm = 0.00000-LVStart.Beta();
		double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

		std::pair<double, double> lora1 = R3BQFSGenerator::Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
		std::pair<double, double> lora2 = R3BQFSGenerator::Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

		P1cm.SetZ(lora1.second);
		P2cm.SetZ(lora2.second);
		//////////////////////////////Kai's Modification//////////////////////////                                                       
		//-------- Rotating back to the beam direction -----------
		//            TVector3 P1L = R3BQFSGenerator::DREHUNG(P1cm, Pa);
		//            TVector3 P2L = R3BQFSGenerator::DREHUNG(P2cm, Pa);
		//            
		///////////////////////////////////////////// 06.08.2024 /////////////////////////////////////
		
		TVector3 direction = LVStart.Vect();	
		direction = direction.Unit();
		P1cm.RotateUz(direction);
		P2cm.RotateUz(direction);

		TVector3 P1L = P1cm;
		TVector3 P2L = P2cm;

		///////////////////////////////////////////// 06.08.2024 /////////////////////////////////////
		
		double Angle_R   = gRandom->Uniform(0., 2.*TMath::Pi());
		vx = gRandom->Gaus(mean_vertex.X(), sigma_vertex.X());
		vy = gRandom->Gaus(mean_vertex.Y(), sigma_vertex.Y());
		vz = mean_vertex.Z() + gRandom->Uniform(-sigma_vertex.Z()/2.,sigma_vertex.Z()/2.);// 

		evt = true;

		LOG(debug) << "R3BQFSGenerator: Sending p2pevt: P1 : " << P1L.Px() << " , " << P1L.Py() << " , "
			<< P1L.Pz() << "\n P2 : " << P2L.Px() << " , " << P2L.Py() << " , " << P2L.Pz() << " ";
		R3BTrackingParticle f_p1(1,vx,vy,vz,P1L.Px()/1000.,P1L.Py()/1000.,P1L.Pz()/1000.,1,1);//P1L.P()/P1L.E(),P1L.M());
		for_ascii.push_back(f_p1);

		if(IS_PN)
		{
			R3BTrackingParticle f_n1(0,vx,vy,vz,P2L.Px()/1000.,P2L.Py()/1000.,P2L.Pz()/1000.,1,1);//P2L.P()/P2L.E(),P2L.M());
			for_ascii.push_back(f_n1);
		}
		else
		{
			R3BTrackingParticle f_p2(1,vx,vy,vz,P2L.Px()/1000.,P2L.Py()/1000.,P2L.Pz()/1000.,1,1);//P2L.P()/P2L.E(),P2L.M());
			for_ascii.push_back(f_p2);
		}	

		auto TotalEnergy = sqrt(PB.Mag2() + MB*MB);
		R3BTrackingParticle f_ion(5,vx,vy,vz,PB.Px()/1000.,PB.Py()/1000.,PB.Pz()/1000.,1,1);//1,PB.M());
		for_ascii.push_back(f_ion);

		if(AddGamma && LVg!= nullptr  )
		{
			TLorentzVector LVfragment;
			LVfragment.SetPxPyPzE(PB.Px(),PB.Py(), PB.Pz(), TotalEnergy);
			//	TLorentzVector gammaLabframe;
			SetLorentzBoost(*LVg, LVfragment);		

		}
		if(LVg != nullptr)
		{
			delete LVg;
			LVg = nullptr;
		}
		
		delete gA;
		delete bA;
	}
	return for_ascii;
}

Bool_t R3BQFSGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
	double E = ENERGY;
	// GenerateEvent(ENERGY,primGen);

	if (E <= 0)
	{
		LOG(error) << "R3BQFSGenerator: E < 0!";
		return kFALSE;
	}
	fRandom.SetSeed(0);
	//  	gRandom = new TRandom3(); //default
	//  	gRandom->SetSeed(0);//using computer time
	//  	TRandom2 r1;
	//  	r1.SetSeed(0);
	//R3BQFSGenerator::Print();
	bool evt = false;


	//-------------- The reaction 11C(p,pn)10C --------------------
	double Tkin = E * A_BEAM; // Total kinetic energy (MeV) of the projectile nucleus for Inverse Kinematics

	if(!INVERSE){Tkin = E;} //For Direct Kinematics

	const double sigma = MOM_SIGMA;
	TLorentzVector LVa; // Target
	TLorentzVector LVi; // Proton Beam 
	TLorentzVector* LVg = nullptr; // For Gamma

	double* gA = nullptr;
	double gi;
	double* bA = nullptr;
	double bi;
	double PA, Pi;
	double EA, Ei;
	double S_first;

	//Note: MA for inverse-kinematics beam, MB = fragment, Mi and Ma are proton or neutron
	if(INVERSE) //Beam Parameters for Inverse Kinematics
	{
		//--------------------- Beam parameters -----------------------------
		PA = sqrt(Tkin * (Tkin + 2 * MA));         // Total 3-momentum of the beam
		EA = sqrt(MA * MA + PA * PA);              // Total energy of the beam
		bA = new double(-PA / EA);                             // Beta of the beam
		gA = new double(1/TMath::Sqrt(1-*bA* *bA));                // Gamma of the beam
		S_first = (EA + Mi) * (EA + Mi) - PA * PA; // Invariant mass (Mandelstam S-variable)
		//	TLorentzVector LVa; // Proton Target
		//	TLorentzVector LVi; // Beam
		LVi.SetPxPyPzE(0., 0.,0.,Mi); //Lorentz Vector of the Beam
		//	TLorentzVector* LVg = nullptr; // For Gamma

		LOG(debug) << "\n****** Beam Parameters for Inverse Kinematics ********";
		LOG(debug) << "\nMA:\t" << MA << " MeV";
		LOG(debug) << "\nTotal momentum:\t" << PA << " MeV";
		LOG(debug) << "\nTotal energy:\t" << EA << " MeV";
		LOG(debug) << "\nBeta (beam):\t" << (-*bA) << "\nGamma (beam):\t" << *gA << "\n\n";
	}

	else //For Direct Kinematics
	{
		//Using the same variables' names here, so the following code is shorter
		Pi = sqrt(Tkin*(Tkin+2*Mi)); // Total 3 momentum of the proton beam
		Ei = sqrt(Mi*Mi + Pi*Pi);    // Total energy of the proton beam
		bi = -Pi/Ei;		  // Beta of the proton beam
		gi = 1/sqrt(1-bi*bi);	  // Gamma of the proton beam
		S_first = (Ei+MA)*(Ei+MA)-Pi*Pi; //Madelstam S-invariant of the proton beam

		//	TLorentzVector LVa; // Target
		//	TLorentzVector LVi; // Proton Beam 
		LVi.SetPxPyPzE(0.,0.,Pi, Ei); //of the Bram

		//	TLorentzVector* LVg = nullptr; //For Gamma Radition from the fragment

		LOG(debug) << "\n******Beam Parameters for Direct Kinematics******"; 
		LOG(debug) << "\nMi: \t" << Mi << "MeV"; //Mass of the proton beam
		LOG(debug) << "\nMA: \t" << MA << "MeV"; //Mass of A target
		LOG(debug) << "\nTotal Momentum:\t" << Pi << "MeV"; //Total Momentum of the proton beam
		LOG(debug) << "\nBeta (beam):\t" << (-bi) << "MeV"; //Beta of the proton beami
	}

	while (!evt) // eventloop
	{

		//------------ Internal momentum of a cluster -------------------
		TVector3 Pa;
		/*Pa.SetX(input_xy->GetRandom()); ### it is for momentum distribution from carlos calculations
		 *input_filename Pa.SetY(input_xy->GetRandom()); Pa.SetZ(input_z->GetRandom()); ### */
		Pa.SetX(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		Pa.SetY(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		Pa.SetZ(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
		//------------ Internal momentum of the residual-----------------
		TVector3 PB;
		PB.SetX(-Pa.X());
		PB.SetY(-Pa.Y());
		PB.SetZ(-Pa.Z());


		// From the energy conservation in the virtual dissociation A->B+a
		double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());

		if (rrtt <= 0)
		{
			// cout<<"\nerror off-shell mass!!";//non-zero and real off-shell mass
			// cout<<"\nP:"<< PBx << "\t" << PBy << "\t" << PBz_rf << "\n";//non-zero and real off-shell mass
			continue;
		}

		// Off-shell mass of the bound cluster
		double Ma_off = sqrt(rrtt);

		// Total energies of "a" and "B" in the restframe of "A"
		double EaL = sqrt(Ma_off * Ma_off + Pa.Mag2());
		double EBL = sqrt(MB * MB + PB.Mag2());



		if(INVERSE)
		{
			//------- Lorentz transformations into laboratory system ---------
			std::pair<double, double> lora = R3BQFSGenerator::Lorentz(*gA,*bA, EaL, Pa.Z());
			EaL = lora.first; // cluster energy in lab
			Pa.SetZ(lora.second);    // cluster Pz in lab
		//	LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaLL); //<<<<<<<<<<<<<<<<<<<<<<<<<<		

			std::pair<double, double> lorB = R3BQFSGenerator::Lorentz(*gA, *bA, EBL, PB.Z());
			EBL = lorB.first; // energy of the residual B in lab
			PB.SetZ(lorB.second);    // Pz of the residual B in lab
			
		}

		/////////////////////////////////////////08.08.2024////////////////////////////////////////
		// Gamma radition from fragment
		if(AddGamma)
		{		
//			cout << "Gamma Energy Check:  " << GammaEnergy << endl;
			std::vector<Double_t> Gamma_vec = GetGamma();
			LVg = new TLorentzVector(Gamma_vec[0], Gamma_vec[1], Gamma_vec[2], Gamma_vec[3]);
		}
		/////////////////////////////////////////08.08.2024/////////////////////////////////////////

		//		std::pair<double, double> lorB = R3BQFSGenerator::Lorentz(gA, bA, EB, PB.Z());
		//		double EBL = lorB.first; // energy of the residual B in lab
		//		PB.SetZ(lorB.second);    // Pz of the residual B in lab
		
	//	Ea = EaL - Ma_off;
		LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaL); 
		//---------- Generating CM scattering process ----------
		//            double S = Ma_off * Ma_off + Mi * Mi + 2 * Mi * EaL; // Mandelstam invariant  //<<<<<<<<<<
		TLorentzVector LVStart = LVa + LVi;
		double S = LVStart.Mag2(); // Manndelstam invariant

		// Now generate CM scattering kinematics
		CM_values CM = R3BQFSGenerator::CENMASS(S, Ma_off, Mi, Ma, ISOTROPIC);
		if (!CM.good)
			continue; // non-physical output

		TVector3 P1cm(0., 0., 1.), P2cm(0., 0., 1.);
		double phi_rand = fRandom.Uniform(-PI, PI);

		P2cm.SetMag(CM.p_clust);
		P2cm.SetTheta(CM.theta_clust);
		P2cm.SetPhi(phi_rand);

		P1cm.SetX(-P2cm.X());
		P1cm.SetY(-P2cm.Y());
		P1cm.SetZ(-P2cm.Z());

		//------- Calculate realtive to the direction of the quasi-particle (cluster) --------
		//            double beta_cm = -Pa.Mag() / (EaL + Mi); // <<<<<<<<<<<<<
		double beta_cm = 0.00000-LVStart.Beta();
		double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

		std::pair<double, double> lora1 = R3BQFSGenerator::Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
		std::pair<double, double> lora2 = R3BQFSGenerator::Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

		P1cm.SetZ(lora1.second);
		P2cm.SetZ(lora2.second);
		//////////////////////////////Kai's Modification//////////////////////////                                                       
		//-------- Rotating back to the beam direction -----------
		//            TVector3 P1L = R3BQFSGenerator::DREHUNG(P1cm, Pa);
		//            TVector3 P2L = R3BQFSGenerator::DREHUNG(P2cm, Pa);
		//            
		///////////////////////////////////////////// 06.08.2024 /////////////////////////////////////
		TVector3 direction = LVStart.Vect();	
		direction = direction.Unit();
		P1cm.RotateUz(direction);
		P2cm.RotateUz(direction);

		TVector3 P1L = P1cm;
		TVector3 P2L = P2cm;

		///////////////////////////////////////////// 06.08.2024 /////////////////////////////////////
		
		double Angle_R   = gRandom->Uniform(0., 2.*TMath::Pi());
		vx = gRandom->Gaus(mean_vertex.X(), sigma_vertex.X());
		vy = gRandom->Gaus(mean_vertex.Y(), sigma_vertex.Y());
		vz = mean_vertex.Z() + gRandom->Uniform(-sigma_vertex.Z()/2.,sigma_vertex.Z()/2.);// 

		evt = true;
		LOG(debug) << "R3BQFSGenerator: Sending p2pevt: P1 : " << P1L.Px() << " , " << P1L.Py() << " , "
			<< P1L.Pz() << "\n P2 : " << P2L.Px() << " , " << P2L.Py() << " , " << P2L.Pz() << " ";
		primGen->AddTrack(2212, P1L.Px() / 1000., P1L.Py() / 1000., P1L.Pz() / 1000., vx, vy, vz);
		if(IS_PN)
		{
		primGen->AddTrack(2112, P2L.Px() / 1000., P2L.Py() / 1000., P2L.Pz() / 1000., vx, vy, vz); //for (p,np)
		}
		else
		{
		primGen->AddTrack(2212, P2L.Px()/1000., P2L.Py()/1000., P2L.Pz()/1000., vx, vy, vz); //for (p,2p)
		}	
	
		TParticlePDG* thisPart = TDatabasePDG::Instance()->GetParticle(fIonFragment->GetName());
		if (!thisPart)
			LOG(fatal) << "FairIonGenerator: Ion " << fIonFragment->GetName() << " not found in database!";
		fIonFragmentPDG = thisPart->PdgCode();
		//cout << "\nFrgagment PDG:" << fIonFragmentPDG << "\n";


		//////////////////////////////Kai's Modification//////////////////////////

		auto TotalEnergy = sqrt(PB.Mag2() + MB*MB);
		//		primGen->AddTrack(fIonFragmentPDG, PB.Px()/1000., PB.Py()/1000., PB.Pz()/1000.,0.,0.,0.,-1,true, TotalEnergy/1000.);
		//primGen->AddTrack(fIonFragmentPDG, PB.Px()/1000., PB.Py()/1000., PB.Pz()/1000.,0.,0.,0.,-1,true, TotalEnergy/1000.);
		primGen->AddTrack(fIonFragmentPDG, PB.Px()/1000., PB.Py()/1000., PB.Pz()/1000.,vx, vy, vz);

		if(AddGamma && LVg!= nullptr  )
		{
			TLorentzVector LVfragment;
			LVfragment.SetPxPyPzE(PB.Px(),PB.Py(), PB.Pz(), TotalEnergy);
		//	TLorentzVector gammaLabframe;
			SetLorentzBoost(*LVg, LVfragment);		
			
			//	primGen->AddTrack(22,Gamma_vec[0], Gamma_vec[1], Gamma_vec[2], 0., 0., 0., -1, true, Gamma_vec[3]); 
			//primGen->AddTrack(22, LVg->Px()/1000., LVg->Py()/1000., LVg->Pz()/1000.,0.,0.,0., -1, true, LVg->E());
			primGen->AddTrack(22, LVg->Px()/1000., LVg->Py()/1000., LVg->Pz()/1000.,vx, vy, vz);
		}

		//Clearing Memories
		//	std::vector<Double_t> Gamma_vec = SetGamma(gammaenergy);
		if(LVg != nullptr)
		{
			delete LVg;
			LVg = nullptr;
		}
		
		delete gA;
		delete bA;

		//////////////////////////////Kai's Modification//////////////////////////
		//--------- Filling in the ROOTTree variables------------
	}


	return kTRUE;
}

void R3BQFSGenerator::SetValues(double E, int A, double MOM, double exe, double gammaenergy , bool invert, bool iso, bool gamma)
{
	ENERGY = E;
	MOM_SIGMA = MOM;
	A_BEAM = A;
	double UNIT = 931.494061;
//	MA = A * UNIT;
	MB =(A - 1) * UNIT + exe;
//	MB = (A - 1) * UNIT;
//	IS_PN = is_pn;
	if(IS_PN)
	{
		Ma = 939.565;
	}
	else
	{
		Ma = 938.727;
	}
	Mi = 938.727;
//	T_LIMIT = 10000000.;
	INVERSE = invert;
	ISOTROPIC = iso;
	AddGamma = gamma;
//	IS_PN = is_pn;
	GammaEnergy = gammaenergy;
	mean_vertex.SetXYZ(0,0,0);
	sigma_vertex.SetXYZ(0,0,0);
	return;
}

TVector3 R3BQFSGenerator::DREHUNG(TVector3 v1, TVector3 v2)
{
	double CT = v2.Z() / v2.Mag(); // cos(theta) of v2 wrt. Z-axis
	double ST = sqrt(1 - CT * CT); // sin(theta)
	double CF = v2.X() / v2.Mag() / ST;
	double SF = v2.Y() / v2.Mag() / ST;

	TVector3 v3;
	double _v3x = v1.X() * CT * CF - v1.Y() * SF + v1.Z() * ST * CF;
	double _v3y = v1.X() * CT * SF + v1.Y() * CF + v1.Z() * ST * SF;
	double _v3z = -v1.X() * ST + v1.Z() * CT;
	v3.SetXYZ(_v3x, _v3y, _v3z);
	return v3;
}

// Kinematical function
double R3BQFSGenerator::CINEMA(double x, double y, double z)
{
	double lambda = x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
	return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
CM_values R3BQFSGenerator::CENMASS(double s, double m2off, double m1, double m2, bool isotropic)
{
	CM_values output;
	output.good = false;
	double X = s;
	double Y = m2off * m2off;
	double Z = m1 * m1;
	double sqrs = sqrt(s);

	// Kinematics before the scattering process
	// (with one off-shell mass)
	double p2_off = sqrt(R3BQFSGenerator::CINEMA(X, Y, Z)) / 2 / sqrs;
	double p1_off = p2_off;
	// CM energies
	double e1_off = (s + Z - Y) / 2 / sqrs;
	double e2_off = (s + Y - Z) / 2 / sqrs;

	// Now take the real masses (after scattering)
	Y = m2 * m2;
	Z = m1 * m1;
	// And check whether the kinematical function is ok
	// for this specific kinematical case
	double error_CI = R3BQFSGenerator::CINEMA(X, Y, Z);
	if (error_CI <= 0.)
	{
		// cout << "\nerror!!! Kinematical function is negative!";
		return output;
	}

	// Kinematics after the scattering process
	// (with all real masses)
	double p2 = sqrt(R3BQFSGenerator::CINEMA(X, Y, Z)) / 2 / sqrs;
	double p1 = p2;
	double e1 = (s + Z - Y) / 2 / sqrs;
	double e2 = (s + Y - Z) / 2 / sqrs;

	// Let's consider momentum transfer <t> from the
	// target particle 1 to the cluster 2
	double tmax = 2 * (m1 * m1 - e1_off * e1 - p1_off * p1); // COSINE=(-1)
	double tmin = 2 * (m1 * m1 - e1_off * e1 + p1_off * p1); // COSINE=(1)
	// cout << "\n\n Tmax = " << tmax;
	// cout << "\n Tmin = " << tmin;
	// cout << "\n Mandels = " << X;

	double t;
	// Generate random momentum transfer for this kinematical case
	if (!isotropic)
	{
		if(!IS_PN)
		{
			t = R3BQFSGenerator::get_T_Cugnon(s, tmax);
		}
		else
		{
			t = R3BQFSGenerator::get_T_Cugnon_PN(s, tmax);
		}
	} // Using parameterized cross sections
	else
	{
		t = fRandom.Uniform(tmax, 0);
	} // Isotropic scattering

	// double COSINE = (t - m2off*m2off - m2*m2 + 2*e2_off*e2)/(2*p2_off*p2);
	double COSINE = (t - 2 * m1 * m1 + 2 * e1_off * e1) / (2 * p1_off * p1);
	if (fabs(COSINE) >= 1)
	{ // momentum transfer out of range
		// cout << "\nerror! Scattering cosine is larger than 1";
		return output;
	}

	// CM scattering angles
	double theta1 = acos(COSINE);
	double theta2 = PI - theta1;

	output.e_clust = e2;
	output.p_clust = p2;
	output.theta_clust = theta2;

	output.e_scat = e1;
	output.p_scat = p1;
	output.theta_scat = theta1;

	output.T = t;
	output.good = true;

	return output;
}

// Calculate 3-momentum in CM system of two particles M1 and M2
// when M1 has kinetic energy TLAB and M2 is at rest
double R3BQFSGenerator::momentum_CM(double TLAB, double M1, double M2)
{
	// Particle M2 is assumed to be in rest
	double PLAB = sqrt(TLAB * (TLAB + 2 * M1));      //  Total 3-momentum of an incident particle in Lab
	double ELAB = sqrt(PLAB * PLAB + M1 * M1);       //  Total energy of an incident particle in lab
	double SLAB = M1 * M1 + M2 * M2 + 2 * M2 * ELAB; // Mandelstam invariant S in lab
	double PCM = PLAB * M2 / sqrt(SLAB);             // Momentum of both particles in CM frame
	return PCM;
}

std::pair<double, double> R3BQFSGenerator::Lorentz(double g, double b, double e, double p)
{
	double eL = g * e - g * b * p;
	double pL = g * p - g * b * e;
	return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c)� units)
// distributed according to the parameterized proton-proton
// invarant cross section. Pass as a parameter "sm" the
// Mandelstam variable S (in MeV�)
// and the maximum possible momentum transfer
double R3BQFSGenerator::get_T(double sm, double max)
{
	// TRandom1 rand;
	// rand.SetSeed(0);

	double Tmax = max * 0.000001; // convert to GeV� units
	// double Tmin = min*0.000001;

	// double Tmax = -2*pCM*pCM*(1 - cos(PI))*0.000001; //in (GeV/c)�
	Double_t rr = fRandom.Uniform(-1., 1.); // to randomize wrt 90 degrees
	double mandels = sm * 0.000001;         // in GeV�
	// cout << "\nMandelstam S = " << mandels << "\t Tmax/2 = " << Tmax/2 << "\t Random: " << rr;

	// Probability function from the parameterization
	TF1* foo = new TF1("foo", "[0]*exp(x*[1])*(1+0.02*exp((-6)*x))", Tmax / 2, 0);
	double c = 0.;
	if (mandels <= 4.79)
		c = -3283.75 + 3064.11 * mandels - 1068.44 * mandels * mandels + 164.844 * pow(mandels, 3) -
			9.48152 * pow(mandels, 4);
	else if (mandels > 4.79)
		c = -776.822 + 586.016 * mandels - 175.347 * mandels * mandels + 26.1823 * pow(mandels, 3) -
			1.94889 * pow(mandels, 4) + 0.0578352 * pow(mandels, 5);

	foo->FixParameter(0, 25.); // normalization constant (could be anything)
	foo->FixParameter(1, c);

	double Trand = foo->GetRandom(Tmax / 2, 0.); // from 90 to 0 degrees
	if (rr > 0)

		Trand = Tmax - Trand; // symmetrization relative to 90 degrees
	// cout << "\n Tmax/2 = " << Tmax/2 *  1000000;
	// cout << "\n Random T = " << Trand*1000000;
	delete foo;
	return (Trand * 1000000); // returning value in MeV�
}

double R3BQFSGenerator::get_T_Cugnon(double sm, double max)
{
	double Tmax = max * 0.000001; // convert to GeV� units
	// double Tmin = min*0.000001;


//	Double_t rr = fRandom.Uniform(-1., 1.); // to randomize wrt 90 degrees
	double mandels = sm * 0.000001;         // in GeV�
	//cout << "\nMandelstam S = " << mandels << "\t Tmax= " << Tmax << endl;
	Double_t Mp = 938.727/1000.; //Converting to GeV to fit in Cugnon Parameterization
	//cout << "\n Proton Mass:%f GeV/c^2" << Mp<< endl;
	// Probability function from the parameterization
	Double_t Bpp;
	Double_t p_lab = TMath::Sqrt(mandels*(mandels-4*Mp*Mp))/(2*Mp); //Convert p_lab into Mandelstam mandels
	//cout << "\n Plab::" << p_lab;	
	double pp = pow(p_lab,8);
	if(p_lab <= 2)
	{
		Bpp =(5.5*pp)/(7.7+pp);
	}
	else
	{
		Bpp = 5.334+0.67*(p_lab-2);
	}
	//cout << "\nBpp: " << Bpp << endl;
	TF1 *foo =  new TF1("foo", "exp(x*[0])",Tmax,0);	
	foo->FixParameter(0,Bpp);
//	cout << "\n foo Value: " << foo << endl;
	double Trand = foo->GetRandom(Tmax,0);//This is the Inverse Transform sampling method, we will random t that satisfies the above expression.
	//Condition Checking
	//cout << "\n Check Trandom" << Trand << endl;	
	//cout << "\nDifferential Cross Section Check: " << TMath::Exp(Trand*Bpp) << endl;
	return (Trand * 1000000); // returning value in MeV�
//	return std::make_pair(t_value* 1000000, TMath::Exp(Bpp*t_value) );
}

double R3BQFSGenerator::get_T_Cugnon_PN(double sm, double max)
{
	double Tmax = max * 0.000001;
	double mandels = sm*0.000001;
	double Mp = 938.727/1000.;
	double Mn = 939.565/1000.;
	//cout << "\n Proton Mass: " << Mp << "\t Neutron Mass: " << Mn;
	double Bpp=0;
	double Bnp=0;
	double a = 0;
	
	double_t Plab = TMath::Sqrt((mandels-pow(Mp-Mn,2))*(mandels-pow(Mp+Mn,2)))/(2*Mp);

	if(Plab<2)
	{
		Bpp = 5.5*pow(Plab,8)/(7.7 + pow(Plab,8));
	}
	else
	{
		Bpp = 5.334+0.67*(Plab-2);
	}

	if(Plab<0.225)
	{
		Bnp = 0;
	}
	if(0.225<=Plab && Plab<0.6)
	{
		Bnp = 16.53*(Plab-0.225);
	}
	if(0.6<= Plab && Plab <1.6)
	{
		Bnp = -1.63*Plab+7.16;
	}
	else
	{
		Bnp = Bpp;
	}

	if(Plab>0.8)
	{
		a = 0.64/pow(Plab,2);
	}
	else
	{
		a = 1;
	}


	TF1 * foo = new TF1("foo","exp(-[0]*x )+[1]*exp([0]*([3]-[2]+x))", 0, Tmax);

	foo->FixParameter(0,Bnp);
	foo->FixParameter(1,a);
	foo->FixParameter(2,mandels);
	foo->FixParameter(3,2*(Mn*Mn+Mp*Mp));
	double Trand = foo->GetRandom(Tmax,0);
	return (Trand*1000000); //return the value in MeV

}

ClassImp(R3BQFSGenerator)
