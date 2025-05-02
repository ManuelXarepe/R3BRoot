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

/*
 * R3BQFSGenerator.h
 *
 *  Created on: Jul 23, 2018
 *      Author: gu92joq
 */

#ifndef R3BQFSGENERATOR_H_
#define R3BQFSGENERATOR_H_

#include <map>
#include <utility>

#include <TRandom3.h>
#include <TVector3.h>

#include <FairGenerator.h>
#include <TLorentzVector.h>
#include "R3BTrackingParticle.h"
/////////////////Kai's Modification//////////////////
#include <FairIon.h>


/////////////////Kai's Modification//////////////////


/**
 * Event generator for p2p events in inverse and direct kinematics, alpha version
 * Supply energy, nucleus A and excitation energy by separate method SetValues
 */

struct CM_values
{
    // Internal cluster
    double e_clust;
    double p_clust;
    double theta_clust;
    // Scattered particle
    double e_scat;
    double p_scat;
    double theta_scat;
    // indicates satisfactory kinematics (i.e. energy & momentum conservation)
    bool good;
    double T;
};

class R3BQFSGenerator : public FairGenerator
{
  public:
    // R3BQFSGenerator() : FairGenerator("R3BQFSGenerator", "R3BQFSGenerator") {}
    // R3BQFSGenerator(const char *name, const char *title) : FairGenerator(name, title) {}
    R3BQFSGenerator();
    ~R3BQFSGenerator(){};

    virtual Bool_t Init();
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);
    void SetProjectile(int A, int Z, double E); // Set A, Z, energy per nucleon
    void SetHeavyFragment(int A, int Z); // Set A, Z
    void SetInverse(bool inverse);    // Set Inverse or direct kinematics
    void SetIsotropic(bool isotropic);          // Set Isotropic or parametric cross section for proton emission
    void SetMomDistrib(double mom);             // Set Sigma of momentum distribution
    void SetExcitation(double exe);             // Set Excitation energy
    void SetLightNucleus(double ma, double mi,bool ); // Set Incoming outgoing particle mass
    void SetAddGamma(double energy);
    void SetVertexDistrib(TVector3 mean, TVector3 sigma);
    std::vector<Double_t> GetGamma();
    void Print();
    void SetLorentzBoost(TLorentzVector& gammaMomentum, const TLorentzVector& fragmentMomentum);
   // double GetSM() const{return S;}		//Extract the values of Mandelstam variable to plot the parameterization
    double get_T_Cugnon(double, double);
    double get_T_Cugnon_PN(double, double);
    std::vector <R3BTrackingParticle> get_ascii_input();
    //void get_ascii_input(std::vector <R3BTrackingParticle> *for_ascii);
  protected:
    /**
     * Set Values for target end beam nucleus, Cross section and kinematics
     */
    double vx = 0;
    double vy = 0;
    double vz = 0;
    TVector3 mean_vertex;
    TVector3 sigma_vertex;
    double ENERGY;
    int A_BEAM;
    double MA;
    double MB;
    double Ma;
    double GammaEnergy;
    double Mi;
 //   double T_LIMIT;
    double MOM_SIGMA;
    bool ISOTROPIC;
    bool INVERSE; 
    bool AddGamma;
    bool IS_PN;
    TRandom3 fRandom;
    double PI = 3.14159265358979323846;
    TVector3 DREHUNG(TVector3 v1, TVector3 v2);
    double CINEMA(double, double, double);
    CM_values CENMASS(double, double, double, double, bool);
    double momentum_CM(double, double, double);
    double get_T(double, double);
   // double get_T_Cugnon(double, double);
    std::pair<double, double> Lorentz(double, double, double, double);
    void SetValues(double E, int A, double MOM, double exe, double gammaenergy, bool invert, bool iso, bool gammas_pn);
    FairIon* fIonBeam; // Pointer to the FairIon
    FairIon* fIonFragment; // Pointer to the FairIon
    Int_t fIonBeamPDG;
    Int_t fIonFragmentPDG;
   // TF1 * foo_Cugnon_DXS;//dSigma/dt parameterization from Cugnon
	
 private:
    //double S; //Mandelstam variable    
    ClassDef(R3BQFSGenerator, 1);
};

#endif
