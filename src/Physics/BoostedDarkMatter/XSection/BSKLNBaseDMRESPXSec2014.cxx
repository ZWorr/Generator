//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Cross Section Calculation

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________
#include <TMath.h>
#include <TSystem.h>

//Check what these are: update "Interactions" for DM
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"//
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"//
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

//_____DM change
#include "Physics/BoostedDarkMatter/XSection/BSKLNBaseDMRESPXSec2014.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMI.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplDM.h"
#include "Framework/ParticleData/PDGLibrary.h"
//______


using namespace genie;
using namespace genie::constants;

//DM updated code
//____________________________________________________________________________
BSKLNBaseDMRESPXSec2014::BSKLNBaseDMRESPXSec2014(string name) :
XSecAlgorithmI(name)
{

}
//____________________________________________________________________________
BSKLNBaseDMRESPXSec2014::BSKLNBaseDMRESPXSec2014(string name, string config) :
XSecAlgorithmI(name, config)
{

}
//____________________________________________________________________________
BSKLNBaseDMRESPXSec2014::~BSKLNBaseDMRESPXSec2014()
{

}


//_________________________________________________________________________

double BSKLNBaseDMRESPXSec2014::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction -> Kine();
  const Target & target = init_state.Tgt();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "Using v^" << fVelMode << " dependence";

//____________________________________________________________________________
//Get Kinematic Parameters:

  //Invarient mass W
  double W  = kinematics.W();
  //q^(2) momentum transfer
  double q2 = kinematics.q2();
  //DM mass: mx
  double mchi   = init_state.GetProbeP4(kRfHitNucRest)->M();
  //E1:
  double E      = init_state.ProbeE(kRfHitNucRest);
  //m:
  double Mnuc   = target.HitNucMass();

//Get Resonance Parameters:

  // Get the input baryon resonance
    Resonance_t resonance = interaction->ExclTag().Resonance();
    string      resname   = utils::res::AsString(resonance);
    bool        is_delta  = utils::res::IsDelta (resonance); //?? don't use?


  // Get the Dark Matter, hit nucleon & DM current
      int  nucpdgc   = target.HitNucPdg();
      int  probepdgc = init_state.ProbePdg();
      bool is_dm     = pdg::IsDarkMatter        (probepdgc);
      bool is_dmbar  = pdg::IsAntiDarkMatter    (probepdgc);
      bool is_p      = pdg::IsProton  (nucpdgc);
      bool is_n      = pdg::IsNeutron (nucpdgc);
      bool is_DM     = proc_info.IsDarkMatterResonant();

  //Need break for is not DM
    if(!is_DM) {
    return 0;
    }
//____________________________________________________________________________


    //Auxillary Kinematic factors
    double W2     = TMath::Power(W,    2);       //W^2
    double Mnuc2  = TMath::Power(Mnuc, 2);       //m^2
    double sqrtq2 = TMath::Sqrt(-q2);            //Sqrt(-q^2)
    double mchi2 = TMath::Power(mchi, 2);        // mx^2
    double mZprime2 = TMath::Power(fMedMass, 2); // mZ'^2

  //Refrence Frame changes
    //Define K:
    //K = (W^2 - m^2) / 2m
    double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
    //Solve for K in RRF; define as K^* and form relation between Lab and RRF
    //K^* = K * (m/W)
    double kstar  = k * (Mnuc/W);
    //Solve for nu and nu^*
    //(W^2 - m^2) / 2m = nu + (q^2 / 2m)
    double v      = k - 0.5 * (q2/Mnuc);
    //(W^2 - m^2) / 2M = nu^* - (q^2 / 2W)
    double vstar  = kstar + 0.5 * (q2/W);

    //nu^2
    double v2     = TMath::Power(v, 2);
    //nu^(*2)
    double v2star = TMath::Power(vstar,2);
  //Use Invarient: q^2 = nu^2 - Q^2 = nu^(*2) - Q^(*2) to solve for Q and Q^*
    //Q^2 = nu^2 - q^2
    double Q2     = v2 - q2;
    //Q^(*2)
    double Q2star = v2star - q2;
    //Q
    double Q      = TMath::Sqrt(Q2);
    //Q^(*)
    double Qstar  = TMath::Sqrt(Q2star);
    //Use nu to solve for E'
    //nu = E - E'
    double Eprime = E - v;

  //Terms useful in both scalar and fermion DM XSec
    //Common Factor:  (1 / (mZ'^2 + Q^2)^2)*(-q^2/Q^2)
    double sigcommon = (1/TMath::Power((mZprime2 + Q2),2)) * (-q2/Q2);

    double E12 = TMath::Power((E + Eprime), 2);          // (E1 + E')^2
    double Xminus = (4.*mchi2/q2) - 1;                   // (2*mx)^2 / q^2 - 1
    double Xplus = (4.*mchi2/q2) + 1;                    // (2*mx)^2 / q^2 + 1
    double denom = 1/(4. * (TMath::Power(E, 2) - mchi2));// 1 / 4*(E^2 - mx^2)

    double E12_Q2 = (E12 - Q2)*denom;
    double E12_Q2_V = (E12 + Q2*Xplus)*denom;
    double E12_Q2_A = (E12 - Q2*Xminus)*denom;
    double E12_Q2_SA = (E12 + Q2*Xminus)*denom;
    double VXA = (4.*Q*(E + Eprime))*denom;

//____________________________________________________________________________

//____________________________________________________________________________
//Same for DM:
//Dealing with the overlap of DIS and RES regimes.

// Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
if(fUsingDisResJoin) {
  if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("BSKLNBaseDMRESPXSec2014", pDEBUG)
      << "RES/DIS Join Scheme: XSec[RES, W=" << W
      << " >= Wcut=" << fWcut << "] = 0";
#endif
    return 0;
  }
}
//____________________________________________________________________________
// Get baryon resonance parameters
  //N:
  int    IR  = utils::res::ResonanceIndex    (resonance);
  //L:
  int    LR  = utils::res::OrbitalAngularMom (resonance);
  //M
  double MR  = utils::res::Mass              (resonance);
  //Width:
  double WR  = utils::res::Width             (resonance);
   double NR  = fNormBW?utils::res::BWNorm    (resonance,fN0ResMaxNWidths,fN2ResMaxNWidths,fGnResMaxNWidths):1;

//Breit-Wigner:
  // Following NeuGEN, avoid problems with underlying unphysical
  // model assumptions by restricting the allowed W phase space
  // around the resonance peak
 if (fNormBW) {
        if      (W > MR + fN0ResMaxNWidths * WR && IR==0) return 0.;
        else if (W > MR + fN2ResMaxNWidths * WR && IR==2) return 0.;
        else if (W > MR + fGnResMaxNWidths * WR)          return 0.;
  }


//_________________________________________________________________________
//Common Adjustment Factor in RS: G[V,A](q^2) = (1 - q^2/4m^2)^(1/2-N) (  1/ (1 - q^2/m[V,A]^2)  )^2
  // Calculate the Feynman-Kislinger-Ravndall parameters

  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);

  if(fGV){

    LOG("BSKLNBaseDMRESPXSec2014",pDEBUG) <<"Using new GV";
    double CV0 =  1./(1-q2/fMv2/4.);
    double CV3 =  2.13 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV4 = -1.51 * CV0 * TMath::Power( 1-q2/fMv2,-2);
    double CV5 =  0.48 * CV0 * TMath::Power( 1-q2/fMv2/0.766, -2);

    double GV3 =  0.5 / TMath::Sqrt(3) * ( CV3 * (W + Mnuc)/Mnuc
                  + CV4 * (W2 + q2 -Mnuc2)/2./Mnuc2
                  + CV5 * (W2 - q2 -Mnuc2)/2./Mnuc2 );

    double GV1 = - 0.5 / TMath::Sqrt(3) * ( CV3 * (Mnuc2 -q2 +Mnuc*W)/W/Mnuc
                 + CV4 * (W2 +q2 - Mnuc2)/2./Mnuc2
                 + CV5 * (W2 -q2 - Mnuc2)/2./Mnuc2 );

    GV = 0.5 * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR)
         * TMath::Sqrt( 3 * GV3*GV3 + GV1*GV1);
  }

  if(fGA){
    LOG("BSKLNBaseDMRESPXSec2014",pDEBUG) << "Using new GA";

    double CA5_0 = 1.2;
    double CA5 = CA5_0 *  TMath::Power( 1./(1-q2/fMa2), 2);
    //  GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5/fZeta;
    GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5;

    LOG("BSKLNBaseDMRESPXSec2014",pINFO) <<"GA= " <<GA << "  C5A= " <<CA5;
  }
  //JN end of new form factors code

//____________________________________________________________________________
//Define FKR kinematic factors:

  //2mg^2 = ((W+m)^2 - q^2)/2M
  double d      = (TMath::Power(W+Mnuc,2.) - q2)/(2. * W);
  //useful note: 2Wd = 4mMg^2 = (W+m)^2 - q^2

  // (2/Omega)^(1/2)
  double sq2omg = TMath::Sqrt(2./fOmega);
  //useful note: 1/sq2omg = Sqrt(Omega/2)

  //N*Omega
  double nomg   = IR * fOmega;

//__________Adjusted for DM____________
//Lamda = (2/Omega)^(1/2) * Q^(*)
  fFKRDM.Lamda  = sq2omg * Qstar;

//T^(V) = [ GV / (3*M) ] * Sqrt( Omega/2 )
  fFKRDM.Tv     = GV / (3.*W*sq2omg); //Same for DM

//T^(A) = [ (GA*Z) / (3*M) ] * Sqrt( Omega/2 ) * [ Q^(*) / 2mg^2 ]
  fFKRDM.Ta     = ((GA * fZeta)/(3. * W)) * (1./sq2omg) * (Qstar / d);

//R^(V) = Sqrt(2) * GV * Q^(*) * [ (M + m) / ((M + m)^2 - q^2 ) ]
  fFKRDM.Rv     = kSqrt2 * GV * Qstar * ((W + Mnuc)/(2. * W * d));

//R^(A) = [ (Sqrt(2) * GA * Z / (6*M) ] * (M + m + ((N * Omega) / 2mg^2) )
  fFKRDM.Ra     = ((kSqrt2 * GA * fZeta)/(6. * W)) * (W + Mnuc + (nomg/d));

//S = [(-q^2) / Q^(*2)] * [ GV / 6*M^2] * (3Mm + q^2 - m^2)
  fFKRDM.S      = (-q2/Q2star) * (GV / (6. * W2)) * (3. * W * Mnuc + q2 - Mnuc2);

//B_s = [(Z*GA) / (3*M)] * Sqrt(Omega/2) * (1 + (nu^(*) / 2mg^2))
  fFKRDM.Bs     = ( (fZeta * GA)/(3. * W * sq2omg) ) * (1 + (vstar / d));

//C_s = [(Z*GA) / (6*M*Q^(*))] * (M^2 - m^2 + N*Omega*( nu^(*) / 2mg^2))
  fFKRDM.Cs     = ( (fZeta * GA)/(6. * W * Qstar) ) * (W2 - Mnuc2 + nomg*(vstar / d));

//B_z = [(Z*GA) / (3*M)] * Sqrt(Omega/2) * (1/Q^(*)) * (M - m)
  fFKRDM.Bz     = ( (fZeta * GA)/(3. * W * Qstar * sq2omg) ) * (W - Mnuc);

//C_z =  [(Z*GA) / (6*M)] * (2*M + ( (3*q^2) / 2mg^2) + ( N*Omega / 2mg^2))
  fFKRDM.Cz     = ( (fZeta * GA)/(6. * W) ) * (2. * W + ( (3.*q2 + nomg)/d ));


//____________________________________________________________________________

  const RSHelicityAmplModelDMI * hamplmod = 0;

//Propogator = ( (gZ'^4)/ (mZ'^2 + Q2)^2 )
//Full cross-section prefactor: (gZ'^4 / 8pi)*(1 / (mZ'^2 + Q^2)^2)*(-q^2/Q^2)
  double gZprime4 = TMath::Power(fgZp, 4);
//Average over polarizations. Compute Hel Amps in RRF. Include (W/Mnuc) for dsig/dWdq2
  double sig0 = (gZprime4/(16.*kPi))*sigcommon*TMath::Power((W/Mnuc),2);

//Prefactor for L and R terms is 1.
//Prefactor for S and Z terms:
  double scSZ  = -Q2star/q2;
  //____________________________________________________________________________
  double sigL =0;
  double sigR =0;
  double sigS =0;
  double sigZ =0;//new for DM
  //_________________________________________________________________________
  if(is_DM) {
      if (is_p) { hamplmod = fHAmplModelDMp;}
      else      { hamplmod = fHAmplModelDMn;}
    }
//____________________________________________________________________________
    // Compute the cross section
     assert(hamplmod);

     const RSHelicityAmplDM & hampl = hamplmod->Compute(resonance, fFKRDM);

     sigL = (hampl.Amp2Minus3 () + hampl.Amp2Minus1 ());
     sigR = (hampl.Amp2Plus3() + hampl.Amp2Plus1());
     sigS = scSZ * (hampl.Amp20Plus () + hampl.Amp20Minus());
     sigZ = scSZ * (hampl.Ampz20Plus () + hampl.Ampz20Minus());
//____________________________________________________________________________

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "sig_{0} = " << sig0;
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "sig_{L} = " << sigL;
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "sig_{R} = " << sigR;
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "sig_{S} = " << sigS;
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG) << "sig_{Z} = " << sigZ;
#endif

//____________________________________________________________________________
//DM Cross Section Calculation:

//XSec calc:
  double xsec = 0.0;
//Fermions:
if (fVelMode == 0) {
  double fQchiA2 = TMath::Power(fQchiA, 2); //QA^2
  double fQchiV2 = TMath::Power(fQchiV, 2); //QV^2
  //additional term from 4th polarization: - ( (mZ'^2 - q^2)^2 / (mZ'^4) )
  double Xterm = -(TMath::Power((mZprime2 - q2), 2)/TMath::Power(mZprime2, 2));

  double L = fQchiA2*E12_Q2_A + fQchiV2*E12_Q2_V + fQchiV*fQchiA*VXA;
  double R = fQchiA2*E12_Q2_A + fQchiV2*E12_Q2_V - fQchiV*fQchiA*VXA;
  double S = 2.*fQchiA2*E12_Q2_SA + 2.*fQchiV2*E12_Q2;
  double Z = fQchiA2*Q2*(Xplus+Xminus)*denom*Xterm;

     if (is_dm) {
         xsec = sig0*(L*sigL + R*sigR + S*sigS + Z*sigZ);
     }
     else
     if (is_dmbar) {
         xsec = sig0*(R*sigR + L*sigL + S*sigS + Z*sigZ);
     }
   }
   else if (fVelMode == 2) {
     double fQchiS2 = TMath::Power(fQchiS,2); //QS^2

    double RL = fQchiS2*E12_Q2_SA;
    double S = fQchiS2*2.*E12*denom;

      xsec = sig0*(RL*sigL + RL*sigR + S*sigS);
     }
  xsec = TMath::Max(0.,xsec);
//____________________________________________________________________________












//____________________________________________________________________________
//Breit-Wigner:Un-changede from resonant neutrino case

  // Check whether the cross section is to be weighted with a Breit-Wigner distribution
  // (default: true)
  double bw = 1.0;
  if ( fWghtBW ) {
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR);
  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseDMRESPXSec2014", pDEBUG)
      << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
#endif
  xsec *= bw;

//____________________________________________________________________________
//____________________________________________________________________________
//Jacobian: Updated for resonant DM case.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BSKLNBaseDMRESPXSec2014", pINFO)
      << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
      << "](W=" << W << ", q2=" << q2 << ", E=" << E << ") = " << xsec;
#endif

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2fE ) {
     double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);          //Jacobian function is calculated in KineUtils.cxx
     xsec *= J;                                                                 //passing "interaction" and "kps" where
  }                                                                             //the calculation goes FROM {W,Q2} TO "kps"}
                                                                                //since the calculation is done in {W,Q2} ps.
  // If requested return the free nucleon xsec even for input nuclear tgt       //kps and kPSWQ2fE should be values held in KinePhaseSpace.h
  if ( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;                 //Interaction is held in Interaction.cxx(h) where we have introduced
                                                                                //RESDM for resonant DM interactions with ScatteringType kScDarkMatterResonant
  int Z = target.Z();                                                           //which is stored in ScatteringType.h
  int A = target.A();                                                           //These Values are then passed to KineUtils to calculate the Jacobian.
  int N = A-Z;                                                                  //Within KineUtils, the limits are determined in KPhaseSpace
                                                                                //KPhaseSpace determines the limits based on the ProcessInfo tags of the scattering types
                                                                                //which, here, is "IsDarkMatterResonant"
  // Take into account the number of scattering centers in the target           //so, we must now update KPhaseSpace.cxx for IsDarkMatterResonant (done)
  int NNucl = (is_p) ? Z : N;                                                   //Then, as KPhaseSpace calculates kinematic limits, KineUtils should be updated to account for those limits
  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

//_________________________________________________________________________
//____________________________________________________________________________
//Pauli Blocking: Un-changed from resonant neutrino case

  if ( fUsePauliBlocking && A!=1 )
  {
    // Calculation of Pauli blocking according references:
    //
    //     [1] S.L. Adler,  S. Nussinov,  and  E.A.  Paschos,  "Nuclear
    //         charge exchange corrections to leptonic pion  production
    //         in  the (3,3) resonance  region,"  Phys. Rev. D 9 (1974)
    //         2125-2143 [Erratum Phys. Rev. D 10 (1974) 1669].
    //     [2] J.Y. Yu, "Neutrino interactions and  nuclear  effects in
    //         oscillation experiments and the  nonperturbative disper-
    //         sive  sector in strong (quasi-)abelian  fields,"  Ph. D.
    //         Thesis, Dortmund U., Dortmund, 2002 (unpublished).
    //     [3] E.A. Paschos, J.Y. Yu,  and  M. Sakuda,  "Neutrino  pro-
    //         duction  of  resonances,"  Phys. Rev. D 69 (2004) 014013
    //         [arXiv: hep-ph/0308130].

    double P_Fermi = 0.0;

    // Maximum value of Fermi momentum of target nucleon (GeV)
    if ( A<6 || ! fUseRFGParametrization )
    {
        // look up the Fermi momentum for this target
        FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
        const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
        P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }

     double FactorPauli_RES = 1.0;

     double k0 = 0., q = 0., q0 = 0.;

     if (P_Fermi > 0.)
     {
        k0 = (W2-Mnuc2-Q2)/(2*W);
        k = TMath::Sqrt(k0*k0+Q2);  // previous value of k is overridden
        q0 = (W2-Mnuc2+kPionMass2)/(2*W);
        q = TMath::Sqrt(q0*q0-kPionMass2);
     }

     if ( 2*P_Fermi < k-q )
        FactorPauli_RES = 1.0;
     if ( 2*P_Fermi >= k+q )
        FactorPauli_RES = ((3*k*k+q*q)/(2*P_Fermi)-(5*TMath::Power(k,4)+TMath::Power(q,4)+10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
     if ( 2*P_Fermi >= k-q && 2*P_Fermi <= k+q )
        FactorPauli_RES = ((q+k)*(q+k)-4*P_Fermi*P_Fermi/5-TMath::Power(k-q, 3)/(2*P_Fermi)+TMath::Power(k-q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);

     xsec *= FactorPauli_RES;
  }
  return xsec;
}
//____________________________________________________________________________






//______________________________________________________________________________

double BSKLNBaseDMRESPXSec2014::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}

//____________________________________________________________________________
bool BSKLNBaseDMRESPXSec2014::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();                      //ProcessInfo.cxx(h) has been updated to include
  const XclsTag &      xcls       = interaction->ExclTag();                       //IsDarkMatterResonant and IsDarkMatterResonant
                                                                                  //by including kScDarkMatterResonant and kScDarkMatter
                                                                                  //in ScatteringType.h which is a bookkeeping file.
  if(!proc_info.IsDarkMatterResonant()) return false;
  if(!xcls.KnownResonance())  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));

  if (!is_pn) return false;


  int  probe   = init_state.ProbePdg();
  bool is_dm   = proc_info.IsDarkMatter();

  if(!is_dm) return false;

  return true;
}
//____________________________________________________________________________



//______________added DM V and A charges
//____________________________________________________________________________
void BSKLNBaseDMRESPXSec2014::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________
void BSKLNBaseDMRESPXSec2014::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________

//____________________________________________________________________________
void BSKLNBaseDMRESPXSec2014::LoadConfig(void)
{
  // dark matter couplings to mediator
  double QchiL, QchiR;
  this->GetParam( "DarkLeftCharge", QchiL ) ;
  this->GetParam( "DarkRightCharge", QchiR ) ;
  this->GetParam( "DarkScalarCharge", fQchiS ) ;
  fQchiV = 0.5*(QchiL + QchiR);
  fQchiA = 0.5*(- QchiL + QchiR);

  // velocity dependence of interaction: fermion/scalar DM
  this->GetParamDef("velocity-mode", fVelMode, 0 );

  // mediator coupling
  this->GetParam("ZpCoupling", fgZp ) ;

  // mediator mass
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  // Load all configuration data or set defaults
//For Kinematic Expressions:
  this->GetParam( "RES-Zeta"   , fZeta  ) ;   ///< FKR parameter Zeta
  this->GetParam( "RES-Omega"  , fOmega ) ;   ///< FKR parameter Omega

//For Form Factors:
  this->GetParam( "minibooneGA", fGA    ) ;   //< axial transition form factor
  this->GetParam( "minibooneGV", fGV    ) ;   //< vector transition form factor
  double ma, mv ;
  this->GetParam( "RES-Ma", ma ) ;
  this->GetParam( "RES-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);    ///< (axial mass)^2
  fMv2 = TMath::Power(mv,2);     ///< (vector mass)^2

//For Breit-Wigner Distributions
  this->GetParamDef( "BreitWignerWeight", fWghtBW, true ) ; //< weight with resonance breit-wigner?
  this->GetParamDef( "BreitWignerNorm",   fNormBW, true);   //< normalize resonance breit-wigner to 1?


//For Pauli-Blocking
  this->GetParam("FermiMomentumTable", fKFTable); //< table of Fermi momentum (kF) constants for various nuclei
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization); //< use parametrization for fermi momentum insted of table?
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);  //< account for Pauli blocking?



//_________________________________________________________________________
// Load all the sub-algorithms needed

//Helicity amplitude calculations
  fHAmplModelDMp    = 0;
  fHAmplModelDMn    = 0;

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelDMp = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMp","Default")); //DM + p
  fHAmplModelDMn = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMn","Default")); //DM + n

  assert( fHAmplModelDMp );
  assert( fHAmplModelDMn );

//_________________________________________________________________________

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ; //< use a DIS/RES joining scheme
  fWcut = 999999;
  if(fUsingDisResJoin) {
    this->GetParam( "Wcut", fWcut ) ; //< apply DIS/RES joining scheme < Wcut
  }

  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model.
  this->GetParamDef( "MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0 ) ; //< limits allowed phase space for n=2 res
  this->GetParamDef( "MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0 ) ; //< limits allowed phase space for n=0 res
  this->GetParamDef( "MaxNWidthForGNRes", fGnResMaxNWidths, 4.0 ) ; //< limits allowed phase space for other res

//_________________________________________________________________________
  // Load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
