//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Changes for Res DM by Zach Orr, Colorado State University
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/BoostedDarkMatter/XSection/DMRESPXSec.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMI.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplDM.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DMRESPXSec::DMRESPXSec() :
XSecAlgorithmI("genie::DMRESPXSec")
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
DMRESPXSec::DMRESPXSec(string config) :
XSecAlgorithmI("genie::DMRESPXSec", config)
{
  fNuTauRdSpl    = 0;
  fNuTauBarRdSpl = 0;
}
//____________________________________________________________________________
DMRESPXSec::~DMRESPXSec()
{
  if(fNuTauRdSpl)    delete fNuTauRdSpl;
  if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;
}
//____________________________________________________________________________
double DMRESPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  const Kinematics & kinematics = interaction -> Kine();
  LOG("DMRESPXSec", pDEBUG) << "Using v^" << fVelMode << " dependence";

  // Get kinematical parameters
  double W  = kinematics.W();
  double q2 = kinematics.q2();

  // Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("DMRes", pDEBUG)
         << "RES/DIS Join Scheme: XSec[RES, W=" << W
         << " >= Wcut=" << fWcut << "] = 0";
#endif
       return 0;
    }
  }

  // Get the input baryon resonance
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  bool        is_delta  = utils::res::IsDelta (resonance);

  // Get the dark matter, hit nucleon & DM current
  int  nucpdgc   = target.HitNucPdg();
  int  probepdgc = init_state.ProbePdg();
  bool is_dm     = pdg::IsDarkMatter         (probepdgc);
  bool is_dmbar  = pdg::IsAntiDarkMatter    (probepdgc);

//  bool is_lplus  = pdg::IsPosChargedLepton (probepdgc);
//  bool is_lminus = pdg::IsNegChargedLepton (probepdgc);

  bool is_p      = pdg::IsProton  (nucpdgc);
  bool is_n      = pdg::IsNeutron (nucpdgc);
  bool is_DM     = proc_info.IsDarkMatterResonant();

  //Need break for is not DM
    if(!is_DM) {
    return 0;
    }

  // Get baryon resonance parameters
  int    IR  = utils::res::ResonanceIndex    (resonance);
  int    LR  = utils::res::OrbitalAngularMom (resonance);
  double MR  = utils::res::Mass              (resonance);
  double WR  = utils::res::Width             (resonance);
  double NR  = fNormBW?utils::res::BWNorm    (resonance,fN0ResMaxNWidths,fN2ResMaxNWidths,fGnResMaxNWidths):1;

  // Following NeuGEN, avoid problems with underlying unphysical
  // model assumptions by restricting the allowed W phase space
  // around the resonance peak
  if (fNormBW) {
	if      (W > MR + fN0ResMaxNWidths * WR && IR==0) return 0.;
	else if (W > MR + fN2ResMaxNWidths * WR && IR==2) return 0.;
	else if (W > MR + fGnResMaxNWidths * WR)          return 0.;
  }

  //Auxillary Kinematic factors
  double W2     = TMath::Power(W,    2);       //W^2
  double Mnuc   = target.HitNucMass();   //Nucleon mass
  double mchi   = init_state.GetProbeP4(kRfHitNucRest)->M(); //DM mass: mx
  double Mnuc2  = TMath::Power(Mnuc, 2);       //m^2
  double sqrtq2 = TMath::Sqrt(-q2);            //Sqrt(-q^2)
  double mchi2 = TMath::Power(mchi, 2);        // mx^2
  double mZprime2 = TMath::Power(fMedMass, 2); // mZ'^2
  double E      = init_state.ProbeE(kRfHitNucRest);   //E1

  // Compute auxiliary & kinematical factors
  double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
  double v      = k - 0.5 * q2/Mnuc;
  double v2     = TMath::Power(v, 2);
  double kstar  = k * (Mnuc/W);
  double vstar  = kstar + 0.5 * (q2/W);
  double v2star = TMath::Power(vstar,2);
  double Q2     = v2 - q2;
  double Q      = TMath::Sqrt(Q2);
  double Q2star = v2star - q2;
  double Qstar  = TMath::Sqrt(Q2star);
  double Eprime = E - v;

  double E12 = TMath::Power((E + Eprime), 2);
  double Xminus = (4.*mchi2/q2) - 1;
  double Xplus = (4.*mchi2/q2) + 1;
  double denom = 1/(4. * (TMath::Power(E, 2) - mchi2));

  double E12_Q2 = (E12 - Q2)*denom;
  double E12_Q2_V = (E12 + Q2*Xplus)*denom;
  double E12_Q2_A = (E12 - Q2*Xminus)*denom;
  double E12_Q2_SA = (E12 + Q2*Xminus)*denom;
  double VXA = (4.*Q*(E + Eprime))*denom;


  // Calculate the Feynman-Kislinger-Ravndall parameters

  double Go  = TMath::Power(1 - 0.25 * q2/Mnuc2, 0.5-IR);
  double GV  = Go * TMath::Power( 1./(1-q2/fMv2), 2);
  double GA  = Go * TMath::Power( 1./(1-q2/fMa2), 2);


  //____________________________________________________________________________
  /*
  if(fGV){

    LOG("DMRESPXSec",pDEBUG) <<"Using new GV";
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
    LOG("DMRESPXSec",pDEBUG) << "Using new GA";

    double CA5_0 = 1.2;
    double CA5 = CA5_0 *  TMath::Power( 1./(1-q2/fMa2), 2);
      GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5/fZeta;
    GA = 0.5 * TMath::Sqrt(3.) * TMath::Power( 1 - q2/(Mnuc + W)/(Mnuc + W), 0.5-IR) * (1- (W2 +q2 -Mnuc2)/8./Mnuc2) * CA5;

    LOG("DMRESPXSec",pINFO) <<"GA= " <<GA << "  C5A= " <<CA5;
  }
  */
  //____________________________________________________________________________



  double d      = (TMath::Power(W+Mnuc,2.) - q2)/(2. * W);
  double sq2omg = TMath::Sqrt(2./fOmega);
  double nomg   = IR * fOmega;


  fFKRDM.Lamda  = sq2omg * Qstar;
  fFKRDM.Tv     = GV / (3.*W*sq2omg);
  fFKRDM.Ta     = ((GA * fZeta)/(3. * W)) * (1./sq2omg) * (Qstar / d);
  fFKRDM.Rv     = kSqrt2 * GV * Qstar * ((W + Mnuc)/(2. * W * d));
  fFKRDM.Ra     = ((kSqrt2 * GA * fZeta)/(6. * W)) * (W + Mnuc + (nomg/d));
  fFKRDM.S      = (-q2/Q2star) * (GV / (6. * W2)) * (3. * W * Mnuc + q2 - Mnuc2);
  fFKRDM.Bs     = ( (fZeta * GA)/(3. * W * sq2omg) ) * (1 + (vstar / d));
  fFKRDM.Cs     = ( (fZeta * GA)/(6. * W * Qstar) ) * (W2 - Mnuc2 + nomg*(vstar / d));
  fFKRDM.Bz     = ( (fZeta * GA)/(3. * W * Qstar * sq2omg) ) * (W - Mnuc);
  fFKRDM.Cz     = ( (fZeta * GA)/(6. * W) ) * (2. * W + ( (3.*q2 + nomg)/d ));


#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("FKR", pDEBUG)
     << "FKR params for RES = " << resname << " : " << fFKR;
#endif

  // Calculate the Rein-Sehgal Helicity Amplitudes

  const RSHelicityAmplModelDMI * hamplmod = 0;
  if(is_DM) {
      if (is_p) { hamplmod = fHAmplModelDMp;}
      else      { hamplmod = fHAmplModelDMn;}
    }
  assert(hamplmod);

  const RSHelicityAmplDM & hampl = hamplmod->Compute(resonance, fFKRDM);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMHAmpl", pDEBUG)
     << "Helicity Amplitudes for DMRES = " << resname << " : " << hampl;
#endif

  // Compute the cross section
  double gZprime4 = TMath::Power(fgZp, 4);
  double sigcommon = (1/TMath::Power((mZprime2 + Q2),2)) * (-q2/Q2);
  double sig0 = (gZprime4/(16.*kPi))*sigcommon*TMath::Power((W/Mnuc),2);
  double scSZ  = -Q2star/q2;

  double sigL =0;
  double sigR =0;
  double sigS =0;
  double sigZ =0;

  sigL = (hampl.Amp2Minus3 () + hampl.Amp2Minus1 ());
  sigR = (hampl.Amp2Plus3() + hampl.Amp2Plus1());
  sigS = scSZ * (hampl.Amp20Plus () + hampl.Amp20Minus());
  sigZ = scSZ * (hampl.Ampz20Plus () + hampl.Ampz20Minus());

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMRes", pDEBUG) << "sig_{0} = " << sig0;
  LOG("DMRes", pDEBUG) << "sig_{L} = " << sigL;
  LOG("DMRes", pDEBUG) << "sig_{R} = " << sigR;
  LOG("DMRes", pDEBUG) << "sig_{S} = " << sigS;
#endif

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
         xsec = sig0*(R*sigL + L*sigR + S*sigS + Z*sigZ);
     }
   }
   else if (fVelMode == 2) {
     double fQchiS2 = TMath::Power(fQchiS,2); //QS^2

    double RL = fQchiS2*E12_Q2_SA;
    double S = fQchiS2*2.*E12*denom;

      xsec = sig0*(RL*sigL + RL*sigR + S*sigS);
     }
  xsec = TMath::Max(0.,xsec);

  // Check whether the cross section is to be weighted with a
  // Breit-Wigner distribution (default: true)
  double bw = 1.0;
  if(fWghtBW) {
     bw = utils::bwfunc::BreitWignerL(W,LR,MR,WR,NR);

  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DMRes", pDEBUG)
       << "BreitWigner(RES=" << resname << ", W=" << W << ") = " << bw;
#endif
  xsec *= bw;



//________________________________________________________________________________________
  // Apply NeuGEN nutau cross section reduction factors
  double rf = 1.0;
  Spline * spl = 0;
  xsec *= rf;


  //Apply given scaling factor
  double xsec_scale = 1.;
  if      (is_DM) { xsec_scale = fXSecScaleDM; }
  xsec *= xsec_scale;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DMRes", pINFO)
    << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
          << "](W=" << W << ", q2=" << q2 << ", E=" << E << ") = " << xsec;
#endif

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;


  int Z = target.Z();
  int A = target.A();
  int N = A-Z;

  // Take into account the number of scattering centers in the target
  int NNucl = (is_p) ? Z : N;

  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)

  if (fUsePauliBlocking && A!=1)
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
     if (A<6 || !fUseRFGParametrization)
     {
         // Look up the Fermi momentum for this target
         FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
         const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
         P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // Define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // Correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }

     double FactorPauli_RES = 1.0;

     double k0 = 0., q = 0., q0 = 0.;

     if (P_Fermi > 0.)
     {
        k0 = (W2-Mnuc2-Q2)/(2*W);
        k = TMath::Sqrt(k0*k0+Q2);                  // previous value of k is overridden
        q0 = (W2-Mnuc2+kPionMass2)/(2*W);
        q = TMath::Sqrt(q0*q0-kPionMass2);
     }

     if (2*P_Fermi < k-q)
        FactorPauli_RES = 1.0;
     if (2*P_Fermi >= k+q)
        FactorPauli_RES = ((3*k*k+q*q)/(2*P_Fermi)-(5*TMath::Power(k,4)+TMath::Power(q,4)+10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
     if (2*P_Fermi >= k-q && 2*P_Fermi <= k+q)
        FactorPauli_RES = ((q+k)*(q+k)-4*P_Fermi*P_Fermi/5-TMath::Power(k-q, 3)/(2*P_Fermi)+TMath::Power(k-q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);

     xsec *= FactorPauli_RES;
  }

  return xsec;
}
//____________________________________________________________________________
double DMRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DMRESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsDarkMatterResonant()) return false;
  if(!xcls.KnownResonance())  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));

  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_DM = proc_info.IsDarkMatterResonant();
  bool is_dm   = proc_info.IsDarkMatter();

  if (!is_dm && !is_DM) return false;

  return true;
}
//____________________________________________________________________________
void DMRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMRESPXSec::LoadConfig(void)
{
  // Cross section scaling factor
  this->GetParam( "RES-DM-XSecScale", fXSecScaleDM ) ;

  //FKR Params
  this->GetParam( "RES-Zeta", fZeta ) ;
  this->GetParam( "RES-Omega", fOmega ) ;

  // velocity dependence of interaction: fermion/scalar DM
  this->GetParamDef("velocity-mode", fVelMode, 0 );
  //DM mass
  double ma, mv ;
  this->GetParam( "RES-Ma", ma ) ;
  this->GetParam( "RES-Mv", mv ) ;
  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);
  //DM Charges
  double QchiL, QchiR;
  this->GetParam( "DarkLeftCharge", QchiL ) ;
  this->GetParam( "DarkRightCharge", QchiR ) ;
  this->GetParam( "DarkScalarCharge", fQchiS ) ;
  fQchiV = 0.5*(QchiL + QchiR);
  fQchiA = 0.5*(- QchiL + QchiR);
  // DM Mediator
  this->GetParam("ZpCoupling", fgZp ) ;
  fMedMass = PDGLibrary::Instance()->Find(kPdgMediator)->Mass();

  this->GetParamDef( "BreitWignerWeight", fWghtBW, true ) ;
  this->GetParamDef( "BreitWignerNorm",   fNormBW, true);

//  this->GetParam( "minibooneGA", fGA    ) ;   //< axial transition form factor
//  this->GetParam( "minibooneGV", fGV    ) ;   //< vector transition form factor


  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

  // Load all the sub-algorithms needed

  fHAmplModelDMp    = 0;
  fHAmplModelDMn    = 0;

  AlgFactory * algf = AlgFactory::Instance();

  fHAmplModelDMp = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMp","Default")); //DM + p
  fHAmplModelDMn = dynamic_cast<const RSHelicityAmplModelDMI *> (
      algf->GetAlgorithm("genie::RSHelicityAmplModelDMn","Default")); //DM + n

  assert( fHAmplModelDMp );
  assert( fHAmplModelDMn );

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;
  fWcut = 999999;
  if(fUsingDisResJoin) {
    this->GetParam( "Wcut", fWcut ) ;
  }

  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model.
  this->GetParamDef( "MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0 ) ;
  this->GetParamDef( "MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0 ) ;
  this->GetParamDef( "MaxNWidthForGNRes", fGnResMaxNWidths, 4.0 ) ;

//______________________________________________________________________________________
  // NeuGEN reduction factors for nu_tau: a gross estimate of the effect of
  // neglected form factors in the R/S model
  this->GetParamDef( "UseNuTauScalingFactors", fUsingNuTauScaling, true ) ;
  if(fUsingNuTauScaling) {
     if(fNuTauRdSpl)    delete fNuTauRdSpl;
     if(fNuTauBarRdSpl) delete fNuTauBarRdSpl;

     assert( std::getenv( "GENIE") );
     string base = std::getenv( "GENIE") ;

     string filename = base + "/data/evgen/rein_sehgal/res/nutau_xsec_scaling_factors.dat";
     LOG("DMRes", pNOTICE)
                << "Loading nu_tau xsec reduction spline from: " << filename;
     fNuTauRdSpl = new Spline(filename);

     filename = base + "/data/evgen/rein_sehgal/res/nutaubar_xsec_scaling_factors.dat";
     LOG("DMRes", pNOTICE)
           << "Loading bar{nu_tau} xsec reduction spline from: " << filename;
     fNuTauBarRdSpl = new Spline(filename);
  }

//______________________________________________________________________________________


  // load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
