//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Cross Section Calculation (header)


Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#ifndef _BSKLN_BASE_RES_DM_PXSEC_2014_H_
#define _BSKLN_BASE_RES_DM_PXSEC_2014_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h" //Resonance_t
#include "Physics/Resonance/XSection/FKRDM.h"//Lambda,T,R,B,C,S

namespace genie {

  class RSHelicityAmplModelDMI;
  class Spline;
  class XSecIntegratorI;

  class BSKLNBaseDMRESPXSec2014: public XSecAlgorithmI {

    public:
      //______________added for DM charges
          Configure(string config);
          //______________
      virtual ~BSKLNBaseDMRESPXSec2014();

      // implement the XSecAlgorithmI interface
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);


//______________added for DM charges
private:
      LoadConfig(void);

      double fQchiV;
      double fQchiA;
      double fQchiS;
      int    fVelMode;
      double fMedMass;
      double fgZp;
      //______________

    protected:

      BSKLNBaseDMRESPXSec2014(string name);
      BSKLNBaseDMRESPXSec2014(string name, string config);

      void LoadConfig (void);

      mutable FKRDM fFKRDM;

      const RSHelicityAmplModelDMI * fHAmplModelDMp;
      const RSHelicityAmplModelDMI * fHAmplModelDMn;

//Update for DM:
//_________________________________________________________________________
      // configuration data
      const XSecIntegratorI * fXSecIntegrator;

      bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
      double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
      double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
      double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
      double   fGnResMaxNWidths;   ///< limits allowed phase space for other res


      bool     fWghtBW;            ///< weight with resonance breit-wigner?
      bool     fNormBW;            ///< normalize resonance breit-wigner to 1?
      double   fZeta;              ///< FKR parameter Zeta
      double   fOmega;             ///< FKR parameter Omega
      double   fMa2;               ///< (axial mass)^2
      double   fMv2;               ///< (vector mass)^2
      double   fSin48w;            ///< sin^4(Weingberg angle)
      double   fVud2;              ///< |Vud|^2(square of magnitude ud-element of CKM-matrix)


      bool     fUsingNuTauScaling; ///< use NeuGEN nutau xsec reduction factors?


      string fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
      bool fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
      bool fUsePauliBlocking;      ///< account for Pauli blocking?

      double   fXSecScaleCC;       ///< external CC xsec scaling factor
      double   fXSecScaleNC;       ///< external NC xsec scaling factor

      bool fKLN;
      bool fBRS;

      bool fGA;
      bool fGV;




  };

}       // genie namespace
//_________________________________________________________________________

#endif  // _BSKLN_BASE_RES_DM_PXSEC_2014_H_
