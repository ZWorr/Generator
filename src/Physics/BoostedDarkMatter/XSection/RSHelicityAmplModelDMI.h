//____________________________________________________________________________
/*
RSHelicityAmplModelDMI: the interface for the Dark Matter helicity Amplitudes (header)

"Physics/Resonance/XSection/RSHelicityAmplModelI.h" -> "Physics/Resonance/XSection/RSHelicityAmplModelDMI.h"

Within RSHelicityAmplModelDMI.h make the following changes from RSHelicityAmplModelI.h:
  #include "Physics/Resonance/XSection/FKR.h" -> #include "Physics/Resonance/XSection/FKRDM.h"
  #include "Physics/Resonance/XSection/RSHelicityAmpl.h" -> #include "Physics/Resonance/XSection/RSHelicityAmplDM.h"
  _REIN_SEHGAL_HELICITY_AMPL_MODEL_I_H_ -> _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_
  RSHelicityAmplModelI -> RSHelicityAmplModelDMI
  RSHelicityAmpl -> RSHelicityAmplDM
  FKR -> FKRDM and fkr -> fkrdm



Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_
#define _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_

#include "Framework/Algorithm/Algorithm.h" //No change?
#include "Framework/ParticleData/BaryonResonance.h" //No change?
#include "Physics/Resonance/XSection/FKRDM.h" //Specific to DM
#include "Physics/Resonance/XSection/RSHelicityAmplDM.h" //Specific to DM

namespace genie {

class RSHelicityAmplModelDMI : public Algorithm
{
public:
//______________________ added for quark charges
    Configure(string config);
//______________________

  virtual ~RSHelicityAmplModelDMI();

  // define the RSHelicityAmplModelI interface
  virtual const RSHelicityAmplDM & Compute(Resonance_t res, const FKRDM & fkrdm) const = 0;



//______________added for quark charges
private:
  LoadConfig(void);

  double fQuV;
  double fQuA;
  double fQdV;
  double fQdA;
//________________


protected:
  RSHelicityAmplModelDMI();
  RSHelicityAmplModelDMI(string name);
  RSHelicityAmplModelDMI(string name, string config);
};

}        // namespace

#endif   // _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_
