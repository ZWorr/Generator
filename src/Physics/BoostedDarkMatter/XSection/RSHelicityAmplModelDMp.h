//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Amplitude Calculation for proton Interactions (header)

Taken from the neutrino case for EM interactions:
  "Physics/Resonance/XSection/RSHelicityAmplModelEMp.h" -> "Physics/Resonance/XSection/RSHelicityAmplModelDMp.h"

Within RSHelicityAmplModelDMp.h make the changes:
  #include "Physics/Resonance/XSection/RSHelicityAmplModelI.h" -> #include "Physics/Resonance/XSection/RSHelicityAmplModelDMI.h"
  Change all "EM" to "DM"
  RSHelicityAmplModelDMI interface implementation: add "DM" to all RSHelicityAmpl
  FKR -> FKRDM and fkr -> fkrdm
  ? fAmpl -> fAmplDM ?


Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________


#ifndef _HELICITY_AMPL_MODEL_DM_P_H_
#define _HELICITY_AMPL_MODEL_DM_P_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelDMI.h" //Specific to DM

namespace genie {

class RSHelicityAmplModelDMp : public RSHelicityAmplModelDMI {

public:
  RSHelicityAmplModelDMp();
  RSHelicityAmplModelDMp(string config);
  virtual ~RSHelicityAmplModelDMp();

  // RSHelicityAmplModelDMI interface implementation
  const RSHelicityAmplDM & Compute(Resonance_t res, const FKRDM & fkrdm) const;

private:
  mutable RSHelicityAmplDM fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_DM_N_H_
