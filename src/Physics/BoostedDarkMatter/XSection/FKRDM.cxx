//____________________________________________________________________________
/*
Dark Matter model parameters calculated in the same form as
the Feynmann-Kislinger-Ravndall (FKR) baryon excitation model parameters.

"Physics/Resonance/XSection/FKR.cxx" -> "Physics/Resonance/XSection/FKRDM.cxx"

Within FKRDM.cxx make the following changes from FKR.cxx:
  #include "Physics/Resonance/XSection/FKR.h" -> #include "Physics/Resonance/XSection/FKRDM.h"
  "FKR" -> "FKRDM"
  C -> Cs and Cz
  B -> Bs and Bz

  (Note: the DM calculation does not
  currently use variables {Rplus,Rminus,Tplus,Tminus,R,T} but I have
  left them defined in the code here just incase)

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/ParticleData/BaryonResUtils.h" //No change?
#include "Framework/Conventions/Constants.h" //No change?
#include "Framework/Messenger/Messenger.h" //No change?
#include "Physics/Resonance/XSection/FKRDM.h" //Specific to DM

using std::endl;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const FKRDM & parameters)
  {
     parameters.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
FKRDM::FKRDM()
{
  this->Reset();
}
//____________________________________________________________________________
FKRDM::~FKRDM()
{

}
//____________________________________________________________________________
void FKRDM::Print(ostream & stream) const
{
  stream << endl;
  stream << " lamda = " << Lamda   << endl;
  stream << " Tv    = " << Tv      << endl;
  stream << " Rv    = " << Rv      << endl;
  stream << " S     = " << S       << endl;
  stream << " Ta    = " << Ta      << endl;
  stream << " Ra    = " << Ra      << endl;
  stream << " Bs    = " << Bs      << endl;
  stream << " Cs    = " << Cs      << endl;
  stream << " Bz    = " << Bz      << endl;
  stream << " Cz    = " << Cz      << endl;
}
//____________________________________________________________________________
void FKRDM::Reset(void)
{
  Lamda        =  0.0;
  Tv           =  0.0;
  Rv           =  0.0;
  S            =  0.0;
  Ta           =  0.0;
  Ra           =  0.0;
  Bs           =  0.0;
  Cs           =  0.0;
  Bz           =  0.0;
  Cz           =  0.0;
}
//____________________________________________________________________________
