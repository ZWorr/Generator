//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Helicity Amplitudes for neutron Interactions

RSHelicityAmplModelEMn.cxx -> RSHelicityAmplModelDMn.cxx

Changes made to RSHelicityAmplModelEMn.cxx:
  #include "Physics/Resonance/XSection/RSHelicityAmplModelEMn.h" -> #include "Physics/Resonance/XSection/RSHelicityAmplModelDMn.h"
  "EM" -> "DM"
  RSHelicityAmplModelI -> RSHelicityAmplModelDMI
  FKR -> FKRDM and fkr -> fkrdm
  New terms specific to DM: fz0Plus and fz0Minus

  From "AhrensDMELPXSec.cxx"
  {QuV,QuA,QdV,QdA} -> {fQuV,fQuA,fQdV,fQdA}




Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h" //No change? <- why just for the n?
#include "Framework/ParticleData/BaryonResUtils.h" //No change? (FKRDM.cxx)
#include "Framework/Conventions/Constants.h" //No change? (FKRDM.cxx)
#include "Framework/Messenger/Messenger.h" //No change? (FKRDM.cxx)
#include "Physics/Resonance/XSection/RSHelicityAmplModelDMn.h"//Specific to DM
#include "Physics/Resonance/XSection/RSHelicityAmplDM.h" //Specific to DM <- why just for n?

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RSHelicityAmplModelEMn::RSHelicityAmplModelEMn() :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMn")
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMn::RSHelicityAmplModelEMn(string config) :
RSHelicityAmplModelI("genie::RSHelicityAmplModelEMn", config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelEMn::~RSHelicityAmplModelEMn()
{

}
//____________________________________________________________________________
const RSHelicityAmpl &
   RSHelicityAmplModelEMn::Compute(
           Resonance_t res, const FKR & fkr) const
{
  switch(res) {

    case (kP33_1232) :
    {
      fAmpl.fMinus3 = kSqrt6 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
      fAmpl.fMinus1 = kSqrt2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
      fAmpl.fPlus1 = kSqrt2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
      fAmpl.fPlus3 = kSqrt6 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
      fAmpl.f0Plus = 2 * kSqrt2 * fkrdm.Cs * (fQuA - fQdA);
      fAmpl.f0Minus = 2 * kSqrt2 * fkrdm.Cs * (fQuA - fQdA);
      fAmpl.fz0Plus = 2 * kSqrt2 * fkrdm.Cz * (fQuA - fQdA);
      fAmpl.fz0Minus = 2 * kSqrt2 * fkrdm.Cz * (fQuA - fQdA);

      break;
    }
    case (kS11_1535) :
    {
    fAmpl.fMinus3 = 0;
    fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQuA + 5 * fQdA) - fkrdm.Rv * (fQuV + 5 * fQdV)) + kSqrt3 * fkrdm.Ta * (fQdA - fQuA) + kSqrt3 * fkrdm.Tv * (fQuV - fQdV);
    fAmpl.fPlus1 = -k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQuA + 5 * fQdA) + fkrdm.Rv * (fQuV + 5 * fQdV)) + kSqrt3 * fkrdm.Ta * (fQuA - fQdA) + kSqrt3 * fkrdm.Tv * (fQuV - fQdV);
    fAmpl.fPlus3 = 0;
    fAmpl.f0Plus = k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQdV - fQuV) - 3 * fkrdm.Bs * (fQuA + 5 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQuA + 5 * fQdA));
    fAmpl.f0Minus = k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQuV - fQdV) - 3 * fkrdm.Bs * (fQuA + 5 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQuA + 5 * fQdA));
    fAmpl.fz0Plus = -k1_Sqrt6 * ((fQuA + 5 * fQdA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz));
    fAmpl.fz0Minus = -k1_Sqrt6 * ((fQuA + 5 * fQdA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz));

      break;
    }
    case (kD13_1520) :
    {
    fAmpl.fMinus3 = k3_Sqrt2 * (fkrdm.Tv * (fQuV - fQdV) + fkrdm.Ta * (fQdA - fQuA));
    fAmpl.fMinus1 = (1/2) * k1_Sqrt3 * (2 * fkrdm.Lamda * (fkrdm.Rv * (5 * fQdV + fQuV) - fkrdm.Ra * (5 * fQdA + fQuA)) + 3 * kSqrt2 * (fkrdm.Tv * (fQuV - fQdV) + fkrdm.Ta * (fQdA - fQuA)));
    fAmpl.fPlus1 = (1/2) * k1_Sqrt3 * (-2 * fkrdm.Lamda * (fkrdm.Ra * fQuA + fkrdm.Rv * (5 * fQdV + fQuV)) + fQdA * (3 * kSqrt2 * fkrdm.Ta - 10 * fkrdm.Lamda * fkrdm.Ra) + 3 * kSqrt2 * fkrdm.Tv * (fQdV - fQuV) - 3 * kSqrt2 * fQuA * fkrdm.Ta);
    fAmpl.fPlus3 = -k3_Sqrt2 * (fkrdm.Tv * (fQdV - fQuV) + fkrdm.Ta * (fQdA - fQuA));
    fAmpl.f0Plus = -fkrdm.Lamda * k1_Sqrt3 * (3 * S * (fQdV - fQuV) + fkrdm.Cs * (fQuA + 5 * fQdA));
    fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt3 * (3 * S * (fQuV - fQdV) + fkrdm.Cs * (fQuA + 5 * fQdA));
    fAmpl.fz0Plus = -fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQuA + 5 * fQdA);
    fAmpl.fz0Minus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQuA + 5 * fQdA);

      break;
    }
    case (kS11_1650) :
    {
    fAmpl.fMinus3 = 0;
    fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda *(fkrdm.Rv * (2 * fQuV + fQdV) - fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.fPlus1 = k1_Sqrt6* fkrdm.Lamda *(fkrdm.Rv * (2 * fQuV + fQdV) + fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.fPlus3 = 0;
    fAmpl.f0Plus = -kSqrt2_3 * (2 * fQuA + fQdA) * (3 * fkrdm.Bs - fkrdm.Lamda * fkrdm.Cs);
    fAmpl.f0Minus = -kSqrt2_3 * (2 * fQuA + fQdA) * (3 * fkrdm.Bs - fkrdm.Lamda * fkrdm.Cs);
    fAmpl.fz0Plus = -kSqrt2_3 * (2 * fQuA + fQdA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);
    fAmpl.fz0Minus = -kSqrt2_3 * (2 * fQuA + fQdA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);

      break;
    }
    case (kD13_1700) :
    {
    fAmpl.fMinus3 = (1/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Rv * (2 * fQuV + fQdV) - fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.fMinus1 = k1_Sqrt30 * fkrdm.Lamda * (fkrdm.Rv * (2 * fQuV + fQdV) - fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.fPlus1 = -k1_Sqrt30 * fkrdm.Lamda * (fkrdm.Rv * (2 * fQuV + fQdV) + fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.fPlus3 = -(1/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Rv * (2 * fQuV + fQdV) + fkrdm.Ra * (2 * fQuA + fQdA));
    fAmpl.f0Plus = (kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.f0Minus = -(kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.fz0Plus = (kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cz * (2 * fQuA + fQdA);
    fAmpl.fz0Minus = -(kSqrt2 * k1_Sqrt15) * fkrdm.Lamda * fkrdm.Cz * (2 * fQuA + fQdA);

      break;
    }
    case (kD15_1675) :
    {
    fAmpl.fMinus3 = kSqrt3_5 * fkrdm.Lamda * (fkrdm.Ra * (2 * fQuA + fQdA) - fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fMinus1 = kSqrt3_10 * fkrdm.Lamda * (fkrdm.Ra * (2 * fQuA + fQdA) - fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fPlus1 = -kSqrt3_10 * fkrdm.Lamda * (fkrdm.Ra * (2 * fQuA + fQdA) + fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fPlus3 = -kSqrt3_5 * fkrdm.Lamda * (fkrdm.Ra * (2 * fQuA + fQdA) + fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.f0Plus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.f0Minus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.fz0Plus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cz * (2 * fQuA + fQdA);
    fAmpl.fz0Minus = -(kSqrt3 * kSqrt2_5) * fkrdm.Lamda * fkrdm.Cz * (2 * fQuA + fQdA);

      break;
    }
    case (kS31_1620) :
    {
    fAmpl.fMinus3 = 0;
    fAmpl.fMinus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV)) + kSqrt3 * fkrdm.Ta * (fQuA - fQdA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
    fAmpl.fPlus1 = k1_Sqrt6 * fkrdm.Lamda * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV)) + kSqrt3 * fkrdm.Ta * (fQdA - fQuA) + kSqrt3 * fkrdm.Tv * (fQdV - fQuV);
    fAmpl.fPlus3 = 0;
    fAmpl.f0Plus = k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQuV - fQdV) + 3 * fkrdm.Bs * (fQuA - fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - fQuA));
    fAmpl.f0Minus = k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQdV - fQuV) + 3 * fkrdm.Bs * (fQuA - fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQdA - fQuA));
    fAmpl.fz0Plus = -k1_Sqrt6 * (fQdA - fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);
    fAmpl.fz0Minus = -k1_Sqrt6 * (fQdA - fQuA) * (3 * fkrdm.Bz - fkrdm.Lamda * fkrdm.Cz);

      break;
    }
    case (kD33_1700) :
    {
    fAmpl.fMinus3 = -3 * k1_Sqrt2 * (fkrdm.Tv * (fQuV - fQdV) + fkrdm.Ta * (fQdA - fQuA));
    fAmpl.fMinus1 = (1/2) * k1_Sqrt3 * ((fQuA - fQdA) * (2 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) + (fQdV - fQuV) * (2 * fkrdm.Lamda * fkrdm.Rv + 3 * kSqrt2 * fkrdm.Tv));
    fAmpl.fPlus1 = (1/2) * k1_Sqrt3 * ((fQuA - fQdA) * (2 * fkrdm.Lamda * fkrdm.Ra + 3 * kSqrt2 * fkrdm.Ta) - (fQdV - fQuV) * (2 * fkrdm.Lamda * fkrdm.Rv + 3 * kSqrt2 * fkrdm.Tv));
    fAmpl.fPlus3 = -3 * k1_Sqrt2 * (fkrdm.Tv * (fQdV - fQuV) + fkrdm.Ta * (fQdA - fQuA));
    fAmpl.f0Plus = fkrdm.Lamda * k1_Sqrt3 * (3 * S * (fQdV - fQuV) + fkrdm.Cs * (fQuA - fQdA));
    fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt3 * (3 * S * (fQdV - fQuV) + fkrdm.Cs * (fQdA - fQuA));
    fAmpl.fz0Plus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQuA - fQdA);
    fAmpl.fz0Minus = fkrdm.Lamda * fkrdm.Cz * k1_Sqrt3 * (fQdA - fQuA);

      break;
    }
    case (kP11_1440) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

      fAmpl.fMinus3 = 0;
      fAmpl.fMinus1 = (1/2) * k1_Sqrt3 * L2 * (fkrdm.Ra * (fQuA - 4 fQdA) - fkrdm.Rv * (fQuV - 4 * fQdV));
      fAmpl.fPlus1 = (1/2) * k1_Sqrt3 * L2 * (fkrdm.Ra * (fQuA - 4 fQdA) + fkrdm.Rv * (fQuV - 4 * fQdV));
      fAmpl.fPlus3 = 0;
      fAmpl.f0Plus = (fkrdm.Lamda/2) * k1_Sqrt3 * (-3 * fkrdm.Lamda * S * (fQuV + 2 * fQdV) - 2 * fkrdm.Bs * (fQuA - 4 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQuA - 4 * fQdA));
      fAmpl.f0Minus = (fkrdm.Lamda/2) * k1_Sqrt3 * (-3 * fkrdm.Lamda * S * (fQuV + 2 * fQdV) + 2 * fkrdm.Bs * (fQuA - 4 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (4 * fQdA - fQuA));
      fAmpl.fz0Plus = (fkrdm.Lamda/2) * k1_Sqrt3 * ((fQuA - 4 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz));
      fAmpl.fz0Minus = -(fkrdm.Lamda/2) * k1_Sqrt3 * ((fQuA - 4 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz));

      break;
    }
    case (kP33_1600) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = k1_Sqrt2 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fMinus1 = k1_Sqrt6 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus1 = k1_Sqrt6 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus3 = k1_Sqrt2 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.f0Plus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2 * fkrdm.Bs);
    fAmpl.f0Minus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 2 * fkrdm.Bs);
    fAmpl.fz0Plus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz);
    fAmpl.fz0Minus = kSqrt2_3 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz);

      break;
    }
    case (kP13_1720) :
    {
    fAmpl.fMinus3 = (1/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Ta * (fQuA + 2 * fQdA) - fkrdm.Tv * (fQuV + 2 * fQdV));
    fAmpl.fMinus1 = k1_Sqrt15 * (fkrdm.Lamda/2) * (2 * fkrdm.Lamda * (fkrdm.Ra * (fQuA - 4 * fQdA) - fkrdm.Rv * (fQuV - 4 * fQdV)) + 9 * kSqrt2 * (fkrdm.Tv * (fQuV + 2 * fQdV) - fkrdm.Ta * (fQuA + 2 * fQdA)));
    fAmpl.fPlus1 = k1_Sqrt15 * (fkrdm.Lamda/2) * (9 * kSqrt2 * (fkrdm.Tv * (fQuV + 2 * fQdV) + fkrdm.Ta * (fQuA + 2 * fQdA)) - 2 * fkrdm.Lamda * (fkrdm.Ra * (fQuA - 4 * fQdA) + fkrdm.Rv * (fQuV - 4 * fQdV)));
    fAmpl.fPlus3 = -(1/kSqrt10) * 3 * fkrdm.Lamda * (fkrdm.Ta * (fQuA + 2 * fQdA) + fkrdm.Tv * (fQuV + 2 * fQdV));
    fAmpl.f0Plus = fkrdm.Lamda * k1_Sqrt15 * (-3 * fkrdm.Lamda * S * (fQuV + 2 * fQdV) - 5 * fkrdm.Bs * (fQuA - 4 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQuA - 4 * fQdA));
    fAmpl.f0Minus = fkrdm.Lamda * k1_Sqrt15 * (3 * fkrdm.Lamda * S * (fQuV + 2 * fQdV) - 5 * fkrdm.Bs * (fQuA - 4 * fQdA) + fkrdm.Lamda * fkrdm.Cs * (fQuA - 4 * fQdA));
    fAmpl.fz0Plus = fkrdm.Lamda * k1_Sqrt15 * ((fQuA - 4 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz));
    fAmpl.fz0Minus = fkrdm.Lamda * k1_Sqrt15 * ((fQuA - 4 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz));

      break;
    }
    case (kF15_1680) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = 3 * fkrdm.Lamda * kSqrt2_5 * (fkrdm.Tv * (fQuV + 2 * fQdV) - fkrdm.Ta * (fQuA + 2 * fQdA));
    fAmpl.fMinus1 = k1_Sqrt5 * (fkrdm.Lamda/2) * (kSqrt2 * fkrdm.Lamda * (fkrdm.Rv * (fQuV - 4 * fQdV) - fkrdm.Ra * (fQuA - 4 * fQdA)) + 6 * (fkrdm.Tv * (fQuV + 2 * fQdV) - fkrdm.Ta * (fQuA + 2 * fQdA)));
    fAmpl.fPlus1 = k1_Sqrt5 * (fkrdm.Lamda/2) * (-kSqrt2 * fkrdm.Lamda * (fkrdm.Rv * (fQuV - 4 * fQdV) + fkrdm.Ra * (fQuA - 4 * fQdA)) - 6 * (fkrdm.Tv * (fQuV + 2 * fQdV) + fkrdm.Ta * (fQuA + 2 * fQdA)));
    fAmpl.fPlus3 = -3 * fkrdm.Lamda * kSqrt2_5 * (fkrdm.Tv * (fQuV + 2 * fQdV) + fkrdm.Ta * (fQuA + 2 * fQdA));
    fAmpl.f0Plus = (1/kSqrt10) * L2 * (3 * S * (fQuV + 2 fQdV) - fkrdm.Cs * (fQuA - 4 * fQdA));
    fAmpl.f0Minus = (1/kSqrt10) * L2 * (3 * S * (fQuV + 2 fQdV) + fkrdm.Cs * (fQuA - 4 * fQdA));
    fAmpl.fz0Plus = -(1/kSqrt10) * L2 * fkrdm.Cz * (fQuA - 4 * fQdA);
    fAmpl.fz0Minus = (1/kSqrt10) * L2 * fkrdm.Cz * (fQuA - 4 * fQdA);

      break;
    }
    case (kP31_1910) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = 0;
    fAmpl.fMinus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fPlus1 = k1_Sqrt15 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus3 = 0;
    fAmpl.f0Plus = -k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
    fAmpl.f0Minus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
    fAmpl.fz0Plus = -k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);
    fAmpl.fz0Minus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);

      break;
    }
    case (kP33_1920) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fMinus1 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus1 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus3 = k1_Sqrt5 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.f0Plus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
    fAmpl.f0Minus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cs - 5 * fkrdm.Bs);
    fAmpl.fz0Plus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);
    fAmpl.fz0Minus = k1_Sqrt15 * 2 * fkrdm.Lamda * (fQdA - fQuA) * (fkrdm.Lamda * fkrdm.Cz - 5 * fkrdm.Bz);

      break;
    }
    case (kF35_1905) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = -3 * L2 * (kSqrt2 * k1_Sqrt35) * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fMinus1 = L2 * k1_Sqrt35 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.fPlus1 = L2 * k1_Sqrt35 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fPlus3 = -3 * L2 * (kSqrt2 * k1_Sqrt35) * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQdV - fQuV));
    fAmpl.f0Plus = 2 * L2 * k1_Sqrt35 * fkrdm.Cs * (fQdA - fQuA);
    fAmpl.f0Minus = 2 * L2 * k1_Sqrt35 * fkrdm.Cs * (fQuA - fQdA);
    fAmpl.fz0Plus = 2 * L2 * k1_Sqrt35 * fkrdm.Cz * (fQdA - fQuA);
    fAmpl.fz0Minus = 2 * L2 * k1_Sqrt35 * fkrdm.Cz * (fQuA - fQdA);

      break;
    }
    case (kF37_1950) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = kSqrt2_7 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fMinus1 = kSqrt6_35 * L2 * (fkrdm.Ra * (fQdA - fQuA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fPlus1 = kSqrt6_35 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.fPlus3 = kSqrt2_7 * L2 * (fkrdm.Ra * (fQuA - fQdA) + fkrdm.Rv * (fQuV - fQdV));
    fAmpl.f0Plus = 2 * kSqrt6_35 * L2 * fkrdm.Cs * (fQuA - fQdA);
    fAmpl.f0Minus = 2 * kSqrt6_35 * L2 * fkrdm.Cs * (fQuA - fQdA);
    fAmpl.fz0Plus = 2 * kSqrt6_35 * L2 * fkrdm.Cz * (fQuA - fQdA);
    fAmpl.fz0Minus = 2 * kSqrt6_35 * L2 * fkrdm.Cz * (fQuA - fQdA);

      break;
    }
    case (kP11_1710) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = 0;
    fAmpl.fMinus1 = (1/2) * k1_Sqrt6 * L2 * (fkrdm.Ra * (fQuA + 5 * fQdA) - fkrdm.Rv * (fQuV + 5 * fQdV));
    fAmpl.fPlus1 = (1/2) * k1_Sqrt6 * L2 * (fkrdm.Ra * (fQuA + 5 * fQdA) + fkrdm.Rv * (fQuV + 5 * fQdV));
    fAmpl.fPlus3 = 0;
    fAmpl.f0Plus = fkrdm.Lamda * (1/2) * k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQdV - fQuV) + (fQuA + 5 * fQdA) * (fkrdm.Lamda * fkrdm.Cs - 2 * fkrdm.Bs));
    fAmpl.f0Minus = -fkrdm.Lamda * (1/2) * k1_Sqrt6 * (3 * fkrdm.Lamda * S * (fQdV - fQuV) + (fQuA + 5 * fQdA) * (fkrdm.Lamda * fkrdm.Cs - 2 * fkrdm.Bs));
    fAmpl.fz0Plus = fkrdm.Lamda * (1/2) * k1_Sqrt6 * (fQuA + 5 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz);
    fAmpl.fz0Minus = -fkrdm.Lamda * (1/2) * k1_Sqrt6 * (fQuA + 5 * fQdA) * (fkrdm.Lamda * fkrdm.Cz - 2 * fkrdm.Bz);

      break;
    }
    case (kF17_1970) :
    {
      double L2  = TMath::Power(fkrdm.Lamda, 2);

    fAmpl.fMinus3 = (k1_Sqrt2 * k1_Sqrt7) * L2 * (fkrdm.Ra * (2 * fQuA + fQdA) - fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fMinus1 = (kSqrt3_2 * k1_Sqrt35) * L2 * (fkrdm.Ra * (2 * fQuA + fQdA) - fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fPlus1 = -(kSqrt3_2 * k1_Sqrt35) * L2 * (fkrdm.Ra * (2 * fQuA + fQdA) + fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.fPlus3 = -(k1_Sqrt2 * k1_Sqrt7) * L2 * (fkrdm.Ra * (2 * fQuA + fQdA) + fkrdm.Rv * (2 * fQuV + fQdV));
    fAmpl.f0Plus = -kSqrt6_35 * L2 * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.f0Minus = -kSqrt6_35 * L2 * fkrdm.Cs * (2 * fQuA + fQdA);
    fAmpl.fz0Plus = -kSqrt6_35 * L2 * fkrdm.Cz * (2 * fQuA + fQdA);
    fAmpl.fz0Minus = -kSqrt6_35 * L2 * fkrdm.Cz * (2 * fQuA + fQdA);

      break;
    }

   default:
   {
     LOG("RSHAmpl", pWARN) << "*** UNRECOGNIZED RESONANCE!";
     fAmpl.fMinus1 = 0.;
     fAmpl.fPlus1  = 0.;
     fAmpl.fMinus3 = 0.;
     fAmpl.fPlus3  = 0.;
     fAmpl.f0Minus = 0.;
     fAmpl.f0Plus  = 0.;
     break;
   }

  }//switch

  return fAmpl;
}
//____________________________________________________________________________
