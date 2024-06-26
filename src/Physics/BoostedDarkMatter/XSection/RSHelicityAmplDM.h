//____________________________________________________________________________
/*
Dark Matter Calculations of |f|^2 terms for the cross section, based on the
procedure of FKR. (header)

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________




//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmpl

\brief    A class holding the Rein-Sehgal's helicity amplitudes.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          RSHelicityAmplModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _RS_HELICITY_AMPL_DM_H_
#define _RS_HELICITY_AMPL_DM_H_

#include <iostream>

#include <TMath.h>

#include "Framework/Interaction/Interaction.h"
#include "Physics/BoostedDarkMatter/XSection/FKRDM.h"

using std::ostream;

namespace genie {

class RSHelicityAmplDM;
ostream & operator<< (ostream & stream, const RSHelicityAmplDM & hamp);

class RSHelicityAmplDM {
//Only two possibilites
friend class RSHelicityAmplModelDMp;
friend class RSHelicityAmplModelDMn;

public:

  RSHelicityAmplDM();
  RSHelicityAmplDM(const RSHelicityAmplDM & hamp);
  ~RSHelicityAmplDM() { }

  //! return helicity amplitude
  double AmpMinus1 (void) const  { return fMinus1; } /* f(-1) */
  double AmpPlus1  (void) const  { return fPlus1;  } /* f(+1) */
  double AmpMinus3 (void) const  { return fMinus3; } /* f(-3) */
  double AmpPlus3  (void) const  { return fPlus3;  } /* f(+3) */
  double Amp0Minus (void) const  { return f0Minus; } /* f(0-) */
  double Amp0Plus  (void) const  { return f0Plus;  } /* f(0+) */
  //New terms specific to DM
  double Ampz0Minus (void) const  { return fz0Minus; } /* fz(0-) */
  double Ampz0Plus  (void) const  { return fz0Plus;  } /* fz(0+) */

  //! return |helicity amplitude|^2
  double Amp2Minus1 (void) const { return TMath::Power(fMinus1, 2.); } /* |f(-1)|^2 */
  double Amp2Plus1  (void) const { return TMath::Power(fPlus1,  2.); } /* |f(+1)|^2 */
  double Amp2Minus3 (void) const { return TMath::Power(fMinus3, 2.); } /* |f(-3)|^2 */
  double Amp2Plus3  (void) const { return TMath::Power(fPlus3,  2.); } /* |f(+3)|^2 */
  double Amp20Minus (void) const { return TMath::Power(f0Minus, 2.); } /* |f(0-)|^2 */
  double Amp20Plus  (void) const { return TMath::Power(f0Plus,  2.); } /* |f(0+)|^2 */
  //New terms specific to DM
  double Ampz20Minus (void) const { return TMath::Power(fz0Minus, 2.); } /* |fz(0-)|^2 */
  double Ampz20Plus  (void) const { return TMath::Power(fz0Plus,  2.); } /* |fz(0+)|^2 */

  friend ostream & operator<< (ostream & stream, const RSHelicityAmplDM & hamp);

  void Print(ostream & stream) const;

private:

  void   Init(void);

  double fMinus1;
  double fPlus1;
  double fMinus3;
  double fPlus3;
  double f0Minus;
  double f0Plus;
  //New terms specific to DM
  double fz0Minus;
  double fz0Plus;
};

}        // genie namespace

#endif   // _RS_HELICITY_AMPL_DM_H_
