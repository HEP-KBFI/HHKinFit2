/*
 * class for fit objects with constant mass
 */

#ifndef HHFitObjectConstM_
#define HHFitObjectConstM_

#include "HHLorentzVector.h"
#include "HHFitObjectE.h"
#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObjectEConstM : public HHFitObjectE {
 public:
  HHFitObjectEConstM(HHLorentzVector const& initial4vector);
  HHLorentzVector constrainEtoMinv(double m, HHLorentzVector const& other4vector) const;
  double calculateEConstrainedToMinv(double m, HHLorentzVector const& other4vector) const;
  HHLorentzVector changeE(double E) const;

  void print() const;

};
}
#endif /* HHFitObjectConstM_ */
