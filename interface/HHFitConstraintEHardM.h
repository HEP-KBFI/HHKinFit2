/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraintEHardM_
#define HHFitConstraintEHardM_

#include "HHFitConstraint.h"
#include "HHFitObject.h"

#include <vector>

namespace HHKinFit2{
class HHFitConstraintEHardM : public HHFitConstraint {
 public:
  HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, double mass);

  void prepare(bool respectLimits=true);
  double getChi2() const;
  double getLikelihood() const;

 private:
  HHFitObject* const m_constrainedobject;
  double const m_mass;
  
};
}

#endif /* HHFitConstraintEHardM_ */
