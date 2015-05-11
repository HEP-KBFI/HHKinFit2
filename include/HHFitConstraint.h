/*
 * class for general fit constraint
 */

#ifndef HHFitConstraint_
#define HHFitConstraint_

#include "HHFitObject.h"

namespace HHKinFit2{
class HHFitConstraint {
 public:
  HHFitConstraint(HHFitObject* fitobject);
  virtual ~HHFitConstraint() {};

  virtual double getChi2() const = 0;
  virtual double getLikelihood() const = 0;
  virtual void printChi2() const {};

 protected:
  HHFitObject* const m_fitobject;

};
}
#endif /* HHFitConstraint_ */
