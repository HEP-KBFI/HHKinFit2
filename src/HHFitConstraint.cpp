#ifdef HHKINFIT2
#include "HHFitConstraint.h"
#else
#include "../interface/HHFitConstraint.h"
#endif

HHKinFit2::HHFitConstraint::HHFitConstraint(HHFitObject* fitobject)
  :m_fitobject(fitobject){

}

void
HHKinFit2::HHFitConstraint::prepare(bool respectLimits){
}
