#include "../interface/HHFitConstraintSoftBoundary.h"
#include "../interface/exceptions/HHCovarianceMatrixException.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include <cmath>
#include <iostream>
#include <sstream>


HHKinFit2::HHFitConstraintSoftBoundary::HHFitConstraintSoftBoundary(HHFitObject* object, double percent)
  : HHFitConstraint(object),
    m_percent(percent)
{
}

double
HHKinFit2::HHFitConstraintSoftBoundary::getChi2() const{
  return(-2*log(getLikelihood()));
}

double
HHKinFit2::HHFitConstraintSoftBoundary::getLikelihood() const{
  double value = 1;
  double efit = m_fitobject->getFit4Vector().E();
  double elimit = m_fitobject->getInitial4Vector().E();
  if (efit<elimit)
    value = TMath::Gaus(efit,elimit,m_percent*elimit);
  return(value);
}
