/*
 * class for composite objects
 */

#ifndef HHFitObjectComposite_
#define HHFitObjectComposite_

#include "HHLorentzVector.h"
#include "HHFitObject.h"
#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObjectComposite : public HHFitObject {
 public:
	  HHFitObjectComposite(std::vector<HHFitObject*> const& subobjects);
	  HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2);
	  HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3);
	  HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3, HHFitObject* subobject4, HHFitObject* subobject5);

 protected:
  HHLorentzVector getFit4Vector() const;
  HHLorentzVector getInitial4Vector() const;
  void setCovMatrix(TMatrixD const cov);
  TMatrixD getCovMatrix() const;

  void setSubobjects(std::vector<HHFitObject*> const& subobjects);
  void addSubobject(HHFitObject* subobject);

  virtual void print() const;

 private:
  std::vector<HHFitObject*> m_subobjects;
  bool m_cov_set;
};
}
#endif /* HHFitObjectComposite_ */
