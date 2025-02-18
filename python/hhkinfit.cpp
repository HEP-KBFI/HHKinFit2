#include <pybind11/pybind11.h>

#include <iostream>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector2.h"
#include "../interface/HHKinFitMasterHeavyHiggs.h"

using namespace HHKinFit2;

class HHFitterCF{
public:
  HHFitterCF();
  void fit();
  double* getFittedB1();
  double* getFittedB2();
  double* getFittedTau1();
  double* getFittedTau2();
  double getMH();
  double getChi2();
  double getFitProb();
  int getConvergence();
  void setB1(double px, double py, double pz, double energy);
  void setB2(double px, double py, double pz, double energy);
  void setTauVis1(double px, double py, double pz, double energy);
  void setTauVis2(double px, double py, double pz, double energy);
  void setMetCov(double xx, double xy, double yx, double yy);
  void setMet(double metx, double mety);
  void setBJetRes(double res1, double res2);
  private:
  TLorentzVector bjet1;
  TLorentzVector bjet2;
  TLorentzVector tauvis1;
  TLorentzVector tauvis2;
  TVector2 met;
  TMatrixD* met_cov;
  HHKinFitMasterHeavyHiggs* heavyhiggsfit;
  TLorentzVector fittedTau1;
  TLorentzVector fittedTau2;
  TLorentzVector fittedBJet1;
  TLorentzVector fittedBJet2;
  double bjetres1;
  double bjetres2;
};


HHFitterCF::HHFitterCF(){
    bjet1=TLorentzVector(-18.0706,28.8731,19.1597,39.3622);
    bjet2=TLorentzVector(14.1252,71.0363,-160.683,176.315);
    tauvis1=TLorentzVector(-20.9251,-24.865,-4.66771,32.8318);
    tauvis2=TLorentzVector(53.6093,-24.6612,-52.1178,78.7299);
    met=TVector2(31.2118,-38.0866);
    met_cov= new TMatrixD(2,2);
    (*met_cov)[0][0]= 200;
    (*met_cov)[0][1]= 0;
    (*met_cov)[1][0]= 0;
    (*met_cov)[1][1]= 200;
    bjetres1 = -1;
    bjetres2 = -1;
    heavyhiggsfit = new HHKinFitMasterHeavyHiggs(bjet1, bjet2, tauvis1, tauvis2, met, *met_cov);
    heavyhiggsfit->fit();
    fittedTau1 = heavyhiggsfit->getFittedTau1();
    fittedTau2 = heavyhiggsfit->getFittedTau2();
    fittedBJet1 = heavyhiggsfit->getFittedBJet1();
    fittedBJet2 = heavyhiggsfit->getFittedBJet2();
  }

  void HHFitterCF::fit(){
    heavyhiggsfit = new HHKinFitMasterHeavyHiggs(bjet1, bjet2, tauvis1, tauvis2, met, *met_cov);
    heavyhiggsfit->fit();
    fittedTau1 = heavyhiggsfit->getFittedTau1();
    fittedTau2 = heavyhiggsfit->getFittedTau2();
    fittedBJet1 = heavyhiggsfit->getFittedBJet1();
    fittedBJet2 = heavyhiggsfit->getFittedBJet2();
  }
  double* HHFitterCF::getFittedB1(){
    double *out = new double[4];
    out[0] = fittedBJet1.Px();
    out[1] = fittedBJet1.Py();
    out[2] = fittedBJet1.Pz();
    out[3] = fittedBJet1.E();
    return out;
  }
  double* HHFitterCF::getFittedB2(){
    double *out = new double[4];
    out[0] = fittedBJet2.Px();
    out[2] = fittedBJet2.Py();
    out[2] = fittedBJet2.Pz();
    out[3] = fittedBJet2.E();
    return out;
  }
  double* HHFitterCF::getFittedTau1(){
    double *out = new double[4];
    out[0] = fittedTau1.Px();
    out[1] = fittedTau1.Py();
    out[2] = fittedTau1.Pz();
    out[3] = fittedTau1.E();
    return out;
  }
  double* HHFitterCF::getFittedTau2(){
    double *out = new double[4];
    out[0] = fittedTau2.Px();
    out[2] = fittedTau2.Py();
    out[2] = fittedTau2.Pz();
    out[3] = fittedTau2.E();
    return out;
  }
  double HHFitterCF::getMH(){ return heavyhiggsfit->getMH(); }
  double HHFitterCF::getChi2(){ return heavyhiggsfit->getChi2(); }
  double HHFitterCF::getFitProb(){ return heavyhiggsfit->getFitProb(); }
  int HHFitterCF::getConvergence(){ return heavyhiggsfit->getConvergence(); }

  void HHFitterCF::setB1(double px, double py, double pz, double energy){
    bjet1.SetPxPyPzE(px, py, pz, energy);
  }
  void HHFitterCF::setB2(double px, double py, double pz, double energy){
    bjet2.SetPxPyPzE(px, py, pz, energy);
  }
  void HHFitterCF::setTauVis1(double px, double py, double pz, double energy){
    tauvis1.SetPxPyPzE(px, py, pz, energy);
  }
  void HHFitterCF::setTauVis2(double px, double py, double pz, double energy){
    tauvis2.SetPxPyPzE(px, py, pz, energy);
  }
  void HHFitterCF::setMet(double metx, double mety){
    met.SetX(metx);
    met.SetY(mety);
  }
  void HHFitterCF::setMetCov(double xx, double xy, double yx, double yy){
    (*met_cov)[0][0]= xx;
    (*met_cov)[0][1]= xy;
    (*met_cov)[1][0]= yx;
    (*met_cov)[1][1]= yy;
  }
  void HHFitterCF::setBJetRes(double res1, double res2){
    bjetres1=res1;
    bjetres2=res2;
  }

float fitSingleEvent2DKinFit(){
  TLorentzVector bjet1(-18.0706,28.8731,19.1597,39.3622);
  TLorentzVector bjet2(14.1252,71.0363,-160.683,176.315);
  TLorentzVector tauvis1(-20.9251,-24.865,-4.66771,32.8318);
  TLorentzVector tauvis2(53.6093,-24.6612,-52.1178,78.7299);
  TVector2 met(31.2118,-38.0866);
  TMatrixD met_cov(2,2);
  met_cov[0][0]= 200;
  met_cov[0][1]= 0;
  met_cov[1][0]= 0;
  met_cov[1][1]= 200;

  HHKinFitMasterHeavyHiggs heavyhiggsfit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov);

  heavyhiggsfit.fit();

  std::cout << "Heavy Higgs fit finished." << std::endl;

  std::cout << "m_H: " << heavyhiggsfit.getMH() << "GeV" << std::endl;
  std::cout << "chi2: " << heavyhiggsfit.getChi2() << std::endl;
  std::cout << "Fit prob: " << heavyhiggsfit.getFitProb() << std::endl;
  std::cout << "Convergence: " << heavyhiggsfit.getConvergence() << std::endl;

  TLorentzVector fittedTau1 = heavyhiggsfit.getFittedTau1();
  TLorentzVector fittedTau2 = heavyhiggsfit.getFittedTau2();
  TLorentzVector fittedBJet1 = heavyhiggsfit.getFittedBJet1();
  TLorentzVector fittedBJet2 = heavyhiggsfit.getFittedBJet2();

  std::cout << "Energy of fitted tau1: " << fittedTau1.E() << std::endl;
  std::cout << "Energy of fitted tau2: " << fittedTau2.E() << std::endl;
  std::cout << "Energy of fitted BJet1: " << fittedBJet1.E() << std::endl;
  std::cout << "Energy of fitted BJet2: " << fittedBJet2.E() << std::endl;

  return (0);


}
namespace py = pybind11;

PYBIND11_MODULE(testbind, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  m.def("fitSingleEvent2DKinFit", &fitSingleEvent2DKinFit, "fittest");

  py::class_<HHFitterCF>(m, "HHFitterCF")
    .def(py::init())
    .def("fit", &HHFitterCF::fit)
    .def("getFittedB1", &HHFitterCF::getFittedB1)
    .def("getFittedB2", &HHFitterCF::getFittedB2)
    .def("getFittedTau1", &HHFitterCF::getFittedTau1)
    .def("getFittedTau2", &HHFitterCF::getFittedTau2)
    .def("getMH", &HHFitterCF::getMH)
    .def("getChi2", &HHFitterCF::getChi2)
    .def("getFitProb", &HHFitterCF::getFitProb)
    .def("getConvergence", &HHFitterCF::getConvergence)
    .def("setB1", &HHFitterCF::setB1)
    .def("setB2", &HHFitterCF::setB2)
    .def("setTauVis1", &HHFitterCF::setTauVis1)
    .def("setTauVis2", &HHFitterCF::setTauVis2)
    .def("setMetCov", &HHFitterCF::setMetCov)
    .def("setMet", &HHFitterCF::setMet)
    .def("setBJetRes", &HHFitterCF::setBJetRes);
}
