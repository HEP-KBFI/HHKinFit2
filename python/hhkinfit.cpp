#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
  std::vector<double> getFittedB1();
  std::vector<double> getFittedB2();
  std::vector<double> getFittedTau1();
  std::vector<double> getFittedTau2();
  double getMH();
  double getChi2();
  double getFitProb();
  int getConvergence();
  void setB1(std::vector<double>);
  void setB2(std::vector<double> p4);
  void setTauVis1(std::vector<double> p4);
  void setTauVis2(std::vector<double> p4);
  void setMetCov(std::vector<double> cov);
  void setMet(std::vector<double> metxy);
  void setBJetRes(std::vector<double> res);
  float getTau1Ratio();
  float getTau2Ratio();

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
    heavyhiggsfit = new HHKinFitMasterHeavyHiggs(bjet1, bjet2, tauvis1, tauvis2, met, *met_cov, bjetres1, bjetres2);
    heavyhiggsfit->fit();
    fittedTau1 = heavyhiggsfit->getFittedTau1();
    fittedTau2 = heavyhiggsfit->getFittedTau2();
    fittedBJet1 = heavyhiggsfit->getFittedBJet1();
    fittedBJet2 = heavyhiggsfit->getFittedBJet2();
  }

void HHFitterCF::fit(){
  heavyhiggsfit = new HHKinFitMasterHeavyHiggs(bjet1, bjet2, tauvis1, tauvis2, met, *met_cov, bjetres1, bjetres2);
  heavyhiggsfit->fit();
  fittedTau1 = heavyhiggsfit->getFittedTau1();
  fittedTau2 = heavyhiggsfit->getFittedTau2();
  fittedBJet1 = heavyhiggsfit->getFittedBJet1();
  fittedBJet2 = heavyhiggsfit->getFittedBJet2();
}
std::vector<double> HHFitterCF::getFittedB1(){
  std::vector<double> out = {fittedBJet1.Pt(), fittedBJet1.Eta(), fittedBJet1.Phi(), fittedBJet1.M()};
  return out;
}
std::vector<double> HHFitterCF::getFittedB2(){
  std::vector<double> out = {fittedBJet2.Pt(), fittedBJet2.Eta(), fittedBJet2.Phi(), fittedBJet2.M()};
  return out;
}
std::vector<double> HHFitterCF::getFittedTau1(){
  std::vector<double> out = {fittedTau1.Pt(), fittedTau1.Eta(), fittedTau1.Phi(), fittedTau1.M()};
  return out;
}
std::vector<double> HHFitterCF::getFittedTau2(){
  std::vector<double> out = {fittedTau2.Pt(), fittedTau2.Eta(), fittedTau2.Phi(), fittedTau2.M()};
  return out;
}
double HHFitterCF::getMH(){ return heavyhiggsfit->getMH(); }
double HHFitterCF::getChi2(){ return heavyhiggsfit->getChi2(); }
double HHFitterCF::getFitProb(){ return heavyhiggsfit->getFitProb(); }
int HHFitterCF::getConvergence(){ return heavyhiggsfit->getConvergence(); }

void HHFitterCF::setB1(std::vector<double> p4){
  double pt = p4[0];
  double eta = p4[1];
  double phi = p4[2];
  double m = p4[3];
  bjet1.SetPtEtaPhiM(pt, eta, phi, m);
}
void HHFitterCF::setB2(std::vector<double> p4){
  double pt = p4[0];
  double eta = p4[1];
  double phi = p4[2];
  double m = p4[3];
  bjet2.SetPtEtaPhiM(pt, eta, phi, m);
}
void HHFitterCF::setTauVis1(std::vector<double> p4){
  double pt = p4[0];
  double eta = p4[1];
  double phi = p4[2];
  double m = p4[3];
  tauvis1.SetPtEtaPhiM(pt, eta, phi, m);
}
void HHFitterCF::setTauVis2(std::vector<double> p4){
  double pt = p4[0];
  double eta = p4[1];
  double phi = p4[2];
  double m = p4[3];
  tauvis2.SetPtEtaPhiM(pt, eta, phi, m);
}
void HHFitterCF::setMet(std::vector<double> metxy){
  // met.SetX(metxy[0]);
  // met.SetY(metxy[1]);
  met.SetMagPhi(metxy[0], metxy[1]);
}
void HHFitterCF::setMetCov(std::vector<double> cov){
  (*met_cov)[0][0]= cov[0];
  (*met_cov)[0][1]= cov[1];
  (*met_cov)[1][0]= cov[1];
  (*met_cov)[1][1]= cov[2];
}
void HHFitterCF::setBJetRes(std::vector<double> res){
  bjetres1=res[0];
  bjetres2=res[1];
}

float HHFitterCF::getTau1Ratio(){
  return fittedTau1.E()/tauvis1.E();
}

float HHFitterCF::getTau2Ratio(){
  return fittedTau2.E()/tauvis2.E();
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

// hh_kinfit.py code
bool isDummy(const std::vector<double>& obj) {
    return obj.size() == 4 && obj[0] == 0.0 && obj[1] == 0.0 && obj[2] == 0.0 && obj[3] == 1.0;
}

std::vector<std::vector<double>> batchFit(
    const std::vector<std::vector<double>>& bjet1_list,
    const std::vector<std::vector<double>>& bjet2_list,
    const std::vector<std::vector<double>>& tauvis1_list,
    const std::vector<std::vector<double>>& tauvis2_list,
    const std::vector<std::vector<double>>& met_list,
    const std::vector<std::vector<double>>& met_cov_list,
    const std::vector<double>& bjet_resolutions
) {
    std::vector<std::vector<double>> results;

    for (size_t i = 0; i < bjet1_list.size(); i++) {
        if (isDummy(bjet1_list[i]) || isDummy(bjet2_list[i]) || isDummy(tauvis1_list[i]) || isDummy(tauvis2_list[i])) {
            results.push_back(std::vector<double>(6 * bjet_resolutions.size(), -99999.0));
            continue;
        }

        std::vector<double> event_results; // results vector initialized
        for (double b_res : bjet_resolutions) { // loop over the b-jet resolution values
            // std::cout << "using bjet resolution: " << b_res << std::endl;
            try {
                HHFitterCF fitter;
                // Set inputs
                fitter.setB1(bjet1_list[i]);
                fitter.setB2(bjet2_list[i]);
                fitter.setTauVis1(tauvis1_list[i]);
                fitter.setTauVis2(tauvis2_list[i]);
                fitter.setMet(met_list[i]);
                fitter.setMetCov(met_cov_list[i]);
                fitter.setBJetRes({b_res, b_res});

                fitter.fit();  // Perform the fit
                // std::cout << "Fit is done for event  " << i << " with b_res = " << b_res << std::endl;

                // Store fit results
                event_results.insert(event_results.end(), {fitter.getMH(), fitter.getTau1Ratio(), fitter.getTau2Ratio(), 
                                                          fitter.getFitProb(), fitter.getChi2(), static_cast<double>(fitter.getConvergence())});
            }
            catch (const std::exception& e) {
                // std::cerr << "Fit failed for event " << i << " with b_res = " << b_res << ": " << e.what() << std::endl;
                event_results.insert(event_results.end(), { -99999.0, -99999.0, -99999.0, -99999.0, -99999.0, -99999.0 });
            }
        }
        results.push_back(event_results);
    }
    return results;
}




namespace py = pybind11;

PYBIND11_MODULE(hhkinfit2, m) {
  m.doc() = "pybind11 version of hhkinfit2"; // optional module docstring

  m.def("fitSingleEvent2DKinFit", &fitSingleEvent2DKinFit, "fittest");

  m.def("batchFit", &batchFit, "Batch process multiple HH KinFit events");

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
