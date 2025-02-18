cmake_minimum_required(VERSION 3.5)

project(HHKINFIT2)

find_package(ROOT REQUIRED COMPONENTS Physics Hist Graf)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(pybind11 REQUIRED)

add_library(libHHKinFit2 SHARED src/HHFitConstraint4VectorBJet.cpp src/HHFitConstraint4Vector.cpp src/HHFitConstraint.cpp src/HHFitConstraintEHardM.cpp src/HHFitConstraintLikelihood.cpp src/HHFitConstraintSoftBoundary.cpp src/HHFitObjectComposite.cpp src/HHFitObject.cpp src/HHFitObjectEConstBeta.cpp src/HHFitObjectEConstM.cpp src/HHFitObjectE.cpp src/HHFitObjectMET.cpp src/HHKinFit.cpp src/HHKinFitMasterHeavyHiggs.cpp src/HHKinFitMasterSingleHiggs.cpp src/HHKinFitMasterSingleHiggsSoftLimits.cpp src/HHLorentzVector.cpp src/PSMath.cpp interface/HHFitConstraint4VectorBJet.h interface/HHFitConstraint4Vector.h interface/HHFitConstraintEHardM.h interface/HHFitConstraint.h interface/HHFitConstraintLikelihood.h interface/HHFitConstraintSoftBoundary.h interface/HHFitObjectComposite.h interface/HHFitObjectEConstBeta.h interface/HHFitObjectEConstM.h interface/HHFitObjectE.h interface/HHFitObject.h interface/HHFitObjectMET.h interface/HHKinFit.h interface/HHKinFitMasterHeavyHiggs.h interface/HHKinFitMasterSingleHiggs.h interface/HHKinFitMasterSingleHiggsSoftLimits.h interface/HHLorentzVector.h interface/PSMath.h interface/exceptions/HHCovarianceMatrixException.h interface/exceptions/HHEnergyConstraintException.h interface/exceptions/HHEnergyRangeException.h interface/exceptions/HHInvMConstraintException.h interface/exceptions/HHLimitSettingException.h)
target_include_directories(libHHKinFit2 PUBLIC include)
target_compile_features(libHHKinFit2 PUBLIC cxx_std_11)
target_link_libraries(libHHKinFit2 PUBLIC ROOT::Physics ROOT::Hist ROOT::Graf)

add_executable(fitSingleEvent2DKinFit examples/fitSingleEvent2DKinFit.C )
target_link_libraries(fitSingleEvent2DKinFit PUBLIC libHHKinFit2 )

find_package(pybind11 REQUIRED)
pybind11_add_module(testbind python/hhkinfit.cpp)
target_link_libraries(testbind PUBLIC libHHKinFit2)
install(TARGETS testbind
  COMPONENT python
  LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
  )