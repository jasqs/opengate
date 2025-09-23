/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateNanodosimetryActor_h
#define GateNanodosimetryActor_h

#include "G4Cache.hh"
#include "G4EmCalculator.hh"
// #include "G4NistManager.hh"
#include "G4VPrimitiveScorer.hh"
// #include "GateHelpersImage.h"
// #include "GateWeightedEdepActor.h"
#include "GateVActor.h"
#include "itkImage.h"
// #include <pybind11/stl.h>

namespace py = pybind11;

class GateNanodosimetryActor : public GateVActor
{

public:
  // Image type needs to be 3D double by default
  typedef itk::Image<double, 3> Image3DType;

  // Constructor
  explicit GateNanodosimetryActor(py::dict &user_info);

  void InitializeUserInfo(py::dict &user_info) override;

  void InitializeCpp() override;

  void BeginOfRunActionMasterThread(int run_id) override;

  void BeginOfEventAction(const G4Event *event) override;

  // Called every time a Run starts (all threads)
  // void BeginOfRunAction(const G4Run *run) override;

  void SteppingAction(G4Step *) override;

  void GetVoxelPosition(G4Step *step, G4ThreeVector &position, bool &isInside, Image3DType::IndexType &index) const;

  inline void SetPhysicalVolumeName(std::string s) { fPhysicalVolumeName = s; }
  
  // double ScoringQuantityFn(G4Step *step, double *secondQuantity) override;

  //   void AddValuesToImages(G4Step *step,itk::Image<double, 3>::IndexType
  //   index) override;

  // Option: which quantities are to be saved
  bool fClusterDoseFlag{}; // cluster dose
  bool fFluenceFlag{};     // fluence, i.e. total track length per voxel volume
  bool fIpFlag{};          // track-length weighted average Ip

  void SetClusterDoseFlag(const bool b) { fClusterDoseFlag = b; }
  bool GetClusterDoseFlag() const { return fClusterDoseFlag; }

  void SetFluenceFlag(const bool b) { fFluenceFlag = b; }
  bool GetFluenceFlag() const { return fFluenceFlag; }

  void SetIpFlag(const bool b) { fIpFlag = b; }
  bool GetIpFlag() const { return fIpFlag; }

  int NbOfEvent = 0;
  int NbOfThreads = 0;

  // The image is accessible on py side (shared by all threads)
  Image3DType::Pointer cpp_clusterDose_image;
  Image3DType::Pointer cpp_fluence_image;
  Image3DType::Pointer cpp_ip_image;
  Image3DType::SizeType size_clusterDose{};
  double fVoxelVolume{};

  std::string fPhysicalVolumeName;
  G4ThreeVector fTranslation;
  std::string fHitType;

  // std::string fAveragingMethod;
};

#endif // GateNanodosimetryActor_h