/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateClusterDoseActor_h
#define GateClusterDoseActor_h

#include "G4Cache.hh"
#include "GateVActor.h"
#include "G4ParticleDefinition.hh"
#include "itkImage.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

class GateClusterDoseActor : public GateVActor {

public:
  using Image3DType = itk::Image<double, 3>;

  explicit GateClusterDoseActor(py::dict &user_info);

  void InitializeUserInfo(py::dict &user_info) override;

  void InitializeCpp() override;

  void BeginOfRunAction(const G4Run *run) override;

  void EndOfRunAction(const G4Run *run) override;

  void SteppingAction(G4Step *) override;

  void BeginOfEventAction(const G4Event *event) override;

  void BeginOfRunActionMasterThread(int run_id) override;

  std::string GetPhysicalVolumeName() const { return fPhysicalVolumeName; }

  void SetPhysicalVolumeName(std::string s) {
    fPhysicalVolumeName = std::move(s);
  }

  void AddProcessedLUT(const std::string &particleName,
                       const std::vector<double> &energyGrid,
                       const std::vector<double> &values);
  void ClearProcessedLUTs() { fIonisationDetailLUTByParticleName.clear(); }

  Image3DType::Pointer cpp_numerator_image;
  Image3DType::Pointer cpp_denominator_image;
  int NbOfEvent = 0;
  long fClampBelowRangeCount = 0;
  long fClampAboveRangeCount = 0;

private:
  struct IonisationDetailLUT {
    std::vector<double> energy;
    std::vector<double> values;
  };

  struct InterpolationResult {
    double value = 0.0;
    bool clampedBelow = false;
    bool clampedAbove = false;
  };

  struct threadLocalT {
    std::unordered_map<const G4ParticleDefinition *,
                       const IonisationDetailLUT *>
        ionisationDetailLUTByParticleDef;
    long clampBelowRangeCount = 0;
    long clampAboveRangeCount = 0;
  };

  InterpolationResult
  InterpolateIonisationDetailValue(const IonisationDetailLUT &lut,
                                   double energy) const;
  void ValidateIonisationDetailLUT(const std::string &particleName,
                                   const std::vector<double> &energyGrid,
                                   const std::vector<double> &values) const;

  std::string fPhysicalVolumeName;
  G4ThreeVector fTranslation;
  std::string fHitType;
  std::string fIonizationParameter;
  std::unordered_map<std::string, IonisationDetailLUT>
      fIonisationDetailLUTByParticleName;
  G4Cache<threadLocalT> fThreadLocalData;
};

#endif // GateClusterDoseActor_h
