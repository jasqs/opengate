/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateClusterDoseActor.h"

#include "GateHelpersDict.h"
#include "GateHelpersImage.h"

#include "G4ParticleTable.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

G4Mutex SetPixelIonisationDetailMutex = G4MUTEX_INITIALIZER;
G4Mutex SetNbEventMutexIonisationDetail = G4MUTEX_INITIALIZER;
G4Mutex IonisationDetailCountersMutex = G4MUTEX_INITIALIZER;

GateClusterDoseActor::GateClusterDoseActor(py::dict &user_info)
    : GateVActor(user_info, true) {
  fActions.insert("SteppingAction");
  fActions.insert("BeginOfRunAction");
  fActions.insert("EndOfRunAction");
}

void GateClusterDoseActor::InitializeUserInfo(py::dict &user_info) {
  GateVActor::InitializeUserInfo(user_info);
  fTranslation = DictGetG4ThreeVector(user_info, "translation");
  fHitType = DictGetStr(user_info, "hit_type");
  fIonizationParameter = DictGetStr(user_info, "ionization_parameter");
}

void GateClusterDoseActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  cpp_numerator_image = Image3DType::New();
  cpp_denominator_image = Image3DType::New();
}

void GateClusterDoseActor::ValidateIonisationDetailLUT(
    const std::string &particleName, const std::vector<double> &energyGrid,
    const std::vector<double> &values) const {
  if (particleName.empty()) {
    throw std::invalid_argument(
        "Ionisation-detail LUT particle name cannot be empty.");
  }
  if (energyGrid.size() != values.size()) {
    throw std::invalid_argument(
        "Ionisation-detail LUT energy and value arrays must have the same size "
        "for particle '" +
        particleName + "'.");
  }
  if (energyGrid.size() < 2) {
    throw std::invalid_argument(
        "Ionisation-detail LUT requires at least two points for particle '" +
        particleName + "'.");
  }
  if (!std::is_sorted(energyGrid.begin(), energyGrid.end())) {
    throw std::invalid_argument(
        "Ionisation-detail LUT energy grid must be sorted for particle '" +
        particleName + "'.");
  }
}

void GateClusterDoseActor::AddProcessedLUT(const std::string &particleName,
                                           const std::vector<double> &energyGrid,
                                           const std::vector<double> &values) {
  ValidateIonisationDetailLUT(particleName, energyGrid, values);
  fIonisationDetailLUTByParticleName[particleName] =
      IonisationDetailLUT{energyGrid, values};
}

void GateClusterDoseActor::BeginOfRunAction(const G4Run * /*run*/) {
  auto &l = fThreadLocalData.Get();
  l.ionisationDetailLUTByParticleDef.clear();
  l.clampBelowRangeCount = 0;
  l.clampAboveRangeCount = 0;

  auto *particleTable = G4ParticleTable::GetParticleTable();
  for (const auto &[particleName, lut] : fIonisationDetailLUTByParticleName) {
    const auto *particleDefinition = particleTable->FindParticle(particleName);
    if (particleDefinition != nullptr) {
      l.ionisationDetailLUTByParticleDef[particleDefinition] = &lut;
    }
  }
}

void GateClusterDoseActor::EndOfRunAction(const G4Run * /*run*/) {
  auto &l = fThreadLocalData.Get();
  G4AutoLock mutex(&IonisationDetailCountersMutex);
  fClampBelowRangeCount += l.clampBelowRangeCount;
  fClampAboveRangeCount += l.clampAboveRangeCount;
}

void GateClusterDoseActor::BeginOfEventAction(const G4Event * /*event*/) {
  G4AutoLock mutex(&SetNbEventMutexIonisationDetail);
  NbOfEvent++;
}

void GateClusterDoseActor::BeginOfRunActionMasterThread(int run_id) {
  AttachImageToVolume<Image3DType>(cpp_numerator_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_denominator_image, fPhysicalVolumeName,
                                   fTranslation);
  NbOfEvent = 0;
  fClampBelowRangeCount = 0;
  fClampAboveRangeCount = 0;
}

GateClusterDoseActor::InterpolationResult
GateClusterDoseActor::InterpolateIonisationDetailValue(
    const IonisationDetailLUT &lut, const double energy) const {
  if (lut.energy.empty() || lut.values.empty() ||
      lut.energy.size() != lut.values.size()) {
    return {};
  }

  if (energy <= lut.energy.front()) {
    return {lut.values.front(), true, false};
  }
  if (energy >= lut.energy.back()) {
    return {lut.values.back(), false, true};
  }

  const auto upper = std::upper_bound(lut.energy.begin(), lut.energy.end(), energy);
  const auto upperIndex = static_cast<size_t>(
      std::distance(lut.energy.begin(), upper));
  const auto lowerIndex = upperIndex - 1;

  const auto e0 = lut.energy[lowerIndex];
  const auto e1 = lut.energy[upperIndex];
  const auto f0 = lut.values[lowerIndex];
  const auto f1 = lut.values[upperIndex];

  if (e1 == e0) {
    return {f0, false, false};
  }

  const auto weight = (energy - e0) / (e1 - e0);
  return {f0 + weight * (f1 - f0), false, false};
}

void GateClusterDoseActor::SteppingAction(G4Step *step) {
  G4ThreeVector position;
  bool isInside = false;
  Image3DType::IndexType index;
  GetStepVoxelPosition<Image3DType>(step, fHitType, cpp_numerator_image,
                                    position, isInside, index);

  if (!isInside) {
    return;
  }

  auto &l = fThreadLocalData.Get();
  const auto *particleDefinition = step->GetTrack()->GetParticleDefinition();
  const auto lutIt =
      l.ionisationDetailLUTByParticleDef.find(particleDefinition);
  if (lutIt == l.ionisationDetailLUTByParticleDef.end()) {
    return;
  }

  const auto preEnergy =
      step->GetPreStepPoint()->GetKineticEnergy() / CLHEP::MeV;
  const auto postEnergy =
      step->GetPostStepPoint()->GetKineticEnergy() / CLHEP::MeV;
  (void)fIonizationParameter;

  const auto preValue =
      InterpolateIonisationDetailValue(*lutIt->second, preEnergy);
  const auto postValue =
      InterpolateIonisationDetailValue(*lutIt->second, postEnergy);
  l.clampBelowRangeCount +=
      static_cast<long>(preValue.clampedBelow) +
      static_cast<long>(postValue.clampedBelow);
  l.clampAboveRangeCount +=
      static_cast<long>(preValue.clampedAbove) +
      static_cast<long>(postValue.clampedAbove);

  const auto value = std::abs(preValue.value - postValue.value);
  if (value <= 0.0) {
    return;
  }

  G4AutoLock mutex(&SetPixelIonisationDetailMutex);
  ImageAddValue<Image3DType>(cpp_numerator_image, index, value);
  ImageAddValue<Image3DType>(cpp_denominator_image, index, 1.0);
}
