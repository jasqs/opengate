/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateNanodosimetryActor.h"
#include "G4Navigator.hh"
#include "G4RandomTools.hh"
#include "G4RunManager.hh"
#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateHelpersImage.h"

#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

// Mutex that will be used by thread to write images
G4Mutex SetNanodosimetryPixelMutex = G4MUTEX_INITIALIZER;
G4Mutex SetNanodosimetryNbEventMutex = G4MUTEX_INITIALIZER;

GateNanodosimetryActor::GateNanodosimetryActor(py::dict &user_info) : GateVActor(user_info) {}

void GateNanodosimetryActor::InitializeUserInfo(py::dict &user_info)
{
  // IMPORTANT: call the base class method
  GateVActor::InitializeUserInfo(user_info);

  // translation
  fTranslation = DictGetG4ThreeVector(user_info, "translation");
  // Hit type (random, pre, post etc)
  fHitType = DictGetStr(user_info, "hit_type");
}

void GateNanodosimetryActor::InitializeCpp()
{
  GateVActor::InitializeCpp();
  NbOfThreads = G4Threading::GetNumberOfRunningWorkerThreads();

  // Create the image pointers
  // (the size and allocation will be performed on the py side)
  if (fClusterDoseFlag)
  {
    cpp_clusterDose_image = Image3DType::New();
  }
  if (fFluenceFlag)
  {
    cpp_fluence_image = Image3DType::New();
  }
  if (fIpFlag)
  {
    cpp_ip_image = Image3DType::New();
  }
}

void GateNanodosimetryActor::BeginOfRunActionMasterThread(int run_id)
{

  // Reset the number of events (per run)
  NbOfEvent = 0;

  // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<Image3DType>(cpp_clusterDose_image, fPhysicalVolumeName, fTranslation);
  auto spacing = cpp_clusterDose_image->GetSpacing();
  fVoxelVolume = spacing[0] * spacing[1] * spacing[2];

  // compute volume of a dose voxel
  Image3DType::RegionType region = cpp_clusterDose_image->GetLargestPossibleRegion();
  size_clusterDose = region.GetSize();

  if (fClusterDoseFlag)
  {
    AttachImageToVolume<Image3DType>(cpp_clusterDose_image, fPhysicalVolumeName, fTranslation);
  }
  if (fFluenceFlag)
  {
    AttachImageToVolume<Image3DType>(cpp_fluence_image, fPhysicalVolumeName, fTranslation);
  }
  if (fIpFlag)
  {
    AttachImageToVolume<Image3DType>(cpp_ip_image, fPhysicalVolumeName, fTranslation);
  }
}

// void GateDoseActor::BeginOfRunAction(const G4Run *run) {
//   int N_voxels = size_edep[0] * size_edep[1] * size_edep[2];
//   if (fEdepSquaredFlag) {
//     PrepareLocalDataForRun(fThreadLocalDataEdep.Get(), N_voxels);
//   }
//   if (fDoseSquaredFlag) {
//     PrepareLocalDataForRun(fThreadLocalDataDose.Get(), N_voxels);
//   }
// }

void GateNanodosimetryActor::BeginOfEventAction(const G4Event *event)
{
  G4AutoLock mutex(&SetNanodosimetryNbEventMutex);
  NbOfEvent++;
}

void GateNanodosimetryActor::GetVoxelPosition(G4Step *step, G4ThreeVector &position, bool &isInside, Image3DType::IndexType &index) const
{
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  // consider random position between pre and post
  if (fHitType == "pre")
  {
    position = preGlobal;
  }
  if (fHitType == "random")
  {
    auto x = G4UniformRand();
    auto direction = postGlobal - preGlobal;
    position = preGlobal + x * direction;
  }
  if (fHitType == "middle")
  {
    auto direction = postGlobal - preGlobal;
    position = preGlobal + 0.5 * direction;
  }

  auto localPosition = touchable->GetHistory()->GetTransform(0).TransformPoint(position);

  // convert G4ThreeVector to itk PointType
  Image3DType::PointType point;
  point[0] = localPosition[0];
  point[1] = localPosition[1];
  point[2] = localPosition[2];

  isInside = cpp_clusterDose_image->TransformPhysicalPointToIndex(point, index);
}

void GateNanodosimetryActor::SteppingAction(G4Step *step)
{
  auto event_id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  // Get the voxel index in which the stepping occurred
  G4ThreeVector position;
  bool isInside;
  Image3DType::IndexType index;
  GetVoxelPosition(step, position, isInside, index);

  G4double trackLength;

  if (isInside)
  {
    auto weight = step->GetTrack()->GetWeight();
    trackLength = step->GetStepLength();

    // all ImageAddValue calls in a mutexed {}-scope
    {
      G4AutoLock mutex(&SetNanodosimetryPixelMutex);
      if (fFluenceFlag)
      {
        // ImageAddValue<Image3DType>(cpp_clusterDose_image, index, clusterDose);
        ImageAddValue<Image3DType>(cpp_fluence_image, index, trackLength);
      }
    } // mutex scope
  }
}

// double GateNanodosimetryActor::ScoringQuantityFn(G4Step *step, double *secondQuantity)
// {

//   auto &l = fThreadLocalData.Get();
//   auto dedx_currstep = l.dedx_currstep;

//   if (fScoreInOtherMaterial)
//   {
//     auto SPR_otherMaterial = GetSPROtherMaterial(step);

//     dedx_currstep *= SPR_otherMaterial;
//   }
//   return dedx_currstep;
// }
