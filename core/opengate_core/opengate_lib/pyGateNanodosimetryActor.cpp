/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateNanodosimetryActor.h"

class PyGateNanodosimetryActor : public GateNanodosimetryActor
{
public:
  // Inherit the constructors
  using GateNanodosimetryActor::GateNanodosimetryActor;

  void BeginOfRunActionMasterThread(int run_id) override
  {
    PYBIND11_OVERLOAD(void, GateNanodosimetryActor, BeginOfRunActionMasterThread, run_id);
  }

  int EndOfRunActionMasterThread(int run_id) override
  {
    PYBIND11_OVERLOAD(int, GateNanodosimetryActor, EndOfRunActionMasterThread, run_id);
  }
};

void init_GateNanodosimetryActor(py::module &m)
{
  py::class_<GateNanodosimetryActor, PyGateNanodosimetryActor, std::unique_ptr<GateNanodosimetryActor, py::nodelete>, GateVActor>(m, "GateNanodosimetryActor")
      .def(py::init<py::dict &>())
      .def("InitializeUserInfo", &GateNanodosimetryActor::InitializeUserInfo)
      .def("BeginOfRunActionMasterThread", &GateNanodosimetryActor::BeginOfRunActionMasterThread)
      // .def("EndOfRunActionMasterThread", &GateNanodosimetryActor::EndOfRunActionMasterThread)

      .def("SetClusterDoseFlag", &GateNanodosimetryActor::SetClusterDoseFlag)
      .def("GetClusterDoseFlag", &GateNanodosimetryActor::GetClusterDoseFlag)
      .def("SetFluenceFlag", &GateNanodosimetryActor::GetClusterDoseFlag)
      .def("GetFluenceFlag", &GateNanodosimetryActor::GetClusterDoseFlag)
      .def("SetIpFlag", &GateNanodosimetryActor::GetClusterDoseFlag)
      .def("GetIpFlag", &GateNanodosimetryActor::GetClusterDoseFlag)
      // .def("GetPhysicalVolumeName", &GateNanodosimetryActor::GetPhysicalVolumeName)
      .def("SetPhysicalVolumeName", &GateNanodosimetryActor::SetPhysicalVolumeName)

      .def_readwrite("NbOfEvent", &GateNanodosimetryActor::NbOfEvent)
      .def_readwrite("cpp_clusterDose_image", &GateNanodosimetryActor::cpp_clusterDose_image)
      .def_readwrite("cpp_fluence_image", &GateNanodosimetryActor::cpp_fluence_image)
      .def_readwrite("cpp_ip_image", &GateNanodosimetryActor::cpp_ip_image)
      .def_readwrite("fPhysicalVolumeName", &GateNanodosimetryActor::fPhysicalVolumeName);
}
