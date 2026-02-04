#pragma once

#include "pybind11/functional.h"
#include "pybind11/operators.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"

namespace py = pybind11;

// Modules
void ParserModule(py::module &);
void VectorModule(py::module &);
void TransformModule(py::module &);
void ShapeModule(py::module &);
void ElementModule(py::module &);
void MaterialModule(py::module &);
void VolumeModule(py::module &);
void WorldModule(py::module &);
