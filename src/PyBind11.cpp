#include <Python.h>
#include <pybind11/pybind11.h>
#include "Context.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc() = "Robosample bindings";

    py::class_<Context>(m, "Context")
        .def(py::init<>())
        .def("initializeFromFile", &Context::initializeFromFile, "Initialize from file")
        .def("Run", (void (Context::*)(void)) &Context::Run, "Run the simulation");
}
