#include "darkprop.hpp"

using RandomNumber = dp::RandomNumber<double>;
using Event = dp::Event<Vector3d, double>;

void bind_Base(py::module_ &m)
{
  py::class_<dp::Target>(m, "Target")
    .def(py::init<>())
    .def(py::init<const std::string&, int, int, double>())
    .def_readwrite("Z", &dp::Target::Z)
    .def_readwrite("A", &dp::Target::A)
    .def_readwrite("m", &dp::Target::m)
    .def_readwrite("name", &dp::Target::name);

  py::class_<RandomNumber>(m, "RandomNumber")
    .def(py::init<int>(), "seed"_a=1)
    .def("uniform_phi", &RandomNumber::uniform_phi)
    .def("uniform_costheta", &RandomNumber::uniform_costheta)
    .def("uniform_xi", &RandomNumber::uniform_xi)
    .def("uniform_ab", &RandomNumber::uniform_ab)
    .def("normal", &RandomNumber::normal);

  py::class_<Event>(m, "Event")
    .def(py::init<double, double, Vector3d, Vector3d, double>())
    .def_readwrite("t", &Event::t)
    .def_readwrite("T", &Event::T)
    .def_readwrite("r", &Event::r)
    .def_readwrite("p3", &Event::p3)
    .def_readwrite("weight", &Event::weight);

  m.def("scatter", &dp::scatter<Vector3d, double>);
  m.def("propagate", &dp::propagate<Vector3d, double>);
}
