#include "darkprop.hpp"
#include <pybind11/operators.h>

void bind_Vector3d(py::module_ &m)
{
  py::class_<Vector3d>(m, "Vector3d")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def(py::init<const Vector3d&>())
    .def(py::init<const std::vector<double>&>())
    .def(py::self + py::self)
    .def(py::self + double())
    .def(double() + py::self)
    .def(py::self - py::self)
    .def(py::self - double())
    .def(double() - py::self)
    .def(-py::self)
    .def(double() * py::self)
    .def(py::self * double())
    .def(py::self / double())
    .def(py::self += py::self)
    .def(py::self += double())
    .def(py::self -= py::self)
    .def(py::self -= double())
    .def(py::self == py::self)
    .def("isApprox", &Vector3d::isApprox)
    .def("norm", &Vector3d::norm)
    .def("normalized", &Vector3d::normalized)
    .def("dot", &Vector3d::dot)
    .def("cross", &Vector3d::cross)
    .def("to_vec", &Vector3d::to_vec)
    .def("__repr__", [](const Vector3d& v) {
      return "Vector3d[" + std::to_string(v.x) + ", "
                         + std::to_string(v.y) + ", "
                         + std::to_string(v.z) + "]";});

    py::implicitly_convertible<py::list, Vector3d>();
    py::implicitly_convertible<py::tuple, Vector3d>();
}
