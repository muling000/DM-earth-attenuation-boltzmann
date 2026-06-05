#include "darkprop.hpp"
#include "darkprop/Utils.hpp"

void bind_Utils(py::module_ &m)
{
  m.def("rotation_axis_angle", &dp::rotation_axis_angle<Vector3d, double>);
  m.def("scatter_fixed_target", &dp::scatter_fixed_target<Vector3d, double>);
  m.def("go_straight_and_check_the_boundary", &dp::go_straight_and_check_the_boundary<Vector3d, double>);
  m.def("reduce_m", &dp::reduce_m<double>);
}
