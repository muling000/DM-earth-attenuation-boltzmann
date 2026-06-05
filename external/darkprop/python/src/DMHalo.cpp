#include "darkprop.hpp"

void bind_DMHalo(py::module_ &m)
{
  m.def("norm_maxwell", &dp::norm_maxwell<double>, "v0"_a=dp::constants::v0);
  m.def("norm_esc", &dp::norm_esc<double>,
        "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0);
  m.def("fv_halo_1d", &dp::fv_halo_1d<double>,
        "v"_a, "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0);
  m.def("fv_halo_1d_earth", &dp::fv_halo_1d_earth<double>,
        "v"_a, "vearth"_a=dp::constants::vearth,
        "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0);
}
