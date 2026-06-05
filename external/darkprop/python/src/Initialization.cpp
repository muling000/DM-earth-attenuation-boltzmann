#include "darkprop.hpp"

void bind_Initialization(py::module_ &m)
{
  m.def("sample_halo_speed", &dp::sample_halo_speed<double>,
        "rn"_a, "vmin"_a=0.0, "vmax"_a=dp::constants::vesc,
        "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0);
  m.def("random_isotropic_vector", &dp::random_isotropic_vector<Vector3d, double>,
        "rn"_a, "len"_a=1.0);
  m.def("sample_halo_velocity_earth_frame",
        &dp::sample_halo_velocity_earth_frame<Vector3d, double>,
        "v_earth"_a, "rn"_a, "vmin"_a=0.0, "vmax"_a=-1.0,
        "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0);
  m.def("random_incident_position", &dp::random_incident_position<Vector3d, double>, "rn"_a,
        "rn"_a, "vi", "r_earth"_a=dp::constants::rEarth, "R"_a=1.1*dp::constants::rEarth);
  m.def("init_halo", &dp::init_halo<Vector3d, double>,
        "p"_a, "t"_a, "vmin"_a, "v_earth"_a, "rn"_a, "R"_a,
        "vmax"_a=-1.0, "vesc"_a=dp::constants::vesc, "v0"_a=dp::constants::v0,
        "radius"_a=dp::constants::rEarth);
  m.def("init_Tbl", &dp::init_Tbl<Vector3d, double>,
        "p"_a, "T"_a, "b"_a, "l"_a, "rn"_a,
        "t"_a=0.0, "R"_a=1.1*dp::constants::rEarth, "radius"_a=dp::constants::rEarth);
  m.def("init_Tbl_vertical", &dp::init_Tbl_vertical<Vector3d, double>,
        "p"_a, "T"_a, "b"_a, "l"_a,
        "t"_a=0.0, "R"_a=dp::constants::rEarth);
  m.def("init_T_isotropic", &dp::init_T_isotropic<Vector3d, double>,
        "p"_a, "T"_a, "rn"_a,
        "t"_a=0.0, "R"_a=1.1*dp::constants::rEarth, "radius"_a=dp::constants::rEarth);
  m.def("init_rT_isotropic", &dp::init_rT_isotropic<Vector3d, double>,
        "p"_a, "r"_a, "T"_a, "rn"_a, "t"_a=0.0);
  m.def("inject", &dp::inject<Vector3d, double>,
        "p"_a, "radius"_a=dp::constants::rEarth);
}
