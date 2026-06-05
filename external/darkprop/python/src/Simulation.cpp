#include "darkprop.hpp"

void bind_Simulation(py::module_ &m)
{
  m.def("simulate_cross_sphere", &dp::simulate_cross_sphere<Vector3d, double>,
        "p"_a, "earth"_a, "Tcut"_a, "r0"_a, "R"_a, "rn"_a,
        "weight0"_a=1.0, "t_long_track"_a=1000000);
  m.def("simulate_cross_depth", &dp::simulate_cross_depth<Vector3d, double>,
        "p"_a, "earth"_a, "Tcut"_a, "depth"_a, "rn"_a,
        "weight0"_a=1.0, "t_long_track"_a=1000000);
  m.def("simulate_track", &dp::simulate_track<Vector3d, double>,
        "p"_a, "earth"_a, "Tcut"_a, "rn"_a,
        "weight0"_a=1.0, "long_track"_a=1000000);
  m.def("analyse_track", &dp::analyse_track<Vector3d, double>);
}
