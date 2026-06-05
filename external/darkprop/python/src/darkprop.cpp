#include "darkprop.hpp"

void bind_Const(py::module_ &);
void bind_Vector3d(py::module_ &);
void bind_DMHalo(py::module_ &);
void bind_Base(py::module_ &);
void bind_Particle(py::module_ &);
void bind_Earth(py::module_ &);
void bind_Utils(py::module_ &);
void bind_Initialization(py::module_ &);
void bind_Simulation(py::module_ &);
void bind_Injector(py::module_ &);

PYBIND11_MODULE(_darkpropy, m)
{
  bind_Const(m);
  bind_Vector3d(m);
  bind_DMHalo(m);
  bind_Base(m);
  bind_Particle(m);
  bind_Earth(m);
  bind_Utils(m);
  bind_Initialization(m);
  bind_Simulation(m);
  bind_Injector(m);
}
