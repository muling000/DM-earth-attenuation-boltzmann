#include "darkprop.hpp"

void bind_Const(py::module_ &mod)
{
  using namespace dp::units;
  auto units = mod.def_submodule("units", "Units");
  units.attr("GeV") = GeV;
  units.attr("eV") = eV;
  units.attr("keV") = keV;
  units.attr("MeV") = MeV;
  units.attr("TeV") = TeV;
  units.attr("PeV") = PeV;

  units.attr("fm") = fm;
  units.attr("mm") = mm;
  units.attr("cm") = cm;
  units.attr("cm2") = cm2;
  units.attr("cm3") = cm3;
  units.attr("m") = m;
  units.attr("au") = au;
  units.attr("km") = km;
  units.attr("pc") = pc;
  units.attr("kpc") = kpc;
  units.attr("Mpc") = Mpc;

  units.attr("barn") = barn;
  units.attr("pb") = pb;

  units.attr("sec") = sec;
  units.attr("s") = sec;
  units.attr("minute") = minute;
  units.attr("hour") = hour;
  units.attr("day") = day;
  units.attr("year") = year;
  units.attr("yr") = year;

  units.attr("kg") = kg;
  units.attr("gram") = gram;
  units.attr("g") = gram;
  units.attr("g_cm3") = g_cm3;

  units.attr("N") = N;
  units.attr("J") = J;
  units.attr("F") = F;
  units.attr("A") = A;
  units.attr("C") = C;

  using namespace dp::constants;
  auto constants = mod.def_submodule("constants", "Physical constants");
  constants.attr("mp") = mp;
  constants.attr("mn") = mn;
  constants.attr("me") = me;
  constants.attr("mu") = mu;

  constants.attr("alpha") = alpha;
  constants.attr("alpha2") = alpha2;

  constants.attr("rEarth") = rEarth;
  constants.attr("rSun") = rSun;

  constants.attr("v0") = v0;
  constants.attr("vesc") = vesc;
  constants.attr("vearth") = vearth;
}
