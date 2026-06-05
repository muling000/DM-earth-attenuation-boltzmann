/**
 * @file
 * Physical constants and units. Refer to PDG 2020.
 * https://pdg.lbl.gov/2021/reviews/contents_sports.html
 */

#pragma once
#include <cmath>

namespace darkprop {

namespace units {
    constexpr double GeV = 1.0;
    constexpr double eV = 1e-9 * GeV;
    constexpr double keV = 1e-6 * GeV;
    constexpr double MeV = 1e-3 * GeV;
    constexpr double TeV = 1e3 * GeV;
    constexpr double PeV = 1e6 * GeV;


    constexpr double fm = 1.0 / (197.3269804 * MeV);
    constexpr double mm = 1e12 * fm;
    constexpr double cm = 10.0 * mm;
    constexpr double cm2 = cm * cm;
    constexpr double cm3 = cm2 * cm;
    constexpr double m = 1e3 * mm;
    constexpr double km = 1e3 * m;
    constexpr double au = 149597870700 * m;
    constexpr double pc = 3.08567758149e16 * m;
    constexpr double kpc = 1e3 * pc;
    constexpr double Mpc = 1e6 * pc;

    constexpr double barn = 1e-24 * cm * cm;
    constexpr double pb = 1e-36 * cm * cm;

    constexpr double sec = 299792458 * m;
    constexpr double minute = 60 * sec;
    constexpr double hour = 60 * minute;

    constexpr double day = 24 * hour;
    // mean sidereal day
    // constexpr double day = 23 * hour + 56 * minute + 4.09053 * sec;
    // constexpr double year = 31558149.8 * sec;  ///< sidereal year
    constexpr double year = 31556925.1 * sec;  ///< tropical year


    // Planck constant, reduced (h bar) = 1.0545718e-34 J sec
    constexpr double kg = sec / (1.0545718e-34 * m * m);
    constexpr double gram = 1e-3 * kg;
    constexpr double g_cm3 = gram / (cm * cm * cm);

    constexpr double N = kg * m / sec / sec;     ///< newton
    constexpr double J = kg * m * m / sec / sec; ///< joule
    constexpr double F = m / 8.8541878128e-12;   ///< farad
    const double A = std::sqrt(1.00000000055 * 4*M_PI * 1e-7 * N);  ///< ampere
    const double C = A * sec;  ///< coulomb
} // namespace darkprop::units

namespace constants {
    constexpr double mp = 938.27208816 * units::MeV;  ///< proton mass
    constexpr double mn = 939.56542052 * units::MeV;  ///< neutron mass
    constexpr double me = 0.51099895 * units::MeV;    ///< electron mass
    constexpr double mu = 931.49410242 * units::MeV;  ///< atomic mass unit

    constexpr double alpha = 7.29735256931e-3;  ///< fine-structure constant
    constexpr double alpha2 = alpha * alpha;

    constexpr double rEarth = 6371.0 * units::km;
    constexpr double rSun = 69634000000 * units::cm;

    // dark matter halo parameters (they are not constants, though)
    constexpr double v0 = 220.0 * units::km / units::sec;
    constexpr double vesc = 544.0 * units::km / units::sec;
    constexpr double vearth = 240.0 * units::km / units::sec;
} // namespace darkprop::constants

} // namespace darkprop
