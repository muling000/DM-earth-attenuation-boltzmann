/*
 * Simulate 50 trajectories.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include "darkprop/core.hpp"

int main(void)
{
    using namespace darkprop::constants;  // using rEarth
    using namespace darkprop::units;      // using cm, GeV, sec

    double sigma = 1e-32 * cm * cm;  // cross section
    double mchi = 10 * GeV;          // DM mass
    double t = 0.0;                  // initial time

    double vcut = 1.0 * cm / sec;               // cutoff velocity
    double Tcut = 0.5 * mchi * vcut * vcut;     // cutoff kinetic energy

    // initial position (constants::rEarth is defined in Const.hpp)
    darkprop::Vector3d init_r(0.0, 0.0, rEarth);

    // initial velocity
    double vi = 220.0 * km / sec;
    darkprop::Vector3d init_v(0.0, 0.0, -vi);

    // instantiating the DM particle and the Earth
    darkprop::SIDM dm(sigma, mchi, t);
    darkprop::HomoEarth earth;

    // need a random number generator
    darkprop::RandomNumber rn;

    std::filesystem::create_directory("out");

    // simulate 50 tracks
    std::size_t number_of_tracks = 50;
    for (std::size_t i = 0; i < number_of_tracks; ++i) {
        dm.t = 0.0;           // time is not used in this case
        dm.setR(init_r);      // reset position
        dm.setV(init_v);      // reset velocity
        dm.in_medium = true;  // reset the flag

        // simulate full trajectory
        auto events = simulate_track(dm, earth, Tcut, rn);

        // output: only store the position of all events in csv files
        std::string filename = "out/trajectory-";
        filename = filename + std::to_string(i) + ".csv";
        std::ofstream off(filename);
        off << std::scientific;
        off.precision(16);
        for (const auto& event : events) {
            off << event.r[0] << " "
                << event.r[1] << " "
                << event.r[2] << std::endl;
        }
    }
    return 0;
}
