# Boltzmann Benchmark

&emsp;&emsp;This fork adds an optional `BoltzmannBenchmark` Earth model for comparison with the deterministic Boltzmann-equation calculation.

## Changes from upstream

- Added the optional configuration setting `earthmodel = "BoltzmannBenchmark"`.
- Added a homogeneous Earth with density `2.7 g/cm^3` and the same eight unnormalized elemental mass fractions as the Boltzmann benchmark.
- Fixed the Earth speed to `240 km/s` in benchmark mode.
- Sampled free paths directly from the homogeneous exponential distribution.
- Updated the `cm`, `gram`, and nucleon-mass constants to the benchmark values.
- Added `bin/boltzmann_benchmark.cfg` for `mchi = 5 GeV`, `sigma_chiN = 5e-32 cm^2`, and a detector depth of `2.4 km`.

## Unchanged behavior

&emsp;&emsp;The original PREM model remains the default, and its tables and propagation algorithm are unchanged. The global constant updates cause a sub-percent numerical shift in legacy PREM runs. Random seeding, nonrelativistic spin-independent collision kinematics, the Analyzer, and native output formats are unchanged.
