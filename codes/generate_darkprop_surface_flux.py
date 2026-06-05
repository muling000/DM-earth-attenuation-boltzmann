from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import iteration_dirac_dm_isoscalar_vector_coupling as transport


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a DarkProp flux input from the benchmark halo surface spectrum "
            "used by the Dirac-DM benchmark iteration."
        )
    )
    parser.add_argument("--mchi-mev", type=float, default=5000.0)
    parser.add_argument("--sigma-chin-cm2", type=float, nargs="+", default=[1.0e-33])
    parser.add_argument("--n-energy", type=int, default=300)
    parser.add_argument("--v-grid-controller", type=float, default=30.0)
    parser.add_argument("--v-min-fraction", type=float, default=1.0e-20)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "external"
        / "darkprop"
        / "examples"
        / "surface-flux",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    m_chi = args.mchi_mev * transport.MEV
    t_grid = transport.build_quasi_uniform_speed_energy_grid(
        m_chi,
        args.n_energy,
        controller=args.v_grid_controller,
        v_min_fraction=args.v_min_fraction,
    )
    g_surface = transport.make_benchmark_surface_spectrum(m_chi)(t_grid)

    # The Boltzmann iteration evolves G = p^2 f. For an isotropic boundary condition,
    # dPhi/(dT dOmega) = G / (2*pi)^3. DarkProp wants GeV^-1 cm^-2 s^-1 sr^-1.
    flux_per_sr_per_gev = g_surface / (2.0 * np.pi) ** 3 * transport.GEV
    t_grid_gev = t_grid / transport.GEV

    flux_columns = [flux_per_sr_per_gev for _ in args.sigma_chin_cm2]
    flux_table = np.column_stack([t_grid_gev, *flux_columns])

    mchi_path = output_dir / "mchi.txt"
    sigma_path = output_dir / "sigma.txt"
    flux_path = output_dir / "surface_flux.txt"

    np.savetxt(mchi_path, np.array([m_chi / transport.GEV]), fmt="%.16e")
    np.savetxt(sigma_path, np.asarray(args.sigma_chin_cm2, dtype=np.float64), fmt="%.16e")
    np.savetxt(flux_path, flux_table, fmt="%.16e")

    total_intensity = float(np.trapezoid(flux_per_sr_per_gev, t_grid_gev))
    total_all_sky_flux = 4.0 * np.pi * total_intensity
    print(f"wrote {mchi_path}")
    print(f"wrote {sigma_path}")
    print(f"wrote {flux_path}")
    print(f"m_chi = {m_chi / transport.GEV:.8e} GeV")
    print(f"sigma_chiN columns = {len(args.sigma_chin_cm2)}")
    print(f"T range = [{t_grid_gev[0]:.8e}, {t_grid_gev[-1]:.8e}] GeV")
    print(f"integral dPhi/(dOmega) = {total_intensity:.8e} cm^-2 s^-1 sr^-1")
    print(f"integral all-sky dPhi = {total_all_sky_flux:.8e} cm^-2 s^-1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
