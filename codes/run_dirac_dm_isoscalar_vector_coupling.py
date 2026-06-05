from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import iteration_dirac_dm_isoscalar_vector_coupling as mod


PRESETS = {
    "smoke": {
        "n_u": 12,
        "n_t": 24,
        "n_x": 4,
        "n_y": 4,
        "n_z": 4,
    },
    "quick": {
        "n_u": 20,
        "n_t": 80,
        "n_x": 4,
        "n_y": 6,
        "n_z": 8,
    },
    "production": {
        "n_u": 40,
        "n_t": 300,
        "n_x": 4,
        "n_y": 8,
        "n_z": 8,
    },
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the Dirac-DM isoscalar-vector benchmark iteration and save detector spectra."
    )
    parser.add_argument("--preset", choices=tuple(PRESETS), default="quick")
    parser.add_argument("--mchi-mev", type=float, default=100.0, help="Dark-matter mass in MeV.")
    parser.add_argument(
        "--sigma-chin-cm2",
        type=float,
        default=1.0e-33,
        help="Non-relativistic DM-nucleon cross section in cm^2.",
    )
    parser.add_argument("--detector-depth-km", type=float, default=2.4)
    parser.add_argument("--n-tail", type=int, default=16, help="Radial tail points from detector_depth + tail_tau*lambda to the Earth center.")
    parser.add_argument(
        "--main-tau-step",
        type=float,
        default=0.2,
        help="Maximum radial spacing in units of lambda up to detector_depth + tail_tau*lambda.",
    )
    parser.add_argument(
        "--tail-tau",
        type=float,
        default=8.0,
        help="Resolve the radial grid uniformly in optical depth through detector_depth + tail_tau*lambda.",
    )
    parser.add_argument(
        "--max-main-step-km",
        type=float,
        default=50.0,
        help="Absolute cap on the fine-region radial spacing, used when lambda is large.",
    )
    parser.add_argument("--tail-power", type=float, default=1.5)
    parser.add_argument("--n-u", type=int, default=None, help="Override preset angular grid points.")
    parser.add_argument("--n-t", type=int, default=None, help="Override preset kinetic-energy grid points.")
    parser.add_argument("--n-x", type=int, default=None, help="Override preset path quadrature order.")
    parser.add_argument("--n-y", type=int, default=None, help="Override preset energy-loss quadrature order.")
    parser.add_argument("--n-z", type=int, default=None, help="Override preset azimuthal quadrature order.")
    parser.add_argument(
        "--tau-segment-edges",
        default="0,0.125,0.25,0.5,1,2,4,8,16,24",
        help="Comma-separated optical-depth segment edges for the exponentially weighted path integral.",
    )
    parser.add_argument(
        "--energy-grid-mode",
        choices=("quasi-uniform-v", "log"),
        default="quasi-uniform-v",
        help="Energy-grid construction: 'quasi-uniform-v' builds an almost uniform speed grid and maps it to T.",
    )
    parser.add_argument(
        "--t-min-fraction",
        type=float,
        default=1.0e-3,
        help="Lowest kinetic energy as a fraction of T_max. Used by the log energy grid.",
    )
    parser.add_argument(
        "--v-grid-controller",
        type=float,
        default=30.0,
        help="Shape parameter for the quasi-uniform-v energy grid. Larger values are closer to uniform spacing in speed.",
    )
    parser.add_argument(
        "--v-min-fraction",
        type=float,
        default=1.0e-20,
        help="Minimum speed as a fraction of vmax in quasi-uniform-v mode.",
    )
    parser.add_argument("--eps", type=float, default=1.0e-4)
    parser.add_argument(
        "--stop-energy-fraction",
        type=float,
        default=1.0e-3,
        help="Only detector-spectrum bins with T >= fraction * T_max enter the stopping metric.",
    )
    parser.add_argument(
        "--stop-consecutive",
        type=int,
        default=3,
        help="Stop only after the metric stays below eps for this many consecutive orders.",
    )
    parser.add_argument(
        "--stop-floor-fraction",
        type=float,
        default=1.0e-12,
        help="Detector-spectrum floor as a fraction of the peak zeroth-order detector spectrum.",
    )
    parser.add_argument("--max-orders", type=int, default=2000)
    parser.add_argument("--strict", action="store_true", help="Enable strict grid-coverage checks.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output .npz path. Defaults to codes/output/...",
    )
    return parser


def build_runtime_configuration(args: argparse.Namespace) -> dict[str, float | int | str]:
    preset = PRESETS[args.preset]
    config: dict[str, float | int | str] = {
        "preset": args.preset,
        "mchi_mev": args.mchi_mev,
        "sigma_chin_cm2": args.sigma_chin_cm2,
        "detector_depth_km": args.detector_depth_km,
        "n_tail": args.n_tail,
        "main_tau_step": args.main_tau_step,
        "tail_tau": args.tail_tau,
        "max_main_step_km": args.max_main_step_km,
        "tail_power": args.tail_power,
        "energy_grid_mode": args.energy_grid_mode,
        "t_min_fraction": args.t_min_fraction,
        "v_grid_controller": args.v_grid_controller,
        "v_min_fraction": args.v_min_fraction,
        "tau_segment_edges": args.tau_segment_edges,
        "eps": args.eps,
        "stop_energy_fraction": args.stop_energy_fraction,
        "stop_consecutive": args.stop_consecutive,
        "stop_floor_fraction": args.stop_floor_fraction,
        "max_orders": args.max_orders,
        "strict": int(args.strict),
    }
    config.update(preset)
    for key in ("n_u", "n_t", "n_x", "n_y", "n_z"):
        value = getattr(args, key)
        if value is not None:
            config[key] = value
    return config


def parse_tau_segment_edges(value: str) -> list[float]:
    edges = [float(part) for part in value.split(",") if part.strip()]
    if len(edges) < 2:
        raise ValueError("--tau-segment-edges must contain at least two comma-separated values.")
    return edges


def default_output_path(args: argparse.Namespace) -> Path:
    output_dir = Path(__file__).resolve().parent / "output"
    mchi_str = f"{args.mchi_mev:.6g}"
    sigma_str = f"{args.sigma_chin_cm2:.3e}"
    if args.energy_grid_mode == "quasi-uniform-v":
        grid_token = f"{args.energy_grid_mode}_vmin{float_token(args.v_min_fraction)}"
    else:
        grid_token = f"{args.energy_grid_mode}_tmin{float_token(args.t_min_fraction)}"
    filename = (
        f"dirac_dm_vector_mchiMeV_{mchi_str}_sigmaChiN_{sigma_str}"
        f"_{grid_token}_{args.preset}_radial_default.npz"
    )
    return output_dir / filename


def float_token(value: float) -> str:
    return f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e").replace("e+", "e")


def build_solver_and_geometry(
    config: dict[str, float | int | str],
) -> tuple[mod.BenchmarkIterationSolver, float]:
    m_chi = float(config["mchi_mev"]) * mod.MEV
    sigma_chi_n = float(config["sigma_chin_cm2"])
    detector_depth = float(config["detector_depth_km"]) * mod.KM

    model = mod.make_benchmark_model_from_sigma(
        m_chi=m_chi,
        sigma_chi_n=sigma_chi_n,
        sigma_in_cm2=True,
    )
    detector_radius = model.earth_radius - detector_depth

    r_grid = mod.build_radial_grid(
        model.earth_radius,
        mean_free_path=mod.mean_free_path_for_model(model),
        detector_depth=detector_depth,
        main_tau_step=float(config["main_tau_step"]),
        tail_tau=float(config["tail_tau"]),
        n_tail=int(config["n_tail"]),
        tail_power=float(config["tail_power"]),
        max_main_step=float(config["max_main_step_km"]) * mod.KM,
    )
    u_grid = mod.build_clenshaw_curtis_u_grid(int(config["n_u"]))
    energy_grid_mode = str(config["energy_grid_mode"])
    if energy_grid_mode == "quasi-uniform-v":
        t_grid = mod.build_quasi_uniform_speed_energy_grid(
            model.m_chi,
            int(config["n_t"]),
            controller=float(config["v_grid_controller"]),
            v_min_fraction=float(config["v_min_fraction"]),
        )
    else:
        t_grid = mod.build_log_energy_grid(
            model.m_chi,
            int(config["n_t"]),
            t_min_fraction=float(config["t_min_fraction"]),
        )
    grid = mod.PhaseSpaceGrid(r=r_grid, u=u_grid, T=t_grid)
    quadrature = mod.QuadratureSpec(
        n_x=int(config["n_x"]),
        n_y=int(config["n_y"]),
        n_z=int(config["n_z"]),
    )
    solver = mod.BenchmarkIterationSolver(
        model=model,
        grid=grid,
        quadrature=quadrature,
        strict=bool(config["strict"]),
        tau_segment_edges=parse_tau_segment_edges(str(config["tau_segment_edges"])),
    )
    return solver, detector_radius


def high_energy_mask(kinetic_grid: np.ndarray, fraction: float) -> np.ndarray:
    threshold = fraction * kinetic_grid[-1]
    mask = kinetic_grid >= threshold
    if not np.any(mask):
        mask[-1] = True
    return mask


def stopping_metric(
    term_spectrum: np.ndarray,
    cumulative_spectrum: np.ndarray,
    monitor_mask: np.ndarray,
    floor_value: float,
) -> float:
    numerator = np.abs(term_spectrum[monitor_mask])
    denominator = np.maximum(np.abs(cumulative_spectrum[monitor_mask]), floor_value)
    return float(np.max(numerator / denominator))


def run_iteration(
    solver: mod.BenchmarkIterationSolver,
    detector_radius: float,
    config: dict[str, float | int | str],
) -> dict[str, np.ndarray | float | int]:
    surface_spectrum = mod.make_benchmark_surface_spectrum(solver.model.m_chi)
    term_field = solver.initial_term(surface_spectrum)
    term_detector = mod.detector_spectrum(term_field, solver.grid.r, solver.grid.u, detector_radius)

    term_spectra = [term_detector.copy()]
    cumulative_detector = term_detector.copy()
    cumulative_spectra = [cumulative_detector.copy()]
    metric_history = [np.nan]
    order_time_sec = [0.0]

    monitor_mask = high_energy_mask(solver.grid.T, float(config["stop_energy_fraction"]))
    floor_value = max(
        float(config["stop_floor_fraction"]) * float(np.max(np.abs(term_detector))),
        1.0e-300,
    )
    below_eps_count = 0
    completed_orders = 0

    print(
        "order=0"
        f" detector_peak={np.max(np.abs(term_detector)):.3e}"
        f" monitored_bins={int(np.count_nonzero(monitor_mask))}"
    )

    for order in range(1, int(config["max_orders"]) + 1):
        start = time.perf_counter()
        term_field = solver.apply_iteration(term_field)
        elapsed = time.perf_counter() - start

        term_detector = mod.detector_spectrum(term_field, solver.grid.r, solver.grid.u, detector_radius)
        cumulative_detector = cumulative_detector + term_detector
        metric = stopping_metric(term_detector, cumulative_detector, monitor_mask, floor_value)

        term_spectra.append(term_detector.copy())
        cumulative_spectra.append(cumulative_detector.copy())
        metric_history.append(metric)
        order_time_sec.append(elapsed)
        completed_orders = order

        print(
            f"order={order}"
            f" time={elapsed:.2f}s"
            f" detector_peak={np.max(np.abs(term_detector)):.3e}"
            f" stop_metric={metric:.3e}"
        )

        if metric < float(config["eps"]):
            below_eps_count += 1
        else:
            below_eps_count = 0

        if below_eps_count >= int(config["stop_consecutive"]):
            print(
                f"stopping: metric < eps for {below_eps_count} consecutive orders "
                f"(eps={float(config['eps']):.3e})"
            )
            break

    return {
        "term_spectra": np.stack(term_spectra, axis=0),
        "cumulative_spectra": np.stack(cumulative_spectra, axis=0),
        "metric_history": np.asarray(metric_history, dtype=np.float64),
        "order_time_sec": np.asarray(order_time_sec, dtype=np.float64),
        "monitor_mask": np.asarray(monitor_mask, dtype=bool),
        "floor_value": float(floor_value),
        "completed_orders": int(completed_orders),
    }


def save_output(
    output_path: Path,
    solver: mod.BenchmarkIterationSolver,
    detector_radius: float,
    config: dict[str, float | int | str],
    run_data: dict[str, np.ndarray | float | int],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    detector_depth = solver.model.earth_radius - detector_radius
    velocity_grid = np.sqrt(2.0 * solver.grid.T / solver.model.m_chi)
    incoming_surface = mod.make_benchmark_surface_spectrum(solver.model.m_chi)(solver.grid.T)

    np.savez_compressed(
        output_path,
        config_json=json.dumps(config, sort_keys=True),
        r_grid=solver.grid.r,
        u_grid=solver.grid.u,
        T_grid=solver.grid.T,
        velocity_grid=velocity_grid,
        detector_radius=detector_radius,
        detector_depth=detector_depth,
        incoming_surface_spectrum=incoming_surface,
        mean_free_path=solver.mean_free_path,
        cutoff_scale=solver.model.cutoff_scale,
        term_spectra=run_data["term_spectra"],
        cumulative_spectra=run_data["cumulative_spectra"],
        metric_history=run_data["metric_history"],
        order_time_sec=run_data["order_time_sec"],
        monitor_mask=run_data["monitor_mask"],
        floor_value=run_data["floor_value"],
        completed_orders=run_data["completed_orders"],
    )


def print_model_diagnostics(
    solver: mod.BenchmarkIterationSolver,
    config: dict[str, float | int | str],
) -> None:
    sigma_chi_n_cm2 = float(config["sigma_chin_cm2"])
    cutoff_scale_gev = solver.model.cutoff_scale / mod.GEV
    mean_free_path_km = solver.mean_free_path / mod.KM
    mean_free_path_cm = solver.mean_free_path / mod.CM

    print(
        "model diagnostics:"
        f" m_chi={solver.model.m_chi / mod.MEV:.6g} MeV"
        f" sigma_chiN={sigma_chi_n_cm2:.6e} cm^2"
        f" cutoff_scale={cutoff_scale_gev:.6e} GeV"
        f" mean_free_path={mean_free_path_km:.6e} km ({mean_free_path_cm:.6e} cm)"
    )
    print("species diagnostics:")
    for species in solver.model.species:
        sigma_a_natural = solver._sigma_nr(species)
        sigma_a_cm2 = mod.sigma_natural_to_cm2(sigma_a_natural)
        inv_mfp_a = species.number_density * sigma_a_natural
        lambda_a_km = (1.0 / inv_mfp_a) / mod.KM
        print(
            f"  {species.name:>2s}"
            f" A={species.mass_number:>2d}"
            f" sigma_chiA={sigma_a_cm2:.6e} cm^2"
            f" lambda_A={lambda_a_km:.6e} km"
        )


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    config = build_runtime_configuration(args)
    output_path = args.output or default_output_path(args)

    solver, detector_radius = build_solver_and_geometry(config)
    print_model_diagnostics(solver, config)
    run_data = run_iteration(solver, detector_radius, config)
    save_output(output_path, solver, detector_radius, config, run_data)

    print(f"saved output to {output_path}")
    print(
        f"grid: Nr={solver.grid.r.size} Nu={solver.grid.u.size} Nt={solver.grid.T.size}; "
        f"quadrature: Nx={solver.quadrature.n_x} Ny={solver.quadrature.n_y} Nz={solver.quadrature.n_z}; "
        f"path=tau-segmented path_nodes={solver._omega_x.shape[2]}"
    )
    print(
        f"detector depth = {(solver.model.earth_radius - detector_radius) / mod.KM:.3f} km; "
        f"completed orders = {int(run_data['completed_orders'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
