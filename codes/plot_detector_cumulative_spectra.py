from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import iteration_dirac_dm_isoscalar_vector_coupling as mod


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot cumulative detector spectra from a run_dirac_dm_isoscalar_vector_coupling output .npz."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input .npz file produced by run_dirac_dm_isoscalar_vector_coupling.py.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output image path. For --x-axis both, _speed and _energy are appended to the stem.",
    )
    parser.add_argument(
        "--x-axis",
        choices=("speed", "energy", "both"),
        default="both",
        help="Horizontal axis for the plot.",
    )
    parser.add_argument(
        "--max-curves",
        type=int,
        default=0,
        help="Maximum number of cumulative-order curves to draw. Use 0 to draw every order.",
    )
    parser.add_argument(
        "--include-final",
        action="store_true",
        help="Always include the final cumulative curve when downsampling orders.",
    )
    parser.add_argument(
        "--show-monitor-threshold",
        action="store_true",
        help="Draw a vertical line at the first monitored high-energy bin if the file contains monitor_mask.",
    )
    return parser.parse_args()


def default_output_path(input_path: Path, x_axis: str) -> Path:
    suffix = f"cumulative_{x_axis}"
    return input_path.with_name(f"{input_path.stem}_{suffix}.png")


def output_path_for_axis(output: Path | None, input_path: Path, x_axis: str) -> Path:
    if output is None:
        return default_output_path(input_path, x_axis)
    if x_axis in ("speed", "energy"):
        return output
    return output.with_name(f"{output.stem}_{x_axis}{output.suffix or '.png'}")


def choose_order_indices(n_orders: int, max_curves: int, include_final: bool) -> np.ndarray:
    if max_curves <= 0 or n_orders <= max_curves:
        return np.arange(n_orders, dtype=np.int64)

    indices = np.linspace(0, n_orders - 1, num=max_curves, dtype=np.int64)
    indices = np.unique(indices)
    if include_final and indices[-1] != n_orders - 1:
        indices = np.unique(np.concatenate([indices, np.array([n_orders - 1], dtype=np.int64)]))
    return indices


def infer_dm_mass_from_grids(t_grid: np.ndarray, v_grid: np.ndarray) -> float:
    positive = v_grid > 0.0
    if not np.any(positive):
        raise ValueError("velocity_grid must contain positive values.")
    inferred = 2.0 * t_grid[positive] / np.square(v_grid[positive])
    return float(np.median(inferred))


def convert_spectrum_for_axis(
    cumulative: np.ndarray,
    t_grid: np.ndarray,
    v_grid: np.ndarray,
    x_axis: str,
) -> tuple[np.ndarray, np.ndarray, str, str]:
    flux_prefactor = 1.0 / (2.0 * np.pi * np.pi)

    if x_axis == "energy":
        return (
            t_grid,
            cumulative * flux_prefactor,
            "Kinetic Energy T [eV]",
            r"Cumulative Flux $d\Phi/dT$ [cm$^{-2}$ s$^{-1}$ eV$^{-1}$]",
        )

    m_chi = infer_dm_mass_from_grids(t_grid, v_grid)
    speed_density = cumulative * flux_prefactor * (m_chi * v_grid[np.newaxis, :])
    return (
        v_grid,
        speed_density,
        "Speed v/c",
        r"Cumulative Flux $d\Phi/dv$ [cm$^{-2}$ s$^{-1}$]",
    )


def plot_single_axis(
    *,
    data: np.lib.npyio.NpzFile,
    x_axis: str,
    output_path: Path,
    order_indices: np.ndarray,
    show_monitor_threshold: bool,
) -> None:
    cumulative = np.asarray(data["cumulative_spectra"], dtype=np.float64)
    t_grid = np.asarray(data["T_grid"], dtype=np.float64)
    v_grid = np.asarray(data["velocity_grid"], dtype=np.float64)
    monitor_mask = np.asarray(data["monitor_mask"], dtype=bool) if "monitor_mask" in data.files else None

    x_values, plotted_spectra, x_label, y_label = convert_spectrum_for_axis(cumulative, t_grid, v_grid, x_axis)
    final_spectrum = plotted_spectra[-1]

    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")

    for rank, order in enumerate(order_indices):
        color = cmap(rank / max(len(order_indices) - 1, 1))
        ax.plot(x_values, plotted_spectra[order], color=color, linewidth=1.0, alpha=0.9)

    first_idx = None
    if monitor_mask is not None and np.any(monitor_mask):
        first_idx = int(np.argmax(monitor_mask))
        y_limit = 1.5 * float(final_spectrum[first_idx])
        if y_limit > 0.0:
            ax.set_ylim(0.0, y_limit)

    if first_idx is not None and show_monitor_threshold:
        ax.axvline(
            x_values[first_idx],
            color="black",
            linewidth=1.0,
            linestyle="--",
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title("Detector Flux by Cumulative Scattering Order")
    ax.grid(True, which="both", alpha=0.25)

    n_orders = cumulative.shape[0]
    completed_orders = int(np.asarray(data["completed_orders"]).item()) if "completed_orders" in data.files else n_orders - 1
    detector_depth_km = None
    if "detector_depth" in data.files:
        detector_depth_km = float(np.asarray(data["detector_depth"]).item()) / mod.KM

    summary_parts = [f"curves={len(order_indices)}", f"completed orders={completed_orders}"]
    if detector_depth_km is not None:
        summary_parts.append(f"depth={detector_depth_km:.3f} km")
    ax.text(
        0.98,
        0.98,
        "\n".join(summary_parts),
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        ha="right",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    print(f"saved plot to {output_path}")


def main() -> int:
    args = parse_args()
    data = np.load(args.input)

    n_orders = int(np.asarray(data["cumulative_spectra"]).shape[0])
    order_indices = choose_order_indices(n_orders, args.max_curves, args.include_final)

    if args.x_axis == "both":
        for axis_name in ("speed", "energy"):
            plot_single_axis(
                data=data,
                x_axis=axis_name,
                output_path=output_path_for_axis(args.output, args.input, axis_name),
                order_indices=order_indices,
                show_monitor_threshold=args.show_monitor_threshold,
            )
    else:
        plot_single_axis(
            data=data,
            x_axis=args.x_axis,
            output_path=output_path_for_axis(args.output, args.input, args.x_axis),
            order_indices=order_indices,
            show_monitor_threshold=args.show_monitor_threshold,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
