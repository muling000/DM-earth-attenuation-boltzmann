from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


C_KM_S = 299792.458
C_CM_S = 2.99792458e10
MCHI_GEV = 5.0


def parse_args() -> argparse.Namespace:
    code_dir = Path(__file__).resolve().parent
    repository_root = code_dir.parent
    parser = argparse.ArgumentParser(
        description="Compare DaMaSCUS, deterministic Boltzmann, and DarkProp dPhi/dv."
    )
    parser.add_argument(
        "--damascus-data",
        type=Path,
        default=repository_root / "external" / "damascus",
    )
    parser.add_argument(
        "--boltzmann",
        type=Path,
        default=code_dir
        / "output"
        / (
            "dirac_dm_vector_mchiMeV_5000_sigmaChiN_5.000e-32_"
            "quasi-uniform-v_vmin1e-20_production_radial_default.npz"
        ),
    )
    parser.add_argument(
        "--darkprop",
        type=Path,
        default=code_dir
        / "output"
        / "darkprop_mc_dphidv_vmin1kms_fine1000.dat",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=code_dir / "output" / "damascus_boltzmann_dphidv.pdf",
    )
    parser.add_argument("--ratio-bins", type=int, default=14)
    return parser.parse_args()


def load_damascus(
    root: Path,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    velocity = np.fromfile(root / "data" / "velocity.0", dtype=np.float64)
    if velocity.size % 3:
        raise ValueError("DaMaSCUS velocity file does not contain complete 3-vectors.")
    speed = np.linalg.norm(velocity.reshape(-1, 3), axis=1)
    weight = np.fromfile(root / "data" / "weights.0", dtype=np.float64)
    if speed.size != weight.size:
        raise ValueError("DaMaSCUS velocity and weight sample counts differ.")

    rho_row = np.atleast_2d(np.loadtxt(root / "results" / "density.rho"))[0]
    rho_gev_cm3 = float(rho_row[1])
    rho_error_gev_cm3 = float(rho_row[2])
    return speed, weight, rho_gev_cm3, rho_error_gev_cm3


def bin_damascus_dphidv(
    speed: np.ndarray,
    weight: np.ndarray,
    rho_gev_cm3: float,
    rho_error_gev_cm3: float,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    number_density = rho_gev_cm3 / MCHI_GEV
    weight_sum = np.sum(weight)
    weighted_speed, _ = np.histogram(speed, bins=edges, weights=weight * speed)
    weighted_speed_sq, _ = np.histogram(
        speed, bins=edges, weights=np.square(weight * speed)
    )
    widths_beta = np.diff(edges)
    centers_beta = np.sqrt(edges[:-1] * edges[1:])

    # DaMaSCUS samples the local number-density distribution. Multiplication
    # by the particle speed converts it to the scalar differential flux.
    dphidbeta = number_density * C_CM_S * weighted_speed / weight_sum / widths_beta
    stat_error = (
        number_density
        * C_CM_S
        * np.sqrt(weighted_speed_sq)
        / weight_sum
        / widths_beta
    )
    density_error = dphidbeta * rho_error_gev_cm3 / rho_gev_cm3
    error_dphidbeta = np.hypot(stat_error, density_error)

    return (
        centers_beta * C_KM_S,
        dphidbeta / C_KM_S,
        error_dphidbeta / C_KM_S,
    )


def load_boltzmann(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, float, int, dict[str, object]]:
    with np.load(path) as data:
        velocity_beta = np.asarray(data["velocity_grid"], dtype=np.float64)
        kinetic_energy = np.asarray(data["T_grid"], dtype=np.float64)
        cumulative = np.asarray(data["cumulative_spectra"], dtype=np.float64)[-1]
        monitor_mask = np.asarray(data["monitor_mask"], dtype=bool)
        completed_orders = int(np.asarray(data["completed_orders"]).item())
        config = json.loads(str(data["config_json"]))

    positive = velocity_beta > 0.0
    mchi_ev = float(
        np.median(2.0 * kinetic_energy[positive] / np.square(velocity_beta[positive]))
    )
    dphidbeta = cumulative * mchi_ev * velocity_beta / (2.0 * np.pi * np.pi)
    convergence_speed = float(velocity_beta[monitor_mask][0]) * C_KM_S
    return (
        velocity_beta * C_KM_S,
        dphidbeta / C_KM_S,
        convergence_speed,
        completed_orders,
        config,
    )


def load_darkprop(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(path)
    return (
        table[:, 0] * C_KM_S,
        table[:, 1] / C_KM_S,
        table[:, 2] / C_KM_S,
    )


def log_interpolate(
    x: np.ndarray,
    y: np.ndarray,
    x_new: np.ndarray,
) -> np.ndarray:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    output = np.full_like(x_new, np.nan, dtype=np.float64)
    valid = (x_new >= np.min(x[mask])) & (x_new <= np.max(x[mask]))
    output[valid] = np.exp(
        np.interp(np.log(x_new[valid]), np.log(x[mask]), np.log(y[mask]))
    )
    return output


def integrate_interval(
    speed: np.ndarray,
    dphidv: np.ndarray,
    lower: float,
    upper: float,
) -> float:
    mask = (speed > lower) & (speed < upper)
    x = np.concatenate(([lower], speed[mask], [upper]))
    y = np.concatenate(
        (
            [np.interp(lower, speed, dphidv)],
            dphidv[mask],
            [np.interp(upper, speed, dphidv)],
        )
    )
    return float(np.trapezoid(y, x))


def integrate_histogram_interval(
    centers: np.ndarray,
    heights: np.ndarray,
    errors: np.ndarray,
    lower: float,
    upper: float,
) -> tuple[float, float]:
    if centers.size < 2:
        raise ValueError("Histogram integration requires at least two bins.")
    width = float(np.median(np.diff(centers)))
    left = centers - 0.5 * width
    right = centers + 0.5 * width
    overlap = np.maximum(0.0, np.minimum(right, upper) - np.maximum(left, lower))
    integral = float(np.sum(heights * overlap))
    error = float(np.sqrt(np.sum(np.square(errors * overlap))))
    return integral, error


def weighted_mean_standard_error(weight: np.ndarray, values: np.ndarray) -> float:
    sample_size = weight.size
    weight_sum = np.sum(weight)
    mean_weight = weight_sum / sample_size
    mean_value = np.sum(weight * values) / weight_sum
    centered_product = weight * values - mean_weight * mean_value
    centered_weight = weight - mean_weight
    sum1 = np.sum(np.square(centered_product))
    sum2 = np.sum(centered_weight * centered_product)
    sum3 = np.sum(np.square(centered_weight))
    variance = (
        sample_size
        / (sample_size - 1)
        / weight_sum**2
        * (
            sum1
            - 2.0 * mean_value * sum2
            + mean_value**2 * sum3
        )
    )
    return float(np.sqrt(max(variance, 0.0)))


def main() -> int:
    args = parse_args()
    speed_beta, weight, rho, rho_error = load_damascus(args.damascus_data)
    (
        boltzmann_speed,
        boltzmann_dphidv,
        convergence_speed,
        completed_orders,
        config,
    ) = load_boltzmann(args.boltzmann)
    darkprop_speed, darkprop_dphidv, darkprop_error = load_darkprop(args.darkprop)

    comparison_vmax = float(np.max(boltzmann_speed))
    ratio_edges = np.geomspace(
        convergence_speed / C_KM_S,
        comparison_vmax / C_KM_S,
        args.ratio_bins + 1,
    )
    ratio_speed, ratio_damascus, ratio_damascus_error = bin_damascus_dphidv(
        speed_beta,
        weight,
        rho,
        rho_error,
        ratio_edges,
    )
    ratio_boltzmann = log_interpolate(
        boltzmann_speed,
        boltzmann_dphidv,
        ratio_speed,
    )
    damascus_ratio = ratio_damascus / ratio_boltzmann
    damascus_ratio_error = ratio_damascus_error / ratio_boltzmann

    darkprop_boltzmann = log_interpolate(
        boltzmann_speed,
        boltzmann_dphidv,
        darkprop_speed,
    )
    darkprop_ratio = darkprop_dphidv / darkprop_boltzmann
    darkprop_ratio_error = darkprop_error / darkprop_boltzmann

    number_density = rho / MCHI_GEV
    high_values = speed_beta * (
        (speed_beta >= convergence_speed / C_KM_S)
        & (speed_beta <= comparison_vmax / C_KM_S)
    )
    high_mean = np.sum(weight * high_values) / np.sum(weight)
    high_mean_error = weighted_mean_standard_error(weight, high_values)
    damascus_high_flux = number_density * C_CM_S * high_mean
    damascus_high_error = damascus_high_flux * np.hypot(
        high_mean_error / high_mean,
        rho_error / rho,
    )
    boltzmann_high_flux = integrate_interval(
        boltzmann_speed,
        boltzmann_dphidv,
        convergence_speed,
        comparison_vmax,
    )
    darkprop_high_flux, darkprop_high_error = integrate_histogram_interval(
        darkprop_speed,
        darkprop_dphidv,
        darkprop_error,
        convergence_speed,
        comparison_vmax,
    )
    damascus_in_range_samples = np.count_nonzero(
        (speed_beta >= convergence_speed / C_KM_S)
        & (speed_beta <= comparison_vmax / C_KM_S)
    )

    fig, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(8.2, 6.8),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.15]},
        constrained_layout=True,
    )

    positive_boltzmann = boltzmann_dphidv > 0.0
    axis.plot(
        boltzmann_speed[positive_boltzmann],
        boltzmann_dphidv[positive_boltzmann],
        color="black",
        linewidth=1.8,
        label="Boltzmann iteration",
        zorder=4,
    )

    positive_darkprop = (
        (darkprop_dphidv > 0.0)
        & (darkprop_speed >= convergence_speed)
        & (darkprop_speed <= comparison_vmax)
    )
    axis.errorbar(
        darkprop_speed[positive_darkprop],
        darkprop_dphidv[positive_darkprop],
        yerr=np.vstack(
            (
                np.minimum(
                    darkprop_error[positive_darkprop],
                    0.8 * darkprop_dphidv[positive_darkprop],
                ),
                darkprop_error[positive_darkprop],
            )
        ),
        fmt="o",
        markersize=2.8,
        markeredgewidth=0.0,
        linewidth=0.7,
        capsize=1.2,
        color="#2673b8",
        alpha=0.9,
        label="DarkProp MC (1.05M in range)",
        zorder=3,
    )

    positive_damascus = ratio_damascus > 0.0
    axis.errorbar(
        ratio_speed[positive_damascus],
        ratio_damascus[positive_damascus],
        yerr=np.vstack(
            (
                np.minimum(
                    ratio_damascus_error[positive_damascus],
                    0.8 * ratio_damascus[positive_damascus],
                ),
                ratio_damascus_error[positive_damascus],
            )
        ),
        fmt="s",
        markersize=3.0,
        markeredgewidth=0.0,
        linewidth=0.8,
        capsize=1.2,
        color="#d95f02",
        alpha=0.9,
        label=(
            f"DaMaSCUS MC "
            f"({damascus_in_range_samples / 1.0e6:.2f}M in range)"
        ),
        zorder=2,
    )

    axis.text(
        0.02,
        0.97,
        (
            rf"Comparison range: $v\geq {convergence_speed:.1f}\,\mathrm{{km/s}}$"
            "\n(Boltzmann stopping criterion enforced)"
        ),
        transform=axis.transAxes,
        fontsize=8.5,
        color="0.3",
        va="top",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.8,
        },
    )
    axis.text(
        0.98,
        0.96,
        (
            rf"$\Phi(v\geq {convergence_speed:.1f}\,\mathrm{{km/s}})"
            rf"/\Phi_{{\rm Boltz}}$"
            "\n"
            rf"DaMaSCUS: ${damascus_high_flux / boltzmann_high_flux:.3f}"
            rf"\pm{damascus_high_error / boltzmann_high_flux:.3f}$"
            "\n"
            rf"DarkProp: ${darkprop_high_flux / boltzmann_high_flux:.3f}"
            rf"\pm{darkprop_high_error / boltzmann_high_flux:.3f}$"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "0.8",
            "alpha": 0.9,
        },
    )

    ratio_axis.axhline(1.0, color="black", linewidth=1.0)
    valid_damascus_ratio = (
        np.isfinite(damascus_ratio)
        & np.isfinite(damascus_ratio_error)
        & (ratio_damascus > 0.0)
    )
    ratio_axis.errorbar(
        ratio_speed[valid_damascus_ratio],
        damascus_ratio[valid_damascus_ratio],
        yerr=np.vstack(
            (
                np.minimum(
                    damascus_ratio_error[valid_damascus_ratio],
                    0.8 * damascus_ratio[valid_damascus_ratio],
                ),
                damascus_ratio_error[valid_damascus_ratio],
            )
        ),
        fmt="s",
        markersize=3.2,
        markeredgewidth=0.0,
        linewidth=0.8,
        capsize=1.3,
        color="#d95f02",
    )
    darkprop_relative_error = np.full_like(darkprop_error, np.inf)
    np.divide(
        darkprop_error,
        darkprop_dphidv,
        out=darkprop_relative_error,
        where=darkprop_dphidv > 0.0,
    )
    valid_darkprop_ratio = (
        (darkprop_speed >= convergence_speed)
        & (darkprop_speed <= comparison_vmax)
        & np.isfinite(darkprop_ratio)
        & np.isfinite(darkprop_ratio_error)
        & (darkprop_dphidv > 0.0)
        & (darkprop_relative_error < 0.5)
    )
    darkprop_ratio_indices = np.flatnonzero(valid_darkprop_ratio)[::2]
    ratio_axis.errorbar(
        darkprop_speed[darkprop_ratio_indices],
        darkprop_ratio[darkprop_ratio_indices],
        yerr=darkprop_ratio_error[darkprop_ratio_indices],
        fmt="o",
        markersize=2.8,
        markeredgewidth=0.0,
        linewidth=0.7,
        capsize=1.1,
        color="#2673b8",
    )

    axis.set_xscale("log")
    axis.set_yscale("log")
    ratio_axis.set_xscale("log")
    axis.set_xlim(convergence_speed, comparison_vmax)
    axis.set_ylim(1.0e-8, 3.0e2)
    ratio_axis.set_ylim(0.0, 2.0)
    axis.set_ylabel(
        r"$d\Phi/dv$ [cm$^{-2}$ s$^{-1}$ (km/s)$^{-1}$]"
    )
    ratio_axis.set_xlabel(r"Speed $v$ [km/s]")
    ratio_axis.set_ylabel("MC / Boltz.")
    axis.grid(True, which="both", alpha=0.22)
    ratio_axis.grid(True, which="both", alpha=0.22)
    axis.legend(loc="lower left", fontsize=8.8)
    fig.suptitle(
        (
            r"$m_\chi=5\,\mathrm{GeV}$, "
            r"$\sigma_{\chi N}=5\times10^{-32}\,\mathrm{cm^2}$, "
            r"depth $=2.4\,\mathrm{km}$"
        ),
        fontsize=14,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    png_output = args.output.with_suffix(".png")
    fig.savefig(png_output, dpi=220)
    plt.close(fig)

    print(f"saved {args.output}")
    print(f"saved {png_output}")
    print(f"DaMaSCUS samples = {speed_beta.size}")
    print(f"DaMaSCUS rho = {rho:.8g} +/- {rho_error:.8g} GeV/cm^3")
    print(f"Boltzmann completed orders = {completed_orders}")
    print(f"Boltzmann comparison threshold = {convergence_speed:.8g} km/s")
    print(
        "high-speed fluxes [cm^-2 s^-1]: "
        f"DaMaSCUS={damascus_high_flux:.8g}+/-{damascus_high_error:.8g}, "
        f"Boltzmann={boltzmann_high_flux:.8g}, "
        f"DarkProp={darkprop_high_flux:.8g}+/-{darkprop_high_error:.8g}"
    )
    print(
        "high-speed integrated ratios: "
        f"DaMaSCUS/Boltzmann={damascus_high_flux / boltzmann_high_flux:.8g}, "
        f"DarkProp/Boltzmann={darkprop_high_flux / boltzmann_high_flux:.8g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
