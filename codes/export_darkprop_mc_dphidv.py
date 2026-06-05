from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def bl_to_xyz(bl: np.ndarray) -> np.ndarray:
    xyz = np.empty((bl.shape[0], 3), dtype=np.float64)
    b = bl[:, 0]
    ell = bl[:, 1]
    xyz[:, 0] = np.cos(b) * np.cos(ell)
    xyz[:, 1] = np.cos(b) * np.sin(ell)
    xyz[:, 2] = np.sin(b)
    return xyz


def reconstruct_dphidv(
    hdf5_path: Path,
    *,
    group_name: str,
    bins: np.ndarray,
    min_abs_cos: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    with h5py.File(hdf5_path, "r") as h5:
        mchi = float(np.asarray(h5["mchi"]))
        r0 = float(np.asarray(h5["R0"]))
        norm = float(np.asarray(h5["norm"]))
        group = h5[group_name]
        rd = float(np.asarray(group["R"]))
        nsim = float(np.asarray(group["Nsim"]))
        nsample = float(np.asarray(group["Nsample"]))
        kinetic_energy = np.asarray(group["T"], dtype=np.float64)
        momentum_direction = np.asarray(group["p"], dtype=np.float64)
        position_direction = np.asarray(group["r"], dtype=np.float64)
        weight = np.asarray(group["weight"], dtype=np.float64)

    finite = np.isfinite(kinetic_energy) & np.isfinite(weight) & (kinetic_energy > 0.0)
    kinetic_energy = kinetic_energy[finite]
    momentum_direction = momentum_direction[finite]
    position_direction = position_direction[finite]
    weight = weight[finite]

    ep = bl_to_xyz(momentum_direction)
    er = bl_to_xyz(position_direction)
    abs_cos = np.abs(np.sum(ep * er, axis=1))
    valid = abs_cos > min_abs_cos
    kinetic_energy = kinetic_energy[valid]
    event_weight = weight[valid] / abs_cos[valid]

    # Match the non-relativistic v grid convention used in the Boltzmann iteration.
    velocity = np.sqrt(2.0 * kinetic_energy / mchi)
    centers = 0.5 * (bins[:-1] + bins[1:])
    widths = np.diff(bins)

    hist_sum, _ = np.histogram(velocity, bins=bins, weights=event_weight)
    hist_sum_sq, _ = np.histogram(velocity, bins=bins, weights=event_weight * event_weight)

    meta = {
        "mchi_GeV": mchi,
        "R0_km": r0,
        "Rd_km": rd,
        "norm_cm-2_s-1_sr-1": norm,
        "Nsample": nsample,
        "Nsim": nsim,
        "valid_events": float(velocity.size),
        "hist_sum": hist_sum,
        "hist_sum_sq": hist_sum_sq,
    }
    return centers, hist_sum, hist_sum_sq, meta


def reconstruct_merged_dphidv(
    hdf5_paths: list[Path],
    *,
    group_name: str,
    nbins: int,
    vmin: float,
    vmax: float,
    min_abs_cos: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    bins = np.linspace(vmin, vmax, nbins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    widths = np.diff(bins)

    total_hist = np.zeros(nbins, dtype=np.float64)
    total_hist_sq = np.zeros(nbins, dtype=np.float64)
    total_nsim = 0.0
    total_nsample = 0.0
    total_valid = 0.0
    reference: dict[str, float] | None = None

    for path in hdf5_paths:
        _, hist, hist_sq, meta = reconstruct_dphidv(
            path,
            group_name=group_name,
            bins=bins,
            min_abs_cos=min_abs_cos,
        )
        if reference is None:
            reference = meta
        else:
            for key in ("mchi_GeV", "R0_km", "Rd_km", "norm_cm-2_s-1_sr-1"):
                if not np.isclose(reference[key], meta[key], rtol=1.0e-6, atol=0.0):
                    raise ValueError(f"{path} has incompatible {key}: {meta[key]} != {reference[key]}")
        total_hist += hist
        total_hist_sq += hist_sq
        total_nsim += meta["Nsim"]
        total_nsample += meta["Nsample"]
        total_valid += meta["valid_events"]

    if reference is None:
        raise ValueError("No HDF5 files were provided.")

    # DarkProp's official reconstruction gives the angular-averaged flux per sr.
    # Multiply by 4*pi to match the iteration script's dPhi/dv [cm^-2 s^-1].
    scale = (
        4.0
        * np.pi
        * reference["R0_km"]
        * reference["R0_km"]
        * reference["norm_cm-2_s-1_sr-1"]
        / (4.0 * reference["Rd_km"] * reference["Rd_km"] * total_nsim)
    )
    dphidv = total_hist / widths * scale
    err = np.sqrt(total_hist_sq) / widths * scale
    meta = {
        "mchi_GeV": reference["mchi_GeV"],
        "R0_km": reference["R0_km"],
        "Rd_km": reference["Rd_km"],
        "norm_cm-2_s-1_sr-1": reference["norm_cm-2_s-1_sr-1"],
        "Nsample": total_nsample,
        "Nsim": total_nsim,
        "valid_events": total_valid,
        "vmin": vmin,
        "vmax": vmax,
        "min_abs_cos": min_abs_cos,
        "n_files": float(len(hdf5_paths)),
    }
    return centers, dphidv, err, meta


def main() -> int:
    parser = argparse.ArgumentParser(description="Export DarkProp MC dPhi/dv as a table.")
    parser.add_argument(
        "--darkprop",
        type=Path,
        nargs="+",
        default=Path(__file__).resolve().parent.parent
        / "external"
        / "darkprop"
        / "examples"
        / "surface-flux"
        / "out"
        / "surface-flux-si-compare"
        / "surface-flux-si-compare_mchi5.000e+00GeV_sigma1.000e-33cm2.hdf5",
    )
    parser.add_argument("--group", default="depth_0")
    parser.add_argument("--nbins", type=int, default=180)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=0.0027)
    parser.add_argument(
        "--min-abs-cos",
        type=float,
        default=1.0e-3,
        help="Drop nearly tangent detector-sphere crossings to reduce the 1/|cos| estimator variance.",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    darkprop_paths = args.darkprop if isinstance(args.darkprop, list) else [args.darkprop]
    output = args.output
    if output is None:
        output = darkprop_paths[0].with_name("darkprop_mc_dphidv.dat")
    output.parent.mkdir(parents=True, exist_ok=True)

    velocity, dphidv, err, meta = reconstruct_merged_dphidv(
        darkprop_paths,
        group_name=args.group,
        nbins=args.nbins,
        vmin=args.vmin,
        vmax=args.vmax,
        min_abs_cos=args.min_abs_cos,
    )
    table = np.column_stack([velocity, dphidv, err])
    header = (
        "v_over_c dPhi_dv_cm^-2_s^-1 stat_error_cm^-2_s^-1\n"
        + " ".join(f"{key}={value:.16e}" for key, value in meta.items())
    )
    np.savetxt(output, table, fmt="%.16e", header=header)

    nonzero = dphidv > 0.0
    print(f"wrote {output}")
    print(f"rows = {table.shape[0]}")
    print(f"Nsample = {meta['Nsample']:.0f}")
    print(f"Nsim = {meta['Nsim']:.0f}")
    if np.any(nonzero):
        print(f"dPhi/dv range = [{np.min(dphidv[nonzero]):.8e}, {np.max(dphidv):.8e}] cm^-2 s^-1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
