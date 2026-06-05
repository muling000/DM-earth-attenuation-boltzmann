from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import iteration_dirac_dm_isoscalar_vector_coupling as wl


def bl_to_xyz(bl: np.ndarray) -> np.ndarray:
    xyz = np.empty((bl.shape[0], 3), dtype=np.float64)
    b = bl[:, 0]
    ell = bl[:, 1]
    xyz[:, 0] = np.cos(b) * np.cos(ell)
    xyz[:, 1] = np.cos(b) * np.sin(ell)
    xyz[:, 2] = np.sin(b)
    return xyz


def log_interp(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    x_new = np.asarray(x_new, dtype=np.float64)
    mask = (x > 0.0) & (y > 0.0)
    out = np.zeros_like(x_new)
    valid = (x_new >= np.min(x[mask])) & (x_new <= np.max(x[mask]))
    out[valid] = np.exp(
        np.interp(np.log(x_new[valid]), np.log(x[mask]), np.log(y[mask]))
    )
    return out


def reconstruct_darkprop_spectrum(
    hdf5_path: Path,
    *,
    group_name: str,
    bins: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    with h5py.File(hdf5_path, "r") as h5:
        r0 = float(np.asarray(h5["R0"]))
        norm = float(np.asarray(h5["norm"]))
        group = h5[group_name]
        rd = float(np.asarray(group["R"]))
        nsim = float(np.asarray(group["Nsim"]))
        nsample = float(np.asarray(group["Nsample"]))
        t = np.asarray(group["T"], dtype=np.float64)
        p = np.asarray(group["p"], dtype=np.float64)
        r = np.asarray(group["r"], dtype=np.float64)
        weight = np.asarray(group["weight"], dtype=np.float64)

    finite = np.isfinite(t) & np.isfinite(weight) & (t > 0.0)
    t = t[finite]
    p = p[finite]
    r = r[finite]
    weight = weight[finite]

    ep = bl_to_xyz(p)
    er = bl_to_xyz(r)
    abs_cos = np.abs(np.sum(ep * er, axis=1))
    valid = abs_cos > 1.0e-12
    t = t[valid]
    event_weight = weight[valid] / abs_cos[valid]

    hist_sum, _ = np.histogram(t, bins=bins, weights=event_weight)
    hist_sum_sq, _ = np.histogram(t, bins=bins, weights=event_weight * event_weight)
    widths = np.diff(bins)
    centers = np.sqrt(bins[:-1] * bins[1:])
    scale = r0 * r0 * norm / (4.0 * rd * rd * nsim)
    flux = hist_sum / widths * scale
    err = np.sqrt(hist_sum_sq) / widths * scale
    meta = {
        "R0_km": r0,
        "Rd_km": rd,
        "norm": norm,
        "Nsample": nsample,
        "Nsim": nsim,
        "valid_events": float(t.size),
    }
    return centers, flux, err, meta


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare a DarkProp WL-flux run with the original iteration result."
    )
    parser.add_argument(
        "--darkprop",
        type=Path,
        default=Path(__file__).resolve().parent
        / "darkprop-v0.3.0"
        / "examples"
        / "wl-flux"
        / "out"
        / "wl-flux-si-compare"
        / "wl-flux-si-compare_mchi5.000e+00GeV_sigma1.000e-33cm2.hdf5",
    )
    parser.add_argument(
        "--iteration",
        type=Path,
        default=Path(__file__).resolve().parent
        / "output"
        / "dirac_dm_vector_mchiMeV_5000_sigmaChiN_1.000e-33_wl-like_vmin1e-20_production.npz",
    )
    parser.add_argument("--group", default="depth_0")
    parser.add_argument("--nbins", type=int, default=80)
    parser.add_argument("--tmin-gev", type=float, default=1.8482279121277079e-10)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    iteration = np.load(args.iteration)
    t_iter_gev = np.asarray(iteration["T_grid"], dtype=np.float64) / wl.GEV
    cumulative = np.asarray(iteration["cumulative_spectra"], dtype=np.float64)[-1]
    iter_flux = cumulative / (2.0 * np.pi) ** 3 * wl.GEV

    tmax = float(t_iter_gev[-1])
    bins = np.geomspace(args.tmin_gev, tmax, args.nbins + 1)
    centers, dp_flux, dp_err, meta = reconstruct_darkprop_spectrum(
        args.darkprop,
        group_name=args.group,
        bins=bins,
    )
    iter_on_bins = log_interp(t_iter_gev, iter_flux, centers)

    ratio = np.full_like(centers, np.nan)
    ratio_err = np.full_like(centers, np.nan)
    positive = (iter_on_bins > 0.0) & (dp_flux > 0.0)
    ratio[positive] = dp_flux[positive] / iter_on_bins[positive]
    ratio_err[positive] = dp_err[positive] / iter_on_bins[positive]

    good = positive & (dp_err / np.maximum(dp_flux, 1.0e-300) < 0.35)
    if np.any(good):
        median_ratio = float(np.median(ratio[good]))
        p16, p84 = np.percentile(ratio[good], [16, 84])
    else:
        median_ratio = float("nan")
        p16 = float("nan")
        p84 = float("nan")

    output = args.output
    if output is None:
        output = args.darkprop.with_name("wl_flux_darkprop_vs_iteration.png")
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(8.0, 6.4),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )
    plot_mask = t_iter_gev >= args.tmin_gev
    ax.plot(
        t_iter_gev[plot_mask],
        t_iter_gev[plot_mask] * iter_flux[plot_mask],
        color="black",
        linewidth=1.8,
        label="Iteration",
    )
    ax.errorbar(
        centers,
        centers * dp_flux,
        yerr=centers * dp_err,
        fmt="o",
        markersize=3.0,
        linewidth=0.8,
        capsize=1.5,
        color="#2f7ed8",
        label="DarkProp MC",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(args.tmin_gev, tmax)
    ax.set_ylabel(r"$T\,d\bar\Phi/(dT\,d\Omega)$ [cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best")
    ax.set_title(r"$m_\chi=5$ GeV, $\sigma_{\chi N}=10^{-33}$ cm$^2$, depth $=2.4$ km")

    rax.axhline(1.0, color="black", linewidth=1.0)
    rax.errorbar(
        centers[positive],
        ratio[positive],
        yerr=ratio_err[positive],
        fmt="o",
        markersize=3.0,
        linewidth=0.8,
        capsize=1.5,
        color="#2f7ed8",
    )
    rax.set_xscale("log")
    rax.set_ylim(0.0, 2.0)
    rax.set_xlabel(r"Kinetic energy $T$ [GeV]")
    rax.set_ylabel("MC / iter.")
    rax.grid(True, which="both", alpha=0.25)

    fig.savefig(output, dpi=200)
    plt.close(fig)

    print(f"wrote {output}")
    print(f"Nsample = {meta['Nsample']:.0f}")
    print(f"Nsim = {meta['Nsim']:.0f}")
    print(f"valid events = {meta['valid_events']:.0f}")
    print(f"norm = {meta['norm']:.8e} cm^-2 s^-1 sr^-1")
    print(f"median ratio = {median_ratio:.6g}")
    print(f"central 68% ratio interval = [{p16:.6g}, {p84:.6g}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
