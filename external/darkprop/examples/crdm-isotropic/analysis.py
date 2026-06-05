#!/usr/bin/env python3

""" Tools to reconstruct the underground DM flux. """
# pylint: disable=invalid-name,too-many-locals

import numpy as np
import h5py
from scipy import interpolate


def sigma_formatter(sig, unit="cm$^2$"):
    """ Scientific format """
    power = int(np.floor(np.log10(sig)))
    return ("$\\sigma_{\\chi p} = "
            f"{sig/10**power:.1f}\\times 10^{{{power:d}}}$ {unit}")


def loginterpolator(x, y, **kwargs):
    """Log-log interpolator.

    Parameters
    ----------
    x : 1-D ndarray
        x array.
    y : 1-D ndarray
        y array.

    Returns
    -------
    wrapper : function
        Interpolator function.

    """
    bounds_error = kwargs.pop("bounds_error", False)
    y = np.where(y <= 0, 1e-300, y)
    fill_value = kwargs.pop("fill_value", -np.inf)
    loglog = interpolate.interp1d(np.log(x), np.log(y),
                                  bounds_error=bounds_error,
                                  fill_value=fill_value, **kwargs)

    def wrapper_loglog(x, inter=loglog):
        x = np.array(x, ndmin=1, dtype=float)
        valid = x > 0.0
        result = np.zeros_like(x, dtype=float)
        result[valid] = np.exp(inter(np.log(x[valid])))
        return result

    return wrapper_loglog


def bl_to_xyz(bl):
    """ Galactic coordinate to Cartitian coordinate."""
    xyz = np.empty((bl.shape[0], 3), dtype=bl.dtype)
    b = bl[:, 0]
    l = bl[:, 1]  # noqa: E741
    xyz[:, 0] = np.cos(b) * np.cos(l)
    xyz[:, 1] = np.cos(b) * np.sin(l)
    xyz[:, 2] = np.sin(b)
    return xyz


def hist_estimate(x, *, weights=None, nbins=100):
    """ Histogram. """
    x = x.astype(np.float128)
    bins = np.geomspace(np.min(x), np.max(x), nbins + 1)
    centers = np.sqrt(bins[:-1] * bins[1:])
    # reverse to prevent the bug of np.histogram in the case of small weights
    bins = -bins[::-1]

    if weights is None:
        weights = np.full_like(x, 1)
    else:
        weights = weights.astype(np.float128)

    bars, bins = np.histogram(-x, bins=bins, density=True, weights=weights)
    return loginterpolator(centers, bars[::-1])


def spectrumz(filename, groupname, nbins=50):
    """ Underground DM flux. """
    with h5py.File(filename, 'r') as f:
        R0 = f['R0'][...]
        norm = f['norm'][...]

        group = f[groupname]
        Rd = group['R'][...]
        Nsim = group["Nsim"][...]

        dset_T = group["T"]
        T = np.empty(dset_T.shape, np.float64)
        dset_T.read_direct(T)

        dset_p = group["p"]
        p = np.empty(dset_p.shape, np.float64)
        dset_p.read_direct(p)

        dset_r = group["r"]
        r = np.empty(dset_r.shape, np.float64)
        dset_r.read_direct(r)

        dset_w = group["weight"]
        w = np.empty(dset_w.shape, np.float64)
        dset_w.read_direct(w)

        valid = ~np.isnan(T)
        T = T[valid]
        p = p[valid, :]
        r = r[valid, :]
        w = w[valid]

        ep = bl_to_xyz(p)
        er = bl_to_xyz(r)
        abscos = np.abs(np.sum(ep * er, axis=1))
        weights = w / abscos
        total_num = np.sum(weights)
        if T.size < 10:
            def wrapper0(E):
                return np.zeros_like(E)
            return wrapper0

        pdf = hist_estimate(T, weights=weights, nbins=nbins)

    def wrapper(E):
        return pdf(E) * total_num * R0**2 * norm / (4 * Rd**2 * Nsim)

    return wrapper


if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt

    basedir = os.path.dirname(__file__)
    mchi = np.loadtxt(os.path.join(basedir, "mchi.txt"))
    sigmas = np.loadtxt(os.path.join(basedir, "sigma.txt"))
    surface_fluxes = np.loadtxt(os.path.join(basedir, "surface_flux.txt"))
    Nbins = 100

    fluxfs = []
    for sigma in sigmas:
        datafile = os.path.join(basedir, (f"out/crdm/crdm_mchi{mchi:.3e}GeV_"
                                          f"sigma{sigma:.3e}cm2.hdf5"))
        # reconstruct the underground flux
        fluxf = spectrumz(datafile, "depth_0", nbins=Nbins)
        fluxfs.append(fluxf)

    Ekin = surface_fluxes[:, 0]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for i, sigma in enumerate(sigmas):
        axes[i].plot(Ekin, Ekin * surface_fluxes[:, i+1], label="surface")
        axes[i].plot(Ekin, Ekin * fluxfs[i](Ekin), label="underground")

        axes[i].set_xscale('log')
        axes[i].set_yscale('log')
        axes[i].set_xlim(1e-4, 2e1)
        if i == 0:
            axes[i].set_ylim(1e-10, 1e-5)
        if i == 1:
            axes[i].set_ylim(1e-8, 1e-3)
        axes[i].set_xlabel(r'$T_\chi$ [GeV]')
        axes[i].set_ylabel(r'$T_\chi$ $\mathrm{d}\bar\Phi_{\chi}/\mathrm{d}T$ '
                           r'[cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        axes[i].legend(loc="lower left")
        text = f"{sigma_formatter(sigma)}\n$m_\\chi = {mchi*1000}$ MeV"
        axes[i].text(0.6, 0.85, text, transform=axes[i].transAxes)
    fig.savefig(os.path.join(basedir, "crdm-flux.png"),
                dpi=200, bbox_inches='tight')
