#!/usr/bin/env python3

""" Reconstruct the SolarDM flux on the Earth. """
# pylint: disable=invalid-name

import numpy as np
import h5py
from scipy import interpolate


def sigma_formatter(sig, unit="cm$^2$"):
    """ Scientific format """
    power = int(np.floor(np.log10(sig)))
    return ("$\\sigma_{\\chi p} = "
            f"{sig/10**power:.1f}\\times 10^{{{power:d}}}$ {unit}")


def hist_estimate(x, *, weights=None, nbins=100):
    """ Histogram. """

    bins = np.linspace(np.min(x), np.max(x), nbins + 1)
    xitp = (bins[:-1] + bins[1:]) / 2
    xitp[0] = bins[0]
    xitp[-1] = bins[-1]
    if weights is None:
        weights = np.full_like(x, 1)
    bars, bins = np.histogram(x, bins=bins, density=True, weights=weights)
    return interpolate.interp1d(xitp, bars, bounds_error=False, fill_value=0)


def spectrum(filename, groupname, nbins=50):
    """Construct spectrum """
    with h5py.File(filename, 'r') as f:
        group = f[groupname]

        dset_T = group["T"]
        T = np.empty(dset_T.shape, np.float64)
        dset_T.read_direct(T)

        dset_w = group["weight"]
        w = np.empty(dset_w.shape, np.float64)
        dset_w.read_direct(w)

        valid = ~np.isnan(T)
        T = T[valid]
        w = w[valid]

        if T.size < 10:
            def wrapper0(E):
                return np.zeros_like(E)
            return wrapper0

        pdf = hist_estimate(T, weights=w, nbins=nbins)

    def wrapper(E):
        au = 14959787070000  # cm
        production_rate = 4.6e26  # sec^-1
        return 1 / (4 * np.pi * au**2) * production_rate * pdf(E)

    return wrapper


if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt

    basedir = os.path.dirname(__file__)
    sigmas = [1e-50, 1e-34]
    mchi = 1.00075e-3
    Nbins = 500

    prefix = "solardm"

    fluxfs = []
    for sigma in sigmas:
        datafile = os.path.join(
            basedir, "out", prefix,
            f"{prefix}_mchi{mchi:.3e}GeV_sigma{sigma:.3e}cm2.hdf5")
        fluxf = spectrum(datafile, "depth_0", nbins=Nbins)
        fluxfs.append(fluxf)

    Ekin = np.linspace(0, 0.004, 1000)

    fig, ax = plt.subplots()
    for i, sigma in enumerate(sigmas):
        label = "w/o attenuation" if i == 0 else "w/ attenuation"
        ax.plot((Ekin + mchi) * 1000, fluxfs[i](Ekin) / 1000, label=label)

    ax.set_yscale('log')
    ax.set_xlim(1, 4.5)
    ax.set_xlabel(r'$E_\chi$ [MeV]')
    ax.set_ylabel(r'$\mathrm{d}\Phi_{\oplus}/\mathrm{d}E_\chi$ '
                  r'[cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
    ax.set_ylim(1e-7, 1e1)
    ax.legend(loc="lower left")
    text = f"{sigma_formatter(sigma)}\n$m_\\chi = {mchi*1000:.1f}$ MeV"
    ax.text(0.6, 0.85, text, transform=ax.transAxes)
    fig.savefig(os.path.join(basedir, "solardm-flux.png"),
                dpi=200, bbox_inches='tight')
