from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
import numpy.typing as npt

try:
    from numba import njit, prange

    _NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency fallback
    _NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore[no-redef]
        def decorator(func):
            return func

        return decorator

    def prange(*args):  # type: ignore[no-redef]
        return range(*args)


FloatArray = npt.NDArray[np.float64]
SpectrumInput = Callable[[FloatArray], FloatArray] | Sequence[float] | FloatArray

# Unit and input conventions extracted from the Mathematica derivation script.
EV = 1.0
GEV = 1.0e9 * EV
MEV = 1.0e6 * EV
KEV = 1.0e3 * EV
CM = 5.0677307165483375e13 / GEV
METER = 100.0 * CM
SEC = 1.5192674479961274e24 / GEV
KG = 5.609588606807705e26 * GEV
GRAM = 1.0e-3 * KG
KM = 1000.0 * METER
# Match DarkProp's constants::mu convention used by SIDM::sigmaEff.
AMU = 0.93149410242 * GEV
NUCLEON_MASS = AMU

BENCHMARK_EARTH_DENSITY = 2.7 * GRAM / CM**3
BENCHMARK_EARTH_RADIUS = 6371.0 * KM
BENCHMARK_DM_DENSITY = 0.3 * GEV / CM**3
BENCHMARK_VRMS = 270.0 / 300000.0
BENCHMARK_VSUN = 220.0 / 300000.0
BENCHMARK_VESC = 544.0 / 300000.0
BENCHMARK_VEARTH = 240.0 / 300000.0
BENCHMARK_VMAX = BENCHMARK_VESC + BENCHMARK_VEARTH
BENCHMARK_CRUST_COMPOSITION = (
    ("O", 16, 0.466),
    ("Si", 28, 0.277),
    ("Al", 27, 0.081),
    ("Fe", 56, 0.050),
    ("Ca", 40, 0.036),
    ("K", 39, 0.028),
    ("Na", 23, 0.026),
    ("Mg", 24, 0.021),
)
JINPING_DEPTH = 2.4 * KM


@dataclass(frozen=True)
class NuclearSpecies:
    name: str
    mass_number: int
    mass: float
    number_density: float

    def __post_init__(self) -> None:
        if self.mass_number <= 0:
            raise ValueError("mass_number must be positive.")
        if self.mass <= 0.0:
            raise ValueError("mass must be positive.")
        if self.number_density <= 0.0:
            raise ValueError("number_density must be positive.")


@dataclass(frozen=True)
class BenchmarkModel:
    m_chi: float
    cutoff_scale: float
    earth_radius: float
    species: tuple[NuclearSpecies, ...]

    def __post_init__(self) -> None:
        if self.m_chi <= 0.0:
            raise ValueError("m_chi must be positive.")
        if self.cutoff_scale <= 0.0:
            raise ValueError("cutoff_scale must be positive.")
        if self.earth_radius <= 0.0:
            raise ValueError("earth_radius must be positive.")
        if not self.species:
            raise ValueError("At least one nuclear species is required.")


@dataclass(frozen=True)
class PhaseSpaceGrid:
    r: FloatArray
    u: FloatArray
    T: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "r", _as_1d_array(self.r, "r"))
        object.__setattr__(self, "u", _as_1d_array(self.u, "u"))
        object.__setattr__(self, "T", _as_1d_array(self.T, "T"))
        if np.any(self.r < 0.0):
            raise ValueError("r grid must be non-negative.")
        if np.any(np.abs(self.u) > 1.0):
            raise ValueError("u grid must lie in [-1, 1].")
        if np.any(self.T <= 0.0):
            raise ValueError("T grid must be strictly positive.")

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.r.size, self.u.size, self.T.size)

    @property
    def size(self) -> int:
        nr, nu, nt = self.shape
        return nr * nu * nt


@dataclass(frozen=True)
class QuadratureSpec:
    n_x: int
    n_y: int
    n_z: int

    def __post_init__(self) -> None:
        if self.n_x <= 0 or self.n_y <= 0 or self.n_z <= 0:
            raise ValueError("Quadrature orders must be positive.")


@dataclass
class _AxisStencil:
    lo: npt.NDArray[np.int64]
    hi: npt.NDArray[np.int64]
    w_lo: FloatArray
    w_hi: FloatArray


@dataclass
class _SpeciesCache:
    species: NuclearSpecies
    kappa: float
    y_nodes: FloatArray
    omega_y: FloatArray
    t_prime: FloatArray
    t_stencil: _AxisStencil
    c_angle: FloatArray
    kernel: FloatArray


class BenchmarkIterationSolver:
    """
    Matrix-free implementation of Eq. `benchmark_iteration_matrix_final`
    for the benchmark vector operator in the non-relativistic benchmark
    subsection. The code keeps the discrete structure of the note:

        F_{i+1}(r_a, u_b, T_c)
          = sum_A n_A sum_k sum_m sum_n
            omega_x omega_y omega_z K_tilde
            M_r M_u M_T F_i .

    Off-grid samples are evaluated with trilinear interpolation on
    (r, u, T). The solver does not auto-run a benchmark when imported.
    """

    def __init__(
        self,
        model: BenchmarkModel,
        grid: PhaseSpaceGrid,
        quadrature: QuadratureSpec,
        *,
        strict: bool = True,
        tau_segment_edges: Sequence[float] | None = None,
    ) -> None:
        self.model = model
        self.grid = grid
        self.quadrature = quadrature
        self.strict = strict
        self.tau_segment_edges = _tau_segment_edges(tau_segment_edges)

        if grid.r[-1] < model.earth_radius:
            raise ValueError("The r grid must reach the Earth's radius.")

        self._x_nodes, self._x_weights = np.polynomial.legendre.leggauss(quadrature.n_x)
        self._y_nodes, self._y_weights = np.polynomial.legendre.leggauss(quadrature.n_y)
        self._z_nodes = np.cos(
            (2.0 * np.arange(1, quadrature.n_z + 1) - 1.0) * math.pi / (2.0 * quadrature.n_z)
        ).astype(np.float64)
        self._omega_z = math.pi / quadrature.n_z

        self.mean_free_path = self._compute_mean_free_path()
        self._build_geometry_cache()
        self._species_cache = tuple(self._build_species_cache(species) for species in model.species)
        self._build_numba_cache()

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.grid.shape

    def initial_term(self, surface_spectrum: SpectrumInput) -> FloatArray:
        spectrum = _evaluate_spectrum(surface_spectrum, self.grid.T)
        attenuation = np.exp(-(self._x_grid + self._x_oplus) / self.mean_free_path)
        return attenuation[:, :, None] * spectrum[None, None, :]

    def apply_iteration(self, field: FloatArray) -> FloatArray:
        field = np.ascontiguousarray(self._validate_field(field))
        return _apply_iteration_kernel(
            field,
            self._omega_x_numba,
            self._r_lo_numba,
            self._r_hi_numba,
            self._rw_lo_numba,
            self._rw_hi_numba,
            self._species_base_weight_numba,
            self._species_t_lo_numba,
            self._species_t_hi_numba,
            self._species_tw_lo_numba,
            self._species_tw_hi_numba,
            self._species_u_lo_numba,
            self._species_u_hi_numba,
            self._species_uw_lo_numba,
            self._species_uw_hi_numba,
        )

    def _build_geometry_cache(self) -> None:
        r_vals = self.grid.r[:, None]
        u_vals = self.grid.u[None, :]

        self._x_grid = r_vals * u_vals
        self._y_grid = r_vals * np.sqrt(np.maximum(0.0, 1.0 - u_vals**2))
        under_root = np.maximum(0.0, self.model.earth_radius**2 - self._y_grid**2)
        self._x_oplus = np.sqrt(under_root)

        eta = np.maximum(0.0, (self._x_grid + self._x_oplus) / self.mean_free_path)
        self._build_tau_segmented_path_quadrature(eta)
        x_shifted = self._x_grid[:, :, None] - self._x_pp
        y_broadcast = self._y_grid[:, :, None]
        self._r_star = np.sqrt(np.maximum(0.0, x_shifted * x_shifted + y_broadcast * y_broadcast))
        self._check_grid_coverage("r*", self._r_star, self.grid.r)
        self._r_stencil = _linear_stencil(self.grid.r, self._r_star)

    def _build_tau_segmented_path_quadrature(self, eta: FloatArray) -> None:
        tau_nodes: list[FloatArray] = []
        omega_nodes: list[FloatArray] = []
        for tau_lo, tau_hi in zip(self.tau_segment_edges[:-1], self.tau_segment_edges[1:]):
            active = eta > tau_lo
            upper = np.minimum(eta, tau_hi)
            width = np.where(active, upper - tau_lo, 0.0)
            center = tau_lo + 0.5 * width
            tau_k = center[:, :, None] + 0.5 * width[:, :, None] * self._x_nodes[None, None, :]
            omega_k = (
                self.mean_free_path
                * 0.5
                * width[:, :, None]
                * self._x_weights[None, None, :]
                * np.exp(-tau_k)
            )
            tau_nodes.append(np.where(active[:, :, None], tau_k, 0.0))
            omega_nodes.append(np.where(active[:, :, None], omega_k, 0.0))

        self._omega_x = np.concatenate(omega_nodes, axis=2)
        self._x_pp = self.mean_free_path * np.concatenate(tau_nodes, axis=2)

    def _build_species_cache(self, species: NuclearSpecies) -> _SpeciesCache:
        kappa = 4.0 * self.model.m_chi * species.mass / (self.model.m_chi + species.mass) ** 2
        y_nodes = 0.5 * kappa * (1.0 + self._y_nodes)
        t_prime = self.grid.T[:, None] / (1.0 - y_nodes[None, :])
        omega_y = (
            0.5
            * kappa
            * self._y_weights[None, :]
            * self.grid.T[:, None]
            / (1.0 - y_nodes[None, :]) ** 2
        )
        self._check_grid_coverage(f"T' ({species.name})", t_prime, self.grid.T)

        t_stencil = _linear_stencil(self.grid.T, t_prime)
        c_angle = self._c_angle(species, self.grid.T[:, None], t_prime)
        kernel = self._tilde_kernel(species, self.grid.T[:, None], t_prime)

        return _SpeciesCache(
            species=species,
            kappa=kappa,
            y_nodes=y_nodes.astype(np.float64),
            omega_y=omega_y.astype(np.float64),
            t_prime=t_prime.astype(np.float64),
            t_stencil=t_stencil,
            c_angle=c_angle.astype(np.float64),
            kernel=kernel.astype(np.float64),
        )

    def _build_numba_cache(self) -> None:
        self._u_grid_numba = np.ascontiguousarray(self.grid.u, dtype=np.float64)
        self._sqrt_one_minus_u_sq_numba = np.ascontiguousarray(
            np.sqrt(np.maximum(0.0, 1.0 - self.grid.u**2)),
            dtype=np.float64,
        )
        self._z_nodes_numba = np.ascontiguousarray(self._z_nodes, dtype=np.float64)
        self._omega_x_numba = np.ascontiguousarray(self._omega_x, dtype=np.float64)
        self._r_lo_numba = np.ascontiguousarray(self._r_stencil.lo, dtype=np.int64)
        self._r_hi_numba = np.ascontiguousarray(self._r_stencil.hi, dtype=np.int64)
        self._rw_lo_numba = np.ascontiguousarray(self._r_stencil.w_lo, dtype=np.float64)
        self._rw_hi_numba = np.ascontiguousarray(self._r_stencil.w_hi, dtype=np.float64)

        self._species_number_density_numba = np.ascontiguousarray(
            np.array([cache.species.number_density for cache in self._species_cache], dtype=np.float64)
        )
        self._species_omega_y_numba = np.ascontiguousarray(
            np.stack([cache.omega_y for cache in self._species_cache], axis=0),
            dtype=np.float64,
        )
        self._species_kernel_numba = np.ascontiguousarray(
            np.stack([cache.kernel for cache in self._species_cache], axis=0),
            dtype=np.float64,
        )
        self._species_c_angle_numba = np.ascontiguousarray(
            np.stack([cache.c_angle for cache in self._species_cache], axis=0),
            dtype=np.float64,
        )
        self._species_t_lo_numba = np.ascontiguousarray(
            np.stack([cache.t_stencil.lo for cache in self._species_cache], axis=0),
            dtype=np.int64,
        )
        self._species_t_hi_numba = np.ascontiguousarray(
            np.stack([cache.t_stencil.hi for cache in self._species_cache], axis=0),
            dtype=np.int64,
        )
        self._species_tw_lo_numba = np.ascontiguousarray(
            np.stack([cache.t_stencil.w_lo for cache in self._species_cache], axis=0),
            dtype=np.float64,
        )
        self._species_tw_hi_numba = np.ascontiguousarray(
            np.stack([cache.t_stencil.w_hi for cache in self._species_cache], axis=0),
            dtype=np.float64,
        )
        self._species_base_weight_numba = np.ascontiguousarray(
            self._species_number_density_numba[:, None, None]
            * self._species_omega_y_numba
            * self._species_kernel_numba
            * self._omega_z,
            dtype=np.float64,
        )
        self._build_u_stencil_cache()

    def _build_u_stencil_cache(self) -> None:
        u_values = self._u_grid_numba[None, :, None, None, None]
        u_prefactor = self._sqrt_one_minus_u_sq_numba[None, :, None, None, None]
        c_angle = self._species_c_angle_numba[:, None, :, :, None]
        z_nodes = self._z_nodes_numba[None, None, None, None, :]
        radius = u_prefactor * np.sqrt(np.maximum(0.0, 1.0 - c_angle * c_angle))
        u_samples = u_values * c_angle + radius * z_nodes

        if self.strict:
            overflow = float(np.max(np.abs(u_samples) - 1.0))
            if overflow > 5.0e-10:
                raise ValueError(
                    "Sampled u' lies outside [-1, 1]. "
                    "This indicates a kinematic or grid inconsistency."
                )

        u_samples = np.clip(u_samples, self.grid.u[0], self.grid.u[-1])
        u_stencil = _linear_stencil(self.grid.u, u_samples)
        self._species_u_lo_numba = np.ascontiguousarray(u_stencil.lo, dtype=np.int64)
        self._species_u_hi_numba = np.ascontiguousarray(u_stencil.hi, dtype=np.int64)
        self._species_uw_lo_numba = np.ascontiguousarray(u_stencil.w_lo, dtype=np.float64)
        self._species_uw_hi_numba = np.ascontiguousarray(u_stencil.w_hi, dtype=np.float64)

    def _compute_mean_free_path(self) -> float:
        return mean_free_path_for_model(self.model)

    def _sigma_nr(self, species: NuclearSpecies) -> float:
        return sigma_nr_for_species(self.model, species)

    def _matrix_element_sq(
        self,
        species: NuclearSpecies,
        kinetic_out: FloatArray,
        kinetic_in: FloatArray,
    ) -> FloatArray:
        e_out = self.model.m_chi + kinetic_out
        e_in = self.model.m_chi + kinetic_in
        prefactor = 8.0 * species.mass_number**2 / self.model.cutoff_scale**4
        bracket = (
            species.mass**2 * (e_in * e_in + e_out * e_out)
            - species.mass * (e_in - e_out) * (self.model.m_chi**2 + species.mass**2)
        )
        return prefactor * bracket

    def _tilde_kernel(
        self,
        species: NuclearSpecies,
        kinetic_out: FloatArray,
        kinetic_in: FloatArray,
    ) -> FloatArray:
        p_in_sq = np.maximum(0.0, kinetic_in * (kinetic_in + 2.0 * self.model.m_chi))
        matrix_sq = self._matrix_element_sq(species, kinetic_out, kinetic_in)
        denominator = 32.0 * math.pi**2 * species.mass * np.maximum(p_in_sq, 1.0e-300)
        return matrix_sq / denominator

    def _c_angle(
        self,
        species: NuclearSpecies,
        kinetic_out: FloatArray,
        kinetic_in: FloatArray,
    ) -> FloatArray:
        e_out = self.model.m_chi + kinetic_out
        e_in = self.model.m_chi + kinetic_in
        p_out = np.sqrt(np.maximum(0.0, kinetic_out * (kinetic_out + 2.0 * self.model.m_chi)))
        p_in = np.sqrt(np.maximum(0.0, kinetic_in * (kinetic_in + 2.0 * self.model.m_chi)))
        numerator = e_out * e_in - self.model.m_chi**2 + species.mass * (e_out - e_in)
        denominator = np.maximum(p_out * p_in, 1.0e-300)
        c_angle = numerator / denominator
        if self.strict:
            overflow = np.max(np.abs(c_angle) - 1.0)
            if overflow > 5.0e-10:
                raise ValueError(
                    f"Kinematics produced |C_A| > 1 for species {species.name}. "
                    "Adjust the grid or disable strict mode to clip small excursions."
                )
        return np.clip(c_angle, -1.0, 1.0)

    def _check_grid_coverage(self, name: str, samples: FloatArray, grid_values: FloatArray) -> None:
        sample_min = float(np.min(samples))
        sample_max = float(np.max(samples))
        grid_min = float(grid_values[0])
        grid_max = float(grid_values[-1])
        tol = 1.0e-12 * max(1.0, abs(grid_min), abs(grid_max))

        if self.strict and (sample_min < grid_min - tol or sample_max > grid_max + tol):
            raise ValueError(
                f"{name} samples leave the grid range [{grid_min}, {grid_max}]. "
                "Enlarge the grid or disable strict mode to clamp."
            )

    def _validate_field(self, field: FloatArray) -> FloatArray:
        arr = np.asarray(field, dtype=np.float64)
        if arr.shape != self.shape:
            raise ValueError(f"field shape must be {self.shape}, got {arr.shape}.")
        return arr


@njit(cache=True, inline="always")
def _trilinear_sample(
    field: FloatArray,
    r_lo: int,
    r_hi: int,
    rw_lo: float,
    rw_hi: float,
    u_lo: int,
    u_hi: int,
    uw_lo: float,
    uw_hi: float,
    t_lo: int,
    t_hi: int,
    tw_lo: float,
    tw_hi: float,
) -> float:
    value_rlo_tlo = uw_lo * field[r_lo, u_lo, t_lo] + uw_hi * field[r_lo, u_hi, t_lo]
    value_rlo_thi = uw_lo * field[r_lo, u_lo, t_hi] + uw_hi * field[r_lo, u_hi, t_hi]
    value_rhi_tlo = uw_lo * field[r_hi, u_lo, t_lo] + uw_hi * field[r_hi, u_hi, t_lo]
    value_rhi_thi = uw_lo * field[r_hi, u_lo, t_hi] + uw_hi * field[r_hi, u_hi, t_hi]

    value_rlo = tw_lo * value_rlo_tlo + tw_hi * value_rlo_thi
    value_rhi = tw_lo * value_rhi_tlo + tw_hi * value_rhi_thi
    return rw_lo * value_rlo + rw_hi * value_rhi


@njit(cache=True, parallel=True)
def _apply_iteration_kernel(
    field: FloatArray,
    omega_x: FloatArray,
    r_lo_all: npt.NDArray[np.int64],
    r_hi_all: npt.NDArray[np.int64],
    rw_lo_all: FloatArray,
    rw_hi_all: FloatArray,
    species_base_weight: FloatArray,
    species_t_lo: npt.NDArray[np.int64],
    species_t_hi: npt.NDArray[np.int64],
    species_tw_lo: FloatArray,
    species_tw_hi: FloatArray,
    species_u_lo: npt.NDArray[np.int64],
    species_u_hi: npt.NDArray[np.int64],
    species_uw_lo: FloatArray,
    species_uw_hi: FloatArray,
) -> FloatArray:
    nr, nu, nt = field.shape
    ns = species_base_weight.shape[0]
    nx = omega_x.shape[2]
    ny = species_base_weight.shape[2]
    nz = species_u_lo.shape[4]

    result = np.zeros_like(field)
    total_size = nr * nu * nt

    for flat_index in prange(total_size):
        a = flat_index // (nu * nt)
        rem = flat_index - a * (nu * nt)
        b = rem // nt
        c = rem - b * nt

        accum = 0.0

        for s in range(ns):
            for k in range(nx):
                weight_x = omega_x[a, b, k]
                if weight_x == 0.0:
                    continue

                r_lo = r_lo_all[a, b, k]
                r_hi = r_hi_all[a, b, k]
                rw_lo = rw_lo_all[a, b, k]
                rw_hi = rw_hi_all[a, b, k]

                for m in range(ny):
                    base_weight = weight_x * species_base_weight[s, c, m]
                    if base_weight == 0.0:
                        continue

                    t_lo = species_t_lo[s, c, m]
                    t_hi = species_t_hi[s, c, m]
                    tw_lo = species_tw_lo[s, c, m]
                    tw_hi = species_tw_hi[s, c, m]

                    sample_sum = 0.0
                    for n in range(nz):
                        sample_sum += _trilinear_sample(
                            field,
                            r_lo,
                            r_hi,
                            rw_lo,
                            rw_hi,
                            species_u_lo[s, b, c, m, n],
                            species_u_hi[s, b, c, m, n],
                            species_uw_lo[s, b, c, m, n],
                            species_uw_hi[s, b, c, m, n],
                            t_lo,
                            t_hi,
                            tw_lo,
                            tw_hi,
                        )

                    accum += base_weight * sample_sum

        result[a, b, c] = accum

    return result


def _as_1d_array(values: Sequence[float] | FloatArray, name: str) -> FloatArray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D array.")
    if arr.size == 0:
        raise ValueError(f"{name} must be non-empty.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr


def _tau_segment_edges(values: Sequence[float] | None) -> FloatArray:
    if values is None:
        values = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 24.0)
    edges = np.asarray(values, dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("tau_segment_edges must contain at least two values.")
    if edges[0] != 0.0:
        raise ValueError("tau_segment_edges must start at 0.")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("tau_segment_edges must be strictly increasing.")
    return np.ascontiguousarray(edges, dtype=np.float64)


def _linear_stencil(grid_values: FloatArray, samples: FloatArray) -> _AxisStencil:
    samples = np.asarray(samples, dtype=np.float64)
    if grid_values.size == 1:
        zeros = np.zeros(samples.shape, dtype=np.int64)
        ones = np.ones(samples.shape, dtype=np.float64)
        return _AxisStencil(lo=zeros, hi=zeros, w_lo=ones, w_hi=np.zeros_like(ones))

    clipped = np.clip(samples, grid_values[0], grid_values[-1])
    hi = np.searchsorted(grid_values, clipped, side="right")
    hi = np.clip(hi, 1, grid_values.size - 1).astype(np.int64)
    lo = hi - 1

    g_lo = grid_values[lo]
    g_hi = grid_values[hi]
    denom = g_hi - g_lo
    frac = np.where(denom > 0.0, (clipped - g_lo) / denom, 0.0)
    w_hi = frac.astype(np.float64)
    w_lo = (1.0 - frac).astype(np.float64)
    return _AxisStencil(lo=lo, hi=hi, w_lo=w_lo, w_hi=w_hi)


def _evaluate_spectrum(surface_spectrum: SpectrumInput, kinetic_grid: FloatArray) -> FloatArray:
    if callable(surface_spectrum):
        try:
            trial = np.asarray(surface_spectrum(kinetic_grid), dtype=np.float64)
        except Exception:
            trial = np.asarray([], dtype=np.float64)
        trial = np.squeeze(trial)
        if trial.shape == kinetic_grid.shape:
            values = trial
        else:
            values = np.asarray(
                [_evaluate_spectrum_point(surface_spectrum, float(t)) for t in kinetic_grid],
                dtype=np.float64,
            )
    else:
        values = np.asarray(surface_spectrum, dtype=np.float64)

    values = np.squeeze(values)
    if values.shape != kinetic_grid.shape:
        raise ValueError(
            "surface_spectrum must return one value per T-grid point, "
            f"expected shape {kinetic_grid.shape}, got {values.shape}."
        )
    return values


def _evaluate_spectrum_point(
    surface_spectrum: Callable[[FloatArray], FloatArray],
    kinetic_energy: float,
) -> float:
    try:
        value = surface_spectrum(np.array([kinetic_energy], dtype=np.float64))
    except Exception:
        value = surface_spectrum(float(kinetic_energy))  # type: ignore[arg-type]
    return float(np.squeeze(np.asarray(value, dtype=np.float64)))


def benchmark_earth_species() -> tuple[NuclearSpecies, ...]:
    species: list[NuclearSpecies] = []
    for name, mass_number, mass_fraction in BENCHMARK_CRUST_COMPOSITION:
        number_density = mass_fraction * BENCHMARK_EARTH_DENSITY / (mass_number * AMU)
        species.append(
            NuclearSpecies(
                name=name,
                mass_number=mass_number,
                mass=mass_number * AMU,
                number_density=number_density,
            )
        )
    return tuple(species)


def benchmark_velocity_distribution(velocity: float | FloatArray) -> float | FloatArray:
    v = np.asarray(velocity, dtype=np.float64)
    nesc = math.pi * BENCHMARK_VSUN**2 * (
        math.sqrt(math.pi) * BENCHMARK_VSUN * math.erf(BENCHMARK_VESC / BENCHMARK_VSUN)
        - 2.0 * BENCHMARK_VESC * math.exp(-(BENCHMARK_VESC / BENCHMARK_VSUN) ** 2)
    )
    exp_esc = math.exp(-(BENCHMARK_VESC / BENCHMARK_VSUN) ** 2)

    term_1 = 2.0 * np.exp(-(v * v + BENCHMARK_VEARTH * BENCHMARK_VEARTH) / (BENCHMARK_VSUN * BENCHMARK_VSUN))
    term_1 *= np.sinh(2.0 * v * BENCHMARK_VEARTH / (BENCHMARK_VSUN * BENCHMARK_VSUN))

    term_2 = (
        np.exp(-((v + BENCHMARK_VEARTH) ** 2) / (BENCHMARK_VSUN * BENCHMARK_VSUN)) - exp_esc
    ) * _heaviside(np.abs(v + BENCHMARK_VEARTH) - BENCHMARK_VESC)
    term_3 = (
        np.exp(-((v - BENCHMARK_VEARTH) ** 2) / (BENCHMARK_VSUN * BENCHMARK_VSUN)) - exp_esc
    ) * _heaviside(np.abs(v - BENCHMARK_VEARTH) - BENCHMARK_VESC)

    prefactor = math.pi * v * BENCHMARK_VSUN**2 / (BENCHMARK_VEARTH * nesc)
    result = prefactor * (term_1 + term_2 - term_3) * _heaviside(BENCHMARK_VEARTH + BENCHMARK_VESC - v)
    if np.isscalar(velocity):
        return float(np.asarray(result))
    return np.asarray(result, dtype=np.float64)


def benchmark_energy_distribution(m_chi: float, kinetic_energy: float | FloatArray) -> float | FloatArray:
    energy = np.asarray(kinetic_energy, dtype=np.float64)
    result = np.zeros_like(energy, dtype=np.float64)

    positive = energy > 0.0
    if np.any(positive):
        velocity = np.sqrt(2.0 * energy[positive] / m_chi)
        result[positive] = np.asarray(benchmark_velocity_distribution(velocity), dtype=np.float64) / np.sqrt(
            2.0 * m_chi * energy[positive]
        )

    if np.isscalar(kinetic_energy):
        return float(result.reshape(-1)[0])
    return result


def benchmark_gincoming(m_chi: float, kinetic_energy: float | FloatArray) -> float | FloatArray:
    energy = np.asarray(kinetic_energy, dtype=np.float64)
    velocity = np.zeros_like(energy, dtype=np.float64)
    positive = energy > 0.0
    if np.any(positive):
        velocity[positive] = np.sqrt(2.0 * energy[positive] / m_chi)

    unit_factor = CM**2 * SEC * EV
    result = unit_factor * 2.0 * math.pi**2 * velocity * BENCHMARK_DM_DENSITY / m_chi
    result *= np.asarray(benchmark_energy_distribution(m_chi, energy), dtype=np.float64)

    if np.isscalar(kinetic_energy):
        return float(result.reshape(-1)[0])
    return result


def make_benchmark_surface_spectrum(m_chi: float) -> Callable[[FloatArray], FloatArray]:
    def surface_spectrum(kinetic_grid: FloatArray) -> FloatArray:
        return np.asarray(benchmark_gincoming(m_chi, kinetic_grid), dtype=np.float64)

    return surface_spectrum


def make_benchmark_model(
    m_chi: float,
    cutoff_scale: float,
    *,
    earth_radius: float = BENCHMARK_EARTH_RADIUS,
) -> BenchmarkModel:
    return BenchmarkModel(
        m_chi=m_chi,
        cutoff_scale=cutoff_scale,
        earth_radius=earth_radius,
        species=benchmark_earth_species(),
    )


def reduced_mass(mass_1: float, mass_2: float) -> float:
    return mass_1 * mass_2 / (mass_1 + mass_2)


def sigma_cm2_to_natural(sigma_cm2: float) -> float:
    return sigma_cm2 * CM**2


def sigma_natural_to_cm2(sigma_natural: float) -> float:
    return sigma_natural / (CM**2)


def cutoff_scale_from_sigma_chi_n(m_chi: float, sigma_chi_n: float) -> float:
    mu_chi_n = reduced_mass(m_chi, NUCLEON_MASS)
    return (mu_chi_n * mu_chi_n / (math.pi * sigma_chi_n)) ** 0.25


def make_benchmark_model_from_sigma(
    m_chi: float,
    sigma_chi_n: float,
    *,
    sigma_in_cm2: bool = True,
    earth_radius: float = BENCHMARK_EARTH_RADIUS,
) -> BenchmarkModel:
    sigma_natural = sigma_cm2_to_natural(sigma_chi_n) if sigma_in_cm2 else sigma_chi_n
    cutoff_scale = cutoff_scale_from_sigma_chi_n(m_chi, sigma_natural)
    return make_benchmark_model(
        m_chi=m_chi,
        cutoff_scale=cutoff_scale,
        earth_radius=earth_radius,
    )


def kinetic_energy_from_speed(m_chi: float, speed: float | FloatArray) -> float | FloatArray:
    speed_arr = np.asarray(speed, dtype=np.float64)
    kinetic_energy = 0.5 * m_chi * speed_arr * speed_arr
    if np.isscalar(speed):
        return float(np.asarray(kinetic_energy))
    return kinetic_energy


def build_log_energy_grid(
    m_chi: float,
    n_t: int,
    *,
    t_min_fraction: float = 1.0e-3,
    vmax: float = BENCHMARK_VMAX,
) -> FloatArray:
    t_max = 0.5 * m_chi * vmax * vmax
    t_min = max(t_min_fraction * t_max, 1.0e-12 * EV)
    return np.geomspace(t_min, t_max, n_t, dtype=np.float64)


def build_quasi_uniform_speed_energy_grid(
    m_chi: float,
    n_t: int,
    *,
    vmax: float = BENCHMARK_VMAX,
    controller: float = 30.0,
    v_min_fraction: float = 1.0e-20,
) -> FloatArray:
    if n_t < 2:
        raise ValueError("n_t must be at least 2.")
    if controller <= 0.0:
        raise ValueError("controller must be positive.")
    if v_min_fraction < 0.0:
        raise ValueError("v_min_fraction must be non-negative.")

    log_points = np.linspace(math.log(controller), math.log(controller + 1.0), n_t, dtype=np.float64)
    v_points = vmax * (np.exp(log_points) - controller + v_min_fraction)
    v_points[0] = vmax * v_min_fraction
    v_points[-1] = vmax
    t_points = 0.5 * m_chi * np.square(v_points)
    return np.ascontiguousarray(t_points, dtype=np.float64)


def build_clenshaw_curtis_u_grid(n_u: int) -> FloatArray:
    if n_u < 2:
        raise ValueError("n_u must be at least 2.")
    theta = np.linspace(0.0, math.pi, n_u)
    return np.ascontiguousarray(-np.cos(theta), dtype=np.float64)


def clenshaw_curtis_weights(n_u: int) -> FloatArray:
    if n_u < 2:
        raise ValueError("n_u must be at least 2.")
    n_intervals = n_u - 1
    if n_intervals == 1:
        return np.ascontiguousarray(np.array([1.0, 1.0], dtype=np.float64))

    theta = math.pi * np.arange(n_u, dtype=np.float64) / n_intervals
    weights = np.zeros(n_u, dtype=np.float64)
    interior = np.arange(1, n_intervals)
    v = np.ones(n_intervals - 1, dtype=np.float64)

    if n_intervals % 2 == 0:
        weights[0] = 1.0 / (n_intervals * n_intervals - 1.0)
        weights[-1] = weights[0]
        for k in range(1, n_intervals // 2):
            v -= 2.0 * np.cos(2.0 * k * theta[interior]) / (4.0 * k * k - 1.0)
        v -= np.cos(n_intervals * theta[interior]) / (n_intervals * n_intervals - 1.0)
    else:
        weights[0] = 1.0 / (n_intervals * n_intervals)
        weights[-1] = weights[0]
        for k in range(1, (n_intervals + 1) // 2):
            v -= 2.0 * np.cos(2.0 * k * theta[interior]) / (4.0 * k * k - 1.0)

    weights[interior] = 2.0 * v / n_intervals
    return np.ascontiguousarray(weights[::-1], dtype=np.float64)


def build_radial_grid(
    earth_radius: float,
    *,
    mean_free_path: float,
    detector_depth: float = JINPING_DEPTH,
    main_tau_step: float = 0.2,
    tail_tau: float = 8.0,
    n_tail: int = 16,
    tail_power: float = 1.5,
    max_main_step: float = 50.0 * KM,
) -> FloatArray:
    if detector_depth <= 0.0 or detector_depth >= earth_radius:
        raise ValueError("detector_depth must lie strictly inside the Earth.")
    if mean_free_path <= 0.0:
        raise ValueError("mean_free_path must be positive.")
    if n_tail < 2:
        raise ValueError("n_tail must be at least 2.")
    if min(main_tau_step, tail_tau, tail_power, max_main_step) <= 0.0:
        raise ValueError("radial optical-depth parameters must be positive.")

    resolved_depth = min(earth_radius, detector_depth + tail_tau * mean_free_path)
    main_step = min(main_tau_step * mean_free_path, max_main_step)
    main_count = int(math.ceil(resolved_depth / main_step)) + 1
    main_depth = np.linspace(0.0, resolved_depth, main_count, endpoint=True, dtype=np.float64)

    depth_segments = [main_depth, np.array([detector_depth], dtype=np.float64)]
    if earth_radius > resolved_depth:
        t = np.linspace(0.0, 1.0, n_tail, endpoint=True, dtype=np.float64)
        depth_segments.append(resolved_depth + (earth_radius - resolved_depth) * t**tail_power)

    depth_grid = np.unique(np.round(np.concatenate(depth_segments), decimals=12))
    grid = np.unique(np.clip(earth_radius - depth_grid, 0.0, earth_radius))
    return np.ascontiguousarray(grid, dtype=np.float64)


def sigma_nr_for_species(model: BenchmarkModel, species: NuclearSpecies) -> float:
    mu = model.m_chi * species.mass / (model.m_chi + species.mass)
    return species.mass_number**2 * mu * mu / (math.pi * model.cutoff_scale**4)


def mean_free_path_for_model(model: BenchmarkModel) -> float:
    inv_mfp = 0.0
    for species in model.species:
        inv_mfp += species.number_density * sigma_nr_for_species(model, species)
    if inv_mfp <= 0.0:
        raise ValueError("The mean free path must be positive.")
    return 1.0 / inv_mfp


def nonuniform_trapezoid_weights(grid_values: Sequence[float] | FloatArray) -> FloatArray:
    grid = np.asarray(grid_values, dtype=np.float64)
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError("grid_values must be a 1D array with at least two points.")
    weights = np.empty_like(grid)
    weights[0] = 0.5 * (grid[1] - grid[0])
    weights[-1] = 0.5 * (grid[-1] - grid[-2])
    if grid.size > 2:
        weights[1:-1] = 0.5 * (grid[2:] - grid[:-2])
    return weights


def angular_average_weights(u_values: Sequence[float] | FloatArray) -> FloatArray:
    grid = _as_1d_array(u_values, "u_values")
    cc_grid = build_clenshaw_curtis_u_grid(grid.size)
    if np.allclose(grid, cc_grid, rtol=1.0e-12, atol=1.0e-14):
        weights = clenshaw_curtis_weights(grid.size)
    else:
        weights = nonuniform_trapezoid_weights(grid)
    return np.ascontiguousarray(0.5 * weights, dtype=np.float64)


def detector_spectrum(
    field: FloatArray,
    r_grid: Sequence[float] | FloatArray,
    u_grid: Sequence[float] | FloatArray,
    detector_radius: float,
) -> FloatArray:
    field_arr = np.asarray(field, dtype=np.float64)
    r_values = np.asarray(r_grid, dtype=np.float64)
    u_values = np.asarray(u_grid, dtype=np.float64)
    if field_arr.ndim != 3:
        raise ValueError("field must have shape (N_r, N_u, N_T).")
    if field_arr.shape[0] != r_values.size or field_arr.shape[1] != u_values.size:
        raise ValueError("field shape does not match r_grid/u_grid.")

    r_stencil = _linear_stencil(r_values, np.array([detector_radius], dtype=np.float64))
    detector_lo = int(r_stencil.lo[0])
    detector_hi = int(r_stencil.hi[0])
    rw_lo = float(r_stencil.w_lo[0])
    rw_hi = float(r_stencil.w_hi[0])

    detector_slice = rw_lo * field_arr[detector_lo] + rw_hi * field_arr[detector_hi]
    u_weights = angular_average_weights(u_values)
    return np.ascontiguousarray(np.tensordot(u_weights, detector_slice, axes=(0, 0)), dtype=np.float64)


def _heaviside(values: FloatArray) -> FloatArray:
    arr = np.asarray(values, dtype=np.float64)
    return np.where(arr >= 0.0, 1.0, 0.0)
