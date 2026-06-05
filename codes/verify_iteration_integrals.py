from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np

import iteration_dirac_dm_isoscalar_vector_coupling as mod


@dataclass(frozen=True)
class ValidationResult:
    name: str
    value: float
    tolerance: float

    @property
    def passed(self) -> bool:
        return self.value <= self.tolerance


def build_validation_solver() -> mod.BenchmarkIterationSolver:
    m_chi = 100.0 * mod.MEV
    cutoff_scale = 1.0e3 * mod.GEV
    model = mod.make_wl_benchmark_model(m_chi=m_chi, cutoff_scale=cutoff_scale)
    grid = mod.PhaseSpaceGrid(
        r=np.linspace(0.0, model.earth_radius, 9),
        u=np.linspace(-1.0, 1.0, 7),
        T=np.geomspace(1.0 * mod.EV, 1.0e5 * mod.EV, 11),
    )
    quadrature = mod.QuadratureSpec(n_x=6, n_y=6, n_z=8)
    return mod.BenchmarkIterationSolver(model, grid, quadrature, strict=False)


def smooth_field_function(
    r: np.ndarray | float,
    u: np.ndarray | float,
    kinetic_energy: np.ndarray | float,
    *,
    earth_radius: float,
    energy_scale: float,
) -> np.ndarray:
    r_arr = np.asarray(r, dtype=np.float64)
    u_arr = np.asarray(u, dtype=np.float64)
    t_arr = np.asarray(kinetic_energy, dtype=np.float64)
    return np.exp(-0.35 * r_arr / earth_radius) * (1.0 + 0.2 * u_arr) * (1.0 + 0.15 * t_arr / energy_scale)


def multilinear_field_function(
    r: np.ndarray | float,
    u: np.ndarray | float,
    kinetic_energy: np.ndarray | float,
    *,
    earth_radius: float,
    energy_scale: float,
) -> np.ndarray:
    rr = np.asarray(r, dtype=np.float64) / earth_radius
    uu = np.asarray(u, dtype=np.float64)
    tt = np.asarray(kinetic_energy, dtype=np.float64) / energy_scale
    return (
        1.0
        + 0.3 * rr
        - 0.25 * uu
        + 0.2 * tt
        + 0.1 * rr * uu
        - 0.08 * rr * tt
        + 0.07 * uu * tt
        + 0.05 * rr * uu * tt
    )


def path_quadrature_normalization(solver: mod.BenchmarkIterationSolver) -> ValidationResult:
    eta = np.maximum(0.0, (solver._x_grid + solver._x_oplus) / solver.mean_free_path)
    eta_resolved = np.minimum(eta, solver.tau_segment_edges[-1])
    exact = solver.mean_free_path * (1.0 - np.exp(-eta_resolved))
    approx = np.sum(solver._omega_x, axis=2)
    scale = np.maximum(np.abs(exact), 1.0)
    rel_error = np.abs(approx - exact) / scale
    return ValidationResult("path quadrature normalization", float(np.max(rel_error)), 3.0e-3)


def y_quadrature_exactness(solver: mod.BenchmarkIterationSolver) -> ValidationResult:
    max_error = 0.0
    max_degree = 2 * solver.quadrature.n_y - 1

    for cache in solver._species_cache:
        T = solver.grid.T[:, None]
        T_prime = cache.t_prime
        y = 1.0 - T / T_prime

        for degree in range(max_degree + 1):
            integrand = (T / (T_prime * T_prime)) * y**degree
            approx = np.sum(cache.omega_y * integrand, axis=1)
            exact = cache.kappa ** (degree + 1) / (degree + 1) * np.ones_like(approx)
            max_error = max(max_error, float(np.max(np.abs(approx - exact))))

    return ValidationResult("y quadrature exactness", max_error, 1.0e-11)


def z_quadrature_exactness(solver: mod.BenchmarkIterationSolver) -> ValidationResult:
    max_error = 0.0
    max_degree = 2 * solver.quadrature.n_z - 1
    z_nodes = solver._z_nodes
    weight = solver._omega_z

    for degree in range(max_degree + 1):
        approx = weight * float(np.sum(z_nodes**degree))
        exact = chebyshev_moment_exact(degree)
        max_error = max(max_error, abs(approx - exact))

    return ValidationResult("z quadrature exactness", max_error, 1.0e-12)


def u_clenshaw_curtis_exactness() -> ValidationResult:
    max_error = 0.0
    for n_u in (8, 20, 40):
        u_nodes = mod.build_clenshaw_curtis_u_grid(n_u)
        weights = mod.clenshaw_curtis_weights(n_u)
        for degree in range(n_u):
            approx = float(np.sum(weights * u_nodes**degree))
            if degree % 2 == 0:
                exact = 2.0 / (degree + 1)
            else:
                exact = 0.0
            max_error = max(max_error, abs(approx - exact))

    return ValidationResult("u Clenshaw-Curtis exactness", max_error, 1.0e-12)


def radial_grid_small_lambda_spacing() -> ValidationResult:
    model = mod.make_wl_benchmark_model_from_sigma(
        m_chi=5000.0 * mod.MEV,
        sigma_chi_n=5.0e-32,
        sigma_in_cm2=True,
    )
    mean_free_path = mod.mean_free_path_for_model(model)
    detector_depth = 2.4 * mod.KM
    main_tau_step = 0.2
    tail_tau = 8.0
    r_grid = mod.build_radial_grid(
        model.earth_radius,
        mean_free_path=mean_free_path,
        detector_depth=detector_depth,
        main_tau_step=main_tau_step,
        tail_tau=tail_tau,
        n_tail=16,
        max_main_step=50.0 * mod.KM,
    )
    depth = np.sort((model.earth_radius - r_grid) / mean_free_path)
    resolved_depth = detector_depth / mean_free_path + tail_tau
    resolved = (depth[:-1] >= -1.0e-12) & (depth[1:] <= resolved_depth + 1.0e-12)
    max_gap = float(np.max(np.diff(depth)[resolved]))
    return ValidationResult("radial grid small-lambda spacing", max(0.0, max_gap - main_tau_step), 1.0e-12)


def radial_grid_large_lambda_spacing() -> ValidationResult:
    earth_radius = mod.WL_EARTH_RADIUS
    detector_depth = 2.4 * mod.KM
    mean_free_path = 10000.0 * mod.KM
    tail_tau = 8.0
    max_main_step = 50.0 * mod.KM
    r_grid = mod.build_radial_grid(
        earth_radius,
        mean_free_path=mean_free_path,
        detector_depth=detector_depth,
        main_tau_step=0.2,
        tail_tau=tail_tau,
        n_tail=16,
        max_main_step=max_main_step,
    )
    depth = np.sort(earth_radius - r_grid)
    resolved_depth = min(earth_radius, detector_depth + tail_tau * mean_free_path)
    resolved = (depth[:-1] >= -1.0e-12) & (depth[1:] <= resolved_depth + 1.0e-12)
    max_gap = float(np.max(np.diff(depth)[resolved]))
    return ValidationResult("radial grid large-lambda spacing", max(0.0, max_gap - max_main_step), 1.0e-8)


def chebyshev_moment_exact(degree: int) -> float:
    if degree % 2 == 1:
        return 0.0
    half_degree = degree // 2
    return math.pi * math.comb(2 * half_degree, half_degree) / (4.0**half_degree)


def direct_discrete_iteration_exact(
    solver: mod.BenchmarkIterationSolver,
    field_function,
) -> np.ndarray:
    nr, nu, nt = solver.shape
    result = np.zeros((nr, nu, nt), dtype=np.float64)

    for a, r_value in enumerate(solver.grid.r):
        for b, u_value in enumerate(solver.grid.u):
            for c, kinetic_energy in enumerate(solver.grid.T):
                accum = 0.0
                for cache in solver._species_cache:
                    for k in range(solver.quadrature.n_x):
                        weight_x = solver._omega_x[a, b, k]
                        if weight_x == 0.0:
                            continue

                        r_star = solver._r_star[a, b, k]
                        for m in range(solver.quadrature.n_y):
                            weight = (
                                cache.species.number_density
                                * weight_x
                                * cache.omega_y[c, m]
                                * cache.kernel[c, m]
                                * solver._omega_z
                            )
                            if weight == 0.0:
                                continue

                            c_val = cache.c_angle[c, m]
                            radius = math.sqrt(max(0.0, 1.0 - u_value * u_value))
                            radius *= math.sqrt(max(0.0, 1.0 - c_val * c_val))
                            t_prime = cache.t_prime[c, m]
                            u_samples = u_value * c_val + radius * solver._z_nodes
                            u_samples = np.clip(u_samples, solver.grid.u[0], solver.grid.u[-1])
                            values = field_function(
                                np.full_like(u_samples, r_star),
                                u_samples,
                                np.full_like(u_samples, t_prime),
                                earth_radius=solver.model.earth_radius,
                                energy_scale=solver.grid.T[-1],
                            )
                            accum += weight * float(np.sum(values))

                result[a, b, c] = accum

    return result


def interpolation_exactness(solver: mod.BenchmarkIterationSolver) -> ValidationResult:
    rr = solver.grid.r[:, None, None]
    uu = solver.grid.u[None, :, None]
    tt = solver.grid.T[None, None, :]
    field = multilinear_field_function(
        rr,
        uu,
        tt,
        earth_radius=solver.model.earth_radius,
        energy_scale=solver.grid.T[-1],
    )
    iterated = solver.apply_iteration(field)
    exact = direct_discrete_iteration_exact(solver, multilinear_field_function)
    error = float(np.max(np.abs(iterated - exact)))
    return ValidationResult("interpolation exactness", error, 1.0e-10)


def triple_integral_reference(solver: mod.BenchmarkIterationSolver) -> ValidationResult:
    a = solver.grid.r.size // 2
    b = solver.grid.u.size // 2
    c = solver.grid.T.size // 2

    discrete = direct_discrete_iteration_point(
        solver,
        a,
        b,
        c,
        smooth_field_function,
    )
    reference = high_order_reference_point(
        solver,
        a,
        b,
        c,
        smooth_field_function,
        n_x_ref=36,
        n_y_ref=36,
        n_z_ref=80,
    )
    abs_error = abs(discrete - reference)
    rel_error = abs_error / max(abs(reference), 1.0e-30)
    return ValidationResult("triple integral reference", max(abs_error, rel_error), 1.0e-6)


def direct_discrete_iteration_point(
    solver: mod.BenchmarkIterationSolver,
    a: int,
    b: int,
    c: int,
    field_function,
) -> float:
    u_value = solver.grid.u[b]
    accum = 0.0

    for cache in solver._species_cache:
        for k in range(solver.quadrature.n_x):
            weight_x = solver._omega_x[a, b, k]
            if weight_x == 0.0:
                continue

            r_star = solver._r_star[a, b, k]
            for m in range(solver.quadrature.n_y):
                weight = (
                    cache.species.number_density
                    * weight_x
                    * cache.omega_y[c, m]
                    * cache.kernel[c, m]
                    * solver._omega_z
                )
                if weight == 0.0:
                    continue

                c_val = cache.c_angle[c, m]
                radius = math.sqrt(max(0.0, 1.0 - u_value * u_value))
                radius *= math.sqrt(max(0.0, 1.0 - c_val * c_val))
                t_prime = cache.t_prime[c, m]
                u_samples = u_value * c_val + radius * solver._z_nodes
                u_samples = np.clip(u_samples, solver.grid.u[0], solver.grid.u[-1])
                values = field_function(
                    np.full_like(u_samples, r_star),
                    u_samples,
                    np.full_like(u_samples, t_prime),
                    earth_radius=solver.model.earth_radius,
                    energy_scale=solver.grid.T[-1],
                )
                accum += weight * float(np.sum(values))

    return accum


def high_order_reference_point(
    solver: mod.BenchmarkIterationSolver,
    a: int,
    b: int,
    c: int,
    field_function,
    *,
    n_x_ref: int,
    n_y_ref: int,
    n_z_ref: int,
) -> float:
    r = solver.grid.r[a]
    u = solver.grid.u[b]
    kinetic_energy = solver.grid.T[c]
    x_coord = r * u
    y_coord = r * math.sqrt(max(0.0, 1.0 - u * u))
    eta = max(0.0, (solver._x_grid[a, b] + solver._x_oplus[a, b]) / solver.mean_free_path)
    eta_resolved = min(float(eta), float(solver.tau_segment_edges[-1]))
    q_max = 1.0 - math.exp(-eta_resolved)

    q_nodes, q_weights = np.polynomial.legendre.leggauss(n_x_ref)
    q_values = 0.5 * q_max * (1.0 + q_nodes)
    x_pp = -solver.mean_free_path * np.log1p(-q_values)
    omega_x = solver.mean_free_path * 0.5 * q_max * q_weights

    z_nodes = np.cos((2.0 * np.arange(1, n_z_ref + 1) - 1.0) * math.pi / (2.0 * n_z_ref))
    omega_z = math.pi / n_z_ref

    total = 0.0
    for species in solver.model.species:
        kappa = 4.0 * solver.model.m_chi * species.mass / (solver.model.m_chi + species.mass) ** 2
        y_nodes, y_weights = np.polynomial.legendre.leggauss(n_y_ref)
        y_values = 0.5 * kappa * (1.0 + y_nodes)
        t_prime = kinetic_energy / (1.0 - y_values)
        omega_y = 0.5 * kappa * y_weights * kinetic_energy / (1.0 - y_values) ** 2
        t_array = np.full_like(t_prime, kinetic_energy)
        kernel = solver._tilde_kernel(species, t_array, t_prime)
        c_angle = solver._c_angle(species, t_array, t_prime)
        radius_prefactor = math.sqrt(max(0.0, 1.0 - u * u))

        for x_weight, x_shift in zip(omega_x, x_pp, strict=True):
            r_star = math.sqrt(max(0.0, (x_coord - x_shift) ** 2 + y_coord * y_coord))
            for m in range(n_y_ref):
                base_weight = species.number_density * x_weight * omega_y[m] * kernel[m] * omega_z
                if base_weight == 0.0:
                    continue

                c_val = float(c_angle[m])
                radius = radius_prefactor * math.sqrt(max(0.0, 1.0 - c_val * c_val))
                u_samples = u * c_val + radius * z_nodes
                u_samples = np.clip(u_samples, solver.grid.u[0], solver.grid.u[-1])
                values = field_function(
                    np.full_like(u_samples, r_star),
                    u_samples,
                    np.full_like(u_samples, t_prime[m]),
                    earth_radius=solver.model.earth_radius,
                    energy_scale=solver.grid.T[-1],
                )
                total += base_weight * float(np.sum(values))

    return total


def print_results(results: list[ValidationResult]) -> None:
    width = max(len(result.name) for result in results)
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(
            f"{status}  {result.name:<{width}}  "
            f"error={result.value:.3e}  tol={result.tolerance:.3e}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the transformed iteration integrals.")
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="Skip the high-order triple-integral reference check.",
    )
    args = parser.parse_args()

    solver = build_validation_solver()
    results = [
        path_quadrature_normalization(solver),
        y_quadrature_exactness(solver),
        z_quadrature_exactness(solver),
        u_clenshaw_curtis_exactness(),
        radial_grid_small_lambda_spacing(),
        radial_grid_large_lambda_spacing(),
        interpolation_exactness(solver),
    ]
    if not args.skip_reference:
        results.append(triple_integral_reference(solver))

    print_results(results)
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
