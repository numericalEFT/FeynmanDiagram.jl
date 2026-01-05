import numpy as np
from scipy.special import ellipk
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib.pyplot as plt


class HubbardModel2D_U0:
    def __init__(self, t=1.0, beta=10.0):
        """
        Initialize the 2D Square Lattice Hubbard Model (U=0) at finite T.

        Args:
            t (float): Hopping parameter
            beta (float): Inverse temperature (1/kT)
        """
        self.t = t
        self.beta = beta
        self.bandwidth = 4 * t  # Energy ranges from -4t to 4t

    def density_of_states(self, E):
        """
        Analytic DOS for 2D square lattice using Elliptic K.
        rho(E) = (1 / 2*pi^2*t) * K(1 - (E/4t)^2)
        """
        if abs(E) >= self.bandwidth:
            return 0.0

        # Avoid exact singularity at E=0 for stability, though quad handles it well
        if abs(E) < 1e-13:
            return np.inf

        term = (E / (4 * self.t)) ** 2
        return (1.0 / (2 * np.pi**2 * self.t)) * ellipk(1 - term)

    def fermi_function(self, E, mu):
        """
        Fermi-Dirac distribution: f(E) = 1 / (exp(beta*(E-mu)) + 1)
        """
        # Numerically stable implementation to avoid overflow
        x = self.beta * (E - mu)
        if x > 100:
            return 0.0
        elif x < -100:
            return 1.0
        else:
            return 1.0 / (np.exp(x) + 1.0)

    def get_n_from_mu(self, mu):
        """
        Calculate density n given chemical potential mu at finite T.
        n = 2 * integral( rho(E) * f(E) ) dE
        """

        # Integrand: DOS(E) * Fermi(E)
        def integrand(E):
            return self.density_of_states(E) * self.fermi_function(E, mu)

        # We integrate over the entire bandwidth [-4t, 4t].
        # For very low T (high beta), the Fermi function is a step function.
        # quad usually handles this well, but providing 'points' at mu helps
        # the integrator locate the step.

        points = [0.0]  # Tell integrator about the DOS singularity at 0
        if -self.bandwidth < mu < self.bandwidth:
            points.append(mu)  # Tell integrator about the Fermi step location

        val, err = quad(integrand, -self.bandwidth, self.bandwidth, points=points)

        return 2.0 * val  # Factor of 2 for spin degeneracy

    def get_mu_from_n(self, target_n):
        """
        Inverse function: Calculate mu given density n.
        Uses Brent's method to find root of n(mu) - target = 0.
        """
        if target_n <= 0.0:
            return -np.inf
        if target_n >= 2.0:
            return np.inf

        # Optimization bounds
        # At T=0, bounds are strictly [-4t, 4t].
        # At finite T, mu can theoretically be outside the band to get n near 0 or 2.
        # We expand search bounds slightly based on T.
        bound_buffer = 10.0 / self.beta if self.beta > 0.1 else 100.0
        lower_bound = -self.bandwidth - bound_buffer
        upper_bound = self.bandwidth + bound_buffer

        def objective(mu):
            return self.get_n_from_mu(mu) - target_n

        try:
            mu_sol = brentq(objective, lower_bound, upper_bound)
            return mu_sol
        except ValueError:
            # Fallback for extreme temperatures or densities
            print("Warning: Standard bounds failed, trying expanded search...")
            return brentq(objective, lower_bound * 2, upper_bound * 2)


# --- Example Usage ---
if __name__ == "__main__":
    beta = 8.0
    model = HubbardModel2D_U0(t=1.0, beta=beta)

    print(f"--- 2D Hubbard Model (U=0, beta={beta}) ---")

    # 1. Forward: mu -> n
    # At mu=0, n should be exactly 1.0 (Particle-hole symmetry)
    mu_test = 0.0
    n_res = model.get_n_from_mu(mu_test)
    print(f"mu = {mu_test:.4f} -> n = {n_res:.6f}")

    # 2. Inverse: n -> mu
    n_target = 0.875
    mu_sol = model.get_mu_from_n(n_target)
    print(f"n = {n_target:.4f} -> mu = {mu_sol:.6f}")

    # Verification
    n_check = model.get_n_from_mu(mu_sol)
    print(f"Check: n(mu_sol) = {n_check:.6f}")

    # --- Plotting ---
    mu_vals = np.linspace(-6, 6, 100)
    n_vals = [model.get_n_from_mu(mu) for mu in mu_vals]

    plt.figure(figsize=(8, 6))
    plt.plot(mu_vals, n_vals, label=f"$\\beta={beta}$")
    plt.xlabel(r"$\mu/t$")
    plt.ylabel(r"$n$")
    plt.title(r"2D Hubbard $U=0$: Density vs Chemical Potential")
    plt.grid(True, alpha=0.3)
    plt.axvline(0, color="k", linestyle="--", alpha=0.3)
    plt.axhline(1, color="k", linestyle="--", alpha=0.3)
    plt.legend()
    plt.show()
