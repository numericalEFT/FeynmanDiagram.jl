import numpy as np
from itertools import product  # Not strictly needed with the new basis generation


# Helper function to count set bits (electrons) in an integer representation
def count_set_bits(n):
    """Counts the number of set bits in an integer n."""
    count = 0
    while n > 0:
        n &= n - 1
        count += 1
    return count


class HubbardEDGrandCanonical:
    """
    Exact Diagonalization for the 1D Hubbard model in the Grand Canonical Ensemble.
    """

    def __init__(self, L, t, U, mu, pbc=True):
        """
        Initializes the Hubbard model.

        Args:
            L (int): Number of sites.
            t (float): Hopping parameter.
            U (float): On-site interaction strength.
            mu (float): Chemical potential.
            pbc (bool): If True, use periodic boundary conditions. Otherwise, open.
        """
        self.L = L
        self.t = t
        self.U = U
        self.mu = mu
        self.pbc = pbc

        self.basis_states = []  # List of (s_up, s_down) tuples
        self.basis_map = {}  # Dictionary: state_tuple -> index
        # Stores (N_up, N_down) for each basis state, corresponding to self.basis_states by index
        self.particle_numbers_for_basis_states = []

        self._generate_basis()

        if not self.basis_states:  # Should only happen if L=0
            print(f"Warning: No basis states generated for L={self.L}.")
            self.H_prime = np.array([])
            self.eigenvalues_prime = np.array([])
        else:
            # H_prime is H_0 - mu * N_tot
            self.H_prime = np.zeros(
                (len(self.basis_states), len(self.basis_states)), dtype=np.float64
            )
            self._build_hamiltonian()
            self.eigenvalues_prime = None  # To be computed by diagonalize()

    def _generate_basis(self):
        """Generates all possible basis states (s_up, s_down) for all particle numbers."""
        idx = 0
        # Iterate over all possible spin-up configurations (0 to 2^L - 1)
        for s_up_config in range(1 << self.L):
            # Iterate over all possible spin-down configurations
            for s_down_config in range(1 << self.L):
                state_tuple = (s_up_config, s_down_config)
                self.basis_states.append(state_tuple)
                self.basis_map[state_tuple] = idx

                n_up = count_set_bits(s_up_config)
                n_down = count_set_bits(s_down_config)
                self.particle_numbers_for_basis_states.append((n_up, n_down))
                idx += 1

    def _get_fermionic_sign(self, config, an_site, cr_site):
        """
        Calculates the fermionic sign for c_cr_site^dag * c_an_site acting on config.
        Sites are 0-indexed.
        """
        if an_site == cr_site:
            return 1

        sign = 1
        mask_before_an = (1 << an_site) - 1
        if (config >> an_site) & 1:
            sign = (-1) ** count_set_bits(config & mask_before_an)
        else:
            return 0

        temp_config = config ^ (1 << an_site)
        mask_before_cr = (1 << cr_site) - 1
        sign *= (-1) ** count_set_bits(temp_config & mask_before_cr)

        return sign

    def _build_hamiltonian(self):
        """Constructs the Grand Canonical Hamiltonian matrix H' = H_0 - mu * N_tot."""
        if not self.basis_states:
            return

        for k_idx_initial, (s_up_initial, s_down_initial) in enumerate(
            self.basis_states
        ):
            N_up_initial, N_down_initial = self.particle_numbers_for_basis_states[
                k_idx_initial
            ]
            N_total_initial = N_up_initial + N_down_initial

            # --- Diagonal Part: On-site interaction (U) and Chemical Potential (mu) ---
            diag_energy = 0
            # U term
            for site in range(self.L):
                is_up_occupied = (s_up_initial >> site) & 1
                is_down_occupied = (s_down_initial >> site) & 1
                if is_up_occupied and is_down_occupied:
                    diag_energy += self.U

            # mu term: H_prime includes -mu * N_total
            diag_energy -= self.mu * N_total_initial

            self.H_prime[k_idx_initial, k_idx_initial] += diag_energy

            # --- Off-diagonal Part: Hopping term (-t) ---
            # Hopping conserves N_up and N_down individually, so it connects states
            # that would have been in the same (N_up, N_down) block.
            for site_i in range(self.L):
                neighbors_j = []
                if self.pbc:
                    neighbors_j.append((site_i - 1 + self.L) % self.L)
                    neighbors_j.append((site_i + 1) % self.L)
                    if self.L == 2:
                        neighbors_j = list(set(neighbors_j))
                else:  # OBC
                    if site_i > 0:
                        neighbors_j.append(site_i - 1)
                    if site_i < self.L - 1:
                        neighbors_j.append(site_i + 1)

                for site_j in neighbors_j:
                    if site_i == site_j:
                        continue

                    # Spin-up hopping
                    if ((s_up_initial >> site_i) & 1) and not (
                        (s_up_initial >> site_j) & 1
                    ):
                        fermionic_sign = self._get_fermionic_sign(
                            s_up_initial, site_i, site_j
                        )
                        s_up_final = s_up_initial ^ (1 << site_i) | (1 << site_j)
                        final_state_tuple = (s_up_final, s_down_initial)
                        k_idx_final = self.basis_map[final_state_tuple]
                        self.H_prime[k_idx_final, k_idx_initial] += (
                            -self.t * fermionic_sign
                        )

                    # Spin-down hopping
                    if ((s_down_initial >> site_i) & 1) and not (
                        (s_down_initial >> site_j) & 1
                    ):
                        fermionic_sign = self._get_fermionic_sign(
                            s_down_initial, site_i, site_j
                        )
                        s_down_final = s_down_initial ^ (1 << site_i) | (1 << site_j)
                        final_state_tuple = (s_up_initial, s_down_final)
                        k_idx_final = self.basis_map[final_state_tuple]
                        self.H_prime[k_idx_final, k_idx_initial] += (
                            -self.t * fermionic_sign
                        )

    def diagonalize(self):
        """Diagonalizes the H_prime Hamiltonian."""
        if len(self.basis_states) == 0:
            self.eigenvalues_prime = np.array([])
            return self.eigenvalues_prime

        self.eigenvalues_prime, _ = np.linalg.eigh(self.H_prime)
        return self.eigenvalues_prime

    def get_grand_potential(self, T, k_B=1.0):
        """
        Calculates the Grand Potential Omega = -k_B T ln(Z_G).
        Z_G = sum(exp(-E'_m / (k_B T))), where E'_m are eigenvalues of H_prime.
        """
        if self.eigenvalues_prime is None:
            self.diagonalize()

        if len(self.eigenvalues_prime) == 0:
            return np.nan

        if T == 0:
            # At T=0, Omega is the minimum eigenvalue of H_prime
            return (
                np.min(self.eigenvalues_prime)
                if self.eigenvalues_prime.size > 0
                else np.nan
            )

        beta = 1.0 / (k_B * T)
        min_E_prime = np.min(self.eigenvalues_prime)
        scaled_eigenvalues_prime = self.eigenvalues_prime - min_E_prime

        partition_function_ZG_scaled_sum = np.sum(
            np.exp(-beta * scaled_eigenvalues_prime)
        )

        if (
            partition_function_ZG_scaled_sum <= 0
        ):  # Handles underflow or all E very large
            return np.inf if T > 0 else min_E_prime  # At T=0, Omega is min E'

        log_ZG_scaled_sum = np.log(partition_function_ZG_scaled_sum)
        Omega = min_E_prime - (1.0 / beta) * log_ZG_scaled_sum

        return Omega

    def get_average_total_particle_number(self, T, k_B=1.0):
        """
        Calculates the average total number of particles <N_tot>.
        Requires diagonalizing H_prime and also its eigenvectors.
        <N_tot> = (1/Z_G) * sum_m <psi_m| N_tot |psi_m> * exp(-beta * E'_m)
        Alternatively, <N_tot> = -d(Omega)/d(mu)
        Here we use the direct summation with eigenvectors.
        """
        if (
            self.eigenvalues_prime is None
            or self.H_prime is None
            or len(self.basis_states) == 0
        ):
            self.diagonalize()  # This only gets eigenvalues currently by default.
            # We need eigenvectors too for this method.

        # Re-diagonalize to get eigenvectors if not already stored
        if len(self.basis_states) == 0:
            return np.nan
        eigenvalues_prime, eigenvectors_prime = np.linalg.eigh(self.H_prime)

        if T == 0:
            ground_state_idx = np.argmin(eigenvalues_prime)
            ground_state_vector = eigenvectors_prime[:, ground_state_idx]
            avg_N = 0
            for i, basis_vec_coeffs in enumerate(ground_state_vector):
                n_up, n_down = self.particle_numbers_for_basis_states[i]
                avg_N += (basis_vec_coeffs**2) * (n_up + n_down)
            return avg_N

        beta = 1.0 / (k_B * T)
        min_E_prime = np.min(eigenvalues_prime)

        Z_G_terms = np.exp(-beta * (eigenvalues_prime - min_E_prime))
        Z_G = np.sum(Z_G_terms)

        if Z_G == 0:
            return np.nan

        avg_N_numerator = 0
        for m_idx in range(len(eigenvalues_prime)):  # Sum over eigenstates of H'
            E_prime_m = eigenvalues_prime[m_idx]
            psi_m = eigenvectors_prime[:, m_idx]  # m-th eigenvector

            N_expectation_for_psi_m = 0
            for k_basis_idx, component_psi_mk in enumerate(psi_m):
                N_up_k, N_down_k = self.particle_numbers_for_basis_states[k_basis_idx]
                N_total_k = N_up_k + N_down_k
                N_expectation_for_psi_m += (
                    component_psi_mk**2
                ) * N_total_k  # Assumes real eigenvectors

            avg_N_numerator += N_expectation_for_psi_m * np.exp(
                -beta * (E_prime_m - min_E_prime)
            )

        return avg_N_numerator / Z_G

    def get_particle_density(self, T, site=None, k_B=1.0):
        """
        Calculates the particle density.

        If site is None, returns the total average density <N_tot>/L.
        If site is specified, returns the local density <n_i> at that site.

        Args:
            T (float): Temperature.
            site (int, optional): The site index (0 to L-1). Defaults to None.
            k_B (float): Boltzmann constant. Defaults to 1.0.

        Returns:
            float: The calculated particle density.
        """
        # Case 1: Calculate total average particle density
        if site is None:
            avg_N = self.get_average_total_particle_number(T, k_B)
            return avg_N / self.L if self.L > 0 else 0

        # Case 2: Calculate local particle density at a specific site
        if not (0 <= site < self.L):
            raise ValueError(f"Site index {site} is out of bounds for L={self.L}.")

        eigenvalues_prime, eigenvectors_prime = self.diagonalize()
        if eigenvalues_prime.size == 0:
            return np.nan

        # T=0: Calculate expectation value for the ground state
        if T == 0:
            gs_vector = eigenvectors_prime[:, np.argmin(eigenvalues_prime)]
            n_i_gs = 0
            for k_basis_idx, vec_comp in enumerate(gs_vector):
                s_up, s_down = self.basis_states[k_basis_idx]
                # Occupation of site i = (n_i_up + n_i_down)
                n_i_k = ((s_up >> site) & 1) + ((s_down >> site) & 1)
                n_i_gs += (vec_comp**2) * n_i_k
            return n_i_gs

        # T>0: Calculate the thermal average
        beta = 1.0 / (k_B * T)
        min_E_prime = np.min(eigenvalues_prime)
        Z_G_terms = np.exp(-beta * (eigenvalues_prime - min_E_prime))
        Z_G = np.sum(Z_G_terms)
        if Z_G == 0:
            return np.nan

        avg_n_i_numerator = 0
        for m_idx in range(len(eigenvalues_prime)):  # Sum over all eigenstates m
            psi_m = eigenvectors_prime[:, m_idx]
            n_i_expect_m = 0  # Expectation value of n_i for eigenstate m
            for k_basis_idx, component in enumerate(psi_m):
                s_up_k, s_down_k = self.basis_states[k_basis_idx]
                n_i_k = ((s_up_k >> site) & 1) + ((s_down_k >> site) & 1)
                n_i_expect_m += (component**2) * n_i_k

            avg_n_i_numerator += n_i_expect_m * Z_G_terms[m_idx]

        return avg_n_i_numerator / Z_G


# --- Example Usage ---
if __name__ == "__main__":
    L_sites = 3
    t_hopping = 1.0
    U_interaction = 1.0
    # mu_chemical_potential = 1.0
    mu_chemical_potential = 0.5
    # mu=U/2 often corresponds to half-filling in large systems. For small L, it's more complex.

    print(f"Hubbard Model Parameters (Grand Canonical):")
    print(f"Number of sites (L): {L_sites}")
    print(f"Hopping (t): {t_hopping}")
    print(f"On-site interaction (U): {U_interaction}")
    print(f"Chemical Potential (mu): {mu_chemical_potential}")
    print("-" * 40)

    hubbard_model_gc = HubbardEDGrandCanonical(
        # L_sites, t_hopping, U_interaction, mu_chemical_potential, pbc=True
        L_sites,
        t_hopping,
        U_interaction,
        mu_chemical_potential,
        # pbc=False,
        pbc=True,
    )

    print(hubbard_model_gc.H_prime)

    if len(hubbard_model_gc.basis_states) > 0:
        print(
            f"Total number of basis states (4^L): {len(hubbard_model_gc.basis_states)}"
        )

        eigenvalues_gc = hubbard_model_gc.diagonalize()
        # print(
        #     f"Lowest few eigenvalues of H' (H_0 - mu*N): {eigenvalues_gc[: min(5, len(eigenvalues_gc))]}"
        # )
        # print(f"Eigenvalues of H' (H_0 - mu*N): {eigenvalues_gc}")

        Temperature = 0.2
        # Temperature = 2.5
        k_Boltzmann = 1.0
        grand_potential = hubbard_model_gc.get_grand_potential(
            T=Temperature, k_B=k_Boltzmann
        )
        # print(f"Grand Potential (Omega) at T={Temperature}: {grand_potential:.4f}")
        print(f"Grand Potential (Omega) at T={Temperature}: {grand_potential}")

        grand_potential_T0 = hubbard_model_gc.get_grand_potential(T=0)
        print(f"Grand Potential (Omega) at T=0: {grand_potential_T0:.4f}")

        avg_N_T1 = hubbard_model_gc.get_average_total_particle_number(T=Temperature)
        print(
            f"Average total particle number <N_tot> at T={Temperature}: {avg_N_T1:.4f}"
        )

        avg_N_T0 = hubbard_model_gc.get_average_total_particle_number(T=0)
        print(
            f"Average total particle number <N_tot> at T=0 (Ground State): {avg_N_T0:.4f}"
        )

    else:
        print("No basis states generated. Cannot proceed.")

    print("-" * 40)

    # # Example: Non-interacting case U=0, mu=0
    # # Expected ground state: E=0 for N=0 (empty state)
    # print("Example: U=0, mu=0 (should favor empty state at T=0 if t > 0)")
    # hubbard_non_interacting = HubbardEDGrandCanonical(
    #     L_sites, t_hopping, 0.0, 0.0, pbc=True
    # )
    # if len(hubbard_non_interacting.basis_states) > 0:
    #     evals_ni = hubbard_non_interacting.diagonalize()
    #     # print(f"Eigenvalues of H' (U=0, mu=0): {np.sort(evals_ni)}")
    #     omega_ni_T0 = hubbard_non_interacting.get_grand_potential(T=0)
    #     print(
    #         f"Grand Potential at T=0 (U=0, mu=0): {omega_ni_T0:.4f}"
    #     )  # Should be 0 if empty state is GS
    #     avg_N_ni_T0 = hubbard_non_interacting.get_average_total_particle_number(T=0)
    #     print(f"Average N_tot at T=0 (U=0, mu=0): {avg_N_ni_T0:.4f}")

    # # Example: mu chosen to favor single particle
    # # Single particle energies for L=3, t=1, PBC: -2*t*cos(0) = -2; -2*t*cos(2pi/3) = 1; -2*t*cos(4pi/3) = 1
    # # So ground state of H_0 for N=1 is -2.
    # # If mu is, say, -1, then E' for N=1 GS is -2 - (-1)*1 = -1.
    # # E' for N=0 GS is 0 - (-1)*0 = 0. So N=1 GS is favored.
    # print("-" * 40)
    # print("Example: U=0, mu=-1.0 (should favor N=1 state at T=0)")
    # hubbard_N1_favored = HubbardEDGrandCanonical(
    #     L_sites, t_hopping, 0.0, -1.0, pbc=True
    # )
    # if len(hubbard_N1_favored.basis_states) > 0:
    #     omega_N1_T0 = hubbard_N1_favored.get_grand_potential(T=0)
    #     print(f"Grand Potential at T=0 (U=0, mu=-1.0): {omega_N1_T0:.4f}")
    #     avg_N_N1_T0 = hubbard_N1_favored.get_average_total_particle_number(T=0)
    #     print(
    #         f"Average N_tot at T=0 (U=0, mu=-1.0): {avg_N_N1_T0:.4f}"
    #     )  # Expect close to 1 (or 2 due to spin degeneracy)
    #     # Actually, 1.0 because specific N_up, N_down sectors are implicitly summed.
    #     # The N_tot operator does not distinguish spin, so it's sum of all particles.
