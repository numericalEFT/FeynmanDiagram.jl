import numpy as np
from math import log, exp

############################################################
# Utility helpers
############################################################


def bit_count(x: int) -> int:
    # Python 3.10+: return x.bit_count()
    # Keep generic for portability:
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


def site_index(x, y, Lx, Ly):
    """Map 2D lattice coordinate -> linear site index [0 .. Lx*Ly-1]."""
    return y * Lx + x  # row-major


def neighbors_2d(Lx, Ly, pbc_x=True, pbc_y=True):
    """
    Precompute nearest-neighbor pairs for hopping on a 2D square lattice.
    Returns a list of lists: nbrs[i] = [j1, j2, ...] where j are NN of i.
    We'll include +x and +y hops, and later add both directions explicitly.
    """
    N = Lx * Ly
    nbrs = [[] for _ in range(N)]
    for y in range(Ly):
        for x in range(Lx):
            i = site_index(x, y, Lx, Ly)

            # +x neighbor
            x_right = x + 1
            if x_right < Lx:
                j = site_index(x_right, y, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)
            elif pbc_x and Lx > 1:
                j = site_index(0, y, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)

            # +y neighbor
            y_up = y + 1
            if y_up < Ly:
                j = site_index(x, y_up, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)
            elif pbc_y and Ly > 1:
                j = site_index(x, 0, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)

    # Remove accidental duplicates (can happen in 1xL or Lx1 cases)
    for i in range(N):
        nbrs[i] = sorted(list(set(nbrs[i])))
    return nbrs


def apply_hop(config: int, i: int, j: int):
    """
    Apply c_j^\dagger c_i to a given spin config bitstring.
    i = origin site, j = destination site.
    Return (new_config, fermion_sign) if hop is allowed, else (None, 0).

    Conventions:
    - If there's an electron at i and NOT at j, we move it.
    - Fermionic sign is (-1)^(#occupied sites strictly between i and j in linear ordering).
      We assume canonical ordering of sites 0..N-1 for the creation/annihilation.
    """
    occ_i = (config >> i) & 1
    occ_j = (config >> j) & 1
    if occ_i == 0 or occ_j == 1:
        return (None, 0.0)  # illegal hop

    # sign from annihilating at i: count how many occupied < i
    mask_before_i = (1 << i) - 1
    sign = (-1) ** (bit_count(config & mask_before_i))

    # remove particle at i
    config_after_ann = config ^ (1 << i)

    # sign from creating at j: count how many occupied < j in the new config
    mask_before_j = (1 << j) - 1
    sign *= (-1) ** (bit_count(config_after_ann & mask_before_j))

    # place particle at j
    new_config = config_after_ann | (1 << j)

    return (new_config, float(sign))


############################################################
# Main ED class
############################################################


class HubbardED2DGrandCanonical:
    """
    Exact Diagonalization for the 2D Hubbard model on an Lx x Ly cluster
    in the grand canonical ensemble.
    """

    def __init__(self, Lx, Ly, t, U, mu, pbc_x=True, pbc_y=True):
        """
        Args:
            Lx, Ly : linear dimensions
            t      : hopping amplitude
            U      : onsite repulsion
            mu     : chemical potential
            pbc_x, pbc_y : periodic boundary conditions in x / y
        """
        self.Lx = Lx
        self.Ly = Ly
        self.N_sites = Lx * Ly
        self.t = t
        self.U = U
        self.mu = mu
        self.pbc_x = pbc_x
        self.pbc_y = pbc_y

        # basis: list of (conf_up, conf_dn)
        self.basis_states = []
        self.basis_index = {}
        self.N_updown_list = []  # store (N_up, N_dn) for each basis idx

        self._build_basis()

        dim = len(self.basis_states)
        self.Hprime = np.zeros((dim, dim), dtype=np.float64)

        # Precompute nearest neighbors on lattice
        self.nbrs = neighbors_2d(Lx, Ly, pbc_x=pbc_x, pbc_y=pbc_y)

        self._build_hamiltonian()

        # Will be filled by diagonalize()
        self.evals = None
        self.evecs = None

    ########################################################
    # Basis
    ########################################################
    def _build_basis(self):
        """
        Full grand-canonical Fock space:
        all 2^(N_sites) configs for spin-up times all 2^(N_sites) configs for spin-down.
        """
        idx = 0
        for up in range(1 << self.N_sites):
            Nu = bit_count(up)
            for dn in range(1 << self.N_sites):
                Nd = bit_count(dn)

                st = (up, dn)
                self.basis_states.append(st)
                self.basis_index[st] = idx
                self.N_updown_list.append((Nu, Nd))
                idx += 1

    ########################################################
    # Hamiltonian
    ########################################################
    def _build_hamiltonian(self):
        """
        Build dense H' = H - mu N.
        H = kinetic + U * sum_i n_up n_dn
        - mu N just shifts diagonal.
        """
        dim = len(self.basis_states)
        H = self.Hprime  # alias

        for k, (up_conf, dn_conf) in enumerate(self.basis_states):
            Nu, Nd = self.N_updown_list[k]
            Ntot = Nu + Nd

            # On-site U term + chemical potential shift
            diag_E = 0.0
            for site in range(self.N_sites):
                occ_up = (up_conf >> site) & 1
                occ_dn = (dn_conf >> site) & 1
                if occ_up and occ_dn:
                    diag_E += self.U
            diag_E -= self.mu * Ntot
            H[k, k] = diag_E

            # Kinetic term: for each spin independently, hop along NN bonds
            # We'll do directed hops i->j. That automatically fills both H[kf,ki]
            # and (later when we loop over other k) H[ki,kf] so we don't need
            # separate Hermitian symmetrization.
            for i_site in range(self.N_sites):
                for j_site in self.nbrs[i_site]:
                    if j_site == i_site:
                        continue

                    # spin-up hopping
                    new_up, sign_up = apply_hop(up_conf, i_site, j_site)
                    if new_up is not None:
                        st_f = (new_up, dn_conf)
                        kf = self.basis_index[st_f]
                        H[kf, k] += -self.t * sign_up

                    # spin-down hopping
                    new_dn, sign_dn = apply_hop(dn_conf, i_site, j_site)
                    if new_dn is not None:
                        st_f = (up_conf, new_dn)
                        kf = self.basis_index[st_f]
                        H[kf, k] += -self.t * sign_dn

        # Done. Note: This double-counts bonds i<->j (because neighbors_2d
        # symmetrizes), but that's FINE because hopping i->j and j->i are distinct
        # matrix elements. We're not adding both directions for the *same* operator;
        # we're literally applying the operator c_j^dag c_i for all ordered pairs.
        # The result is still Hermitian.

    ########################################################
    # Diagonalization
    ########################################################
    def diagonalize(self):
        """
        Diagonalize H' and cache results.
        """
        if self.evals is None or self.evecs is None:
            self.evals, self.evecs = np.linalg.eigh(self.Hprime)
        return self.evals, self.evecs

    ########################################################
    # Thermodynamics helpers
    ########################################################
    def _thermal_weights(self, T, kB=1.0):
        """
        Return normalized Boltzmann weights w_m \propto exp(-beta E_m')
        and the partition function Z_G.
        Also returns beta and E0 shift to stabilize numerics.
        """
        evals, _ = self.diagonalize()
        if T == 0:
            # Pure ground state projector
            gs_idx = np.argmin(evals)
            w = np.zeros_like(evals)
            w[gs_idx] = 1.0
            ZG = 1.0
            return w, ZG, np.inf, 0.0  # beta=inf conceptually, shift=0

        beta = 1.0 / (kB * T)
        Emin = np.min(evals)
        shifted = evals - Emin
        boltz_unnorm = np.exp(-beta * shifted)
        ZG = np.sum(boltz_unnorm)
        w = boltz_unnorm / ZG
        return w, ZG, beta, Emin

    def grand_potential(self, T, kB=1.0):
        """
        Omega = -kB T ln Z_G, with H' eigenvalues.
        At T=0, Omega = min(E').
        """
        evals, _ = self.diagonalize()
        if len(evals) == 0:
            return np.nan
        if T == 0:
            return float(np.min(evals))

        w, ZG, beta, Emin = self._thermal_weights(T, kB)
        # Z_G = sum exp(-beta (E - Emin)) = ZG
        # ln Z_G = ln(ZG) - beta * (-Emin)?? -> careful:
        # We defined: shifted = E - Emin
        # boltz = exp(-beta * shifted)
        # ZG = sum boltz = sum exp(-beta(E - Emin))
        # actual Z_G = sum exp(-beta E) = exp(-beta*Emin) * ZG
        # so ln Z_G = -beta*Emin + ln(ZG)
        lnZ_actual = (-beta * Emin) + log(ZG)
        return -(1.0 / beta) * lnZ_actual

    def average_total_particle_number(self, T, kB=1.0):
        """
        <N_tot> = sum_m w_m <m| N |m>.
        N = sum_i (n_up + n_dn).
        """
        evals, evecs = self.diagonalize()
        dim = len(self.basis_states)
        w, ZG, beta, Emin = self._thermal_weights(T, kB)

        # Precompute N_tot for each basis state (diagonal operator)
        N_basis = np.empty(dim, dtype=np.float64)
        for idx, (Nu, Nd) in enumerate(self.N_updown_list):
            N_basis[idx] = Nu + Nd

        # Expectation value in each eigenstate m:
        N_expect_m = np.einsum(
            "im,i,im->m", evecs, N_basis, evecs
        )  # sum_i |psi_i|^2 N_i
        return float(np.dot(w, N_expect_m))

    def filling(self, T, kB=1.0):
        """n = <N_tot>/N_sites."""
        return self.average_total_particle_number(T, kB) / self.N_sites

    def local_density(self, T, site, kB=1.0):
        """
        <n_i> = <n_i_up + n_i_dn>.
        site: int in [0, N_sites-1]
        """
        evals, evecs = self.diagonalize()
        w, ZG, beta, Emin = self._thermal_weights(T, kB)

        dim = len(self.basis_states)

        # n_i is diagonal in Fock basis
        n_i_basis = np.empty(dim, dtype=np.float64)
        for idx, (up_conf, dn_conf) in enumerate(self.basis_states):
            occ = ((up_conf >> site) & 1) + ((dn_conf >> site) & 1)
            n_i_basis[idx] = occ

        n_i_expect_m = np.einsum("im,i,im->m", evecs, n_i_basis, evecs)
        return float(np.dot(w, n_i_expect_m))

    def average_density(self, T, kB=1.0):
        """
        (1/N_sites) sum_i <n_i>.
        For translationally invariant clusters, this ~ filling(),
        but for small clusters with broken symmetry it's nice to measure directly.
        """
        return np.mean([self.local_density(T, i, kB=kB) for i in range(self.N_sites)])

    def double_occupancy(self, T, kB=1.0):
        """
        D = (1/N_sites) sum_i <n_{i up} n_{i dn}>.
        Diagonal in Fock basis.
        """
        evals, evecs = self.diagonalize()
        w, ZG, beta, Emin = self._thermal_weights(T, kB)

        dim = len(self.basis_states)
        d_i_total_expect_m = np.zeros_like(evals)

        for site in range(self.N_sites):
            d_i_basis = np.empty(dim, dtype=np.float64)
            for idx, (up_conf, dn_conf) in enumerate(self.basis_states):
                occ_up = (up_conf >> site) & 1
                occ_dn = (dn_conf >> site) & 1
                d_i_basis[idx] = occ_up * occ_dn

            d_i_expect_m = np.einsum("im,i,im->m", evecs, d_i_basis, evecs)
            d_i_total_expect_m += d_i_expect_m

        D_expect = np.dot(w, d_i_total_expect_m) / self.N_sites
        return float(D_expect)

    def kinetic_energy(self, T, kB=1.0):
        """
        K = < -t sum_<ij>,sigma (c_j^\dagger c_i + h.c.) >
        We'll compute <psi_m | H_kin | psi_m> and thermal average.
        We already have the *full* H' matrix and we know:
            H' = H_kin + H_U - mu N.
        So:
            H_kin = H' - H_U + mu N.
        We'll build diagonal operator (H_U - mu N) and subtract.
        """
        evals, evecs = self.diagonalize()
        w, ZG, beta, Emin = self._thermal_weights(T, kB)

        dim = len(self.basis_states)

        # Build diag array for (H_U - mu N) in basis, because those are diagonal in Fock basis.
        diag_HUint_minus_muN = np.zeros(dim, dtype=np.float64)
        for idx, (up_conf, dn_conf) in enumerate(self.basis_states):
            Nu, Nd = self.N_updown_list[idx]
            Ntot = Nu + Nd
            HU = 0.0
            for site in range(self.N_sites):
                if ((up_conf >> site) & 1) and ((dn_conf >> site) & 1):
                    HU += self.U
            diag_HUint_minus_muN[idx] = HU - self.mu * Ntot

        # <psi_m|H'|psi_m> is just evals[m]
        # so <psi_m|H_kin|psi_m> = evals[m] - <psi_m|diag_HUint_minus_muN|psi_m>
        diag_expect_m = np.einsum("im,i,im->m", evecs, diag_HUint_minus_muN, evecs)
        K_expect_m = evals - diag_expect_m  # vector over m
        return float(np.dot(w, K_expect_m))


############################################################
# Example usage / quick test
############################################################
if __name__ == "__main__":
    Lx, Ly = 3, 2  # 2x2 plaquette (16 spin states per spin -> 256 total states)
    # Lx, Ly = 2, 2  # 2x2 plaquette (16 spin states per spin -> 256 total states)
    # Lx, Ly = 4, 4  # 2x2 plaquette (16 spin states per spin -> 256 total states)
    t = 1.0

    # for U in [10.0, 20.0]:
    for U in [5.0]:
        # for U in [5.0, 10.0, 15.0, 20.0, 30.0]:
        # U = 5.0
        # U = 10.0
        mu = U / 2
        # mu = U
        T = 1.0
        # for mu in [U / 2 - 3.0, U / 2 - 2.0, U / 2 - 1.0, U / 2]:
        model = HubbardED2DGrandCanonical(Lx, Ly, t, U, mu, pbc_x=True, pbc_y=True)
        print(
            "Hubbard model (Lx={}, Ly={}, t={}, U={}, mu={}):".format(Lx, Ly, t, U, mu)
        )
        print("dim(H) =", len(model.basis_states))
        print("Omega(T) =", model.grand_potential(T))
        print("<N> =", model.average_total_particle_number(T))
        print("filling n =", model.filling(T))
        print("Double occupancy D =", model.double_occupancy(T))
        print("Avg density (check) =", model.average_density(T))
        print(
            "Local n(i) per site =",
            [model.local_density(T, i) for i in range(model.N_sites)],
        )
        print("Kinetic energy =", model.kinetic_energy(T))
        print("========================================\n")

        # for dmu in [-0.5, -0.1, 0.0, 0.1, 0.5]:
        # for dmu in [-0.1, 0.0, 0.1]:
        for dmu in [0.0]:
            # mu = U / 2  # often near half-filling in bipartite lattices
            # mu = 1.0

            # free_model = HubbardED2DGrandCanonical(Lx, Ly, t, 0.0, mu, pbc_x=True, pbc_y=True)
            free_model = HubbardED2DGrandCanonical(
                Lx,
                Ly,
                0.0,
                U,
                mu + dmu,
                pbc_x=True,
                pbc_y=True,
            )

            print(f"Non-hopping model (U={U}, T={T}, mu={mu} + {dmu}):")
            print("Omega0(T)        =", free_model.grand_potential(T))
            print("<N>(T)           =", free_model.average_total_particle_number(T))
            print("filling n(T)     =", free_model.filling(T))
            print("D0(T)            =", free_model.double_occupancy(T))
            print(
                "Local n(i) per site =",
                [free_model.local_density(T, i) for i in range(free_model.N_sites)],
                "\n",
            )
            print("========================================\n")

            print(
                "(Omega(T) - Omega0(T)) / V =",
                (model.grand_potential(T) - free_model.grand_potential(T)) / (Lx * Ly),
            )
