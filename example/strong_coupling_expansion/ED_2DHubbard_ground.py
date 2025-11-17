import numpy as np
from itertools import combinations
from math import comb
from scipy.sparse import dok_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


def site_index(x, y, Lx, Ly):
    return y * Lx + x  # row-major


def make_neighbors(Lx, Ly, pbc_x=True, pbc_y=True):
    """
    neighbors[i] = list of NN sites j for site i
    on an Lx x Ly square lattice with optional periodic BCs.
    """
    N = Lx * Ly
    nbrs = [[] for _ in range(N)]
    for y in range(Ly):
        for x in range(Lx):
            i = site_index(x, y, Lx, Ly)

            # +x neighbor
            x_r = x + 1
            if x_r < Lx:
                j = site_index(x_r, y, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)
            elif pbc_x and Lx > 1:
                j = site_index(0, y, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)

            # +y neighbor
            y_u = y + 1
            if y_u < Ly:
                j = site_index(x, y_u, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)
            elif pbc_y and Ly > 1:
                j = site_index(x, 0, Lx, Ly)
                nbrs[i].append(j)
                nbrs[j].append(i)

    # unique + sorted
    for i in range(N):
        nbrs[i] = sorted(list(set(nbrs[i])))
    return nbrs


def bit_count(x: int) -> int:
    # portable version of x.bit_count()
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


def apply_hop(config: int, i: int, j: int):
    """
    Apply c_j^† c_i to a spin sector bitstring 'config'.
    Returns (new_config, sign) or (None, 0) if forbidden.

    Sign rule:
      c_j^† c_i |state>
      = (-1)^(#occupied sites < i) (-1)^(#occupied sites < j AFTER removing at i)
        |state with particle moved i->j|
    """
    occ_i = (config >> i) & 1
    occ_j = (config >> j) & 1
    if occ_i == 0 or occ_j == 1:
        return (None, 0.0)

    # fermion sign from annihilating at i
    mask_before_i = (1 << i) - 1
    sign = (-1) ** (bit_count(config & mask_before_i))

    # remove from i
    cfg2 = config ^ (1 << i)

    # fermion sign from creating at j in cfg2
    mask_before_j = (1 << j) - 1
    sign *= (-1) ** (bit_count(cfg2 & mask_before_j))

    # put at j
    cfg_new = cfg2 | (1 << j)
    return (cfg_new, float(sign))


class HubbardEDSectorSparse:
    """
    Exact diagonalization *in a fixed particle-number sector*
    for the 2D Hubbard model on an Lx x Ly cluster, using sparse storage.
    We build H' = H - mu N. (grand-canonical shifted operator)
    We then get the lowest k eigenvalues/eigenvectors with eigsh.
    Observables are evaluated from those eigenvectors (T=0 approx).
    """

    def __init__(self, Lx, Ly, t, U, mu, Nup, Ndn, pbc_x=True, pbc_y=True):
        """
        Args:
          Lx, Ly : lattice size
          t, U, mu : Hubbard params
          Nup, Ndn : fixed number of up/down electrons in this sector
          pbc_x, pbc_y : boundary conditions
        """
        self.Lx = Lx
        self.Ly = Ly
        self.Ns = Lx * Ly
        self.t = t
        self.U = U
        self.mu = mu
        self.Nup = Nup
        self.Ndn = Ndn
        self.pbc_x = pbc_x
        self.pbc_y = pbc_y

        # Precompute nearest-neighbor list
        self.nbrs = make_neighbors(Lx, Ly, pbc_x, pbc_y)

        # Build sector basis
        # basis_up: all bitstrings with exactly Nup bits set
        # basis_dn: all bitstrings with exactly Ndn bits set
        self.basis_up = list(self._fixed_pop_bitstrings(self.Ns, Nup))
        self.basis_dn = list(self._fixed_pop_bitstrings(self.Ns, Ndn))

        # Cartesian product for full basis in this sector
        # index_map[(up_bits, dn_bits)] = basis index
        self.basis_states = []
        self.index_map = {}
        idx = 0
        for up_cfg in self.basis_up:
            for dn_cfg in self.basis_dn:
                self.basis_states.append((up_cfg, dn_cfg))
                self.index_map[(up_cfg, dn_cfg)] = idx
                idx += 1

        self.dim = len(self.basis_states)
        print(f"[info] Sector dimension dim = {self.dim}")

        # We'll assemble the Hamiltonian in DOK then convert to CSR
        self.Hprime = dok_matrix((self.dim, self.dim), dtype=np.float64)
        self._build_sparse_hprime()
        self.Hprime = self.Hprime.tocsr()

        # place-holders for eigensystem (lowest k)
        self.evals = None
        self.evecs = None

    def _fixed_pop_bitstrings(self, L, N):
        """
        Generate all bitstrings of length L with exactly N ones.
        Yield as integers.
        """
        for occ_sites in combinations(range(L), N):
            cfg = 0
            for s in occ_sites:
                cfg |= 1 << s
            yield cfg

    def _build_sparse_hprime(self):
        """
        Construct H' = H - mu N in this fixed-(Nup,Ndn) sector.

        Because Nup, Ndn are fixed, -mu N is just a constant shift
        -mu*(Nup+Ndn) on the diagonal — easy.

        H has two parts:
          - Kinetic: -t sum_{<ij>,σ} (c†_{jσ} c_{iσ} + h.c.)
          - Interaction: U sum_i n_{i↑} n_{i↓}  (diagonal in Fock basis)
        """
        const_mu_shift = -self.mu * (self.Nup + self.Ndn)

        for bra_idx, (up_cfg, dn_cfg) in enumerate(self.basis_states):
            # --- diagonal contribution: U * sum_i n_up(i) n_dn(i) + const_mu_shift
            double_occ = 0.0
            for site in range(self.Ns):
                if ((up_cfg >> site) & 1) and ((dn_cfg >> site) & 1):
                    double_occ += 1.0
            diag_E = self.U * double_occ + const_mu_shift
            self.Hprime[bra_idx, bra_idx] += diag_E

            # --- off-diagonal kinetic terms
            # spin-up hopping
            for i_site in range(self.Ns):
                if ((up_cfg >> i_site) & 1) == 0:
                    continue  # no up electron at i_site
                for j_site in self.nbrs[i_site]:
                    new_up, sign_up = apply_hop(up_cfg, i_site, j_site)
                    if new_up is not None:
                        ket_state = (new_up, dn_cfg)
                        ket_idx = self.index_map.get(ket_state, None)
                        if ket_idx is not None:
                            self.Hprime[ket_idx, bra_idx] += -self.t * sign_up

            # spin-down hopping
            for i_site in range(self.Ns):
                if ((dn_cfg >> i_site) & 1) == 0:
                    continue
                for j_site in self.nbrs[i_site]:
                    new_dn, sign_dn = apply_hop(dn_cfg, i_site, j_site)
                    if new_dn is not None:
                        ket_state = (up_cfg, new_dn)
                        ket_idx = self.index_map.get(ket_state, None)
                        if ket_idx is not None:
                            self.Hprime[ket_idx, bra_idx] += -self.t * sign_dn

    def solve_low_energy(self, k=10):
        """
        Compute the lowest k eigenvalues/eigenvectors of H' using ARPACK.
        Store them in self.evals (shape k) and self.evecs (dim x k).
        """
        # 'SA' = smallest algebraic eigenvalues
        evals, evecs = eigsh(self.Hprime, k=k, which="SA")
        # Sort ascending just in case
        order = np.argsort(evals)
        self.evals = evals[order]
        self.evecs = evecs[:, order]
        return self.evals, self.evecs

    ########################################################
    # Ground-state observables (T = 0)
    ########################################################
    def ground_state_energy(self):
        """
        Return E0' (the lowest eigenvalue of H').
        """
        if self.evals is None:
            raise RuntimeError("Call solve_low_energy() first.")
        return float(self.evals[0])

    def gs_vector(self):
        """
        Ground-state eigenvector in the chosen sector.
        """
        if self.evecs is None:
            raise RuntimeError("Call solve_low_energy() first.")
        return self.evecs[:, 0]

    def gs_double_occupancy(self):
        """
        D = (1/Ns) * sum_i <n_i_up n_i_dn> in GS.
        Diagonal operator in the Fock basis.
        """
        psi0 = self.gs_vector()
        prob = psi0.conjugate() * psi0  # |psi_i|^2 in basis
        D_sum = 0.0
        for basis_idx, (up_cfg, dn_cfg) in enumerate(self.basis_states):
            docc_local = 0
            for site in range(self.Ns):
                occu = (up_cfg >> site) & 1
                occd = (dn_cfg >> site) & 1
                docc_local += occu & occd
            D_sum += prob[basis_idx] * docc_local
        return float(D_sum / self.Ns)

    def gs_density_per_site(self):
        """
        n_i = <n_i_up + n_i_dn> for each site in GS.
        Returns array of length Ns.
        """
        psi0 = self.gs_vector()
        prob = psi0.conjugate() * psi0

        n_site = np.zeros(self.Ns, dtype=np.float64)
        for basis_idx, (up_cfg, dn_cfg) in enumerate(self.basis_states):
            for site in range(self.Ns):
                occ = ((up_cfg >> site) & 1) + ((dn_cfg >> site) & 1)
                n_site[site] += prob[basis_idx] * occ
        return n_site  # length Ns

    def gs_total_particle_number(self):
        """
        <N> in GS. In a fixed-(Nup,Ndn) sector this is just Nup+Ndn,
        but we compute it via expectation for completeness.
        """
        n_site = self.gs_density_per_site()
        return float(np.sum(n_site))

    def gs_kinetic_energy(self):
        """
        K = <psi0| H_kin |psi0>, where
         H' = H_kin + U*sum_i n_i_up n_i_dn - mu N.
        So:
         H_kin = H' - (U * double_occ_op) + mu N.
        We'll evaluate all three expectation values term by term.
        """
        psi0 = self.gs_vector()
        prob = psi0.conjugate() * psi0

        # <psi0|H'|psi0> = E0'
        E0p = self.ground_state_energy()

        # <double_occ_op>
        D_abs = 0.0
        for basis_idx, (up_cfg, dn_cfg) in enumerate(self.basis_states):
            docc_local = 0
            for site in range(self.Ns):
                occu = (up_cfg >> site) & 1
                occd = (dn_cfg >> site) & 1
                docc_local += occu & occd
            D_abs += prob[basis_idx] * docc_local
        # <N>
        N_tot = 0.0
        for basis_idx, (up_cfg, dn_cfg) in enumerate(self.basis_states):
            Nu = bit_count(up_cfg)
            Nd = bit_count(dn_cfg)
            N_tot += prob[basis_idx] * (Nu + Nd)

        K_expect = E0p - (self.U * D_abs) + self.mu * N_tot
        return float(K_expect)


############################################################
# Example usage for 4x4 half filling
############################################################
if __name__ == "__main__":
    Lx, Ly = 4, 4  # 4x4 cluster → 16 sites
    t = 1.0
    U = 10.0
    mu = U / 2  # chemical potential near half filling for bipartite lattice
    Nup = 8
    Ndn = 8

    model = HubbardEDSectorSparse(Lx, Ly, t, U, mu, Nup, Ndn, pbc_x=True, pbc_y=True)

    # Ask ARPACK for, say, the 20 lowest eigenstates in this sector.
    # NOTE: this can still be very expensive in RAM/CPU for 4x4 at half filling.
    evals, evecs = model.solve_low_energy(k=20)

    print("Lowest energies (shifted H'):", evals[:5])
    print("Ground state grand potential Omega0 = E0' =", model.ground_state_energy())
    print("<N> in GS =", model.gs_total_particle_number())
    print("Double occupancy D =", model.gs_double_occupancy())
    dens = model.gs_density_per_site()
    print("Average density per site =", np.mean(dens))
    print("Kinetic energy K =", model.gs_kinetic_energy())
