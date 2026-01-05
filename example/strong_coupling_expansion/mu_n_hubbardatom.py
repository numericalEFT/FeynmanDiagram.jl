import numpy as np
import matplotlib.pyplot as plt


class HubbardAtom:
    def __init__(self, beta=10.0, U=2.0, h=0.0):
        """
        Initialize the Single-Site Hubbard Model (Hubbard Atom).

        Args:
            beta (float): Inverse temperature (1/kT)
            U (float): On-site interaction strength
            h (float): Magnetic field (Zeeman splitting term)
        """
        self.beta = beta
        self.U = U
        self.h = h

    def get_n_from_mu(self, mu):
        """
        Calculate density n given chemical potential mu.
        Uses the exact partition function sum.
        """
        # Exponentials for the Boltzmann weights
        # Using np.exp with standard overflow protection logic usually not strictly
        # needed for simple 0D unless beta*mu is massive, but good for clarity.

        # Terms: Empty (1), Singly Occupied (exp1, exp2), Doubly Occupied (exp_double)
        # Factor out exp(beta*mu) for numerical stability if needed,
        # but here we stick to the raw definition for readability.

        term_up = np.exp(self.beta * (mu + self.h))
        term_down = np.exp(self.beta * (mu - self.h))
        term_double = np.exp(self.beta * (2 * mu - self.U))

        Z = 1.0 + term_up + term_down + term_double

        # n = <N> / Z
        # Occupation numbers: 0, 1, 1, 2
        numerator = (0 * 1.0) + (1 * term_up) + (1 * term_down) + (2 * term_double)

        return numerator / Z

    def get_mu_from_n(self, target_n):
        """
        Calculate chemical potential mu given density n.
        Uses a numerically stable quadratic solver.
        """
        # 1. 边界情况处理
        if target_n <= 1e-14:
            return -np.inf
        if target_n >= 2.0 - 1e-14:
            return np.inf

        # 2. 计算系数: Ax^2 + Bx + C = 0, x = exp(beta * mu)
        # A = (2-n) * exp(-beta * U)
        # B = (1-n) * 2 * cosh(beta * h)
        # C = -n

        # 为了防止 A 下溢 (Underflow) 导致除以零，这里做个保护
        # 如果 beta*U 太大，exp(-beta*U) 可能在 float64 下变成 0
        try:
            factor_double = np.exp(-self.beta * self.U)
        except FloatingPointError:
            factor_double = 0.0

        factor_single = 2.0 * np.cosh(self.beta * self.h)

        a = (2.0 - target_n) * factor_double
        b = (1.0 - target_n) * factor_single
        c = -target_n

        # 3. 数值稳定的求根公式
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            # 理论上不可能发生，除非数值误差极大
            raise ValueError(f"Discriminant < 0 for n={target_n}")

        sqrt_disc = np.sqrt(discriminant)

        # 核心修复：避免两个相近的大数相减
        # 如果 b > 0 (对应 n < 1)，标准公式 -b + sqrt(...) 会发生消减。
        # 此时应使用变形公式: x = -2c / (b + sqrt(...))
        if b > 0:
            numerator = -2 * c
            denominator = b + sqrt_disc
            x = numerator / denominator
        else:
            # 如果 b <= 0 (对应 n >= 1)，标准公式无消减风险，或者两数同号相加
            # x = (-b + sqrt(...)) / 2a
            numerator = -b + sqrt_disc
            denominator = 2 * a

            # 特殊情况：如果 A 极其接近 0 (U 非常大)，退化为线性方程 Bx + C = 0
            if abs(denominator) < 1e-100:
                # x = -C / B
                x = -c / b
            else:
                x = numerator / denominator

        if x <= 0:
            # 如果仍然出错，说明参数极端到了浮点数无法表示的地步
            raise ValueError(
                f"Solver failed: x={x}, likely due to extreme beta/U values."
            )

        return np.log(x) / self.beta


# --- Example Usage ---
if __name__ == "__main__":
    # Parameters
    beta = 5.0
    U = 4.0
    # h = 0.5  # Non-zero magnetic field breaks spin symmetry
    h = 0.0

    atom = HubbardAtom(beta=beta, U=U, h=h)

    print(f"--- Hubbard Atom (beta={beta}, U={U}, h={h}) ---")

    # 1. Forward: mu -> n
    test_mu = 2.0  # Half-filling should be at mu = U/2 = 2.0
    n_calc = atom.get_n_from_mu(test_mu)
    print(f"Given mu = {test_mu}, calculated n = {n_calc:.6f}")

    # 2. Inverse: n -> mu
    test_n = 0.875
    mu_calc = atom.get_mu_from_n(test_n)
    print(f"Given n  = {test_n}, calculated mu = {mu_calc:.6f}")
    # Check consistency
    n_check = atom.get_n_from_mu(mu_calc)
    print(f"Consistency check: n(mu={mu_calc}) = {n_check}")

    atom1 = HubbardAtom(beta=beta / 2, U=U, h=h)
    mu_calc = atom1.get_mu_from_n(test_n)
    print(f"Given n  = {test_n}, calculated mu = {mu_calc:.6f}")
    n_check = atom1.get_n_from_mu(mu_calc)
    print(f"Consistency check: n(mu={mu_calc}) = {n_check}")

    atom1 = HubbardAtom(beta=beta * 2, U=U, h=h)
    mu_calc = atom1.get_mu_from_n(test_n)
    print(f"Given n  = {test_n}, calculated mu = {mu_calc:.6f}")
    n_check = atom1.get_n_from_mu(mu_calc)
    print(f"Consistency check: n(mu={mu_calc}) = {n_check}")

    atom1 = HubbardAtom(beta=200.0, U=U, h=h)
    mu_calc = atom1.get_mu_from_n(test_n)
    print(f"Given n  = {test_n}, calculated mu = {mu_calc:.6f}")
    n_check = atom1.get_n_from_mu(mu_calc)
    print(f"Consistency check: n(mu={mu_calc}) = {n_check}")

    # --- Plotting ---
    mu_vals = np.linspace(-2, 6, 200)
    n_vals = [atom.get_n_from_mu(mu) for mu in mu_vals]

    plt.figure(figsize=(8, 6))
    plt.plot(
        mu_vals, n_vals, label=f"Density n vs $\mu$\n($\\beta={beta}, U={U}, h={h}$)"
    )

    # Mark Half-filling
    mu_half = atom.get_mu_from_n(1.0)
    plt.scatter(
        [mu_half],
        [1.0],
        color="red",
        zorder=5,
        label=f"Half-filling ($\mu={mu_half:.2f}$)",
    )

    plt.xlabel(r"Chemical Potential $\mu$")
    plt.ylabel(r"Density $n$")
    plt.title("Hubbard Atom: Density vs Chemical Potential")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
