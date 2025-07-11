import sympy as sp
import numpy as np
from scipy.integrate import dblquad

# Note: This script now requires the SciPy library.
# You can install it by running: pip install scipy

# Original Mathematica Expression, for reference:
#
# g[\[Tau]_] := -(Exp[\[Mu] \[Tau]] +
#      Exp[\[Beta] \[Mu]] Exp[(\[Mu] -
#            U) \[Tau]]) HeavisideTheta[\[Tau]] + (Exp[\[Mu] (\[Tau] + \[Beta])] +
#      Exp[\[Beta] (2 \[Mu] - U)] Exp[(\[Mu] -
#            U) \[Tau]]) HeavisideTheta[-\[Tau]]
# Z = (1 + 2 Exp[\[Beta] \[Mu]] + Exp[\[Beta] (2 \[Mu] - U)])^3
# F3 = Integrate[
#    g[\[Tau]] g[t - \[Tau]] g[\[Beta] - t], {\[Tau], 0, \[Beta]}, {t,
#     0, \[Beta]}]*2/Z/\[Beta]


def create_symbolic_components():
    """
    Sets up the symbolic components of the F3 expression using sympy.
    Instead of creating Integral objects, this function returns the raw
    symbolic integrands and other components.
    """
    # Define symbolic variables.
    beta, mu, U, tau, t = sp.symbols("beta mu U tau t", real=True)
    beta = sp.Symbol("beta", real=True, positive=True)
    x = sp.Symbol("x", real=True)

    # Expression for g(x) when x > 0
    g_pos = -(sp.exp(mu * x) + sp.exp(beta * mu) * sp.exp((mu - U) * x))

    # Expression for g(x) when x < 0
    g_neg = sp.exp(mu * (x + beta)) + sp.exp(beta * (2 * mu - U)) * sp.exp((mu - U) * x)

    # --- Construct the integrand components ---
    g1 = g_pos.subs(x, tau)  # g(τ), where τ >= 0
    g3 = g_pos.subs(x, beta - t)  # g(β - t), where t <= β

    # Integrand for the region where t > τ, so g(t - τ) uses g_pos
    integrand_pos = g1 * g_pos.subs(x, t - tau) * g3

    # Integrand for the region where t < τ, so g(t - τ) uses g_neg
    integrand_neg = g1 * g_neg.subs(x, t - tau) * g3

    # Normalization factor Z
    Z = (1 + 2 * sp.exp(beta * mu) + sp.exp(beta * (2 * mu - U))) ** 3

    # Return all the symbolic components and the variables
    components = {
        "integrand_pos": integrand_pos,
        "integrand_neg": integrand_neg,
        "Z": Z,
    }
    variables = {"beta": beta, "mu": mu, "U": U, "t": t, "tau": tau}

    return components, variables


def main():
    """
    Main function to perform numerical evaluation using SciPy.
    """
    print("Setting up symbolic expressions with SymPy...")
    components, variables = create_symbolic_components()

    # Unpack variables for easier access
    beta, mu, U = variables["beta"], variables["mu"], variables["U"]
    t, tau = variables["t"], variables["tau"]

    # --- Numerical Calculation using SciPy ---
    print("\nStarting numerical evaluation with SciPy...")

    # Define the numerical values for the parameters.
    num_vals = {beta: 0.2, mu: 1.0, U: 4.0}
    print(
        f"\nUsing numerical values: {{'beta': {num_vals[beta]}, 'mu': {num_vals[mu]}, 'U': {num_vals[U]}}}"
    )

    # Substitute numerical values into the symbolic components
    integrand_pos_num = components["integrand_pos"].subs(num_vals)
    integrand_neg_num = components["integrand_neg"].subs(num_vals)
    Z_num = components["Z"].subs(num_vals)
    beta_val = num_vals[beta]

    # Convert the symbolic integrands into fast numerical functions
    # The order [t, tau] is important for dblquad, which expects func(y, x)
    # where y is the inner variable of integration (t) and x is the outer (tau).
    f_pos = sp.lambdify([t, tau], integrand_pos_num, "numpy")
    f_neg = sp.lambdify([t, tau], integrand_neg_num, "numpy")

    try:
        # Perform the numerical double integration using SciPy's dblquad
        # dblquad(func, x_min, x_max, y_min_func(x), y_max_func(x))
        # Here, x is tau, and y is t.

        print("Evaluating integral for t > tau...")
        # For this part, t ranges from tau to beta
        integral_val_pos, err_pos = dblquad(
            f_pos, 0, beta_val, lambda tau_val: tau_val, lambda tau_val: beta_val
        )

        print("Evaluating integral for t < tau...")
        # For this part, t ranges from 0 to tau
        integral_val_neg, err_neg = dblquad(
            f_neg, 0, beta_val, lambda tau_val: 0, lambda tau_val: tau_val
        )

        # --- Calculate and Display Final Result ---
        integral_total = integral_val_pos + integral_val_neg
        Z_val = Z_num.evalf()  # .evalf() is fine for this non-integral part

        # The final F3 expression
        F3_val = (2 / Z_val) * integral_total

        print("\n" + "=" * 40)
        print("      NUMERICAL RESULTS")
        print("=" * 40)
        print(f"Value of Z: {Z_val:.6f}")
        print(f"Total integral value (numerical): {integral_total:.10f}")
        print(f"Error estimate: +/- {err_pos + err_neg:.2e}")
        print("-" * 40)
        print(f"Final Numerical Result for F3: {F3_val:.10f}")
        print("=" * 40)

    except Exception as e:
        print(f"\nAn error occurred during the numerical evaluation: {e}")
        print("Please ensure you have SciPy installed ('pip install scipy').")


if __name__ == "__main__":
    main()
