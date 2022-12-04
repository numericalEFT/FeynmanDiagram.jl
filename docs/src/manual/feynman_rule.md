# Feynman Rules

In general, we follow the convention in the textbook "Quantum Many-particle Systems" by J. Negele and H. Orland, Page 95,

## Fourier Transform

```math
G(\tau) = \frac{1}{\beta} \sum_n G(i\omega_n) \text{e}^{-i\omega_n \tau}
```

```math
G(i\omega_n) = \int_0^\beta G(\tau) \text{e}^{i\omega_n \tau} d\tau
```

where the Matsubara-frequency ``\omega_n=2\pi n/\beta`` for boson and ``\omega_n = 2\pi (n+1)/\beta`` for fermion.

## Action and Partition Sum

The partition sum associates with a generic action,

```math
Z = \int \mathcal{D}\bar{\psi}\mathcal{D}\psi \exp\left(-S\right),
```

where the action takes a generic form,

```math
S = \bar{\psi}_1\left(\frac{\partial}{\partial \tau} +\epsilon_k \right)\psi_1 + V_{1234}\bar{\psi}_1\bar{\psi}_2\psi_3\psi_4,
```

where ``\bar{\psi}`` and ``\psi`` are Grassman fields.

In the Matsubara-frequency domain, the action is,

```math
S = \bar{\psi}_1\left(-i\omega_n +\epsilon_k \right)\psi_1 + V_{1234}\bar{\psi}_1\bar{\psi}_2\psi_3\psi_4,
```

## Bare Propagator

- Imaginary time

```math
g(\tau, k) = \left<\mathcal{T} \psi(k, \tau) \bar{\psi}(k, 0) \right>_0= \frac{e^{-\epsilon_k \tau}}{1+e^{-\epsilon_k \beta}}\theta(\tau)+\xi \frac{e^{-\epsilon_k (\beta+\tau)}}{1+e^{-\epsilon_k \beta}}\theta(-\tau),
```

where ``\xi`` is ``+1`` for boson and ``-1`` for fermion.

- Matusbara frequency

```math
g(i\omega_n, k) = -\frac{1}{i\omega_n-\epsilon_k},
```

Then the action takes a simple form,

```math
S = \bar{\psi}_1g_{12}^{-1}\psi_2 + V_{1234}\bar{\psi}_1\bar{\psi}_2\psi_3\psi_4,
```

## Dressed Propagator and Self-energy

The dressed propagator is given by,

```math
G(\tau, k) = \left<\mathcal{T} \psi(k, \tau) \bar{\psi}(k, 0) \right>,
```

and we define the self-energy ``\Sigma`` as the one-particle irreducible vertex function,

```math
G^{-1} = g^{-1} + \Sigma,
```

so that

```math
G = g - g\Sigma g + g\Sigma g \Sigma g - ...
```

## Perturbative Expansion of the Green's Function

![Sign rule for the Wick contractions.](../../assets/diagrams/green0.svg#green0)

![Diagrammatic expansion of the Green's function.](../../assets/diagrams/green.svg#green)

The sign of a Green's function diagram is given by ``(-1)^{n_v} \xi^{n_F}``, where

1. ``n_v`` is the number of interactions.
2. ``n_F`` is the number of the fermionic loops.

## Feynman Rules for the Self-energy

From the Green's function diagrams, one can derive the __negative__ self-energy diagram,

![Diagrammatic expansion of the self-energy.](../../assets/diagrams/sigma.svg#sigma)

```math
\begin{aligned}
-\Sigma = & (-1) \xi V_{34} g_{44}+(-1) V_{34} g_{34} \\
+&(-1)^2 \xi V_{34} V_{56} g_{46} g_{64} g_{43}+(-1)^2 V_{34} V_{56} g_{35} g_{54} g_{42}+\cdots
\end{aligned}
```

The sign of a __negative__ self-energy ``-\Sigma`` diagram is given by ``(-1)^{n_v} \xi^{n_F}``, where

1. ``n_v`` is the number of interactions.
2. ``n_F`` is the number of the fermionic loops.

## Feynman Rules for the 3-point Vertex Function

The self-energy is related to the 3-point vertex function through an equation,

```math
-\left(\Sigma_{3, x} -\Sigma^{Hartree}_{3, x}\right) = G_{3,y} \cdot \left(-V_{3, 4}\right) \cdot \Gamma^3_{4,y,x},
```

where the indices $x, y$ could be different from diagrams to diagrams, and $\Gamma_3$ is the inproper three-vertex function. Eliminate the additional sign, one derives,

```math
\Sigma_{3, x} -\Sigma^{Hartree}_{3, x} = G_{3,y} \cdot V_{3, 4} \cdot \Gamma^3_{4,y,x},
```

![Diagrammatic expansion of the 3-point vertex function.](../../assets/diagrams/gamma3.svg#gamma3)

The diagram weights are given by,

```math
\begin{aligned}
\Gamma^{(3)}= & 1 + (-1) \xi V_{56} g_{46} g_{64} + (-1) V_{56} g_{54} g_{46}\\
+&(-1)^2 \xi^2 V_{56} V_{78} g_{46} g_{64} g_{58} g_{85}+(-1)^2\xi V_{56} V_{78} g_{74} g_{46}+\cdots
\end{aligned}
```

The sign of ``\Gamma^{(3)}`` diagram is given by ``(-1)^{n_v} \xi^{n_F}``.

## Feynman Rules for the 4-point Vertex Function

The 4-point vertex function is related to the 3-point vertex function through an equation,

```math
\Gamma^{(3)}_{4,y,x} = \xi \cdot G_{4,s} \cdot G_{t, 4} \cdot \Gamma^{(4)}_{s, t, y, x},
```

where the indices $x, y, s, t$ could be different from diagrams to diagrams.

![Diagrammatic expansion of the 4-point vertex function.](../../assets/diagrams/gamma4.svg#gamma4)

The diagram weights are given by,

```math
\begin{aligned}
\Gamma^{(4)}= & (-1) V_{56}^{\text{direct}} + (-1)\xi V_{56}^{exchange}\\
+&(-1)^2 \xi V_{56} V_{78} g_{58} g_{85}+(-1)^2 V_{56} V_{78}+\cdots,
\end{aligned}
```

where we used the identity ``\xi^2 = 1``.

The sign of ``\Gamma^{(4)}`` diagram is given by ``(-1)^{n_v} \xi^{n_F}`` multiplied with a sign from the permutation of the external legs.

## Feynman Rules for the Susceptibility

The susceptibility can be derived from ``\Gamma^{(4)}``.

```math
\chi_{1,2} \equiv \left<\mathcal{T} n_1 n_2\right>_{\text{connected}} = \xi G_{1,2} G_{2, 1} + \xi G_{1,s} G_{t, 1} \Gamma^{(4)}_{s, t, y, x} G_{2,y} G_{x, 2}
```

![Diagrammatic expansion of the susceptibility.](../../assets/diagrams/susceptibility.svg#susceptibility)

We define the polarization ``P`` as the one-interaction irreducible (or proper) vertex function,

```math
\chi^{-1} = P^{-1} + V,
```
