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
Z = \int \mathcal{D}\bar{\psi}\psi \exp\left(-S\right),
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
g(\tau, k) = \left<\mathcal{T} \psi(k, \tau) \bar{\psi}(k, 0) \right>_0= \frac{e^{-\epsilon_k \tau}}{1+e^{-\epsilon_k \beta}}\theta(\tau)+\xi \frac{e^{-\epsilon_k (\beta+\tau)}}{1+e^{-\epsilon_k \beta}}\theta(-\tau)
```

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

![Sign rule for the Wick contractions.](../assets/diagrams/grassman/green0.svg#green0)

![Diagrammatic expansion of the Green's function.](../assets/diagrams/grassman/green.svg#green)

The sign of a Green's function diagram is given by ``(-1)^{n_v} \xi^{n_F}``, where
1. ``n_v`` is the number of interactions.
2. ``n_F`` is the number of the fermionic loops.

## Feynman Rules for the Self-energy

![Diagrammatic expansion of the self-energy.](../assets/diagrams/grassman/sigma.svg#sigma)

## Feynman Rules for the 3-point Vertex Function

## Feynman Rules for the 4-point Vertex Function

## Feynman Rules for the Polarization
