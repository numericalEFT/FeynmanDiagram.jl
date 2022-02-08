# Interaction/Scattering-Amplitude Convention

In general, the interaction (or scattering amplitude) between two spin-$1/2$ particles is controlled by only two parameters. Here we briefly review some of the common conventions of the parameters in literature, and show how to map between different conventions.

We assume a particle from the left collides with a particle from the right. Due to the interaction between the two particles, one of them scatters to the left, and the other to the right. 

We will use the index $\alpha,\beta,\gamma,\delta$ to label the left incoming, the left outgoing, the right incoming and the right outgoing particles. We assume four legs associate with four momentum-frequency vectors $k_1, k_2, k_3, k_4$.

There are two possible scatterings: a direct scattering where the two particles don't permutate after the collision, and an exchange scattering where the two particles permutate. Therefore, the total scattering amplitude has two contributions,
```math
v_{\alpha\beta\gamma\delta}(12;34) = v^d_{\alpha\beta\gamma\delta}(12;34)+v^e_{\alpha\beta\gamma\delta}(12;34)
```
where the exchange contribution is different from the direct counterpart up to an exchange of two outgoing legs and an overall particle statistic sign $\xi =\pm 1$,
```math
v^e_{\alpha\beta\gamma\delta}(12;34) = \xi v^d_{\alpha\delta\gamma\beta}(14;32)
```
where we have abbreviated the momentum-frequency vectors as indices.

## 1. Spin symmetric and asymmetric convention:

A interaction is split into a term that is independent of the spin (e.g., Coulomb interaction between electrons) and a term with spin dependence (e.g., ferromagnetic/antiferromagnetic interaction).
This convention is commonly used in textbooks of Fermi liquid theory.

```math
v^d_{\alpha\beta\gamma\delta}(12;34) \equiv v^{s}_{12;34}\delta_{\alpha\beta}\delta_{\gamma\delta} + v^a_{12;34} \vec{\sigma}_{\alpha\beta}\cdot \vec{\sigma}_{\gamma\delta}
```
where the superscript $d$ means the scattering is direct.

If two particles permute after the collision, then one needs to exchange $\beta \leftrightarrow \gamma$ and other internal variables,
```math
v^e_{\alpha\beta\gamma\delta}(12;34) = \xi v^s_{14;32}\delta_{\alpha\delta}\delta_{\gamma\beta} + \xi v^a_{14;32} \vec{\sigma}_{\alpha\delta}\cdot \vec{\sigma}_{\gamma\beta}
```

```math
v^e_{\alpha\beta\gamma\delta}(12;34) = \xi\frac{v^s_{14;32}+3v^a_{14;32}}{2}\delta_{\alpha\beta}\delta_{\gamma\delta} + \xi \frac{v^s_{14;32}-v^a_{14;32}}{2} \vec{\sigma}_{\alpha\beta}\cdot \vec{\sigma}_{\gamma\delta}
```

To derive the above equations, we use the identity $\vec{\sigma}_{\alpha\beta}\cdot \vec{\sigma}_{\gamma\delta}=2 \delta_{\alpha \delta} \delta_{\beta \gamma}-\delta_{\alpha \beta} \delta_{\gamma \delta}$, which gives the following equations,
```math
\delta_{\alpha\delta}\delta_{\gamma\beta}=\frac{1}{2}\delta_{\alpha\beta}\delta_{\gamma\delta}+\frac{1}{2}\vec{\sigma}_{\alpha\beta}\cdot \vec{\sigma}_{\gamma\delta},
```
and,
```math
\vec{\sigma}_{\alpha\delta}\cdot \vec{\sigma}_{\gamma\beta}=\frac{3}{2}\delta_{\alpha\beta}\delta_{\gamma\delta}-\frac{1}{2}\vec{\sigma}_{\alpha\beta}\cdot \vec{\sigma}_{\gamma\delta}
```

## 2. Spin $\uparrow\uparrow$ and $\uparrow\downarrow$ convention:

An alternative parameterization of the interaction matrix is by specifiying the spin components $v_{\uparrow\uparrow} \equiv v_{\uparrow\uparrow\uparrow\uparrow}$ and $v_{\uparrow\downarrow} \equiv v_{\uparrow\uparrow\downarrow\downarrow}$. They can be derived from $v_s$ and $v_a$ by the simple relation,
```math
v^d_{\uparrow\uparrow}(12;34) = v^s_{12;34}+v^a_{12;34}
```
```math
v^d_{\uparrow\downarrow}(12;34) = v^s_{12;34}-v^a_{12;34}
```

If two particles permute after the collision, the exchange interaction $v^e$,
```math
v^e_{\uparrow\uparrow}(12;34) = \xi\frac{v^s_{14;32}+3v^a_{14;32}}{2}+\xi\frac{v^s_{14;32}-v^s_{14;32}}{2} = \xi v^s_{14;32}+ \xi v^a_{14;32}
```
```math
v^e_{\uparrow\downarrow}(12;34) = \xi\frac{v^s_{14;32}+3v^a_{14;32}}{2}-\xi\frac{v^s_{14;32}-v^s_{14;32}}{2} = 2 \xi v^a_{14;32}
```