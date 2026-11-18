# CommunityMatrices.jl

Interaction matrices of ecological communities, and what they say about stability.

"Interaction matrix" names several different objects: Lotka–Volterra coefficients, the Jacobian
at an equilibrium, the per-capita community matrix, feeding fluxes and sign patterns. Software
written for one silently accepts another. Here each is its own type, conversions between them
are explicit functions that ask for what they need, and every matrix carries a metadata record
of where it came from.

- **Typed matrices:** [`AlphaMatrix`](@ref), [`Jacobian`](@ref), [`CommunityMatrix`](@ref),
  [`FluxMatrix`](@ref), [`SignedDigraph`](@ref) and [`ProbabilisticNetwork`](@ref), with
  [`MatrixMetadata`](@ref).
- **Lotka–Volterra assembly:** equilibria and feasibility, subcommunities, invasion rates,
  steady states reached from a species pool, and trajectories under a forcing such as an
  antibiotic pulse.
- **Structural stability:** the size of the feasibility domain and the angle of the growth
  rates from its centre.
- **Sign patterns:** quasi-sign-stability and qualitative press responses.
- **Time series:** local Jacobians by the S-map, and temporal invariability.
- **Empirical data:** readers and bundled data for the Kinlock (2019) and Adler et al. (2018)
  plant interaction matrices and the Stein et al. (2013) gut-microbiota model.

The resilience metrics come from
[ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/), which defines them on
operators of every kind (equilibria, projection matrices, periodic and Markov operators) and is
re-exported here. Its documentation also lists the resilience tutorials of the other ecorecipes
packages.

```julia
using CommunityMatrices

a = AlphaMatrix([-1.0 -0.4; -0.3 -1.0], [1.0, 0.8])   # Lotka–Volterra: dN/dt = N (r + A N)
isfeasible(a), equilibrium(a)
feasibility_domain(a)                                   # (Ω, Ω^(1/S), se)
J = jacobian(a)                                         # a typed Jacobian at the equilibrium
resilience(J), reactivity(J), max_amplification(J)
quasi_sign_stability(SignedDigraph(J))
```

```@docs
CommunityMatrices
```
