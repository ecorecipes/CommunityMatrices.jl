# CommunityMatrices.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/CommunityMatrices.jl/dev/)

Interaction matrices of ecological communities, and what they say about stability: typed
Lotka–Volterra, Jacobian, community, flux and sign matrices with provenance metadata;
Lotka–Volterra assembly; structural stability; sign-pattern analyses; S-map Jacobians from time
series; and readers, with the data, for published plant and gut-microbiota interaction
matrices.

The resilience metrics come from
[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl) and are re-exported.

```julia
using CommunityMatrices

a = AlphaMatrix([-1.0 -0.4; -0.3 -1.0], [1.0, 0.8])   # Lotka–Volterra: dN/dt = N (r + A N)
isfeasible(a), equilibrium(a)
feasibility_domain(a)                                   # (Ω, Ω^(1/S), se)
J = jacobian(a)                                         # a typed Jacobian at the equilibrium
resilience(J), reactivity(J), max_amplification(J)
quasi_sign_stability(SignedDigraph(J))
```

Tutorials (Quarto, in `vignettes/`):
1. structural stability and sign patterns of empirical plant communities (Kinlock 2019; Adler
   et al. 2018);
2. a gut microbiota under antibiotics (Stein et al. 2013): equilibria, invasion, the
   clindamycin pulse, and what S-maps can and cannot recover from the data.

`data/` holds the three datasets, each with a README giving its source, licence (CC0, or CC BY
4.0 and MIT for Stein) and checksums.

The package was split from EcologicalResilience.jl (commit `e84f286`), whose reproduction of
Domínguez-García, Dakos & Kéfi (2019) continues as DominguezGarcia2019.jl.

## Status

Not registered; the ecorecipes dependencies are fetched by URL through `[sources]`
(Julia ≥ 1.11). Tests: `julia --project -e 'using Pkg; Pkg.test()'`.
