# Interaction matrices and Lotka–Volterra

Typed matrices make explicit which object a number belongs to, and conversions between them
ask for what they need.

```@docs
MatrixMetadata
InteractionMatrix
matrix
matrix_metadata
taxa
nspecies
AlphaMatrix
glv_rhs!
equilibrium
isfeasible
Jacobian
Jacobian(::ContinuousJacobian)
CommunityMatrix
per_capita
FluxMatrix
jacobian(::AlphaMatrix)
jacobian(::FluxMatrix)
SignedDigraph
ProbabilisticNetwork
```

## Assembly and trajectories

```@docs
subcommunity
invasion_rates
steady_state
glv_trajectory
```

## Simulated perturbations

A [`GLVCommunity`](@ref) implements ResilienceCore's community interface, so its mortality
presses, tolerances and removals, and `stability_profile` (all 27 metrics of
Domínguez-García et al. 2019), run on it. The baseline mortality must be supplied.

```@docs
GLVCommunity
```
