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
jacobian_metrics
```

## Assembly and trajectories

```@docs
subcommunity
invasion_rates
steady_state
glv_trajectory
```
