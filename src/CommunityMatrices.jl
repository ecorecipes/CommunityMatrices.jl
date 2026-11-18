"""
    CommunityMatrices

Interaction matrices of ecological communities, and what can be said about stability from them.

- **Typed matrices.** "Interaction matrix" names several objects: Lotka–Volterra coefficients
  ([`AlphaMatrix`](@ref)), the Jacobian at an equilibrium ([`Jacobian`](@ref)), per-capita
  community matrices ([`CommunityMatrix`](@ref)), feeding fluxes ([`FluxMatrix`](@ref)), sign
  patterns ([`SignedDigraph`](@ref)) and link probabilities ([`ProbabilisticNetwork`](@ref)).
  Each is its own type, conversions between them ask for what they need, and every matrix
  carries a [`MatrixMetadata`](@ref) record of where it came from.
- **Lotka–Volterra assembly:** equilibria, feasibility, subcommunities, invasion rates, steady
  states reached from a species pool, and forced trajectories.
- **Structural stability:** the feasibility domain and the angle from its centre (Saavedra et
  al. 2017).
- **Sign patterns:** quasi-sign-stability and qualitative press responses.
- **Time series:** local Jacobians by the S-map, and temporal invariability.
- **Empirical data:** readers for the Kinlock (2019) and Adler et al. (2018) plant interaction
  matrices and the Stein et al. (2013) gut-microbiota model, with the data.

The resilience metrics themselves (`resilience`, `reactivity`, `max_amplification`, the median
metrics, `stochastic_invariability`, the sensitivities) come from
[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl) and are re-exported; they
accept a [`Jacobian`](@ref) directly, through its `ContinuousJacobian`.
"""
module CommunityMatrices

using CSV: CSV
using DataFrames: DataFrame, groupby, rename!
using LinearAlgebra: I, Diagonal, diag, diagind, dot, lu, mul!, norm
using Statistics: mean, std
using OrdinaryDiffEqTsit5: Tsit5
using Random: Random, AbstractRNG
using SciMLBase: SciMLBase, ODEProblem, DiscreteCallback
using ResilienceCore: ResilienceCore, ContinuousTime, EquilibriumJacobian, ContinuousJacobian, resilience, isstable,
                      reactivity, amplification_envelope, propagator, max_amplification, pulse_covariance,
                      median_reactivity, median_displacement, median_max_amplification, median_resilience,
                      sensitivity_matrix, press_sensitivity, total_sensitivity, stochastic_invariability,
                      applicable_metrics, linearise
import ResilienceCore: jacobian, SignPattern, quasi_sign_stability, press_sign_determinacy

# re-exported from ResilienceCore
export ContinuousJacobian, resilience, isstable, reactivity, amplification_envelope, max_amplification,
       pulse_covariance, median_reactivity, median_displacement, median_max_amplification, median_resilience,
       sensitivity_matrix, press_sensitivity, total_sensitivity, stochastic_invariability, applicable_metrics,
       linearise, jacobian, SignPattern, quasi_sign_stability, press_sign_determinacy
# this package
export MatrixMetadata, InteractionMatrix, matrix, matrix_metadata, taxa, nspecies, AlphaMatrix, glv_rhs!, equilibrium,
       isfeasible, Jacobian, CommunityMatrix, per_capita, FluxMatrix, SignedDigraph, ProbabilisticNetwork,
       jacobian_metrics, subcommunity, invasion_rates, steady_state, glv_trajectory
export feasibility_domain, structural_vector, feasibility_angle
export smap_jacobians, invariability
export read_kinlock, alpha_from_rii, read_adler, read_stein

include("matrices.jl")
include("structural.jl")
include("timeseries.jl")
include("readers.jl")

end
