# Metrics from ResilienceCore

The resilience and stability metrics are defined in
[ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/), on operators whose
type says what they are, and re-exported here: `resilience`, `isstable`, `reactivity`,
`amplification_envelope`, `propagator`, `max_amplification`, `pulse_covariance`,
`median_reactivity`, `median_displacement`, `median_max_amplification`, `median_resilience`,
`sensitivity_matrix`, `press_sensitivity`, `total_sensitivity`, `stochastic_invariability` and
`applicable_metrics`.

Each accepts a [`Jacobian`](@ref) directly, by way of its `ContinuousJacobian(j)`; the median
metrics then use pulses proportional to the Jacobian's equilibrium. A plain matrix `J` is
analysed as `ContinuousJacobian(J)`, or `ContinuousJacobian(J, equilibrium)` for the median
metrics. [`jacobian_metrics`](@ref) collects the Jacobian-based metrics in one record.
