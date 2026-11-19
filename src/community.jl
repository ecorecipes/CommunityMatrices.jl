# A Lotka–Volterra community that implements ResilienceCore's community interface, so that the
# simulated perturbation experiments (mortality presses and tolerances, removals) and
# `stability_profile` run on it.

"""
    GLVCommunity(a::AlphaMatrix, mortality; N0 = fill(1e-3, nspecies(a)), steady = (;))

A Lotka–Volterra community `dNᵢ/dt = Nᵢ(rᵢ + Σⱼ Aᵢⱼ Nⱼ)` assembled from `N0` with
[`steady_state`](@ref), for ResilienceCore's simulation experiments and `stability_profile`.

The net growth rate does not say how much of it is mortality, so the caller supplies the
baseline per-capita mortality `mortality` (`dᵢ ≥ 0`, with `rᵢ = bᵢ − dᵢ`); there is no
default. Mortality multiplied by `fᵢ` (a press, as in `with_mortality`) changes the growth rate
to `rᵢ − (fᵢ − 1) dᵢ`. A species with `dᵢ = 0` is unaffected by a press on its own mortality.

After each perturbation the community returns to steady state with `steady_state` (keyword
arguments `steady`, e.g. `(tmax = 1e5,)`). Species absent before a perturbation are not
reintroduced. If an equilibrium has negative abundances, those species are removed and the
rest re-solved.
"""
struct GLVCommunity
    alpha::AlphaMatrix
    mortality::Vector{Float64}
    factors::Vector{Float64}
    N::Vector{Float64}
    steady::NamedTuple
end

function GLVCommunity(a::AlphaMatrix, mortality::AbstractVector; N0::AbstractVector = fill(1e-3, nspecies(a)),
                      steady::NamedTuple = (;))
    n = nspecies(a)
    length(mortality) == n || throw(ArgumentError("mortality has $(length(mortality)) entries for $n species"))
    all(≥(0), mortality) || throw(ArgumentError("mortality must be non-negative"))
    length(N0) == n || throw(ArgumentError("N0 has $(length(N0)) entries for $n species"))
    f = ones(n)
    return GLVCommunity(a, collect(Float64, mortality), f, settle(a, N0, steady), steady)
end

pressed(c::GLVCommunity, f) = AlphaMatrix(c.alpha.A, c.alpha.r .- (f .- 1) .* c.mortality; meta = c.alpha.meta)

function settle(a::AlphaMatrix, N0::AbstractVector, steady::NamedTuple)
    x = collect(Float64, N0)
    for _ in 0:length(x)
        N = steady_state(a, x; steady...).N
        any(<(0), N) || return N
        x = max.(N, 0.0)
    end
    return max.(steady_state(a, x; steady...).N, 0.0)
end

steadykw(c::GLVCommunity, steady) = steady === nothing ? c.steady : merge(c.steady, steady)

ResilienceCore.state(c::GLVCommunity) = c.N
ResilienceCore.present_species(c::GLVCommunity) = findall(>(0), c.N)

function ResilienceCore.with_mortality(c::GLVCommunity, factors; newton::Bool = true, steady = nothing)
    f = collect(Float64, factors)
    N = settle(pressed(c, f), c.N, steadykw(c, steady))
    return GLVCommunity(c.alpha, c.mortality, f, N, c.steady)
end

function ResilienceCore.without_species(c::GLVCommunity, i; steady = nothing)
    x = copy(c.N)
    x[i] = 0.0
    N = settle(pressed(c, c.factors), x, steadykw(c, steady))
    return GLVCommunity(c.alpha, c.mortality, c.factors, N, c.steady)
end

function ResilienceCore.jacobian(c::GLVCommunity)
    idx = present_species(c)
    return matrix(jacobian(pressed(c, c.factors); N = c.N))[idx, idx]
end
