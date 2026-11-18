# A gut microbiota under antibiotics: stability from a fitted model and from data
Simon Frost

- [Overview](#overview)
- [The fitted interaction matrix](#the-fitted-interaction-matrix)
- [The community the model predicts](#the-community-the-model-predicts)
- [Stability of the predicted
  community](#stability-of-the-predicted-community)
- [The clindamycin pulse in the
  model](#the-clindamycin-pulse-in-the-model)
- [Model against observation](#model-against-observation)
- [Stability from the time series
  alone](#stability-from-the-time-series-alone)
- [What these data can and cannot
  validate](#what-these-data-can-and-cannot-validate)
- [Caveats](#caveats)
- [References](#references)

## Overview

Stein et al. (2013) fitted a generalised Lotka–Volterra (gLV) model to
the intestinal microbiota of mice. They used 11 genera from 9 mice: some
were given a single dose of clindamycin, and some were also given
*Clostridium difficile* spores. The fit includes each genus’s
susceptibility to the antibiotic. This makes it a rare empirical system
with both an interaction matrix and a recorded pulse perturbation.

This vignette uses the published fit to show the package’s
interaction-matrix layer:

- the fitted [`AlphaMatrix`](@ref) and why its full equilibrium is not
  feasible;
- assembling the community the model predicts (`steady_state`,
  `subcommunity`, `invasion_rates`);
- the stability metrics of that community (`jacobian_metrics`);
- the clindamycin pulse, simulated with `glv_trajectory` and compared
  with the mice;
- what can be estimated from the time series alone (`smap_jacobians`,
  `invariability`).

The conclusion is mostly negative. These data are too short, and the fit
too loose, to validate resilience metrics. The vignette shows why, and
what would be needed.

``` julia
using CommunityMatrices
using DataFrames, Statistics, LinearAlgebra, Random, Printf
using Plots
st = read_stein(joinpath(pkgdir(CommunityMatrices), "data", "stein2013"))
a, ε, obs = st.alpha, st.susceptibility, st.data
genera = taxa(a)
short(x) = replace(x, "undefined_genus_of_unclassified_" => "u.g. ", "undefined_genus_of_" => "u.g. ",
                   "unclassified_" => "uncl. ", "Clostridium_difficile" => "C. difficile", "_" => " ")
(genera = nspecies(a), samples = nrow(unique(obs[:, [:mouse, :time]])), mice = length(unique(obs.mouse)))
```

    (genera = 11, samples = 77, mice = 9)

The data are bundled in `data/stein2013/` with their provenance: the
published fit (Dataset S1 of the article, CC BY 4.0) and the time series
(from the cLV repository, MIT). `read_stein` returns the fit as an
`AlphaMatrix`, with the units and source recorded in its metadata:

``` julia
matrix_metadata(a).derivation, matrix_metadata(a).units
```

    ("Stein et al. 2013, Dataset S1 sheet MmuE (published regularised fit)", "per day")

## The fitted interaction matrix

In the model, dNᵢ/dt = Nᵢ (rᵢ + Σⱼ Aᵢⱼ Nⱼ + εᵢ u(t)), where u(t) is the
clindamycin signal. Abundances are in 10¹¹ 16S copies per cm³ and time
is in days.

``` julia
A = matrix(a)
off = A[.!I(nspecies(a))]
(self_regulated = all(<(0), diag(A)), positive = count(>(0), off), negative = count(<(0), off))
```

    (self_regulated = true, positive = 45, negative = 65)

``` julia
lim = maximum(abs, A)
heatmap(short.(genera), short.(genera), A; yflip = true, c = :RdBu, clims = (-lim, lim),
        xrotation = 45, size = (720, 600), colorbar_title = "Aᵢⱼ")
```

![The published interaction matrix: the effect of the column genus on
the per-capita growth of the row genus. Rows are sorted by clindamycin
susceptibility, most harmed at the
top.](02_gut_microbiota_files/figure-commonmark/cell-5-output-1.svg)

``` julia
DataFrame(genus = short.(genera), growth = round.(a.r; digits = 2), susceptibility = round.(ε; digits = 2))
```

<div><div style = "float: left;">11×3 DataFrame</div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | genus                   |  growth | susceptibility |
|----:|:------------------------|--------:|---------------:|
|     | String                  | Float64 |        Float64 |
|   1 | Barnesiella             |    0.37 |          -3.29 |
|   2 | u.g. Lachnospiraceae    |    0.31 |          -3.04 |
|   3 | uncl. Lachnospiraceae   |    0.36 |          -2.09 |
|   4 | Other                   |    0.54 |          -1.94 |
|   5 | Blautia                 |    0.71 |          -1.35 |
|   6 | u.g. Mollicutes         |    0.47 |           -1.1 |
|   7 | Akkermansia             |    0.23 |          -0.92 |
|   8 | Coprobacillus           |    0.83 |          -0.79 |
|   9 | C. difficile            |    0.39 |          -0.31 |
|  10 | Enterococcus            |    0.29 |           1.07 |
|  11 | u.g. Enterobacteriaceae |    0.32 |            3.7 |

```@raw html
</div>
```

- **Every genus is self-regulated** (Aᵢᵢ \< 0).
- **The interactions are mostly negative:** 65 of the 110 off-diagonal
  effects.
- **Clindamycin harms nine genera and helps two:** *Enterococcus* and
  the Enterobacteriaceae genus.

With every genus present, the equilibrium N\* = −A⁻¹r is not feasible:

``` julia
(feasible = isfeasible(a), negative_abundances = count(<(0), equilibrium(a)))
```

    (feasible = false, negative_abundances = 5)

## The community the model predicts

The model’s resting community comes from `steady_state`:

- it integrates from a starting state;
- removes the genera that fall below a threshold;
- solves for the survivors’ equilibrium;
- reintroduces any genus of the regional pool whose invasion rate is
  positive.

The start is the mean day-0 community of the three clindamycin-only
mice. By default, the pool is the genera present at the start. These
mice never received *C. difficile*.

``` julia
day0 = [mean(obs[(obs.group .== 2) .& (obs.time .== 0) .& (obs.genus .== g), :abundance]) for g in genera]
ss = steady_state(a, day0)
(survivors = short.(genera[ss.survivors]), stable = ss.stable, uninvadable = ss.uninvadable)
```

    (survivors = ["Other", "Blautia", "u.g. Mollicutes", "Akkermansia", "Coprobacillus", "u.g. Enterobacteriaceae"], stable = true, uninvadable = true)

The surviving community has 6 genera. It is stable, and no absent genus
of the pool can invade it. `invasion_rates` gives the per-capita growth
rate of each absent genus at low density:

``` julia
gr = invasion_rates(a, ss.N)
absent = setdiff(eachindex(genera), ss.survivors)
DataFrame(genus = short.(genera[absent]), invasion_rate = round.(gr[absent]; digits = 3),
          present_at_day_0 = day0[absent] .> 0)
```

<div><div style = "float: left;">5×3 DataFrame</div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | genus                 | invasion_rate | present_at_day_0 |
|----:|:----------------------|--------------:|-----------------:|
|     | String                |       Float64 |             Bool |
|   1 | Barnesiella           |        -0.181 |             true |
|   2 | u.g. Lachnospiraceae  |        -0.191 |             true |
|   3 | uncl. Lachnospiraceae |        -0.169 |             true |
|   4 | C. difficile          |         0.238 |            false |
|   5 | Enterococcus          |         -0.01 |             true |

```@raw html
</div>
```

- **The excluded genera** (*Barnesiella*, the two Lachnospiraceae genera
  and *Enterococcus*) have negative invasion rates.
- ***C. difficile* is the exception.** Its invasion rate is positive, so
  it would invade. With all 11 genera in the pool, it does:

``` julia
short.(genera[steady_state(a, day0; pool = eachindex(genera)).survivors])
```

    7-element Vector{String}:
     "Other"
     "Blautia"
     "u.g. Mollicutes"
     "Akkermansia"
     "Coprobacillus"
     "C. difficile"
     "u.g. Enterobacteriaceae"

This is a check on the fit, not only a result. The mice of population 1
were given spores without antibiotics, and *C. difficile* did not
establish in them. At the state the fitted model assembles, the model
shows no colonisation resistance. At the observed day-0 states of those
mice, the model gives a mixed answer:

``` julia
icd = findfirst(==("Clostridium_difficile"), genera)
[(mouse = m, invasion_rate = round(invasion_rates(a, [obs[(obs.mouse .== m) .& (obs.time .== 0) .& (obs.genus .== g),
                                                        :abundance][1] for g in genera])[icd]; digits = 3))
 for m in sort(unique(obs.mouse[obs.group .== 1]))]
```

    3-element Vector{@NamedTuple{mouse::Int64, invasion_rate::Float64}}:
     (mouse = 2, invasion_rate = 0.114)
     (mouse = 5, invasion_rate = -0.021)
     (mouse = 8, invasion_rate = 0.012)

- **Two of the three mice** start where *C. difficile* could invade
  (rates 0.114 and 0.012 per day).
- **Mouse 5** starts where it could not, but only just (−0.021).
- **Mouse 8** carried *C. difficile* at a low level (about 10⁻³)
  throughout. Mice 2 and 5 had none.

The fit is therefore at best weakly consistent with the colonisation
resistance these mice showed.

## Stability of the predicted community

The survivors’ Jacobian comes with its equilibrium, so every metric
applies, including the median metrics of Domínguez-García et al. (2019):

``` julia
sub = subcommunity(a, ss.survivors)
J = jacobian(sub)
m = jacobian_metrics(J)
(Rinf = m.Rinf, return_time = 1 / m.Rinf, R0 = m.R0, Amax = m.Amax, tmax = m.tmax,
 MRinf = m.MRinf, MR0 = m.MR0, MAmax = m.MAmax)
```

    (Rinf = 0.06305951806270148, return_time = 15.858034293977283, R0 = 0.1856312718485403, Amax = 1.4519462534759764, tmax = 4.1, MRinf = 0.09108035781409002, MR0 = -0.45335049726082705, MAmax = 1.0)

- **Slow return:** the asymptotic resilience is 0.063 per day, a return
  time of about 16 days.
- **Reactive:** worst-case perturbations grow at first (R0 \> 0) and are
  amplified up to 1.45× after about 4 days.
- **Median pulses do not amplify** (MR0 \< 0, MAmax = 1). Here a median
  pulse is one scaled by abundance.

## The clindamycin pulse in the model

The antibiotic is a pulse: u = 1 for one day, which adds εᵢ to each
genus’s per-capita growth rate. Here the pulse is applied to the model’s
own equilibrium.

``` julia
times = collect(0:0.25:120.0)
pulse(t) = t < 1 ? 1.0 : 0.0
traj = glv_trajectory(a, ss.N, times; forcing = pulse, susceptibility = ε, tstops = [1.0])
lg(x) = log10.(max.(x, 1e-12))
dist = [norm(lg(traj[k, ss.survivors]) .- lg(ss.N[ss.survivors])) for k in eachindex(times)]
p1 = plot(times, lg(traj[:, ss.survivors]); label = permutedims(short.(genera[ss.survivors])), xlabel = "day",
          ylabel = "log10 abundance", legend = :bottomright, legendfontsize = 6)
p2 = plot(times, dist; label = "", xlabel = "day", ylabel = "distance (log10 units)")
vline!(p2, [1 / m.Rinf]; ls = :dash, label = "1/Rinf")
plot(p1, p2; layout = (1, 2), size = (860, 340))
```

![A one-day clindamycin pulse applied to the model’s equilibrium. Left:
abundances of the six genera. Right: distance from the equilibrium in
log10 units, with the return time 1/Rinf
marked.](02_gut_microbiota_files/figure-commonmark/cell-13-output-1.svg)

``` julia
at(t) = dist[findfirst(==(t), times)]
(peak = maximum(dist), peak_day = times[argmax(dist)], day_23 = at(23.0), day_60 = at(60.0), day_120 = at(120.0))
```

    (peak = 2.7541835353664137, peak_day = 8.5, day_23 = 2.130166972847733, day_60 = 0.8851348149703921, day_120 = 0.04312846280235173)

- **The displacement keeps growing after the pulse ends.** It peaks at
  about 2.8 log10 units around day 8, the transient amplification that
  R0 \> 0 allows.
- **Recovery takes months.** The community is still 2.1 log10 units away
  on day 23, the last day observed in the mice, and returns only over
  the following months.

## Model against observation

For the three clindamycin-only mice, the model is started from each
mouse’s observed day-0 state and given the same one-day pulse. The
prediction is then compared with the observed abundances on days 2–23.

``` julia
function mouse_series(mse)
    s = obs[obs.mouse .== mse, :]
    ts = sort(unique(s.time))
    X = [s[(s.time .== t) .& (s.genus .== g), :abundance][1] for t in ts, g in genera]
    return ts, X
end
fl(x; floor = 1e-6) = log10.(max.(x, floor))
plt = plot(; xlabel = "day", ylabel = "distance from day 0", legend = :bottomright)
rows = []
for (k, mse) in enumerate(sort(unique(obs.mouse[obs.group .== 2])))
    ts, X = mouse_series(mse)
    P = glv_trajectory(a, X[1, :], ts; forcing = pulse, susceptibility = ε, tstops = [1.0])
    dO = [norm(fl(X[i, :]) .- fl(X[1, :])) for i in eachindex(ts)]
    dP = [norm(fl(P[i, :]) .- fl(X[1, :])) for i in eachindex(ts)]
    plot!(plt, ts, dO; color = k, marker = :circle, label = "mouse $mse")
    plot!(plt, ts, dP; color = k, ls = :dash, label = "")
    late = ts .>= 2
    o, p = vec(X[late, :]), vec(P[late, :])
    both = (o .> 0) .& (p .> 1e-12)
    push!(rows, (mouse = mse, cor_floor_1e4 = cor(fl(o; floor = 1e-4), fl(p; floor = 1e-4)),
                 cor_floor_1e6 = cor(fl(o), fl(p)), cor_floor_1e8 = cor(fl(o; floor = 1e-8), fl(p; floor = 1e-8)),
                 cor_present_only = cor(log10.(o[both]), log10.(p[both]))))
end
plt
```

![Distance from each mouse’s own day-0 state (log10 units, zeros floored
at 10⁻⁶) for the three clindamycin-only mice: observed (solid) and
predicted by the fitted model
(dashed).](02_gut_microbiota_files/figure-commonmark/cell-15-output-1.svg)

``` julia
DataFrame(rows)
```

<div><div style = "float: left;">3×5 DataFrame</div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | mouse | cor_floor_1e4 | cor_floor_1e6 | cor_floor_1e8 | cor_present_only |
|----:|------:|--------------:|--------------:|--------------:|-----------------:|
|     | Int64 |       Float64 |       Float64 |       Float64 |          Float64 |
|   1 |     1 |      0.153497 |      0.416137 |      0.632802 |       -0.0655451 |
|   2 |     4 |      0.188612 |      0.443774 |      0.653775 |        0.0325784 |
|   3 |     7 |      0.443673 |       0.58727 |      0.723572 |         0.389298 |

```@raw html
</div>
```

- **The observed shift is far larger than predicted.** The mice move 5–8
  log10 units from their day-0 state by day 2, while the model moves
  about 2.7. The observed shift is driven by genera that drop to zero.
- **None recovers by day 23.** The observed communities stay far from
  day 0. The model drifts further, towards its own equilibrium: the
  observed day-0 state is not the model’s equilibrium.
- **The correlation between predicted and observed log-abundances
  measures the zeros, not the fit.**
  - It rises from about 0.15–0.44 to 0.63–0.72 as the floor for zero
    abundances moves from 10⁻⁴ to 10⁻⁸.
  - Over the genera present in both observation and prediction, it is
    close to zero for two of the three mice.

## Stability from the time series alone

`smap_jacobians` estimates local Jacobians from a multivariate time
series (Deyle et al. 2016), and `invariability` measures temporal
stability as mean/sd. Each needs enough data.

**What the S-map needs.**

- **Equally spaced samples.** The mice were sampled on days 0, 2, 3, 4,
  5, 6, 7, 9, 12, 16 and 23.
- **More time points than genera plus two.** Each mouse has at most 11
  points, for 11 genera.

So the method cannot be applied to these series directly. To see what it
would need, apply it to series simulated from the fitted 6-genus
community. Each day the community follows the model, and is then
multiplied by log-normal process noise (5% on the log scale).

``` julia
function noisy(rng, sub, N0, T; σ = 0.05)
    X = zeros(T, length(N0))
    x = copy(N0)
    for t in 1:T
        X[t, :] = x
        x = glv_trajectory(sub, x, [0.0, 1.0])[end, :] .* exp.(σ .* randn(rng, length(x)))
    end
    return X
end
Xsim = noisy(Xoshiro(1), sub, ss.N[ss.survivors], 500)
Js = smap_jacobians(Xsim; theta = 1.0, times = 101:20:481)
long = [resilience(j) for j in Js]
# error of the whole matrix, against the model's Jacobian at the same state
jac_error = [opnorm(j.J - jacobian(sub; N = j.equilibrium).J) / opnorm(jacobian(sub; N = j.equilibrium).J) for j in Js]
short_windows = [resilience(smap_jacobians(Xsim[s:(s + 10), :]; theta = 0.0, times = [10])[1]) for s in 1:10:480]
(true_Rinf = m.Rinf, from_500_days = median(long), relative_error_of_J = median(jac_error),
 from_11_days_median = median(short_windows), from_11_days_range = extrema(short_windows))
```

    (true_Rinf = 0.06305951806270148, from_500_days = 0.06811233948466101, relative_error_of_J = 0.9268124981741833, from_11_days_median = -0.007026155025189483, from_11_days_range = (-0.5754211286062723, 0.41343949650662937))

``` julia
plot([fill(1, length(long)); fill(2, length(short_windows))], [long; short_windows]; seriestype = :scatter,
     xticks = ([1, 2], ["500 daily points", "11 daily points"]), xlims = (0.5, 2.5), label = "",
     ylabel = "estimated Rinf (per day)", alpha = 0.6, size = (520, 380))
hline!([m.Rinf]; ls = :dash, label = "true Rinf")
```

![Resilience estimated by the S-map from simulated series: 20 target
times in one 500-day series (θ = 1), and 48 separate 11-day windows (θ =
0). The dashed line is the true resilience of the model’s
equilibrium.](02_gut_microbiota_files/figure-commonmark/cell-18-output-1.svg)

- **From 500 daily points**, the S-map recovers the dominant rate of
  return to within about 10%.
- **The full Jacobian is recovered much less well.** The error in its
  norm is about as large as the matrix itself: small noise excites some
  directions too little for their entries to be estimated.
- **From 11 daily points**, the estimates scatter widely around zero:
  they cannot even tell whether the community is stable.

Invariability needs no model, so it can be computed on the mice
directly:

``` julia
invtab = combine(groupby(obs, [:mouse, :group]), [:time, :abundance] => ((t, x) -> begin
                  tot = [sum(x[t .== s]) for s in sort(unique(t))]
                  (points = length(tot), invariability = invariability(tot))
              end) => AsTable)
sort!(invtab, [:group, :mouse])
invtab.treatment = getindex.(Ref(Dict(1 => "spores only", 2 => "clindamycin", 3 => "clindamycin + spores")), invtab.group)
invtab
```

<div><div style = "float: left;">9×5 DataFrame</div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | mouse | group | points | invariability | treatment            |
|----:|------:|------:|-------:|--------------:|:---------------------|
|     | Int64 | Int64 |  Int64 |       Float64 | String               |
|   1 |     2 |     1 |      4 |       5.19644 | spores only          |
|   2 |     5 |     1 |      4 |       4.70595 | spores only          |
|   3 |     8 |     1 |      4 |       2.88673 | spores only          |
|   4 |     1 |     2 |     11 |      0.774032 | clindamycin          |
|   5 |     4 |     2 |     11 |       1.01702 | clindamycin          |
|   6 |     7 |     2 |     11 |       1.17662 | clindamycin          |
|   7 |     3 |     3 |     10 |      0.840661 | clindamycin + spores |
|   8 |     6 |     3 |     11 |       1.11215 | clindamycin + spores |
|   9 |     9 |     3 |     11 |      0.975551 | clindamycin + spores |

```@raw html
</div>
```

- **Untreated mice are steadier.** The total abundance of the mice that
  received no antibiotic varies less (invariability 2.9–5.2) than that
  of treated mice (0.8–1.2).
- **But this measures the perturbation, not an intrinsic property.** The
  treated mice were knocked far from their state, and the untreated ones
  have only four samples each.
- **In the simulation, process noise alone gives invariability near
  30.**

## What these data can and cannot validate

- **They can check the fit’s qualitative claims.** For example, the
  fitted model allows *C. difficile* to invade the community it
  assembles, which the mice of population 1 do not support.
- **They cannot validate return times.** The predicted return time
  (about 16 days, longer after the amplification) is as long as the
  whole series. The observed communities had not recovered when sampling
  stopped.
- **They cannot support data-driven Jacobians.** That would need long,
  evenly spaced series with many more samples than taxa. Long
  perturbation studies such as MDSINE2 (Bucci et al.
  2016) come closer.

## Caveats

- **Zeros.** About 9% of the abundances are exactly zero. Log-scale
  distances and correlations depend on how zeros are floored; this
  vignette uses 10⁻⁶ and shows the sensitivity.
- **One published fit.** The fit is regularised (Stein et al. 2013,
  Dataset S1), and its uncertainty is not used here. Refitting is
  outside the scope of this package.
- **The pulse shape.** The clindamycin signal is 1 on the day of the
  dose and 0 afterwards. A one-day pulse is the simplest reading of that
  signal.

## References

- Bucci V. et al. (2016) MDSINE: Microbial Dynamical Systems INference
  Engine for microbiome time-series analyses. *Genome Biology* 17:121.
- Deyle E.R., May R.M., Munch S.B. & Sugihara G. (2016) Tracking and
  forecasting ecosystem interactions in real time. *Proc. R. Soc. B*
  283:20152258.
- Domínguez-García V., Dakos V. & Kéfi S. (2019) Unveiling dimensions of
  stability in complex ecological networks. *PNAS* 116:25714–25720.
- Stein R.R. et al. (2013) Ecological modeling from time-series
  inference: insight into dynamics and stability of intestinal
  microbiota. *PLoS Comput Biol* 9:e1003388.
