# `jacobian_metrics` now computes every metric with ResilienceCore. These are the values the
# plain-matrix functions of EcologicalResilience.jl (commit e84f286), from which this package
# was split, gave on the same Jacobians: the two implementations agree.
const ER_REFERENCE = Dict(
    "glv4" => (Rinf = 0.013821682180493894, R0 = 0.03310836581127131, Is = 0.010475710823308523,
               Amax = 1.0573755236218165, tmax = 3.5, S = 147.171203924115, absumS = 156.24419071475407,
               MRinf = 0.07411093997251655, MR0 = -1.2570277506921945, MAmax = 1.0, Mtmax = 0.0),
    "nonnormal" => (Rinf = 1.0, R0 = 3.524937810560445, Is = 0.057053391992240785, Amax = 2.5601828481660673,
                    tmax = 0.7, S = 6.5, absumS = 6.5, MRinf = 0.7802163735066893, MR0 = -1.2,
                    MAmax = 1.2187644298044127, Mtmax = 0.6),
    "oscillating" => (Rinf = 0.09999999999999998, R0 = 0.22015621187164242, Is = 0.06941129752148878,
                      Amax = 1.2840903521672857, tmax = 1.8, S = 4.042553191489361, absumS = 4.042553191489361,
                      MRinf = 0.09696156596864723, MR0 = -0.21999999999999997, MAmax = 1.1066054101567444,
                      Mtmax = 2.1),
    "stein_ss" => (Rinf = 0.02362480949010252, R0 = 0.778255751277281, Is = 0.0027256391681037026,
                   Amax = 2.9106718760558734, tmax = 5.5, S = 371.53399883166554, absumS = 419.41915409490286,
                   MRinf = 0.06812179764816151, MR0 = -0.8051535632153662, MAmax = 1.0, Mtmax = 0.0),
    "smap" => (Rinf = 3.2642925842401536e-14, R0 = 1.286882965681846, Is = 5.483593572106729e-15,
               Amax = 4.866452175757379, tmax = 2.9, S = 23.871993729905228, absumS = 79.98753845675598,
               MRinf = NaN, MR0 = -0.6118286425526492, MAmax = 2.428224863529894, Mtmax = 3.2))

# rtol for ordinary values; atol for values that are numerically zero (the S-map Jacobian has a
# nearly undamped mode, and its Rinf of order 1e-14 differs between platforms)
same(x, y) = (isnan(x) && isnan(y)) || isapprox(x, y; rtol = 1e-9, atol = 1e-12)

@testset "jacobian_metrics agrees with EcologicalResilience's plain-matrix metrics" begin
    A = [-1.0 -0.3 0.1 -0.2; -0.2 -1.0 -0.4 0.0; 0.3 -0.1 -1.0 -0.25; -0.1 0.2 -0.3 -1.0]
    st = read_stein(joinpath(@__DIR__, "..", "data", "stein2013"))
    ss = steady_state(st.alpha, fill(1e-3, nspecies(st.alpha)))
    X = [exp(-0.1t) * [1.0, 0.5, 0.25][i] + 0.01sin(t + i) for t in 0:40, i in 1:3]
    cases = Dict("glv4" => jacobian(AlphaMatrix(A, [1.0, 0.8, 1.2, 0.9])),
                 "nonnormal" => Jacobian([-1.0 10.0; 0.0 -2.0], [1.0, 0.5]),
                 "oscillating" => Jacobian([0.1 1.0; -0.5 -0.3], [1.0, 2.0]),
                 "stein_ss" => jacobian(subcommunity(st.alpha, ss.survivors)),
                 "smap" => smap_jacobians(X; theta = 1.0, times = [20])[1])
    for (name, j) in cases
        m, ref = jacobian_metrics(j), ER_REFERENCE[name]
        @test keys(m) == keys(ref)
        for k in keys(ref)
            @test same(m[k], ref[k])
        end
    end
    # the one change: an unstable equilibrium has unbounded amplification (the old scan stopped
    # at its time cap with a finite value)
    m = jacobian_metrics(Jacobian([0.1 0.0; 0.0 -1.0], [1.0, 1.0]))
    @test m.Rinf < 0 && m.Amax == Inf && m.tmax == Inf && isnan(m.MRinf)
end
