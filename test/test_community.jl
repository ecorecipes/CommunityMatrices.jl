using CommunityMatrices.ResilienceCore: state, present_species, PerturbationOptions, global_mortality_tolerance,
                                         local_mortality_tolerances

@testset "GLVCommunity: ResilienceCore's simulated experiments on Lotka–Volterra" begin
    # independent logistic species, dNᵢ/dt = Nᵢ(bᵢ − dᵢ − bᵢNᵢ/Kᵢ): Nᵢ* = Kᵢ(1 − fᵢdᵢ/bᵢ) under a
    # press fᵢ on mortality, so every experiment has a closed form
    b, K, d = [1.0, 1.0, 1.0], [2.0, 1.0, 0.5], [0.099, 0.199, 0.249]
    a = AlphaMatrix(Diagonal(-b ./ K), b .- d)
    c = GLVCommunity(a, d; N0 = fill(0.1, 3))
    Nstar = K .* (1 .- d ./ b)
    @test state(c) ≈ Nstar && present_species(c) == 1:3
    @test_throws ArgumentError GLVCommunity(a, d[1:2])
    @test_throws ArgumentError GLVCommunity(a, -d)

    opts = PerturbationOptions(sequences = 2)
    kexp = findfirst(k -> (opts.global_start + opts.step * k) * maximum(d ./ b) ≥ 1, 0:2000) - 1
    @test global_mortality_tolerance(c; options = opts) ≈ opts.global_start + opts.step * kexp
    @test local_mortality_tolerances(c; options = opts) ≈
          [opts.step * findfirst(k -> (1 + opts.step * k) * d[i] ≥ b[i], 1:2000) for i in 1:3]

    prof = stability_profile(c; options = opts)
    drop = K .* 0.1 .* d ./ b
    B = sum(Nstar)
    @test length(prof) == 27
    @test prof.Rinf ≈ minimum(b .* Nstar ./ K)
    @test prof.RM ≈ sum(drop) / B && prof.max_RM ≈ maximum(drop) / B
    @test prof.mean_CE == 0 && prof.TE ≈ 2 / 3 && prof.mean_RE ≈ 1 / 3

    # an interacting community: removals cause secondary extinctions, and the Jacobian part of
    # the profile is the Jacobian of the survivors
    A = [-1.0 -0.5 0.0; 0.4 -1.0 -0.5; 0.0 0.6 -1.0]
    g = GLVCommunity(AlphaMatrix(A, [1.0, -0.1, -0.1]), [0.0, 0.2, 0.2])
    p = stability_profile(Xoshiro(1), g; options = opts)
    idx = present_species(g)
    jp = jacobian_profile(linearise(g))
    @test p[keys(jp)] == jp
    @test isapprox(jacobian(g), matrix(jacobian(subcommunity(g.alpha, idx))); rtol = 1e-8)
    @test p.max_CE ≥ 1     # removing the basal species takes its consumers with it
end
