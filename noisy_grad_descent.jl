function Ave_Grad(_f, _x, _ℓ, _ρ)
    sum = Compute_Gradient(_f, _x)

    for i = 1:_ρ
        vec = (rand(length(_x)) .- 0.5) * _ℓ + _x
        sum += Compute_Gradient(_f, vec)

        push!(noiseX, vec[1])
        push!(noiseY, vec[2])
    end

    return sum/(_ρ + 1)
end


function Ave_Grad_Descent(_f, _x, _α, _β, _ϵ, _η, _κ, _ℓ, _ρ, maxIt)
    it = 0;
    t = 1
    ∇f = Compute_Gradient(_f, _x)
    normGrad = norm(∇f)

    val = _f(_x)

    push!(xPlot, _x[1])
    push!(yPlot, _x[2])

    while normGrad > _η && it < maxIt
        ∇f = Compute_Gradient(_f, _x)
        ∇f_ave = Ave_Grad(_f, _x, _ℓ, _ρ)
        normGrad = norm(∇f)

        bls = Back_Line_Search(_x, _f, val, ∇f, -∇f_ave, _α, _β, _κ)
        t = bls[1]
        _x -= t*∇f_ave
        val = _f(_x)

        push!(xPlot, _x[1])
        push!(yPlot, _x[2])

        it += 1
    end

    println("\nIterations: ", it)

    push!(solPlotX, _x[1])
    push!(solPlotY, _x[2])

    sol = Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    return sol
end
