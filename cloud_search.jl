include("convex_algorithms.jl")


function Generate_ℓ_Vector(_x, _ℓ, _k, _ρ)
    return (rand(length(_x)) .- 0.5) * _ℓ
end


function Search(local_sol, _f, _ℓ, _γ, _η, _ρ)

    ℓ_start = _ℓ

    for k in 1:_ρ


        _ℓ = ℓ_start
        vec = Generate_ℓ_Vector(local_sol[1], _ℓ, k, _ρ)

        while _f(local_sol[1] + vec) >= local_sol[2] + _η && _f(local_sol[1] - vec) >= local_sol[2] + _η && _ℓ > _η
            #_ℓ *= _γ
            vec = Generate_ℓ_Vector(local_sol[1], _ℓ, k, _ρ)


            push!(searchX, (local_sol[1] + vec)[1])
            push!(searchY, (local_sol[1] + vec)[2])
            push!(searchX, (local_sol[1] - vec)[1])
            push!(searchY, (local_sol[1] - vec)[2])
        end


        if _f(local_sol[1] + vec) < local_sol[2] - _η
            return local_sol[1] + vec, _ℓ
        elseif _f(local_sol[1] - vec) < local_sol[2] - _η
            return local_sol[1] - vec, _ℓ
        end
    end

    return local_sol[1], _ℓ
end


function Cloud_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, _ρ, width, maxIt)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
    s = Search(sol, _f, _ℓ, _γ, _η, _ρ)
    x_prev = s[1]
    it = 0

    while abs(norm(_x) - norm(x_prev)) > _η && it < maxIt

        s = Search(sol, _f, _ℓ, _γ, _η, _ρ)
        x_prev = s[1]

        sol = Grad_Descent(_f, x_prev, _α, _β, _η, _κ)
        s = Search(sol, _f, _ℓ, _γ, _η, _ρ)
        _x = s[1]

        it +=1
    end

    #sol = Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    return sol
end
